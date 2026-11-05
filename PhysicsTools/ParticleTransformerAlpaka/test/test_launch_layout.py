"""Design contracts of the single-launch, shared-memory ParticleTransformer.

These are pure-Python checks of the indexing arithmetic and of the invariants
the sources are expected to keep.  The numerical agreement with the scalar
reference is covered by test/compile_kernel_int8.cc.
"""

import math
import pathlib

MAX_PARTICLES = 16
EMBED_DIM = 128
QKV_DIM = 384
FFN_DIM = 512
HEADS = 8
PAIR_HIDDEN = 64
TOKEN_BLOCK = 4
FFN_TOKEN_TILE = 8
PAIR_TILE = 32
BIG_SCRATCH = MAX_PARTICLES * QKV_DIM
NORM_SCRATCH = FFN_TOKEN_TILE * EMBED_DIM
QUANT_SCRATCH = FFN_TOKEN_TILE * FFN_DIM  # bytes
MAX_PAIRS = MAX_PARTICLES * (MAX_PARTICLES + 1) // 2


def triangular_base(row, size):
    return row * size - (row * (row - 1)) // 2


def pair_slot(i, j, size):
    low, high = min(i, j), max(i, j)
    return triangular_base(low, size) + (high - low)


def pair_from_slot(slot, size):
    """The kernel's closed-form inverse, including its two correction loops."""
    b = 2.0 * size + 1.0
    row = int((b - math.sqrt(b * b - 8.0 * slot)) * 0.5)
    row = max(0, min(row, size - 1))
    while row > 0 and triangular_base(row, size) > slot:
        row -= 1
    while row < size - 1 and triangular_base(row + 1, size) <= slot:
        row += 1
    return row, row + (slot - triangular_base(row, size))


def test_upper_triangular_round_trip():
    for size in range(1, MAX_PARTICLES + 1):
        pairs = size * (size + 1) // 2
        seen = set()
        for slot in range(pairs):
            i, j = pair_from_slot(slot, size)
            assert 0 <= i <= j < size
            assert pair_slot(i, j, size) == slot
            assert pair_slot(j, i, size) == slot  # the features are symmetric
            seen.add((i, j))
        assert len(seen) == pairs
        assert seen == {(i, j) for i in range(size) for j in range(i, size)}


def test_triangular_storage_is_smaller_than_the_full_matrix():
    assert MAX_PAIRS == 136
    assert MAX_PAIRS < MAX_PARTICLES * MAX_PARTICLES


def test_token_padding():
    """`tokens` pads `active` up to a whole register-blocked tile."""
    for active in range(0, MAX_PARTICLES + 1):
        tokens = min(-(-active // TOKEN_BLOCK) * TOKEN_BLOCK, MAX_PARTICLES)
        assert tokens >= active
        assert tokens % TOKEN_BLOCK == 0 or tokens == MAX_PARTICLES
        assert tokens - active < TOKEN_BLOCK
        assert tokens <= MAX_PARTICLES


def test_scratch_is_large_enough_for_every_consumer():
    assert BIG_SCRATCH >= MAX_PARTICLES * QKV_DIM                       # QKV projection
    assert BIG_SCRATCH >= FFN_TOKEN_TILE * (FFN_DIM + EMBED_DIM)        # feed-forward tile
    assert BIG_SCRATCH >= 2 * PAIR_TILE * PAIR_HIDDEN                   # pair ping-pong
    assert BIG_SCRATCH >= MAX_PARTICLES * EMBED_DIM + FFN_TOKEN_TILE * FFN_DIM
    class_scratch = (MAX_PARTICLES + 1) * 2 * EMBED_DIM + 3 * EMBED_DIM + HEADS * (MAX_PARTICLES + 1)
    assert BIG_SCRATCH >= class_scratch
    assert NORM_SCRATCH >= MAX_PARTICLES * 15                           # embedding inputs


def test_shared_footprint_fits_both_back_ends():
    floats = (
        MAX_PARTICLES * EMBED_DIM      # x
        + NORM_SCRATCH                 # na
        + QUANT_SCRATCH // 4           # quant, INT8
        + BIG_SCRATCH                  # big
        + HEADS * MAX_PAIRS            # pair
        + EMBED_DIM                    # cls
        + EMBED_DIM                    # clsNorm
        + 128 + 128                    # reduceA, reduceB
        + (MAX_PARTICLES + 4) * 2      # mean, rstd
        + MAX_PARTICLES * 4            # p4
        + MAX_PARTICLES                # mask
        + 4                            # logits
    )
    raw = floats * 4 + 3 * 4  # active, tokens, pairs
    total = -(-raw // 16) * 16  # the struct is alignas(16)
    assert raw == 47868
    assert total == 47872
    assert total <= 47 * 1024        # the header's own budget
    assert total <= 47 * 1024        # Alpaka CPU block-shared arena default
    assert total <= 48 * 1024        # CUDA static shared-memory limit


def test_source_contract():
    package = pathlib.Path(__file__).parents[1]
    algo = (package / "plugins/alpaka/ParticleTransformerAlgo.dev.cc").read_text()
    header = (package / "plugins/alpaka/ParticleTransformerAlgo.h").read_text()
    producer = (package / "plugins/alpaka/ParticleTransformerProducer.cc").read_text()
    kernel = (package / "interface/alpaka/ParticleTransformerFastKernel.h").read_text()
    transform = (package / "plugins/alpaka/TransformKernel.dev.cc").read_text()
    config = (package / "python/particleTransformerProducer_cfi.py").read_text()

    # One launch for the whole event, one block per jet, no global workspace.
    assert "for (uint32_t jet = block; jet < numberOfJets; jet += blocks)" in algo
    assert "declareSharedVar<::part::JetShared" in algo
    assert "maxBlocks == 0 ? numberOfJets" in algo
    assert "kThreadsPerJet" in algo
    assert "workspace" not in algo.lower()
    assert "workspace" not in producer.lower()
    assert not (package / "interface/alpaka/WorkspaceLayout.h").exists()

    # The API revision guard must agree between the header and the source.
    assert "PARTICLE_TRANSFORMER_ALGO_API_REVISION 5" in header
    assert "PARTICLE_TRANSFORMER_ALGO_API_REVISION != 5" in algo

    # The model is transposed on the host and uploaded once per device.
    assert "transposeWeights" in producer
    assert "std::int8_t[]" in producer  # the weight blob is uploaded as INT8
    assert "modelRegistry" in producer
    assert "std::lock_guard<std::mutex>" in producer
    assert "entry.device == device" in producer
    assert "alpaka::wait(queue)" in producer

    # No barrier may sit inside a divergent loop in the input transform.
    assert "independent_group_elements(acc" not in transform  # only named in a comment
    assert "bunchCrossingOf" in transform
    assert "for (uint32_t jet = block; jet < numberOfJets; jet += blocks)" in transform

    # Only the surviving tokens and the upper triangle are evaluated.
    # Activations feeding the heavy GEMMs are token-contiguous and vector-loaded.
    assert "struct alignas(16) Float4" in kernel
    assert "struct alignas(16) Int4" in kernel
    # nvcc derives the execution space of a lambda and rejects an explicit one.
    assert "ALPAKA_FN_ACC {" not in kernel
    assert "reinterpret_cast<Int4 const*>(inputWords" in kernel
    assert "alignas(16) float na[kNormScratch]" in kernel

    assert "shared.tokens" in kernel
    assert "kMaxPairs" in kernel
    assert "-3.4e38f" in kernel

    # Integer arithmetic: dp4a on INT8 activations, INT32 accumulators, and
    # both scales applied once per output.
    assert "__dp4a(a, b, c)" in kernel
    assert "int sum[TB][CB];" in kernel
    assert "(activationScale * scale[channel]) + bias[channel]" in kernel
    assert "static_assert(KDIM % 16 == 0" in kernel
    assert kernel.count("gemmInt8<") >= 12
    # Only the two odd-width projections keep the FP32 path.
    assert kernel.count("gemmTiled<") == 2

    assert "maxBlocks=cms.uint32(0)" in config
    assert "batchSize" not in config
    assert "synchronizeForTiming=cms.bool(False)" in config
