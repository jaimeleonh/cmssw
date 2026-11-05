"""End-to-end checks of the exporter and the NumPy interpreter.

The expensive part -- agreement with the production kernel -- lives in
test/compile_graph_kernel.cc, which needs a C++ compiler.  What is checked here
is everything that can go wrong before a byte reaches the device: the allocator
handing two live buffers the same memory, the writer disagreeing with the
reader, a shape that was assumed rather than read.
"""

import pathlib
import sys

import numpy as np
import pytest

PACKAGE = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(PACKAGE / "scripts"))

from nnexport import allocator, reference, writer  # noqa: E402
from nnexport.allocator import Interval, _liveness  # noqa: E402
from nnexport.frontends import legacy_generated  # noqa: E402
from nnexport.frontends import particle_transformer as frontend  # noqa: E402
from nnexport.graph import ExportError  # noqa: E402

MODEL = PACKAGE / "data/b_kinadd.nnm"
JET = PACKAGE / "test/printed_jet.txt"

# The value the hand-written scalar kernel produces for the printed jet.  The
# interpreter reaches it through an entirely different route -- a tiled,
# reallocating, op-table walk -- so agreement is a real check, not a tautology.
PRINTED_JET_CLASS0 = 0.228496164


@pytest.fixture(scope="module")
def built():
    """The graph, rebuilt from the shipped quantized model."""
    tensors, scales, constants = legacy_generated.load(
        PACKAGE / "interface/alpaka/GeneratedModel.h",
        PACKAGE / "data/model_weights.int8.bin",
        PACKAGE / "data/model_params.fp32.bin")
    graph, meta = frontend.build(
        tensors, scales, max_particles=16,
        pair_final_gelu=constants["pair_final_gelu"],
        transpose_heads=constants["transpose_scaled_heads"])
    placement = allocator.allocate(graph, float_capacity=10560, int8_capacity=4096)
    return graph, meta, placement


def test_shapes_are_read_from_the_model_not_assumed(built):
    _, meta, _ = built
    assert meta["embed_dim"] == 128
    assert meta["heads"] == 8
    assert meta["head_dim"] == 16
    assert meta["ffn_dim"] == 512
    assert meta["particle_blocks"] == 8
    assert meta["class_blocks"] == 2
    assert meta["input_features"] == 15
    assert meta["classes"] == 2
    assert meta["pair_features"] == 4


def test_allocator_never_overlaps_two_live_buffers(built):
    """The property the whole design rests on.

    A buffer may share memory with another only when their lifetimes are
    disjoint.  Getting this wrong produces plausible-looking numbers, so it is
    checked directly rather than inferred from the output.
    """
    graph, _, _ = built
    live = _liveness(graph)
    roots = [b for b in graph.buffers if b.parent is None and b.id in live]
    for i, first in enumerate(roots):
        for second in roots[i + 1:]:
            if first.arena != second.arena:
                continue
            if not live[first.id].overlaps(live[second.id]):
                continue
            a0, a1 = first.offset, first.offset + first.elements
            b0, b1 = second.offset, second.offset + second.elements
            assert a1 <= b0 or b1 <= a0, (
                f"'{first.name}' [{a0},{a1}) and '{second.name}' [{b0},{b1}) share memory "
                f"while both are live")


def test_views_resolve_inside_their_root(built):
    graph, _, _ = built
    for buffer in graph.buffers:
        if buffer.parent is None:
            continue
        root = buffer.root
        assert buffer.offset >= root.offset
        assert buffer.offset + buffer.cols <= root.offset + root.elements


def test_tiled_regions_are_balanced(built):
    graph, _, _ = built
    from nnexport import format as fmt
    depth = 0
    for op in graph.ops:
        if op.opcode == fmt.Op.TILE_BEGIN:
            depth += 1
            assert depth == 1, "tiled regions must not nest"
        elif op.opcode == fmt.Op.TILE_END:
            depth -= 1
            assert depth >= 0, "a region closed without opening"
    assert depth == 0, "a region was left open"


def test_the_model_fits_with_headroom(built):
    _, _, placement = built
    assert placement.float_elements <= 10560
    assert placement.int8_bytes <= 4096


def test_allocator_reports_the_offending_buffer_when_a_model_does_not_fit(built):
    graph, _, _ = built
    with pytest.raises(ExportError) as raised:
        allocator.allocate(graph, float_capacity=4096, int8_capacity=4096)
    message = str(raised.value)
    assert "does not fit" in message
    assert "'" in message, "the error should name the buffer it could not place"
    assert "--arena-floats" in message or "kArenaFloats" in message


def test_serialised_file_round_trips(built):
    graph, _, placement = built
    image = writer.serialise(graph, placement, arena_floats=10560, arena_int8=4096)
    model = reference.load(image)  # validates both payload hashes
    assert model.header["op_count"] == len(graph.ops)
    assert model.header["buffer_count"] == len(graph.buffers)
    assert model.header["max_tokens"] == 16
    assert model.header["output_size"] == 2


@pytest.mark.skipif(not MODEL.exists(), reason="data/b_kinadd.nnm has not been exported")
def test_shipped_model_reproduces_the_printed_jet():
    sys.path.insert(0, str(PACKAGE / "scripts"))
    from verify_model import read_printed_jet

    model = reference.load(MODEL)
    features, vectors, mask = read_printed_jet(JET, model.header["max_tokens"],
                                               model.header["input_features"])
    output = reference.Interpreter(model).run(features, vectors, mask)
    assert output[0] == pytest.approx(PRINTED_JET_CLASS0, abs=1e-5)
    assert output.sum() == pytest.approx(1.0, abs=1e-6)


@pytest.mark.skipif(not MODEL.exists(), reason="data/b_kinadd.nnm has not been exported")
def test_a_single_particle_jet_runs():
    """The smallest sequence: one pair, one token, every tile a partial one."""
    model = reference.load(MODEL)
    tokens = model.header["max_tokens"]
    features = np.zeros((tokens, model.header["input_features"]), np.float32)
    vectors = np.zeros((tokens, 4), np.float32)
    mask = np.zeros(tokens, np.float32)
    features[0] = 0.5
    vectors[0] = [10.0, 5.0, 3.0, np.sqrt(10.0**2 + 5.0**2 + 3.0**2)]
    mask[0] = 1.0
    output = reference.Interpreter(model).run(features, vectors, mask)
    assert np.all(np.isfinite(output))
    assert output.sum() == pytest.approx(1.0, abs=1e-6)


@pytest.mark.skipif(not MODEL.exists(), reason="data/b_kinadd.nnm has not been exported")
def test_an_empty_jet_produces_zeros():
    model = reference.load(MODEL)
    tokens = model.header["max_tokens"]
    output = reference.Interpreter(model).run(
        np.zeros((tokens, model.header["input_features"]), np.float32),
        np.zeros((tokens, 4), np.float32),
        np.zeros(tokens, np.float32))
    assert np.all(output == 0.0)


@pytest.mark.skipif(not MODEL.exists(), reason="data/b_kinadd.nnm has not been exported")
def test_a_hole_in_the_mask_is_not_the_same_as_a_shorter_jet():
    """`active` is the last unmasked index plus one, so an interior gap and a
    truncated tail are different cases; a masked particle must still be
    excluded from every other token's attention."""
    model = reference.load(MODEL)
    tokens = model.header["max_tokens"]
    rng = np.random.default_rng(12345)
    features = rng.uniform(-2.0, 2.0, (tokens, model.header["input_features"])).astype(np.float32)
    vectors = np.zeros((tokens, 4), np.float32)
    for particle in range(tokens):
        px, py, pz = rng.uniform(-30.0, 30.0, 3)
        vectors[particle] = [px, py, pz, np.sqrt(px * px + py * py + pz * pz)]

    mask = np.zeros(tokens, np.float32)
    mask[:6] = 1.0
    full = reference.Interpreter(model).run(features, vectors, mask)

    holed = mask.copy()
    holed[3] = 0.0
    gapped = reference.Interpreter(model).run(features, vectors, holed)

    assert np.all(np.isfinite(gapped))
    assert gapped.sum() == pytest.approx(1.0, abs=1e-6)
    assert not np.allclose(full, gapped), "masking a particle should change the answer"


def test_pair_features_are_symmetric_and_float32():
    """The upper triangle is only evaluated because the features are exactly
    symmetric under i <-> j, and the diagonal only survives because the
    computation is FP32 throughout."""
    a = np.array([12.0, -4.0, 7.0, np.sqrt(144 + 16 + 49)], np.float32)
    b = np.array([-3.0, 9.0, 2.0, np.sqrt(9 + 81 + 4)], np.float32)
    assert np.array_equal(reference._lorentz_features(a, b),
                          reference._lorentz_features(b, a))
    diagonal = reference._lorentz_features(a, a)
    assert np.all(np.isfinite(diagonal))
    assert diagonal.dtype == np.float32


def test_triangular_indexing_is_a_bijection():
    for size in (1, 2, 5, 16):
        seen = {}
        for i in range(size):
            for j in range(i, size):
                slot = reference._pair_slot(i, j, size)
                assert reference._pair_slot(j, i, size) == slot
                assert slot not in seen, f"slot {slot} used twice at size {size}"
                seen[slot] = (i, j)
                assert reference._pair_from_slot(slot, size) == (i, j)
        assert sorted(seen) == list(range(size * (size + 1) // 2))
