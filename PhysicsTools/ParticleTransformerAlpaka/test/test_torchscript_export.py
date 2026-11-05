"""The TorchScript converter, checked against PyTorch itself.

Models are built and scripted here rather than committed, so the test says
what it converts instead of relying on a fixture nobody can read.  Each one is
converted twice: once in FP32, where agreement with PyTorch should be at the
level of float rounding and anything worse is a conversion bug, and once in
INT8, where a few 1e-4 is the expected quantization error.  Separating the two
is the whole point -- a wrong graph and a coarse quantization look alike if you
only ever measure one number.

Skipped entirely when PyTorch is not importable, which is the usual case
outside a CMSSW area.
"""

import pathlib
import sys

import numpy as np
import pytest

PACKAGE = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(PACKAGE / "scripts"))

torch = pytest.importorskip("torch", reason="PyTorch is needed to write a TorchScript archive")

from nnexport import allocator, reference, writer  # noqa: E402
from nnexport.frontends import torchscript  # noqa: E402
from nnexport.graph import ExportError  # noqa: E402
# The CMSSW torch build has no NumPy bridge; see nnexport/torch_numpy.py.
from nnexport.torch_numpy import to_numpy, to_torch  # noqa: E402

ARENA_FLOATS = 10560
ARENA_INT8 = 4096


def convert(path, *, features, sequence, quantize, max_tokens=16):
    # The arena sizes are what lets the frontend tile a model that would not
    # otherwise fit, so the exporter's own numbers are passed through here too.
    graph, meta = torchscript.build(path, max_tokens=max_tokens, input_features=features,
                                    sequence=sequence, quantize=quantize,
                                    arena_floats=ARENA_FLOATS, arena_int8=ARENA_INT8)
    placement = allocator.allocate(graph, float_capacity=ARENA_FLOATS, int8_capacity=ARENA_INT8)
    image = writer.serialise(graph, placement, arena_floats=ARENA_FLOATS, arena_int8=ARENA_INT8)
    return reference.load(image), meta


def worst_difference(module, model, *, features, sequence, trials=16, seed=99):
    module.eval()
    interpreter = reference.Interpreter(model)
    tokens = model.header["max_tokens"]
    generator = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        multiplicity = int(generator.integers(1, tokens + 1)) if sequence else 1
        live = generator.normal(0.0, 1.0, (multiplicity, features)).astype(np.float32)
        staged = np.zeros((tokens, features), np.float32)
        staged[:multiplicity] = live
        mask = np.zeros(tokens, np.float32)
        mask[:multiplicity] = 1.0
        with torch.no_grad():
            expected = to_numpy(module(to_torch(torch, live))).reshape(-1)
        actual = interpreter.run(staged, np.zeros((tokens, 4), np.float32), mask)
        worst = max(worst, float(np.max(np.abs(actual[: expected.size] - expected))))
    return worst


def save(module, example, tmp_path, name):
    path = tmp_path / f"{name}.pt"
    module.eval()
    torch.jit.save(torch.jit.trace(module, example), str(path))
    return path


# --------------------------------------------------------------------------- #
# models
# --------------------------------------------------------------------------- #


class DeepSets(torch.nn.Module):
    """Per-particle stack, mean pool over particles, jet-level head."""

    def __init__(self, features=15, hidden=64, classes=2):
        super().__init__()
        self.phi = torch.nn.Sequential(
            torch.nn.Linear(features, hidden), torch.nn.BatchNorm1d(hidden), torch.nn.ReLU(),
            torch.nn.Linear(hidden, hidden), torch.nn.LayerNorm(hidden), torch.nn.GELU(),
        )
        self.rho = torch.nn.Sequential(
            torch.nn.Linear(hidden, hidden), torch.nn.ReLU(), torch.nn.Linear(hidden, classes))

    def forward(self, x):
        y = self.phi(x).mean(dim=0)
        return torch.softmax(self.rho(y.unsqueeze(0)), dim=-1)


class ResidualMLP(torch.nn.Module):
    """A jet-level model with residual sums and an interior LayerNorm."""

    def __init__(self, features=10, hidden=32, classes=3):
        super().__init__()
        self.inp = torch.nn.Linear(features, hidden)
        self.a = torch.nn.Linear(hidden, hidden)
        self.b = torch.nn.Linear(hidden, hidden)
        self.norm = torch.nn.LayerNorm(hidden)
        self.out = torch.nn.Linear(hidden, classes)

    def forward(self, x):
        y = torch.relu(self.inp(x))
        y = y + torch.relu(self.a(y))
        y = self.norm(y + self.b(y))
        return torch.softmax(self.out(y), dim=-1)


class ScalarNorm(torch.nn.Module):
    """`alpha * (x - mean) / (std + eps) + bias`, with alpha and bias scalars.

    Not an `nn.LayerNorm`: it emits mean/std/sub/div/add rather than
    `aten::layer_norm`, and torch.std is unbiased, so the interpreter has to
    know both conventions.  Models written from a transformer tutorial normalize
    this way, including the L1 scouting tau identifier.
    """

    def __init__(self, eps=1e-6):
        super().__init__()
        self.alpha = torch.nn.Parameter(torch.ones(1))
        self.bias = torch.nn.Parameter(torch.zeros(1))
        self.eps = eps

    def forward(self, x):
        mean = x.mean(dim=-1, keepdim=True)
        std = x.std(dim=-1, keepdim=True)
        return self.alpha * (x - mean) / (std + self.eps) + self.bias


class Encoder(torch.nn.Module):
    """A transformer encoder with attention pooling and two jet-level heads.

    Everything the straight-line converter could not do: `nn.MultiheadAttention`
    with a padding mask, a hand-written normalization, a softmax over the
    particles feeding a weighted sum, and a tuple of outputs.
    """

    def __init__(self, features=10, embed=32, heads=4, ff=64):
        super().__init__()
        self.embed = torch.nn.Sequential(torch.nn.LayerNorm(features),
                                         torch.nn.Linear(features, embed), torch.nn.GELU())
        self.attn = torch.nn.MultiheadAttention(embed, heads)
        self.n1, self.n2, self.n3 = ScalarNorm(), ScalarNorm(), ScalarNorm()
        self.ff = torch.nn.Sequential(torch.nn.Linear(embed, ff), torch.nn.ReLU(),
                                      torch.nn.Linear(ff, embed))
        self.score = torch.nn.Sequential(torch.nn.Linear(embed, 16), torch.nn.Tanh(),
                                         torch.nn.Linear(16, 1))
        self.head = torch.nn.Linear(embed, 8)
        self.a, self.b = torch.nn.Linear(8, 1), torch.nn.Linear(8, 1)

    def forward(self, x, pad_mask):
        y = self.embed(x).transpose(1, 0)
        y, _ = self.attn(y, y, y, key_padding_mask=pad_mask, need_weights=False)
        y = y.transpose(1, 0)
        y = y + self.n1(y)
        y = y + self.ff(self.n2(y))
        y = self.n3(y)
        weights = torch.softmax(self.score(y), dim=1)
        pooled = (weights * y).sum(dim=1)
        h = torch.relu(self.head(pooled))
        return torch.sigmoid(self.a(h)), self.b(h)

class MaxPool(torch.nn.Module):
    def __init__(self, features=8, hidden=16, classes=2):
        super().__init__()
        self.phi = torch.nn.Sequential(torch.nn.Linear(features, hidden), torch.nn.SiLU())
        self.out = torch.nn.Linear(hidden, classes)

    def forward(self, x):
        y = torch.amax(self.phi(x), dim=0)
        return torch.softmax(self.out(y.unsqueeze(0)), dim=-1)


# --------------------------------------------------------------------------- #
# tests
# --------------------------------------------------------------------------- #


@pytest.mark.parametrize("factory,features,sequence,name", [
    (DeepSets, 15, True, "deepsets"),
    (ResidualMLP, 10, False, "residual"),
    (MaxPool, 8, True, "maxpool"),
])
def test_float_conversion_is_exact(factory, features, sequence, name, tmp_path):
    """In FP32 the interpreter should reproduce PyTorch to float rounding.

    This is the test that catches a mis-converted graph; the INT8 one below
    cannot, because quantization error would hide it.
    """
    torch.manual_seed(11)
    module = factory()
    example = torch.randn(6 if sequence else 1, features)
    path = save(module, example, tmp_path, name)

    model, meta = convert(path, features=features, sequence=sequence, quantize=False)
    assert meta["precision"] == "float32"
    difference = worst_difference(module, model, features=features, sequence=sequence)
    assert difference < 1e-5, f"FP32 conversion differs from PyTorch by {difference:.3e}"


@pytest.mark.parametrize("factory,features,sequence,name", [
    (DeepSets, 15, True, "deepsets"),
    (ResidualMLP, 10, False, "residual"),
    (MaxPool, 8, True, "maxpool"),
])
def test_int8_conversion_stays_in_the_quantization_band(factory, features, sequence, name,
                                                        tmp_path):
    torch.manual_seed(11)
    module = factory()
    example = torch.randn(6 if sequence else 1, features)
    path = save(module, example, tmp_path, name)

    model, meta = convert(path, features=features, sequence=sequence, quantize=True)
    assert meta["precision"] == "int8_weights_fp32_activations"
    difference = worst_difference(module, model, features=features, sequence=sequence)
    assert difference < 5e-3, f"INT8 conversion differs from PyTorch by {difference:.3e}"


def test_pooling_respects_the_mask(tmp_path):
    """Padding must not reach the pool.

    A mean over 16 padded rows instead of the 3 real ones would still produce
    finite, plausible numbers, so this compares against PyTorch run on the real
    particles alone.
    """
    torch.manual_seed(3)
    module = DeepSets()
    path = save(module, torch.randn(4, 15), tmp_path, "masked")
    model, _ = convert(path, features=15, sequence=True, quantize=False)

    tokens = model.header["max_tokens"]
    generator = np.random.default_rng(7)
    live = generator.normal(0.0, 1.0, (3, 15)).astype(np.float32)
    staged = np.zeros((tokens, 15), np.float32)
    staged[:3] = live
    # Deliberately fill the padding with something large, so a pool that
    # included it could not possibly agree.
    staged[3:] = 50.0
    mask = np.zeros(tokens, np.float32)
    mask[:3] = 1.0

    with torch.no_grad():
        expected = to_numpy(module(to_torch(torch, live))).reshape(-1)
    actual = reference.Interpreter(model).run(staged, np.zeros((tokens, 4), np.float32), mask)
    assert np.max(np.abs(actual[: expected.size] - expected)) < 1e-5


def test_a_particle_transformer_is_handed_to_its_own_frontend(tmp_path):
    """Recognition is by parameter name, so a ParT gets fused attention rather
    than a pile of unsupported matmuls."""
    names = ["blocks.0.attn.in_proj.weight", "cls_token", "norm.weight"]
    assert torchscript.looks_like_particle_transformer(names)
    assert torchscript.looks_like_particle_transformer([n.replace("blocks", "mod.blocks", 1)
                                                        if n.startswith("blocks") else n
                                                        for n in names])
    assert not torchscript.looks_like_particle_transformer(
        ["phi.0.weight", "rho.0.weight", "rho.2.bias"])


def test_an_unsupported_op_is_refused_with_its_location(tmp_path):
    """A refusal has to say what it could not convert, or it is useless."""

    class Unsupported(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.a = torch.nn.Linear(8, 8)

        def forward(self, x):
            # An elementwise product of two activations: not in the op set.
            y = self.a(x)
            return torch.softmax(y * torch.sigmoid(y), dim=-1)

    path = save(Unsupported(), torch.randn(1, 8), tmp_path, "unsupported")
    with pytest.raises(ExportError) as raised:
        convert(path, features=8, sequence=False, quantize=False)
    message = str(raised.value)
    assert "cannot map" in message
    assert "aten::mul" in message


def test_a_model_that_never_pools_is_refused(tmp_path):
    """The producer writes one row per jet, so a per-particle output has to be
    caught at export time with an explanation rather than silently truncated."""

    class PerParticle(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.a = torch.nn.Linear(8, 2)

        def forward(self, x):
            return torch.softmax(self.a(x), dim=-1)

    path = save(PerParticle(), torch.randn(5, 8), tmp_path, "perparticle")
    with pytest.raises(ExportError) as raised:
        convert(path, features=8, sequence=True, quantize=False)
    assert "pooling" in str(raised.value)


def test_gelu_tanh_approximation_is_refused(tmp_path):
    """The kernel implements the exact erf form; the tanh approximation is a
    different function and would be a silent 1e-3 shift."""

    class TanhGelu(torch.nn.Module):
        def __init__(self):
            super().__init__()
            self.a = torch.nn.Linear(8, 8)
            self.g = torch.nn.GELU(approximate="tanh")
            self.out = torch.nn.Linear(8, 2)

        def forward(self, x):
            return torch.softmax(self.out(self.g(self.a(x))), dim=-1)

    path = save(TanhGelu(), torch.randn(1, 8), tmp_path, "tanhgelu")
    with pytest.raises(ExportError) as raised:
        convert(path, features=8, sequence=False, quantize=False)
    assert "erf" in str(raised.value)


# --------------------------------------------------------------------------- #
# the real model
# --------------------------------------------------------------------------- #

REAL_MODEL = (PACKAGE.parents[1] / "L1TriggerScouting/JetTagging/data"
              / "softjet_part_unique_cpu.pt")
PRINTED_JET = PACKAGE / "test/printed_jet.txt"


@pytest.mark.skipif(not REAL_MODEL.exists(), reason=f"{REAL_MODEL.name} is not in this checkout")
def test_frozen_particle_transformer_is_recognised_and_its_conventions_read():
    """The shipped L1 scouting ParT is a *frozen* archive: no state_dict, every
    parameter an unnamed constant.  Both graph conventions have to be read out
    of it rather than assumed -- and for this model both are the opposite of
    the QAT checkpoint's, so a default would be silently wrong."""
    from nnexport.frontends import torchscript_part

    module = torch.jit.load(str(REAL_MODEL), map_location="cpu")
    module.eval()
    assert not module.state_dict(), "this test is about the frozen case"

    recovered = torchscript_part.extract(torch, module)
    assert recovered is not None, "the frozen ParT was not recognised"
    tensors, conventions = recovered

    assert conventions["embed_dim"] == 128
    assert conventions["heads"] == 8
    assert conventions["particle_blocks"] == 8
    assert conventions["class_blocks"] == 2
    # Read from the graph: a GELU does follow the last pair convolution here,
    # and the head-scaling einsum is 'tbhd,h->tbdh'.
    assert conventions["pair_final_gelu"] is True
    assert conventions["transpose_scaled_heads"] is True
    assert conventions["pair_norms_folded"] is True

    # Every tensor the builder needs, and nothing left dangling.
    for name in ("cls_token", "norm.weight", "fc.0.weight",
                 "embed.input_bn.running_var", "pair_embed.embed.0.running_mean",
                 "blocks.7.attn.in_proj.weight", "blocks.7.w_resid",
                 "cls_blocks.1.attn.in_proj.weight", "cls_blocks.1.c_attn"):
        assert name in tensors, f"{name} was not recovered"
    assert tensors["blocks.0.attn.in_proj.weight"].shape == (384, 128)
    assert tensors["cls_blocks.0.attn.in_proj.weight"].shape == (384, 128)
    assert tensors["blocks.0.c_attn"].shape == (8,)
    assert tensors["blocks.0.w_resid"].shape == (128,)


@pytest.mark.skipif(not REAL_MODEL.exists(), reason=f"{REAL_MODEL.name} is not in this checkout")
def test_frozen_particle_transformer_matches_torch_on_a_real_jet():
    """Converted in FP32, the graph must reproduce the module more closely than
    the module reproduces itself under a float32-epsilon nudge.

    That floor is the honest bar for a ten-block transformer: this one
    amplifies a 1e-7 relative input change into ~1e-4 of output, so demanding
    1e-7 agreement would be demanding more than float32 can deliver.
    """
    model, meta = convert(REAL_MODEL, features=15, sequence=True, quantize=False)
    assert meta["particle_blocks"] == 8 and meta["class_blocks"] == 2

    module = torch.jit.load(str(REAL_MODEL), map_location="cpu")
    module.eval()

    sys.path.insert(0, str(PACKAGE / "scripts"))
    from verify_model import read_printed_jet

    features, vectors, mask = read_printed_jet(PRINTED_JET, 16, 15)
    live = int(mask.sum())

    def call(f, v, k):
        with torch.no_grad():
            return float(to_numpy(module(torch.zeros(1, 16, 2),
                                         to_torch(torch, f[None].copy()),
                                         to_torch(torch, v[None].copy()),
                                         to_torch(torch, k[None].copy()))).reshape(-1)[0])

    interpreter = reference.Interpreter(model)
    generator = np.random.default_rng(5)
    worst = 0.0
    floor = 0.0
    for multiplicity in (live, 12, 8, 4, 1):
        f = np.zeros((16, 15), np.float32)
        v = np.zeros((16, 4), np.float32)
        k = np.zeros(16, np.float32)
        pick = [i % live for i in range(multiplicity)]
        f[:multiplicity] = features[pick]
        v[:multiplicity] = vectors[pick]
        k[:multiplicity] = 1.0

        expected = call(f, v, k)
        actual = float(np.asarray(interpreter.run(f, v, k))[0])
        worst = max(worst, abs(actual - expected))

        nudged = call((f * (1 + generator.normal(0, 1e-7, f.shape))).astype(np.float32),
                      (v * (1 + generator.normal(0, 1e-7, v.shape))).astype(np.float32), k)
        floor = max(floor, abs(nudged - expected))

    assert worst < max(1e-4, 3 * floor), (
        f"converted graph differs from torch by {worst:.3e}, against a float32 noise "
        f"floor of {floor:.3e}")


# --------------------------------------------------------------------------- #
# attention, hand-written norms, pooling, several outputs
# --------------------------------------------------------------------------- #


def save_masked(module, tmp_path, name, rows=6, features=10):
    path = tmp_path / f"{name}.pt"
    module.eval()
    example = (torch.randn(1, rows, features), torch.zeros(1, rows))
    torch.jit.save(torch.jit.trace(module, example), str(path))
    return path


def masked_difference(module, model, *, features=10, trials=16, seed=5):
    """Worst disagreement over jets of every multiplicity.

    The module is fed its full padded length with an additive -inf mask, which
    is what a fixed-length conversion reproduces: this model pools with a
    softmax that has no mask, so the padded rows do reach the answer and
    "the same jet" means the same padded tensor.
    """
    module.eval()
    interpreter = reference.Interpreter(model)
    tokens = model.header["max_tokens"]
    generator = np.random.default_rng(seed)
    worst = 0.0
    for _ in range(trials):
        multiplicity = int(generator.integers(1, tokens + 1))
        staged = np.zeros((tokens, features), np.float32)
        staged[:multiplicity] = generator.normal(0.0, 1.0, (multiplicity, features))
        mask = np.zeros(tokens, np.float32)
        mask[:multiplicity] = 1.0
        additive = np.where(mask > 0.5, 0.0, -np.inf).astype(np.float32)
        with torch.no_grad():
            result = module(to_torch(torch, staged[None].copy()),
                            to_torch(torch, additive[None].copy()))
        expected = np.concatenate([to_numpy(value).reshape(-1) for value in result])
        actual = interpreter.run(staged, np.zeros((tokens, 4), np.float32), mask)
        worst = max(worst, float(np.max(np.abs(actual[: expected.size] - expected))))
    return worst


def test_a_transformer_encoder_converts_and_matches_torch(tmp_path):
    """The whole new path at once: fused attention, a std()-based norm,
    attention pooling and a tuple of heads."""
    torch.manual_seed(17)
    path = save_masked(Encoder(), tmp_path, "encoder")
    model, meta = convert(path, features=10, sequence=True, quantize=False)

    assert meta["mask_argument"] is True
    assert meta["padding"] == "fixed length"
    assert meta["outputs"] == [1, 1]
    kinds = [layer["op"] for layer in meta["layers"]]
    assert "attention" in kinds
    assert any(layer["op"] == "reduce" and layer["weighted"] for layer in meta["layers"])
    assert any(layer["op"] == "softmax" and layer["over"] == "token"
               for layer in meta["layers"])
    assert any(layer["op"] == "layer_norm" and layer["unbiased"] for layer in meta["layers"])

    difference = masked_difference(torch.jit.load(str(path)), model)
    assert difference < 1e-5, f"FP32 conversion differs from PyTorch by {difference:.3e}"


def test_the_attention_head_count_is_read_from_the_graph(tmp_path):
    """Heads are not guessed from the score scale when the graph says what they
    are -- getting this wrong mixes channels between heads and is silent."""
    from nnexport.frontends import torchscript_patterns as patterns

    torch.manual_seed(4)
    path = save_masked(Encoder(embed=48, heads=6, ff=64), tmp_path, "sixheads")
    frozen = torch.jit.freeze(torch.jit.load(str(path)).eval())
    inner = frozen.graph
    torch._C._jit_pass_inline(inner)
    matches = patterns.find_attention(inner, patterns.uses_of(inner))
    assert len(matches) == 1
    assert matches[0].heads == 6
    assert matches[0].mask is not None


def test_the_unbiased_and_outside_epsilon_conventions_are_not_defaults(tmp_path):
    """A std()-based norm differs from nn.LayerNorm by 0.4% over 128 channels.

    Converting it as an ordinary LayerNorm would still produce plausible
    numbers, so this pins that both flags actually reach the file.
    """
    torch.manual_seed(21)
    path = save_masked(Encoder(), tmp_path, "norms")
    model, meta = convert(path, features=10, sequence=True, quantize=False)
    norms = [layer for layer in meta["layers"] if layer["op"] == "layer_norm"]
    # One nn.LayerNorm in the embedding, three hand-written ones after it.
    assert sum(1 for n in norms if not n["unbiased"]) == 1
    assert sum(1 for n in norms if n["unbiased"] and n["eps_outside"]) == 3


def test_a_model_that_does_not_fit_is_tiled_rather_than_refused(tmp_path):
    """A 16-token, 512-wide feed-forward state does not fit beside the residual
    stream it came from; the exporter schedules it a few rows at a time."""
    torch.manual_seed(2)
    path = save_masked(Encoder(embed=128, heads=8, ff=512), tmp_path, "wide")
    graph, meta = torchscript.build(path, max_tokens=16, input_features=10, quantize=False,
                                    arena_floats=ARENA_FLOATS, arena_int8=ARENA_INT8)
    assert meta["token_tile"] < 16, "this model only fits tiled"
    placement = allocator.allocate(graph, float_capacity=ARENA_FLOATS,
                                   int8_capacity=ARENA_INT8)
    assert placement.float_elements <= ARENA_FLOATS


# --------------------------------------------------------------------------- #
# the tau identifier
# --------------------------------------------------------------------------- #

TAU_MODEL = (PACKAGE.parents[1] / "L1TriggerScouting/JetTagging/data"
             / "softtauid_sigmoid_col.pt")


@pytest.mark.skipif(not TAU_MODEL.exists(), reason=f"{TAU_MODEL.name} is not in this checkout")
def test_the_tau_identifier_converts_with_the_shape_it_was_written_for():
    """It takes a padding mask, so it is converted at its full padded length.

    That is not a stylistic choice: its pooling softmax has no mask, so the
    padded rows contribute, and the module's own answer for the same particles
    depends on how many rows it was given.
    """
    graph, meta = torchscript.build(TAU_MODEL, max_tokens=16, input_features=10,
                                    quantize=False, arena_floats=ARENA_FLOATS,
                                    arena_int8=ARENA_INT8)
    assert meta["batched"] is True
    assert meta["mask_argument"] is True
    assert meta["padding"] == "fixed length"
    assert meta["outputs"] == [1, 1, 1, 1]
    attention = [layer for layer in meta["layers"] if layer["op"] == "attention"]
    assert len(attention) == 1 and attention[0]["heads"] == 8


@pytest.mark.skipif(not TAU_MODEL.exists(), reason=f"{TAU_MODEL.name} is not in this checkout")
def test_the_tau_identifier_matches_torch_in_float():
    """Converted in FP32 it should reproduce the module to float rounding."""
    model, _ = convert(TAU_MODEL, features=10, sequence=True, quantize=False)
    module = torch.jit.load(str(TAU_MODEL), map_location="cpu")
    difference = masked_difference(module, model, trials=24, seed=31)
    assert difference < 1e-5, f"FP32 conversion differs from PyTorch by {difference:.3e}"


@pytest.mark.skipif(not TAU_MODEL.exists(), reason=f"{TAU_MODEL.name} is not in this checkout")
def test_the_shipped_tau_model_file_is_the_current_conversion():
    """data/softtauid_sigmoid.nnm is FP32 on purpose.

    Weight-only INT8 costs this model a few percent of each output's spread --
    it has no quantization-aware training and narrow heads -- and at 1.4 MB
    there is nothing to save.
    """
    shipped = PACKAGE / "data/softtauid_sigmoid.nnm"
    if not shipped.exists():
        pytest.skip("the converted model is not in this checkout")
    model = reference.load(shipped.read_bytes())
    assert model.header["output_size"] == 4
    module = torch.jit.load(str(TAU_MODEL), map_location="cpu")
    assert masked_difference(module, model, trials=8, seed=77) < 1e-5


# --------------------------------------------------------------------------- #
# inferring the feature width


@pytest.mark.parametrize("features", [8, 15, 32])
def test_the_input_width_is_inferred_from_the_model(tmp_path, features):
    """`--input-features` is discoverable, so it should not have to be supplied.

    A wrong width is refused by the first dense layer, so trying the widths the
    weights suggest identifies the right one exactly.
    """
    torch.manual_seed(5)
    module = DeepSets(features=features)
    path = save(module, torch.randn(4, features), tmp_path, f"width{features}")
    loaded = torch.jit.load(str(path), map_location="cpu")
    loaded.eval()
    assert torchscript.infer_input_features(torch, loaded, max_tokens=16) == features


def test_an_unreadable_input_width_is_reported_rather_than_guessed(tmp_path):
    """Inference must fail loudly; a wrong guess would convert the wrong graph."""

    class NoDenseLayer(torch.nn.Module):
        def forward(self, x):
            return x.sum(0, keepdim=True)

    path = save(NoDenseLayer(), torch.randn(4, 15), tmp_path, "nodense")
    loaded = torch.jit.load(str(path), map_location="cpu")
    loaded.eval()
    with pytest.raises(ExportError):
        torchscript.infer_input_features(torch, loaded, max_tokens=16)
