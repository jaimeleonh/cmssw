"""The exporter and the interpreter must agree on the file format.

scripts/nnexport/format.py is a hand-written mirror of
interface/alpaka/graph/ModelFormat.h.  A divergence between them is the one
bug the format cannot survive -- it would not fail to build, and it would not
fail to load; it would quietly evaluate the wrong graph.  So the header is
parsed here and compared with the mirror, field by field.
"""

import pathlib
import re
import sys

import pytest

PACKAGE = pathlib.Path(__file__).parents[1]
sys.path.insert(0, str(PACKAGE / "scripts"))

from nnexport import format as fmt  # noqa: E402

HEADER = (PACKAGE / "interface/alpaka/graph/ModelFormat.h").read_text()
SHARED = (PACKAGE / "interface/alpaka/graph/GraphShared.h").read_text()


def enum_values(text, name):
    """The `name = value` pairs of one C++ enum body.

    Comments are stripped from the whole body before splitting, because the
    prose in this header contains commas.
    """
    body = re.search(rf"enum(?: class)? {name}\s*:?[^{{]*\{{(.*?)\}};", text, re.S)
    assert body is not None, f"enum {name} not found"
    stripped = re.sub(r"//[^\n]*", "", body.group(1))
    values = {}
    running = 0
    for entry in stripped.split(","):
        entry = entry.strip()
        if not entry:
            continue
        match = re.fullmatch(r"(\w+)\s*(?:=\s*(-?\w+))?", entry)
        assert match is not None, f"cannot parse enumerator {entry!r}"
        if match.group(2):
            running = int(match.group(2).rstrip("uU"), 0)
        values[match.group(1)] = running
        running += 1
    return values


def constant(text, name):
    """The value of an `inline constexpr`, following one level of #define.

    Some of these constants are written as an overridable macro so a build can
    be sized to a particular model; the invariant they carry still has to be
    checked, so resolve the macro rather than skip the constant."""
    match = re.search(rf"inline constexpr [\w:]+ {name} = (-?[0-9a-fA-FxX]+)[uU]*;", text)
    if match is not None:
        return int(match.group(1), 0)

    match = re.search(rf"inline constexpr [\w:]+ {name} = (\w+);", text)
    assert match is not None, f"constant {name} not found"
    macro = re.search(rf"#define {match.group(1)} (-?[0-9a-fA-FxX]+)[uU]*\b", text)
    assert macro is not None, f"constant {name} is {match.group(1)}, which is not #defined here"
    return int(macro.group(1), 0)


def test_magic_and_version_match():
    assert constant(HEADER, "kFormatMagic") == fmt.MAGIC
    assert constant(HEADER, "kFormatVersion") == fmt.VERSION


def test_record_arities_match():
    assert constant(HEADER, "kMaxOpInputs") == fmt.MAX_OP_INPUTS
    assert constant(HEADER, "kMaxOpTensors") == fmt.MAX_OP_TENSORS
    assert constant(HEADER, "kMaxOpInts") == fmt.MAX_OP_INTS
    assert constant(HEADER, "kMaxOpFloats") == fmt.MAX_OP_FLOATS


def test_record_sizes_match_the_static_asserts():
    # The header asserts these at compile time; the mirror asserts them at
    # import time.  Both numbers appear here so a change to either is caught.
    assert fmt.OP_STRUCT.size == 64
    assert fmt.BUFFER_STRUCT.size == 24
    assert fmt.HEADER_STRUCT.size == 128
    assert "sizeof(OpRecord) == 64" in HEADER
    assert "sizeof(BufferRecord) == 24" in HEADER
    assert "sizeof(FileHeader) == 128" in HEADER


def test_opcodes_match():
    cpp = enum_values(HEADER, "OpCode")
    cpp.pop("kOpCodeCount", None)
    expected = {
        "Nop": fmt.Op.NOP, "TileBegin": fmt.Op.TILE_BEGIN, "TileEnd": fmt.Op.TILE_END,
        "Linear": fmt.Op.LINEAR, "LayerNorm": fmt.Op.LAYER_NORM, "BatchNorm": fmt.Op.BATCH_NORM,
        "Activate": fmt.Op.ACTIVATE, "Residual": fmt.Op.RESIDUAL, "Quantize": fmt.Op.QUANTIZE,
        "RoundToGrid": fmt.Op.ROUND_TO_GRID, "Attention": fmt.Op.ATTENTION,
        "HeadScale": fmt.Op.HEAD_SCALE, "PairFeatures": fmt.Op.PAIR_FEATURES,
        "PairScatter": fmt.Op.PAIR_SCATTER, "LoadConst": fmt.Op.LOAD_CONST,
        "Copy": fmt.Op.COPY, "Softmax": fmt.Op.SOFTMAX, "LoadInput": fmt.Op.LOAD_INPUT,
        "Output": fmt.Op.OUTPUT, "Reduce": fmt.Op.REDUCE,
    }
    assert set(cpp) == set(expected), "the two sides list different opcodes"
    for name, value in expected.items():
        assert cpp[name] == int(value), f"opcode {name} differs"


def test_activations_match():
    cpp = enum_values(HEADER, "Activation")
    cpp.pop("kActivationCount", None)
    expected = {"Identity": fmt.Act.IDENTITY, "Gelu": fmt.Act.GELU, "Relu": fmt.Act.RELU,
                "Silu": fmt.Act.SILU, "Tanh": fmt.Act.TANH, "Sigmoid": fmt.Act.SIGMOID}
    assert set(cpp) == set(expected)
    for name, value in expected.items():
        assert cpp[name] == int(value), f"activation {name} differs"


def test_reductions_match():
    cpp = enum_values(HEADER, "Reduction")
    cpp.pop("kReductionCount", None)
    expected = {"Sum": fmt.Red.SUM, "Mean": fmt.Red.MEAN, "Max": fmt.Red.MAX}
    assert set(cpp) == set(expected)
    for name, value in expected.items():
        assert cpp[name] == int(value), f"reduction {name} differs"


def test_row_specs_match():
    cpp = enum_values(HEADER, "RowSpec")
    expected = {"kRowsTokens": fmt.ROWS_TOKENS, "kRowsTokensPlusOne": fmt.ROWS_TOKENS_PLUS_ONE,
                "kRowsPairs": fmt.ROWS_PAIRS, "kRowsActive": fmt.ROWS_ACTIVE,
                "kRowsOne": fmt.ROWS_ONE}
    assert set(cpp) == set(expected)
    for name, value in expected.items():
        assert cpp[name] == value, f"row spec {name} differs"


def test_layout_and_arena_enums_match():
    assert enum_values(HEADER, "Arena") == {"kFloatArena": fmt.ARENA_FLOAT,
                                            "kInt8Arena": fmt.ARENA_INT8}
    assert enum_values(HEADER, "Layout") == {"kRowMajor": fmt.ROW_MAJOR,
                                             "kColumnMajor": fmt.COLUMN_MAJOR}
    from nnexport import graph
    assert enum_values(HEADER, "WeightLayout") == {"kWeightOutIn": graph.WEIGHT_OUT_IN,
                                                   "kWeightInOut": graph.WEIGHT_IN_OUT,
                                                   "kWeightDp4a": graph.WEIGHT_DP4A}


@pytest.mark.parametrize("name,value", [
    ("kFlagInt8Weights", fmt.FLAG_INT8_WEIGHTS),
    ("kFlagInt8Input", fmt.FLAG_INT8_INPUT),
    ("kFlagPairBias", fmt.FLAG_PAIR_BIAS),
    ("kFlagMaskKeys", fmt.FLAG_MASK_KEYS),
    ("kFlagKeyOffsetOne", fmt.FLAG_KEY_OFFSET_ONE),
    ("kFlagTransposeHeads", fmt.FLAG_TRANSPOSE_HEADS),
    ("kFlagMaskRows", fmt.FLAG_MASK_ROWS),
    ("kFlagInPlace", fmt.FLAG_IN_PLACE),
    ("kFlagNormUnbiased", fmt.FLAG_NORM_UNBIASED),
    ("kFlagNormEpsOutside", fmt.FLAG_NORM_EPS_OUTSIDE),
    ("kFlagSoftmaxRows", fmt.FLAG_SOFTMAX_ROWS),
    ("kFlagWeightedRows", fmt.FLAG_WEIGHTED_ROWS),
])
def test_flag_bits_match(name, value):
    match = re.search(rf"inline constexpr std::uint16_t {name} = 1u << (\d+);", HEADER)
    assert match is not None, f"flag {name} not found in the header"
    assert (1 << int(match.group(1))) == value


def test_fnv_matches_a_known_vector():
    # The classic FNV-1a test vector.  The C++ side had a transcribed offset
    # basis that was one digit short, which only showed up as a hash mismatch
    # on a real file; this pins both constants to a published value.
    assert fmt.fnv1a(b"") == 0xCBF29CE484222325
    assert fmt.fnv1a(b"a") == 0xAF63DC4C8601EC8C
    assert fmt.fnv1a(b"foobar") == 0x85944171F73967E8
    assert "14695981039346656037ull" in HEADER
    assert "1099511628211ull" in HEADER


def test_exporter_arena_defaults_match_the_interpreter():
    """A model packed for a bigger arena than the build provides is refused at
    load time, so the exporter's default has to be the build's capacity."""
    import export_model
    assert constant(SHARED, "kArenaFloats") == export_model.ARENA_FLOATS
    assert constant(SHARED, "kArenaInt8") == export_model.ARENA_INT8
