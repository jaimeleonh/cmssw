"""Source-level invariants of the graph interpreter that a compiler will not check.

Each of these guards a mistake that compiles cleanly and is wrong at run time.
The first one was made while adding warp mode: a bulk rewrite of
`alpaka::syncBlockThreads(acc)` into `barrier(acc)` also rewrote the body of
`barrier` itself, into a call to itself.  Infinite recursion with no side
effects is undefined behaviour, so nvcc deleted it -- the block kernel lost 94
of its 96 barriers, raced, and still built and ran.  The GPU comparison caught
it; nothing here would have.
"""

import pathlib
import re

GRAPH = pathlib.Path(__file__).resolve().parents[1] / "interface" / "alpaka" / "graph"
LAYERS = (GRAPH / "GraphLayers.h").read_text()
KERNEL = (GRAPH / "GraphKernel.h").read_text()


def definitions(text, name):
    """(signature, body) of each non-template-argument overload of `name`."""
    pattern = re.compile(r"inline [\w:<>&\s]*?\b" + name + r"\(([^)]*)\)\s*\{(.*?)\n  \}", re.S)
    return [(m.group(1).strip(), m.group(2)) for m in pattern.finditer(text)]


def test_the_lane_and_barrier_helpers_do_not_call_themselves():
    for name in ("laneOf", "laneCount", "inner", "barrier"):
        found = definitions(LAYERS, name)
        assert found, f"nn::{name} is not defined in GraphLayers.h"
        for signature, body in found:
            assert not re.search(r"\b" + name + r"\(", body), (
                f"nn::{name}({signature}) calls itself; recursion with no side effects is "
                "undefined behaviour and the compiler is free to delete it")


def test_the_block_barrier_is_a_real_block_barrier():
    generic = [body for signature, body in definitions(LAYERS, "barrier")
               if "SubGroup" not in signature]
    assert len(generic) == 1, "expected exactly one generic nn::barrier"
    assert "alpaka::syncBlockThreads(acc)" in generic[0]


def test_every_barrier_in_the_interpreter_goes_through_nn_barrier():
    """A direct block barrier in op code would deadlock warp mode: the warps of
    a block are different jets with different control flow."""
    assert "alpaka::syncBlockThreads" not in KERNEL
    # In GraphLayers.h the only block barriers are inside the two helpers.
    outside = LAYERS
    for signature, body in definitions(LAYERS, "barrier"):
        outside = outside.replace(body, "")
    assert "alpaka::syncBlockThreads" not in outside


def test_math_calls_are_given_the_real_accelerator():
    """A SubGroup is not an Alpaka accelerator; the math functions need inner(acc)."""
    for text in (LAYERS, KERNEL):
        assert not re.search(r"alpaka::math::\w+\(\s*acc\s*,", text)


def test_the_fp32_gemm_tiles_stay_out_of_line():
    """Inlining them pushed nvcc's inliner over budget and it outlined ~57 other
    functions into ParT's hot path, costing it ~3% for code it never runs."""
    assert re.search(r"__attribute__\(\(noinline\)\)\s+void\s+gemmPlainTile\(", LAYERS)
