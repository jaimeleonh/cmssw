# The graph format: running a model CMSSW was not built for

The original package compiled one network into the release. `GeneratedModel.h`
held the tensor descriptors, `ModelConstants.h` held the shapes, and
`ParticleTransformerFastKernel.h` hard-coded the layer sequence. Retraining
meant re-running the exporter and rebuilding; changing the width or the depth
meant editing C++.

This replaces all three with a data file. CMSSW ships an interpreter that knows
how to evaluate a set of layers; a model file says which of them to run, over
what shapes, and where to put the intermediate results. A new network is a new
file.

The original path is still present and unchanged — `ParticleTransformerProducer`
and its kernel are untouched — so the two can be compared on the same event.

## What lives where

```
interface/alpaka/graph/
  ModelFormat.h    the on-disk contract: opcodes, records, flags.  No Alpaka,
                   no CMSSW, no C++ library: host, device and tools all include it
  GraphShared.h    the fixed shared-memory arena the interpreter provides
  GraphLayers.h    block-cooperative Linear, LayerNorm, BatchNorm, attention,
                   pair features -- the algorithms, with shapes as arguments
  GraphKernel.h    the interpreter: buffer addressing, tiled regions, op dispatch
  ModelLoader.h    host-side validation; nothing reaches the device unchecked

plugins/alpaka/
  GraphInferenceAlgo.{h,dev.cc}   one launch per event, one block per jet
  GraphInferenceProducer.cc       loads a .nnm, uploads once per device

scripts/
  export_model.py       CLI: checkpoint -> .nnm
  export_torchscript.py CLI: TorchScript -> .nnm, with --check against torch
  verify_model.py       runs an exported file and checks its output
  nnexport/
    format.py           mirror of ModelFormat.h (kept honest by a test)
    graph.py            the IR: buffers, views, ops, tiled regions
    allocator.py        liveness analysis and arena packing
    writer.py           serialisation
    reference.py        a NumPy interpreter of an exported file
    quantize.py         symmetric INT8, per output channel
    frontends/
      particle_transformer.py   Weaver ParT -> IR
      legacy_generated.py       an already-exported model -> IR (migration)
      torchscript.py            a scripted/traced torch module -> IR
```

`test/run_graph.cc` builds a small command-line tool that runs any `.nnm`
through the C++ interpreter and prints its outputs. It is how a converted model
gets compared against the framework it came from without a GPU, and the
quickest way to see what a model file actually does.

## The file

One blob, every section 16-byte aligned and in device order, so the loader
hands back pointers into the same bytes the producer uploads:

```
[FileHeader]      magic, version, section extents, payload hashes, limits
[OpRecord   x N]  the program, in execution order
[BufferRecord x M] the arena layout, precomputed by the exporter
[int8 weights]    dense tensors, already in the layout their kernel reads
[float params]    scales, biases, normalisation parameters, FP32 weights
[strings]         op names, diagnostics only
```

Nothing in the file is a pointer. A tensor is an element offset into a blob and
a buffer is an index into the buffer table, so the same bytes are valid on the
host and on every device.

## The op set

Twenty opcodes, listed with their argument layouts in `ModelFormat.h`:
`Linear`, `LayerNorm`, `BatchNorm`, `Activate`, `Residual`, `Quantize`,
`RoundToGrid`, `Attention`, `HeadScale`, `PairFeatures`, `PairScatter`,
`Reduce`, `Copy`, `Softmax`, `LoadInput`, `LoadConst`, `Output`, `Nop`, and the
region markers `TileBegin`/`TileEnd`.

Two are worth singling out. `Attention` is fused rather than decomposed into
matmul and softmax, because the alternative is materialising a
[tokens][tokens][heads] score tensor in shared memory; one lane owns one
(query, head), keeps the score row in registers, and writes the result back
over its own slice of Q. `Reduce` (sum, mean or max over the row axis) is what
turns a per-particle stack into a jet-level one — the only op whose output has
fewer rows than its input, and it skips masked rows so padding cannot reach the
pool and a mean divides by the number of real particles.

Several ops carry conventions in their flags rather than in separate opcodes,
because the difference is a line of arithmetic and getting it wrong is silent:

| flag | op | what it changes |
|---|---|---|
| `kFlagNormUnbiased` | `LayerNorm` | divide the sum of squares by `cols - 1`, which is what `torch.std()` does by default — over 128 channels that is 0.4% |
| `kFlagNormEpsOutside` | `LayerNorm` | scale by `1 / (sqrt(var) + eps)` instead of `1 / sqrt(var + eps)` |
| `kFlagSoftmaxRows` | `Softmax` | normalize down each column over the particles, which is the softmax an attention-pooling layer takes |
| `kFlagWeightedRows` | `Reduce` | multiply row *r* by `in[1][r][0]` first, so `(softmax(score) * x).sum(particles)` is one pass |

## Three ideas carry the design

### Buffers are named, not allocated

A frontend says "I need somewhere to put 16 tokens of 384 channels". Where that
lands is the allocator's problem. A **view** is a buffer record that aliases
another with its own offset, columns and stride, which costs nothing and makes
aliasing explicit: slicing a fused QKV projection into Q, K and V is three
views, and a transformer block reusing its dead K columns for the output
projection is another.

### Tiling is part of the program

A 16-token, 512-channel feed-forward state is 32 kB and does not fit beside
everything else it has to coexist with. The hand-written kernel solved this by
tiling the feed-forward loop to 8 tokens; the IR makes that a construct:

```python
with graph.tile(fmt.ROWS_TOKENS, 8):
    ...   # buffers marked tiled hold 8 rows and are indexed from 0
          # buffers not marked tiled span the sequence and are indexed
          # from the tile base
```

The interpreter runs the region once per tile. Without this the allocator would
simply report that the model does not fit.

`nnexport/tiling.py` does this automatically for a converted model. It finds
maximal runs of ops whose output row *r* depends only on input row *r*, works
out which buffers are born and die inside such a run, and shrinks exactly
those. Two cuts keep a region from costing more than it saves: a region pins
every buffer that crosses its boundary for the whole region, so a run is split
after the last read of something older than it (the fused QKV projection,
read once by the output projection and then dead) and before the write of
something that outlives it. The exporter then tries 16, 8, 4, 2 and 1 rows at
a time and keeps the largest that fits.

### The arena is fixed, so the allocator runs at export time

Block shared memory is sized when the kernel is compiled, so the interpreter
cannot grow an arena to fit whatever it is handed. It publishes one
(`kArenaFloats`, `kArenaInt8` in `GraphShared.h`) and the exporter packs into
it by liveness analysis — two buffers may share storage when their lifetimes
are disjoint — failing at export time, naming the buffer it could not place:

```
export failed: the model does not fit in the float arena: placing buffer
'blocks.0.qkv' (6144 floats) at 8192 needs 14336 floats but the interpreter
provides 10560.
Reduce a tile size in the frontend, or raise kArenaFloats in GraphShared.h
and rebuild.
```

A model file therefore either runs or never gets written. There is no runtime
allocation failure to handle, and the interpreter does no bounds checking in
the inner loop.

## What is still compile-time, and why

The register tile of a GEMM. Each dense layer picks a `TB x CB` accumulator
tile from its runtime shapes and jumps into one of twelve instantiations, so
the accumulator stays in registers and the innermost loops still unroll. The
same is true of attention: the score row is kept in registers, with two
instantiations so a 17-key model does not pay for the longest sequence the
interpreter allows.

What was lost relative to the hand-written kernel is that `KDIM` and `NOUT` are
no longer compile-time, so the k loop has a runtime bound. The 16-byte vector
load of activations survives behind a runtime width check. On an L4 that costs
**9%** against the compiled kernel — see "What the interpreter costs on a GPU"
below for the measurement and what it is made of.

## Adding a layer type

1. Add an opcode at the end of `OpCode` in `ModelFormat.h`, with a comment
   saying what each of `in`, `out`, `rows`, `cols`, `i[]`, `t[]` and `f[]`
   means for it. Appending does not need a format version bump: the loader
   rejects an opcode it does not know, which is what you want for a model file
   that uses a layer an older release cannot evaluate.
2. Mirror it in `scripts/nnexport/format.py`. `test_model_format.py` compares
   the two and fails if they drift.
3. Implement it in `GraphLayers.h` if it is a real algorithm, or inline in
   `executeOp` if it is elementwise. It is called by every lane of the block
   and must end with a barrier.
4. Add an emitter to `Graph` in `graph.py`, and a `_op_<name>` handler to
   `reference.py`.
5. Add whatever the loader needs to check about it.

The NumPy interpreter is the specification: where it and the kernel could
differ, it does what the kernel does, not what NumPy would find natural.

## Converting a TorchScript model

`scripts/export_torchscript.py` takes a module saved with `torch.jit.save` and
writes a `.nnm`. A TorchScript archive carries the code as well as the weights,
so the layer sequence is read from the model rather than assumed:

```bash
python3 scripts/export_torchscript.py model.pt --output data/model.nnm \
    --max-particles 16 --check      # --input-features is inferred when omitted
```

It covers straight-line stacks — `Linear`, `LayerNorm`, `BatchNorm1d`,
ReLU/GELU/SiLU/Tanh/Sigmoid, residual sums, pooling over the particle axis and
softmax along either axis — plus three idioms that TorchScript leaves as a
spray of primitives and `torchscript_patterns.py` puts back together:

* **`nn.MultiheadAttention`**, which becomes the fused `Attention` op. Both
  spellings are read: the single `aten::scaled_dot_product_attention` recent
  torch emits, and the `baddbmm`/`softmax`/`bmm` decomposition older releases
  (and the archives already in this repository) emit. The head count comes from
  the `embed_dim // num_heads` the graph computes, falling back on the score
  scale; a `key_padding_mask` is recognised and dropped, since the interpreter
  masks from the jet's own mask.
* **a normalization written by hand** as
  `alpha * (x - x.mean(-1)) / (x.std(-1) + eps) + bias`, which emits no
  `aten::layer_norm` at all but is one — with the two convention flags above.
* **attention pooling**, `(softmax(score, particles) * x).sum(particles)`,
  which becomes one weighted `Reduce`.

That covers a transformer encoder with a pooled head, including one that
returns several quantities at once: a tuple of jet-level tensors becomes
consecutive outputs in the order the model lists them.

A Weaver ParticleTransformer is still routed to `particle_transformer.py`
instead — its pairwise interaction bias is not something to rediscover from
`matmul` nodes, and that path is the verified one.

Every activation carries a label per torch axis — particle, channel, or a
width-1 axis that only exists for shape reasons — so a `transpose` or a `view`
is checked rather than guessed at, and a reshape is accepted exactly when the
particle and channel extents come out unchanged. The shape arithmetic torch
writes around attention (`size`, `Int`, `NumToTensor`) is evaluated
symbolically, so the views that split and merge the heads are matched rather
than assumed.

With `--sequence` (the default) the module's `forward()` takes
`[particles, features]` and must pool down to one row before returning. With
`--no-sequence` it takes `[1, features]` of jet-level quantities. Whether
there is also a leading batch axis is not assumed: the exporter calls the
module both ways and uses the one that runs, preferring the unbatched form,
since a model written for `[particles, features]` usually also accepts a batch
axis and then pools over the wrong one.

### Fixed length or the jet's own multiplicity

A module whose `forward()` takes a padding mask was written to be fed a fixed
number of rows, and its answer generally depends on how many. The L1 scouting
tau identifier is the case in point: its pooling softmax has **no** mask, so
the padded rows contribute to the jet-level result. The same 7 particles give
0.164 unpadded and 0.137 padded to 16 — a 16% difference that belongs to the
model, not to the conversion.

So such a model is converted at its full padded length (`--padding fixed`, the
default when a mask argument is present): every per-particle op runs over
`--max-particles` rows and the pooling does not skip padding, which reproduces
the module exactly. Attention still masks the padded keys, which is what the
module's `key_padding_mask` did. A module without a mask argument is converted
at the jet's own multiplicity (`--padding variable`).

A TorchScript export carries no activation ranges, so the dp4a path is not
available: weights are quantized per output channel and activations stay FP32,
which is the weight-only scheme this package used before it had a QAT
checkpoint. `--float` keeps FP32 weights too.

### A frozen ParticleTransformer

A module saved after `torch.jit.freeze` has no `state_dict` at all: every
parameter is an unnamed constant in the graph. `torchscript_part.py` recovers
them by walking the graph in execution order and assigning by position,
asserting the expected shape at every step. Two conventions are *read from the
graph* rather than defaulted, and for the shipped L1 scouting model both are
the opposite of the QAT checkpoint's:

* `pair_final_gelu`, from whether a GELU follows the last pair convolution;
* `transpose_scaled_heads`, from the head-scaling einsum's equation
  (`tbhd,h->tbdh` means the permuted convention).

Getting either wrong is silent and large, so neither has a default here.
Freezing also folds each pair-embedding BatchNorm into the convolution before
it; the builder skips the absent norms.

### Checking a conversion

**`--check` is the part worth using.** It runs the exported graph against the
original module and reports the largest disagreement — and, alongside it, how
far the module's *own* output moves when its input is nudged by one part in
1e7, which is about float32 epsilon.

That second number matters more than it sounds. A ten-block transformer
amplifies input rounding: the shipped L1 scouting ParT turns a 1e-7 relative
input change into ~1.6e-4 of output. **No reimplementation can agree with a
model more closely than the model agrees with itself**, so `--check` treats
three times that floor as a legitimate tolerance and prints the ratio.

Run it twice, and read the ratio rather than the raw number:

| | shallow stack | one-block tau encoder | ten-block ParT |
|---|---|---|---|
| `--check --float` | ~1e-7 | 2.4e-6, i.e. **1.0x** the noise floor | ~2e-5, i.e. **0.1x** the noise floor |
| `--check` (INT8) | a few 1e-4 | a few 1e-2 (see below) | ~2e-3, i.e. 13x the floor |

Separating the two is the whole point: a wrong graph and a coarse quantization
produce similar-looking numbers if you only ever measure one. `--check-jet`
draws real particles from a printed-jet file, which is strongly recommended for
a deep model — random features are outside the training distribution, where
the amplification is an order of magnitude worse and the number stops meaning
anything.

Unsupported nodes are listed with their source location, so a refusal says
which line of the original model to change.

### When weight-only INT8 is the wrong default

The tau identifier is the counter-example to "quantize by default". It had no
quantization-aware training and its heads are narrow — a per-output-channel
scale over 32 weights is coarse — so INT8 weights cost a few percent of each
output's spread:

| output | max &#124;difference&#124; from torch | as a fraction of that output's spread |
|---|---|---|
| tau probability | 0.0029 | 4.6% |
| vertex z | 0.039 | 4.3% |
| pt | 0.037 | 2.8% |
| charge probability | 0.0019 | 2.2% |

(Random inputs, so an over-estimate: prior work here showed features outside
the training distribution inflate the difference by about an order of
magnitude. But there is no real-input sample for this model in the repository
to measure with.) At 1.4 MB in FP32 against 0.37 MB in INT8 there is nothing
worth saving, so `data/softtauid_sigmoid.nnm` is the `--float` conversion.

## Adding a model family

Write a frontend. It takes a checkpoint and returns a `Graph`; nothing
downstream knows where the ops came from. `particle_transformer.py` is the
worked example — note that it *reads* the widths, depths, head count and class
count from the tensors rather than assuming them, so a resized ParT needs no
change to it at all.

## Exporting and checking

```bash
# from a quantization-aware Weaver checkpoint
python3 scripts/export_model.py particle_transformer net_best_epoch_state_q.pt \
    --output data/b_kinadd.nnm --report data/b_kinadd.export.json

# from the already-exported model, needing no checkpoint
python3 scripts/export_model.py legacy_generated interface/alpaka/GeneratedModel.h \
    --legacy-weights data/model_weights.int8.bin \
    --legacy-params  data/model_params.fp32.bin \
    --output data/b_kinadd.nnm --report data/b_kinadd.export.json

# from a TorchScript archive, checked against torch itself
python3 scripts/export_torchscript.py \
    ../../L1TriggerScouting/JetTagging/data/softtauid_sigmoid_col.pt \
    --output data/softtauid_sigmoid.nnm --report data/softtauid_sigmoid.export.json \
    --max-particles 16 --float --check --check-tolerance 1e-5

# run the exported file in NumPy and check its output
python3 scripts/verify_model.py data/b_kinadd.nnm \
    --jet test/printed_jet.txt --expect 0.228496164

# the exporter and the interpreter agree on the format; the allocator never
# overlaps two live buffers
python3 -m pytest test/ -q

# the interpreter agrees with the hand-written kernel, on the printed jet and
# on 32 randomised ones at four block sizes
g++ -O2 -std=c++20 -Itest/stubs -I<src-root> -pthread \
    -o graphtest test/compile_graph_kernel.cc
./graphtest data/b_kinadd.nnm data/model_weights.int8.bin \
            data/model_params.fp32.bin test/printed_jet.txt

# any .nnm through the C++ interpreter, at several block sizes.  This is the
# only check that compares the C++ and NumPy interpreters against each other,
# which is where a disagreement between them shows up
g++ -O2 -std=c++20 -Itest/stubs -I<src-root> -pthread \
    -o run_graph test/run_graph.cc
./run_graph data/softtauid_sigmoid.nnm jets.txt 256
```

The export report records the checkpoint hash, the per-tensor quantization
error, every activation scale, and the full arena allocation with each buffer's
offset, size and live range.

## Measured

Interpreter against the hand-written kernel, same quantized model:

| | class 0 |
|---|---|
| hand-written kernel, printed jet | 0.228496164 |
| graph interpreter, printed jet | 0.228496164 |
| NumPy interpreter, printed jet | 0.228496194 |

Over 129 evaluations — the printed jet plus 32 randomised ones covering every
multiplicity from 1 to 16, half with a hole punched in the mask, at 1, 32, 96
and 256 threads per block — the largest difference from the hand-written kernel
is **3.7e-4**. That is the ordinary INT8 instability: a few 1e-7 before rounding
occasionally moves an activation across a quantization boundary. For scale, the
existing kernel differs from its own scalar reference by up to 1.4e-3 on the
same jets.

Shared memory per jet is 47,840 B against the hand-written kernel's 47,872 B,
so occupancy is unchanged.

For a converted TorchScript DeepSets model, the C++ interpreter and the NumPy
one agree to 6.0e-08, and both differ from PyTorch by 4.3e-04, which is the
INT8 weight-only quantization alone.

The shipped L1 scouting ParticleTransformer
(`L1TriggerScouting/JetTagging/data/softjet_part_unique_cpu.pt`), converted and
run on real particles from `test/printed_jet.txt`:

| | max difference from PyTorch |
|---|---|
| FP32 weights (conversion fidelity) | **2.6e-05** |
| INT8 weights (production default) | 2.7e-04 |
| the model's own float32 noise floor | 1.6e-04 |
| C++ interpreter vs NumPy interpreter | 1.2e-07 |

The FP32 conversion is an order of magnitude *inside* the model's own float32
reproducibility, which is as close as the question can be asked.

The L1 scouting tau identifier
(`L1TriggerScouting/JetTagging/data/softtauid_sigmoid_col.pt`) — a one-block
encoder with `nn.MultiheadAttention`, hand-written normalizations, attention
pooling and four heads — over 32 jets covering every multiplicity from 1 to 16,
each fed as a 16-row padded tensor:

| | max difference from PyTorch |
|---|---|
| NumPy interpreter, FP32 weights | 2.9e-06 |
| C++ interpreter, 256 threads | 3.3e-06 |
| C++ interpreter, 32 threads | 4.8e-06 |
| C++ interpreter, 1 thread | 4.3e-06 |
| the model's own float32 noise floor | 2.4e-06 |
| C++ interpreter vs NumPy interpreter | 5.7e-06 |

Its largest output is ~9, so 3e-6 is 4e-7 relative — the conversion sits at the
model's own float32 noise floor.

Converting it needed a tiled schedule (4 rows at a time) to fit its 512-wide
feed-forward state into the arena, which the exporter found by itself; it uses
10,240 of the 10,560 available floats.

## A numerical trap worth knowing about

The pairwise features must be computed in FP32, in the kernel's order. For the
diagonal pair the invariant mass is `(2E)^2 - (2p)^2`, a cancellation of two
nearly equal squares, so what survives is rounding. Computing it in float64
gives a different residue, the clamp at `1e-8` lands on the other side, and the
diagonal of the attention bias comes out visibly different — enough to move the
classifier output by 1e-4. The first NumPy interpreter had exactly this bug and
it looked like ordinary quantization noise until the intermediates were
compared against the C++ kernel tensor by tensor.

## What the interpreter costs on a GPU

Measured on an NVIDIA L4, 50,000 jets, one block of 256 threads per jet, with
`test/bench_graph.cu` and `test/bench_compiled.cu` — the same jets from the
same seed through both kernels.

**Read this table with the occupancy caveat below; the headline number depends
entirely on how full the jets are.** With multiplicity drawn uniformly from
[1, 16], which is what `bench_graph.cu` generates by default:

| | per 50k jets | vs the compiled kernel |
|---|---|---|
| hand-written compiled kernel | 339 ms | — |
| graph interpreter, before this round | 386 ms | +13.8% |
| graph interpreter, now | 371 ms | **+9.3%** |

### The +9.3% was measured at one occupancy, and both kernels move with it

`bench_graph.cu` draws multiplicity uniformly from [1, 16], so the jets are
half empty on average. That matters because the interpreter sizes every op by
the *live* token count (`resolveRows` -> `shared.tokens`, and `kRowsPairs` is
quadratic in it). The compiled kernel has a weaker version of the same
dependence. `--minmult` raises the floor of the distribution; at 56,052 jets
(the real jets-per-event, see below) on an L4 with `b_kinadd.nnm`:

| min multiplicity | graph | compiled | ratio |
|---|---|---|---|
| 1 | 410 ms | 377 ms | 1.09 |
| 8 | 540 ms | 471 ms | 1.15 |
| 16 | 651 ms | 567 ms | 1.15 |

Both columns rise. The ratio is 1.09-1.15 and fairly stable, so the headline
"+9%" is not badly wrong -- but always say which occupancy a number came from,
and rebuild *both* benchmarks before sweeping. An earlier pass of this file
claimed 1.70x from a sweep in which only `bench_graph` had been rebuilt with
`--minmult`; the stale `bench_compiled` ignored the flag, held its column
artificially flat, and manufactured the result. Both binaries now accept it.

### In the framework the gap is wider than the bench predicts

A full `cmsRun` job over the Phase-2 scouting sample (100 events, one thread,
one stream, `cuda_async`, `batchSize=0`, both producers on the same `b_kinadd`
weights), where inference is ~99% of the event:

| | event loop Real/event | `SoftJetIdSC4` |
|---|---|---|
| `ParticleTransformerProducer` (compiled) | 0.370 s | 0.366 s |
| `GraphInferenceProducer` (graph) | 0.481 s | 0.476 s |

`nsys` pins this on the kernels and nothing else: 56,052 jets per event, one
block of 256 threads each, and the inference kernel is 99.5% of GPU time at
**353.7 ms compiled against 465.8 ms graph, a ratio of 1.32**. Neither producer
carries measurable host-side overhead. So the framework gap (1.32) is wider
than the standalone one (1.15) on the same jet count.

The difference is in code generation, not in the graph. For the same source,
`cuobjdump -res-usage` reports:

| build | GraphJetKernel | JetInferenceKernel |
|---|---|---|
| plain `nvcc` (`test/bench_*.cu`) | REG:128 STACK:136 | REG:128-ish |
| CMSSW plugin | **REG:254 STACK:504** | REG:255 STACK:1904 |

Three candidates were tested and ruled out: CMSSW's CUDA flag set reproduces
REG:128 and the same runtime when applied to the benchmark; removing
`__launch_bounds__` from the benchmark changes nothing (nvcc picks 128 anyway);
and the CUDA build uses no `-rdc`. What is left is the Alpaka
`alpaka::detail::gpuKernel` wrapper, which carries no launch bounds, so ptxas
has no occupancy target to hit and spends registers freely. That is the place
to look next, and it would speed up **both** producers.

Separately, and relevant to the ROCm backend which *does* compile with
`-fgpu-rdc`: relocatable device code is expensive here. At 56,052 jets,
`-rdc=true` takes the graph kernel from 414 to 599 ms and the compiled kernel
from 379 to 513 ms -- 35-45% on both.

On the same 50,000 jets the two kernels are **bit-identical on 74.9%** of them
and agree to better than 1e-6 on 99.4%, with a largest difference of 2.4e-3 and
a mean of 1.8e-6. That tail is the same INT8 instability the CPU comparison
shows: a few 1e-7 before rounding occasionally moves an activation across a
quantization boundary. The GPU agreement is much tighter than the CPU's
because both kernels there take the same fused-multiply-add contractions, so
most jets never reach a boundary at all.

That makes the compiled kernel a usable oracle for interpreter changes — but
the bar is "75% exact, tail to a few 1e-3", not exactness. A real regression
shows up immediately as a difference on most jets.

### Where the gap was, and where it was not

The first hypothesis was barriers: the interpreter ends every op with one, 372
per jet against the hand-written kernel's 29 sites. **That was wrong.** Removing
every barrier — which gives the wrong answer, but bounds the gain — made the
kernel *slower*, 412 ms against 386 ms, and `ncu` puts the barrier stall higher
in the compiled kernel (23.6%) than in the interpreter (20.1%). Warps that
arrive early at a barrier are not idle when sixteen of them are resident.

What the profile actually showed was **+50% instructions and +89% global load
requests** for the same network. Three changes addressed that:

*   **The op table is read by value, not by reference.** A `OpRecord const&`
    leaves each field as its own load from global memory; a 64-byte
    `alignas(16)` copy is four vector loads that every lane broadcasts. Same
    for the `BufferRecord`s. Global load requests fell 10%.
*   **No integer division per element.** An SM has no integer divide unit, so
    `index / cols` with a width read from the file expands to roughly twenty
    instructions, and the elementwise ops were doing it once per element. The
    `Lattice` walk in `GraphLayers.h` steps the row and column forward instead.
*   **Pointers walked forward in both GEMMs** rather than `k * stride + column`
    recomputed per access, which was a multiply-add per weight word — as many
    instructions as the dp4a it feeds.

Together: 58.2M cycles to 54.6M, and 386 ms to 371 ms.

Three further attempts did **not** pay and were reverted: `#pragma unroll` on
the k loop (no change), and specialising the attention on its head width (26
instructions saved, but the extra instantiation cost more in instruction cache
than it gained, and it perturbed the aliasing between the attention output and
its own queries).

### What is left

The remaining 9% is diffuse rather than one thing. The largest single item the
profiler flags is a shared-memory bank conflict in the attention dot product:
one lane owns one (query, head), adjacent lanes are 16 floats apart, and with
32 banks that is a 16-way conflict. **The hand-written kernel has exactly the
same conflict** — it is a property of the 384-wide fused QKV layout, not of the
interpreter — so fixing it means padding the stride at export time, and would
speed up both paths. That is a bigger and more interesting job than closing the
remaining gap.

## Warp layout: one jet per warp, for small networks

A block of 256 lanes per jet is right for ParT, where one jet carries eight
encoder blocks of 128-wide attention. It is wrong for a network like the L1
scouting DeepSets (15 -> 32 -> 32, pool, 32 -> 32 -> 1): each op has 16 x 32
elements of work between block-wide barriers, and `ncu` put the barrier stall
at 6.07 per issue with 33% of warps active. TorchScript, which batches every
particle of the event into one cuBLAS GEMM, was 2x faster.

### What it is

`GraphInferenceProducer` has a `layout` parameter:

| value | meaning |
|---|---|
| `auto` (default) | `warp` if the model fits the per-warp arena **and** the back-end is CUDA, otherwise `block` |
| `block` | one jet per 256-lane block, the full 10,560-float arena |
| `warp` | one jet per warp, eight per block, 1,152 floats each; an error if the model does not fit or the back-end is not CUDA |

Every shipped ParT and the tau model use 10,240-10,304 floats, so `auto`
leaves them in block layout. The DeepSets uses 1,024 and goes to warp layout.
`ModelLoader` measures a model's actual high-water mark (`floatsUsed`), so this
needs no re-export.

The op code is shared, not duplicated. It never asks Alpaka directly which lane
it is, how many there are, or to wait: it calls `nn::laneOf`, `nn::laneCount`,
`nn::barrier` and `nn::inner`, and an `nn::SubGroup` answers them for a warp
(`barrier` becomes `__syncwarp`). The ops are templated on the state type,
`GraphShared` or the small `GraphWarpState`.

**The warp barrier is not an optimisation, it is required for correctness.**
The warps of a block are different jets, with different control flow: an op
with no rows returns before its barrier, and a tiled region runs once per tile
of *that* jet. A block barrier would deadlock as soon as two jets in a block
disagreed. That is also why warp layout is CUDA-only: a ROCm wavefront is 64
lanes and would need its own group size and sync, and the CPU back-ends run
one thread per block.

### The FP32 GEMM had to be tiled too

Warp layout alone took the DeepSets from 13.65 to 10.05 ms. With the barrier
stall gone, `ncu` showed the next wall: memory-instruction throttling
(`mio_throttle` 15.3, `lg_throttle` 10.2). `gemmPlain`, the GEMM for FP32
weights, computed one output per lane with a global weight load and a shared
input load **per multiply-add**. For ParT that is one small classifier layer
and never mattered. For a `--float` export it is nearly all the arithmetic.

It now uses the same register tile as the INT8 GEMMs, so each weight load
serves TB rows and each input load serves CB columns. Each accumulator starts
at its bias and walks k in order, which is the arithmetic the old version did.
On the CPU, where nothing is contracted into FMAs, the interpreter is
**byte-identical** to before on all four models at 32 and 256 lanes.

Tiling alone gives 13.65 -> 10.16 ms; warp layout alone 13.65 -> 10.05. They
compound, because each removes the bottleneck the other exposes.

### Keeping it off ParT's hot path

The first tiled version cost ParT ~3%, although ParT only runs its untiled
case. It was not instruction-cache capacity: one extra instantiation cost as
much as eleven. Adding **any** tiled instantiation pushed nvcc's inliner over
its budget for this very large kernel, and it responded by outlining ~57
*other* functions (`RET` count 3 -> 60; stack frame 136 -> 264 bytes). That put
real calls and register spills on paths ParT does run.

`gemmPlainTile` is therefore `__attribute__((noinline))`. With that the stack is
back to 136 bytes, only the twelve tiles are out of line, and ParT is back to
baseline: per-round ratio 0.997 (`b_kinadd`) and 0.998 (`softjet_part_unique`)
over five interleaved rounds, both ranges straddling 1.0. A call per Linear op
is nothing next to the GEMM it performs. `test_graph_kernel_source.py` fails if
the attribute is removed.

### Measured (L4, identical jets)

DeepSets:

| jets | block, before | **warp + tiled** | TorchScript | hand-fused prototype |
|---|---|---|---|---|
| 512 | 0.135 ms | **0.041 ms** | 0.150 ms | 0.040 ms |
| 8,192 | 2.00 ms | **0.336 ms** | 0.491 ms | 0.282 ms |
| 56,052 | 13.65 ms | **2.24 ms** | 6.67 ms | 1.83 ms |

At the real jets-per-event that is **6.1x faster than before and 3.0x faster
than TorchScript**, within 1.2x of a kernel hand-written for this one network.

The tiled GEMM also helps block layout: the tau model (a `--float` export,
every layer FP32) goes from **2782 to ~1100 ms** at 56,052 jets, 2.5x.

In the framework (scouting config, 56,052 jets/event, 40 events, one stream),
`SoftJetIdSC4` per event:

| | per event |
|---|---|
| DeepSets, block layout | 18.26 ms |
| DeepSets, warp layout (`auto`) | **5.41 ms** |
| ParT `b_kinadd`, `auto` (block) | 469 ms, unchanged (476 ms before) |

The framework gain (3.4x) is smaller than the kernel-only one (6.1x) because
the module time includes the input transform, which is the same for both
layouts. On 241 real jets warp and block layout agree **241/241 bit-for-bit**.

### Checking it

- `bench_graph --warp` runs the producer's warp kernel, including its
  pre-staged inputs. `--warp` fails if the model does not fit.
- The warp path is checked against the block path **at the same lane count**.
  Build with `-DNN_THREADS_PER_JET=32 -DNN_WARP_ARENA_FLOATS=10560
  -DNN_WARP_ARENA_INT8=4096`, and one warp per block gets the full arena. That
  runs ParT's whole op set (attention, LayerNorm, pair embedding, INT8 GEMMs)
  through `SubGroup`. With `-fmad=false` the two paths are **bit-identical on
  all four models**.
- With FMA contraction on, the models that read four-vectors differ at the ulp
  level (median 5.9e-8). The block and warp kernels are different template
  instantiations, and nvcc contracts the pair-feature arithmetic differently
  in each. `b_kinadd` amplifies that through INT8 rounding to ~1e-3, the same
  tail the compiled-vs-graph comparison has always had.
- Comparing at *different* lane counts also moves LayerNorm, whose partial
  sums are grouped by lane count.

## Open questions

1. **DP4A is confirmed in the release plugin**: 23,630 `IDP` instructions in
   `pluginPhysicsToolsParticleTransformerAlpakaPluginsCudaAsync.so`, which is
   built for sm_60/70/75/80/89/90, so there is no PTX JIT on an L4 either.
2. The producer builds, links and runs in CMSSW. Two defects had to be fixed
   first, both showing that neither `*_cfi.py` had ever been used: the module
   type was written `"alpaka/GraphInferenceProducer@alpaka"` rather than
   `"GraphInferenceProducer@alpaka"`, and the resident model was allocated
   through the *caching* `make_device_buffer(queue, ...)` while the registry
   holding it is deliberately leaked -- so the caching allocator's
   `liveBlocks_.empty()` assertion aborted the job at exit. It now uses the
   plain device overload.
3. The input transform (`TransformKernel.dev.cc`) is still hard-coded to the
   `b_kinadd` feature list. Generalising it means reading the Weaver data
   config, which is a separate job from this one.
4. `kMaxTokens` is 32 and `kMaxOutputs` is 16. Both are checked at load time
   with a clear message; raising them costs shared memory the arena currently
   uses.
5. **`GraphInferenceAlgo` publishes only `shared.output[0]`**, into a
   single-float output column, and stages 15 features plus a four-vector from
   the `SoftJetInputTensorSoA` layout. The interpreter computes all four of the
   tau model's outputs and the model file declares them, but running that model
   in CMSSW needs a producer variant with a 10-feature input and a four-wide
   output — a data-format change, not a graph one.
6. There is no real-input sample for the tau model in this repository, so its
   numbers above are on random features. Prior work here showed those are far
   outside the training distribution and inflate a difference by about an order
   of magnitude, which cuts both ways: the FP32 agreement is if anything better
   than quoted, and the INT8 gap smaller.
