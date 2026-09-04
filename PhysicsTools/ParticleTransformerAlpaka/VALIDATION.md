# Numerical validation

The included model was exported from `ParT_checkpoint.pt` with symmetric
per-output-channel INT8 quantization of the 48 dense tensors. Biases,
normalisation parameters, activations and accumulators remain FP32.

This document covers the *implementation* of the quantized graph. The
quantization error itself is a property of the export and is unchanged by the
rewrite; `interface/alpaka/GeneratedModel.quantization.json` records the
per-tensor reconstruction error.

## Provenance

The shipped model comes from the quantization-aware checkpoint
`net_best_epoch_state_q.pt` (sha256 `ff19bd5c...`), exported by
`scripts/export_qat_checkpoint.py`.  That checkpoint differs from the ones the
older exporter handled in three ways, all absorbed by the new script:

* every quantized layer carries an extra `.wrapped.` level in its parameter
  names plus `.observer.*` buffers;
* attention is stored as separate `q_proj` / `k_proj` / `v_proj`, which the
  generated graph wants concatenated into one `in_proj`;
* the activation ranges learned during training travel with the model, so a
  future kernel can quantize activations instead of dequantizing weights.

Two graph conventions were set from the model source rather than inherited:
`pair_final_gelu=false`, because `use_pre_activation_pair` strips the GELU
after the last pair-embedding BatchNorm, and `transpose_scaled_heads=false`,
because `c_attn` is applied as `(batch, token, head, head_dim)`.

## Printed 14-particle jet

The two missing particles in the fixed 16-particle input are zero padded and
masked.  The class-0 probability:

| implementation | class 0 |
|---|---:|
| PyTorch, FP32 | 0.227766678 |
| PyTorch, weights INT8 and FP32 activations | 0.227739170 |
| Alpaka scalar reference kernel, integer arithmetic | 0.228496164 |
| Alpaka block-cooperative kernel, dp4a path | 0.228496164 |
| PyTorch, full QAT (INT8 activations too) | 0.229284480 |

The third and fourth rows are what the CMSSW package computes, and they sit
`7.9e-04` from the QAT simulation the model was trained against.  That residual is entirely the first pair-embedding
convolution: quantization-aware training kept it in FP32, and the export
quantizes it because the graph's `pairLinear` accessor has a single return
type.  Restoring it costs a separate descriptor for that one layer; with the
pairwise bias switched off, the exported graph and PyTorch agree to `5.3e-08`.

The classifier had the same problem and was fixed, because its effect was an
order of magnitude larger: quantizing 256 weights sitting directly under the
softmax shifted the output by `7.3e-04`, more than every other layer combined.
It is now stored as FP32 in the parameter blob and described by `FloatLinear`.

For reference, the quantization the model was trained for moves this jet by
`1.5e-03` (row 1 against row 5), so both residuals are well inside it.

## Independent reference

`scripts/` also carries a NumPy reconstruction of the exported graph, built
only from the emitted header and binary blobs.  It reproduces the scalar C++
reference kernel to `1e-09` on this jet, and reproduces PyTorch layer by layer
to a relative `2e-07` through the embedding, all eight particle blocks and both
class blocks.

## Execution layout

The CMSSW path evaluates the whole event in a single launch of
`JetInferenceKernel`, one 256-thread Alpaka block per jet. All activations live
in `part::JetShared`, a 43,776-byte block-shared structure; there is no
device workspace buffer. The footprint fits under the 48 kB CUDA static
shared-memory limit and under the 47 kB default Alpaka CPU block-shared arena,
with about 4 kB of margin.

The dense weights are transposed on the host into the `[in][out]` layout that
the kernel reads (`ModelTranspose.h`). The transposition visits all 48 dense
tensors and covers exactly 2,108,032 INT8 values, verified at startup by the
producer, which throws if the count disagrees with `GeneratedModel.h`. Only the
ordering of the stored bytes changes: no value is requantized and the
per-output-channel scales are untouched, so the model on the device is
bit-for-bit the one the exporter produced.

## Integer arithmetic

The dense projections are evaluated with INT8 activations and INT32
accumulators, so the multiply-accumulate runs on `dp4a` (four per instruction)
instead of the FP32 pipe (one per lane per clock).

A layer takes the integer path when quantization-aware training gave it an
activation scale *and* its input width is a multiple of 16, which is what lets
the kernel fetch sixteen INT8 activations in one 16-byte load.  That covers
99.8% of the multiply-accumulates; the two exceptions are the 15-input token
embedding and the 4-input pair embedding, which keep FP32 arithmetic but still
round their activations onto the INT8 grid so the result matches what training
simulated.

`part::usesIntegerArithmetic` and `part::quantizeOne` live in ModelDefinition.h
and are shared by the host packer, the scalar reference and the block kernel, so
all three necessarily agree on which layers are integer and on how a value is
quantized.  `quantizeOne` rounds half-to-even, matching `torch.round`.

Weights for integer layers are packed `[k/4][out][4]`: four consecutive inputs
in one 32-bit word per output, so a warp's weight fetch is a single coalesced
transaction and each word feeds one `dp4a`.  Float layers keep the `[in][out]`
transposed layout.  The accumulator cannot overflow: with a 512-wide input the
largest partial sum is 512 x 127 x 127, about 8.3e6.

Both scales come out of the sum, since `s_x` is a property of the tensor and
`s_w` of the output channel, so the kernel computes
`bias + s_x * s_w[o] * sum_k q(x_k) q(w_ok)` with one multiply per output.  That
is algebraically what quantization-aware training simulated.

The three attention projections observed their inputs separately.  In a particle
block all three see the same normalised sequence and the scales coincide, but in
a class block the query sees the raw class token: for `cls_blocks.0` the query
scale is 0.000793 against 0.045836 for key and value, a factor of 58.  The
exporter therefore stores three scales per layer, and the class block quantizes
its query separately.

### Expected effect

| | instructions per multiply-accumulate |
|---|---:|
| FP32 with vector activation loads | 1.375 |
| INT8 with dp4a | 0.344 |

Four times fewer instructions on 99.8% of the arithmetic, which against the
6.1% of the kernel that is not GEMM gives about **3.4x**, or roughly 210 ms for
the 50,000-jet event that currently takes 700 ms on an L4.  This is an estimate
from the instruction mix, not a measurement; no GPU was available here.

### Reproducibility

Quantized inference is not bitwise stable across reduction orders.  The block
kernel's LayerNorm sums in a different order from the scalar reference, and a
difference of a few 1e-7 before rounding occasionally moves an activation across
an INT8 boundary.  Each flip shifts one input by a whole quantization step, so
agreement between the two implementations is a few 1e-3 rather than the 1e-7 of
the FP32 build.  Roughly one or two activations flip per jet out of about
190,000.  For a fixed block size the kernel is deterministic; the spread is
between block sizes, and it is well inside the 1.5e-3 that quantization itself
moves the output.

## Register tiling## Register tiling

Each dense layer picks its own `TB x CB` register tile at compile time
(`tokenBlockFor` / `columnBlockFor`). The task space of a GEMM is
`(token tiles) x (column groups)`, and the tile is the largest one that still
hands work to all 256 lanes, preferring a large `TB` first because each weight
element loaded is reused for `TB` tokens.

The first version distributed only over output columns, so a 128-wide
projection gave work to 32 of 256 lanes. Averaged over the graph by
multiply-accumulate count, lane utilisation was 45%; it is now 93%, and the
weight loads per jet fall from 7.25 M to 5.07 M. A large `TB` matters more here
than in the FP32 build, because it also amortises the integer-to-float
conversion of each weight over `TB` tokens.

At 5.07 M INT8 loads per jet the weight traffic is 253 GB per 50,000-jet event,
against roughly 1 TB for the same kernel over FP32 weights, and the 2.1 MB
model stays resident in the L2 of any current device.

| layer | MMAC/jet | tile | lanes used |
|---|---:|---|---:|
| block fc1 128->512 | 8.39 | 8x2 | 100% |
| block fc2 512->128 | 8.39 | 4x1 | 100% |
| block qkv 128->384 | 6.29 | 8x1 | 75% |
| block out 128->128 | 2.10 | 8x1 | 100% |
| pair 64->64 | 1.11 | 8x1 | 100% |
| class-block KV 128->256 | 1.05 | 8x1 | 100% |

The single-token projections of the class blocks cannot fill the block at any
tile; they are 0.5% of the arithmetic.

## Skipping masked tokens is exact

Masked keys are forced to `-3.4e38` before the softmax. `exp` of that minus any
finite maximum underflows to exactly zero, so a masked token cannot contribute
to any surviving token, and evaluating only the first `active` tokens is an
identity rather than an approximation. `active` is defined as the last unmasked
index plus one, which handles the holes that `ComputeInputsKernel` leaves when
it skips candidates with non-positive pt. Stale scratch rows are zeroed at
kernel start so that a zero softmax weight can never multiply an uninitialised
value into a NaN.

The pairwise features are exactly symmetric under `i <-> j` — all four raw
features and the mask — so only the 136 entries of the upper triangle are
evaluated and stored instead of all 256.

## Block-shape independence

`test/compile_kernel_int8.cc` runs the printed 14-particle jet and 34
randomised jets covering every multiplicity from 1 to 16, half of them with a
hole punched in the mask, against the scalar single-thread reference kernel.
Each case is repeated at block sizes of 1, 32, 96 and 256 threads, so the
result must be independent of how the cooperative loops are partitioned.

```text
shared bytes per jet: 47872
printed jet reference 0.228496164 0.771503806
printed jet fast      0.228496164 0.771503866
largest difference from the scalar reference: 0.00136914849
```

The residual `3e-7` comes from the FP32 accumulations being performed in a
different order and from the hoisted scale. It is far below the quantization
error already present in the model.

`test/test_launch_layout.py` verifies the triangular pair indexing round-trip
for every sequence length, the token padding rule, the scratch sizing against
each of its consumers, the 43,776-byte shared footprint, and the source-level
invariants of the design (single launch, no workspace, per-device model upload,
no barrier inside a divergent loop).

## Input transform

`ComputeInputsKernel` now uses one block of 32 threads per jet and locates the
bunch crossing by binary search into the jet lookup. Per-particle four-vectors
are computed in parallel, but the jet four-vector sum stays serial and in index
order, and candidates with non-positive pt contribute exact zeros, so the sum
and every feature derived from it are bit-identical to the previous serial
transform. The old code called `syncBlockThreads` inside a divergent
`independent_group_elements` loop, which is undefined behaviour on the CUDA
back-end whenever the jet count of a bunch crossing is not a multiple of the
block size.

## Corrected graph semantics

The generated implementation follows the model graph in these details:

1. GELU is applied after all three feature-embedding linear layers.
2. The legacy-MHA model applies GELU after the last pair convolution.
3. Learned head scaling is applied after the attention output projection.
4. The legacy-MHA model performs the `tbhd -> tbdh` channel permutation.
5. Class-attention Q is projected from the unnormalized class token; K and V
   are projected from the normalized concatenated class and particle tokens.

## Older checkpoint

`export_weaver_checkpoint.py` also detects the older custom-attention
checkpoint from its original `attn.in_proj.weight` key spelling. It selects no
final pair GELU and no head-channel permutation. `test/test_checkpoint_quantization.py`
compares the quantized and unquantized NumPy graphs for whichever checkpoints
are present.

These local checks validate the generated C++ graph, the block-cooperative
kernel against the scalar reference, and the exported model bytes. They do not
replace a complete `scram b` build for each enabled CMSSW Alpaka backend.
