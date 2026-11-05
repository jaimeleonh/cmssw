# Numerical validation

The included model was exported from `ParT_checkpoint.pt` without quantization.
All learned tensors are stored as their original IEEE-754 float32 bytes.

## Printed 14-particle jet

The two missing particles in the fixed 16-particle input are zero padded and
masked. The class probabilities are:

| implementation | class 0 | class 1 |
|---|---:|---:|
| frozen TorchScript graph, independent NumPy reconstruction | 0.2534404695 | 0.7465595603 |
| generated Alpaka C++ kernel (serial API stub) | 0.25344044 | 0.74655962 |

The absolute class-0 difference is approximately `3e-8`.

## Execution layout

The CMSSW path evaluates the whole event in a single launch of
`JetInferenceKernel`, one 256-thread Alpaka block per jet. All activations live
in `partfp32::JetShared`, a 43,776-byte block-shared structure; there is no
device workspace buffer. The footprint fits under the 48 kB CUDA static
shared-memory limit and under the 47 kB default Alpaka CPU block-shared arena,
with about 4 kB of margin.

The dense weights are transposed on the host into the `[in][out]` layout that
the kernel reads (`ModelTranspose.h`). The transposition visits all 48 dense
tensors and covers exactly 2,108,288 floats, verified at startup by the
producer, which throws if the count disagrees with `GeneratedModel.h`. Because
only the ordering of the stored floats changes, the FP32 values themselves are
bit-identical to the checkpoint.

## Activation layout and instruction issue

An SM retires four instructions per clock and reaches peak FP32 only when all
four are fused multiply-adds. In the first version of the register-tiled GEMM
each `k` iteration issued `TB` scalar shared-memory loads for the activations,
one per token, alongside `CB` weight loads and `TB * CB` multiply-adds. At
`TB = 8, CB = 1` that is 8 LDS + 1 LDG per 8 FFMA: an FMA share of 47%, which
caps the kernel at 47% of the FP32 pipe before occupancy, barriers or memory
are considered at all. The shared-memory *bandwidth* was never the problem --
those reads are broadcasts -- but each one still costs an issue slot.

The activations that feed the heavy GEMMs are therefore stored
token-contiguous: element `(t, k)` lives at `[k * stride + t]` rather than
`[t * stride + k]`. A thread's `TB` tokens for one `k` are then adjacent, and
`TB / 4` sixteen-byte loads replace `TB` scalar ones. `inputStride` and the
tile's first token are both multiples of four and the buffers are declared
`alignas(16)`, so the vector loads are always aligned.

LayerNorm produces these buffers, so it gained the same layout flags and writes
transposed directly; no separate transpose pass exists. The layers converted
are the input embedding, both embedding feed-forward matrices, the block QKV
projection, both block feed-forward matrices and the class-block key/value
projection -- 87% of all multiply-adds. The attention output projection, the
pair encoder and the single-token class projections keep the row-major path,
which the same template still supports.

| | FMA share of issued instructions | ceiling |
|---|---:|---:|
| scalar activation loads | 50% | 50% of FP32 peak |
| vector activation loads | 68% | 68% of FP32 peak |

On an L4 (58 SM, 128 FP32/clk/SM, 2.04 GHz) the arithmetic floor for a
50,000-jet event is 95 ms, so this moves the issue-limited ceiling from about
200 ms to about 143 ms. The result is bit-for-bit unchanged: the regression
test reports the same worst-case deviation as before the change.

## Register tiling

Each dense layer picks its own `TB x CB` register tile at compile time
(`tokenBlockFor` / `columnBlockFor`). The task space of a GEMM is
`(token tiles) x (column groups)`, and the tile is the largest one that still
hands work to all 256 lanes, preferring a large `TB` first because each weight
element loaded is reused for `TB` tokens.

The first version distributed only over output columns, so a 128-wide
projection gave work to 32 of 256 lanes. Averaged over the graph by
multiply-accumulate count, lane utilisation was 45%; it is now 93%, and the
weight loads per jet fall from 7.25 M to 5.06 M.

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

`test/compile_kernel_fp32.cc` runs the printed 14-particle jet and 34
randomised jets covering every multiplicity from 1 to 16, half of them with a
hole punched in the mask, against the scalar single-thread reference kernel.
Each case is repeated at block sizes of 1, 32, 96 and 256 threads, so the
result must be independent of how the cooperative loops are partitioned.

```text
shared bytes per jet: 43776
printed jet reference 0.25344044 0.74655962
printed jet fast      0.25344047 0.74655956
largest difference from the scalar reference: 2.38418579e-07
```

The residual `3e-7` is the last bit of an FP32 accumulation performed in a
different order; it is far below the FP32 epsilon accumulated over the roughly
32 million multiply-adds each jet requires.

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

`export_weaver_checkpoint_fp32.py` also detects the older custom-attention
checkpoint from its original `attn.in_proj.weight` key spelling. It selects no
final pair GELU and no head-channel permutation. On the same printed jet, the
independent NumPy implementation gives `0.22407684` and the generated C++
kernel gives `0.22407664` for class 0.

These local checks validate the generated C++ graph, the block-cooperative
kernel against the scalar reference, and the checkpoint bytes. They do not replace a complete `scram b` build for each
enabled CMSSW Alpaka backend.
