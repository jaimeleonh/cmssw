# ParticleTransformerAlpakaFP32

Full-precision Alpaka translation of the supplied Weaver ParticleTransformer
for CMSSW. This package is independent from `ParticleTransformerAlpaka` and can
be installed alongside it. Its plugin is
`alpaka/ParticleTransformerFP32Producer@alpaka`.

## Precision and storage

No quantization is performed. Every learned tensor is preserved in the
checkpoint's original IEEE-754 float32 representation. Linear and convolution
weights are stored in `data/model_weights.fp32.bin`; biases, normalization
parameters, running statistics, class token, and learned attention/residual
scales are stored in `data/model_params.fp32.bin`. The producer checks both
files' exact byte sizes before uploading them as float buffers.

The generated model contains:

- 2,108,288 FP32 matrix-weight values (8,433,152 bytes);
- 34,340 FP32 remaining parameter values (137,360 bytes);
- FP32 activations and FP32 accumulation.

The binary blobs avoid compiling millions of numeric literals for each CMSSW
backend. `GeneratedModel.h` contains only offsets, tensor descriptors, sizes,
and SHA-256 provenance.

## Regeneration

The exporter accepts both supported Weaver QKV naming conventions:
`attn.in_proj_weight` / `attn.in_proj_bias`, and the older
`attn.in_proj.weight` / `attn.in_proj.bias`.

It uses that original spelling to select the matching Block semantics:

- `legacy_mha`: `torch.nn.MultiheadAttention`, post-projection head scaling
  with the legacy `tbhd -> tbdh` permutation, and a final pair GELU;
- `modern_custom`: Weaver custom attention, post-projection head scaling
  without that permutation, and no final pair GELU.

Both variants use a GELU after every feature-embedding linear layer and use the
unnormalized class token for the class-attention query. The concatenated class
and particle tokens are normalized for the keys and values.

```bash
python3 scripts/export_weaver_checkpoint_fp32.py ParT_checkpoint.pt \
  interface/alpaka/GeneratedModel.h \
  --weights data/model_weights.fp32.bin \
  --params data/model_params.fp32.bin \
  --report data/model_fp32.json
```

Normally `--model-variant auto --pair-final-activation auto` is sufficient.
The choices can be made explicit for reproducibility, for example:

```bash
python3 scripts/export_weaver_checkpoint_fp32.py ParT_checkpoint.pt \
  interface/alpaka/GeneratedModel.h \
  --weights data/model_weights.fp32.bin \
  --params data/model_params.fp32.bin \
  --report data/model_fp32.json \
  --model-variant legacy_mha \
  --pair-final-activation gelu
```

The included data were generated from checkpoint SHA-256
`f655cdf03c0c010f9ff3a7034a8b8c46492ce26d8adaa19e87f71bd27a8524cc`.

## CMSSW installation and build

Place the directory at:

```text
$CMSSW_BASE/src/PhysicsTools/ParticleTransformerAlpakaFP32
```

Then build with:

```bash
cd "$CMSSW_BASE/src"
scram b -j 1
```

Configure the five L1 scouting input tags in
`python/particleTransformerFP32Producer_cfi.py`. The producer creates the same
15-feature, four-vector, and mask inputs as the INT8 package, runs one FP32
inference per jet, and emits `l1sc::SoftJetOutputDeviceTensor` containing the
class-0 probability.

## How the inference kernel is organised

The whole event is evaluated by **one kernel launch**, with **one Alpaka block
of 256 threads per jet** and **no global scratch buffer at all**. Every
activation of the graph — token states, QKV projections, attention scores, the
pair bias, the feed-forward hidden state — lives in 43,776 bytes of block
shared memory (`partfp32::JetShared`), which fits both the 48 kB CUDA static
shared limit and the 47 kB default Alpaka CPU block-shared arena. During
inference the only global traffic is the streaming read of the model weights,
which every resident block shares through L1 and L2.

Five design points carry the speedup:

1. **One launch, one block per jet.** An event with 50,000 jets launches 50,000
   blocks. The previous implementation chopped the event into `batchSize`
   batches of 32 jets and issued about 1,563 dependent launches of 32 blocks
   each, so the device was under 3% occupied and the launches were serialised
   by the queue.
2. **Shared memory instead of a global workspace.** The old kernel kept
   22,336 floats per jet in device memory (up to 349 MiB per stream) and every
   layer round-tripped through it. Shared memory removes that buffer, its
   allocation logic, and the bandwidth it consumed.
3. **Transposed weights.** The checkpoint stores dense weights row-major as
   `[out][in]`, so a thread owning one output column strides through memory.
   `ModelTranspose.h` rewrites them to `[in][out]` on the host at startup, so
   consecutive threads read consecutive columns and each weight fetch is a
   fully coalesced transaction. The transposition is a pure relabelling: no
   value is modified or rounded.
4. **Register-blocked GEMMs.** Every dense layer picks its own `TB x CB`
   accumulator tile at compile time. The task space is
   `(token tiles) x (column groups)`, and the tile is the largest one that
   still gives work to all 256 lanes, preferring a large `TB` because each
   weight element loaded is reused for `TB` tokens. Averaged over the graph,
   93% of lanes are busy inside a GEMM.
5. **Token-contiguous activations.** The activations feeding the heavy GEMMs
   are stored `[channel][token]`, so a thread's `TB` tokens for one `k` are
   adjacent and load in `TB / 4` vector instructions instead of `TB` scalar
   ones. Instruction issue, not shared-memory bandwidth, is what limits the
   inner loop: an SM reaches peak FP32 only when all four instructions it
   retires per clock are multiply-adds, and the scalar form spent half its
   slots on activation loads. See VALIDATION.md for the measured mix.
6. **Only the real work.** Tokens beyond the last valid particle are skipped,
   and only the `i <= j` half of the pairwise interaction matrix is evaluated
   (136 pairs instead of 256), because the four pairwise features and the pair
   mask are all exactly symmetric under `i <-> j`. Skipping masked tokens is an
   identity rather than an approximation: masked keys are forced to `-3.4e38`
   before the softmax, so their weight underflows to exactly zero and they
   cannot influence any surviving token.

The model is transposed and uploaded **once per device** and shared by every
CMSSW stream running on it, instead of once per stream. With several streams
per GPU this keeps a single 8.4 MB copy of the weights in the L2 working set
rather than one copy per stream.

The `batchSize` parameter is gone. Its replacement, `maxBlocks`, caps the grid
and makes each block loop over several jets; it defaults to `0`, meaning one
block per jet, which is the fastest configuration. It exists for occupancy
studies only — there is no memory reason to reduce it.

The input transform was rewritten as well. It now runs one block of 32 threads
per jet with a binary search into the jet-to-bunch-crossing lookup, computes
the per-particle four-vectors in parallel, and keeps only the jet four-vector
sum serial and in index order so the results stay bit-identical. The previous
version looped over bunch crossings and called `syncBlockThreads` inside a
divergent `independent_group_elements` loop, which is undefined behaviour on
the GPU back-ends.

The same code runs on the serial CPU backend, where Alpaka maps the block to
its supported execution model and the cooperative loops cover all elements
sequentially. The kernel is written for any block size: the regression test
exercises 1, 32, 96 and 256 threads per block.

### Occupancy note

43,776 bytes of shared memory per block is deliberate: it buys the removal of
all workspace traffic, at the cost of how many blocks an SM can hold at once.
On a device with 164 kB of shared memory per SM three blocks (768 threads) are
resident per SM, which is enough to hide the weight-stream latency because all
resident blocks read the same weights. If a future device or a CUDA
`cudaFuncAttributeMaxDynamicSharedMemorySize` change makes a different trade
attractive, the knobs are `kFfnTokenTile` and `kPairTile` in
`interface/alpaka/ModelConstants.h`: lowering either shrinks `JetShared` at the
cost of more passes over the weights.

## Local validation

`test/test_fp32_export.py` checks byte-for-byte equality between checkpoint
tensors and both binary model files and evaluates the independent NumPy graph.
`test/compile_kernel_fp32.cc` runs the printed 14-particle jet plus 34
randomised jets — every multiplicity from 1 to 16, half of them with a hole
punched in the mask — through the block-cooperative kernel against a minimal
Alpaka API stub, and compares each of them with the scalar reference kernel at
block sizes of 1, 32, 96 and 256 threads. `test/test_launch_layout.py` checks
the triangular pair indexing, the shared-memory budget, and the source-level
invariants of the new design. For the included model:

```text
NumPy: 0.2534404695 0.7465595603
C++:   0.25344044   0.74655962
```

### Interpreting FastTimerService cleanup time

GPU launches are asynchronous. If this producer is the last GPU module in a
path, CMSSW may perform the first queue synchronization during an event cleanup
transition. FastTimerService can therefore attribute inference execution to
`cleanup`, even though nothing is being freed there. Set
`synchronizeForTiming=True` temporarily to wait immediately after inference;
the time should move from cleanup to this producer. This option is diagnostic
and should remain false for production because it removes CPU/GPU overlap.

Across all 35 jets and all four block shapes the largest deviation from the
scalar reference is `2.38e-07`, i.e. the last bits of an FP32 accumulation
performed in a different order.

The older custom-attention checkpoint is also accepted automatically and gives
NumPy/C++ agreement within `2e-7` on the same input. A full SCRAM build is still
required in a CMSSW environment.
