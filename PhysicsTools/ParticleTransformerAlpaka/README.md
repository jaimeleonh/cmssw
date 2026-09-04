# ParticleTransformerAlpaka — b_kinadd

Weight-only INT8 translation of the supplied Weaver ParticleTransformer checkpoint, plus Alpaka input transformation, inference, and a CMSSW producer.

## Updating an existing CMSSW checkout

Replace the whole `PhysicsTools/ParticleTransformerAlpaka` directory when
installing a new package revision. Do not overlay individual files: the
`ParticleTransformerAlgo` declaration and definition form one versioned API.
After replacing the directory, remove the package's SCRAM intermediates or run
`scram b clean` before rebuilding.

## Resolved model and data contract

The exporter reads the PyTorch zip `state_dict` without importing PyTorch. The supplied checkpoint resolves to:

* input: 15 PF features, 4-vector and mask, fixed length 16;
* embedding: 15 → 128 → 512 → 128;
* pair encoder: four `pp` Lorentz features, 4 → 64 → 64 → 64 → 8;
* 8 self-attention blocks, 8 heads, hidden size 128, FFN size 512;
* 2 class-attention blocks and a 128 → 2 classifier;
* LayerNorm, BatchNorm, GELU, post-projection normalization, head scaling and residual scaling.

Both supported Weaver serialization variants are accepted automatically.
Their original QKV key spelling selects the exact Block graph:

- PyTorch `MultiheadAttention` (`attn.in_proj_weight`) uses post-projection
  head scaling with the legacy `tbhd -> tbdh` permutation and a final pair
  GELU;
- custom attention (`attn.in_proj.weight`) uses post-projection head scaling
  without that permutation and no final pair GELU.

Both paths apply GELU after all feature-embedding layers. Class-attention Q is
projected from the unnormalized class token, while K and V are projected from
the normalized concatenated class and particle tokens. The automatic choices
can be recorded explicitly with `--model-variant` and
`--pair-final-activation`.

The YAML preprocessing order is: transformed `log(pt)`, `log(E)`, `log(pt/jet_pt)`, `log(E/jet_energy)`, and `deltaR`; followed by `deta`, `dphi`, `z0`, `dxy`, `pweight`; then five particle-ID indicators. See `data/b_kinadd.yaml` for exact constants and ordering.

## Quantization

All 48 matrix/1×1-convolution tensors use symmetric INT8 quantization independently per output channel. Biases, normalization parameters, activations and accumulators remain FP32. `GeneratedModel.h` contains only model sizes, hashes, offsets, and device-safe accessors; the model values are stored in `data/model_weights.int8.bin` and `data/model_params.fp32.bin`. The supplied checkpoint SHA-256 is:

`f655cdf03c0c010f9ff3a7034a8b8c46492ce26d8adaa19e87f71bd27a8524cc`

Regenerate it with:

```bash
python3 scripts/export_weaver_checkpoint.py ParT_checkpoint.pt \
  interface/alpaka/GeneratedModel.h \
  --weights data/model_weights.int8.bin \
  --params data/model_params.fp32.bin \
  --report interface/alpaka/GeneratedModel.quantization.json
```


## Model provenance

The shipped model is exported from the quantization-aware checkpoint
`net_best_epoch_state_q.pt` by `scripts/export_qat_checkpoint.py`.  Dense
weights are symmetric INT8 per output channel; the classifier stays FP32
because training kept it that way and quantizing 256 weights under the softmax
moves the output more than every other layer combined.

The checkpoint also carries the activation ranges learned during training.  The
exporter writes them into the parameter blob -- every `Linear` descriptor has an
`activation` pointer holding three scales, for the key, query and value inputs
-- and **the kernel uses them**: activations are quantized to INT8 and the dense
projections run on `dp4a` with INT32 accumulators.  See VALIDATION.md.

## How the inference kernel is organised

The whole event is evaluated by **one kernel launch**, with **one Alpaka block
of 256 threads per jet** and **no global scratch buffer at all**. Every
activation of the graph — token states, QKV projections, attention scores, the
pair bias, the feed-forward hidden state — lives in 43,776 bytes of block
shared memory (`part::JetShared`), which fits both the 48 kB CUDA static
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
   layer round-tripped through it. That FP32 scratch dwarfed the INT8 model it
   was there to serve. Shared memory removes the buffer, its allocation logic,
   and the bandwidth it consumed.
3. **Transposed weights.** The exporter stores dense weights row-major as
   `[out][in]`, so a thread owning one output column strides through memory.
   `ModelTranspose.h` rewrites them to `[in][out]` on the host at startup, so
   consecutive threads read consecutive columns and a warp's weight fetch is a
   single coalesced transaction. The transposition is a pure relabelling: no
   value is modified, rounded or requantized, and the per-channel scales are
   untouched.
4. **Integer arithmetic.** 99.8% of the multiply-accumulates run as `dp4a`,
   four per instruction, against one per lane per clock on the FP32 pipe.
   Weights for those layers are packed `[k/4][out][4]` so a warp's fetch is one
   coalesced transaction and each 32-bit word feeds one `dp4a`.  The two
   projections whose input width is not a multiple of 16 keep FP32 arithmetic
   but still round their activations onto the INT8 grid, so the whole graph
   reproduces what quantization-aware training simulated.

5. **Register-blocked GEMMs, with the scale hoisted.** Every dense layer picks
   its own `TB x CB` accumulator tile at compile time. The task space is
   `(token tiles) x (column groups)`, and the tile is the largest one that
   still gives work to all 256 lanes, preferring a large `TB` because each
   weight element is loaded, converted from INT8 to float, and then reused for
   `TB` tokens. Averaged over the graph, 93% of lanes are busy inside a GEMM.
   The per-output-channel scale is a property of the column, so it comes out of
   the `k` sum entirely: the inner loop accumulates against the raw integers
   and the scale is applied once per output rather than once per multiply-add
   as the scalar reference does. That removes one FP32 multiply from every
   multiply-accumulate in the network.
6. **Token-contiguous activations.** The activations feeding the heavy GEMMs
   are stored `[channel][token]`, so a thread's `TB` tokens for one `k` are
   adjacent and load in `TB / 4` vector instructions instead of `TB` scalar
   ones. Instruction issue, not shared-memory bandwidth, is what limits the
   inner loop: an SM reaches peak FP32 only when all four instructions it
   retires per clock are multiply-adds, and the scalar form spent half its
   slots on activation loads. See VALIDATION.md for the measured mix.
7. **Only the real work.** Tokens beyond the last valid particle are skipped,
   and only the `i <= j` half of the pairwise interaction matrix is evaluated
   (136 pairs instead of 256), because the four pairwise features and the pair
   mask are all exactly symmetric under `i <-> j`. Skipping masked tokens is an
   identity rather than an approximation: masked keys are forced to `-3.4e38`
   before the softmax, so their weight underflows to exactly zero and they
   cannot influence any surviving token.

The model is transposed and uploaded **once per device** and shared by every
CMSSW stream running on it, instead of once per stream. With several streams
per GPU this keeps a single 2.1 MB copy of the weights in the L2 working set
rather than one copy per stream.

Being INT8 helps most here. Weight reads are 253 GB per 50,000-jet event,
against roughly 1 TB for the same design over FP32 weights, and the whole model
is small enough to stay resident in the L2 of any current device. The cost is
one integer-to-float conversion per weight load, which is why the tile chooser
prefers a large `TB`: at `TB = 8` a single conversion is amortised over eight
tokens.

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

## CMSSW installation and build

Place the directory at:

```text
$CMSSW_BASE/src/PhysicsTools/ParticleTransformerAlpaka
```

Then build with:

```bash
cd "$CMSSW_BASE/src"
scram b -j 1
```

Configure the five L1 scouting input tags in
`python/particleTransformerProducer_cfi.py`. The producer builds the
15-feature, four-vector and mask inputs, runs one inference per jet, and emits
`l1sc::SoftJetOutputDeviceTensor` containing the class-0 probability.

## Local validation

`test/test_quantization.py` checks the per-channel INT8 error bound and the
byte sizes and hashes of both model files. `test/compile_kernel_int8.cc` runs
the printed 14-particle jet plus 34 randomised jets — every multiplicity from 1
to 16, half of them with a hole punched in the mask — through the
block-cooperative kernel against a minimal Alpaka API stub, and compares each
of them with the scalar reference kernel at block sizes of 1, 32, 96 and 256
threads. `test/test_launch_layout.py` checks the triangular pair indexing, the
shared-memory budget, and the source-level invariants of the design.

```text
shared bytes per jet: 43776
printed jet reference 0.252941519 0.747058511
printed jet fast      0.252941519 0.747058511
largest difference from the scalar reference: 5.06639481e-07
```

Across all 35 jets and all four block shapes the largest deviation from the
INT8 scalar reference is `5.07e-07`, i.e. the last bits of an FP32 accumulation
performed in a different order. The quantization error itself is unchanged by
this rewrite and is much larger: the same jet gives `0.25344044` through the
FP32 build of the same graph.

### Interpreting FastTimerService cleanup time

GPU launches are asynchronous. If this producer is the last GPU module in a
path, CMSSW may perform the first queue synchronization during an event cleanup
transition. FastTimerService can therefore attribute inference execution to
`cleanup`, even though nothing is being freed there. Set
`synchronizeForTiming=True` temporarily to wait immediately after inference;
the time should move from cleanup to this producer. This option is diagnostic
and should remain false for production because it removes CPU/GPU overlap.

A full SCRAM build is still required in a CMSSW environment.
