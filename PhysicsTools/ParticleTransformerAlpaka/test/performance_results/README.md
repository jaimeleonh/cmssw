# Measured runtime comparison — 15 October 2026

The subsequent [complete CMSSW producer measurements](cmssw/README.md) include actual preprocessing, PyTorchAlpaka adapters, event cleanup, and stream tuning. The results below retain the earlier standalone inference scope.

These are new measurements on **one NVIDIA L4** in `cms-p2scout-ngt-01`, using the CMSSW `CMSSW_16_1_0_pre3` runtime, **Torch 2.6.0a0, CUDA 12.9, and driver 590.48.01**. The other L4 was idle. Each tested event below has **56,052 jets, all with 16 valid particles**. GPU inputs are shared by the two implementations and use the CMSSW-style component-major SoA layout, including 128-byte column alignment. Results for 512 and 8,192 jets and initial contiguous-input experiments are also retained in the JSON files.

The measured Torch side is the **C++ TorchScript inference engine used underneath PyTorchAlpaka**, not the complete `AlpakaModel`/`TensorCollection` adapter or a `cmsRun` producer. Both sides receive pre-existing GPU buffers; input transformation, transfer, loading/export, Python dispatch of individual model operations, and CMSSW event handling are excluded. Native input staging and Torch operations' internal layout conversions are included. The custom side runs the package's actual graph-interpreter headers through a CUDA Alpaka shim.

## Main results: tuned batch sizes, aligned SoA input

Times are median milliseconds per complete event, with five warmup evaluations and five measurement samples for ParT or seven for DeepSets. Each DeepSets large-event sample averages three complete evaluations. Best batch means **best of those tested**, not a proof that every possible kernel/precision configuration was exhausted.

| Matched checkpoint / precision | Graph ms/event | TorchScript ms/event | Best tested Torch batch | Winner | Graph / Torch jets/s |
|---|---:|---:|---:|---|---:|
| DeepSets, FP32 | 2.517 | 4.059 | 4,096 | Graph 1.61× | 22,272,459 / 13,808,850 |
| ParT, FP32 | 3752.792 | 1257.795 | 256 | Torch 2.98× | 14,936 / 44,564 |
| ParT, graph INT8 weights / Torch FP32 | 1470.411 | 1250.980 | 256 | Torch 1.18× | 38,120 / 44,806 |

The graph uses warp layout for DeepSets and block layout for ParT. DeepSets block mode takes 11.641 ms, so its result is strongly dependent on selecting the correct layout.

With **TF32 allowed**, the best tested TorchScript ParT result is **966.302 ms/event**, batch 256, or **58,007 jets/s**. This is a separately labeled precision policy; enabling the flag does not prove every operation used Tensor Cores. DeepSets with TF32 allowed reaches **3.019 ms**, batch 4,096, versus **2.417 ms** for graph warp in the same run: graph remains about **1.25× faster**. The TF32-versus-strict-Torch maximum probability difference at that best batch is **0.00063473**. FP16, TensorRT, AOT compilation, CUDA Graph capture, multiple concurrent inference streams, and other GPUs were not benchmarked.

**Precision clarification from the subsequent CMSSW tests:** the frozen ParT archive embeds `allow_tf32=True` in its convolution operations. Here, a disabled TF32 setting describes the global Torch flags, not proof that every serialized operation forbids TF32. The separate CMSSW runs using `NVIDIA_TF32_OVERRIDE=0` test the stricter policy. Historical raw JSON flag values retain their original meaning.

## Memory measurements

The additional Torch live-tensor peak at the fastest baseline-FP32 configuration is **17.0 MiB for DeepSets** (batch 4,096) and approximately **35.0 MiB for ParT** (batch 256). The ParT repeat paired with INT8 graph inference observes 34.1 MiB; allocator state and workspace reuse account for small peak differences between processes. These are increments above the warmed baseline, not complete device VRAM usage. With tuned batches, the observed activation saving is tens of MiB; the multi-GiB contrast belongs to the slower full-event Torch batch.

Making Torch process all 56,052 ParT jets in **one batch** takes **1920.034 ms** and raises its additional live-tensor peak to **7461.4 MiB** (about 7.29 GiB). It is slower than the tuned batch-256 configuration. Its reserved allocator peak is **11622.0 MiB**; reserved memory is not identical to live activations. The full-event DeepSets batch similarly uses 221 MiB of extra live tensors and takes 6.854 ms.

The graph has **no per-event global activation arena allocation**: the uploaded model/program is 11,744 bytes for DeepSets, 8,583,264 bytes for FP32 ParT, or 2,309,120 bytes for weight-only INT8 ParT. It shows zero forward increment in the Torch allocator and unchanged device-used memory between the pre/post samples of warmed graph evaluations. This does **not** mean zero total VRAM or zero global-memory accesses: model data and inputs/outputs are global, and compiler local/stack storage and runtime memory are outside Torch's allocator.

CUDA reports **47,840 shared bytes per block** for block mode, **43,776 per eight-jet warp block**, **128 registers and 136 local bytes per thread**, and a maximum of **two resident blocks per SM** for these compiled bridge kernels. The L4 has 58 SMs. Queuing 56,052 jets therefore provides work for many waves; it does not reserve 56,052 copies of shared state in VRAM. Common event I/O is excluded from the incremental activation comparison, but remains allocated in both implementations.

Torch model loading initially allocates 12,800 bytes for this DeepSets and 8,581,632 bytes for this ParT in its allocator. Those figures are close to the corresponding graph FP32 model storage; the main FP32 memory difference is activation materialization. INT8 further reduces graph model storage. Raw JSON records include baseline/peak allocated/reserved values and, for the later experiments, device-used memory snapshots. Do not compare absolute process memory from these same-process experiments as if only one engine were loaded.

## Numerical checks and input sample

Both engines load the **same source checkpoint** for each matched comparison. The graph is exported from that checkpoint in the requested precision. Every Torch batch configuration is compared with the graph's output zero over the entire event; all reported outputs are finite.

| Configuration | Maximum absolute probability difference | Mean absolute probability difference |
|---|---:|---:|
| DeepSets FP32 | 8.04663e-07 | 4.60389e-08 |
| ParT FP32 | 0.00183278 | 7.75772e-05 |
| ParT weight-only INT8 | 0.00578082 | 0.000543672 |

For TF32 allowed versus baseline Torch FP32 at its best batch, the maximum difference is **0.000183046**, mean **1.8057e-05**. Floating arithmetic and pair-feature cancellation mean that the FP32 graph/Torch results are not bitwise identical; the ParT maximum is about 0.18 percentage points. The INT8 maximum is about 0.58 percentage points. These are measured conversion/rounding differences, not a validation of tagging accuracy on a held-out physics sample.

Inputs use seed `20261015`, sampling particle feature/four-vector rows from the package's 14-particle printed fixture with replacement to construct 16-particle records. The original preprocessed feature values are retained, rather than recomputing relative features for a new reconstructed jet. This is a reproducible near-physical numerical/timing sample, not a new collision-data dataset. Full 16-particle occupancy removes any timing advantage from dropping a masked tail. `--multiplicity mixed` is supported by the harness but is not the primary dataset reported here.

A separate 128-jet `-fmad=false` diagnostic remains in `part_no_fma_check.json`. It still shows a maximum difference of `4.09e-4` on 128 jets, so disabling contraction alone does not make the two implementations identical. Its timing is not substituted into the production-flag comparison.

## Compiled QAT context: a different checkpoint

The package's existing benchmark programs were also compiled and executed on the L4, with 56,052 jets, 16 valid particles each, seven timed launches, and seed `20260921`:

| Shipped QAT model implementation | Best ms/event | Mean ms/event |
|---|---:|---:|
| Compiled `runParticleTransformer` | 577.662 | 581.964 |
| Graph interpreter, `b_kinadd.nnm` | 673.081 | 675.661 |

These two rows use the same shipped QAT weights as each other, with INT8 activations and DP4A. They use the existing benchmarks' seeded random inputs, not the resampled fixture dataset above. **They are not a matched comparison with the frozen TorchScript ParT checkpoint**. The original QAT TorchScript model is not supplied here, so dividing these times by the frozen-model Torch timings would confound checkpoint and arithmetic choices. Raw output is retained in `compiled_qat_context.txt` and `graph_qat_context.txt`.

## Reproduction and files

From the package directory in a configured CMSSW environment with GPU device access:

```bash
eval "$(scram runtime -sh)"
bash test/build_runtime_benchmark.sh /tmp
python3 test/benchmark_runtime_comparison.py --model deepsets \
  --soa --batches 512,2048,4096,8192 --repeats 7 \
  --output test/performance_results/deepsets_fp32_soa.json
python3 test/benchmark_runtime_comparison.py --model part --soa \
  --sizes 8192,56052 --batches 128,256,512,4096 --omit-full-batch \
  --repeats 5 --output test/performance_results/part_fp32_soa.json
python3 test/benchmark_runtime_comparison.py --model part --soa \
  --sizes 56052 --batches 256 --omit-full-batch --quantize \
  --repeats 5 --output test/performance_results/part_int8_soa.json
python3 test/benchmark_runtime_comparison.py --model part --soa \
  --sizes 56052 --batches 256,512,1024 --omit-full-batch --tf32 \
  --skip-native-timing --repeats 5 \
  --output test/performance_results/part_tf32_soa.json
```

The build targets `sm_89` (L4) and uses `-O3`, C++20 for CUDA, default FMA contraction, no relocatable device code, and the installed Torch C++ ABI. Source/runtime fingerprints are in [environment.json](environment.json); NVIDIA driver/toolkit output is in [gpu_environment.txt](gpu_environment.txt). GPU clocks were not locked and background interference was not exhaustively controlled; raw sample ranges are retained. Small percent-level timing differences should not be treated as decisive.

[summary.json](summary.json) collects the main rows. The `*_soa.json` files contain the primary data, `*_full.json`/`*_tuning.json` retain earlier contiguous-input sweeps, and `*_context.txt` records the separate shipped-QAT timings. Exported `.nnm` files and compiled shared libraries are generated in `/tmp`, with hashes in the reports; rerun export/build to recreate them.

- [benchmark_runtime_comparison.py](../benchmark_runtime_comparison.py): `digest` hashes source artifacts; `libraries` declares the ctypes interfaces; `export` builds/allocates/serializes a graph from the selected TorchScript file; `timed` warms and measures CUDA-event samples plus allocator peaks; `main` selects models/options, constructs common inputs, runs graph and Torch batch sweeps, checks outputs, and writes JSON incrementally. Nested `checkpoint` writes results, `as_soa` creates aligned component-major views, and `run_native`/`run_torch` call the native bridges.
- [benchmark_runtime_bridge.cu](../benchmark_runtime_bridge.cu): empty accelerator tag `Acc`; `GraphHandle` owns uploaded op/buffer/weight/parameter arrays and a view; `check` turns CUDA errors into exceptions; `upload` copies typed arrays; `stage` reads contiguous or aligned SoA input into a jet's arena; `graphBlock` and `graphWarp` call the production interpreter at block/warp scope and publish output zero. C entry points `graph_create`/`graph_destroy` manage the model, `graph_error` exposes errors, `graph_bytes` and `graph_fits_warp` expose metadata, `graph_launch` dispatches, and `graph_resources` reports CUDA function attributes/occupancy limits.
- [benchmark_torch_bridge.cc](../benchmark_torch_bridge.cc): `TorchHandle` owns a loaded TorchScript module and nonowning input/output tensor views; its constructor moves in the module. `ts_create` loads on the selected CUDA device and enables evaluation behavior; `ts_bind` wraps contiguous or aligned SoA device pointers without copying; `ts_run` loops over the selected batches under C++ inference mode and copies each probability to the common output; `ts_destroy` releases state and `ts_error` reports failures. The current-stream/default-stream contract must be preserved if this harness is extended.
- [build_runtime_benchmark.sh](../build_runtime_benchmark.sh): builds both bridges using the installed Torch include/library paths and ABI.

The measured inference results favor graph warp execution for this DeepSets, C++ TorchScript for this FP32/weight-only ParT, and the graph design for avoiding batch-proportional activation storage. The separate CMSSW report establishes the measured producer ordering on the baseline raw workload. A same-checkpoint QAT PyTorch baseline remains unavailable.
