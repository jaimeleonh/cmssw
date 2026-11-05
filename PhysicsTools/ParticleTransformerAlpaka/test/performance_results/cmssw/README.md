# Measurements inside CMSSW — 15 October 2026

These measurements run **`cmsRun` with the actual CMSSW Alpaka CUDA backend and the actual PyTorchAlpaka `AlpakaModel` / `TensorCollection` adapter**. They use the raw input, unpacking, vertices, seeded-cone clustering, and event selection from the user's [run_benchmark.sh](../../../../../L1TriggerScouting/Phase2/test/run_benchmark.sh) and its [configuration](../../../../../L1TriggerScouting/Phase2/test/runScoutingPhase2PF_cfg.py). The baseline files were not edited. The active baseline model was replaced by the same DeepSets and frozen ParT checkpoints used in the earlier standalone comparison, so the graph/Torch pairs have matching weights and output meaning.

The hardware is **one NVIDIA L4 (GPU 0) on `cms-p2scout-ngt-01`**, with `CMSSW_16_1_0_pre3`, Torch 2.6.0a0, CUDA 12.9, and driver 590.48.01. GPU 1 was unused. Primary asynchronous runs execute without another GPU benchmark or build running. GPU clocks were not locked.

## Completed-event throughput and GPU footprint

The following are **one-thread, one-framework-stream, normal asynchronous** results. Each measured orbit includes source reading, PF/vertex unpacking, event selection, clustering, input preparation, inference, output publication, and completion of framework cleanup. Loading/startup is excluded. There is no NanoAOD writing. They are complete-pipeline results, not isolated inference latency. Each inspected orbit has **56,052 jets**.

GPU memory is the **sampled device-wide peak above the idle baseline**, in MiB, from isolated processes; it includes common pipeline buffers, context/code/library allocations, caches, and model/intermediate storage. It is not an incremental Torch activation counter. See the memory definitions below.

| Configuration, one framework stream | Pipeline ms/orbit | Jets/s | GPU MiB above idle |
|---|---:|---:|---:|
| No tagging (unpacking + clustering) | 2.250 | 24,913,655 | 280.0 |
| DeepSets graph FP32, warp layout | 7.182 | 7,804,111 | 408.0 |
| DeepSets PyTorchAlpaka FP32, batch 4096 | 6.657 | 8,420,608 | 502.0 |
| ParT graph FP32 | 1841.702 | 30,435 | 418.0 |
| ParT graph weight-only INT8 | 694.198 | 80,743 | 412.0 |
| ParT PyTorchAlpaka FP32/global TF32 off, batch 256 | 1243.892 | 45,062 | 546.0 |
| ParT PyTorchAlpaka, NVIDIA_TF32_OVERRIDE=0, batch 256 | 1280.919 | 43,759 | 546.0 |
| ParT PyTorchAlpaka, global TF32 allowed, batch 256 | 958.973 | 58,450 | 546.0 |

On the raw-data workload, baseline FP32 Torch ParT is **1.48× faster** than graph FP32. Weight-only INT8 graph is **1.79× faster** than that Torch baseline; this is explicitly a different weight-precision comparison. Its numerical differences are reported below. The stricter NVIDIA TF32 override takes **1280.9 ms/orbit**, and enabling the global TF32 policy takes **959.0 ms/orbit**. These results are specific to the supplied sparse workload and these models.

## Batch and framework-concurrency tuning

Torch maximum batch size and the number of framework streams were tuned separately. A larger number does not guarantee more useful GPU work. Each orbit already supplies tens of thousands of jets, far above the resident-jet capacity of an L4. `maxBlocks=0` submits every graph/compiled jet block; a capped native grid was not used for the primary comparison.

| Torch model | Tested maximum batches | Best batch in synchronized sweep | Producer ms/orbit after five warmups |
|---|---|---:|---:|
| DeepSets | 128, 512, 2048, 4096, 8192, 65536 | 4096 | 4.900 |
| ParT | 128, 256, 512, 1024, 4096, 65536 | 256 | 1235.253 |

These diagnostic producer numbers include preprocessing and an explicit framework synchronization on the tagger. FastTimer's printed per-event values have 0.1 ms resolution. The asynchronous completed-event table above is the main throughput result; it avoids charging GPU waits inconsistently to the producer versus cleanup.

| Model/runtime | Threads = streams | Asynchronous pipeline ms/orbit | GPU MiB above idle |
|---|---:|---:|---:|
| DeepSets graph | 1 | 7.182 | 408.0 |
| DeepSets graph | 4 | 6.210 | 866.0 |
| DeepSets graph | 8 | 6.197 | 1,218.0 |
| DeepSets Torch batch 4096 | 1 | 6.657 | 502.0 |
| DeepSets Torch batch 4096 | 4 | 6.850 | 1,464.0 |
| DeepSets Torch batch 4096 | 8 | 7.099 | 2,712.0 |
| ParT graph FP32 | 1 | 1841.702 | 418.0 |
| ParT graph FP32 | 4 | 1862.635 | 876.0 |
| ParT Torch batch 256 | 1 | 1243.892 | 546.0 |
| ParT Torch batch 256 | 4 | 1342.378 | 1,402.0 |

For DeepSets, repeated fresh-process runs at the graph plateau (four streams) give **6.210, 6.304, 6.246 ms/orbit**, versus **6.657, 6.643, 6.647 ms/orbit** for Torch's best tested concurrency (one stream). The process medians are **6.246 / 6.647 ms**, respectively: a modest **1.06×** observed graph advantage after concurrency tuning. Eight graph streams add memory without a material throughput gain. At one stream, Torch is slightly faster. Treat this small ordering as configuration dependent, particularly with unlocked GPU clocks.

For ParT, compare the one/four-stream results directly: more in-flight orbits increase memory substantially and do not necessarily increase completed-jet throughput. The tested stream count is a throughput parameter, not a measure of how many jet states are simultaneously resident in shared memory.

## Memory interpretation

NVML v1's `used` field includes driver-reserved device memory. Here its idle value is **468.875 MiB**; an idle v2 query reports that amount as `reserved`, with v2 `used=0`. This is why raw v1 peaks are larger than an otherwise idle `nvidia-smi` application footprint. The tables subtract the **own-process run's recorded idle baseline**. See [nvml_memory_definition.json](nvml_memory_definition.json). Other GPU applications would contaminate this device-wide measurement, so the primary asynchronous cases ran in isolation.

At equal one-stream concurrency, graph uses less sampled GPU memory for these workloads: approximately **408 versus 502 MiB** for DeepSets, and **418 versus 546 MiB** for FP32 ParT. INT8 graph uses about **412 MiB**. These include a common no-tagger pipeline footprint of about **280 MiB**. They are peak process/workload comparisons, not a measurement that every byte of the difference is an activation tensor.

The faster DeepSets graph configuration uses four streams and therefore more total memory (**866 MiB**) than single-stream Torch (**502 MiB**). A low-memory inference schedule does not guarantee that the fastest chosen framework configuration has the smallest total footprint.

A whole-orbit Torch ParT batch (`maxBatchSize=65536`, actual size 56,052) takes about **2072 ms in the synchronized producer** and reaches **20,002 MiB above idle** (about 19.53 GiB). It is slower and far larger than tuned batch 256. This total CMSSW footprint must not be equated to the earlier standalone 7.29 GiB incremental live-tensor peak: they measure different scopes, including allocator/queue caches and full producer/event buffers here.

An early synchronized FP32 graph trace (`part_graph_fp32_sync`) overlapped the SCRAM test-plugin build and contains two brief extra device-memory peaks. It is retained and annotated, but its approximately 607 MiB above-idle maximum is **excluded from isolated memory comparisons**. The later isolated graph run is approximately 418 MiB, consistent with its separate validation run. The earlier commentary's comparison based on the transient peak was corrected.

The graph still has no batch-sized global activation arena, but compiler stack/local backing, CUDA context/code, model arrays, SoA inputs/outputs, and allocator pools consume global memory. Host RSS traces are retained in every case JSON. RSS includes loaded library/code pages and host/raw/pinned buffers, not just model storage; it is not peak CUDA memory.

## Why CMSSW changes the standalone ordering

The native CMSSW device binary is compiled differently from the standalone non-RDC CUDA shim. `cuobjdump` of the actual `sm_89` binary reports:

| Production kernel | Registers/thread | Static stack bytes/thread | Shared bytes/block |
|---|---:|---:|---:|
| Graph block | 246 | 496 | 47,840 |
| Graph warp | 254 | 392 | 43,776 |
| Compiled QAT | 255 | 1,824 | 47,872 |

The L4 has **65,536 registers and 102,400 shared bytes per SM**, across **58 SMs**. With 256 threads/block, these register counts permit only **one resident block per SM**; the earlier bridge had 128 registers/thread and could keep two blocks resident. The graph warp implementation therefore has a register bound of eight resident jets per SM in this CMSSW build, versus sixteen in the bridge. Stack figures describe the compiled frame, not an allocation for every queued jet, and `LOCAL:0` in `cuobjdump` does not mean that the nonzero stack or device-call backing is absent. Raw extracted attributes are in [kernel_resources.json](kernel_resources.json) and [device_attributes.json](device_attributes.json).

This is direct evidence that the standalone kernel's resource/occupancy results cannot simply be transferred to the CMSSW plugin. Lower occupancy and different device-call/inlining behavior are plausible contributors to the different timing; these measurements do not isolate the cost of each compiler choice. Producer-level preprocessing/allocation, Torch adapter/batch setup, queue integration, and cleanup also contribute. No Tensor Core instruction counts or per-kernel profiler attribution are claimed here.

The supplied raw workload also has much less useful particle content than the earlier full-16 synthetic test. The graph uses the active extent of each jet and its triangular pair evaluation, whereas the frozen Torch archive retains a full 16-row tensor and square pair representation. This favors native inference particularly at sparse occupancy. INT8 graph therefore wins this raw-data timing comparison even though it lost the earlier full-occupancy inference comparison. Comparing different multiplicities or build modes as though only the adapter changed would give a misleading explanation.

## Input and numerical checks

The first ten inspected orbits each have 56,052 jets, mean **2.65896 valid particles/jet**. **53.18%** have one valid particle and **0.578%** have all sixteen. Their score vectors are identical across those ten events, so they are a repeated benchmark workload, not ten statistically independent physics samples. The complete mask-population histogram is in [summary.json](summary.json). `nJets=16` seeds per BX and the other clustering options follow the baseline. The input file hashes are in [environment.json](environment.json).

Separate validation jobs inspect every **published score**, using identical upstream data and production preprocessing. All scores were finite and nonzero. Each comparison below covers 560,520 paired jets (ten repeated orbits):

| Output pair | Maximum absolute probability difference | Mean absolute probability difference |
|---|---:|---:|
| DeepSets graph FP32 / Torch | 2.39e-07 | 2.823102e-08 |
| ParT graph FP32 / Torch | 0.001567363 | 2.188813e-05 |
| ParT graph weight-only INT8 / Torch | 0.004104614 | 0.0005079263 |
| Graph QAT / compiled QAT | 0.000349611 | 2.428645e-06 |

The ParT INT8 maximum is about 0.41 percentage points on this sample. These checks quantify runtime/conversion differences; they do not establish physics accuracy on an independent validation set. The graph/compiled QAT comparison uses the same shipped QAT checkpoint, which differs from the frozen Torch checkpoint.

**Precision policy:** the baseline Torch cases use FP32 inputs/weights with global matmul and cuDNN TF32 flags disabled. The frozen ParT archive itself contains `_convolution(..., allow_tf32=True)` calls. The additional `--force-ieee` cases set `NVIDIA_TF32_OVERRIDE=0` in the entire process, instead of assuming a global Torch flag overrides every serialized operation. `--tf32` explicitly enables the global policy. Enabling a policy does not prove which matrix instructions every operation used. NVIDIA documents the override in its [cuBLAS 12.9 guide](https://docs.nvidia.com/cuda/archive/12.9.1/cublas/index.html#cublascomputetype-t) and [TF32 library discussion](https://developer.nvidia.com/blog/accelerating-ai-training-with-tf32-tensor-cores/). The numerical table above compares the original baseline policy; the override/TF32 timings are separately labeled and are not accompanied by a new per-jet physics validation.

## Separate shipped-QAT context

| QAT producer, one asynchronous stream | Pipeline ms/orbit | Jets/s | GPU MiB above idle |
|---|---:|---:|---:|
| Compiled ParticleTransformerProducer | 378.313 | 148,163 | 480.0 |
| GraphInferenceProducer, b_kinadd.nnm | 485.982 | 115,338 | 412.0 |

These two use the same QAT blobs as each other. The original matching QAT TorchScript archive is unavailable, so these are **not a same-checkpoint speedup claim against the frozen Torch model**.

The first compiled run processed all 100 events but exited with `SIGABRT` during static-cache destruction. [failed_compiled_cache_teardown.log](failed_compiled_cache_teardown.log) retains that failure. The local [cache fix](compiled_cache_fix.patch) applies the strategy already used by the graph producer: use device-scoped, non-caching model allocations and retain their registry for process lifetime, avoiding late destruction after CUDA/allocator teardown. The numerical kernels and model data are unchanged. The repeated 100-event job, ten-event score validation, and asynchronous job all completed with status zero after the fix. Its bounded model storage is reclaimed when the process exits.

The first attempt to load the frozen CPU model through the default adapter constructor then `Module::to()` also failed because embedded constants stayed on CPU. [failed_frozen_cpu_constants.log](failed_frozen_cpu_constants.log) records it. The benchmark producer uses the already-supported `AlpakaModel(path, queue)` constructor to load-map constants onto the device under `QueueGuard`; it does not modify PyTorchAlpaka.

## Method, source and reproduction

The completed campaign contains **46 successful `cmsRun` jobs and 2,910 framework events/orbits**, plus the two explicitly retained failed diagnostics. All seven score-inspection files passed finite/nonzero checks; final source/runtime hashes and device-binary resource counts were verified. Python syntax, shell syntax, and the cache-fix diff checks passed.

The unchanged production graph/compiled producers are compared with a **test-only PyTorch producer using the same production transform and the real PyTorchAlpaka API**. The baseline's older `SoftJetIdML` transform differs, including barriers inside a divergent nested work loop, so it is not used for the matched main comparison. This is a benchmark producer adaptation, not a claim that the baseline's originally configured 10k model or its different TorchScript input format was measured unchanged.

The [harness reference](../../cmssw_benchmark/README.md) explains every added class/function and provides reproduction commands. [environment.json](environment.json) hashes the baseline scripts, checkpoints, exported graphs, final harness sources, production producer, plugin binaries and raw inputs. [summary.json](summary.json) collects successful cases and comparisons.

Synchronized graph/compiled and DeepSets sweeps use 100 events; ParT Torch tuning uses 20 events, and its whole-event batch uses ten. Separate correctness jobs use ten events. Main asynchronous DeepSets jobs use 100 events with 1/4/8 threads and streams; ParT uses 30 single-stream events and 60 events at four streams. The completion calculation removes the first `max(5,5*streams)` completions and divides the remaining completion count by the elapsed steady-clock time, measured **after event cleanup**. Small DeepSets differences are repeated in fresh processes. Per-event FastTimer text is disabled for completion-timed jobs; ordinary framework reporting remains. The GPU/RSS polling interval is nominally 20 ms, so very short memory spikes can be missed.

Each case retains `.json` settings, command, process exit status, full memory samples, completion throughput, and FastTimer totals; `.log` contains framework output; `.timer.json` is the original FastTimer report; `.completion.jsonl` retains individual completion timestamps. `validate_*.scores.jsonl.gz` contains every inspected score and mask histogram, compressed after the runs. Failed attempts are explicitly prefixed `failed_` and excluded from the successful-case summary.

The source of the framework pipeline is a local RAM-disk raw file. There is no remote storage throughput, NanoAOD output, FP16 inference, CPU/ROCm comparison, TensorRT/AOT/CUDA Graph optimization, or matching Torch-QAT model in these measurements. Clocks were not locked; the figures describe this hardware, source/build revision, model, and multiplicity. Use the raw samples and measured scope rather than treating any one speed ratio as universal.
