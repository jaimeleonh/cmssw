# CMSSW producer benchmark harness

This harness derives its source, unpacking, vertex decoding, seeded-cone clustering, and event selection from `L1TriggerScouting/Phase2/test/runScoutingPhase2PF_cfg.py`, the configuration invoked by the user's `run_benchmark.sh`. It replaces only the scheduled tagger and the reporting destinations. The baseline files are not edited. It runs on GPU 0, using the run-41 files already present under `/mnt/ramdisk/jleonhol/raw` and a separate FU directory under `/tmp`.

Measured results, numerical validation, and limitations are in [the CMSSW results report](../performance_results/cmssw/README.md). These complete producer measurements are separate from the earlier standalone inference results.

## Files and functions

[run_cfg.py](run_cfg.py) executes the baseline configuration using `runpy.run_path`, keeps its parsed `options`, and selects `BENCH_ENGINE`. `graph` instantiates the unchanged production `GraphInferenceProducer`; `compiled` instantiates `ParticleTransformerProducer`; `torch` instantiates the test-only `BenchmarkPyTorchProducer`; `none` omits tagging. `torch_original` is an optional diagnostic using the existing `SoftJetIdML`, whose preprocessing differs and is excluded from the matched comparison. The script rebuilds the path and schedule, overrides FastTimerService reporting, optionally adds a score inspector, and optionally enables the completion timing service. Observer runs suppress per-event FastTimer text to avoid measuring that log I/O. It defines no functions or classes.

[BenchmarkProducers.cc](../alpaka/BenchmarkProducers.cc) contains:

- `inputDescription(d)`: registers the five common upstream `InputTag` parameters.
- `Inputs`: groups the typed PF, association, jet-BX, vertex, and vertex-BX consumption tokens.
- `BenchmarkPyTorchProducer::BenchmarkPyTorchProducer(p)`: registers consumes/produces, retains the checkpoint path, validates a positive maximum batch size, and explicitly configures strict FP32 or TF32 policy.
- `BenchmarkPyTorchProducer::fillDescriptions`: supplies the common tags, TorchScript file, batch-size limit, and TF32 switch.
- `BenchmarkPyTorchProducer::produce`: obtains the same products as the production graph producer, calls its unchanged `makeParticleTransformerInputs`, allocates/zeros the output collection, and registers mutable SoA records in `TensorCollection` batches. Input order is points, features, four-vectors, mask. It enables inference mode and uses `QueueGuard` for model loading as well as execution. On its first nonempty event it constructs the actual `AlpakaModel(path, queue)`, then calls the actual `forward(queue, inputs, outputs)` for every batch and publishes the device score collection. This tests the real adapter and its output-copy path, not a direct substitute `torch::jit` call. The queue constructor maps frozen tensor constants to the device; loading on CPU and later calling `Module::to()` does not move those constants in the supplied frozen archives. Loading is excluded by warmup. The class owns one model per CMSSW stream, as a stream producer normally does.
- `BenchmarkScoreInspector::BenchmarkScoreInspector(p)`: registers the score input, common upstream inputs, and output JSONL destination.
- `BenchmarkScoreInspector::fillDescriptions`: declares those inputs and destination.
- `BenchmarkScoreInspector::produce`: independently repeats the production transform, copies the input and produced scores to pinned host collections, waits for those copies, then records every score, event identifier, jet count, valid-particle histogram, nonfinite/zero counts, and score sum. It runs only in separate single-stream validation jobs. Its copy, synchronization, transform, and text-writing overhead are excluded from timing jobs.

[BenchmarkInputs.dev.cc](../alpaka/BenchmarkInputs.dev.cc) includes the production `TransformKernel.dev.cc` unchanged so the test library can link the same transform. It introduces no new algorithm or named definitions.

[BenchmarkTimingService.cc](../BenchmarkTimingService.cc) supplies completion timestamps for asynchronous and concurrent runs:

- `Completion`: stores event ID, stream index, and steady-clock nanoseconds.
- The constructor registers `complete` on `watchPostClearEvent` and `write` on `watchPostEndJob`, retaining the destination path.
- `complete(context)` timestamps the event after cleanup, protecting the vector with a mutex. This includes asynchronous work waited for during framework cleanup. A producer's host return alone is not a completed-event measurement.
- `write()` writes the retained records after the event loop, avoiding per-event file I/O in the measured interval.

[test/BuildFile.xml](../BuildFile.xml) builds two test plugin libraries: the Alpaka CUDA/serial producers plus unchanged transform, and the ordinary framework timing service. It adds Torch/PyTorchAlpaka dependencies to the test library; it does not change the production plugin dependencies.

[run.py](run.py) drives each isolated process:

- `Memory`: the ctypes layout of NVML's device memory record (`total`, `free`, `used` bytes).
- `parse_events(log)`: extracts per-event CPU and real milliseconds from FastTimerService's module table, separating it from the following path table and run/job averages. Printed values have 0.1 ms resolution.
- `main()`: parses one case, sets `BENCH_*` variables, launches `cmsRun`, and polls NVML plus `/proc/PID/status` at a nominal 20 ms interval. It enforces a timeout, retains logs and FastTimer JSON, excludes five initial events from diagnostic averages, records the complete memory trace and process status, and verifies the completed-event count. `--force-ieee` sets `NVIDIA_TF32_OVERRIDE=0`; this additionally controls the supplied archive’s convolution calls with embedded `allow_tf32=True`, which a global Torch flag alone may not override. With `--observer`, per-event FastTimer logging is disabled; throughput is calculated from the completion timestamp of warmup event `max(5, 5*streams)` to the final completion. The number of intervening completions divided by elapsed time gives sustained orbit throughput. Startup/model loading is outside that interval. It also retains aggregate FastTimer data, whose module totals include startup.

[run_matrix.sh](run_matrix.sh) runs the synchronized batch sweep, production graph/compiled QAT context, and independent validation cases. `run_case(name, args...)` launches a named case and prints its summary, skipping an already recorded case. [run_followup.sh](run_followup.sh) selects the best measured batch using `best_batch(model)`, then tests normal asynchronous execution, several framework stream counts, TF32, and the separate QAT paths. Its `run_case` enables the completion observer. It also repeats the small DeepSets throughput difference in fresh processes, using four streams for graph and one for Torch. Existing result files should be moved aside before deliberately repeating a case under the same name.

[summarize.py](summarize.py) writes a compact summary. `scores(path)` reads plain or gzip JSONL inspection records; `compare(directory, left, right)` verifies event/jet/mask alignment and finite scores, then computes maximum/mean absolute error and RMSE; `main()` gathers successful process cases, converts byte counters to MiB, summarizes multiplicity, pairs the matched score files, and writes `summary.json`.

## Reproduction

In the package directory with CMSSW and the L4 available:

```bash
eval "$(scram runtime -sh)"
scram b -j 8
# Produce the matched FP32 and weight-only INT8 .nnm files in /tmp if needed:
bash test/build_runtime_benchmark.sh /tmp
python3 test/benchmark_runtime_comparison.py --model deepsets --soa \
  --sizes 512 --batches 512 --repeats 1 --output /tmp/export_deepsets.json
python3 test/benchmark_runtime_comparison.py --model part --soa \
  --sizes 512 --batches 256 --repeats 1 --output /tmp/export_part.json
python3 test/benchmark_runtime_comparison.py --model part --soa --quantize \
  --sizes 512 --batches 256 --repeats 1 --output /tmp/export_part_int8.json
bash test/cmssw_benchmark/run_matrix.sh
bash test/cmssw_benchmark/run_followup.sh
```

Example single asynchronous run:

```bash
python3 test/cmssw_benchmark/run.py --engine torch --model part --batch 256 \
  --events 100 --observer --output /tmp/part_torch_cmssw
```

`maxBlocks=0` submits all jet blocks for graph/compiled paths. Torch's `maxBatchSize` is separately tuned; setting it to at least 56,052 processes the whole supplied orbit in one batch. Framework streams are another independent dimension of concurrency. A single orbit already provides many more jets than the L4 can keep resident at once.
