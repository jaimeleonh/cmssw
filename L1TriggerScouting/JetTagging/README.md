# Jet tagging at L1 Scouting — ONNX Runtime implementation

This package runs b-jet tagging (ParticleTransformer) on soft jets at L1
scouting, with two implementations living side by side:

| | PyTorch (`SoftJetIdML`) | ONNX Runtime (`SoftJetIdONNX`) |
|---|---|---|
| Model format | TorchScript (`.pt`) | ONNX (`.onnx`) |
| Inference | `PhysicsTools/PyTorchAlpaka` (`torch::AlpakaModel`) | `PhysicsTools/ONNXRuntimeAlpaka` (`ort::AlpakaSession`) |
| Input tensor SoA | `SoftJetInputTensorSoA` (particle-major) | `SoftJetONNXInputTensorSoA` (feature-major) |
| Output | `SoftJetOutputDeviceTensor` | `SoftJetOutputDeviceTensor` (identical) |

Both are alpaka modules (serial + cuda backends) and produce the same output
data format, so they are drop-in interchangeable in the reconstruction
sequences and in `ScPhase2SoftJetOutputTensorToOrbitFlatTable`.

## `SoftJetIdONNX`

Replicates `SoftJetIdML` using ONNXRuntimeAlpaka:

- `plugins/alpaka/TransformKernelONNX.{h,dev.cc}` — builds the
  ParticleTransformer input tensors directly on device from
  `PFCandidateDeviceCollection`, `AssociationMapDevice`, jet/vertex
  `BxLookupDevice` and `VertexDeviceCollection`. The two-pass algorithm and
  the feature arithmetic are identical to the PyTorch
  `ComputeParticleTransformerFeaturesKernel`; the only difference is the
  memory layout of the per-jet tensors (see below).
- `plugins/alpaka/SoftJetIdONNX.cc` — drives the transform, wraps the tensors
  in `cms::Ort::alpakatools::TensorCollection`, runs the session bound to the
  module queue (`AlpakaSession` does the SoA→row-major gather through
  `PackKernels` and the queue synchronisation), extracts the b-jet score, and
  puts a `SoftJetOutputDeviceTensor` into the event.

### Why a transposed SoA

The exported ONNX models expect feature-major tensors, e.g.
`pf_features [batch, 15, n_particles]`. In `SoftJetInputTensorSoA` the eigen
blocks are `[16 particles, 15 features]`, which would arrive at ONNX Runtime
transposed. `DataFormats/L1ScoutingSoA/interface/SoftJetONNXTensorSoA.h`
defines the same content with the eigen blocks transposed
(`[15 features, 16 particles]`, `[4, 16]`, `[1, 16]`), so after the
`PackKernels` gather the contiguous buffer is exactly what the model expects —
no transpose pass, no points column (the model derives pairwise coordinates
from the vectors).

The `softmax [batch, 2]` output is staged in a local two-column SoA
(`SoftJetONNXScoresSoA` in `TransformKernelONNX.h`) and `extractScores`
copies the signal-class score into the final tensor. Class 0 is the b score
(the training uses `labels: [jet_isB, jet_isnotB]`).

### Batch size = 1

The current exports (``model_example_ParticleTransformer.onnx``) contain a
`pair_embed/Reshape` node with a hardcoded `[8, 136]` target shape
(136 = 16·17/2 pairs) that only matches batch size 1, even though the batch
dimension is declared dynamic. `SoftJetIdONNX` therefore defaults to
`maxBatchSize = 1` (one jet per inference call). The batching machinery is
generic — re-exporting the model without that reshape and raising
`maxBatchSize` re-enables batched inference with no code change.

### Usage

```python
from L1TriggerScouting.JetTagging.modules import l1sc_SoftJetIdONNX_alpaka

process.ONNXService = cms.Service("ONNXService")   # required by ONNXRuntimeAlpaka

process.SoftJetIdSC4 = l1sc_SoftJetIdONNX_alpaka(
    alpaka = cms.untracked.PSet(backend = cms.untracked.string(options.backend)),
    model = cms.FileInPath("L1TriggerScouting/JetTagging/data/btagging/softjet_part_example.onnx"),
    pf = cms.InputTag("scPhase2PFRawToDigiAlpaka", "candidates"),
    clusters = cms.InputTag("scPhase2SC4PFAlpaka"),
    jetBxLookup = cms.InputTag("scPhase2SC4PFAlpaka"),
    vertices = cms.InputTag("scPhase2VertexRawToDigiAlpaka"),
    vertexBxLookup = cms.InputTag("scPhase2VertexRawToDigiAlpaka", "bxLookup"),
    maxBatchSize = cms.uint32(1),
)
```

The parameter names match `SoftJetIdML`, so it can be swapped into the
existing `runScoutingPhase2*` configurations in place of
`l1sc_SoftJetIdML_alpaka` / `l1sc_SoftJetIdMLTorchScript_alpaka`.

## Test

`test/runSoftJetTaggingONNXTest.py` is a self-contained end-to-end test:

1. `SoftJetONNXTestSource` builds a synthetic event from explicit parameters:
   5 jets over 3 BX crossings (one empty, one >16-constituent jet exercising
   the truncation) plus one vertex per used BX.
2. `SoftJetTagONNX` runs the transform + ONNX inference.
3. `SoftJetONNXTestSink` copies the scores to the host and compares them with
   the reference values baked into the cfg (tolerance 1e-4); any mismatch
   fails the job.

Run it with

```bash
cmsRun src/L1TriggerScouting/JetTagging/test/runSoftJetTaggingONNXTest.py backend=serial_sync
# or, on a machine with a GPU:
cmsRun src/L1TriggerScouting/JetTagging/test/runSoftJetTaggingONNXTest.py backend=cuda_async
```

The reference scores are produced by `test/generateSoftJetONNXTest.py`, which
replicates the kernel preprocessing step by step in float32 and scores the
same tensors with the Python ONNX Runtime. Re-run it and paste its output
block into the cfg whenever the model or the synthetic event changes:

```bash
python3 src/L1TriggerScouting/JetTagging/test/generateSoftJetONNXTest.py
```

The model used by the test is
`data/btagging/softjet_part_example.onnx` (example ParticleTransformer,
inputs `pf_features [N,15,n_pf]`, `pf_vectors [N,4,n_pf]`, `pf_mask [N,1,n_pf]`,
output `softmax [N,2]`).
