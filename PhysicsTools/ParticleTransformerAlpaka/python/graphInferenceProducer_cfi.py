import FWCore.ParameterSet.Config as cms

# The model-agnostic producer.  `modelFile` is the whole network -- layer
# sequence, shapes, buffer layout and weights -- so pointing this at a
# different .nnm swaps the model without rebuilding anything.
#
# Produce one with:
#   python3 scripts/export_model.py particle_transformer <checkpoint>.pt \
#       --output data/<name>.nnm --report data/<name>.export.json
# The type is the plugin name the producer registers, which is
# `DEFINE_FWK_ALPAKA_MODULE(GraphInferenceProducer)` in the accelerator
# namespace -- not a path.  `edmPluginDump | grep GraphInference` is the way to
# check it if this ever stops loading.
graphInference = cms.EDProducer(
    "GraphInferenceProducer@alpaka",
    pf=cms.InputTag(""),
    clusters=cms.InputTag(""),
    jetBxLookup=cms.InputTag(""),
    vertices=cms.InputTag(""),
    vertexBxLookup=cms.InputTag(""),
    modelFile=cms.string("PhysicsTools/ParticleTransformerAlpaka/data/b_kinadd.nnm"),
    # Zero launches one block per jet, which is the fastest configuration.
    # A non-zero value caps the grid and makes each block loop over several
    # jets; it is only useful for occupancy studies.
    maxBlocks=cms.uint32(0),
    # How jets map onto the device.  "auto" runs a small network (one whose
    # activations fit the per-warp arena, e.g. a DeepSets) one jet per warp,
    # eight per block, on CUDA -- several times faster -- and everything else,
    # including every ParT, one jet per block.  "block" and "warp" force one;
    # "warp" is an error if the model does not fit or the back-end is not CUDA.
    layout=cms.string("auto"),
    # Diagnostic only: moves asynchronous inference time from framework
    # cleanup into this producer in FastTimerService reports.
    synchronizeForTiming=cms.bool(False),
)
