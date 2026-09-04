import FWCore.ParameterSet.Config as cms

particleTransformer = cms.EDProducer(
    "alpaka/ParticleTransformerProducer@alpaka",
    pf=cms.InputTag(""),
    clusters=cms.InputTag(""),
    jetBxLookup=cms.InputTag(""),
    vertices=cms.InputTag(""),
    vertexBxLookup=cms.InputTag(""),
    weightsFile=cms.string("PhysicsTools/ParticleTransformerAlpaka/data/model_weights.int8.bin"),
    paramsFile=cms.string("PhysicsTools/ParticleTransformerAlpaka/data/model_params.fp32.bin"),
    # Zero launches one block per jet, which is the fastest configuration.
    # A non-zero value caps the grid and makes each block loop over several
    # jets; it is only useful for occupancy studies.
    maxBlocks=cms.uint32(0),
    # Diagnostic only: moves asynchronous inference time from framework
    # cleanup into this producer in FastTimerService reports.
    synchronizeForTiming=cms.bool(False),
)
