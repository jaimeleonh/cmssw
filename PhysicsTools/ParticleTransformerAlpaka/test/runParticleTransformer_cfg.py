import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.alpaka_cff import alpaka

process = cms.Process("PART", alpaka)
process.load("PhysicsTools.ParticleTransformerAlpaka.particleTransformerProducer_cfi")
process.particleTransformer.pf = "REPLACE_WITH_PF_CANDIDATES"
process.particleTransformer.clusters = "REPLACE_WITH_JET_CANDIDATE_ASSOCIATION"
process.particleTransformer.jetBxLookup = "REPLACE_WITH_JET_BX_LOOKUP"
process.particleTransformer.vertices = "REPLACE_WITH_VERTICES"
process.particleTransformer.vertexBxLookup = "REPLACE_WITH_VERTEX_BX_LOOKUP"
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))
process.path = cms.Path(process.particleTransformer)
