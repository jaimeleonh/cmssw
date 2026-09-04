import FWCore.ParameterSet.Config as cms
from Configuration.ProcessModifiers.alpaka_cff import alpaka

process = cms.Process("PARTFP32", alpaka)
process.load("PhysicsTools.ParticleTransformerAlpakaFP32.particleTransformerFP32Producer_cfi")
process.particleTransformerFP32.pf = "REPLACE_WITH_PF_CANDIDATES"
process.particleTransformerFP32.clusters = "REPLACE_WITH_JET_CANDIDATE_ASSOCIATION"
process.particleTransformerFP32.jetBxLookup = "REPLACE_WITH_JET_BX_LOOKUP"
process.particleTransformerFP32.vertices = "REPLACE_WITH_VERTICES"
process.particleTransformerFP32.vertexBxLookup = "REPLACE_WITH_VERTEX_BX_LOOKUP"
process.source = cms.Source("PoolSource", fileNames=cms.untracked.vstring())
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))
process.path = cms.Path(process.particleTransformerFP32)
