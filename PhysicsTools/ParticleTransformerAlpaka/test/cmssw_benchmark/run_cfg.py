"""Derivative of Phase2/test/runScoutingPhase2PF_cfg.py; selected by BENCH_* env."""
import os
import runpy
import FWCore.ParameterSet.Config as cms

base = os.environ["CMSSW_BASE"] + "/src/L1TriggerScouting/Phase2/test/runScoutingPhase2PF_cfg.py"
ns = runpy.run_path(base)
process, options = ns["process"], ns["options"]
engine = os.environ.get("BENCH_ENGINE", "graph")
model = os.environ.get("BENCH_MODEL", "deepsets")
sync = os.environ.get("BENCH_SYNC", "0") == "1"
common = dict(
    alpaka=cms.untracked.PSet(backend=cms.untracked.string("cuda_async"),
                             synchronize=cms.untracked.bool(sync)),
    pf=cms.InputTag("scPhase2PFRawToDigiAlpaka", "candidates"),
    clusters=cms.InputTag("scPhase2SC4PFAlpaka"),
    jetBxLookup=cms.InputTag("scPhase2SC4PFAlpaka"),
    vertices=cms.InputTag("scPhase2VertexRawToDigiAlpaka", "candidates"),
    vertexBxLookup=cms.InputTag("scPhase2VertexRawToDigiAlpaka", "bxLookup"))
if engine == "graph":
    process.SoftJetIdSC4 = cms.EDProducer("GraphInferenceProducer@alpaka", **common,
        modelFile=cms.string(os.environ.get("BENCH_GRAPH", "/tmp/benchmark_%s_float.nnm" % model)),
        maxBlocks=cms.uint32(options.batchSize), layout=cms.string("auto"))
elif engine == "compiled":
    process.SoftJetIdSC4 = cms.EDProducer("ParticleTransformerProducer@alpaka", **common,
        maxBlocks=cms.uint32(options.batchSize))
elif engine in ("torch", "torch_original"):
    name = "l1sc::SoftJetIdML@alpaka" if engine == "torch_original" else "BenchmarkPyTorchProducer@alpaka"
    extras = {} if engine == "torch_original" else dict(tf32=cms.bool(os.environ.get("BENCH_TF32", "0") == "1"))
    modelname = "softjet_deepsets_cpu.pt" if model == "deepsets" else "softjet_part_unique_cpu.pt"
    process.SoftJetIdSC4 = cms.EDProducer(name, **common,
        model=cms.FileInPath("L1TriggerScouting/JetTagging/data/" + modelname),
        maxBatchSize=cms.uint32(options.batchSize), **extras)
elif engine != "none":
    raise ValueError(engine)
# Rebuild the scheduled path so the selected producer replaces the original reference.
process.p_sc4Alpaka = cms.Path(process.scPhase2PFRawToDigiAlpaka +
    process.scPhase2VertexRawToDigiAlpaka + process.goodOrbitsByNBX + process.scPhase2SC4PFAlpaka)
if engine != "none":
    process.p_sc4Alpaka += process.SoftJetIdSC4
if os.environ.get("BENCH_VALIDATE", "0") == "1":
    process.validateScores = cms.EDProducer("BenchmarkScoreInspector@alpaka",
        alpaka=cms.untracked.PSet(backend=cms.untracked.string("cuda_async")),
        src=cms.InputTag("SoftJetIdSC4"), clusters=common["clusters"],
        pf=common["pf"], jetBxLookup=common["jetBxLookup"],
        vertices=common["vertices"], vertexBxLookup=common["vertexBxLookup"],
        output=cms.string(os.environ["BENCH_OUTPUT"] + ".scores.jsonl"))
    process.p_sc4Alpaka += process.validateScores
process.schedule = cms.Schedule(process.p_sc4Alpaka)
process.FastTimerService.jsonFileName = cms.untracked.string(os.environ["BENCH_OUTPUT"] + ".timer.json")
process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.options.accelerators = cms.untracked.vstring("gpu-nvidia")

process.FastTimerService.printEventSummary = cms.untracked.bool(True)
process.FastTimerService.enableDQM = cms.untracked.bool(False)
process.MessageLogger.cerr.FastReport = cms.untracked.PSet(limit=cms.untracked.int32(1000000))
if os.environ.get("BENCH_OBSERVER", "0") == "1":
    process.FastTimerService.printEventSummary = cms.untracked.bool(False)
    process.BenchmarkTimingService = cms.Service("BenchmarkTimingService",
        output=cms.string(os.environ["BENCH_OUTPUT"] + ".completion.jsonl"))
