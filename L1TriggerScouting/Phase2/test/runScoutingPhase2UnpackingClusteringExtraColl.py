import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("TestUnpacking", eras.Phase2C17I13M9)

from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
def LazyVar(expr, valtype, doc=None, precision=-1):
    return Var(expr, valtype, doc, precision, lazyEval=True)

# from L1TriggerScouting.Phase2.options_cff import options, VarParsing
# options.register ('njets',
#                   16, 
#                   VarParsing.VarParsing.multiplicity.singleton,
#                   VarParsing.VarParsing.varType.int,         
#                   'Number of jet seeds to reconstruct with SeededCone'
# )
# options.register ('minSeedPt',
#                   0.0, 
#                   VarParsing.VarParsing.multiplicity.singleton,
#                   VarParsing.VarParsing.varType.float,
#                   'Minimum pt cut for seeded-cone jet seeds')
# options.register ('jetR',
#                   0.4, 
#                   VarParsing.VarParsing.multiplicity.singleton,
#                   VarParsing.VarParsing.varType.float,
#                   'Jet radius')
# options.register ('jetReFitR',
#                   0.0, 
#                   VarParsing.VarParsing.multiplicity.singleton,
#                   VarParsing.VarParsing.varType.float,
#                   'Jet radius (after ReFit)')
# options.register ('dumpClusters',
#                   False, 
#                   VarParsing.VarParsing.multiplicity.singleton,
#                   VarParsing.VarParsing.varType.bool,         
#                   'Dump clusters to options.outFile')

# options.parseArguments()

# process.options = cms.untracked.PSet(
#     numberOfThreads = cms.untracked.uint32(options.numThreads),
#     numberOfStreams = cms.untracked.uint32(options.numFwkStreams),
#     numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
#     wantSummary = cms.untracked.bool(True)
# )

# logging configuration
process.load('Configuration.StandardSequences.Services_cff')
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("Configuration.StandardSequences.Accelerators_cff")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.options.wantSummary = cms.untracked.bool(True)
process.maxEvents.input = cms.untracked.int32(10)

process.load('Configuration.Geometry.GeometryExtendedRun4D110Reco_cff')
process.load('Configuration.Geometry.GeometryExtendedRun4D110_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff') # needed to read HCal TPs
process.load('SimCalorimetry.HGCalSimProducers.hgcalDigitizer_cfi') # needed for HGCAL_noise_fC
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '141X_mcRun4_realistic_v3', '')

process.load('L1Trigger.Phase2L1ParticleFlow.l1ctLayer1_cff')
process.load('L1Trigger.L1TTrackMatch.l1tGTTInputProducer_cfi')
process.load('L1Trigger.L1TTrackMatch.l1tTrackSelectionProducer_cfi')
process.l1tTrackSelectionProducer.processSimulatedTracks = False # these would need stubs, and are not used anyway
process.load('L1Trigger.VertexFinder.l1tVertexProducer_cfi')
from L1Trigger.Configuration.SimL1Emulator_cff import l1tSAMuonsGmt
process.l1tSAMuonsGmt = l1tSAMuonsGmt.clone()
from L1Trigger.L1CaloTrigger.l1tPhase2L1CaloEGammaEmulator_cfi import l1tPhase2L1CaloEGammaEmulator
process.l1tPhase2L1CaloEGammaEmulator = l1tPhase2L1CaloEGammaEmulator.clone()
from L1Trigger.L1CaloTrigger.l1tPhase2CaloPFClusterEmulator_cfi import l1tPhase2CaloPFClusterEmulator
process.l1tPhase2CaloPFClusterEmulator = l1tPhase2CaloPFClusterEmulator.clone()
from L1Trigger.L1CaloTrigger.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator_cfi import l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator = l1tPhase2GCTBarrelToCorrelatorLayer1Emulator.clone()

# process.L1TInputTask = cms.Task(
#    process.l1tSAMuonsGmt,
#    process.l1tPhase2L1CaloEGammaEmulator,
#    process.l1tPhase2CaloPFClusterEmulator,
#    process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
# )

from IOPool.Input.modules import PoolSource
process.source = PoolSource(
    fileNames = [
        #"file:" + cms.FileInPath("L1TriggerScouting/TauTagging/data/pfOnly.root").value()
        #'file:/afs/cern.ch/work/g/gpetrucc/l1p2/l1scout-cmssw/CMSSW_16_0_0_pre1/src/L1TriggerScouting/TauTagging/data/pfOnly.root'
        #'/store/cmst3/group/l1tr/vcamagni/L1TauID/DATA/FPinputs/m40/4STEPS/142Xv0/inputs140X_7099346_0.root'
        '/store/cmst3/group/l1tr/FastPUPPI/15_1_X/fpinputs_140X/v1/caseC_m220_67/4STEPS/151Xv0/inputs151X_14682953_1699.root'
    ]
)

process.load('Configuration.StandardSequences.SimL1Emulator_cff')

from L1Trigger.Phase2L1ParticleFlow.L1BJetProducer_cff import l1tBJetProducerPuppiCorrectedEmulator
process.bjetsNoMinPt = l1tBJetProducerPuppiCorrectedEmulator.clone(maxJets = 16, minPt = 0)

from L1Trigger.Phase2L1ParticleFlow.l1tSC4NGJetProducer_cfi import l1tSC4NGJetProducer
process.NGJetNoMinPt = l1tSC4NGJetProducer.clone(minPt = cms.double(0.))

process.runPF = cms.Task(
    process.l1tTkMuonsGmt,
    process.l1tSAMuonsGmt,
    process.l1tPhase2L1CaloEGammaEmulator,
    process.l1tPhase2CaloPFClusterEmulator,
    process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator,
    process.l1tGTTInputProducer,
    process.l1tTrackSelectionProducer,
    process.l1tVertexFinderEmulator, 
    process.L1TLayer1TaskInputsTask,
    process.L1TLayer1Task,
    process.l1tLayer2EG,
    process.l1tLayer1BarrelExtended,
    process.l1tLayer1HGCalExtended,
    process.l1tLayer1HGCalNoTK,
    process.l1tLayer1HF,
    process.l1tLayer1Extended,
    process.L1TPFJetsEmulationTask,
    process.L1TPFJetsExtendedTask,
    process.L1TBJetsTask,
    process.NGJetNoMinPt,
    process.bjetsNoMinPt
)

process.p_pf = cms.Path(
    process.runPF
)


# cms.Path( 
#         process.l1tGTTInputProducer +
#         process.l1tTrackSelectionProducer +
#         process.l1tVertexFinderEmulator +
#         process.l1tLayer1BarrelExtended + 
#         process.l1tLayer1HGCalExtended +
#         process.l1tLayer1HGCalNoTK +
#         process.l1tLayer1HF +
#         process.l1tLayer1Extended +
#         process.L1TPFJetsEmulationTask,
#         process.L1TPFJetsExtendedTask,
#         process.L1TBJetsTask
# )
# process.runPF.associate(process.L1TLayer1TaskInputsTask)
# process.runPF.associate(process.L1TInputTask)

process.pfBarrelTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
                name = cms.string("L1PFBarrel"),
                src = cms.InputTag("l1tLayer1BarrelExtended:PF"),
                cut = cms.string(""),
                doc = cms.string(""),
                singleton = cms.bool(False), # the number of entries is variable
                extension = cms.bool(False), # this is the main table
                variables = cms.PSet(
                    pt  = Var("pt",  float, precision=16),
                    eta  = Var("eta", float, precision=16),
                    phi = Var("phi", float, precision=16),
                    pdgId = Var("pdgId", int,),
                    z0 = Var("vz", float, precision=16),
                    dxy = LazyVar("dxy", float, precision=16),
                    quality = LazyVar("hwQual", int),
                    puppiw = LazyVar("puppiWeight", float, precision=16),
                )
)
process.pfEndcapTable = process.pfBarrelTable.clone(
                name = cms.string("L1PFEndcap"),
                src = cms.InputTag("l1tLayer1HGCalExtended:PF"),
)
process.puppiTable = process.pfBarrelTable.clone(
                name = cms.string("L1Puppi"),
                src = cms.InputTag("l1tLayer1Extended:Puppi"),
)


process.packPFBarrel = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1BarrelExtended:PF"),
    fedIDs = cms.vuint32(1,2,3),
    splitFactor = cms.uint32(3),
    scoutingHeader = cms.bool(True)
)

process.packPFEndcap = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1HGCalExtended:PF"),
    fedIDs = cms.vuint32(4,5),
    splitFactor = cms.uint32(2),
    scoutingHeader = cms.bool(True)
)

process.packPuppi = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1Extended:Puppi"),
    fedIDs = cms.vuint32(0),
    splitFactor = cms.uint32(1),
    scoutingHeader = cms.bool(True)
)

from L1TriggerScouting.Phase2.modules import (
    l1sc_L1TScPhase2PuppiRawToDigi_alpaka, l1sc_L1TScPhase2SCJets_alpaka
)

process.unpackPFBarrelAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
    #alpaka = cms.untracked.PSet( backend = cms.untracked.string("serial_sync") ),
    streams = process.packPFBarrel.fedIDs,
    splitFactor = process.packPFBarrel.splitFactor,
    src = cms.InputTag('packPFBarrel'),
)
process.unpackPFEndcapAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
    #alpaka = cms.untracked.PSet( backend = cms.untracked.string("serial_sync") ),
    streams = process.packPFEndcap.fedIDs,
    splitFactor = process.packPFEndcap.splitFactor,
    src = cms.InputTag('packPFEndcap'),
)

process.unpackPFBarrel = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPFBarrel.fedIDs,
    splitFactor = process.packPFBarrel.splitFactor,
    src = cms.InputTag('packPFBarrel'),
)
process.unpackPFEndcap = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPFEndcap.fedIDs,
    splitFactor = process.packPFEndcap.splitFactor,
    src = cms.InputTag('packPFEndcap'),
)
process.unpackPuppi = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPuppi.fedIDs,
    splitFactor = process.packPuppi.splitFactor,
    src = cms.InputTag('packPuppi'),
)

process.scPFBarrelStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPFBarrel"),
    name = cms.string("L1PFUnpackBarrel"),
    doc = cms.string("L1PF candidates from Barrel, unpacked to OrbitCollection"),
)
process.scPFEndcapStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPFEndcap"),
    name = cms.string("L1PFUnpackEndcap"),
    doc = cms.string("L1PF candidates from Endcap, unpacked to OrbitCollection"),
)
process.scPFEndcapStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPFEndcap"),
    name = cms.string("L1PFUnpackEndcap"),
    doc = cms.string("L1PF candidates from Endcap, unpacked to OrbitCollection"),
)
process.scPuppiStructToTable = cms.EDProducer("ScPuppiToOrbitFlatTable",
    src = cms.InputTag("unpackPuppi"),
    name = cms.string("L1PuppiUnpack"),
    doc = cms.string("L1Puppi candidates, unpacked to OrbitCollection"),
)
process.scPFBarrelSoAToTable = cms.EDProducer("PFSoAToOrbitFlatTable",
    srcBx = cms.InputTag("unpackPFBarrelAlpaka"),
    srcPF = cms.InputTag("unpackPFBarrelAlpaka"),
    name = cms.string("L1PFUnpackAlpakaBarrel"),
    doc = cms.string("L1PF candidates from Barrel, unpacked to Alpaka SoA"),
)
process.scPFEndcapSoAToTable = cms.EDProducer("PFSoAToOrbitFlatTable",
    srcBx = cms.InputTag("unpackPFEndcapAlpaka"),
    srcPF = cms.InputTag("unpackPFEndcapAlpaka"),
    name = cms.string("L1PFUnpackAlpakaEndcap"),
    doc = cms.string("L1PF candidates from Endcap, unpacked to Alpaka SoA"),
)

process.genJetsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("ak4GenJetsNoNu"),
        # cut  = cms.string("pt > 10 && abs(eta) < 5.0"),
        cut  = cms.string(""),
        name = cms.string("GenJet"),
        doc = cms.string("GenJets (ak4)"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float),
            phi = Var("phi", float),
            eta  = Var("eta", float),
            mass  = Var("mass", float),
        )
)

process.genPartTable = cms.EDProducer("SimpleGenParticleFlatTableProducer",
    src = cms.InputTag("genParticles"),
    cut = cms.string("pt > 0.5 && abs(eta) < 5.5"),
    #cut = cms.string("status == 1 && pt > 0.5 && abs(eta) < 5.5"),
    name = cms.string("GenPart"),
    doc = cms.string("gen particles (pt > 0.5 GeV, |eta| < 5.5)"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table
    variables = cms.PSet(
        pt  = Var("pt",  float),
        phi = Var("phi", float),
        eta  = Var("eta", float),
        mass  = Var("mass", float),
        status  = Var("status", int),
        z0   = Var("vz",  float, doc="Production point along the beam axis"),
        dxy   = Var("vertex.Rho",  float, doc="transverse distance of production point from the beam axis"),
        charge  = Var("charge", int, doc="electric charge"),
        pdgId  = Var("pdgId", int, doc="pdgId code"),
        isPrompt  = Var("statusFlags().isPrompt()", int, doc="Prompt"),
        genPartIdxMother = Var("?numberOfMothers>0?motherRef(0).key():-1", "int16", doc="index of the mother particle"),
    )
)

process.p_gen = cms.Path(
    process.genJetsTable +
    process.genPartTable
)

# Needs FastPUPPI/NtupleProducer
#process.l1VertexTable = cms.EDProducer("VertexWordFlatTableProducer",
#        name = cms.string("L1Vtx"),
#        cut  = cms.string(""),
#        src = cms.InputTag("l1tVertexFinderEmulator","L1VerticesEmulation"),
#        doc = cms.string("Primary vertices reconstructed by L1T"),
#        singleton = cms.bool(False), # the number of entries is variable
#        extension = cms.bool(False), # this is the main table
#        variables = cms.PSet(
#            sumpt = Var("pt",  float,precision=10),
#            z     = Var("z0",  float,precision=16),
#        )
#)
#process.p_vtx = cms.Path(
#    process.l1VertexTable
#)

process.puppiJetsTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        src = cms.InputTag("l1tSC4PFL1PuppiExtendedCorrectedEmulator"),
        cut  = cms.string("pt > 0"),
        name = cms.string("PuppiJet"),
        doc = cms.string("Puppi Jets reconstructed by L1T (extended tracking)"),
        singleton = cms.bool(False), # the number of entries is variable
        extension = cms.bool(False), # this is the main table
        variables = cms.PSet(
            pt  = Var("pt",  float),
            phi = Var("phi", float),
            eta  = Var("eta", float),
            mass  = Var("mass", float),
            numberOfDaughters = Var("numberOfDaughters()", int, doc="number of constituents"),
        ),
        externalVariables = cms.PSet(
            btagScore = ExtVar(cms.InputTag("bjetsNoMinPt", "L1PFBJets"), float, doc="b-tag score (L1 model)"),
        )
)

process.p_puppijet = cms.Path(
    process.puppiJetsTable
)

process.p = cms.Path(
    process.packPuppi +
    #process.packPFBarrel + process.packPFEndcap + process.packPuppi +
    process.unpackPuppi +
    #process.unpackPFBarrel + process.unpackPFEndcap + process.unpackPuppi +
    #process.unpackPFBarrelAlpaka + process.unpackPFEndcapAlpaka +
    process.puppiTable# +
    #process.pfBarrelTable + process.pfEndcapTable + process.puppiTable +
    #process.scPuppiStructToTable #+
    #process.scPFBarrelStructToTable + process.scPFEndcapStructToTable + process.scPuppiStructToTable #+
    #process.scPFBarrelSoAToTable + process.scPFEndcapSoAToTable
)

process.unpackPuppiAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
    #alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
    streams = process.packPuppi.fedIDs,
    splitFactor = process.packPuppi.splitFactor,
    src = cms.InputTag('packPuppi'),
    environment = cms.untracked.int32(0),
)

process.scPhase2SC4PFAlpaka = l1sc_L1TScPhase2SCJets_alpaka(
    #alpaka = cms.untracked.PSet( backend = cms.untracked.string(options.backend) ),
    src = cms.InputTag("unpackPuppiAlpaka"),
    rParam = cms.double(0.4),
    nJets = cms.uint32(999),
)
process.p_sc4Alpaka = cms.Path(
    process.unpackPuppiAlpaka +
    process.scPhase2SC4PFAlpaka
)

process.dumpJetsTable = cms.EDProducer("ClusterObjSoAToNanoaodFlatTable",
    srcClusters = cms.InputTag("scPhase2SC4PFAlpaka"),
    name = cms.string("SC4AllIterJet"),
    doc = cms.string(""),
)

process.clusterJetIndexTable = cms.EDProducer("ClusterMapperSoAToNanoaodFlatTable",
    srcClusters = cms.InputTag("scPhase2SC4PFAlpaka"),
    name = process.puppiTable.name,
    clustering_name = process.dumpJetsTable.name,
    doc = cms.string(""),
)

process.p_dump = cms.Path(process.dumpJetsTable + process.clusterJetIndexTable)

process.sc4NGJetTable = cms.EDProducer("SimpleTriggerL1PFJetFlatTableProducer",
    src = cms.InputTag('NGJetNoMinPt', 'l1tSC4NGJets'),
    name = cms.string("L1puppiJetSC4NG"),
    doc = cms.string("NextGeneration SeededCone 0.4 Puppi jet with multi-class tagging, origin: Correlator"),
    variables = cms.PSet(
        pt  = Var("pt",  float),
        phi = Var("phi", float),
        eta  = Var("eta", float),
        mass  = Var("mass", float),
        udsTagScore = Var('getTagScore("uds")', float),
        gTagScore = Var('getTagScore("g")', float),
        bTagScore = Var('getTagScore("b")', float),
        cTagScore = Var('getTagScore("c")', float),
        tau_pTagScore = Var('getTagScore("tau_p")', float),
        tau_nTagScore = Var('getTagScore("tau_n")', float),
        eTagScore = Var('getTagScore("e")', float),
        muTagScore = Var('getTagScore("mu")', float),
        numberOfDaughters = Var("numberOfDaughters()", int, doc="number of constituents"),
    )
)

process.p_ng_dump = cms.Path(process.sc4NGJetTable)

process.outStandard = cms.OutputModule("NanoAODOutputModule",
    fileName = cms.untracked.string("l1Nano.root"),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
    outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string("ZLIB"),
)
# process.outScout = cms.OutputModule("OrbitNanoAODOutputModule",
#     fileName = cms.untracked.string("scoutingNano.root"),
#     SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
#     skipEmptyBXs = cms.bool(True),
#     outputCommands = cms.untracked.vstring("drop *", "keep l1ScoutingRun3OrbitFlatTable_*_*_*"),
#     compressionLevel = cms.untracked.int32(4),
#     compressionAlgorithm = cms.untracked.string("ZLIB"),
# )
process.e = cms.EndPath(process.outStandard)# + process.outScout)
