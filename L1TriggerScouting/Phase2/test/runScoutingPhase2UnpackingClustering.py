import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("TestUnpacking", eras.Phase2C17I13M9)

from PhysicsTools.NanoAOD.common_cff import Var, ExtVar
def LazyVar(expr, valtype, doc=None, precision=-1):
    return Var(expr, valtype, doc, precision, lazyEval=True)

from L1TriggerScouting.Phase2.options_cff import options, VarParsing
options.register ('njets',
                  16, 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.int,         
                  'Number of jet seeds to reconstruct with SeededCone'
)
options.register ('minSeedPt',
                  0.0, 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.float,
                  'Minimum pt cut for seeded-cone jet seeds')
options.register ('jetR',
                  0.4, 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.float,
                  'Jet radius')
options.register ('dumpClusters',
                  False, 
                  VarParsing.VarParsing.multiplicity.singleton,
                  VarParsing.VarParsing.varType.bool,
                  'Dump clusters to options.outFile')

# Additional algorithms
# canonical names:
#   SCGreedy                = iterative greedy seeded cone
#   SCNMS                   = old non-iterative seeded cone with split radii:
#                             RSeed for seed finding + old centroid accumulation,
#                             RClu for final nearest-axis assignment
#   SCNMSWeighted           = same old non-iterative seeded cone logic, but final
#                             assignment uses old weighted metric
#   SCNMSWeightedMultiIter  = newer weighted NMS-style seeded cone with explicit
#                             centroid iterations and separate RCen
#   LinkTree                = link-based clustering

options.register('jetAlgo',
                 'SCGreedy',
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.string,
                 'Jet algorithm: SCGreedy | SCNMS | SCNMSWeighted | SCNMSWeightedMultiIter | LinkTree')
options.register('RSeed',
                 0.3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Seed-finding radius; for SCNMS/SCNMSWeighted also used for old centroid accumulation')
options.register('RCen',
                 0.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Centroid radius (used only by SCNMSWeightedMultiIter)')
options.register('RClu',
                 0.4,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Final particle-assignment radius for SCNMS / SCNMSWeighted / SCNMSWeightedMultiIter')
options.register('RLink',
                 0.3,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'Link radius for LinkTree')
options.register('alphaSeed',
                 2.0,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.float,
                 'alpha in score = dr2 / pt_seed^alpha (used by SCNMSWeightedMultiIter)')
options.register('nCentroidIters',
                 1,
                 VarParsing.VarParsing.multiplicity.singleton,
                 VarParsing.VarParsing.varType.int,
                 'Number of centroid iterations for SCNMSWeightedMultiIter: 0=use seed axis, >0 refine axis')

options.parseArguments()

process.options = cms.untracked.PSet(
    numberOfThreads = cms.untracked.uint32(options.numThreads),
    numberOfStreams = cms.untracked.uint32(options.numFwkStreams),
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    wantSummary = cms.untracked.bool(True)
)

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

process.L1TInputTask = cms.Task(
    process.l1tSAMuonsGmt,
    process.l1tPhase2L1CaloEGammaEmulator,
    process.l1tPhase2CaloPFClusterEmulator,
    process.l1tPhase2GCTBarrelToCorrelatorLayer1Emulator
)

from IOPool.Input.modules import PoolSource
process.source = PoolSource(
    fileNames = [
        #"file:" + cms.FileInPath("L1TriggerScouting/TauTagging/data/pfOnly.root").value()
        #'file:/afs/cern.ch/work/g/gpetrucc/l1p2/l1scout-cmssw/CMSSW_16_0_0_pre1/src/L1TriggerScouting/TauTagging/data/pfOnly.root'
        #'/store/cmst3/group/l1tr/vcamagni/L1TauID/DATA/FPinputs/m40/4STEPS/142Xv0/inputs140X_7099346_0.root'
        '/store/cmst3/group/l1tr/FastPUPPI/15_1_X/fpinputs_140X/v1/caseC_m220_67/4STEPS/151Xv0/inputs151X_14682953_1699.root'
    ]
)

process.runPF = cms.Path( 
        process.l1tGTTInputProducer +
        process.l1tTrackSelectionProducer +
        process.l1tVertexFinderEmulator +
        process.l1tLayer1BarrelExtended + 
        process.l1tLayer1HGCalExtended +
        process.l1tLayer1HGCalNoTK +
        process.l1tLayer1HF +
        process.l1tLayer1Extended 
)
process.runPF.associate(process.L1TLayer1TaskInputsTask)
process.runPF.associate(process.L1TInputTask)

process.packPuppi = cms.EDProducer("ScPhase2PuppiPacker",
    src = cms.InputTag("l1tLayer1Extended:Puppi"),
    fedIDs = cms.vuint32(0),
    splitFactor = cms.uint32(1),
    scoutingHeader = cms.bool(True)
)

from L1TriggerScouting.Phase2.modules import (
    l1sc_L1TScPhase2PuppiRawToDigi_alpaka, l1sc_L1TScPhase2SCJets_alpaka
)
process.unpackPuppi = cms.EDProducer('ScPhase2PuppiRawToDigi',
    fedIDs = process.packPuppi.fedIDs,
    splitFactor = process.packPuppi.splitFactor,
    src = cms.InputTag('packPuppi'),
)

process.p = cms.Path(
    process.packPuppi +
    process.unpackPuppi
)

if options.dumpClusters:
    process.puppiTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
        name = cms.string("L1Puppi"),
        src = cms.InputTag("l1tLayer1Extended:Puppi"),
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
    process.p_dumpclusters = cms.Path(process.puppiTable)

if options.run in ("sc4Alpaka",):
    process.unpackPuppiAlpaka = l1sc_L1TScPhase2PuppiRawToDigi_alpaka(
        streams = process.packPuppi.fedIDs,
        splitFactor = process.packPuppi.splitFactor,
        src = cms.InputTag('packPuppi'),
        environment = cms.untracked.int32(options.environment),
    )

    sc4Alpaka_kwargs = {
        "src": cms.InputTag("unpackPuppiAlpaka"),
        "algo": cms.string(options.jetAlgo)
    }

    assert options.jetAlgo in [
        "SCGreedy", "SCNMS", "SCNMSWeighted", "SCNMSWeightedMultiIter", "LinkTree"
    ]

    if options.jetAlgo == "SCGreedy":
        sc4Alpaka_kwargs.update({
            "rParam": cms.double(options.jetR),
            "nJets": cms.uint32(options.njets),
        })
    elif options.jetAlgo == "SCNMS" or options.jetAlgo == "SCNMSWeighted":
        sc4Alpaka_kwargs.update({
            "RSeed": cms.double(options.RSeed),
            "RClu": cms.double(options.RClu),
        })
    elif options.jetAlgo == "SCNMSWeightedMultiIter":
        sc4Alpaka_kwargs.update({
            "RSeed": cms.double(options.RSeed),
            "RCen": cms.double(options.RCen),
            "RClu": cms.double(options.RClu),
            "alphaSeed": cms.double(options.alphaSeed),
            "minSeedPt": cms.double(options.minSeedPt),
            "nCentroidIters": cms.uint32(options.nCentroidIters),
        })
    elif options.jetAlgo == "LinkTree":
        sc4Alpaka_kwargs.update({
            "RLink": cms.double(options.RLink),
            "minSeedPt": cms.double(options.minSeedPt),
        })

    process.scPhase2SC4PFAlpaka = l1sc_L1TScPhase2SCJets_alpaka(
        **sc4Alpaka_kwargs
    )

    process.p_sc4Alpaka = cms.Path(
        process.unpackPuppiAlpaka +
        process.scPhase2SC4PFAlpaka
    )

    if options.dumpClusters:
        process.dumpJetsTable = cms.EDProducer("ClusterObjSoAToNanoaodFlatTable",
            srcClusters = cms.InputTag("scPhase2SC4PFAlpaka"),
            name = cms.string(process.scPhase2SC4PFAlpaka.algo.value() + "Jet"),
            doc = cms.string(""),
        )
        process.clusterJetIndexTable = cms.EDProducer("ClusterMapperSoAToNanoaodFlatTable",
            srcClusters = cms.InputTag("scPhase2SC4PFAlpaka"),
            name = process.puppiTable.name,
            clustering_name = process.dumpJetsTable.name,
            doc = cms.string(""),
        )
        process.p_dump = cms.Path(process.dumpJetsTable + process.clusterJetIndexTable)

if options.dumpClusters:
    process.outStandard = cms.OutputModule("NanoAODOutputModule",
        fileName = cms.untracked.string("plainNano.root"),
        SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring()),
        outputCommands = cms.untracked.vstring("drop *", "keep nanoaodFlatTable_*Table_*_*"),
        compressionLevel = cms.untracked.int32(4),
        compressionAlgorithm = cms.untracked.string("ZLIB"),
    )
    process.e = cms.EndPath(process.outStandard)
