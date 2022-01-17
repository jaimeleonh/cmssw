import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_cff import Run3
import os
import re

options = VarParsing('analysis')
options.register('sampleType', 'Run3MC', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Indicates the sample type: Run3MC, Run2MC or Run2Data")
options.register('fileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of root files to process.")
options.register('fileNamePrefix', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "Prefix to add to input file names.")
options.register('lumiFile', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "JSON file with lumi mask.")
options.register('eventList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string,
                 "List of events to process.")
options.register('dumpPython', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "Dump full config into stdout.")
options.register('Summary', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool,
                 "want summary")

options.parseArguments()
dataDict = {"Run3MC" : False, "Run2Data" : True}
isData = dataDict[options.sampleType]
isMC = not(isData)
print(isMC)
processName = 'MYHLT'
process = cms.Process(processName)
if(options.sampleType=='Run3MC'):
    process = cms.Process(processName, Run3)
# import of standard configurations

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('HLTrigger.Configuration.HLT_diTauJet_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

if isMC:
    process.load('SimGeneral.MixingModule.mixNoPU_cfi')

process.source = cms.Source('PoolSource', fileNames = cms.untracked.vstring())

# Input source
def readFileList(fileList, inputFileName, fileNamePrefix):
    with open(inputFileName, 'r') as inputFile:
        for name in inputFile.readlines():
            if len(name) > 0 and name[0] != '#':
                fileList.append(fileNamePrefix + name)

def addFilesToList(fileList, inputFiles, fileNamePrefix):
    """read intput file list from a another list"""
    for name in inputFiles:
        if len(name) > 0 and name[0] != '#':
            fileList.append(fileNamePrefix + name)

if len(options.fileList) > 0:
    readFileList(process.source.fileNames, options.fileList, options.fileNamePrefix)
elif len(options.inputFiles) > 0:
    addFilesToList(process.source.fileNames, options.inputFiles, options.fileNamePrefix)
else:
    #print(options.inputFiles)
    #print(options.fileList)

    #process.source.fileNames = cms.untracked.vstring('file:/store/mc/Run3Winter20DRPremixMiniAOD/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v1/20000/1E71AB01-1FF5-F049-94D4-642325AFF937.root')  # 5000 evts

    #process.source.fileNames = cms.untracked.vstring('/store/mc/Run3Winter20DRPremixMiniAOD/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v1/20000/A95D23CB-88AC-6849-8980-5D12C8417007.root') # 11000 evts
    # process.source.fileNames = cms.untracked.vstring('file:/eos/home-v/vdamante/A95D23CB-88AC-6849-8980-5D12C8417007.root') # 11000 evts
    # process.source.fileNames = cms.untracked.vstring('/store/mc/Run3Winter21DRMiniAOD/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/270000/005b56c1-0107-46b3-9740-1c6efc559295.root') #1000 evts
    #process.source.fileNames = cms.untracked.vstring('file:/store/mc/Run3Winter20DRPremixMiniAOD/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v1/20000/F2760E46-A3DB-DA4F-A6EC-525C10EDCBC7.root') # 1000 evts
    #process.source.fileNames = cms.untracked.vstring('file:/store/mc/Run3Winter20DRPremixMiniAOD/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v1/20000/657F58E7-64E3-AA4D-B505-6D0F39997487.root')
    # process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_tau/amascell/005b56c1-0107-46b3-9740-1c6efc559295.root')
    # process.source.fileNames = cms.untracked.vstring('file:/eos/cms/store/group/phys_tau/TauTrigger/store/mc/Run3Winter20DRPremixMiniAOD/VBFHToTauTau_M125_TuneCUETP8M1_14TeV_powheg_pythia8/GEN-SIM-RAW/110X_mcRun3_2021_realistic_v6-v1/20000/28DC45D2-CB46-6D48-A100-379B076BFE1D.root')
    process.source.fileNames=cms.untracked.vstring('/store/mc/Run3Winter21DRMiniAOD/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/270000/005b56c1-0107-46b3-9740-1c6efc559295.root')
    # process.source.fileNames = cms.untracked.vstring('/store/mc/Run3Winter21DRMiniAOD/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/GEN-SIM-DIGI-RAW/FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/270000/434ff3bc-f045-4662-a793-2206c25de018.root')


if len(options.lumiFile) > 0:
    import FWCore.PythonUtilities.LumiList as LumiList
    process.source.lumisToProcess = LumiList.LumiList(filename = options.lumiFile).getVLuminosityBlockRange()

if options.eventList != '':
    process.source.eventsToProcess = cms.untracked.VEventRange(re.split(',', options.eventList))

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents),
    output = cms.optional.untracked.allowed(cms.int32,cms.PSet)
)


process.options = cms.untracked.PSet(
    FailPath = cms.untracked.vstring(),
    IgnoreCompletely = cms.untracked.vstring(),
    Rethrow = cms.untracked.vstring(),
    #SkipEvent = cms.untracked.vstring('ProductNotFound'),
    SkipEvent = cms.untracked.vstring(),
    allowUnscheduled = cms.obsolete.untracked.bool,
    canDeleteEarly = cms.untracked.vstring(),
    emptyRunLumiMode = cms.obsolete.untracked.string,
    eventSetup = cms.untracked.PSet(
        forceNumberOfConcurrentIOVs = cms.untracked.PSet(
            allowAnyLabel_=cms.required.untracked.uint32
        ),
        numberOfConcurrentIOVs = cms.untracked.uint32(1)
    ),
    fileMode = cms.untracked.string('FULLMERGE'),
    forceEventSetupCacheClearOnNewRun = cms.untracked.bool(False),
    makeTriggerResults = cms.obsolete.untracked.bool,
    numberOfConcurrentLuminosityBlocks = cms.untracked.uint32(1),
    numberOfConcurrentRuns = cms.untracked.uint32(1),
    numberOfStreams = cms.untracked.uint32(0),
    numberOfThreads = cms.untracked.uint32(1),
    printDependencies = cms.untracked.bool(False),
    sizeOfStackForThreadsInKB = cms.optional.untracked.uint32,
    throwIfIllegalParameter = cms.untracked.bool(True),
    wantSummary = cms.untracked.bool(False)
)


# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('tau_hlt nevts:-1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag import GlobalTag

GlobalTagName = 'auto:run3_data_GRun' if isData else 'auto:run3_mc_GRun'

print(GlobalTagName)
process.GlobalTag = GlobalTag(process.GlobalTag, GlobalTagName, '')

# Path and EndPath definitions
#process.L1simulation_step = cms.Path(cms.Sequence(cms.Task(process.SimL1TGlobalTask, process.simBscDigis, process.simCastorTechTrigDigis, process.simCscTriggerPrimitiveDigis, process.simCsctfDigis, process.simCsctfTrackDigis, process.simDtTriggerPrimitiveDigis, process.simDttfDigis, process.simGctDigis, process.simGmtDigis, process.simHcalTechTrigDigis, process.simRctDigis, process.simRpcTechTrigDigis, process.simRpcTriggerDigis)))
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule()
#process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step])
from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)


# customisation of the process.
process.hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60 = cms.EDFilter( "HLT3DoublePFTauPFJet",
    saveTags = cms.bool( True ),
    originTag1 = cms.VInputTag( 'hltHpsSelectedPFTausMediumDitauWPDeepTauDoubleTauJet' ),
    originTag2 = cms.VInputTag( 'hltHpsSelectedPFTausMediumDitauWPDeepTauDoubleTauJet' ),
    originTag3 = cms.VInputTag( 'hltPFJetsCorrectedMatchedToL1' ),
    inputTag1 = cms.InputTag( 'hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet' ),
    inputTag2 = cms.InputTag( 'hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet' ),
    inputTag3 = cms.InputTag( 'hltPFJets60L1HLTMatched' ),
    triggerType1 = cms.int32( 84 ),
    triggerType2 = cms.int32( 84 ),
    triggerType3 = cms.int32( 85 ),
    MinDphi = cms.double( 0.0 ),
    MaxDphi = cms.double( -1.0 ),
    MinDeta = cms.double( 1.0 ),
    MaxDeta = cms.double( -1.0 ),
    MinMinv = cms.double( 0.0 ),
    MaxMinv = cms.double( -1.0 ),
    MinDelR = cms.double( 0.5 ),
    MaxDelR = cms.double( 99999.0 ),
    MinPt = cms.double( 1.0 ),
    MaxPt = cms.double( -1.0 ),
    MinN = cms.int32( 1 )
)

# Automatic addition of the customisation function

if isData:
    # Customisation from command line
    from HLTrigger.Configuration.customizeHLTforCMSSW import customiseFor2018Input
    process = customiseFor2018Input(process)
else:
    from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC
    process = customizeHLTforMC(process)

# # customise tracking for Run3
# from HLTrigger.Configuration.customizeHLTforRun3Tracking import customizeHLTforRun3Tracking
# process = customizeHLTforRun3Tracking(process)
#from HLTrigger.Configuration.customizeHLTforPatatrack import customizeHLTforPatatrackTriplets
#process = customizeHLTforPatatrackTriplets(process)
#process.HLTL2TauTagNNSequenceDoubleTauJet._tasks.add(process.HLTDoLocalPixelTask)
#process.HLTL2TauTagNNSequenceDoubleTauJet._tasks.add(process.HLTRecoPixelTracksTask)
#process.HLTL2TauTagNNSequenceDoubleTauJet._tasks.add(process.HLTRecopixelvertexingTask)

# from TauMLTools.Studies.diTauJet_Full_optimization.createTauNNTagTuple_diTauJet_v2 import *
# process = updateDiTauJets(process)
# process = add_counters(process, options.outputFile, isMC)

# #from applyL2TauTag import update
# from TauTriggerTools.HLTProducers.applyL2TauTag import update as update_L2
# process = update_L2(process)

# from TauTriggerTools.HLTProducers.deepTauAtHLT import update as update_deepTau
# process = update_deepTau(process, useReg=False, resetWP=True, addCounters=True, isData=isData)

#from TauTriggerTools.HLTProducers.deepTauAtHLT import add_counters
#process = add_counters(process, isData=isData, resetWP=True)

# process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseChargedIsoPFTauHPS30_eta2p1_CrossL1_v1, process.HLT_IsoMu20_eta2p1_LooseChargedIsoPFTauHPS27_eta2p1_CrossL1_v4, process.HLT_DoubleMediumChargedIsoPFTauHPS35_Trk1_eta2p1_Reg_v4, process.HLT_MediumChargedIsoPFTau50_Trk30_eta2p1_1pr_MET100_v12, process.HLT_MediumChargedIsoPFTau180HighPtRelaxedIso_Trk50_eta2p1_v12, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

updatedTauName = "slimmedTausNewID" #name of pat::Tau collection with new tau-Ids
import RecoTauTag.RecoTau.tools.runTauIdMVA as tauIdConfig
tauIdEmbedder = tauIdConfig.TauIDEmbedder(process, cms, #debug = False,
                    updatedTauName = updatedTauName,
                    toKeep = ["deepTau2017v2p1", #deepTau TauIDs
                               ])
tauIdEmbedder.runTauID()

process.demo = cms.EDAnalyzer('NtupleMaker'
     #, tracks = cms.untracked.InputTag('ctfWithMaterialTracks')
     , triggerResults = cms.untracked.InputTag("TriggerResults", "", "MYHLT") #change to HLT instead of MYHLT if not rerun
     , triggerEvent = cms.untracked.InputTag("hltTriggerSummaryAOD", "", "MYHLT")
     , triggerEventWithRefs = cms.untracked.InputTag("hltTriggerSummaryRAW", "", "MYHLT")
     #, PFTauCollection = cms.untracked.InputTag("hltHpsPFTauProducer","","MYHLT")
     #, hltHpsTracks = cms.untracked.InputTag("hltHpsPFTauTrack", "","MYHLT")

     , SkipEvent = cms.untracked.vstring('ProductNotFound')

     , fillingTriggers = cms.untracked.bool(True)
     , fillingEventInfo = cms.untracked.bool(True)
     , fillingL1 = cms.untracked.bool(True)    # False for data
     , fillingTaus = cms.untracked.bool(True)  # False for data
     , fillingJets = cms.untracked.bool(True)  # False for data
     , development = cms.untracked.bool(False)
     , doGenParticles = cms.untracked.bool(False)

     , jetTriggerPrimitives = cms.untracked.InputTag("caloStage2Digis", "Jet", "RECO")
     , tauTriggerPrimitives = cms.untracked.InputTag("caloStage2Digis", "Tau", "RECO")

     , genParticleSrc = cms.untracked.InputTag("prunedGenParticles")
     , VtxLabel = cms.untracked.InputTag("offlineSlimmedPrimaryVertices")
     , rhoLabel = cms.untracked.InputTag("fixedGridRhoFastjetAll")
     , tauSrc = cms.untracked.InputTag("slimmedTausNewID")
     , ak4JetSrc = cms.untracked.InputTag("slimmedJets")
)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("histoAOD.root"),
    closeFileFast = cms.untracked.bool(True)
)



process.demoPath = cms.EndPath(
        process.rerunMvaIsolationSequence *  # Comment for data
        getattr(process,updatedTauName) *    # Comment for data
        process.demo
)


process.schedule = cms.Schedule( [*(process.HLTriggerFirstPath,process.HLT_DoubleMediumDeepTauIsoPFTauHPS30_L2NN_eta2p1_PFJet60_v2, process.HLTAnalyzerEndpath, process.HLTriggerFinalPath, process.endjob_step ), process.demoPath], tasks=[process.patAlgosToolsTask])

#process.schedule = cms.Schedule(*[ process.HLTriggerFirstPath, process.HLT_DoubleMediumDeepTauIsoPFTauHPS35_L2NN_eta2p1_v1, process.HLT_Ele24_eta2p1_WPTight_Gsf_LooseDeepTauPFTauHPS30_eta2p1_CrossL1_v1, process.HLT_IsoMu20_eta2p1_LooseDeepTauPFTauHPS27_eta2p1_CrossL1_v1, process.HLT_LooseDeepTauPFTauHPS180_L2NN_eta2p1_v1, process.HLT_LooseDeepTauPFTauHPS50_L2NN_eta2p1_MET100_v1, process.HLTriggerFinalPath, process.endjob_step ], tasks=[process.patAlgosToolsTask])

# End of customisation functions

# Customisation from command line

process.load('FWCore.MessageLogger.MessageLogger_cfi')
x = process.maxEvents.input.value()
x = x if x >= 0 else 10000
process.MessageLogger.cerr.FwkReport.reportEvery = max(1, min(1000, x // 10))

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)

if options.dumpPython:
    print (process.dumpPython())

process.options.wantSummary = cms.untracked.bool(options.Summary)
