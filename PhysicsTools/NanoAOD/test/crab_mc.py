from CRABClient.UserUtilities import config
config = config()

#config.General.requestName = 'nano6'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True


config.JobType.maxMemoryMB = 3000
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'NANO_NANO.py'
#config.JobType.allowUndistributedCMSSW = True
#config.JobType.outputFiles = ['lzma.root']

#config.Data.inputDataset = '/GluGluToHHTo2B2Tau_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
# config.Data.inputDataset = '/Muon/Run2022C-PromptReco-v1/MINIAOD'
#config.Data.inputDataset = '/SingleMuon/Run2022C-PromptReco-v1/MINIAOD'
config.Data.inputDataset = '/Muon/Run2022F-PromptReco-v1/MINIAOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.unitsPerJob = 100
#config.Data.totalUnits = 1000000
config.Data.outLFNDirBase = '/store/user/jleonhol/L1_Run2022F/'
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoTestFull'
config.Data.lumiMask = '/afs/cern.ch/work/j/jleonhol/private/L1/nanoaod/CMSSW_12_4_8/src/Cert_Collisions2022_355100_362167_13p6TeV_DCSOnly_TkPx.json'
#config.Data.lumiMask = '/afs/cern.ch/work/j/jleonhol/private/L1/nanoaod/CMSSW_12_4_8/src/Cert_Collisions2022_eraF_360390_361580_Golden.json'
# config.Data.runRange = '357078'
# config.Data.runRange = '357438,357440,357479,357542,357610,357611,357612,357688,357696,357697,357769,357770,357771,357802'
# config.Data.runRange = '356077,356076,356075,356071,356005,356003,355933,355921,355872,355680'
config.Data.runRange = '362153,362154'


config.Site.storageSite = 'T3_CH_CERNBOX'
