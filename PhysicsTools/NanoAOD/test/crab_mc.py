from CRABClient.UserUtilities import config
config = config()

config.General.requestName = 'nano6'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'NANO_NANO.py'
config.JobType.allowUndistributedCMSSW = True
#config.JobType.outputFiles = ['lzma.root']

#config.Data.inputDataset = '/GluGluToHHTo2B2Tau_node_cHHH1_TuneCP5_PSWeights_13TeV-powheg-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
#config.Data.inputDataset = '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'
#config.Data.inputDataset = '/GluGluToHHTo2B2Tau_node_SM_TuneCP5_PSWeights_13TeV-madgraph-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM'
config.Data.inputDataset = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'
config.Data.inputDBS = 'global'
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
#config.Data.totalUnits = 1000000
config.Data.outLFNDirBase = '/store/user/jleonhol/NanoTest/' 
config.Data.publication = False
config.Data.outputDatasetTag = 'NanoTestFullDis'



config.Site.storageSite = 'T2_ES_CIEMAT'
