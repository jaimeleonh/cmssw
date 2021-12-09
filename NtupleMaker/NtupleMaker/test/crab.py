#Kyungwook's example crab file

import sys

from CRABClient.UserUtilities import config
config = config()

config.General.workArea = 'test'
config.General.transferOutputs = True
config.General.transferLogs = True

# config.JobType.numCores = 4


config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'hlt_full.py'#'AODConfFile.py'#'zeroBiasHLT.py'#'HLTandRAW4520Taus.py'#'ditau_and_vbf.py'
# config.JobType.maxJobRuntimeMin = 120

config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 1
config.Data.outLFNDirBase = '/store/user/jleonhol/hltNtuples/' #'/store/user/knam'

# config.Data.ignoreLocality = True
#config.Site.whitelist = ['T2_US_*', 'T2_RU_JINR', 'T1_RU_JINR']
#config.Site.whitelist = ['T3_US_FNALLPC']
#config.Site.whitelist = ['T2_US_*']
# config.Site.whitelist = ['T1_US_FNAL', 'T2_FR_GRIF_LLR', 'T2_HU_Budapest']
#config.Site.ignoreGlobalBlacklist = True
#config.JobType.maxJobRuntimeMin = 2000
# config.JobType.maxMemoryMB = 4000
#config.JobType.numCores = 4
#config.JobType.inputFiles = ['L1Menu_Collisions2018_v1_0_0-d1_fixed.xml']
config.JobType.inputFiles = ['L1Menu_Collisions2022_v0_1_1_modified.xml']

config.Data.inputDataset = '/GluGluToHHTo2B2Tau_node_cHHH1_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'
#config.Data.inputDataset = '/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_TuneCP5_14TeV-madgraph-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/MINIAODSIM'
#config.Data.inputDataset = '/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/MINIAODSIM'
#config.Data.inputDataset = '/GluGluHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/MINIAODSIM'
config.Data.secondaryInputDataset = '/GluGluToHHTo2B2Tau_node_cHHH1_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
#config.Data.secondaryInputDataset = '/VBFHHTo2B2Tau_CV_1_C2V_1_C3_1_TuneCP5_14TeV-madgraph-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v2/GEN-SIM-DIGI-RAW'
#config.Data.secondaryInputDataset = '/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/GEN-SIM-DIGI-RAW'
#config.Data.secondaryInputDataset = '/GluGluHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/GEN-SIM-DIGI-RAW'

#config.Data.runRange = '321755,323725,323755,323790,323841,323940,323976,323978,324021,324077,324201,324237,324245,324293,324315,324420,324747,324785,324835,324897,324970,324980,324997,325022,325057,325097-325099'

config.Site.storageSite = 'T2_ES_CIEMAT' #'T3_KR_KNU'

# if __name__ == '__main__':

    # from CRABAPI.RawCommand import crabCommand

    #for i in range(1,9): 
      #config.General.requestName = 'ZB_selectedRuns_HLT_CorrectedRmvOl_1_EphemeralZeroBias{}'.format(i)
      #config.Data.inputDataset = '/EphemeralZeroBias{}/Run2018D-v1/RAW'.format(i)
      #crabCommand('submit', config = config)
    
    # newest ggH
    # config.Data.inputDataset = '/GluGluHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/MINIAODSIM'
    # config.Data.secondaryInputDataset = '/GluGluHToTauTau_M-125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/GEN-SIM-DIGI-RAW'
    # newest VBF
    #config.Data.inputDataset = '/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/MINIAODSIM'
    #config.Data.secondaryInputDataset = '/VBFHToTauTau_M125_TuneCP5_14TeV-powheg-pythia8/Run3Winter21DRMiniAOD-FlatPU30to80FEVT_112X_mcRun3_2021_realistic_v16-v1/GEN-SIM-DIGI-RAW'

    # older samples
    #config.Data.inputDataset = '/GluGluHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM'
    #config.Data.inputDataset = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext1-v1/MINIAODSIM'
    #config.Data.inputDataset = '/EphemeralZeroBias8/Run2018D-PromptReco-v2/MINIAOD'
    #config.Data.secondaryInputDataset = '/EphemeralZeroBias8/Run2018D-v1/RAW' 
    #config.Data.inputDataset = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISpring18MiniAOD-NZSPU28to70_100X_upgrade2018_realistic_v10-v1/MINIAODSIM'
    #config.Data.secondaryInputDataset = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISpring18DR-NZSPU28to70_100X_upgrade2018_realistic_v10-v1/GEN-SIM-RAW'
    #config.Data.inputDataset = '/VBFHToTauTau_M125_13TeV_powheg_pythia8/RunIISpring18DR-NZSPU28to70_100X_upgrade2018_realistic_v10-v1/GEN-SIM-RAW'
 
    # crabCommand('submit', config = config)
