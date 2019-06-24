import FWCore.ParameterSet.Config as cms

doAging = True

process = cms.Process("L1DTTrigPhase2Prod")

process.load("Geometry.CMSCommonData.cmsIdealGeometryXML_cff")
process.load("Geometry.DTGeometry.dtGeometry_cfi")
process.DTGeometryESModule.applyAlignment = False

process.load("L1Trigger.DTPhase2Trigger.dtTriggerPhase2PrimitiveDigis_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = "80X_dataRun2_2016SeptRepro_v7"

if doAging:
  process.GlobalTag = GlobalTag(process.GlobalTag, '100X_upgrade2023_realistic_Candidate_2018_07_26_10_13_24', 'MuonSystemAging_3000fbm1,MuonSystemAgingRcd,sqlite_file:./MuonSystemAging.db')
  process.load("Phase2L1Trigger.CalibratedDigis.CalibratedDigis_cfi")

process.load("L1Trigger.DTPhase2Trigger.dtTriggerPhase2PrimitiveDigis_cfi")

process.source = cms.Source("PoolSource",fileNames = cms.untracked.vstring('file:/pool/ciencias/userstorage/folguera/P2L1T/digis_segments_Run2016BSingleMuonRAW-RECO_camilo.root'))
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(1000))

process.out = cms.OutputModule("PoolOutputModule",
       outputCommands = cms.untracked.vstring('keep *'),
       fileName = cms.untracked.string('DTTriggerPhase2Primitives.root')
)


from Phase2L1Trigger.CalibratedDigis.CalibratedDigis_cfi import CalibratedDigis as _CalibratedDigis

if doAging:
  from SimMuon.DTDigitizer.dtChamberMasker_cfi import dtChamberMasker as _dtChamberMasker
  print("[appendDTChamberMasker] : Found CalibratedDigis, appending producer for aged DTs and Our TriggerPrimitives producer dtPhase2Emulator\n")
  process.simAgedDtTriggerDigis = _dtChamberMasker.clone()
  process.simAgedDtTriggerDigis.digiTag = "muonDTDigis"
  process.simCalibratedDigis = _CalibratedDigis.clone()
  process.simCalibratedDigis.dtDigiTag = "simAgedDtTriggerDigis"
  process.withAgedDtTriggerSequence = cms.Sequence(process.simAgedDtTriggerDigis + process.simCalibratedDigis)
  process.RandomNumberGeneratorService = cms.Service("RandomNumberGeneratorService",simAgedDtTriggerDigis = cms.PSet(initialSeed = cms.untracked.uint32(789342)))    

else:
  process.simCalibratedDigis = _CalibratedDigis.clone()
  process.withAgedDtTriggerSequence = cms.Sequence(process.simCalibratedDigis)

process.dtTriggerPhase2PrimitiveDigis.digiTag = "simCalibratedDigis"

process.p = cms.Path(process.withAgedDtTriggerSequence*process.dtTriggerPhase2PrimitiveDigis)
process.this_is_the_end = cms.EndPath(process.out)
