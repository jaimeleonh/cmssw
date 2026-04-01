import FWCore.ParameterSet.Config as cms

dimuStruct = cms.EDProducer("ScPhase2TrackerMuonDiMuDemo",
    src = cms.InputTag("scPhase2TrackerMuonRawToDigiStruct"),
)

zdeeStruct = cms.EDProducer("ScPhase2TkEmDarkPhotonDiEle",
    src = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    ptMin = cms.vdouble(5.0, 4.0),
    idScore = cms.vdouble(0.5, 0.5),  # in GT format, i.e. range is [0, 1] and not [-1, 1]
    etaMax = cms.double(1.479),
    dzMax = cms.double(1.0),
    relIso = cms.double(0.4),
)

