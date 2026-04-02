import FWCore.ParameterSet.Config as cms


hjpsigammaRecMesonStruct = cms.EDProducer("ScPhase2BosonToRecMesonGamma",
    srcMeson = cms.InputTag("puppiRecMesonStruct", "jpsi"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
    minMassBoson = cms.double(60), # relax from 100 to include the Z
    maxMassBoson = cms.double(200),
    minPtQ = cms.double(15), # relax cuts from 30 to get more events in the scouting stream for now
    minPtGamma = cms.double(20), # ditto
    maxRelIsoQ = cms.double(0.25),
    analysisName = cms.string("HJPsiGamma")
)
hjpsigammaTTrackRecMesonStruct = hjpsigammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:jpsi",
    analysisName = "HJPsiGammaTTrack",
)
hjpsigammaMuMuRecMesonStruct = hjpsigammaRecMesonStruct.clone(
    srcMeson = "puppiRecMesonStruct:jpsi",
    analysisName = "HJPsiGammaMuMu",
)
hjpsigammaEERecMesonStruct = hjpsigammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:jpsi",
    analysisName = "HJPsiGammaEE",
)

h2phiRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "phi"),
    minMassBoson = cms.double(60), # relax from 100 to include the Z
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(10),
    maxIso = cms.double(0.25),
    analysisName = cms.string("H2Phi")
)
h2phiTTrackRecMesonStruct = h2phiRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:phi",
    srcMeson2 = "ttrackRecMesonStruct:phi",
    analysisName = "H2PhiTTrack",
)   

hphigammaRecMesonStruct = cms.EDProducer("ScPhase2BosonToRecMesonGamma",
    srcMeson = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
    minMassBoson = cms.double(60), # relax from 100 to include the Z
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    minPtGamma = cms.double(30),
    maxRelIsoQ = cms.double(0.25),
    analysisName = cms.string("HPhiGamma")
)
hphigammaTTrackRecMesonStruct = hphigammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:phi",
    analysisName = "HPhiGammaTTrack",
)


hphijpsiRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "phi"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "jpsi"),
    sameDaughersCollection = cms.bool(True),
    minMassBoson = cms.double(60), # relax from 100 to allow the Z
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    maxIso = cms.double(0.25),
    analysisName = cms.string("HPhiJPsi")
)
hphijpsiTTrackRecMesonStruct = hphijpsiRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:phi",
    srcMeson2 = "ttrackRecMesonStruct:jpsi",
    analysisName = "HPhiJPsiTTrack",
)
hphijpsiMuMuRecMesonStruct = hphijpsiRecMesonStruct.clone(
    srcMeson1 = "puppiRecMesonStruct:phi",
    srcMeson2 = "tkMuonRecMesonStruct:jpsi",
    sameDaughersCollection = False,
    analysisName = "HPhiJPsiMuMu",
)
hphijpsiEERecMesonStruct = hphijpsiRecMesonStruct.clone(
    srcMeson1 = "puppiRecMesonStruct:phi",
    srcMeson2 = "tkEleRecMesonStruct:jpsi",
    sameDaughersCollection = False,
    analysisName = "HPhiJPsiEE",
)


h2rhoRecMesonStruct = cms.EDProducer("ScPhase2BosonTo2RecMeson",
    srcMeson1 = cms.InputTag("puppiRecMesonStruct", "rho"),
    srcMeson2 = cms.InputTag("puppiRecMesonStruct", "rho"),
    minMassBoson = cms.double(60), # relax from 100 to allow the Z
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(3),
    maxIso = cms.double(0.25),
    analysisName = cms.string("H2Rho")
)
h2rhoTTrackRecMesonStruct = h2rhoRecMesonStruct.clone(
    srcMeson1 = "ttrackRecMesonStruct:rho",
    srcMeson2 = "ttrackRecMesonStruct:rho",
    analysisName = "H2RhoTTrack",
)

hrhogammaRecMesonStruct = cms.EDProducer("ScPhase2BosonToRecMesonGamma",
    srcMeson = cms.InputTag("puppiRecMesonStruct", "rho"),
    srcGamma = cms.InputTag("recIsoTkEmStruct"),
    minMassBoson = cms.double(60), # relax from 100 to allow the Z
    maxMassBoson = cms.double(150),
    minPtQ = cms.double(30),
    minPtGamma = cms.double(30),
    maxRelIsoQ = cms.double(0.25),
    analysisName = cms.string("HRhoGamma")
)
hrhogammaTTrackRecMesonStruct = hrhogammaRecMesonStruct.clone(
    srcMeson = "ttrackRecMesonStruct:rho",
    analysisName = "HRhoGammaTTrack",
)

w3piStruct = cms.EDProducer("ScPhase2PuppiW3PiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)

wdsgStruct = cms.EDProducer("ScPhase2PuppiWDsGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    ptHad = cms.vdouble(10., 4., 3.),
    ptTkEm = cms.double(25.),
    maxDeltaRHad = cms.double(0.15),
    minMassDs = cms.double(1.2), # tight is 1.75
    maxMassDs = cms.double(3.0), # tight is 2.3
    minMass = cms.double(60.),
    maxMass = cms.double(100.),
    maxDeltaRDsTkEm = cms.double(3.5),
    minDeltaPhiDsTkEm = cms.double(2.5),
    relIsoDs = cms.double(0.5), # tight cut is 0.45
    relIsoTkEm = cms.double(0.5), # tight cut is 0.25
    isolationMinDeltaRDs = cms.double(0.00),
    isolationMaxDeltaRDs = cms.double(0.50),
    isolationMinDeltaRTkEm = cms.double(0.02),
    isolationMaxDeltaRTkEm = cms.double(0.50),
)


wpigStruct = cms.EDProducer("ScPhase2PuppiWPiGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
    ptPi = cms.double(20), # tight cut is 25
    ptTkEm = cms.double(20),
    minMass = cms.double(60.),
    maxMass = cms.double(100.),
    minDeltaRPiTkEm = cms.double(0.5),
    relIsoPi = cms.double(0.50), # tight cut is 0.3
    relIsoTkEm = cms.double(0.50), # tight cut is 0.3
    isolationMinDeltaRPi = cms.double(0.00),
    isolationMaxDeltaRPi = cms.double(0.50),
    isolationMinDeltaRTkEm = cms.double(0.02),
    isolationMaxDeltaRTkEm = cms.double(0.50)
)

hrhogStruct = cms.EDProducer("ScPhase2PuppiHRhoGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

hphigStruct = cms.EDProducer("ScPhase2PuppiHPhiGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

hjpsigStruct = cms.EDProducer("ScPhase2PuppiHJPsiGammaDemo",
    srcPuppi = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
    srcTkEm = cms.InputTag("scPhase2TkEmRawToDigiStruct"),
)

h2rhoStruct = cms.EDProducer("ScPhase2PuppiH2RhoDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)

h2phiStruct = cms.EDProducer("ScPhase2PuppiH2PhiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)

hphijpsiStruct = cms.EDProducer("ScPhase2PuppiHPhiJPsiDemo",
    src = cms.InputTag("scPhase2PuppiRawToDigiStruct"),
)