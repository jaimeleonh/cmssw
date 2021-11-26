// Updated by Abdollah Mohammadi (KSU)  [18 May 2015]
// abdollah.mohammadi@cern.ch

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "NtupleMaker/NtupleMaker/plugins/NtupleMaker.h"

using namespace std;

// (local) variables associated with tree branches
Int_t          nTau_;

int runNumber;
int lumiBlock;
unsigned long long eventNumberID;

// decay mode discriminators

// Tau L1 Primitives
vector<float> 	tauL1PrimitivesPt_;
vector<float> 	tauL1PrimitivesEta_;
vector<float> 	tauL1PrimitivesPhi_;
vector<float> 	tauL1PrimitivesEnergy_;
vector<int> 	tauL1PrimitivesIso_;

//Tau Id & Isolation
vector<bool>   taupfTausDiscriminationByDecayModeFinding_;
vector<bool>   taupfTausDiscriminationByDecayModeFindingNewDMs_;

vector<bool>   tauByMVA6VLooseElectronRejection_;
vector<bool>   tauByMVA6LooseElectronRejection_;
vector<bool>   tauByMVA6MediumElectronRejection_;
vector<bool>   tauByMVA6TightElectronRejection_;
vector<bool>   tauByMVA6VTightElectronRejection_;

vector<bool>   tauByLooseMuonRejection3_;
vector<bool>   tauByTightMuonRejection3_;

vector<bool>   tauByLooseCombinedIsolationDeltaBetaCorr3Hits_;
vector<bool>   tauByMediumCombinedIsolationDeltaBetaCorr3Hits_;
vector<bool>   tauByTightCombinedIsolationDeltaBetaCorr3Hits_;
vector<float>   tauCombinedIsolationDeltaBetaCorrRaw3Hits_;

vector<bool>    tauByVVVLooseDeepTau2017v2p1VSjet_;
vector<bool>    tauByVVLooseDeepTau2017v2p1VSjet_;
vector<bool>    tauByVLooseDeepTau2017v2p1VSjet_;
vector<bool>    tauByLooseDeepTau2017v2p1VSjet_;
vector<bool>    tauByMediumDeepTau2017v2p1VSjet_;
vector<bool>    tauByTightDeepTau2017v2p1VSjet_;
vector<bool>    tauByVTightDeepTau2017v2p1VSjet_;
vector<bool>    tauByVVTightDeepTau2017v2p1VSjet_;
vector<bool>    tauByVVVLooseDeepTau2017v2p1VSe_;
vector<bool>    tauByVVLooseDeepTau2017v2p1VSe_;
vector<bool>    tauByVLooseDeepTau2017v2p1VSe_;
vector<bool>    tauByLooseDeepTau2017v2p1VSe_;
vector<bool>    tauByMediumDeepTau2017v2p1VSe_;
vector<bool>    tauByTightDeepTau2017v2p1VSe_;
vector<bool>    tauByVTightDeepTau2017v2p1VSe_;
vector<bool>    tauByVVTightDeepTau2017v2p1VSe_;
//vector<bool>    tauByVVVLooseDeepTau2017v2p1VSmu_;
//vector<bool>    tauByVVLooseDeepTau2017v2p1VSmu_;
vector<bool>    tauByVLooseDeepTau2017v2p1VSmu_;
vector<bool>    tauByLooseDeepTau2017v2p1VSmu_;
vector<bool>    tauByMediumDeepTau2017v2p1VSmu_;
vector<bool>    tauByTightDeepTau2017v2p1VSmu_;
//vector<bool>    tauByVTightDeepTau2017v2p1VSmu_;
//vector<bool>    tauByVVTightDeepTau2017v2p1VSmu_;

vector<float>        tauByIsolationMVArun2v1DBnewDMwLTraw_;
vector<float>        tauByIsolationMVArun2v1DBoldDMwLTraw_;
vector<float>        tauByIsolationMVArun2v1PWnewDMwLTraw_;
vector<float>        tauByIsolationMVArun2v1PWoldDMwLTraw_;
vector<bool>        tauByVTightIsolationMVArun2v1DBnewDMwLT_;
vector<bool>        tauByVTightIsolationMVArun2v1DBoldDMwLT_;
vector<bool>        tauByVTightIsolationMVArun2v1PWnewDMwLT_;
vector<bool>        tauByVTightIsolationMVArun2v1PWoldDMwLT_;
vector<bool>        tauByTightIsolationMVArun2v1DBnewDMwLT_;
vector<bool>        tauByTightIsolationMVArun2v1DBoldDMwLT_;
vector<bool>        tauByTightIsolationMVArun2v1PWnewDMwLT_;
vector<bool>        tauByTightIsolationMVArun2v1PWoldDMwLT_;
vector<bool>        tauByMediumIsolationMVArun2v1DBnewDMwLT_;
vector<bool>        tauByMediumIsolationMVArun2v1DBoldDMwLT_;
vector<bool>        tauByMediumIsolationMVArun2v1PWnewDMwLT_;
vector<bool>        tauByMediumIsolationMVArun2v1PWoldDMwLT_;
vector<bool>        tauByLooseIsolationMVArun2v1DBnewDMwLT_;
vector<bool>        tauByLooseIsolationMVArun2v1DBoldDMwLT_;
vector<bool>        tauByLooseIsolationMVArun2v1PWnewDMwLT_;
vector<bool>        tauByLooseIsolationMVArun2v1PWoldDMwLT_;
vector<bool>        tauByVLooseIsolationMVArun2v1DBnewDMwLT_;
vector<bool>        tauByVLooseIsolationMVArun2v1DBoldDMwLT_;
vector<bool>        tauByVLooseIsolationMVArun2v1PWnewDMwLT_;
vector<bool>        tauByVLooseIsolationMVArun2v1PWoldDMwLT_;


//Tau Kinematics
vector<float> tauEta_;
vector<float> tauPhi_;
vector<float> tauPt_;
vector<float> tauEt_;
vector<float> tauCharge_;
vector<int>   tauDecayMode_;
vector<float> tauP_;
vector<float> tauPx_;
vector<float> tauPy_;
vector<float> tauPz_;
vector<float> tauVz_;
vector<float> tauEnergy_;
vector<float> tauMass_;
vector<float> tauDxy_;
vector<float> tauZImpact_;

//Tau Ingredients
vector<float> tauChargedIsoPtSum_;
vector<float> tauNeutralIsoPtSum_;
vector<float> tauPuCorrPtSum_;
vector<int> tauNumSignalPFChargedHadrCands_;
vector<int> tauNumSignalPFNeutrHadrCands_;
vector<int> tauNumSignalPFGammaCands_;
vector<int> tauNumSignalPFCands_;
vector<int> tauNumIsolationPFChargedHadrCands_;
vector<int> tauNumIsolationPFNeutrHadrCands_;
vector<int> tauNumIsolationPFGammaCands_;
vector<int> tauNumIsolationPFCands_;
vector<bool>  tauLeadChargedHadronExists_;
vector<float> tauLeadChargedHadronEta_;
vector<float> tauLeadChargedHadronPhi_;
vector<float> tauLeadChargedHadronPt_;
vector<float> tauneutralIsoPtSumWeight_;
vector<float> taufootprintCorrection_;
vector<float> tauphotonPtSumOutsideSignalCone_;
vector<float> taudz_;
vector<float> taudxy_;

void NtupleMaker::branchesEventInfo(TTree* tree)
{
    tree->Branch("runNumber", &runNumber);
    tree->Branch("lumiBlock", &lumiBlock);
    tree->Branch("eventNumberID", &eventNumberID);
}

void NtupleMaker::branchesL1Taus(TTree* tree) 
{
    tree->Branch("tauL1PrimitivesPt", &tauL1PrimitivesPt_);
    tree->Branch("tauL1PrimitivesEta", &tauL1PrimitivesEta_);
    tree->Branch("tauL1PrimitivesPhi", &tauL1PrimitivesPhi_);
    tree->Branch("tauL1PrimitivesEnergy", &tauL1PrimitivesEnergy_);
    tree->Branch("tauL1PrimitivesIso", &tauL1PrimitivesIso_);
}

void NtupleMaker::branchesTaus(TTree* tree)
{

//    tree->Branch("runNumber", &runNumber);
//    tree->Branch("lumiBlock", &lumiBlock);
//    tree->Branch("eventNumberID", &eventNumberID);

    tree->Branch("nTau", &nTau_);
    
    // Tau L1 Primitives
    //Tau Id & Isolation
    tree->Branch("taupfTausDiscriminationByDecayModeFinding", &taupfTausDiscriminationByDecayModeFinding_);
    tree->Branch("taupfTausDiscriminationByDecayModeFindingNewDMs", &taupfTausDiscriminationByDecayModeFindingNewDMs_);
    
    tree->Branch("tauByMVA6VLooseElectronRejection", &tauByMVA6VLooseElectronRejection_);
    tree->Branch("tauByMVA6LooseElectronRejection", &tauByMVA6LooseElectronRejection_);
    tree->Branch("tauByMVA6MediumElectronRejection", &tauByMVA6MediumElectronRejection_);
    tree->Branch("tauByMVA6TightElectronRejection", &tauByMVA6TightElectronRejection_);
    tree->Branch("tauByMVA6VTightElectronRejection", &tauByMVA6VTightElectronRejection_);
    
    tree->Branch("tauByLooseMuonRejection3", &tauByLooseMuonRejection3_);
    tree->Branch("tauByTightMuonRejection3", &tauByTightMuonRejection3_);
    
    tree->Branch("tauByLooseCombinedIsolationDeltaBetaCorr3Hits", &tauByLooseCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauByMediumCombinedIsolationDeltaBetaCorr3Hits", &tauByMediumCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauByTightCombinedIsolationDeltaBetaCorr3Hits", &tauByTightCombinedIsolationDeltaBetaCorr3Hits_);
    tree->Branch("tauCombinedIsolationDeltaBetaCorrRaw3Hits", &tauCombinedIsolationDeltaBetaCorrRaw3Hits_);
   
    tree->Branch("tauByVVVLooseDeepTau2017v2p1VSjet", &tauByVVVLooseDeepTau2017v2p1VSjet_);
    tree->Branch("tauByVVLooseDeepTau2017v2p1VSjet", &tauByVVLooseDeepTau2017v2p1VSjet_);
    tree->Branch("tauByVLooseDeepTau2017v2p1VSjet", &tauByVLooseDeepTau2017v2p1VSjet_);
    tree->Branch("tauByLooseDeepTau2017v2p1VSjet", &tauByLooseDeepTau2017v2p1VSjet_);
    tree->Branch("tauByMediumDeepTau2017v2p1VSjet", &tauByMediumDeepTau2017v2p1VSjet_);
    tree->Branch("tauByTightDeepTau2017v2p1VSjet", &tauByTightDeepTau2017v2p1VSjet_);
    tree->Branch("tauByVTightDeepTau2017v2p1VSjet", &tauByVTightDeepTau2017v2p1VSjet_);
    tree->Branch("tauByVVTightDeepTau2017v2p1VSjet", &tauByVVTightDeepTau2017v2p1VSjet_);
    tree->Branch("tauByVVVLooseDeepTau2017v2p1VSe", &tauByVVVLooseDeepTau2017v2p1VSe_);
    tree->Branch("tauByVVLooseDeepTau2017v2p1VSe", &tauByVVLooseDeepTau2017v2p1VSe_);
    tree->Branch("tauByVLooseDeepTau2017v2p1VSe", &tauByVLooseDeepTau2017v2p1VSe_);
    tree->Branch("tauByLooseDeepTau2017v2p1VSe", &tauByLooseDeepTau2017v2p1VSe_);
    tree->Branch("tauByMediumDeepTau2017v2p1VSe", &tauByMediumDeepTau2017v2p1VSe_);
    tree->Branch("tauByTightDeepTau2017v2p1VSe", &tauByTightDeepTau2017v2p1VSe_);
    tree->Branch("tauByVTightDeepTau2017v2p1VSe", &tauByVTightDeepTau2017v2p1VSe_);
    tree->Branch("tauByVVTightDeepTau2017v2p1VSe", &tauByVVTightDeepTau2017v2p1VSe_);
    //tree->Branch("tauByVVVLooseDeepTau2017v2p1VSmu", &tauByVVVLooseDeepTau2017v2p1VSmu_);
    //tree->Branch("tauByVVLooseDeepTau2017v2p1VSmu", &tauByVVLooseDeepTau2017v2p1VSmu_);
    tree->Branch("tauByVLooseDeepTau2017v2p1VSmu", &tauByVLooseDeepTau2017v2p1VSmu_);
    tree->Branch("tauByLooseDeepTau2017v2p1VSmu", &tauByLooseDeepTau2017v2p1VSmu_);
    tree->Branch("tauByMediumDeepTau2017v2p1VSmu", &tauByMediumDeepTau2017v2p1VSmu_);
    tree->Branch("tauByTightDeepTau2017v2p1VSmu", &tauByTightDeepTau2017v2p1VSmu_);
    //tree->Branch("tauByVTightDeepTau2017v2p1VSmu", &tauByVTightDeepTau2017v2p1VSmu_);
    //tree->Branch("tauByVVTightDeepTau2017v2p1VSmu", &tauByVVTightDeepTau2017v2p1VSmu_);
 
    tree->Branch("tauByIsolationMVArun2v1DBnewDMwLTraw", &tauByIsolationMVArun2v1DBnewDMwLTraw_);
    tree->Branch("tauByIsolationMVArun2v1DBoldDMwLTraw", &tauByIsolationMVArun2v1DBoldDMwLTraw_);
    tree->Branch("tauByIsolationMVArun2v1PWnewDMwLTraw", &tauByIsolationMVArun2v1PWnewDMwLTraw_);
    tree->Branch("tauByIsolationMVArun2v1PWoldDMwLTraw", &tauByIsolationMVArun2v1PWoldDMwLTraw_);
    tree->Branch("tauByVTightIsolationMVArun2v1DBnewDMwLT", &tauByVTightIsolationMVArun2v1DBnewDMwLT_);
    tree->Branch("tauByVTightIsolationMVArun2v1DBoldDMwLT", &tauByVTightIsolationMVArun2v1DBoldDMwLT_);
    tree->Branch("tauByVTightIsolationMVArun2v1PWnewDMwLT", &tauByVTightIsolationMVArun2v1PWnewDMwLT_);
    tree->Branch("tauByVTightIsolationMVArun2v1PWoldDMwLT", &tauByVTightIsolationMVArun2v1PWoldDMwLT_);
    tree->Branch("tauByTightIsolationMVArun2v1DBnewDMwLT", &tauByTightIsolationMVArun2v1DBnewDMwLT_);
    tree->Branch("tauByTightIsolationMVArun2v1DBoldDMwLT", &tauByTightIsolationMVArun2v1DBoldDMwLT_);
    tree->Branch("tauByTightIsolationMVArun2v1PWnewDMwLT", &tauByTightIsolationMVArun2v1PWnewDMwLT_);
    tree->Branch("tauByTightIsolationMVArun2v1PWoldDMwLT", &tauByTightIsolationMVArun2v1PWoldDMwLT_);
    tree->Branch("tauByMediumIsolationMVArun2v1DBnewDMwLT", &tauByMediumIsolationMVArun2v1DBnewDMwLT_);
    tree->Branch("tauByMediumIsolationMVArun2v1DBoldDMwLT", &tauByMediumIsolationMVArun2v1DBoldDMwLT_);
    tree->Branch("tauByMediumIsolationMVArun2v1PWnewDMwLT", &tauByMediumIsolationMVArun2v1PWnewDMwLT_);
    tree->Branch("tauByMediumIsolationMVArun2v1PWoldDMwLT", &tauByMediumIsolationMVArun2v1PWoldDMwLT_);
    tree->Branch("tauByLooseIsolationMVArun2v1DBnewDMwLT", &tauByLooseIsolationMVArun2v1DBnewDMwLT_);
    tree->Branch("tauByLooseIsolationMVArun2v1DBoldDMwLT", &tauByLooseIsolationMVArun2v1DBoldDMwLT_);
    tree->Branch("tauByLooseIsolationMVArun2v1PWnewDMwLT", &tauByLooseIsolationMVArun2v1PWnewDMwLT_);
    tree->Branch("tauByLooseIsolationMVArun2v1PWoldDMwLT", &tauByLooseIsolationMVArun2v1PWoldDMwLT_);
    tree->Branch("tauByVLooseIsolationMVArun2v1DBnewDMwLT", &tauByVLooseIsolationMVArun2v1DBnewDMwLT_);
    tree->Branch("tauByVLooseIsolationMVArun2v1DBoldDMwLT", &tauByVLooseIsolationMVArun2v1DBoldDMwLT_);
    tree->Branch("tauByVLooseIsolationMVArun2v1PWnewDMwLT", &tauByVLooseIsolationMVArun2v1PWnewDMwLT_);
    tree->Branch("tauByVLooseIsolationMVArun2v1PWoldDMwLT", &tauByVLooseIsolationMVArun2v1PWoldDMwLT_);
    
    //Tau Kinematics
    tree->Branch("tauEta"  ,&tauEta_);
    tree->Branch("tauPhi"  ,&tauPhi_);
    tree->Branch("tauPt"  ,&tauPt_);
    tree->Branch("tauEt"  ,&tauEt_);
    tree->Branch("tauCharge"  ,&tauCharge_);
    tree->Branch("tauP"  ,&tauP_);
    tree->Branch("tauPx"  ,&tauPx_);
    tree->Branch("tauPy"  ,&tauPy_);
    tree->Branch("tauPz"  ,&tauPz_);
    tree->Branch("tauVz"  ,&tauVz_);
    tree->Branch("tauEnergy"  ,&tauEnergy_);
    tree->Branch("tauMass"  ,&tauMass_);
    tree->Branch("tauDxy"  ,&tauDxy_);
    tree->Branch("tauZImpact"  ,&tauZImpact_);
    
    // Tau Ingredients
    tree->Branch("tauDecayMode"  ,&tauDecayMode_);
    tree->Branch("tauLeadChargedHadronExists"  ,&tauLeadChargedHadronExists_);
    tree->Branch("tauLeadChargedHadronEta"  ,&tauLeadChargedHadronEta_);
    tree->Branch("tauLeadChargedHadronPhi"  ,&tauLeadChargedHadronPhi_);
    tree->Branch("tauLeadChargedHadronPt"  ,&tauLeadChargedHadronPt_);
    tree->Branch("tauChargedIsoPtSum"  ,&tauChargedIsoPtSum_);
    tree->Branch("tauNeutralIsoPtSum"  ,&tauNeutralIsoPtSum_);
    tree->Branch("tauPuCorrPtSum"  ,&tauPuCorrPtSum_);
    tree->Branch("tauNumSignalPFChargedHadrCands"  ,&tauNumSignalPFChargedHadrCands_);
    tree->Branch("tauNumSignalPFNeutrHadrCands"  ,&tauNumSignalPFNeutrHadrCands_);
    tree->Branch("tauNumSignalPFGammaCands"  ,&tauNumSignalPFGammaCands_);
    tree->Branch("tauNumSignalPFCands"  ,&tauNumSignalPFCands_);
    tree->Branch("tauNumIsolationPFChargedHadrCands"  ,&tauNumIsolationPFChargedHadrCands_);
    tree->Branch("tauNumIsolationPFNeutrHadrCands"  ,&tauNumIsolationPFNeutrHadrCands_);
    tree->Branch("tauNumIsolationPFGammaCands"  ,&tauNumIsolationPFGammaCands_);
    tree->Branch("tauNumIsolationPFCands"  ,&tauNumIsolationPFCands_);
    
    
    tree->Branch("taufootprintCorrection"  ,&taufootprintCorrection_);
    tree->Branch("tauphotonPtSumOutsideSignalCone"  ,&tauphotonPtSumOutsideSignalCone_);
    tree->Branch("taudz"  ,&taudz_);
    tree->Branch("taudxy"  ,&taudxy_);
    
    
}

void NtupleMaker::fillEventInfo(const edm::Event& iEvent)
{
    runNumber = iEvent.id().run();
    lumiBlock = iEvent.id().luminosityBlock();
    eventNumberID = iEvent.id().event();
}

void NtupleMaker::fillL1Taus(const edm::Event& iEvent)
{
    tauL1PrimitivesPt_.clear();
    tauL1PrimitivesEta_.clear();
    tauL1PrimitivesPhi_.clear();
    tauL1PrimitivesEnergy_.clear();
    tauL1PrimitivesIso_.clear();

    edm::Handle<BXVector<l1t::Tau>> tauL1Handle;
    iEvent.getByToken(tauTriggerPrimitives_, tauL1Handle);

    for(BXVector<l1t::Tau>::const_iterator itau = tauL1Handle->begin(0); itau != tauL1Handle->end(0); ++itau) {
        tauL1PrimitivesPt_.push_back(itau->pt());
        tauL1PrimitivesEta_.push_back(itau->eta());
        tauL1PrimitivesPhi_.push_back(itau->phi());
        tauL1PrimitivesEnergy_.push_back(itau->energy());
        tauL1PrimitivesIso_.push_back(itau->hwIso());
    }

}

void NtupleMaker::fillTaus(const edm::Event& e)
{
/***
    // Tau L1 primitives
    tauL1PrimitivesPt_.clear(); 
    tauL1PrimitivesEta_.clear(); 
    tauL1PrimitivesPhi_.clear(); 
    tauL1PrimitivesEnergy_.clear(); 
    tauL1PrimitivesIso_.clear(); 
   ***/
    // Tau Id & Isolation
    
    taupfTausDiscriminationByDecayModeFinding_.clear();
    taupfTausDiscriminationByDecayModeFindingNewDMs_.clear();
    
    tauByMVA6VLooseElectronRejection_.clear();
    tauByMVA6LooseElectronRejection_.clear();
    tauByMVA6MediumElectronRejection_.clear();
    tauByMVA6TightElectronRejection_.clear();
    tauByMVA6VTightElectronRejection_.clear();
    
    tauByLooseMuonRejection3_.clear();
    tauByTightMuonRejection3_.clear();

    tauByLooseCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauByMediumCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauByTightCombinedIsolationDeltaBetaCorr3Hits_.clear();
    tauCombinedIsolationDeltaBetaCorrRaw3Hits_.clear();

    tauByVVVLooseDeepTau2017v2p1VSjet_.clear();
    tauByVVLooseDeepTau2017v2p1VSjet_.clear();
    tauByVLooseDeepTau2017v2p1VSjet_.clear();
    tauByLooseDeepTau2017v2p1VSjet_.clear();
    tauByMediumDeepTau2017v2p1VSjet_.clear();
    tauByTightDeepTau2017v2p1VSjet_.clear();
    tauByVTightDeepTau2017v2p1VSjet_.clear();
    tauByVVTightDeepTau2017v2p1VSjet_.clear();
    tauByVVVLooseDeepTau2017v2p1VSe_.clear();
    tauByVVLooseDeepTau2017v2p1VSe_.clear();
    tauByVLooseDeepTau2017v2p1VSe_.clear();
    tauByLooseDeepTau2017v2p1VSe_.clear();
    tauByMediumDeepTau2017v2p1VSe_.clear();
    tauByTightDeepTau2017v2p1VSe_.clear();
    tauByVTightDeepTau2017v2p1VSe_.clear();
    tauByVVTightDeepTau2017v2p1VSe_.clear();
    //tauByVVVLooseDeepTau2017v2p1VSmu_.clear();
    //tauByVVLooseDeepTau2017v2p1VSmu_.clear();
    tauByVLooseDeepTau2017v2p1VSmu_.clear();
    tauByLooseDeepTau2017v2p1VSmu_.clear();
    tauByMediumDeepTau2017v2p1VSmu_.clear();
    tauByTightDeepTau2017v2p1VSmu_.clear();
    //tauByVTightDeepTau2017v2p1VSmu_.clear();
    //tauByVVTightDeepTau2017v2p1VSmu_.clear();
    
    tauByIsolationMVArun2v1DBnewDMwLTraw_.clear();
    tauByIsolationMVArun2v1DBoldDMwLTraw_.clear();
    tauByIsolationMVArun2v1PWnewDMwLTraw_.clear();
    tauByIsolationMVArun2v1PWoldDMwLTraw_.clear();
    tauByVTightIsolationMVArun2v1DBnewDMwLT_.clear();
    tauByVTightIsolationMVArun2v1DBoldDMwLT_.clear();
    tauByVTightIsolationMVArun2v1PWnewDMwLT_.clear();
    tauByVTightIsolationMVArun2v1PWoldDMwLT_.clear();
    tauByTightIsolationMVArun2v1DBnewDMwLT_.clear();
    tauByTightIsolationMVArun2v1DBoldDMwLT_.clear();
    tauByTightIsolationMVArun2v1PWnewDMwLT_.clear();
    tauByTightIsolationMVArun2v1PWoldDMwLT_.clear();
    tauByMediumIsolationMVArun2v1DBnewDMwLT_.clear();
    tauByMediumIsolationMVArun2v1DBoldDMwLT_.clear();
    tauByMediumIsolationMVArun2v1PWnewDMwLT_.clear();
    tauByMediumIsolationMVArun2v1PWoldDMwLT_.clear();
    tauByLooseIsolationMVArun2v1DBnewDMwLT_.clear();
    tauByLooseIsolationMVArun2v1DBoldDMwLT_.clear();
    tauByLooseIsolationMVArun2v1PWnewDMwLT_.clear();
    tauByLooseIsolationMVArun2v1PWoldDMwLT_.clear();
    tauByVLooseIsolationMVArun2v1DBnewDMwLT_.clear();
    tauByVLooseIsolationMVArun2v1DBoldDMwLT_.clear();
    tauByVLooseIsolationMVArun2v1PWnewDMwLT_.clear();
    tauByVLooseIsolationMVArun2v1PWoldDMwLT_.clear();
    
    //Tau Kinematics
    tauEta_.clear();
    tauPhi_.clear();
    tauPt_.clear();
    tauEt_.clear();
    tauCharge_.clear();
    tauP_.clear();
    tauPx_.clear();
    tauPy_.clear();
    tauPz_.clear();
    tauVz_.clear();
    tauEnergy_.clear();
    tauMass_.clear();
    tauDxy_.clear();
    tauZImpact_.clear();
    
    // Tau Ingredients
    tauDecayMode_.clear();
    tauLeadChargedHadronExists_.clear();
    tauLeadChargedHadronEta_.clear();
    tauLeadChargedHadronPhi_.clear();
    tauLeadChargedHadronPt_.clear();
    tauChargedIsoPtSum_.clear();
    tauNeutralIsoPtSum_.clear();
    tauPuCorrPtSum_.clear();
    tauNumSignalPFChargedHadrCands_.clear();
    tauNumSignalPFNeutrHadrCands_.clear();
    tauNumSignalPFGammaCands_.clear();
    tauNumSignalPFCands_.clear();
    tauNumIsolationPFChargedHadrCands_.clear();
    tauNumIsolationPFNeutrHadrCands_.clear();
    tauNumIsolationPFGammaCands_.clear();
    tauNumIsolationPFCands_.clear();
    
    tauneutralIsoPtSumWeight_.clear();
    taufootprintCorrection_.clear();
    tauphotonPtSumOutsideSignalCone_.clear();
    taudz_.clear();
    taudxy_.clear();
    
    
    nTau_ = 0;
    
    //startTaus L1
/***
    edm::Handle<BXVector<l1t::Tau>> tauL1Handle;
    e.getByToken(tauTriggerPrimitives_, tauL1Handle);

    for(BXVector<l1t::Tau>::const_iterator itau = tauL1Handle->begin(); itau != tauL1Handle->end(); ++itau) {
	tauL1PrimitivesPt_.push_back(itau->pt());
	tauL1PrimitivesEta_.push_back(itau->eta());
	tauL1PrimitivesPhi_.push_back(itau->phi());
	tauL1PrimitivesEnergy_.push_back(itau->energy());
	tauL1PrimitivesIso_.push_back(itau->hwIso());
    }
***/
    //original Tau lines

    edm::Handle<vector<pat::Tau> > tauHandle;
    e.getByToken(tauCollection_, tauHandle);
    
    if (!tauHandle.isValid()) {
        edm::LogWarning("ggNtuplizer") << "no pat::Tau in event";
        return;
    }
    
    //startTaus Lvdp
    for(vector<pat::Tau>::const_iterator itau = tauHandle->begin(); itau != tauHandle->end(); ++itau) {
        
        // Tau Id & Isolation
        taupfTausDiscriminationByDecayModeFinding_.push_back(itau->tauID("decayModeFinding"));
        taupfTausDiscriminationByDecayModeFindingNewDMs_.push_back(itau->tauID("decayModeFindingNewDMs"));
        
        tauByMVA6VLooseElectronRejection_.push_back(itau->tauID("againstElectronVLooseMVA6"));
        tauByMVA6LooseElectronRejection_.push_back(itau->tauID("againstElectronLooseMVA6"));
        tauByMVA6MediumElectronRejection_.push_back(itau->tauID("againstElectronMediumMVA6"));
        tauByMVA6TightElectronRejection_.push_back(itau->tauID("againstElectronTightMVA6"));
        tauByMVA6VTightElectronRejection_.push_back(itau->tauID("againstElectronVTightMVA6"));
        
        tauByLooseMuonRejection3_.push_back(itau->tauID("againstMuonLoose3"));
        tauByTightMuonRejection3_.push_back(itau->tauID("againstMuonTight3"));
        
        tauByLooseCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byLooseCombinedIsolationDeltaBetaCorr3Hits"));
        tauByMediumCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byMediumCombinedIsolationDeltaBetaCorr3Hits"));
        tauByTightCombinedIsolationDeltaBetaCorr3Hits_.push_back(itau->tauID("byTightCombinedIsolationDeltaBetaCorr3Hits"));
        tauCombinedIsolationDeltaBetaCorrRaw3Hits_.push_back(itau->tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"));
        
        tauByVVVLooseDeepTau2017v2p1VSjet_.push_back(itau->tauID("byVVVLooseDeepTau2017v2p1VSjet"));
        tauByVVLooseDeepTau2017v2p1VSjet_.push_back(itau->tauID("byVVLooseDeepTau2017v2p1VSjet"));
        tauByVLooseDeepTau2017v2p1VSjet_.push_back(itau->tauID("byVLooseDeepTau2017v2p1VSjet"));
        tauByLooseDeepTau2017v2p1VSjet_.push_back(itau->tauID("byLooseDeepTau2017v2p1VSjet"));
        tauByMediumDeepTau2017v2p1VSjet_.push_back(itau->tauID("byMediumDeepTau2017v2p1VSjet"));
        tauByTightDeepTau2017v2p1VSjet_.push_back(itau->tauID("byTightDeepTau2017v2p1VSjet"));
        tauByVTightDeepTau2017v2p1VSjet_.push_back(itau->tauID("byVTightDeepTau2017v2p1VSjet"));
        tauByVVTightDeepTau2017v2p1VSjet_.push_back(itau->tauID("byVVTightDeepTau2017v2p1VSjet"));
        tauByVVVLooseDeepTau2017v2p1VSe_.push_back(itau->tauID("byVVVLooseDeepTau2017v2p1VSe"));
        tauByVVLooseDeepTau2017v2p1VSe_.push_back(itau->tauID("byVVLooseDeepTau2017v2p1VSe"));
        tauByVLooseDeepTau2017v2p1VSe_.push_back(itau->tauID("byVLooseDeepTau2017v2p1VSe"));
        tauByLooseDeepTau2017v2p1VSe_.push_back(itau->tauID("byLooseDeepTau2017v2p1VSe"));
        tauByMediumDeepTau2017v2p1VSe_.push_back(itau->tauID("byMediumDeepTau2017v2p1VSe"));
        tauByTightDeepTau2017v2p1VSe_.push_back(itau->tauID("byTightDeepTau2017v2p1VSe"));
        tauByVTightDeepTau2017v2p1VSe_.push_back(itau->tauID("byVTightDeepTau2017v2p1VSe"));
        tauByVVTightDeepTau2017v2p1VSe_.push_back(itau->tauID("byVVTightDeepTau2017v2p1VSe"));
        //tauByVVVLooseDeepTau2017v2p1VSmu_.push_back(itau->tauID("byVVVLooseDeepTau2017v2p1VSmu"));
        //tauByVVLooseDeepTau2017v2p1VSmu_.push_back(itau->tauID("byVVLooseDeepTau2017v2p1VSmu"));
        tauByVLooseDeepTau2017v2p1VSmu_.push_back(itau->tauID("byVLooseDeepTau2017v2p1VSmu"));
        tauByLooseDeepTau2017v2p1VSmu_.push_back(itau->tauID("byLooseDeepTau2017v2p1VSmu"));
        tauByMediumDeepTau2017v2p1VSmu_.push_back(itau->tauID("byMediumDeepTau2017v2p1VSmu"));
        tauByTightDeepTau2017v2p1VSmu_.push_back(itau->tauID("byTightDeepTau2017v2p1VSmu"));
        //tauByVTightDeepTau2017v2p1VSmu_.push_back(itau->tauID("byVTightDeepTau2017v2p1VSmu"));
        //tauByVVTightDeepTau2017v2p1VSmu_.push_back(itau->tauID("byVVTightDeepTau2017v2p1VSmu"));
   
 
        tauByIsolationMVArun2v1DBnewDMwLTraw_.push_back(itau->tauID("byIsolationMVArun2v1DBnewDMwLTraw"));
        tauByIsolationMVArun2v1DBoldDMwLTraw_.push_back(itau->tauID("byIsolationMVArun2v1DBoldDMwLTraw"));
        //tauByIsolationMVArun2v1PWnewDMwLTraw_.push_back(itau->tauID("byIsolationMVArun2v1PWnewDMwLTraw"));
        //tauByIsolationMVArun2v1PWoldDMwLTraw_.push_back(itau->tauID("byIsolationMVArun2v1PWoldDMwLTraw"));
        tauByVLooseIsolationMVArun2v1DBnewDMwLT_.push_back(itau->tauID("byVLooseIsolationMVArun2v1DBnewDMwLT"));
        tauByVLooseIsolationMVArun2v1DBoldDMwLT_.push_back(itau->tauID("byVLooseIsolationMVArun2v1DBoldDMwLT"));
        //tauByVLooseIsolationMVArun2v1PWnewDMwLT_.push_back(itau->tauID("byVLooseIsolationMVArun2v1PWnewDMwLT"));
        //tauByVLooseIsolationMVArun2v1PWoldDMwLT_.push_back(itau->tauID("byVLooseIsolationMVArun2v1PWoldDMwLT"));
        tauByLooseIsolationMVArun2v1DBnewDMwLT_.push_back(itau->tauID("byLooseIsolationMVArun2v1DBnewDMwLT"));
        tauByLooseIsolationMVArun2v1DBoldDMwLT_.push_back(itau->tauID("byLooseIsolationMVArun2v1DBoldDMwLT"));
        //tauByLooseIsolationMVArun2v1PWnewDMwLT_.push_back(itau->tauID("byLooseIsolationMVArun2v1PWnewDMwLT"));
        //tauByLooseIsolationMVArun2v1PWoldDMwLT_.push_back(itau->tauID("byLooseIsolationMVArun2v1PWoldDMwLT"));
        tauByMediumIsolationMVArun2v1DBnewDMwLT_.push_back(itau->tauID("byMediumIsolationMVArun2v1DBnewDMwLT"));
        tauByMediumIsolationMVArun2v1DBoldDMwLT_.push_back(itau->tauID("byMediumIsolationMVArun2v1DBoldDMwLT"));
        //tauByMediumIsolationMVArun2v1PWnewDMwLT_.push_back(itau->tauID("byMediumIsolationMVArun2v1PWnewDMwLT"));
        //tauByMediumIsolationMVArun2v1PWoldDMwLT_.push_back(itau->tauID("byMediumIsolationMVArun2v1PWoldDMwLT"));
        tauByTightIsolationMVArun2v1DBnewDMwLT_.push_back(itau->tauID("byTightIsolationMVArun2v1DBnewDMwLT"));
        tauByTightIsolationMVArun2v1DBoldDMwLT_.push_back(itau->tauID("byTightIsolationMVArun2v1DBoldDMwLT"));
        //tauByTightIsolationMVArun2v1PWnewDMwLT_.push_back(itau->tauID("byTightIsolationMVArun2v1PWnewDMwLT"));
        //tauByTightIsolationMVArun2v1PWoldDMwLT_.push_back(itau->tauID("byTightIsolationMVArun2v1PWoldDMwLT"));
        tauByVTightIsolationMVArun2v1DBnewDMwLT_.push_back(itau->tauID("byVTightIsolationMVArun2v1DBnewDMwLT"));
        tauByVTightIsolationMVArun2v1DBoldDMwLT_.push_back(itau->tauID("byVTightIsolationMVArun2v1DBoldDMwLT"));
        //tauByVTightIsolationMVArun2v1PWnewDMwLT_.push_back(itau->tauID("byVTightIsolationMVArun2v1PWnewDMwLT"));
        //tauByVTightIsolationMVArun2v1PWoldDMwLT_.push_back(itau->tauID("byVTightIsolationMVArun2v1PWoldDMwLT"));
    
        
        //Tau Kinematics
        tauEta_.push_back(itau->eta());
        tauPhi_.push_back(itau->phi());
        tauPt_.push_back(itau->pt());
        tauEt_.push_back(itau->et());
        tauCharge_.push_back(itau->charge());
        tauP_.push_back(itau->p() );
        tauPx_.push_back(itau->px() );
        tauPy_.push_back(itau->py() );
        tauPz_.push_back(itau->pz() );
        tauVz_.push_back(itau->vz() );
        tauEnergy_.push_back(itau->energy() );
        tauMass_.push_back(itau->mass());
        tauDxy_.push_back(itau->dxy() );
        tauZImpact_.push_back(itau->vertex().z() + 130./tan(itau->theta()));
        
        
        // Tau Ingredients
        tauDecayMode_.push_back(itau->decayMode());
        tauChargedIsoPtSum_.push_back(itau->tauID("chargedIsoPtSum") );
        tauNeutralIsoPtSum_.push_back(itau->tauID("neutralIsoPtSum")  );
        tauPuCorrPtSum_.push_back(itau->tauID("puCorrPtSum")  );
        tauneutralIsoPtSumWeight_.push_back(itau->tauID("neutralIsoPtSumWeight"));
        taufootprintCorrection_.push_back(itau->tauID("footprintCorrection"));
        tauphotonPtSumOutsideSignalCone_.push_back(itau->tauID("photonPtSumOutsideSignalCone"));
        
        tauNumSignalPFChargedHadrCands_.push_back(itau->signalChargedHadrCands().size());
        tauNumSignalPFNeutrHadrCands_.push_back(itau->signalNeutrHadrCands().size());
        tauNumSignalPFGammaCands_.push_back(itau->signalGammaCands().size());
        tauNumSignalPFCands_.push_back(itau->signalCands().size());
        
        tauNumIsolationPFChargedHadrCands_.push_back(itau->isolationChargedHadrCands().size());
        tauNumIsolationPFNeutrHadrCands_.push_back(itau->isolationNeutrHadrCands().size());
        tauNumIsolationPFGammaCands_.push_back(itau->isolationGammaCands().size());
        tauNumIsolationPFCands_.push_back(itau->isolationCands().size());
        
        edm::Handle<reco::VertexCollection> vertexs;
        e.getByToken(vtxLabel_, vertexs);
        
        if (vertexs->size()>0) {
            pat::PackedCandidate const* packedLeadTauCand = dynamic_cast<pat::PackedCandidate const*>(itau->leadChargedHadrCand().get());
            taudz_.push_back(packedLeadTauCand->dz());
            taudxy_.push_back(packedLeadTauCand->dxy());
            tauLeadChargedHadronExists_.push_back(true);
            tauLeadChargedHadronEta_.push_back(packedLeadTauCand->eta());
            tauLeadChargedHadronPhi_.push_back(packedLeadTauCand->phi());
            tauLeadChargedHadronPt_.push_back(packedLeadTauCand->pt());
        }
        
        
        
        
        ++nTau_;
        
    } // loop over tau candidates
    
}
