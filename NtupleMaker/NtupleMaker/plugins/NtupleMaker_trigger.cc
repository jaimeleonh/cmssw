// Braden Allmond Nov 22nd 2021
// NtupleMaker
// Takes a dataset and makes a tree of branches filled with data
// In this case, our NtupleMaker is for an HLT study,
// so we have branches for multiple HLT filters. 
// Here's a diagram of "triggers" and "filters"
//   Trigger (HLT means High Level Trigger)
// ---------------------------------------------------------------------------
// |L1 Decision| Middle Filter 1| Middle Filter 2| Final Filter| HLT Decision|
// ---------------------------------------------------------------------------
// Some filters/modules are shared between HLT paths,
// meaning sometimes those filters are only triggered
// by one path and sometimes they're triggered by both.
// From one filter's information alone, it's not possible to 
// tell which path the filter was triggered by. If you want to 
// know which filter is triggered by which path, you have to
// daisy-chain the filter decisions for a path explcitly, 
// meaning you check each filter decision in a path before
// the one you care about. If filter in the path before the
// one you care about is passed, then the filter you're looking at
// was triggered in the path you're studying in. That would look
// like this.
//    SomeHLTPath
// -----------------
// |1|1|1|1|0|0|
// -----------------
// Above, we can see that the fourth filter was triggered by
// this path, because each filter before that was triggered
// by this path as well. If we somehow find something like 
// the following in our analysis
//   SomeOtherHLTPath
// -----------------
// |0|0|0|1|0|1|0|0|
// -----------------
// We can safely those filters were not triggered by the path
// we're looking at, and were instead triggered by a different
// path with shared modules. 
//
// Note: I changed all branch names to use exact module names from HLT.
// This seemed more straightforward than coming up with good variable names.
// I'll make a table/sheet of the module names and what they do.
// I'll also update all macros/trigger_trees/and analyzers that are
// affected by this branch name rewriting.
//
// InclusiveVBF = Old VBF = L1_DoubleJet_110_35_DoubleJet35_Mass_Min620
//   it's called inclusive VBF bc a two jet L1 seed includes VBF events with any final lepton state
// VBFPlusTwoTau = New VBF/Proposed VBF = L1_DoubleJet35_Mass_Min420_IsoTau45er2p1_RmvOl
//   includes VBF events with two hadronic taus in final lepton state
// VBFPlusOneTau = same L1 as above
//   includes VBF events with two hadronic taus, one hadronic tau one muon, or one hadronic tau one electron in final state

#include "NtupleMaker/NtupleMaker/plugins/NtupleMaker.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/Common/interface/TriggerResults.h"

using namespace std; // I think best practice is to include <vector> explicitly at the top of the file

// int passDiTau35HLT;
// int passInclusiveVBFHLT;
// int passVBFPlusTwoTauHLT;
// int passVBFPlusOneTauHLT;
int passDoubleTauJetHLT;
int passDoubleTau32HLT;
int passDoubleTau35HLT;

// std::vector <std::string> hlt_paths = {
// };

std::map <std::string, int> pass_hlt;

// 

float 	pt_;
float 	eta_;
float 	phi_;
float 	energy_;

int	nEvents;
//full trigger filter name is commented above each of its branches
//hltL1sDoubleTauBigOR
// int 		passhltL1sDoubleTauBigOR;
// vector<float>	hltL1sDoubleTauBigOR_pt;
// vector<float>	hltL1sDoubleTauBigOR_eta;
// vector<float>	hltL1sDoubleTauBigOR_phi;
// vector<float>	hltL1sDoubleTauBigOR_energy;
// //hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg
// int		passhltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg;
// vector<float>	hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_pt;
// vector<float>	hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_eta;
// vector<float>	hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_phi;
// vector<float>	hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_energy;
// //hltL1VBFDiJetOR
// int	passhltL1VBFDiJetOR;
// vector<float>	hltL1VBFDiJetOR_pt;
// vector<float>	hltL1VBFDiJetOR_eta;
// vector<float>	hltL1VBFDiJetOR_phi;
// vector<float>	hltL1VBFDiJetOR_energy;
// //hltL1VBFDiJetIsoTau (previously called hltL1NewVBFDiJet) 
// int	passhltL1VBFDiJetIsoTau;
// vector<int>	hltL1VBFDiJetIsoTau_nJets;
// vector<float>	hltL1VBFDiJetIsoTau_jetPt;
// vector<float>	hltL1VBFDiJetIsoTau_jetEta;
// vector<float>	hltL1VBFDiJetIsoTau_jetPhi;
// vector<float>	hltL1VBFDiJetIsoTau_jetEnergy;
// vector<int>	hltL1VBFDiJetIsoTau_nTaus;
// vector<float>	hltL1VBFDiJetIsoTau_tauPt;
// vector<float>	hltL1VBFDiJetIsoTau_tauEta;
// vector<float>	hltL1VBFDiJetIsoTau_tauPhi;
// vector<float>	hltL1VBFDiJetIsoTau_tauEnergy;

//hltL1DoubleTau32
int	passhltL1DoubleTau32;
int	hltL1DoubleTau32_nTaus;
vector<float>	hltL1DoubleTau32_tauPt;
vector<float>	hltL1DoubleTau32_tauEta;
vector<float>	hltL1DoubleTau32_tauPhi;
vector<float>	hltL1DoubleTau32_tauEnergy;

//hltL1DoubleTau35
int	passhltL1DoubleTau35;
int	hltL1DoubleTau35_nTaus;
vector<float>	hltL1DoubleTau35_tauPt;
vector<float>	hltL1DoubleTau35_tauEta;
vector<float>	hltL1DoubleTau35_tauPhi;
vector<float>	hltL1DoubleTau35_tauEnergy;

//hltL1DoubleTauJet
int	passhltL1DoubleTauJet;
int	hltL1DoubleTauJet_nJets;
vector<float>	hltL1DoubleTauJet_jetPt;
vector<float>	hltL1DoubleTauJet_jetEta;
vector<float>	hltL1DoubleTauJet_jetPhi;
vector<float>	hltL1DoubleTauJet_jetEnergy;
int	hltL1DoubleTauJet_nTaus;
vector<float>	hltL1DoubleTauJet_tauPt;
vector<float>	hltL1DoubleTauJet_tauEta;
vector<float>	hltL1DoubleTauJet_tauPhi;
vector<float>	hltL1DoubleTauJet_tauEnergy;

int passhltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet;
int hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_nObjs;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_pt;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_eta;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_phi;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_energy;

int passhltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet;
int hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_nObjs;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_pt;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_eta;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_phi;
vector<float> hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_energy;

int passhltHpsPFTauTrack;
int hltHpsPFTauTrack_nObjs;
vector<float> hltHpsPFTauTrack_pt;
vector<float> hltHpsPFTauTrack_eta;
vector<float> hltHpsPFTauTrack_phi;
vector<float> hltHpsPFTauTrack_energy;

int passhltL2DoubleTauTagNNFilterDoubleTauJet;
int hltL2DoubleTauTagNNFilterDoubleTauJet_nObjs;
vector<float> hltL2DoubleTauTagNNFilterDoubleTauJet_pt;
vector<float> hltL2DoubleTauTagNNFilterDoubleTauJet_eta;
vector<float> hltL2DoubleTauTagNNFilterDoubleTauJet_phi;
vector<float> hltL2DoubleTauTagNNFilterDoubleTauJet_energy;

int passhltPFJets60L1HLTMatched;
int hltPFJets60L1HLTMatched_nObjs;
vector<float> hltPFJets60L1HLTMatched_pt;
vector<float> hltPFJets60L1HLTMatched_eta;
vector<float> hltPFJets60L1HLTMatched_phi;
vector<float> hltPFJets60L1HLTMatched_energy;

int passhltHpsOverlapFilterDeepTauDoublePFTau30PFJet60;
int hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_nObjs;
vector<float> hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_pt;
vector<float> hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_eta;
vector<float> hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_phi;
vector<float> hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_energy;




void NtupleMaker::branchesTriggers(TTree* tree){
    
    tree->Branch("nEvents", &nEvents);
    
    // tree->Branch("passDiTau35HLT", &passDiTau35HLT);
    
    // tree->Branch("passInclusiveVBFHLT", &passInclusiveVBFHLT);    
    
    // tree->Branch("passVBFPlusTwoTauHLT", &passVBFPlusTwoTauHLT);
    
    tree->Branch("passDoubleTauJetHLT", &passDoubleTauJetHLT);
    tree->Branch("passDoubleTau32HLT", &passDoubleTau32HLT);
    tree->Branch("passDoubleTau35HLT", &passDoubleTau35HLT);
    
    // for (auto trigger: hlt_paths) {
        // tree->Branch(trigger, &int(pass_hlt[trigger]));
    // }
    // for (auto trigger: hlt_paths) {
        // pass_hlt[trigger] = 0;
        // // tree->Branch(trigger, &pass_hlt[trigger]);
    // }
    
    // std::map<std::string, int>::iterator it;
    // for (it = pass_hlt.begin(); it != pass_hlt.end(); it++){
         
    // }

    // for (auto const &trigger: pass_hlt) {
        // tree->Branch(trigger.first.c_str(), trigger.second);
    // }
    
    // tree->Branch("passVBFPlusOneTauHLT", &passVBFPlusOneTauHLT);
    
    // tree->Branch("passhltL1sDoubleTauBigOR", &passhltL1sDoubleTauBigOR);
    // tree->Branch("hltL1sDoubleTauBigOR_pt", &hltL1sDoubleTauBigOR_pt);
    // tree->Branch("hltL1sDoubleTauBigOR_eta", &hltL1sDoubleTauBigOR_eta);
    // tree->Branch("hltL1sDoubleTauBigOR_phi", &hltL1sDoubleTauBigOR_phi);
    // tree->Branch("hltL1sDoubleTauBigOR_energy", &hltL1sDoubleTauBigOR_energy);
    
    // // don't need the cutflow, so I'm only storing the final filter of DiTau35HLT which is used for HLT-AOD matching purposes
    // tree->Branch("passhltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg", &passhltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg);  
    // tree->Branch("hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_pt", &hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_pt);  
    // tree->Branch("hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_eta", &hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_eta);  
    // tree->Branch("hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_phi", &hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_phi);  
    // tree->Branch("hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_energy", &hltHpsDoublePFTau35TrackPt1TightChargedIsolationL1HLTMatchedReg_energy);  
    
    // tree->Branch("passhltL1VBFDiJetOR", &passhltL1VBFDiJetOR);
    // tree->Branch("hltL1VBFDiJetOR_pt", &hltL1VBFDiJetOR_pt);
    // tree->Branch("hltL1VBFDiJetOR_eta", &hltL1VBFDiJetOR_eta);
    // tree->Branch("hltL1VBFDiJetOR_phi", &hltL1VBFDiJetOR_phi);
    // tree->Branch("hltL1VBFDiJetOR_energy", &hltL1VBFDiJetOR_energy);
    
    tree->Branch("passhltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet", &passhltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_nObjs", &hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_nObjs);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_pt", &hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_pt);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_eta", &hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_eta);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_phi", &hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_phi);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_energy", &hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_energy);

    tree->Branch("passhltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet", &passhltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_nObjs", &hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_nObjs);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_pt", &hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_pt);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_eta", &hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_eta);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_phi", &hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_phi);
    tree->Branch("hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_energy", &hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_energy);

    tree->Branch("passhltHpsPFTauTrack", &passhltHpsPFTauTrack);
    tree->Branch("hltHpsPFTauTrack_nObjs", &hltHpsPFTauTrack_nObjs);
    tree->Branch("hltHpsPFTauTrack_pt", &hltHpsPFTauTrack_pt);
    tree->Branch("hltHpsPFTauTrack_eta", &hltHpsPFTauTrack_eta);
    tree->Branch("hltHpsPFTauTrack_phi", &hltHpsPFTauTrack_phi);
    tree->Branch("hltHpsPFTauTrack_energy", &hltHpsPFTauTrack_energy);

    tree->Branch("passhltL2DoubleTauTagNNFilterDoubleTauJet", &passhltL2DoubleTauTagNNFilterDoubleTauJet);
    tree->Branch("hltL2DoubleTauTagNNFilterDoubleTauJet_nObjs", &hltL2DoubleTauTagNNFilterDoubleTauJet_nObjs);
    tree->Branch("hltL2DoubleTauTagNNFilterDoubleTauJet_pt", &hltL2DoubleTauTagNNFilterDoubleTauJet_pt);
    tree->Branch("hltL2DoubleTauTagNNFilterDoubleTauJet_eta", &hltL2DoubleTauTagNNFilterDoubleTauJet_eta);
    tree->Branch("hltL2DoubleTauTagNNFilterDoubleTauJet_phi", &hltL2DoubleTauTagNNFilterDoubleTauJet_phi);
    tree->Branch("hltL2DoubleTauTagNNFilterDoubleTauJet_energy", &hltL2DoubleTauTagNNFilterDoubleTauJet_energy);

    tree->Branch("passhltPFJets60L1HLTMatched", &passhltPFJets60L1HLTMatched);
    tree->Branch("hltPFJets60L1HLTMatched_nObjs", &hltPFJets60L1HLTMatched_nObjs);
    tree->Branch("hltPFJets60L1HLTMatched_pt", &hltPFJets60L1HLTMatched_pt);
    tree->Branch("hltPFJets60L1HLTMatched_eta", &hltPFJets60L1HLTMatched_eta);
    tree->Branch("hltPFJets60L1HLTMatched_phi", &hltPFJets60L1HLTMatched_phi);
    tree->Branch("hltPFJets60L1HLTMatched_energy", &hltPFJets60L1HLTMatched_energy);

    tree->Branch("passhltHpsOverlapFilterDeepTauDoublePFTau30PFJet60", &passhltHpsOverlapFilterDeepTauDoublePFTau30PFJet60);
    tree->Branch("hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_nObjs", &hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_nObjs);
    tree->Branch("hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_pt", &hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_pt);
    tree->Branch("hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_eta", &hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_eta);
    tree->Branch("hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_phi", &hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_phi);
    tree->Branch("hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_energy", &hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_energy);


}

void NtupleMaker::fillTriggers(const edm::Event& iEvent){
    
    using namespace edm;
    
    // std::cout << "Filling triggers" << std::endl;
    
    // clearing vectors and initializing flags at the start of every event 
    nEvents = 0;
    
    // passDiTau35HLT = 0;
    // passInclusiveVBFHLT = 0; 
    // passVBFPlusTwoTauHLT = 0;
    // passVBFPlusOneTauHLT = 0;
    passDoubleTauJetHLT = 0;
    passDoubleTau32HLT = 0;
    passDoubleTau35HLT = 0;
    
    // other HLT paths
    // for (auto trigger: hlt_paths) {
        // pass_hlt[trigger] = 0;
    // }
    
    //L1 branches 
    passhltL1DoubleTau32 = 0;
    hltL1DoubleTau32_nTaus = 0;
    hltL1DoubleTau32_tauPt.clear();
    hltL1DoubleTau32_tauEta.clear();
    hltL1DoubleTau32_tauPhi.clear();
    hltL1DoubleTau32_tauEnergy.clear();

    //hltL1DoubleTau35
    passhltL1DoubleTau35 = 0;
    hltL1DoubleTau35_nTaus = 0;
    hltL1DoubleTau35_tauPt.clear();
    hltL1DoubleTau35_tauEta.clear();
    hltL1DoubleTau35_tauPhi.clear();
    hltL1DoubleTau35_tauEnergy.clear();

    //hltL1DoubleTauJet
    passhltL1DoubleTauJet = 0;
    hltL1DoubleTauJet_nJets = 0;
    hltL1DoubleTauJet_jetPt.clear();
    hltL1DoubleTauJet_jetEta.clear();
    hltL1DoubleTauJet_jetPhi.clear();
    hltL1DoubleTauJet_jetEnergy.clear();
    hltL1DoubleTauJet_nTaus = 0;
    hltL1DoubleTauJet_tauPt.clear();
    hltL1DoubleTauJet_tauEta.clear();
    hltL1DoubleTauJet_tauPhi.clear();
    hltL1DoubleTauJet_tauEnergy.clear();

    passhltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet = 0;
    hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_nObjs = 0;
    hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_pt.clear();
    hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_eta.clear();
    hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_phi.clear();
    hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_energy.clear();

    passhltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet = 0;
    hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_nObjs = 0;
    hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_pt.clear();
    hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_eta.clear();
    hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_phi.clear();
    hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_energy.clear();

    passhltHpsPFTauTrack = 0;
    hltHpsPFTauTrack_nObjs = 0;
    hltHpsPFTauTrack_pt.clear();
    hltHpsPFTauTrack_eta.clear();
    hltHpsPFTauTrack_phi.clear();
    hltHpsPFTauTrack_energy.clear();

    passhltL2DoubleTauTagNNFilterDoubleTauJet = 0;
    hltL2DoubleTauTagNNFilterDoubleTauJet_nObjs = 0;
    hltL2DoubleTauTagNNFilterDoubleTauJet_pt.clear();
    hltL2DoubleTauTagNNFilterDoubleTauJet_eta.clear();
    hltL2DoubleTauTagNNFilterDoubleTauJet_phi.clear();
    hltL2DoubleTauTagNNFilterDoubleTauJet_energy.clear();

    passhltPFJets60L1HLTMatched = 0;
    hltPFJets60L1HLTMatched_nObjs = 0;
    hltPFJets60L1HLTMatched_pt.clear();
    hltPFJets60L1HLTMatched_eta.clear();
    hltPFJets60L1HLTMatched_phi.clear();
    hltPFJets60L1HLTMatched_energy.clear();

    passhltHpsOverlapFilterDeepTauDoublePFTau30PFJet60 = 0;
    hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_nObjs = 0;
    hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_pt.clear();
    hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_eta.clear();
    hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_phi.clear();
    hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_energy.clear();


    
    
    // getting trigger results, following this page
    // https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideHLTAnalysis
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByToken(triggerResultToken_, triggerResults);
    edm::Handle<trigger::TriggerEvent> triggerEvent;
    iEvent.getByToken(triggerEventToken_, triggerEvent);
    const edm::TriggerNames triggerNames_ = iEvent.triggerNames(*triggerResults);
    
    // std::cout << "Get by token" << std::endl;
    
    // // saving trigger results to respective branches
    
    // 2 Tau + Jet
    std::string pathNameDoubleTauJet = "HLT_DoubleMediumDeepTauIsoPFTauHPS30_L2NN_eta2p1_PFJet60_v2";
    passDoubleTauJetHLT = triggerResults->accept(triggerNames_.triggerIndex(pathNameDoubleTauJet));
    
    // // 2 Tau, L1 32
    // std::string pathNameDoubleTau32 = "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_v1";
    // passDoubleTau32HLT = triggerResults->accept(triggerNames_.triggerIndex(pathNameDoubleTau32));
    
    // // 2 Tau, L1 35
    // std::string pathNameDoubleTau35 = "HLT_DoubleMediumChargedIsoPFTauHPS40_Trk1_TightID_eta2p1_Reg_L135";
    // passDoubleTau35HLT = triggerResults->accept(triggerNames_.triggerIndex(pathNameDoubleTau35));
    
    // filling branches with triggerObjs information, hltL1VBFDiJetIsoTau object info filled separately since it's a weird L1
    
    // getting trigger refs for hltL1DoubleTauJet filter's tau/jet splitting
    edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
    iEvent.getByToken(triggerEventWithRefsToken_, triggerEventWithRefsHandle_);
    const unsigned int filterIndex(triggerEventWithRefsHandle_->filterIndex(InputTag("hltL1sDoubleTauJet", "", "MYHLT")));
    //making jet object and filling vector
    l1t::JetVectorRef jetCandRefVec;
    trigger::Vids jvids;
    triggerEventWithRefsHandle_->getObjects(filterIndex, jvids, jetCandRefVec);
    //making tau object and filling vector
    l1t::TauVectorRef tauCandRefVec;
    trigger::Vids tvids;
    triggerEventWithRefsHandle_->getObjects(filterIndex, tvids, tauCandRefVec);

    const unsigned int nJets(jvids.size());
    hltL1DoubleTauJet_nJets = nJets;
    if (nJets > 0) {
        for (unsigned int i = 0; i != nJets; ++i) {
            hltL1DoubleTauJet_jetPt.push_back(jetCandRefVec[i]->pt());
            hltL1DoubleTauJet_jetEta.push_back(jetCandRefVec[i]->eta());
            hltL1DoubleTauJet_jetPhi.push_back(jetCandRefVec[i]->phi());
            hltL1DoubleTauJet_jetEnergy.push_back(jetCandRefVec[i]->energy());
        }
    }
    const unsigned int nTaus(tvids.size());
    hltL1DoubleTauJet_nTaus = nTaus;
    if (nTaus > 0) {
        for (unsigned int i = 0; i != nTaus; ++i) {
            hltL1DoubleTauJet_tauPt.push_back(tauCandRefVec[i]->pt());
            hltL1DoubleTauJet_tauEta.push_back(tauCandRefVec[i]->eta());
            hltL1DoubleTauJet_tauPhi.push_back(tauCandRefVec[i]->phi());
            hltL1DoubleTauJet_tauEnergy.push_back(tauCandRefVec[i]->energy());
        }
    }
    
    
    // // getting trigger refs for hltL1sDoubleTauBigOR filter's tau/jet splitting
    // const unsigned int filterIndex32(triggerEventWithRefsHandle_->filterIndex(InputTag("hltL1sDoubleTauBigOR", "", "MYHLT")));
    // //making tau object and filling vector
    // triggerEventWithRefsHandle_->getObjects(filterIndex32, tvids, tauCandRefVec);

    // const unsigned int nTaus32(tvids.size());
    // hltL1DoubleTau32_nTaus = nTaus32;
    // if (nTaus > 0) {
        // for (unsigned int i = 0; i != nTaus32; ++i) {
            // hltL1DoubleTau32_tauPt.push_back(tauCandRefVec[i]->pt());
            // hltL1DoubleTau32_tauEta.push_back(tauCandRefVec[i]->eta());
            // hltL1DoubleTau32_tauPhi.push_back(tauCandRefVec[i]->phi());
            // hltL1DoubleTau32_tauEnergy.push_back(tauCandRefVec[i]->energy());
        // }
    // }
    
    // // getting trigger refs for hltL1DoubleTau35 filter's tau/jet splitting
    // const unsigned int filterIndex35(triggerEventWithRefsHandle_->filterIndex(InputTag("hltL1sDoubleTau35", "", "MYHLT")));
    // //making tau object and filling vector
    // triggerEventWithRefsHandle_->getObjects(filterIndex35, tvids, tauCandRefVec);

    // const unsigned int nTaus35(tvids.size());
    // hltL1DoubleTau35_nTaus = nTaus35;
    // if (nTaus > 0) {
        // for (unsigned int i = 0; i != nTaus35; ++i) {
            // hltL1DoubleTau35_tauPt.push_back(tauCandRefVec[i]->pt());
            // hltL1DoubleTau35_tauEta.push_back(tauCandRefVec[i]->eta());
            // hltL1DoubleTau35_tauPhi.push_back(tauCandRefVec[i]->phi());
            // hltL1DoubleTau35_tauEnergy.push_back(tauCandRefVec[i]->energy());
        // }
    // }    

    //filling the rest, as well as passFilter branches


    const trigger::size_type nFilters(triggerEvent->sizeFilters());
    // const trigger::size_type nObjects(triggerEvent->sizeObjects());
    // const trigger::size_type nCollections(triggerEvent->sizeCollections());
    // std::cout << nFilters << " " << nObjects << " " << nCollections << " " << std::endl;
    
    // std::string hltL1sDoubleTauBigOR_Tag = "hltL1sDoubleTauBigOR::MYHLT"; // ditau L1
    // std::string hltL1VBFDiJetOR_Tag = "hltL1VBFDiJetOR::MYHLT";		  // inclusive L1
    // std::string hltL1VBFDiJetIsoTau_Tag = "hltL1VBFDiJetIsoTau::MYHLT";	  // proposed L1
    std::string hltL1DoubleTauJet_Tag = "hltL1sDoubleTauJet::MYHLT";	  // proposed L1
    std::string hltL1DoubleTau32_Tag = "hltL1sDoubleTauBigOR::MYHLT";	  // ditau32 L1
    std::string hltL1DoubleTau35_Tag = "hltL1sDoubleTau35::MYHLT";	  // ditau35 L1
    
    std::string hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_Tag = "hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet::MYHLT";
    std::string hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_Tag = "hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet::MYHLT";
    std::string hltHpsPFTauTrack_Tag = "hltHpsPFTauTrack::MYHLT";
    std::string hltL2DoubleTauTagNNFilterDoubleTauJet_Tag = "hltL2DoubleTauTagNNFilterDoubleTauJet::MYHLT";
    std::string hltPFJets60L1HLTMatched_Tag = "hltPFJets60L1HLTMatched::MYHLT";
    std::string hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_Tag = "hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60::MYHLT";
    
    // accepted filters per event
    for(trigger::size_type iFilter=0; iFilter!=nFilters; ++iFilter) {
        std::string filterTag = triggerEvent->filterTag(iFilter).encode();
        // std::cout << filterTag << std::endl;
        trigger::Keys objectKeys = triggerEvent->filterKeys(iFilter);
        
        const trigger::TriggerObjectCollection& triggerObjects(triggerEvent->getObjects());
        // fill "pass filter" branches
        int nObjKeys = objectKeys.size();
        // if (filterTag == hltL1VBFDiJetOR_Tag && nObjKeys >= 0) nEvents = 1; // accept pass or fail condition to fill nEvents
        if (filterTag == hltL1DoubleTauJet_Tag && nObjKeys >= 0) nEvents = 1; // accept pass or fail condition to fill nEvents
        
        // L1s
        // if (filterTag == hltL1sDoubleTauBigOR_Tag && nObjKeys >= 2) passhltL1sDoubleTauBigOR = 1;
        // if (filterTag == hltL1VBFDiJetOR_Tag && nObjKeys >= 2) passhltL1VBFDiJetOR = 1;
        // if (filterTag == hltL1VBFDiJetIsoTau_Tag && hltL1VBFDiJetIsoTau_tauPt.size() >= 1
        // && hltL1VBFDiJetIsoTau_jetPt.size() >= 2) passhltL1VBFDiJetIsoTau = 1;
        if (filterTag == hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_Tag && nObjKeys >= 2) {
            passhltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet = 1;
            hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_nObjs = nObjKeys;
        }
        if (filterTag == hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_Tag && nObjKeys >= 2) {
            passhltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet = 1;
            hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_nObjs = nObjKeys;
        }
        if (filterTag == hltHpsPFTauTrack_Tag && nObjKeys >= 2) {
            passhltHpsPFTauTrack = 1;
            hltHpsPFTauTrack_nObjs = nObjKeys;
        }
        if (filterTag == hltL2DoubleTauTagNNFilterDoubleTauJet_Tag && nObjKeys >= 2) {
            passhltL2DoubleTauTagNNFilterDoubleTauJet = 1;
            hltL2DoubleTauTagNNFilterDoubleTauJet_nObjs = nObjKeys;
        }
        if (filterTag == hltPFJets60L1HLTMatched_Tag && nObjKeys >= 2) {
            passhltPFJets60L1HLTMatched = 1;
            hltPFJets60L1HLTMatched_nObjs = nObjKeys;
        }
        if (filterTag == hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_Tag && nObjKeys >= 2) {
            passhltHpsOverlapFilterDeepTauDoublePFTau30PFJet60 = 1;
            hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_nObjs = nObjKeys;
        }


        //loop over trigger objects and store their kinematics to the proper filter branches
        for(trigger::size_type iKey=0; iKey < nObjKeys; ++iKey){
            trigger::size_type objKey = objectKeys.at(iKey);
            const trigger::TriggerObject& triggerObj(triggerObjects[objKey]);
            pt_ = triggerObj.pt();
            eta_ = triggerObj.eta();
            phi_ = triggerObj.phi();
            energy_ = triggerObj.energy();

            if (filterTag == hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_Tag
                    && passhltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet && pt_>0) {
                hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_pt.push_back(pt_);
                hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_eta.push_back(eta_);
                hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_phi.push_back(phi_);
                hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet_energy.push_back(energy_);
            }
            if (filterTag == hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_Tag
                    && passhltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet && pt_>0) {
                hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_pt.push_back(pt_);
                hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_eta.push_back(eta_);
                hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_phi.push_back(phi_);
                hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet_energy.push_back(energy_);
            }
            if (filterTag == hltHpsPFTauTrack_Tag
                    && passhltHpsPFTauTrack && pt_>0) {
                hltHpsPFTauTrack_pt.push_back(pt_);
                hltHpsPFTauTrack_eta.push_back(eta_);
                hltHpsPFTauTrack_phi.push_back(phi_);
                hltHpsPFTauTrack_energy.push_back(energy_);
            }
            if (filterTag == hltL2DoubleTauTagNNFilterDoubleTauJet_Tag
                    && passhltL2DoubleTauTagNNFilterDoubleTauJet && pt_>0) {
                hltL2DoubleTauTagNNFilterDoubleTauJet_pt.push_back(pt_);
                hltL2DoubleTauTagNNFilterDoubleTauJet_eta.push_back(eta_);
                hltL2DoubleTauTagNNFilterDoubleTauJet_phi.push_back(phi_);
                hltL2DoubleTauTagNNFilterDoubleTauJet_energy.push_back(energy_);
            }
            if (filterTag == hltPFJets60L1HLTMatched_Tag
                    && passhltPFJets60L1HLTMatched && pt_>0) {
                hltPFJets60L1HLTMatched_pt.push_back(pt_);
                hltPFJets60L1HLTMatched_eta.push_back(eta_);
                hltPFJets60L1HLTMatched_phi.push_back(phi_);
                hltPFJets60L1HLTMatched_energy.push_back(energy_);
            }
            if (filterTag == hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_Tag
                    && passhltHpsOverlapFilterDeepTauDoublePFTau30PFJet60 && pt_>0) {
                hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_pt.push_back(pt_);
                hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_eta.push_back(eta_);
                hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_phi.push_back(phi_);
                hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60_energy.push_back(energy_);
            }

        } // end loop over trigger object keys
    } // end loop over nfilters
} // end function
