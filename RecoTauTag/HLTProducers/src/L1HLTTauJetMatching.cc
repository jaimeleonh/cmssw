#include "RecoTauTag/HLTProducers/interface/L1HLTTauJetMatching.h"

//
// class declaration
//

L1HLTTauJetMatching::L1HLTTauJetMatching(const edm::ParameterSet& iConfig)
    : tauSrc_(consumes<reco::PFTauCollection>(iConfig.getParameter<edm::InputTag>("TauSrc"))),
      jetSrc_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("JetSrc"))),
      L1JetSrc_(consumes<trigger::TriggerFilterObjectWithRefs>(iConfig.getParameter<edm::InputTag>("L1JetSrc"))),
      minTauPt_(iConfig.getParameter<double>("minTauPt")),
      minJetPt_(iConfig.getParameter<double>("minJetPt")),
      matchingL1HLTR_(iConfig.getParameter<double>("matchingL1HLTR")),
      matchingTauJetR_(iConfig.getParameter<double>("matchingTauJetR")),
      matchingL1HLTR2_(matchingL1HLTR_ * matchingL1HLTR_),
      matchingTauJetR2_(matchingTauJetR_ * matchingTauJetR_) {
  produces<reco::PFTauCollection>("taus");
  produces<reco::PFJetCollection>("jets");
}

L1HLTTauJetMatching::~L1HLTTauJetMatching() {}

void L1HLTTauJetMatching::produce(edm::Event& iEvent, const edm::EventSetup& iES) {
  std::unique_ptr<reco::PFTauCollection> L1TmatchedPFTau(new reco::PFTauCollection);
  std::unique_ptr<reco::PFJetCollection> L1TmatchedPFJet(new reco::PFJetCollection);

  edm::Handle<reco::PFTauCollection> taus;
  iEvent.getByToken(tauSrc_, taus);

  edm::Handle<reco::PFJetCollection> jets;
  iEvent.getByToken(jetSrc_, jets);

  edm::Handle<trigger::TriggerFilterObjectWithRefs> L1Jets;
  iEvent.getByToken(L1JetSrc_, L1Jets);

  l1t::JetVectorRef jetCandRefVec;
  L1Jets->getObjects(trigger::TriggerL1Jet, jetCandRefVec);

  /* Loop over taus that must pass a certain minTauPt_ cut */
  /* then loop over L1T jets and check whether they match, */
  /* if yes -> include the 2 highest pt taus in */
  /* the new L1T matched PFTau collection */

  std::vector<int> iMatchedTaus = {-1, -1};
  std::vector<double> ptMatchedTaus = {-1., -1.};

  for (unsigned int iTau = 0; iTau < taus->size(); iTau++) {
    bool isMatched = false;
    if ((*taus)[iTau].pt() > minTauPt_) {
      for (unsigned int iL1Jet = 0; iL1Jet < jetCandRefVec.size(); iL1Jet++) {
        if (reco::deltaR2((*taus)[iTau].p4(), jetCandRefVec[iL1Jet]->p4()) < matchingL1HLTR2_) {
          isMatched = true;
          break;
        }
      }
    }
    if (isMatched) {
      if ((*taus)[iTau].pt() > ptMatchedTaus[0]) {
        ptMatchedTaus[1] = ptMatchedTaus[0];
        iMatchedTaus[1] = iMatchedTaus[0];
        ptMatchedTaus[0] = (*taus)[iTau].pt() ;
        iMatchedTaus[0] = iTau;
      } else if ((*taus)[iTau].pt() > ptMatchedTaus[1]) {
        ptMatchedTaus[1] = (*taus)[iTau].pt() ;
        iMatchedTaus[1] = iTau;
      }
    }
  }
  
  if (iMatchedTaus[0] != -1)
    L1TmatchedPFTau->push_back((*taus)[iMatchedTaus[0]]);
  if (iMatchedTaus[1] != -1)
    L1TmatchedPFTau->push_back((*taus)[iMatchedTaus[1]]);

  /* Loop over jets that must pass a certain minJetPt_ cut */
  /* then loop over L1T jets and check whether they match, */
  /* if yes -> check if they match with the two taus stored previously*/
  /* if not -> include the highest pt jet in the new L1T matched PFJet collection */
     
  int iMatchedJet = -1;
  double ptMatchedJet = -1.;
  for (unsigned int iJet = 0; iJet < jets->size(); iJet++) {
    bool isMatched = false;
    if ((*jets)[iJet].pt() > minJetPt_) {
      for (unsigned int iL1Jet = 0; iL1Jet < jetCandRefVec.size(); iL1Jet++) {
        if (reco::deltaR2((*jets)[iJet].p4(), jetCandRefVec[iL1Jet]->p4()) < matchingL1HLTR2_) {
          for (unsigned int iMatchedTau = 0; iMatchedTau < L1TmatchedPFTau->size(); iMatchedTau++) {
            if (reco::deltaR2((*jets)[iJet].p4(), (L1TmatchedPFTau->at(iMatchedTau)).p4()) > matchingTauJetR2_) {
              isMatched = true;
              break;
            }
          }
        }
      }
    }
    if (isMatched) {
      if ((*jets)[iJet].pt() > ptMatchedJet) {
        ptMatchedJet = (*jets)[iJet].pt();
        iMatchedJet = iJet;
      }
    }
  }
  if (iMatchedJet != -1)
    L1TmatchedPFJet->push_back((*jets)[iMatchedJet]);

  iEvent.put(std::move(L1TmatchedPFTau), "taus");
  iEvent.put(std::move(L1TmatchedPFJet), "jets");
}

void L1HLTTauJetMatching::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("L1JetSrc", edm::InputTag("whatever"))
      ->setComment("Input filter objects passing L1 seed");
  desc.add<edm::InputTag>("TauSrc", edm::InputTag("whatever"))
      ->setComment("Input collection of PFTaus");
  desc.add<edm::InputTag>("JetSrc", edm::InputTag("whatever"))
      ->setComment("Input collection of PFJets");
  desc.add<double>("minTauPt", 40.0)->setComment("Minimal pT1 of PFTaus to match");
  desc.add<double>("minJetPt", 55.0)->setComment("Minimal pT2 of PFJets to match");
  desc.add<double>("matchingL1HLTR", 0.5)->setComment("dR value used for matching between HLT and L1 objects");
  desc.add<double>("matchingTauJetR", 0.5)->setComment("dR value used for matching between PFTaus and PFJets");
  descriptions.setComment(
    "This module produces a collection of PFTaus and a collection of PFJets matched to the L1 jets.");
  descriptions.add("L1HLTTauJetMatching", desc);
}
//define this as a plug-in
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(L1HLTTauJetMatching);
