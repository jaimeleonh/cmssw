// -*- C++ -*-
//
// Package:    RecoTauTag/HLTProducers
// Class:      L1HLTTauJetMatching
//
/**\class L1HLTTauJetMatching L1HLTTauJetMatching.h 
 RecoTauTag/HLTProducers/interface/L1HLTTauJetMatching.h
 Description: 
 Matching L1 to PF/Calo Jets. Used for HLT_VBF paths.
	*Matches PF/Calo Jets to L1 jets from the dedicated seed
	*Adds selection criteria to the leading/subleading jets as well as the maximum dijet mass
	*Separates collections of PF/Calo jets into two categories
 
 
*/
//
// Original Author:  Jaime Leon Holgado (CIEMAT)
//         Created:  Tue, 28 Sep 2021 13:00:00 GMT
//
//

#ifndef RecoTauTag_HLTProducers_L1HLTTauJetsMatching_h
#define RecoTauTag_HLTProducers_L1HLTTauJetsMatching_h

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TauReco/interface/PFTauFwd.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"

#include "Math/GenVector/VectorUtil.h"
#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/TauReco/interface/PFTau.h"

#include "HLTrigger/HLTcore/interface/defaultModuleLabel.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <map>
#include <vector>

class L1HLTTauJetMatching : public edm::EDProducer {
public:
  explicit L1HLTTauJetMatching(const edm::ParameterSet&);
  ~L1HLTTauJetMatching() override;
  void produce(edm::Event&, const edm::EventSetup&) override;
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  const edm::EDGetTokenT<reco::PFTauCollection> tauSrc_;
  const edm::EDGetTokenT<reco::PFJetCollection> jetSrc_;
  const edm::EDGetTokenT<trigger::TriggerFilterObjectWithRefs> L1JetSrc_;
  const double minTauPt_;
  const double minJetPt_;
  const double matchingL1HLTR_;
  const double matchingTauJetR_;
  const double matchingL1HLTR2_;
  const double matchingTauJetR2_;
};


#endif
