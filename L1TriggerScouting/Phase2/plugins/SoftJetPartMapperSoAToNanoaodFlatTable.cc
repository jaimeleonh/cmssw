#include "FWCore/Framework/interface/MakerMacros.h"

#include <fstream>
#include <iomanip>
#include <memory>
#include <string>
#include <cmath>

#include "FWCore/Framework/interface/global/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDGetToken.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/MessageLogger/interface/MessageDrop.h"

#include "DataFormats/NanoAOD/interface/FlatTable.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"


class SoftJetPartMapperSoAToNanoaodFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit SoftJetPartMapperSoAToNanoaodFlatTable(const edm::ParameterSet&);
  ~SoftJetPartMapperSoAToNanoaodFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::SoftJetOutputHostTensor> srcSoftJetParT_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

SoftJetPartMapperSoAToNanoaodFlatTable::SoftJetPartMapperSoAToNanoaodFlatTable(const edm::ParameterSet& iConfig)
    : srcSoftJetParT_(consumes<l1sc::SoftJetOutputHostTensor>(iConfig.getParameter<edm::InputTag>("src"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")){
  produces<nanoaod::FlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void SoftJetPartMapperSoAToNanoaodFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::SoftJetOutputHostTensor> srcSoftJetParT;
  iEvent.getByToken(srcSoftJetParT_, srcSoftJetParT);

  const auto *SoftJetParT = srcSoftJetParT->const_view().output().data();
  const unsigned int nclusters = srcSoftJetParT->const_view().metadata().size();
  std::vector<float> softJetPartScores{SoftJetParT, SoftJetParT + nclusters};
  // std::vector<int32_t> is_seed{seed, seed + nclusters};

  auto out = std::make_unique<nanoaod::FlatTable>(nclusters, name_, false, true);
  out->setDoc(doc_);
  out->addColumn<float>("softJetParTScore", softJetPartScores, "SoftJet ParT score");
  iEvent.put(std::move(out));
}

void SoftJetPartMapperSoAToNanoaodFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(SoftJetPartMapperSoAToNanoaodFlatTable);
