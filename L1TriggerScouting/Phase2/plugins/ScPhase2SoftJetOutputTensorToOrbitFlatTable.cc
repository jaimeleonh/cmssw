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

#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHost.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"

class ScPhase2SoftJetOutputTensorToOrbitFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ScPhase2SoftJetOutputTensorToOrbitFlatTable(const edm::ParameterSet&);
  ~ScPhase2SoftJetOutputTensorToOrbitFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::BxLookupHost> srcBxClustersMap_;
  edm::EDGetTokenT<l1sc::SoftJetOutputHostTensor> srcOutput_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ScPhase2SoftJetOutputTensorToOrbitFlatTable::ScPhase2SoftJetOutputTensorToOrbitFlatTable(const edm::ParameterSet& iConfig) :
      srcBxClustersMap_(consumes<l1sc::BxLookupHost>(iConfig.getParameter<edm::InputTag>("srcBxClustersMap"))),
      srcOutput_(consumes<l1sc::SoftJetOutputHostTensor>(iConfig.getParameter<edm::InputTag>("srcOutput"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<l1ScoutingRun3::OrbitFlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ScPhase2SoftJetOutputTensorToOrbitFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::BxLookupHost> srcBxClustersMap;
  iEvent.getByToken(srcBxClustersMap_, srcBxClustersMap);
  edm::Handle<l1sc::SoftJetOutputHostTensor> srcOutput;
  iEvent.getByToken(srcOutput_, srcOutput);

  const auto nbx = srcBxClustersMap->const_view().bx().metadata().size();

  // bx->clusters map
  const auto *bxc_offsets = srcBxClustersMap->const_view().offset().offset().data();

  // soft Jet output tensor
  const auto *cls = srcOutput->const_view().output().data();
  const auto nclusters = srcOutput->const_view().metadata().size();

  // offsets for dividing model outputs for each event
  std::vector<unsigned int> bxOffsets;
  bxOffsets.push_back(0);
  bxOffsets.insert(bxOffsets.end(), bxc_offsets, bxc_offsets + nbx + 1);

  // table with output logits
  std::vector<float> cls_vec{cls, cls + nclusters};

  auto out_table = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, name_);
  out_table->setDoc(doc_);
  out_table->addColumn<float>("softJetPart", cls_vec, "Soft Jet ParT output score");
  iEvent.put(std::move(out_table));
}

void ScPhase2SoftJetOutputTensorToOrbitFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcBxClustersMap");
  desc.add<edm::InputTag>("srcOutput");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2SoftJetOutputTensorToOrbitFlatTable);
