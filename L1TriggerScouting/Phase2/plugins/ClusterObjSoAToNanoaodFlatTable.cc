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
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"


class ClusterObjSoAToNanoaodFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ClusterObjSoAToNanoaodFlatTable(const edm::ParameterSet&);
  ~ClusterObjSoAToNanoaodFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::ClusterObjHostCollection> srcClusters_;

  std::string name_, doc_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ClusterObjSoAToNanoaodFlatTable::ClusterObjSoAToNanoaodFlatTable(const edm::ParameterSet& iConfig)
    : srcClusters_(consumes<l1sc::ClusterObjHostCollection>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")) {
  produces<nanoaod::FlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ClusterObjSoAToNanoaodFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::ClusterObjHostCollection> srcClusters;
  iEvent.getByToken(srcClusters_, srcClusters);

  const unsigned int nclusters = srcClusters->const_view().metadata().size();
  const auto *pt = srcClusters->const_view().pt().data();
  const auto *eta = srcClusters->const_view().eta().data();
  const auto *phi = srcClusters->const_view().phi().data();
  const auto *cluster = srcClusters->const_view().cluster().data();
  const auto *ndau = srcClusters->const_view().numberOfDaughters().data();

  std::vector<float> pts, etas, phis;
  std::vector<int32_t> clusters, ndaus;

  for (unsigned int i = 0; i < nclusters; i++) {
    if (pt[i] > 0) {
      pts.push_back(pt[i]);
      etas.push_back(eta[i]);
      phis.push_back(phi[i]);
      clusters.push_back(cluster[i]);
      ndaus.push_back(ndau[i]);
    }
  }

  auto out = std::make_unique<nanoaod::FlatTable>(pts.size(), name_, false, false);
  out->setDoc(doc_);
  out->addColumn<float>("pt", pts, "cluster pt");
  out->addColumn<float>("eta", etas, "cluster eta");
  out->addColumn<float>("phi", phis, "cluster phi");
  out->addColumn<int32_t>("cluster", clusters, "cluster index");
  out->addColumn<int32_t>("ndau", ndaus, "number of daughters");
  iEvent.put(std::move(out));
}

void ClusterObjSoAToNanoaodFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcClusters");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ClusterObjSoAToNanoaodFlatTable);
