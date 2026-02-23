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
#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"


class ClusterMapperSoAToNanoaodFlatTable : public edm::global::EDProducer<> {
public:
  // constructor and destructor
  explicit ClusterMapperSoAToNanoaodFlatTable(const edm::ParameterSet&);
  ~ClusterMapperSoAToNanoaodFlatTable() override {};

  void produce(edm::StreamID, edm::Event&, edm::EventSetup const&) const override;

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  // the tokens to access the data
  edm::EDGetTokenT<l1sc::ClustersHostCollection> srcClusters_;
  edm::EDGetTokenT<l1sc::AssociationMapHost> srcMapCands_;

  std::string name_, doc_, varName_;
};
// -----------------------------------------------------------------------------

// -------------------------------- constructor  -------------------------------

ClusterMapperSoAToNanoaodFlatTable::ClusterMapperSoAToNanoaodFlatTable(const edm::ParameterSet& iConfig)
    : srcClusters_(consumes<l1sc::ClustersHostCollection>(iConfig.getParameter<edm::InputTag>("srcClusters"))),
      srcMapCands_(consumes<l1sc::AssociationMapHost>(iConfig.getParameter<edm::InputTag>("srcMapCands"))),
      name_(iConfig.getParameter<std::string>("name")),
      doc_(iConfig.getParameter<std::string>("doc")),
      varName_(iConfig.getParameter<std::string>("varName")){
  produces<nanoaod::FlatTable>();
}
// -----------------------------------------------------------------------------

// ----------------------- method called for each orbit  -----------------------
void ClusterMapperSoAToNanoaodFlatTable::produce(edm::StreamID, edm::Event& iEvent, edm::EventSetup const&) const {
  edm::Handle<l1sc::ClustersHostCollection> srcClusters;
  iEvent.getByToken(srcClusters_, srcClusters);
  edm::Handle<l1sc::AssociationMapHost> srcMapCands;
  iEvent.getByToken(srcMapCands_, srcMapCands);

  const unsigned int nObj = srcClusters->const_view().metadata().size();
  const unsigned int nClusters = srcMapCands->const_view<l1sc::OffsetsSoA>().metadata().size();
  const auto *cand_indexes = srcMapCands->const_view<l1sc::IndexSoA>().indexes().data();
  const auto *cand_offsets = srcMapCands->const_view<l1sc::OffsetsSoA>().offsets().data();

  std::vector<int32_t> cluster_indexes;

  std::map<size_t, size_t> candidate_cluster_map;
  for (size_t icluster = 0; icluster < nClusters; icluster++) {
    for (size_t icand = cand_offsets[icluster]; icand < cand_offsets[icluster + 1]; icand++) {
      candidate_cluster_map[cand_indexes[icand]] = icluster;
    }
  }

  for (size_t i = 0; i < nObj; i++) {
    if (candidate_cluster_map.find(i) != candidate_cluster_map.end()) {
      cluster_indexes.push_back(static_cast<int32_t>(candidate_cluster_map[i]));
    } else {
      cluster_indexes.push_back(-1);
    }
  }

  auto out = std::make_unique<nanoaod::FlatTable>(cluster_indexes.size(), name_, false, true);
  out->setDoc(doc_);
  out->addColumn<int32_t>(varName_, cluster_indexes, "associated cluster");
  iEvent.put(std::move(out));
}

void ClusterMapperSoAToNanoaodFlatTable::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcClusters");
  desc.add<edm::InputTag>("srcMapCands");
  desc.add<std::string>("name");
  desc.add<std::string>("doc");
  desc.add<std::string>("varName");
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ClusterMapperSoAToNanoaodFlatTable);
