#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/Math/interface/deltaR.h"

#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"
#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"

#include "fastjet/ClusterSequence.hh"

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiAKJetsSoA : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiAKJetsSoA(const edm::ParameterSet &);
  ~ScPhase2PuppiAKJetsSoA() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;

  edm::EDGetTokenT<l1sc::PuppiHostCollection> src_candidates_token_;
  edm::EDGetTokenT<l1sc::BxLookupHostCollection> bx_lookup_token_;
  double R_;
};

ScPhase2PuppiAKJetsSoA::ScPhase2PuppiAKJetsSoA(const edm::ParameterSet &iConfig)
    : src_candidates_token_(consumes<l1sc::PuppiHostCollection>(iConfig.getParameter<edm::InputTag>("src"))),
    bx_lookup_token_(consumes<l1sc::BxLookupHostCollection>(iConfig.getParameter<edm::InputTag>("src"))),
      R_(iConfig.getParameter<double>("rParam")) {
    produces<l1sc::AssociationMapHost>();
    produces<l1sc::BxLookupHostCollection>();
    produces<l1sc::ClustersHostCollection>();
    produces<l1sc::ClusterObjHostCollection>();
}

ScPhase2PuppiAKJetsSoA::~ScPhase2PuppiAKJetsSoA() {};

void ScPhase2PuppiAKJetsSoA::beginStream(edm::StreamID) {}

void ScPhase2PuppiAKJetsSoA::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  using namespace fastjet;

  edm::Handle<l1sc::PuppiHostCollection> src;
  iEvent.getByToken(src_candidates_token_, src);
  l1sc::PuppiHostCollection::ConstView candidates = src->const_view();

   edm::Handle<l1sc::BxLookupHostCollection> bxLookup;
  iEvent.getByToken(bx_lookup_token_, bxLookup); 
  unsigned int nbx = bxLookup->view<l1sc::OffsetsSoA>().metadata().size() - 1;
  auto srcOffsets = bxLookup->const_view<l1sc::OffsetsSoA>().offsets();
  auto srcBx = bxLookup->const_view<l1sc::BxIndexSoA>().bx();

  std::array<int32_t, 2> bxLookupSizes{{int32_t(nbx), int32_t(nbx+1)}};
  auto jetBxLookup = std::make_unique<l1sc::BxLookupHostCollection>(bxLookupSizes, cms::alpakatools::host());
  auto jetOffsets = jetBxLookup->view<l1sc::OffsetsSoA>().offsets();
  auto jetBx = jetBxLookup->view<l1sc::BxIndexSoA>().bx();

  // containers for output products
  // jet 3 momenta
  std::vector<float> pt;
  std::vector<float> eta;
  std::vector<float> phi;
  // clustering information
  auto clusterInfo = std::make_unique<l1sc::ClustersHostCollection>(candidates.metadata().size(), cms::alpakatools::host());
  auto cluster = clusterInfo->view().cluster();
  auto is_seed = clusterInfo->view().is_seed();

  // choose a jet definition
  JetDefinition jet_def(antikt_algorithm, R_);
  std::vector<PseudoJet> particles;
  std::vector<unsigned int> ndaughters;
  std::vector<unsigned int> indices;
  indices.reserve(candidates.metadata().size());

  for (unsigned int block_idx = 0; block_idx < nbx; ++block_idx) {
    jetOffsets[block_idx] = pt.size();
    jetBx[block_idx] = srcBx[block_idx];
    uint32_t begin = srcOffsets[block_idx];
    uint32_t end = srcOffsets[block_idx + 1];
    if (end <= begin)
        continue;
    uint32_t block_dim = end - begin;

    particles.clear();
    for (unsigned i = 0; i < block_dim; i++) {
        unsigned int idx = begin + i;
        cluster[idx] = -1;
        is_seed[idx] = 0;
        float mass = 0.13;
        ROOT::Math::PtEtaPhiMVector p4(candidates.pt()[idx], candidates.eta()[idx], candidates.phi()[idx], mass); 
        particles.emplace_back(p4.px(), p4.py(), p4.pz(), p4.energy());
        particles.back().set_user_index(i);
    }

    // run the clustering, extract the jets
    ClusterSequence cs(particles, jet_def);
    auto jets = sorted_by_pt(cs.inclusive_jets());

    unsigned int icluster = 0;
    for (auto j : jets) {
      pt.push_back(j.pt());
      eta.push_back(j.eta());
      phi.push_back(reco::reducePhiRange(j.phi()));
      const auto & constituents = j.constituents();
      ndaughters.push_back(constituents.size());
      for (const auto & dau : constituents) {
        unsigned int idx = dau.user_index() + begin;
        cluster[idx] = icluster;
        indices.push_back(idx);
      } 
      icluster++;
    }

  }  // loop over BX
  // finalize last bx offset
  jetOffsets[nbx] = pt.size();

  int32_t nclustered = indices.size();
  int32_t njets = pt.size();

  std::array<int32_t, 2> clusterInfoSizes{{nclustered, njets+1}};
  auto map = std::make_unique<l1sc::AssociationMapHost>(clusterInfoSizes, cms::alpakatools::host());
  auto jets = std::make_unique<l1sc::ClusterObjHostCollection>(njets, cms::alpakatools::host());
  auto jetsOut = jets->view();
  auto mapIndices = map->view<l1sc::IndexSoA>().indexes();
  for (int32_t i = 0; i < nclustered; ++i) {
    mapIndices[i] = indices[i];
  }
  auto mapOffsets = map->view<l1sc::OffsetsSoA>().offsets();
  mapOffsets[0] = 0;
  for (int32_t i = 0; i < njets; ++i) {
    mapOffsets[i + 1] = mapOffsets[i] + ndaughters[i];
    jetsOut.pt()[i] = pt[i];
    jetsOut.eta()[i] = eta[i];
    jetsOut.phi()[i] = phi[i];
    jetsOut.numberOfDaughters()[i] = ndaughters[i];
    jetsOut.cluster()[i] = 0;
  }

  iEvent.put(std::move(clusterInfo));
  iEvent.put(std::move(jetBxLookup));
  iEvent.put(std::move(jets));
  iEvent.put(std::move(map));
}

void ScPhase2PuppiAKJetsSoA::endStream() {}

void ScPhase2PuppiAKJetsSoA::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<double>("rParam", 0.4);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiAKJetsSoA);