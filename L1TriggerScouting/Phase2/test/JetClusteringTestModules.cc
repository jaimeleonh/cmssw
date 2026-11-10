#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <vector>

#include "DataFormats/L1ScoutingSoA/interface/PuppiHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHost.h"
#include "DataFormats/L1ScoutingSoA/interface/ClustersHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/Exception.h"

namespace {
  constexpr std::array<uint16_t, 4> bxIds{1, 7, 55, 99};
  void check(bool condition, char const* message) {
    if (!condition)
      throw cms::Exception("JetClusteringTest") << message;
  }
}

// The serial Alpaka backend consumes these same host product types. The second
// event is empty; the first exercises phi wrapping, tied pT and the >128 fallback.
class Phase2JetTestInput : public edm::stream::EDProducer<> {
public:
  explicit Phase2JetTestInput(edm::ParameterSet const&)
      : particles_{produces<l1sc::PuppiHostCollection>()}, bx_{produces<l1sc::BxLookupHost>()} {}

  void produce(edm::Event& event, edm::EventSetup const&) override {
    bool empty = event.id().event() == 2;
    auto particles = l1sc::PuppiHostCollection(empty ? 0 : 137);
    auto bx = l1sc::BxLookupHost(4, 5);
    std::array<uint32_t, 5> offsets{0, 4, 4, 136, 137};
    for (unsigned int i = 0; i < bxIds.size(); ++i)
      bx.view().bx()[i].bx() = bxIds[i];
    for (unsigned int i = 0; i < offsets.size(); ++i)
      bx.view().offset()[i].offset() = empty ? 0 : offsets[i];
    auto view = particles.view();
    auto set = [&](int i, float pt, float eta, float phi) {
      view.pt()[i] = pt;
      view.eta()[i] = eta;
      view.phi()[i] = phi;
      view.z0()[i] = view.dxy()[i] = 0.f;
      view.puppiw()[i] = 1.f;
      view.quality()[i] = 0;
      view.pdgid()[i] = 211;
    };
    if (!empty) {
      set(0, 40.f, 0.f, 3.12f);
      set(1, 20.f, 0.f, -3.12f);
      set(2, 15.f, 2.f, 0.f);
      set(3, 5.f, 2.03f, .01f);
      for (int i = 4; i < 136; ++i)
        set(i, 1.f, .5f + .0001f * (i % 4), .5f);
      set(136, 8.f, -2.f, .4f);
    }
    event.emplace(particles_, std::move(particles));
    event.emplace(bx_, std::move(bx));
  }

private:
  edm::EDPutTokenT<l1sc::PuppiHostCollection> particles_;
  edm::EDPutTokenT<l1sc::BxLookupHost> bx_;
};

class Phase2JetTestCheck : public edm::one::EDAnalyzer<> {
public:
  explicit Phase2JetTestCheck(edm::ParameterSet const& p)
      : particles_{consumes<l1sc::PuppiHostCollection>(edm::InputTag("input"))},
        bx_{consumes<l1sc::BxLookupHost>(edm::InputTag("input"))},
        jets_{consumes<l1sc::ClusterObjHostCollection>(p.getParameter<edm::InputTag>("src"))},
        jetBx_{consumes<l1sc::BxLookupHost>(p.getParameter<edm::InputTag>("src"))},
        clusters_{consumes<l1sc::ClustersHostCollection>(p.getParameter<edm::InputTag>("src"))},
        map_{consumes<l1sc::AssociationMapHost>(p.getParameter<edm::InputTag>("src"))},
        minConstituents_{p.getParameter<unsigned int>("minConstituents")},
        capped_{p.getParameter<bool>("capped")} {}

  void analyze(edm::Event const& event, edm::EventSetup const&) override {
    auto particles = event.get(particles_).const_view();
    auto bx = event.get(bx_).const_view();
    auto jets = event.get(jets_).const_view();
    auto jetBx = event.get(jetBx_).const_view();
    auto clusters = event.get(clusters_).const_view();
    auto map = event.get(map_).const_view();
    bool empty = event.id().event() == 2;
    std::array<unsigned int, 4> counts{capped_ ? 1u : 2u, 0, 1, 1};
    if (minConstituents_ >= 2)
      counts[3] = 0;
    if (minConstituents_ >= 3)
      counts[0] = 0;
    if (empty)
      counts.fill(0);
    check(jetBx.offset()[0].offset() == 0, "Jet offsets must start at zero");
    check(map.offset()[0].offset() == 0, "Constituent offsets must start at zero");
    unsigned int totalJets = 0, totalConstituents = 0;
    std::vector<bool> seen(particles.metadata().size(), false);
    for (unsigned int b = 0; b < bxIds.size(); ++b) {
      check(jetBx.bx()[b].bx() == bxIds[b], "BX IDs were not preserved");
      auto first = jetBx.offset()[b].offset(), end = jetBx.offset()[b + 1].offset();
      check(first == totalJets && end - first == counts[b], "Unexpected number of jets in BX");
      for (unsigned int j = first; j < end; ++j) {
        check(jets.cluster()[j] == int(j - first), "Jet IDs are not contiguous within BX");
        auto from = map.offset()[j].offset(), to = map.offset()[j + 1].offset();
        check(to >= from && to - from == jets.numberOfDaughters()[j], "Constituent count mismatch");
        check(to <= unsigned(map.index().metadata().size()), "Constituent offset outside map");
        check(to > from, "Empty jet");
        float pt = 0.f, eta = 0.f, dphi = 0.f;
        unsigned int previous = 0;
        auto ref = map.index()[from].index();
        for (unsigned int k = from; k < to; ++k) {
          auto i = map.index()[k].index();
          check(i >= bx.offset()[b].offset() && i < bx.offset()[b + 1].offset(), "Cross-BX association");
          check(!seen[i], "Candidate assigned more than once");
          seen[i] = true;
          check(clusters.cluster()[i] == int(j - first), "Cluster label and association disagree");
          check(clusters.is_seed()[i] == 0 || clusters.is_seed()[i] == 1, "Invalid seed flag");
          if (k > from)
            check(particles.pt()[previous] > particles.pt()[i] ||
                      (particles.pt()[previous] == particles.pt()[i] && previous < i),
                  "Constituents must be ordered by descending pT and input index for ties");
          pt += particles.pt()[i];
          eta += particles.pt()[i] * particles.eta()[i];
          dphi += particles.pt()[i] * std::remainder(particles.phi()[i] - particles.phi()[ref], 2.f * float(M_PI));
          previous = i;
        }
        check(std::abs(pt - jets.pt()[j]) < .001f, "Jet pT is not the sum of constituent pT");
        check(std::abs(eta / pt - jets.eta()[j]) < .001f, "Incorrect jet eta centroid");
        check(std::abs(std::remainder(particles.phi()[ref] + dphi / pt - jets.phi()[j], 2.f * float(M_PI))) < .001f,
              "Incorrect jet phi centroid across branch cut");
        check(to - from >= minConstituents_, "Constituent filter was not applied");
        if (b == 0)
          check(std::abs(pt - 60.f) < .001f || (!capped_ && std::abs(pt - 20.f) < .001f), "Incorrect small-BX jet pT");
        if (b == 2)
          check(std::abs(pt - 132.f) < .001f, "Incorrect large-BX jet pT");
        totalConstituents += to - from;
      }
      totalJets = end;
    }
    check(totalJets == unsigned(jets.metadata().size()), "Jet storage size mismatch");
    check(totalConstituents == unsigned(map.index().metadata().size()), "Association map storage size mismatch");
    for (int i = 0; i < particles.metadata().size(); ++i)
      check(seen[i] == (clusters.cluster()[i] >= 0), "Unassigned candidate or missing association");
  }

private:
  edm::EDGetTokenT<l1sc::PuppiHostCollection> particles_;
  edm::EDGetTokenT<l1sc::BxLookupHost> bx_;
  edm::EDGetTokenT<l1sc::ClusterObjHostCollection> jets_;
  edm::EDGetTokenT<l1sc::BxLookupHost> jetBx_;
  edm::EDGetTokenT<l1sc::ClustersHostCollection> clusters_;
  edm::EDGetTokenT<l1sc::AssociationMapHost> map_;
  unsigned int minConstituents_;
  bool capped_;
};

DEFINE_FWK_MODULE(Phase2JetTestInput);
DEFINE_FWK_MODULE(Phase2JetTestCheck);
