#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;

  class L1TScPhase2SCJets : public stream::EDProducer<> {
  public:
    L1TScPhase2SCJets(const edm::ParameterSet &params)
        : EDProducer<>(params),
          src_candidates_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("bxLookup"))},
          clusters_token_{produces()},
          jetBXs_token_{produces()},
          jets_token_{produces()},
          map_token_{produces()},
          R2_{std::pow(params.getParameter<double>("rParam"), 2)},
          nJets_{params.getParameter<unsigned int>("nJets")} {}

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      const auto &src = event.get(src_candidates_token_);
      const auto &bx_lookup = event.get(bx_lookup_token_);

      const auto nsrc = src.const_view().metadata().size();
      auto clusters = ClustersDeviceCollection(event.queue(), int(nsrc));

      if (nJets_ == 0) {
        auto [jetBXs, jets, map] = kernels_.run(event.queue(), src, bx_lookup, R2_, clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));
      } else {
        auto [jetBXs, jets, map] = kernels_.run(event.queue(), src, bx_lookup, R2_, nJets_, clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));
      }

      event.emplace(clusters_token_, std::move(clusters));
    };

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src", edm::InputTag("l1tScPhase2PuppiRawToDigi", "candidates"));
      desc.add<edm::InputTag>("bxLookup", edm::InputTag("l1tScPhase2PuppiRawToDigi", "bxLookup"));
      desc.add<double>("rParam", 0.4);
      desc.add<unsigned int>("nJets", 0);
      descriptions.addWithDefaultLabel(desc);
    };

  private:
    const device::EDGetToken<PuppiDeviceCollection> src_candidates_token_;
    const device::EDGetToken<BxLookupDevice> bx_lookup_token_;

    const device::EDPutToken<ClustersDeviceCollection> clusters_token_;
    const device::EDPutToken<BxLookupDevice> jetBXs_token_;
    const device::EDPutToken<ClusterObjDeviceCollection> jets_token_;
    const device::EDPutToken<AssociationMapDevice> map_token_;

    kernels::L1TScPhase2SCJetsKernels kernels_;

    double R2_;
    unsigned int nJets_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2SCJets);
