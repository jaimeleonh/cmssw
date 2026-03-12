// Alpaka EDProducer wrapper around clustering kernels

// includes for types, function/class declarations, utilities/macros
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/Phase2/interface/L1TScPhase2Common.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"

#include <cmath>
#include <string>

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc { // place module in backend-specific Alpaka namespace

  using namespace ::l1sc;

  class L1TScPhase2SCJets : public stream::EDProducer<> { // new class inheriting from stream::EDProducer -> construct once per stream, call produce() each event/orbit
  public:
    // constructor (declare consumed & produced products; read config param at job start & store in member variables)
    L1TScPhase2SCJets(const edm::ParameterSet& params)
          // init parent base class
        : EDProducer<>(params),

          // input tokens: declare what input product module will read (e.g. PuppiDeviceCollection from input tag "src")
          src_candidates_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("src"))},

          // output tokens: declare products written by this module (tokens later used with event.emplace()) 
          clusters_token_{produces()},
          jetBXs_token_{produces()},
          jets_token_{produces()},
          map_token_{produces()},

          // module parameters (square Radii once in constructor)
          R2_{std::pow(params.getParameter<double>("rParam"), 2)},
          nJets_{params.getParameter<unsigned int>("nJets")},
          algo_{params.getParameter<std::string>("algo")},
          RSeed2_{std::pow(params.getParameter<double>("RSeed"), 2)},
          RCen2_{std::pow(params.getParameter<double>("RCen"), 2)},
          RClu2_{std::pow(params.getParameter<double>("RClu"), 2)},
          RLink2_{std::pow(params.getParameter<double>("RLink"), 2)},
          alphaSeed_{params.getParameter<double>("alphaSeed")},
          minSeedPt_{params.getParameter<double>("minSeedPt")},
          nCentroidIters_(params.getParameter<unsigned int>("nCentroidIters")) {}


    void produce(device::Event& event, const device::EventSetup& event_setup) override {
      // fetch inputs using token (corresponding to consumes() declarations)
      const auto& src = event.get(src_candidates_token_);
      const auto& bx_lookup = event.get(bx_lookup_token_);

      // create collection object clusters; allocate storage (on device/queue ssociated with this event -> event.queue) for nsrc entries
      const auto nsrc = src.const_view().metadata().size(); //number of PF cands
      auto clusters = ClustersDeviceCollection(nsrc, event.queue());

      // Decide algorithm ("algo_" string to boolean)
      const bool autoMode = (algo_ == "auto");
      const bool wantWeighted = (algo_ == "seededConeNMSWeighted");
      const bool wantSeeded = (algo_ == "seededCone");
      const bool wantIter = (algo_ == "iterative");
      const bool wantLinkTree = (algo_ == "linkTree");

      // validate algo string
      if (!autoMode && !wantWeighted && !wantSeeded && !wantIter && !wantLinkTree) {
        throw cms::Exception("Configuration")
            << "L1TScPhase2SCJets: unknown algo='" << algo_
            << "'. Allowed: auto | seededCone | iterative | seededConeNMSWeighted | linkTree";
      }

      if (wantLinkTree) {
        auto [jetBXs, jets, map] = kernels_.runLinkTree(event.queue(),
                                                        src,
                                                        bx_lookup,
                                                        float(RLink2_),
                                                        float(minSeedPt_),
                                                        clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));
      } else if (wantWeighted) {
        auto [jetBXs, jets, map] = kernels_.runSeededConeNMSWeighted(event.queue(),
                                                                     src,
                                                                     bx_lookup,
                                                                     float(RSeed2_),
                                                                     float(RCen2_),
                                                                     float(RClu2_),
                                                                     float(alphaSeed_),
                                                                     float(minSeedPt_),
                                                                     nCentroidIters_,
                                                                     clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));
      } else {
        // old behaviour
        const bool doIter = wantIter || (autoMode && nJets_ != 0);

        if (doIter) {
          auto [jetBXs, jets, map] = kernels_.run(event.queue(), src, bx_lookup, float(R2_), nJets_, clusters);
          event.emplace(jetBXs_token_, std::move(jetBXs));
          event.emplace(jets_token_, std::move(jets));
          event.emplace(map_token_, std::move(map));
        } else {
          auto [jetBXs, jets, map] = kernels_.run(event.queue(), src, bx_lookup, float(R2_), clusters);
          event.emplace(jetBXs_token_, std::move(jetBXs));
          event.emplace(jets_token_, std::move(jets));
          event.emplace(map_token_, std::move(map));
        }
      }

      // move clustering results to event storage
      event.emplace(clusters_token_, std::move(clusters));
    }

    // define config schema for module (which par accepted, types, defaults)
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) { // called during config validation before job start
      // create object holding definition of one config set
      edm::ParameterSetDescription desc;

      // define input parameters
      desc.add<edm::InputTag>("src");
      desc.add<double>("rParam", 0.4);
      desc.add<unsigned int>("nJets", 0);
      desc.add<std::string>("algo", "auto");  // auto | seededCone | iterative | seededConeNMSWeighted | linkTree
      desc.add<double>("RSeed", 0.3);
      desc.add<double>("RCen", 0.4);
      desc.add<double>("RClu", 0.4);
      desc.add<double>("RLink", 0.3);
      desc.add<double>("alphaSeed", 2.0);
      desc.add<double>("minSeedPt", 0.0);
      desc.add<unsigned int>("nCentroidIters", 1);

      // register this config schema for the default instance of the module -> python can now use these params
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // get device pf data
    const device::EDGetToken<PuppiDeviceCollection> src_candidates_token_;
    // get association map if runScouting=False
    const device::EDGetToken<BxLookupDeviceCollection> bx_lookup_token_;
    // put device clustering data
    const device::EDPutToken<ClustersDeviceCollection> clusters_token_;
    const device::EDPutToken<BxLookupDeviceCollection> jetBXs_token_;
    const device::EDPutToken<ClusterObjDeviceCollection> jets_token_;
    const device::EDPutToken<AssociationMapDevice> map_token_;

    // kernels
    kernels::L1TScPhase2SCJetsKernels kernels_;

    // params
    double R2_;
    unsigned int nJets_;

    std::string algo_;
    double RSeed2_, RCen2_, RClu2_, RLink2_;
    double alphaSeed_;
    double minSeedPt_;
    unsigned int nCentroidIters_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

// register class as Alpaka EDProducer module
DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2SCJets);