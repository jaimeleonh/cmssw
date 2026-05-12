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

#include <algorithm>
#include <cmath>
#include <string>

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc { // place module in backend-specific Alpaka namespace

  using namespace ::l1sc;

  class L1TScPhase2SCJets : public stream::EDProducer<> { // new class inheriting from stream::EDProducer -> construct once per stream, call produce() each event/orbit
  public:
    // constructor (declare consumed & produced products; read config param at job start & store in member variables)
   L1TScPhase2SCJets(const edm::ParameterSet& params)
    : EDProducer<>(params),

      src_candidates_token_{consumes(params.getParameter<edm::InputTag>("src"))},
      bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("src"))},

      clusters_token_{produces()},
      jetBXs_token_{produces()},
      jets_token_{produces()},
      map_token_{produces()},

      algo_{params.getParameter<std::string>("algo")},

      // safe defaults; real values are filled below depending on algo_
      R2_{0.0},
      nJets_{0},
      RSeed2_{0.0},
      RCen2_{0.0},
      RClu2_{0.0},
      RLink2_{0.0},
      alphaSeed_{0.0},
      minSeedPt_{0.0},
      nCentroidIters_{0} {
      if (algo_ == "SCGreedy") {
        R2_ = std::pow(params.getParameter<double>("rParam"), 2);
        nJets_ = params.getParameter<unsigned int>("nJets");

      } else if (algo_ == "SCNMS" || algo_ == "SCNMSWeighted") {
        RSeed2_ = std::pow(params.getParameter<double>("RSeed"), 2);
        RClu2_ = std::pow(params.getParameter<double>("RClu"), 2);

      } else if (algo_ == "SCNMSWeightedMultiIter") {
        RSeed2_ = std::pow(params.getParameter<double>("RSeed"), 2);
        RCen2_ = std::pow(params.getParameter<double>("RCen"), 2);
        RClu2_ = std::pow(params.getParameter<double>("RClu"), 2);
        alphaSeed_ = params.getParameter<double>("alphaSeed");
        minSeedPt_ = params.getParameter<double>("minSeedPt");
        nCentroidIters_ = params.getParameter<unsigned int>("nCentroidIters");

      } else if (algo_ == "LinkTree") {
        RLink2_ = std::pow(params.getParameter<double>("RLink"), 2);
        minSeedPt_ = params.getParameter<double>("minSeedPt");

      } else {
        throw cms::Exception("Configuration")
            << "Unsupported algo '" << algo_ << "'";
      }
    }

    void produce(device::Event& event, const device::EventSetup& event_setup) override {
      // fetch inputs using token (corresponding to consumes() declarations)
      const auto& src = event.get(src_candidates_token_);
      const auto& bx_lookup = event.get(bx_lookup_token_);

      // create collection object clusters; allocate storage (on device/queue associated with this event -> event.queue) for nsrc entries
      const auto nsrc = src.const_view().metadata().size(); // number of PF cands
      auto clusters = ClustersDeviceCollection(nsrc, event.queue());

      // canonical names
      const bool wantSCGreedy               = (algo_ == "SCGreedy");
      const bool wantSCNMS                  = (algo_ == "SCNMS");
      const bool wantSCNMSWeighted          = (algo_ == "SCNMSWeighted");
      const bool wantSCNMSWeightedMultiIter = (algo_ == "SCNMSWeightedMultiIter");
      const bool wantLinkTree               = (algo_ == "LinkTree");

      // validate algo string
      if (!wantSCGreedy &&
          !wantSCNMS &&
          !wantSCNMSWeighted &&
          !wantSCNMSWeightedMultiIter &&
          !wantLinkTree) {
        throw cms::Exception("Configuration")
            << "L1TScPhase2SCJets: unknown algo='" << algo_ << "'. "
            << "Allowed: SCGreedy | SCNMS | SCNMSWeighted | SCNMSWeightedMultiIter | LinkTree ";
      }

      // ----------------------------------------------------------------
      // dispatch
      //
      // Algorithm mapping:
      //
      //   SCNMS:
      //     old kernel logic, but with split radii
      //       RSeed -> seed finding + old centroid accumulation
      //       RClu  -> final nearest-axis assignment
      //
      //   SCNMSWeighted:
      //     same as SCNMS, but final assignment uses the old weighted metric
      //
      //   SCNMSWeightedMultiIter:
      //     newer generic NMS-style implementation with explicit
      //     seed/centroid/assignment radii and forced >= 2 centroid iters
      //
      // So for SCNMS / SCNMSWeighted, RCen is intentionally ignored.
      // ----------------------------------------------------------------
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

      } else if (wantSCNMSWeightedMultiIter) {
        // force multiple centroid iterations here even if config forgot to set it
        unsigned int nIters = std::max(2u, nCentroidIters_);

        auto [jetBXs, jets, map] = kernels_.runSCNMSWeightedMultiIter(event.queue(),
                                                                      src,
                                                                      bx_lookup,
                                                                      float(RSeed2_),
                                                                      float(RCen2_),
                                                                      float(RClu2_),
                                                                      float(alphaSeed_),
                                                                      float(minSeedPt_),
                                                                      nIters,
                                                                      clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));

      } else if (wantSCNMSWeighted) {
        auto [jetBXs, jets, map] = kernels_.runSCNMSWeighted(event.queue(),
                                                             src,
                                                             bx_lookup,
                                                             float(RSeed2_),
                                                             float(RClu2_),
                                                             clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));

      } else if (wantSCNMS) {
        auto [jetBXs, jets, map] = kernels_.runSCNMS(event.queue(),
                                                     src,
                                                     bx_lookup,
                                                     float(RSeed2_),
                                                     float(RClu2_),
                                                     clusters);
        event.emplace(jetBXs_token_, std::move(jetBXs));
        event.emplace(jets_token_, std::move(jets));
        event.emplace(map_token_, std::move(map));

      } else {
        // old behaviour for SCGreedy:
        // if nJets != 0 -> iterative SCGreedy
        // if nJets == 0 -> legacy single-radius non-iterative seeded cone
        const bool doIter = (nJets_ != 0);

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

    // define config schema for module (which params accepted, types, defaults)
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) { // called during config validation before job start
      // create object holding definition of one config set
      edm::ParameterSetDescription desc;

      // define input parameters
      desc.add<edm::InputTag>("src");

      // canonical names:
      // SCGreedy | SCNMS | SCNMSWeighted | SCNMSWeightedMultiIter | LinkTree
      desc.ifValue(edm::ParameterDescription<std::string>("algo", "None", true,
          edm::Comment("Allowed values: SCGreedy, SCNMS, SCNMSWeighted, SCNMSWeightedMultiIter, LinkTree")),
        "SCGreedy" >> (edm::ParameterDescription<double>("rParam", 0.3, true) and
          edm::ParameterDescription<unsigned int>("nJets", 0, true)
        ) or
        
        "SCNMS" >> (edm::ParameterDescription<double>("RSeed", 0.3, true) and
          edm::ParameterDescription<double>("RClu", 0.4, true)
        ) or
        
        "SCNMSWeighted" >> (edm::ParameterDescription<double>("RSeed", 0.3, true) and
          edm::ParameterDescription<double>("RClu", 0.4, true)
        ) or

        "SCNMSWeightedMultiIter" >> (edm::ParameterDescription<double>("RSeed", 0.3, true) and
          edm::ParameterDescription<double>("RCen", 0.4, true) and
          edm::ParameterDescription<double>("RClu", 0.4, true) and
          edm::ParameterDescription<double>("alphaSeed", 2.0, true) and
          edm::ParameterDescription<double>("minSeedPt", 0.0, true) and
          edm::ParameterDescription<unsigned int>("nCentroidIters", 1, true)
        ) or
        
        "LinkTree" >> (edm::ParameterDescription<double>("RLink", 0.3, true) and
          edm::ParameterDescription<double>("minSeedPt", 0.0, true)
        ) or 

        "None" >> edm::EmptyGroupDescription()
      );

      // register this config schema for the default instance of the module -> python can now use these params
      descriptions.addWithDefaultLabel(desc);
    }

  private:
    // get device pf data
    const device::EDGetToken<PuppiDeviceCollection> src_candidates_token_;
    // get BX lookup from same input tag
    const device::EDGetToken<BxLookupDeviceCollection> bx_lookup_token_;

    // put device clustering data
    const device::EDPutToken<ClustersDeviceCollection> clusters_token_;
    const device::EDPutToken<BxLookupDeviceCollection> jetBXs_token_;
    const device::EDPutToken<ClusterObjDeviceCollection> jets_token_;
    const device::EDPutToken<AssociationMapDevice> map_token_;

    // kernels
    kernels::L1TScPhase2SCJetsKernels kernels_;

    // params
    std::string algo_;
    double R2_;
    unsigned int nJets_;

    double RSeed2_, RCen2_, RClu2_, RLink2_;
    double alphaSeed_;
    double minSeedPt_;
    unsigned int nCentroidIters_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

// register class as Alpaka EDProducer module
DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2SCJets);
