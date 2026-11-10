#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
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
#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"
#include <cmath>
#include <string>
namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {
  using namespace ::l1sc;
  class L1TScPhase2SCJets : public stream::EDProducer<> {
  public:
    explicit L1TScPhase2SCJets(edm::ParameterSet const& p)
        : EDProducer<>(p),
          src_{consumes(p.getParameter<edm::InputTag>("src"))},
          bx_{consumes(p.getParameter<edm::InputTag>("src"))},
          clusters_{produces()},
          jetBx_{produces()},
          jets_{produces()},
          map_{produces()},
          kernels_{p.getParameter<unsigned int>("minConstituents")},
          algo_{p.getParameter<std::string>("algo")},
          nJets_{p.getParameter<unsigned int>("nJets")},
          r2_{std::pow(p.getParameter<double>("rParam"), 2)},
          rSeed2_{std::pow(p.getParameter<double>("RSeed"), 2)},
          rCen2_{std::pow(p.getParameter<double>("RCen"), 2)},
          rClu2_{std::pow(p.getParameter<double>("RClu"), 2)},
          rLink2_{std::pow(p.getParameter<double>("RLink"), 2)},
          alpha_{p.getParameter<double>("alphaSeed")},
          minPt_{p.getParameter<double>("minSeedPt")},
          maxRounds_{p.getParameter<int>("maxRounds")} {
      if (algo_ != "SeededCone" && algo_ != "SCNMSW" && algo_ != "LinkTree")
        throw cms::Exception("Configuration") << "algo must be SeededCone, SCNMSW or LinkTree";
      if (algo_ == "SeededCone" && !nJets_)
        throw cms::Exception("Configuration") << "SeededCone needs a positive nJets iteration limit";
      if (maxRounds_ == 0 || maxRounds_ < -1)
        throw cms::Exception("Configuration") << "maxRounds must be -1 or a positive integer";
      for (auto radius : {"rParam", "RSeed", "RCen", "RClu", "RLink"}) {
        double value = p.getParameter<double>(radius);
        if (!(value > 0.) || !std::isfinite(value) || !std::isfinite(float(value * value)))
          throw cms::Exception("Configuration") << radius << " must be finite and positive";
      }
      if (!std::isfinite(float(alpha_)) || alpha_ < 0. || !std::isfinite(float(minPt_)) || minPt_ < 0.)
        throw cms::Exception("Configuration") << "alphaSeed and minSeedPt must be finite and nonnegative";
    }
    void produce(device::Event& event, device::EventSetup const&) override {
      auto const& src = event.get(src_);
      auto const& bx = event.get(bx_);
      auto clusters = ClustersDeviceCollection(event.queue(), src.const_view().metadata().size());
      auto result = [&]() {
        if (algo_ == "SeededCone")
          return kernels_.runSeededCone(event.queue(), src, bx, float(r2_), nJets_, clusters);
        if (algo_ == "LinkTree")
          return kernels_.runLinkTree(event.queue(), src, bx, float(rLink2_), float(minPt_), clusters);
        return kernels_.runSCNMSW(
            event.queue(), src, bx, float(rSeed2_), float(rCen2_), float(rClu2_), float(alpha_), maxRounds_, clusters);
      }();
      auto [jetBx, jets, map] = std::move(result);
      event.emplace(jetBx_, std::move(jetBx));
      event.emplace(jets_, std::move(jets));
      event.emplace(map_, std::move(map));
      event.emplace(clusters_, std::move(clusters));
    }
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription d;
      d.add<edm::InputTag>("src");
      d.add<std::string>("algo", "SeededCone");
      d.add<double>("rParam", .4);
      d.add<unsigned int>("nJets", 16);
      d.add<double>("RSeed", .4);
      d.add<double>("RCen", .4);
      d.add<double>("RClu", .4);
      d.add<double>("RLink", .36);
      d.add<double>("alphaSeed", 2.);
      d.add<double>("minSeedPt", 0.);
      d.add<int>("maxRounds", 1);
      d.add<unsigned int>("minConstituents", 0);
      descriptions.addWithDefaultLabel(d);
    }

  private:
    device::EDGetToken<PuppiDeviceCollection> src_;
    device::EDGetToken<BxLookupDevice> bx_;
    device::EDPutToken<ClustersDeviceCollection> clusters_;
    device::EDPutToken<BxLookupDevice> jetBx_;
    device::EDPutToken<ClusterObjDeviceCollection> jets_;
    device::EDPutToken<AssociationMapDevice> map_;
    kernels::L1TScPhase2SCJetsKernels kernels_;
    std::string algo_;
    unsigned int nJets_;
    double r2_, rSeed2_, rCen2_, rClu2_, rLink2_, alpha_, minPt_;
    int maxRounds_;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc
DEFINE_FWK_ALPAKA_MODULE(l1sc::L1TScPhase2SCJets);
