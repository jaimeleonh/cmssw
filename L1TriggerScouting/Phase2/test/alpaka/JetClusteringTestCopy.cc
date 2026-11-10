#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {
  class Phase2JetTestCopy : public stream::EDProducer<> {
  public:
    explicit Phase2JetTestCopy(edm::ParameterSet const& p)
        : EDProducer<>(p), hostParticles_{consumes(p.getParameter<edm::InputTag>("src"))},
          hostBx_{consumes(p.getParameter<edm::InputTag>("src"))}, particles_{produces()}, bx_{produces()} {}
    void produce(device::Event& event, device::EventSetup const&) override {
      auto const& hostParticles = event.get(hostParticles_);
      auto const& hostBx = event.get(hostBx_);
      auto particles = PuppiDeviceCollection(event.queue(), hostParticles.const_view().metadata().size());
      auto bx = BxLookupDevice(event.queue(), hostBx.const_view().bx().metadata().size(),
                               hostBx.const_view().offset().metadata().size());
      alpaka::memcpy(event.queue(), particles.buffer(), hostParticles.buffer());
      alpaka::memcpy(event.queue(), bx.buffer(), hostBx.buffer());
      alpaka::wait(event.queue());
      event.emplace(particles_, std::move(particles));
      event.emplace(bx_, std::move(bx));
    }
    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription d;
      d.add<edm::InputTag>("src");
      descriptions.addWithDefaultLabel(d);
    }
  private:
    edm::EDGetTokenT<::l1sc::PuppiHostCollection> hostParticles_;
    edm::EDGetTokenT<::l1sc::BxLookupHost> hostBx_;
    device::EDPutToken<PuppiDeviceCollection> particles_;
    device::EDPutToken<BxLookupDevice> bx_;
  };
}
DEFINE_FWK_ALPAKA_MODULE(l1sc::Phase2JetTestCopy);
