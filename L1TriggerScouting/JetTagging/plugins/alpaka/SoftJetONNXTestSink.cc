#include <cmath>
#include <cstdlib>
#include <sstream>
#include <vector>

#include <Eigen/Core>

#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/FixedQueueEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  // Copy the SoftJetOutputDeviceTensor back to the host and compare the jet
  // scores with the expected values (see runSoftJetTaggingONNXTest.py).
  // Any mismatch above the tolerance fails the job.
  class SoftJetONNXTestSink : public stream::FixedQueueEDProducer<> {
  public:
    SoftJetONNXTestSink(const edm::ParameterSet& params)
        : FixedQueueEDProducer<>(params),
          src_token_{consumes(params.getParameter<edm::InputTag>("src"))},
          expected_{params.getParameter<std::vector<double>>("expectedScores")},
          tolerance_{params.getParameter<double>("tolerance")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("src");
      desc.add<std::vector<double>>("expectedScores");
      desc.add<double>("tolerance", 1e-4);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, const device::EventSetup& event_setup) override {
      const auto& device_scores = event.get(src_token_);
      const auto size = device_scores.const_view().metadata().size();

      // copy the scores to the host and wait for the copy to complete
      auto host_scores = SoftJetOutputHostTensor(event.queue(), size);
      alpaka::memcpy(event.queue(), host_scores.buffer(), device_scores.buffer());
      alpaka::wait(event.queue());
      auto view = host_scores.const_view();

      if (static_cast<size_t>(size) != expected_.size()) {
        throw cms::Exception("TestFailure")
            << "SoftJetONNXTestSink: got " << size << " jet scores, expected " << expected_.size() << ".";
      }

      bool failed = false;
      std::ostringstream report;
      for (uint32_t i = 0; i < static_cast<uint32_t>(size); ++i) {
        const float score = view.output()[i];
        const double expected = expected_[i];
        const bool ok = std::abs(score - expected) <= tolerance_ + tolerance_ * std::abs(expected);
        failed |= not ok;
        report << "\n  jet " << i << ": score " << score << " expected " << expected << (ok ? "" : "  <-- MISMATCH");
      }
      if (failed) {
        throw cms::Exception("TestFailure") << "SoftJetONNXTestSink: score mismatch:" << report.str();
      }
      edm::LogPrint("SoftJetONNXTestSink") << "All " << size << " jet scores within " << tolerance_ << " of the"
                                           << " expected values:" << report.str();
    }

  private:
    const device::EDGetToken<SoftJetOutputDeviceTensor> src_token_;
    const std::vector<double> expected_;
    const double tolerance_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftJetONNXTestSink);
