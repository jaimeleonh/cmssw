#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernel.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  struct BatchIO {
    cms::torch::alpakatools::TensorCollection<Queue> inputs;
    cms::torch::alpakatools::TensorCollection<Queue> outputs;
  };

  class SoftJetIdML : public stream::EDProducer<> {
  public:
    SoftJetIdML(const edm::ParameterSet &params)
        : EDProducer<>(params),
          pf_token_(consumes(params.getParameter<edm::InputTag>("pf"))),
          association_map_token_{consumes(params.getParameter<edm::InputTag>("clusters"))},
          soft_jet_token_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          max_batch_size_{params.getParameter<uint32_t>("maxBatchSize")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model");
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("clusters");
      desc.add<uint32_t>("maxBatchSize", std::numeric_limits<uint32_t>::max());
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // in/out collections
      const auto &pf = event.get(pf_token_);
      const auto &association_map = event.get(association_map_token_);
      SoftJetInputDeviceTensor input_tensor = kernels::transform(event.queue(), pf, association_map);

      const auto job_size = association_map.const_view().offset().metadata().size() - 1;
      auto output_tensor = SoftJetOutputDeviceTensor(event.queue(), job_size);
      output_tensor.zeroInitialise(event.queue());

      const auto batch_size = std::min<uint32_t>(job_size, max_batch_size_);
      auto num_batches = (job_size + batch_size - 1) / batch_size;
      // printf("Batch ratios %i %i %i\n", job_size, num_batches, num_batches * batch_size);

      // records
      auto input_records = input_tensor.view().records();
      auto output_records = output_tensor.view().records();

      std::deque<BatchIO> batches;
      for (auto batch_idx = 0; batch_idx < num_batches; ++batch_idx) {
        // std::cout << "Batch " << batch_idx << std::endl;
        BatchIO batch{cms::torch::alpakatools::TensorCollection<Queue>(batch_size, job_size),
                      cms::torch::alpakatools::TensorCollection<Queue>(batch_size, job_size)};

        batch.inputs.add<SoftJetInputTensorSoA>("pf_points", batch_idx, input_records.points());
        batch.inputs.add<SoftJetInputTensorSoA>("pf_features", batch_idx, input_records.features());
        batch.inputs.add<SoftJetInputTensorSoA>("pf_vectors", batch_idx, input_records.vectors());
        batch.inputs.add<SoftJetInputTensorSoA>("pf_mask", batch_idx, input_records.mask());

        batch.outputs.add<SoftJetOutputTensorSoA>("softmax", batch_idx, output_records.output());

        batches.push_back(std::move(batch));
      }

      for (auto &batch : batches) {
        c10::InferenceMode guard(true);
        model_.to(event.queue());
        model_.forward(event.queue(), batch.inputs, batch.outputs);
      }

      // auto out = output_tensor.view();
      // auto output = out.output();

      // for (uint32_t i = 0; i <  2 * batch_size; ++i) {
      //   std::cout << "jet " << i
      //             << " score=" << output[i]
      //             << std::endl;
      // }

      // // put device-side product into event
      event.emplace(soft_jet_token_, std::move(output_tensor));
    }

  private:
    // event query tokens
    const device::EDGetToken<PFCandidateDeviceCollection> pf_token_;
    // clustering output
    const device::EDGetToken<ClustersDeviceCollection> clusters_token_;
    const device::EDGetToken<AssociationMapDevice> association_map_token_;
    // put ml output into event
    const device::EDPutToken<SoftJetOutputDeviceTensor> soft_jet_token_;
    // model
    torch::AlpakaModel model_;
    // scouting switch
    const uint32_t max_batch_size_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftJetIdML);
