#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/FixedQueueEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernel.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"


namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  struct BatchIO {
    cms::torch::alpakatools::TensorCollection<Queue> inputs;
    cms::torch::alpakatools::TensorCollection<Queue> outputs;
  };

  class SoftJetIdML : public stream::FixedQueueEDProducer<> {
  public:
    SoftJetIdML(const edm::ParameterSet &params)
        : FixedQueueEDProducer<>(params),
          pf_token_{consumes(params.getParameter<edm::InputTag>("pf"))},
          association_map_token_{consumes(params.getParameter<edm::InputTag>("clusters"))},
          jet_bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("jetBxLookup"))},
          vertex_token_{consumes(params.getParameter<edm::InputTag>("vertices"))},
          vertex_bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("vertexBxLookup"))},
          soft_jet_token_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          max_batch_size_{params.getParameter<uint32_t>("maxBatchSize")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model");
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("clusters");
      desc.add<edm::InputTag>("jetBxLookup");
      desc.add<edm::InputTag>("vertices");
      desc.add<edm::InputTag>("vertexBxLookup");
      desc.add<uint32_t>("maxBatchSize", std::numeric_limits<uint32_t>::max());
      descriptions.addWithDefaultLabel(desc);
    }

    void beginStream(edm::StreamID, Queue queue) {
        // Warming up the model with dummy data. This allows for fast inference since the first iteration
        // and should ensure a more homogeneous memory balance

        SoftJetInputDeviceTensor input_tensor(queue, max_batch_size_);
        SoftJetOutputDeviceTensor output_tensor(queue, max_batch_size_);
        auto input_records = input_tensor.view().records();
        auto output_records = output_tensor.view().records();

        for (auto ii = 0; ii < warmupIterations_; ++ii) {
          cms::torch::alpakatools::TensorCollection<Queue> inputs_d, outputs_d;

          // register inputs
          inputs_d.add<SoftJetInputTensorSoA>("pf_points", input_records.points());
          inputs_d.add<SoftJetInputTensorSoA>("pf_features", input_records.features());
          inputs_d.add<SoftJetInputTensorSoA>("pf_vectors", input_records.vectors());
          inputs_d.add<SoftJetInputTensorSoA>("pf_mask", input_records.mask());

          // register outputs
          outputs_d.add<SoftJetOutputTensorSoA>("softmax", output_records.output());
      
          // inference on dummy data
          // model_.forward(queue, inputs_d, outputs_d, ::torch::kHalf);
          model_.forward(queue, inputs_d, outputs_d);
        }
      }

    void produce(device::Event &event, const device::EventSetup &event_setup) override {
      // in/out collections
      const auto &pf = event.get(pf_token_);
      const auto &association_map = event.get(association_map_token_);
      const auto &jetBxLookup = event.get(jet_bx_lookup_token_);
      const auto &vertices = event.get(vertex_token_);
      const auto &vertexBxLookup = event.get(vertex_bx_lookup_token_);
      SoftJetInputDeviceTensor input_tensor = kernels::transform(event.queue(), pf, association_map, jetBxLookup, vertices, vertexBxLookup);

      const auto job_size = association_map.const_view().offset().metadata().size() - 1;
      auto output_tensor = SoftJetOutputDeviceTensor(event.queue(), job_size);
      output_tensor.zeroInitialise(event.queue());

      const auto batch_size = std::min<uint32_t>(job_size, max_batch_size_);
      // printf("Batch ratios %i %i %i\n", job_size, num_batches, num_batches * batch_size);

      // records
      auto input_records = input_tensor.view().records();
      auto output_records = output_tensor.view().records();

      // batched tensors
      cms::torch::alpakatools::BatchedTensorCollection<Queue> inputs, outputs;

      // register inputs
      inputs.addBatched<SoftJetInputTensorSoA>("pf_points", batch_size, input_records.points());
      inputs.addBatched<SoftJetInputTensorSoA>("pf_features", batch_size, input_records.features());
      inputs.addBatched<SoftJetInputTensorSoA>("pf_vectors", batch_size, input_records.vectors());
      inputs.addBatched<SoftJetInputTensorSoA>("pf_mask", batch_size, input_records.mask());

      // register outputs
      outputs.addBatched<SoftJetOutputTensorSoA>("softmax", batch_size, output_records.output());

      // do inference
      // model_.forward(event.queue(), inputs, outputs, ::torch::kHalf);
      model_.forward(event.queue(), inputs, outputs);

      // // put device-side product into event
      event.emplace(soft_jet_token_, std::move(output_tensor));
    }

  private:
    // event query tokens
    const device::EDGetToken<PFCandidateDeviceCollection> pf_token_;
    // clustering output
    const device::EDGetToken<ClustersDeviceCollection> clusters_token_;
    const device::EDGetToken<AssociationMapDevice> association_map_token_;
    const device::EDGetToken<BxLookupDevice> jet_bx_lookup_token_;
    const device::EDGetToken<VertexDeviceCollection> vertex_token_;
    const device::EDGetToken<BxLookupDevice> vertex_bx_lookup_token_;
    // put ml output into event
    const device::EDPutToken<SoftJetOutputDeviceTensor> soft_jet_token_;
    // model
    torch::AlpakaModel model_;
    // scouting switch
    const uint32_t max_batch_size_;
    // warmup iterations
    const int warmupIterations_ = 10;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::torchtest

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftJetIdML);
