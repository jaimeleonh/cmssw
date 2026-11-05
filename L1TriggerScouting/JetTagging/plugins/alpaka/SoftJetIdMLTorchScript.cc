#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/FixedQueueEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernelTorchScript.h"
#include "PhysicsTools/PyTorchAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/PyTorchAlpaka/interface/alpaka/AlpakaModel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  struct BatchIO {
    cms::torch::alpakatools::TensorCollection<Queue> inputs;
    cms::torch::alpakatools::TensorCollection<Queue> outputs;
  };

  class SoftJetIdMLTorchScript : public stream::FixedQueueEDProducer<> {
  public:
    SoftJetIdMLTorchScript(const edm::ParameterSet& params)
        : FixedQueueEDProducer<>(params),
          pf_token_{consumes(params.getParameter<edm::InputTag>("pf"))},
          jets_token_{consumes(params.getParameter<edm::InputTag>("jets"))},
          association_map_token_{consumes(params.getParameter<edm::InputTag>("clusters"))},
          jet_bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("jetBxLookup"))},
          vertex_token_{consumes(params.getParameter<edm::InputTag>("vertices"))},
          vertex_bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("vertexBxLookup"))},
          soft_jet_token_{produces()},
          model_(params.getParameter<edm::FileInPath>("model").fullPath()),
          max_batch_size_{params.getParameter<uint32_t>("maxBatchSize")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model");
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("jets");
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

        SoftJetNewInputDeviceTensor input_tensor(queue, max_batch_size_);
        SoftJetOutputDeviceTensor output_tensor(queue, max_batch_size_);
        auto input_records = input_tensor.view().records();
        auto output_records = output_tensor.view().records();

        for (auto ii = 0; ii < warmupIterations_; ++ii) {
          cms::torch::alpakatools::TensorCollection<Queue> inputs_d, outputs_d;

          // register inputs
          inputs_d.add<SoftJetNewInputTensorSoA>("pf_points", input_records.points());
          inputs_d.add<SoftJetNewInputTensorSoA>("pf_features", input_records.features());
          inputs_d.add<SoftJetNewInputTensorSoA>("pf_mask", input_records.mask());
          inputs_d.add<SoftJetNewInputTensorSoA>("jet_features", input_records.jet_features());

          // register outputs
          outputs_d.add<SoftJetOutputTensorSoA>("softmax", output_records.output());
      
          // inference on dummy data
          model_.forward(queue, inputs_d, outputs_d);
          // model_.forward(queue, inputs_d, outputs_d);
        }
      }

    void produce(device::Event& event, const device::EventSetup& event_setup) override {
      const auto& pf = event.get(pf_token_);
      const auto& jets = event.get(jets_token_);
      const auto& association_map = event.get(association_map_token_);
      const auto& jetBxLookup = event.get(jet_bx_lookup_token_);
      const auto& vertices = event.get(vertex_token_);
      const auto& vertexBxLookup = event.get(vertex_bx_lookup_token_);
      SoftJetNewInputDeviceTensor input_tensor =
          kernels::torchscript::transform(event.queue(), pf, jets, association_map, jetBxLookup, vertices, vertexBxLookup);

      const auto job_size = association_map.const_view().offset().metadata().size() - 1;
      auto output_tensor = SoftJetOutputDeviceTensor(event.queue(), job_size);
      output_tensor.zeroInitialise(event.queue());

      // The part_10k.pt model is batched and consumes precomputed features:
      //   pf_points [B,48,3] (pt, deta, dphi), pf_features [B,48,14], pf_mask [B,48],
      //   jet_features [B,4] (log(1+pt), eta, log(1+mass), log1p(nconst)) -> [B,1] class-1 prob.
      // The inputs are registered directly from the SoA records (no manual copy), and the
      // model runs batched over the jets via the TensorCollection / AlpakaModel interface.
      const auto batch_size = std::min<uint32_t>(job_size, max_batch_size_);

      // records
      auto input_records = input_tensor.view().records();
      auto output_records = output_tensor.view().records();

      cms::torch::alpakatools::BatchedTensorCollection<Queue> inputs, outputs;

      inputs.addBatched<SoftJetNewInputTensorSoA>("pf_points", batch_size, input_records.points());
      inputs.addBatched<SoftJetNewInputTensorSoA>("pf_features", batch_size, input_records.features());
      inputs.addBatched<SoftJetNewInputTensorSoA>("pf_mask", batch_size, input_records.mask());
      inputs.addBatched<SoftJetNewInputTensorSoA>("jet_features", batch_size, input_records.jet_features());

      outputs.addBatched<SoftJetOutputTensorSoA>("softmax", batch_size, output_records.output());

      model_.forward(event.queue(), inputs, outputs);

      // // put device-side product into event
      event.emplace(soft_jet_token_, std::move(output_tensor));

    }

  private:
    const device::EDGetToken<PFCandidateDeviceCollection> pf_token_;
    const device::EDGetToken<ClusterObjDeviceCollection> jets_token_;
    const device::EDGetToken<AssociationMapDevice> association_map_token_;
    const device::EDGetToken<BxLookupDevice> jet_bx_lookup_token_;
    const device::EDGetToken<VertexDeviceCollection> vertex_token_;
    const device::EDGetToken<BxLookupDevice> vertex_bx_lookup_token_;
    const device::EDPutToken<SoftJetOutputDeviceTensor> soft_jet_token_;
    torch::AlpakaModel model_;
    const uint32_t max_batch_size_;
    // warmup iterations
    const int warmupIterations_ = 10;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftJetIdMLTorchScript);
