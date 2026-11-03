#include <deque>
#include <limits>

#include <Eigen/Core>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetONNXDeviceTensor.h"
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
#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernelONNX.h"
#include "PhysicsTools/ONNXRuntimeAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/ONNXRuntimeAlpaka/interface/alpaka/AlpakaSession.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  struct BatchIO {
    cms::Ort::alpakatools::TensorCollection<Queue> inputs;
    cms::Ort::alpakatools::TensorCollection<Queue> outputs;
  };

  // ONNX Runtime counterpart of SoftJetIdML: jets' constituents are transformed
  // on device into the ParticleTransformer input tensors, scored with ONNX
  // Runtime through PhysicsTools/ONNXRuntimeAlpaka, and the b-jet score is
  // written to the same SoftJetOutputDeviceTensor produced by the PyTorch module.
  //
  // Legacy ParticleTransformer graphs exported with advanced indexed assignment
  // require maxBatchSize=1. Re-export with batch-preserving pair scatters to
  // support larger batches and a smaller final batch.
  class SoftJetIdONNX : public stream::FixedQueueEDProducer<> {
  public:
    SoftJetIdONNX(const edm::ParameterSet& params)
        : FixedQueueEDProducer<>(params),
          pf_token_{consumes(params.getParameter<edm::InputTag>("pf"))},
          association_map_token_{consumes(params.getParameter<edm::InputTag>("clusters"))},
          jet_bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("jetBxLookup"))},
          vertex_token_{consumes(params.getParameter<edm::InputTag>("vertices"))},
          vertex_bx_lookup_token_{consumes(params.getParameter<edm::InputTag>("vertexBxLookup"))},
          soft_jet_token_{produces()},
          session_(params.getParameter<edm::FileInPath>("model").fullPath()),
          max_batch_size_{params.getParameter<uint32_t>("maxBatchSize")} {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::FileInPath>("model");
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("clusters");
      desc.add<edm::InputTag>("jetBxLookup");
      desc.add<edm::InputTag>("vertices");
      desc.add<edm::InputTag>("vertexBxLookup");
      desc.add<uint32_t>("maxBatchSize", 1);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, const device::EventSetup& event_setup) override {
      // input collections
      const auto& pf = event.get(pf_token_);
      const auto& association_map = event.get(association_map_token_);
      const auto& jetBxLookup = event.get(jet_bx_lookup_token_);
      const auto& vertices = event.get(vertex_token_);
      const auto& vertexBxLookup = event.get(vertex_bx_lookup_token_);
      SoftJetONNXInputDeviceTensor input_tensor =
          kernels::transformONNX(event.queue(), pf, association_map, jetBxLookup, vertices, vertexBxLookup);

      const auto job_size = association_map.const_view().offset().metadata().size() - 1;
      auto output_tensor = SoftJetOutputDeviceTensor(event.queue(), job_size);
      output_tensor.zeroInitialise(event.queue());

      if (job_size > 0) {
        const auto batch_size = std::min<uint32_t>(job_size, max_batch_size_);
        const auto num_batches = (job_size + batch_size - 1) / batch_size;

        // raw softmax scores, one [batch, 2] tensor per batch
        auto scores = SoftJetONNXScoresDevice(event.queue(), job_size);

        auto input_records = input_tensor.view().records();
        auto output_records = scores.view().records();

        std::deque<BatchIO> batches;
        for (uint32_t batch_idx = 0; batch_idx < num_batches; ++batch_idx) {
          BatchIO batch{cms::Ort::alpakatools::TensorCollection<Queue>(batch_size, job_size),
                        cms::Ort::alpakatools::TensorCollection<Queue>(batch_size, job_size)};

          batch.inputs.add<SoftJetONNXInputTensorSoA>("pf_features", batch_idx, input_records.features());
          batch.inputs.add<SoftJetONNXInputTensorSoA>("pf_vectors", batch_idx, input_records.vectors());
          batch.inputs.add<SoftJetONNXInputTensorSoA>("pf_mask", batch_idx, input_records.mask());

          batch.outputs.add<SoftJetONNXScoresSoA>("softmax", batch_idx, output_records.b(), output_records.notB());

          batches.push_back(std::move(batch));
        }

        for (auto& batch : batches) {
          // ONNX Runtime writes the softmax scores to the staged [batch, 2]
          // buffers, and unpacks them into the scores SoA columns
          session_.forward(event.queue(), batch.inputs, batch.outputs);
        }

        // extract the b-jet score into the output tensor
        kernels::extractScores(event.queue(), scores, output_tensor);
      }

      // put device-side product into event
      event.emplace(soft_jet_token_, std::move(output_tensor));
    }

  private:
    // event query tokens
    const device::EDGetToken<PFCandidateDeviceCollection> pf_token_;
    // clustering output
    const device::EDGetToken<AssociationMapDevice> association_map_token_;
    const device::EDGetToken<BxLookupDevice> jet_bx_lookup_token_;
    const device::EDGetToken<VertexDeviceCollection> vertex_token_;
    const device::EDGetToken<BxLookupDevice> vertex_bx_lookup_token_;
    // put ml output into event
    const device::EDPutToken<SoftJetOutputDeviceTensor> soft_jet_token_;
    // ONNX Runtime session bound to the module queue
    ort::AlpakaSession session_;
    // batch size for the inference
    const uint32_t max_batch_size_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftJetIdONNX);
