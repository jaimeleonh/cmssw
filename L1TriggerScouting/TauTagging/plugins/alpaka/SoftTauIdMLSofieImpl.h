#ifndef L1TriggerScouting_TauTagging_SoftTauIdMLSofieImpl_h
#define L1TriggerScouting_TauTagging_SoftTauIdMLSofieImpl_h

#include <Eigen/Core>
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftTauDeviceTensor.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/FixedQueueEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "L1TriggerScouting/TauTagging/interface/SoftTauIdMLCommon.h"
#include "L1TriggerScouting/TauTagging/plugins/alpaka/TransformKernel.h"
#include "L1TriggerScouting/TauTagging/plugins/alpaka/STTSofieKernels.h"
#include "PerfTools/Perfetto/interface/CMSSWPerfettoTrace.h"
#include "PhysicsTools/SofieAlpaka/interface/TensorCollection.h"
#include "PhysicsTools/SofieAlpaka/interface/alpaka/AlpakaModel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

class SoftTauIdMLSofie : public stream::FixedQueueEDProducer<> {
 public:
  SoftTauIdMLSofie(const edm::ParameterSet& params)
      : FixedQueueEDProducer<>(params),
        pf_candidates_token_(
            consumes(params.getParameter<edm::InputTag>("srcCandidates"))),
        bx_clusters_map_token_(
            consumes(params.getParameter<edm::InputTag>("srcBxClustersMap"))),
        cluster_cands_map_token_(consumes(
            params.getParameter<edm::InputTag>("srcClustersCandsMap"))),
        clusters_token_(
            consumes(params.getParameter<edm::InputTag>("srcClusters"))),
        cluster_cands_map_sorted_token_(produces("clusterCandsMapSorted")),
        model_(params.getParameter<edm::FileInPath>("model").fullPath()),
        substep_(parseSubstep(params.getParameter<std::string>("substep"))),
        batch_size_(params.getParameter<uint32_t>("batchSize")) {
    if (substep_ >= Substep::Reshaping) {
      soft_tau_input_token_ = produces("softTauInputs");
    }
    if (substep_ >= Substep::Tagging) {
      soft_tau_output_token_ = produces("softTauOutputs");
    }
  }

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    edm::ParameterSetDescription desc;
    desc.add<edm::InputTag>("srcCandidates");
    desc.add<edm::InputTag>("srcBxClustersMap");
    desc.add<edm::InputTag>("srcClustersCandsMap");
    desc.add<edm::InputTag>("srcClusters");
    desc.add<edm::FileInPath>("model");
    desc.add<std::string>("substep", "tagging");
    desc.add<uint32_t>("batchSize", 32u);
    descriptions.addWithDefaultLabel(desc);
  }

  void produce(device::Event& event,
               const device::EventSetup& event_setup) override {
    CMS_PERFETTO_SCOPE("SoftTauIdMLSofie produce()");
    const auto& pf = event.get(pf_candidates_token_);  // pf collection
    const auto& bx_clusters_map =
        event.get(bx_clusters_map_token_);  // bx -> clusters map
    const auto& cluster_cands_map =
        event.get(cluster_cands_map_token_);  // clusters->candidates map
    const auto& clusters =
        event.get(clusters_token_);  // cluster ID for each candidate

    // sort the clusters->candidates association map by pt
    auto cluster_cands_map_sorted = kernels::sortClustersCandsMap(
        event.queue(), pf, bx_clusters_map, cluster_cands_map, clusters);

    if (substep_ >= Substep::Reshaping) {
      // get filled input tensor
      SoftTauInputDeviceTensor input_tensor =
          kernels::transform(event.queue(), pf, cluster_cands_map_sorted);

      if (substep_ >= Substep::Tagging) {
        CMS_PERFETTO_SCOPE("Run Batched STT Inference");

        // initialize output tensor
        const auto job_size = cluster_cands_map.const_view().offset().metadata().size() - 1;
        auto output_tensor = SoftTauOutputDeviceTensor(event.queue(), job_size);
        output_tensor.zeroInitialise(event.queue());

        // set batch size
        assert(batch_size_ > 0 &&
               "batch_size_ is expected to be greater than zero, as unbatched "
               "inference will probably make the device go out of memory");

        // records
        auto input_records =
            input_tensor.view()
                .records();  // pay A LOT OF attention here to constness
        auto output_records = output_tensor.view().records();

        // batched inference
        auto num_batches = (job_size + batch_size_ - 1) / batch_size_;
        
        for (unsigned int batch_idx = 0; batch_idx < (num_batches - 1); ++batch_idx) {
            // batched tensors
            cms::sofie::alpakatools::TensorCollection<Queue> inputs(batch_size_, job_size), outputs(batch_size_, job_size);

            // add inputs
            inputs.add<SoftTauInputTensorSoA>("features", batch_idx, input_records.features());
            inputs.add<SoftTauInputTensorSoA>("pad_mask", batch_idx, input_records.pad_mask());
            
            // add outputs
            outputs.add<SoftTauOutputTensorSoA>("cls", batch_idx, output_records.cls());
            outputs.add<SoftTauOutputTensorSoA>("pt", batch_idx, output_records.pt());
            outputs.add<SoftTauOutputTensorSoA>("charge", batch_idx, output_records.charge());

            // run inference
            model_.forward(event.queue(), inputs, outputs);
        }   

        // emplace output
        event.emplace(soft_tau_output_token_, std::move(output_tensor));
      }

      // emplace input
      event.emplace(soft_tau_input_token_, std::move(input_tensor));
    }

    // emplace sorted map
    event.emplace(cluster_cands_map_sorted_token_,
                  std::move(cluster_cands_map_sorted));
  }

 private:
  // event query tokens
  const device::EDGetToken<PFCandidateDeviceCollection> pf_candidates_token_;
  // input bx -> clusters map
  const device::EDGetToken<BxLookupDevice> bx_clusters_map_token_;
  // input clusters -> candidates map that must be sorted by pt
  const device::EDGetToken<AssociationMapDevice> cluster_cands_map_token_;
  // cluster ID for each candidate
  const device::EDGetToken<ClustersDeviceCollection> clusters_token_;
  // sorted clusters -> candidates map that is emplaced in the event
  const device::EDPutToken<AssociationMapDevice>
      cluster_cands_map_sorted_token_;
  // put ml input into event
  device::EDPutToken<SoftTauInputDeviceTensor> soft_tau_input_token_;
  // put ml output into event
  device::EDPutToken<SoftTauOutputDeviceTensor> soft_tau_output_token_;
  // model
  sofie::Model<SOFIE_STT_ort_relu_fp32_b2048_s16::Session<alpaka::AccToTag<Acc1D>>> model_;
  // processing step
  Substep substep_;
  // scouting switch
  const uint32_t batch_size_;
  // warmup iterations
  const int warmupIterations_ = 10;
};

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftTauIdMLSofie);
#endif
