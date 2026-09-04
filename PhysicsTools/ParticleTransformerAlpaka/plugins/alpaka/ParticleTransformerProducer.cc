#include <Eigen/Core>
#include <alpaka/alpaka.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <mutex>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/EDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelDefinition.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelTranspose.h"
#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/ParticleTransformerAlgo.h"
#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/TransformKernel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace {

    using DeviceWeights = decltype(cms::alpakatools::make_device_buffer<std::int8_t[]>(
        std::declval<Queue const&>(), std::declval<std::size_t>()));
    using DeviceFloats = decltype(
        cms::alpakatools::make_device_buffer<float[]>(std::declval<Queue const&>(), std::declval<std::size_t>()));

    // The transposed model is 2.1 MB of INT8 weights plus 237 kB of FP32
    // scales, biases and normalisation parameters, and it is read by every
    // block of every inference launch.  One copy per CMSSW stream would
    // multiply the L2 working set by the number of streams sharing a GPU, so
    // the buffers are owned per *device* and shared by all streams on it.
    struct ResidentModel {
      DeviceWeights weights;
      DeviceFloats params;
    };

    struct ModelRegistryEntry {
      Device device;
      std::string weightsFile;
      std::string paramsFile;
      std::shared_ptr<ResidentModel> model;
    };

    std::mutex& modelRegistryMutex() {
      static std::mutex mutex;
      return mutex;
    }

    std::vector<ModelRegistryEntry>& modelRegistry() {
      static std::vector<ModelRegistryEntry> registry;
      return registry;
    }

    template <typename T>
    void readExact(edm::FileInPath const& file, T* destination, std::size_t count) {
      std::ifstream input(file.fullPath(), std::ios::binary | std::ios::ate);
      auto const expectedBytes = count * sizeof(T);
      if (!input || static_cast<std::size_t>(input.tellg()) != expectedBytes) {
        throw cms::Exception("ParticleTransformerModel")
            << "Model file " << file.fullPath() << " has the wrong size; expected " << expectedBytes << " bytes";
      }
      input.seekg(0);
      input.read(reinterpret_cast<char*>(destination), static_cast<std::streamsize>(expectedBytes));
      if (!input) {
        throw cms::Exception("ParticleTransformerModel") << "Could not read model file " << file.fullPath();
      }
    }

    // Reads the exported checkpoint, converts the dense weights to the
    // [in][out] layout the kernel expects and uploads the result, at most once
    // per (device, file pair).
    std::shared_ptr<ResidentModel> residentModel(Queue& queue,
                                                 edm::FileInPath const& weightsFile,
                                                 edm::FileInPath const& paramsFile) {
      auto const device = alpaka::getDev(queue);
      std::lock_guard<std::mutex> guard(modelRegistryMutex());
      for (auto const& entry : modelRegistry()) {
        if (entry.device == device && entry.weightsFile == weightsFile.fullPath() &&
            entry.paramsFile == paramsFile.fullPath())
          return entry.model;
      }

      auto rowMajor = cms::alpakatools::make_host_buffer<std::int8_t[]>(queue, ::part::generated::weights_size);
      auto columnMajor =
          cms::alpakatools::make_host_buffer<std::int8_t[]>(queue, ::part::generated::weights_size);
      auto parameters = cms::alpakatools::make_host_buffer<float[]>(queue, ::part::generated::params_size);
      readExact(weightsFile, rowMajor.data(), ::part::generated::weights_size);
      readExact(paramsFile, parameters.data(), ::part::generated::params_size);

      // The weight blob holds nothing but dense quantized tensors, but a
      // verbatim copy first keeps the buffer defined even if a future export
      // adds padding between them.
      std::copy_n(rowMajor.data(), ::part::generated::weights_size, columnMajor.data());
      auto const covered =
          ::part::transposeWeights(rowMajor.data(), parameters.data(), columnMajor.data());
      if (covered != ::part::generated::weights_size) {
        throw cms::Exception("ParticleTransformerModel")
            << "The transposed model covers " << covered << " of " << ::part::generated::weights_size
            << " floats; the weight file does not match GeneratedModel.h";
      }

      auto model = std::make_shared<ResidentModel>(ResidentModel{
          cms::alpakatools::make_device_buffer<std::int8_t[]>(queue, ::part::generated::weights_size),
          cms::alpakatools::make_device_buffer<float[]>(queue, ::part::generated::params_size)});
      alpaka::memcpy(queue, model->weights, columnMajor);
      alpaka::memcpy(queue, model->params, parameters);
      // Complete the upload while still holding the lock: another stream that
      // finds this entry must be able to launch against it immediately, and
      // the host staging buffers have to stay alive until the copy is done.
      alpaka::wait(queue);

      modelRegistry().push_back(
          ModelRegistryEntry{device, weightsFile.fullPath(), paramsFile.fullPath(), model});
      return model;
    }

  }  // namespace

  class ParticleTransformerProducer : public stream::EDProducer<> {
  public:
    explicit ParticleTransformerProducer(edm::ParameterSet const& params)
        : EDProducer<>(params),
          pfToken_(consumes(params.getParameter<edm::InputTag>("pf"))),
          associationMapToken_(consumes(params.getParameter<edm::InputTag>("clusters"))),
          jetBxLookupToken_(consumes(params.getParameter<edm::InputTag>("jetBxLookup"))),
          vertexToken_(consumes(params.getParameter<edm::InputTag>("vertices"))),
          vertexBxLookupToken_(consumes(params.getParameter<edm::InputTag>("vertexBxLookup"))),
          weightsFile_(params.getParameter<std::string>("weightsFile")),
          paramsFile_(params.getParameter<std::string>("paramsFile")),
          maxBlocks_(params.getParameter<unsigned int>("maxBlocks")),
          synchronizeForTiming_(params.getParameter<bool>("synchronizeForTiming")),
          outputToken_(produces()) {}

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("clusters");
      desc.add<edm::InputTag>("jetBxLookup");
      desc.add<edm::InputTag>("vertices");
      desc.add<edm::InputTag>("vertexBxLookup");
      desc.add<std::string>("weightsFile", "PhysicsTools/ParticleTransformerAlpaka/data/model_weights.int8.bin");
      desc.add<std::string>("paramsFile", "PhysicsTools/ParticleTransformerAlpaka/data/model_params.fp32.bin");
      // Zero launches one block per jet, which is what saturates the device.
      // A non-zero value caps the grid and makes every block loop over several
      // jets; it exists for occupancy studies, not for memory reasons, since
      // every activation of the graph lives in block shared memory and the
      // kernel needs no per-jet global scratch at all.
      desc.add<unsigned int>("maxBlocks", 0);
      desc.add<bool>("synchronizeForTiming", false);
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, device::EventSetup const&) override {
      auto const& pf = event.get(pfToken_);
      auto const& associationMap = event.get(associationMapToken_);
      auto const& jetBxLookup = event.get(jetBxLookupToken_);
      auto const& vertices = event.get(vertexToken_);
      auto const& vertexBxLookup = event.get(vertexBxLookupToken_);
      auto& queue = event.queue();

      if (!model_)
        model_ = residentModel(queue, weightsFile_, paramsFile_);

      auto inputs = part::makeParticleTransformerInputs(
          queue, pf, associationMap, jetBxLookup, vertices, vertexBxLookup);

      auto const numberOfOffsets = associationMap.const_view().offset().metadata().size();
      auto const numberOfJets = numberOfOffsets == 0 ? 0 : numberOfOffsets - 1;
      l1sc::SoftJetOutputDeviceTensor outputs(queue, numberOfJets);
      outputs.zeroInitialise(queue);
      if (numberOfJets != 0) {
        algo_.run(queue,
                  model_->weights.data(),
                  model_->params.data(),
                  inputs,
                  outputs,
                  static_cast<uint32_t>(numberOfJets),
                  maxBlocks_);
        // This option is for FastTimerService diagnosis only.  Normally the
        // queue remains asynchronous; without this wait the framework can
        // account the outstanding inference time to its cleanup transition.
        if (synchronizeForTiming_)
          alpaka::wait(queue);
      }
      event.emplace(outputToken_, std::move(outputs));
    }

  private:
    device::EDGetToken<l1sc::PFCandidateDeviceCollection> const pfToken_;
    device::EDGetToken<l1sc::AssociationMapDevice> const associationMapToken_;
    device::EDGetToken<l1sc::BxLookupDevice> const jetBxLookupToken_;
    device::EDGetToken<l1sc::VertexDeviceCollection> const vertexToken_;
    device::EDGetToken<l1sc::BxLookupDevice> const vertexBxLookupToken_;
    edm::FileInPath const weightsFile_;
    edm::FileInPath const paramsFile_;
    uint32_t const maxBlocks_;
    bool const synchronizeForTiming_;
    device::EDPutToken<l1sc::SoftJetOutputDeviceTensor> const outputToken_;

    std::shared_ptr<ResidentModel> model_;
    ParticleTransformerAlgo algo_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(ParticleTransformerProducer);
