// A producer that runs whatever model it is pointed at.
//
// It differs from ParticleTransformerProducer in one way that matters: the
// network is a data file, not a compiled-in graph.  Swapping models -- a
// retrained checkpoint, a different width, a different number of blocks -- is
// a configuration change, and nothing in this plugin or the kernel it launches
// has to be rebuilt.
//
// The upload is per *device* and shared by every CMSSW stream on it, as in the
// original producer: one copy of a 2 MB model in the L2 working set rather than
// one per stream.

#include <Eigen/Core>
#include <alpaka/alpaka.hpp>

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <utility>
#include <vector>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/Utilities/interface/FileInPath.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/FixedQueueEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelLoader.h"
#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/GraphInferenceAlgo.h"
#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/TransformKernel.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace {

    // Allocated against the *device*, not a queue.  The queue overload of
    // make_device_buffer goes through the caching allocator, and this registry
    // is deliberately never freed -- so those blocks would still be live when
    // the caching allocator is destroyed, and its `liveBlocks_.empty()`
    // assertion aborts the job at exit.  The device overload is a plain
    // allocation, which is what a leak wants.
    template <typename T>
    using DeviceArray = decltype(cms::alpakatools::make_device_buffer<T[]>(
        std::declval<Device const&>(), std::declval<std::size_t>()));

    // The four sections the kernel reads, uploaded once per device.
    struct ResidentGraph {
      DeviceArray<nn::OpRecord> ops;
      DeviceArray<nn::BufferRecord> buffers;
      DeviceArray<std::int8_t> weights;
      DeviceArray<float> params;
      nn::ModelView view;
      // Whether the model's activations fit the per-warp arena, and how much
      // of the arena they actually use -- for choosing, and explaining, the
      // layout.
      bool fitsWarpLayout = false;
      int floatsUsed = 0;
    };

    struct GraphRegistryEntry {
      Device device;
      std::string file;
      std::shared_ptr<ResidentGraph> graph;
    };

    std::mutex& graphRegistryMutex() {
      static std::mutex mutex;
      return mutex;
    }

    // Deliberately leaked, and device-scoped rather than queue-scoped: CMSSW
    // destroys the CUDA context before function-local statics, so a normal
    // destructor would call into a dead context at the end of every job.
    std::vector<GraphRegistryEntry>& graphRegistry() {
      static std::vector<GraphRegistryEntry>* registry = new std::vector<GraphRegistryEntry>();
      return *registry;
    }

    std::vector<char> readWholeFile(edm::FileInPath const& file) {
      std::ifstream input(file.fullPath(), std::ios::binary | std::ios::ate);
      if (!input)
        throw cms::Exception("GraphModel") << "Cannot open model file " << file.fullPath();
      auto const bytes = static_cast<std::size_t>(input.tellg());
      input.seekg(0);
      std::vector<char> image(bytes);
      input.read(image.data(), static_cast<std::streamsize>(bytes));
      if (!input)
        throw cms::Exception("GraphModel") << "Cannot read model file " << file.fullPath();
      return image;
    }

    std::shared_ptr<ResidentGraph> residentGraph(Queue& queue, edm::FileInPath const& file) {
      auto const device = alpaka::getDev(queue);
      std::lock_guard<std::mutex> guard(graphRegistryMutex());
      for (auto const& entry : graphRegistry())
        if (entry.device == device && entry.file == file.fullPath())
          return entry.graph;

      auto const image = readWholeFile(file);
      nn::LoadedModel loaded;
      std::string error;
      if (!nn::loadModel(image.data(), image.size(), loaded, error))
        throw cms::Exception("GraphModel")
            << "Model file " << file.fullPath() << " was rejected: " << error;

      auto const& header = *loaded.header;
      // Every model this producer drives consumes the per-particle features.
      // The four-vectors are only needed by a model with a pairwise
      // interaction embedding, so a converted MLP or DeepSets network simply
      // has no input for them and none are staged.
      if (!loaded.features.valid())
        throw cms::Exception("GraphModel")
            << "Model file " << file.fullPath()
            << " declares no feature input, so there is nothing for this producer to feed it";

      // Staging buffers, so the sections reach the device in one copy each.
      auto hostOps = cms::alpakatools::make_host_buffer<nn::OpRecord[]>(queue, header.opCount);
      auto hostBuffers =
          cms::alpakatools::make_host_buffer<nn::BufferRecord[]>(queue, header.bufferCount);
      auto hostWeights =
          cms::alpakatools::make_host_buffer<std::int8_t[]>(queue, header.weightsCount);
      auto hostParams = cms::alpakatools::make_host_buffer<float[]>(queue, header.paramsCount);
      std::copy_n(loaded.ops, header.opCount, hostOps.data());
      std::copy_n(loaded.buffers, header.bufferCount, hostBuffers.data());
      std::copy_n(loaded.weights, header.weightsCount, hostWeights.data());
      std::copy_n(loaded.params, header.paramsCount, hostParams.data());

      auto graph = std::make_shared<ResidentGraph>(ResidentGraph{
          cms::alpakatools::make_device_buffer<nn::OpRecord[]>(device, header.opCount),
          cms::alpakatools::make_device_buffer<nn::BufferRecord[]>(device, header.bufferCount),
          cms::alpakatools::make_device_buffer<std::int8_t[]>(device, header.weightsCount),
          cms::alpakatools::make_device_buffer<float[]>(device, header.paramsCount),
          nn::ModelView{},
          loaded.fitsWarpMode(),
          loaded.floatsUsed});
      alpaka::memcpy(queue, graph->ops, hostOps);
      alpaka::memcpy(queue, graph->buffers, hostBuffers);
      alpaka::memcpy(queue, graph->weights, hostWeights);
      alpaka::memcpy(queue, graph->params, hostParams);
      // Complete the upload while still holding the lock: another stream that
      // finds this entry must be able to launch against it immediately, and
      // the host staging buffers have to stay alive until the copy is done.
      alpaka::wait(queue);

      graph->view = loaded.view(graph->weights.data(), graph->params.data(), graph->ops.data(),
                                graph->buffers.data());

      graphRegistry().push_back(GraphRegistryEntry{device, file.fullPath(), graph});
      return graph;
    }

  }  // namespace

  class GraphInferenceProducer : public stream::FixedQueueEDProducer<> {
  public:
    explicit GraphInferenceProducer(edm::ParameterSet const& params)
        : FixedQueueEDProducer<>(params),
          pfToken_(consumes(params.getParameter<edm::InputTag>("pf"))),
          associationMapToken_(consumes(params.getParameter<edm::InputTag>("clusters"))),
          jetBxLookupToken_(consumes(params.getParameter<edm::InputTag>("jetBxLookup"))),
          vertexToken_(consumes(params.getParameter<edm::InputTag>("vertices"))),
          vertexBxLookupToken_(consumes(params.getParameter<edm::InputTag>("vertexBxLookup"))),
          modelFile_(params.getParameter<std::string>("modelFile")),
          maxBlocks_(params.getParameter<unsigned int>("maxBlocks")),
          layoutRequest_(params.getParameter<std::string>("layout")),
          synchronizeForTiming_(params.getParameter<bool>("synchronizeForTiming")),
          outputToken_(produces()) {
      if (layoutRequest_ != "auto" && layoutRequest_ != "block" && layoutRequest_ != "warp")
        throw cms::Exception("Configuration")
            << "layout must be \"auto\", \"block\" or \"warp\", not \"" << layoutRequest_ << "\"";
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      desc.add<edm::InputTag>("pf");
      desc.add<edm::InputTag>("clusters");
      desc.add<edm::InputTag>("jetBxLookup");
      desc.add<edm::InputTag>("vertices");
      desc.add<edm::InputTag>("vertexBxLookup");
      // The whole network: op table, buffer layout and weights.  Produced by
      // scripts/export_model.py.
      desc.add<std::string>("modelFile", "PhysicsTools/ParticleTransformerAlpaka/data/b_kinadd.nnm");
      // Zero launches as many blocks as the jets need, which is what saturates
      // the device.
      desc.add<unsigned int>("maxBlocks", 0);
      // How jets map onto the device.  "block" runs one jet per 256-lane block;
      // "warp" runs one per warp, eight per block, which is several times
      // faster for a small network such as a DeepSets.  "auto" picks warp
      // whenever the model's activations fit the per-warp arena and the
      // back-end supports it (CUDA), and block otherwise -- so a ParT, which
      // needs nearly the whole arena, always runs one jet per block.
      desc.add<std::string>("layout", "auto");
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

      if (!graph_) {
        graph_ = residentGraph(queue, modelFile_);
        layout_ = chooseLayout(*graph_);
        // Say what "auto" decided: the two layouts differ several-fold in
        // speed for a small network, and nothing else in the log shows it.
        edm::LogInfo("GraphInferenceProducer")
            << modelFile_.relativePath() << " uses " << graph_->floatsUsed
            << " arena floats; running one jet per "
            << (layout_ == GraphInferenceAlgo::Layout::Warp ? "warp" : "block") << " (layout = \""
            << layoutRequest_ << "\")";
      }

      auto inputs = part::makeParticleTransformerInputs(
          queue, pf, associationMap, jetBxLookup, vertices, vertexBxLookup);

      auto const numberOfOffsets = associationMap.const_view().offset().metadata().size();
      auto const numberOfJets = numberOfOffsets == 0 ? 0 : numberOfOffsets - 1;
      l1sc::SoftJetOutputDeviceTensor outputs(queue, numberOfJets);
      outputs.zeroInitialise(queue);
      if (numberOfJets != 0) {
        algo_.run(queue, graph_->view, inputs, outputs, static_cast<uint32_t>(numberOfJets),
                  maxBlocks_, layout_);

        // alpaka::wait(queue);
      }
      event.emplace(outputToken_, std::move(outputs));
    }

  private:
    device::EDGetToken<l1sc::PFCandidateDeviceCollection> const pfToken_;
    device::EDGetToken<l1sc::AssociationMapDevice> const associationMapToken_;
    device::EDGetToken<l1sc::BxLookupDevice> const jetBxLookupToken_;
    device::EDGetToken<l1sc::VertexDeviceCollection> const vertexToken_;
    device::EDGetToken<l1sc::BxLookupDevice> const vertexBxLookupToken_;
    edm::FileInPath const modelFile_;
    uint32_t const maxBlocks_;
    std::string const layoutRequest_;
    bool const synchronizeForTiming_;
    device::EDPutToken<l1sc::SoftJetOutputDeviceTensor> const outputToken_;

    std::shared_ptr<ResidentGraph> graph_;
    GraphInferenceAlgo::Layout layout_ = GraphInferenceAlgo::Layout::Block;
    GraphInferenceAlgo algo_;

    // "auto" is warp whenever it can be; an explicit "warp" that cannot be
    // honoured is a configuration error, reported rather than quietly
    // downgraded, since the person asking for it is measuring something.
    GraphInferenceAlgo::Layout chooseLayout(ResidentGraph const& graph) const {
      using Layout = GraphInferenceAlgo::Layout;
      bool const available = GraphInferenceAlgo::kWarpLayoutAvailable;
      bool const fits = graph.fitsWarpLayout;

      if (layoutRequest_ == "block")
        return Layout::Block;
      if (layoutRequest_ == "warp") {
        if (!available)
          throw cms::Exception("Configuration")
              << "layout = \"warp\" is only implemented on the CUDA back-end";
        if (!fits)
          throw cms::Exception("Configuration")
              << "layout = \"warp\" needs a model whose activations fit " << nn::kWarpArenaFloats
              << " floats and " << nn::kWarpArenaInt8 << " int8 bytes per jet; "
              << modelFile_.fullPath() << " uses " << graph.floatsUsed
              << " floats.  Use \"block\" or \"auto\".";
        return Layout::Warp;
      }
      return available && fits ? Layout::Warp : Layout::Block;
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

DEFINE_FWK_ALPAKA_MODULE(GraphInferenceProducer);
