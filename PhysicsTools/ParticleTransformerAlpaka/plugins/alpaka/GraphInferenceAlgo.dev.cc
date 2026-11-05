#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/GraphInferenceAlgo.h"

#if !defined(GRAPH_INFERENCE_ALGO_API_REVISION) || GRAPH_INFERENCE_ALGO_API_REVISION != 2
#error "GraphInferenceAlgo.h is stale: replace the complete ParticleTransformerAlpaka package before rebuilding"
#endif

#include <algorithm>
#include <type_traits>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace {

    // One block evaluates one jet.  The grid-stride loop only matters when the
    // caller caps the grid; with the default configuration each block runs a
    // single iteration.
    class GraphJetKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    nn::ModelView model,
                                    l1sc::SoftJetInputDeviceTensor::ConstView inputs,
                                    l1sc::SoftJetOutputDeviceTensor::View outputs,
                                    uint32_t numberOfJets) const {
        auto& shared = alpaka::declareSharedVar<nn::GraphShared, __COUNTER__>(acc);
        auto const lane = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const lanes = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const block = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        auto const particles = static_cast<uint32_t>(model.maxTokens);
        auto const features = static_cast<uint32_t>(model.features.cols);
        // A model without a pairwise embedding has no four-vector input.
        auto const components =
            model.vectors.valid() ? static_cast<uint32_t>(model.vectors.cols) : 0u;

        for (uint32_t jet = block; jet < numberOfJets; jet += blocks) {
          auto const input = inputs[jet];

          // Written straight into the arena rather than through a staging
          // buffer: the arena is what the model was packed for, and a second
          // copy of the inputs would cost shared memory it needs.  runGraph is
          // then told there is nothing left to load.  The Eigen-backed SoA has
          // to be addressed by matrix coordinates; it cannot be reinterpreted
          // as the kernel's row-major arrays.
          for (uint32_t index = lane; index < particles * features; index += lanes) {
            auto const particle = index / features;
            auto const feature = index - particle * features;
            shared.arena[model.features.offset + particle * model.features.stride + feature] =
                input.features()(particle, feature);
          }
          for (uint32_t index = lane; index < particles * components; index += lanes) {
            auto const particle = index / components;
            auto const component = index - particle * components;
            shared.arena[model.vectors.offset + particle * model.vectors.stride + component] =
                input.vectors()(particle, component);
          }
          for (uint32_t particle = lane; particle < particles; particle += lanes)
            shared.mask[particle] = input.mask()(particle);
          alpaka::syncBlockThreads(acc);

          nn::runGraph(acc, model, shared, nn::GraphInputs{nullptr, nullptr});

          if (lane == 0)
            outputs.output()[jet] = shared.output[0];  // YAML class jet_isB
          alpaka::syncBlockThreads(acc);
        }
      }
    };

#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
    // One jet per warp, nn::kWarpJetsPerBlock jets per block.  The op code is
    // exactly the block kernel's, reached through nn::SubGroup: each warp sees
    // itself as a 32-lane "block" with its own small arena and its own
    // warp-level barrier, so the warps of a block never wait for each other.
    // That independence is the point -- the jets in a block have different
    // particle counts, and so different control flow.
    class GraphWarpJetKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    nn::ModelView model,
                                    l1sc::SoftJetInputDeviceTensor::ConstView inputs,
                                    l1sc::SoftJetOutputDeviceTensor::View outputs,
                                    uint32_t numberOfJets) const {
        auto& block = alpaka::declareSharedVar<nn::GraphWarpBlock, __COUNTER__>(acc);
        auto const thread = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const blockIndex = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        uint32_t const warp = thread / nn::kWarpLanes;
        uint32_t const lane = thread % nn::kWarpLanes;
        nn::SubGroup<TAcc> const group{acc, static_cast<int>(lane), nn::kWarpLanes};
        nn::GraphWarpState& shared = block.jet[warp];

        auto const particles = static_cast<uint32_t>(model.maxTokens);
        auto const features = static_cast<uint32_t>(model.features.cols);
        auto const components =
            model.vectors.valid() ? static_cast<uint32_t>(model.vectors.cols) : 0u;
        uint32_t const stride = blocks * nn::kWarpJetsPerBlock;

        // The condition is uniform across the warp, which is all a warp-level
        // barrier needs; different warps may leave the loop at different times.
        for (uint32_t jet = blockIndex * nn::kWarpJetsPerBlock + warp; jet < numberOfJets;
             jet += stride) {
          auto const input = inputs[jet];

          // Staged straight into this warp's arena, as the block kernel does
          // into the block's; runGraph leaves those buffers alone.
          for (uint32_t index = lane; index < particles * features; index += nn::kWarpLanes) {
            auto const particle = index / features;
            auto const feature = index - particle * features;
            shared.arena[model.features.offset + particle * model.features.stride + feature] =
                input.features()(particle, feature);
          }
          for (uint32_t index = lane; index < particles * components; index += nn::kWarpLanes) {
            auto const particle = index / components;
            auto const component = index - particle * components;
            shared.arena[model.vectors.offset + particle * model.vectors.stride + component] =
                input.vectors()(particle, component);
          }
          for (uint32_t particle = lane; particle < particles; particle += nn::kWarpLanes)
            shared.mask[particle] = input.mask()(particle);
          nn::barrier(group);

          nn::runGraph(group, model, shared, nn::GraphInputs{nullptr, nullptr});

          if (lane == 0)
            outputs.output()[jet] = shared.output[0];  // YAML class jet_isB
          nn::barrier(group);
        }
      }
    };
#endif

  }  // namespace

  void GraphInferenceAlgo::run(Queue& queue,
                               nn::ModelView const& model,
                               l1sc::SoftJetInputDeviceTensor const& inputs,
                               l1sc::SoftJetOutputDeviceTensor& outputs,
                               uint32_t numberOfJets,
                               uint32_t maxBlocks,
                               [[maybe_unused]] Layout layout) const {
    if (numberOfJets == 0)
      return;

#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
    if (layout == Layout::Warp) {
      uint32_t const needed = (numberOfJets + nn::kWarpJetsPerBlock - 1) / nn::kWarpJetsPerBlock;
      auto const blocks = maxBlocks == 0 ? needed : std::min(maxBlocks, needed);
      auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, nn::kThreadsPerJet);
      alpaka::exec<Acc1D>(queue, workDiv, GraphWarpJetKernel{}, model, inputs.const_view(),
                          outputs.view(), numberOfJets);
      return;
    }
#endif

    auto const blocks = maxBlocks == 0 ? numberOfJets : std::min(maxBlocks, numberOfJets);
    auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, nn::kThreadsPerJet);
    alpaka::exec<Acc1D>(queue, workDiv, GraphJetKernel{}, model, inputs.const_view(),
                        outputs.view(), numberOfJets);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
