#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/DeepSetsAlgo.h"

#if !defined(DEEPSETS_ALGO_API_REVISION) || DEEPSETS_ALGO_API_REVISION != 1
#error "DeepSetsAlgo.h is stale: replace the complete ParticleTransformerAlpaka package before rebuilding"
#endif

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/DeepSetsFastKernel.h"

#include <algorithm>
#include <type_traits>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace {

    // One block evaluates one jet.  The grid-stride loop only matters when the
    // caller caps the grid; with the default configuration each block runs a
    // single iteration.
    class DeepSetsInferenceKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::part::generated::ModelView model,
                                    l1sc::SoftJetInputDeviceTensor::ConstView inputs,
                                    l1sc::SoftJetOutputDeviceTensor::View outputs,
                                    uint32_t numberOfJets) const {
        auto& shared = alpaka::declareSharedVar<::part::deepsets::DeepSetsShared, __COUNTER__>(acc);
        auto const lane = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const lanes = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const block = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        constexpr uint32_t particles = ::part::deepsets::kMaxParticles;
        constexpr uint32_t features = ::part::generated::input_features;

        for (uint32_t jet = block; jet < numberOfJets; jet += blocks) {
          auto const input = inputs[jet];

          // The Eigen-backed SoA must be addressed by matrix coordinates; it
          // cannot be reinterpreted as the kernel's row-major arrays.  The
          // DeepSets graph reads no four-vector input.
          for (uint32_t index = lane; index < particles * features; index += lanes) {
            auto const particle = index / features;
            shared.features[index] = input.features()(particle, index - particle * features);
          }
          for (uint32_t particle = lane; particle < particles; particle += lanes)
            shared.mask[particle] = input.mask()(particle);
          alpaka::syncBlockThreads(acc);

          ::part::deepsets::runDeepSets(acc, model, shared);

          if (lane == 0)
            outputs.output()[jet] = shared.logits[0];  // class 0 of the two-class head
          alpaka::syncBlockThreads(acc);
        }
      }
    };

#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
    // One jet per warp, part::deepsets::kWarpJetsPerBlock jets per block.  The
    // graph is ~25 kMACs per jet, so a whole block per jet spends its time on
    // barrier latency; one warp per jet pays a warp-level sync instead and
    // hides every staging load behind the other warps' computations -- the
    // same shape the graph interpreter uses for small networks.  The graph
    // code itself is exactly runDeepSets's, reached through WarpAcc.
    class DeepSetsWarpInferenceKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::part::generated::ModelView model,
                                    l1sc::SoftJetInputDeviceTensor::ConstView inputs,
                                    l1sc::SoftJetOutputDeviceTensor::View outputs,
                                    uint32_t numberOfJets) const {
        auto& block = alpaka::declareSharedVar<::part::deepsets::DeepSetsWarpBlock, __COUNTER__>(acc);
        auto const thread = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const blockIndex = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        using namespace ::part::deepsets;
        uint32_t const warp = thread / kWarpLanes;
        uint32_t const lane = thread % kWarpLanes;
        WarpAcc<TAcc> const group{acc, static_cast<int>(lane)};
        DeepSetsShared& shared = block.jet[warp];

        constexpr uint32_t particles = kMaxParticles;
        constexpr uint32_t features = ::part::generated::input_features;
        uint32_t const stride = blocks * kWarpJetsPerBlock;

        // The condition is uniform across the warp, which is all a warp-level
        // sync needs; different warps may leave the loop at different times.
        for (uint32_t jet = blockIndex * kWarpJetsPerBlock + warp; jet < numberOfJets; jet += stride) {
          auto const input = inputs[jet];

          for (uint32_t index = lane; index < particles * features; index += kWarpLanes) {
            auto const particle = index / features;
            shared.features[index] = input.features()(particle, index - particle * features);
          }
          for (uint32_t particle = lane; particle < particles; particle += kWarpLanes)
            shared.mask[particle] = input.mask()(particle);
          groupBarrier(group);

          runDeepSets(group, model, shared);

          if (lane == 0)
            outputs.output()[jet] = shared.logits[0];  // class 0 of the two-class head
          groupBarrier(group);
        }
      }
    };
#endif

  }  // namespace

  void DeepSetsAlgo::run(Queue& queue,
                         std::int8_t const* weights,
                         float const* params,
                         l1sc::SoftJetInputDeviceTensor const& inputs,
                         l1sc::SoftJetOutputDeviceTensor& outputs,
                         uint32_t numberOfJets,
                         uint32_t maxBlocks) const {
    if (numberOfJets == 0)
      return;
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
    // Warp layout: kWarpJetsPerBlock jets per block, warp-level synchronisation.
    uint32_t const jetsPerBlock = ::part::deepsets::kWarpJetsPerBlock;
    auto const wanted = (numberOfJets + jetsPerBlock - 1) / jetsPerBlock;
    auto const blocks = maxBlocks == 0 ? wanted : std::min<uint32_t>(maxBlocks, wanted);
    auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, ::part::deepsets::kThreadsPerJet);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        DeepSetsWarpInferenceKernel{},
                        ::part::generated::ModelView{weights, params},
                        inputs.const_view(),
                        outputs.view(),
                        numberOfJets);
#else
    auto const blocks = maxBlocks == 0 ? numberOfJets : std::min(maxBlocks, numberOfJets);
    auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, ::part::deepsets::kThreadsPerJet);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        DeepSetsInferenceKernel{},
                        ::part::generated::ModelView{weights, params},
                        inputs.const_view(),
                        outputs.view(),
                        numberOfJets);
#endif
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
