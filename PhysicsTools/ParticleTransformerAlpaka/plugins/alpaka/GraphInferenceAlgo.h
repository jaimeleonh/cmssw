#ifndef PhysicsTools_ParticleTransformerAlpaka_GraphInferenceAlgo_h
#define PhysicsTools_ParticleTransformerAlpaka_GraphInferenceAlgo_h

#include <Eigen/Core>

#include <cstdint>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/GraphKernel.h"

#define GRAPH_INFERENCE_ALGO_API_REVISION 2

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  // The model-agnostic counterpart of ParticleTransformerAlgo: it takes a
  // view of an uploaded graph rather than a pair of weight blobs, so the same
  // build runs any model the exporter accepted.
  class GraphInferenceAlgo {
  public:
    // How the jets are laid out over the device.
    //   Block: one jet per 256-lane block, with the full shared arena.  What
    //          every network that needs most of the arena -- ParT, tau -- uses.
    //   Warp:  one jet per warp, eight jets per block, each with a small arena.
    //          Much faster for a small network, whose time otherwise goes on
    //          the block-wide barrier that ends every op; see GraphShared.h.
    enum class Layout { Block, Warp };

    // Warp mode is implemented for CUDA only.  Its barrier is __syncwarp, and
    // a ROCm wavefront is 64 lanes, not 32, so it would need its own group
    // size and sync; the CPU back-ends run one thread per block, where there is
    // nothing to split.  Everywhere else every model runs in block mode.
    static constexpr bool kWarpLayoutAvailable =
#if defined(ALPAKA_ACC_GPU_CUDA_ENABLED)
        true;
#else
        false;
#endif

    // Evaluates every jet of the event in a single kernel launch.  `maxBlocks`
    // caps the grid; zero means as many blocks as the jets need, which is the
    // fastest configuration on every device tested.  Asking for Layout::Warp
    // where it is not available is a logic error the producer rules out.
    void run(Queue& queue,
             nn::ModelView const& model,
             l1sc::SoftJetInputDeviceTensor const& inputs,
             l1sc::SoftJetOutputDeviceTensor& outputs,
             uint32_t numberOfJets,
             uint32_t maxBlocks,
             Layout layout = Layout::Block) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
