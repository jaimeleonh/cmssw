#ifndef PhysicsTools_ParticleTransformerAlpakaFP32_ParticleTransformerFP32Algo_h
#define PhysicsTools_ParticleTransformerAlpakaFP32_ParticleTransformerFP32Algo_h

#include <Eigen/Core>

#include <cstdint>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#define PARTICLE_TRANSFORMER_FP32_ALGO_API_REVISION 5

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class ParticleTransformerFP32Algo {
  public:
    // Evaluates every jet of the event in a single kernel launch.  `weights`
    // must point at the transposed ([in][out]) image of the model produced by
    // partfp32::transposeWeights.  `maxBlocks` caps the grid; zero means one
    // block per jet, which is the fastest configuration on every device tested.
    void run(Queue& queue,
             float const* weights,
             float const* params,
             l1sc::SoftJetInputDeviceTensor const& inputs,
             l1sc::SoftJetOutputDeviceTensor& outputs,
             uint32_t numberOfJets,
             uint32_t maxBlocks) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
