#ifndef PhysicsTools_ParticleTransformerAlpaka_ParticleTransformerAlgo_h
#define PhysicsTools_ParticleTransformerAlpaka_ParticleTransformerAlgo_h

#include <Eigen/Core>

#include <cstdint>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#define PARTICLE_TRANSFORMER_ALGO_API_REVISION 5

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class ParticleTransformerAlgo {
  public:
    // Evaluates every jet of the event in a single kernel launch.  `weights`
    // must point at the transposed ([in][out]) image of the model produced by
    // part::transposeWeights.  `maxBlocks` caps the grid; zero means one
    // block per jet, which is the fastest configuration on every device tested.
    void run(Queue& queue,
             std::int8_t const* weights,
             float const* params,
             l1sc::SoftJetInputDeviceTensor const& inputs,
             l1sc::SoftJetOutputDeviceTensor& outputs,
             uint32_t numberOfJets,
             uint32_t maxBlocks) const;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

#endif
