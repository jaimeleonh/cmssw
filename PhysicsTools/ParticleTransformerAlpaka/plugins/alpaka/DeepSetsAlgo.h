#ifndef PhysicsTools_ParticleTransformerAlpaka_DeepSetsAlgo_h
#define PhysicsTools_ParticleTransformerAlpaka_DeepSetsAlgo_h

#include <Eigen/Core>

#include <cstdint>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

#define DEEPSETS_ALGO_API_REVISION 1

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  class DeepSetsAlgo {
  public:
    // Evaluates every jet of the event in a single kernel launch.  Unlike the
    // ParticleTransformer path, `weights` must point at the exported row-major
    // [out][in] image unchanged: the DeepSets kernel reads the 2048 INT8
    // values directly, with no transposed copy.  `maxBlocks` caps the grid;
    // zero means one block per jet.
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
