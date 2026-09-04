#ifndef PhysicsTools_ParticleTransformerAlpaka_TransformKernel_h
#define PhysicsTools_ParticleTransformerAlpaka_TransformKernel_h

#include <Eigen/Core>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::part {

  // Build the four Weaver inputs used by b_kinadd.yaml. The returned tensor has
  // one record per association-map entry and fixed shapes [16,2], [16,15],
  // [16,4], and [16] for points, features, vectors, and mask respectively.
  l1sc::SoftJetInputDeviceTensor makeParticleTransformerInputs(
      Queue& queue,
      l1sc::PFCandidateDeviceCollection const& pf,
      l1sc::AssociationMapDevice const& associationMap,
      l1sc::BxLookupDevice const& jetBxLookup,
      l1sc::VertexDeviceCollection const& vertices,
      l1sc::BxLookupDevice const& vertexBxLookup);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::part

#endif
