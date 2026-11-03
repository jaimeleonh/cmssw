#ifndef L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernelTorchScript_h
#define L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernelTorchScript_h

#include <alpaka/alpaka.hpp>

#include <Eigen/Core>

#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels::torchscript {

  SoftJetNewInputDeviceTensor transform(
      Queue& queue,
      const PFCandidateDeviceCollection& pf,
      const ClusterObjDeviceCollection& jets,
      const AssociationMapDevice& association_map,
      const BxLookupDevice& jetBxLookup,
      const VertexDeviceCollection& vertices,
      const BxLookupDevice& vertexBxLookup);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels::torchscript

#endif  // L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernelTorchScript_h
