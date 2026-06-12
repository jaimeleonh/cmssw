#ifndef L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernel_h
#define L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernel_h

#include <alpaka/alpaka.hpp>

#include <cstdio>
#include <limits>
#include <fmt/core.h> 

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/Portable/interface/PortableHostObject.h"
#include "DataFormats/Portable/interface/alpaka/PortableObject.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// // These definitions are not stored in DataFormats/L1ScoutingSoA/
// // since are designed to be used as helper types in the kernels
// // and not meant to be stored in the FW event struct
namespace l1sc {

  struct PortableCounter {
    int value;
  };
  using CounterHost = PortableHostObject<PortableCounter>;

}  // namespace l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;
  using CounterDevice = PortableObject<PortableCounter>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace ::l1sc;

  // SoftJetInputDeviceTensor transform(Queue& queue, 
  //                const PFCandidateDeviceCollection& pf, 
  //                const BxLookupDeviceCollection& bx_lookup, 
  //                const ClustersDeviceCollection& clusters);
  // SoftJetInputDeviceTensor transform(Queue& queue, 
  //                const PFCandidateDeviceCollection& pf, 
  //                const ClustersDeviceCollection& clusters);
  SoftJetInputDeviceTensor transform(Queue& queue, 
                 const PFCandidateDeviceCollection& pf, 
                 const AssociationMapDevice& association_map,
                 const BxLookupDevice& jetBxLookup,
                 const VertexDeviceCollection& vertices,
                 const BxLookupDevice& vertexBxLookup);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernel_h