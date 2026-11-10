#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h
#include <alpaka/alpaka.hpp>
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/CounterDevice.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"
namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {
  using namespace ::l1sc;
  struct SCNMSWCounts {
    unsigned int jets;
    unsigned int clustered;
  };
  class L1TScPhase2SCJetsKernels {
  public:
    explicit L1TScPhase2SCJetsKernels(unsigned int minConstituents = 0) : minConstituents_{minConstituents} {}
    using return_type = std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice>;
    return_type runSeededCone(Queue&,
                              const PuppiDeviceCollection&,
                              const BxLookupDevice&,
                              float R2,
                              unsigned int nJets,
                              ClustersDeviceCollection&) const;
    return_type runSCNMSW(Queue&,
                          const PuppiDeviceCollection&,
                          const BxLookupDevice&,
                          float RSeed2,
                          float RCen2,
                          float RClu2,
                          float alpha,
                          int maxRounds,
                          ClustersDeviceCollection&) const;
    return_type runLinkTree(Queue&,
                            const PuppiDeviceCollection&,
                            const BxLookupDevice&,
                            float RLink2,
                            float ptMin,
                            ClustersDeviceCollection&) const;
    // Output helpers are public because NVCC extended lambdas require it.
    return_type finalize(Queue&,
                         const PuppiDeviceCollection&,
                         const BxLookupDevice&,
                         ClustersDeviceCollection&,
                         CounterDevice&,
                         ClusterObjDeviceCollection&,
                         BxLookupDevice&) const;
    return_type finalizeSCNMSW(Queue&,
                               const PuppiDeviceCollection&,
                               const BxLookupDevice&,
                               const ClustersDeviceCollection&,
                               const PortableObject<SCNMSWCounts>&,
                               const ClusterObjDeviceCollection&,
                               BxLookupDevice&,
                               uint32_t*,
                               uint32_t*) const;
  private:
    unsigned int minConstituents_;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
#endif
