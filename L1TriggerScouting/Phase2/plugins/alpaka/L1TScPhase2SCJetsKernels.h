// header file = public interface for kernel helper class; declares what functions exist (and what inputs/outputs of those, what namespace they live)
// include guard: prevent header being included multiple times
#ifndef L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h
#define L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h

#include <alpaka/alpaka.hpp> // main Alpaka library header
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h" // jet to const mapping
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDeviceCollection.h" // offset for BX
#include "DataFormats/L1ScoutingSoA/interface/alpaka/ClustersDeviceCollection.h" // per-part cluster  and isSeed
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PuppiDeviceCollection.h" // per-part pt,eta,... on device
#include "DataFormats/L1ScoutingSoA/interface/alpaka/CounterDevice.h" // global counts (njets, clustered part) on device memory
#include "HeterogeneousCore/AlpakaInterface/interface/config.h" // Aliases like Queue, Acc1D,...

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels { // backend-specific (Alpaka) namespace

  using namespace ::l1sc;

  // algorithm dispatcher API used by producer; wrapper with methods LinkTree,...
  class L1TScPhase2SCJetsKernels { 
  public:
    //constructor
    explicit L1TScPhase2SCJetsKernels(); 
    // short name for common return type off allgorithms (BX lookup/offsets: where each BX's jets are in flat array; final jet collection; const-jet-map)
    typedef std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice> return_type;

    // SC algorithm
    return_type run(Queue& queue,
                    const PuppiDeviceCollection& src,
                    const BxLookupDeviceCollection& bx_lookup,
                    float R2,
                    ClustersDeviceCollection& clusters) const;

    // iterative
    return_type run(Queue& queue,
                    const PuppiDeviceCollection& src,
                    const BxLookupDeviceCollection& bx_lookup,
                    float R2,
                    unsigned int nJets,
                    ClustersDeviceCollection& clusters) const;

    return_type runSeededConeNMSWeighted(Queue& queue,
                                         const PuppiDeviceCollection& src,
                                         const BxLookupDeviceCollection& bxLookup,
                                         float RSeed2,
                                         float RCen2,
                                         float RClu2,
                                         float alphaSeed,
                                         float minSeedPt,
                                         unsigned int nCentroidIters,
                                         ClustersDeviceCollection& clusters) const;

    return_type runLinkTree(Queue& queue,
                            const PuppiDeviceCollection& src,
                            const BxLookupDeviceCollection& bxLookup,
                            float RLink2,
                            float ptMin,
                            ClustersDeviceCollection& clusters) const;

    return_type finalize(Queue& queue,
                         const PuppiDeviceCollection& src,
                         const BxLookupDeviceCollection& bx_lookup,
                         const ClustersDeviceCollection& clusters,
                         const CounterDevice& nJetsTotalDevice,
                         const ClusterObjDeviceCollection& jetsNonZS,
                         BxLookupDeviceCollection& jetBxLookup) const;
  };
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_Phase2_plugins_alpaka_L1TScPhase2SCJetsKernels_h