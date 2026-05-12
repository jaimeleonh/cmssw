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

  // algorithm dispatcher API used by producer; wrapper with methods SCGreedy, SCNMS, ...
  class L1TScPhase2SCJetsKernels {
  public:
    // constructor
    explicit L1TScPhase2SCJetsKernels();

    // short name for common return type of all algorithms
    // (BX lookup/offsets: where each BX's jets are in flat array; final jet collection; const-jet-map)
    typedef std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice> return_type;

    // ------------------------------------------------------------------
    // legacy single-radius non-iterative seeded cone
    // kept mainly for backward compatibility / auto mode
    // (closest-axis assignment after seed finding)
    // ------------------------------------------------------------------
    return_type run(Queue& queue,
                    const PuppiDeviceCollection& src,
                    const BxLookupDeviceCollection& bx_lookup,
                    float R2,
                    ClustersDeviceCollection& clusters) const;

    // ------------------------------------------------------------------
    // SCGreedy = iterative seeded cone (greedy / L1T-style)
    // ------------------------------------------------------------------
    return_type run(Queue& queue,
                    const PuppiDeviceCollection& src,
                    const BxLookupDeviceCollection& bx_lookup,
                    float R2,
                    unsigned int nJets,
                    ClustersDeviceCollection& clusters) const;

    // ------------------------------------------------------------------
    // SCNMS = old non-iterative seeded cone with split radii
    //
    // This preserves the old kernel logic exactly, except that:
    //   RSeed controls seed finding + old centroid/precluster accumulation
    //   RClu  controls final particle assignment to accepted seed axes
    //
    // RCen is intentionally not used here; it belongs only to the newer
    // multi-iteration NMS-style implementation.
    // ------------------------------------------------------------------
    return_type runSCNMS(Queue& queue,
                         const PuppiDeviceCollection& src,
                         const BxLookupDeviceCollection& bxLookup,
                         float RSeed2,
                         float RClu2,
                         ClustersDeviceCollection& clusters) const;

    // ------------------------------------------------------------------
    // SCNMSWeighted = old weighted non-iterative seeded cone with split radii
    //
    // Same as SCNMS above, except final assignment uses the old weighted
    // metric:
    //   metric = dr2 / pt_jet^2
    //
    // alphaSeed and RCen are intentionally not used here; this preserves
    // the old weighted kernel logic while allowing different seed and
    // assignment radii.
    // ------------------------------------------------------------------
    return_type runSCNMSWeighted(Queue& queue,
                                 const PuppiDeviceCollection& src,
                                 const BxLookupDeviceCollection& bxLookup,
                                 float RSeed2,
                                 float RClu2,
                                 ClustersDeviceCollection& clusters) const;

    // ------------------------------------------------------------------
    // SCNMSWeightedMultiIter = newer NMS-style weighted seeded cone
    // explicit radii for seed finding / centroiding / assignment, plus
    // configurable centroid iterations and weighted assignment exponent
    // ------------------------------------------------------------------
    return_type runSCNMSWeightedMultiIter(Queue& queue,
                                          const PuppiDeviceCollection& src,
                                          const BxLookupDeviceCollection& bxLookup,
                                          float RSeed2,
                                          float RCen2,
                                          float RClu2,
                                          float alphaSeed,
                                          float minSeedPt,
                                          unsigned int nCentroidIters,
                                          ClustersDeviceCollection& clusters) const;

    // ------------------------------------------------------------------
    // LinkTree
    // ------------------------------------------------------------------
    return_type runLinkTree(Queue& queue,
                            const PuppiDeviceCollection& src,
                            const BxLookupDeviceCollection& bxLookup,
                            float RLink2,
                            float ptMin,
                            ClustersDeviceCollection& clusters) const;

    // common post-processing:
    // compact jets, build per-BX offsets, and create jet->const association map
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
