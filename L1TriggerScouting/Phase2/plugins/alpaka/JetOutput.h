#ifndef L1TriggerScouting_Phase2_plugins_alpaka_JetOutput_h
#define L1TriggerScouting_Phase2_plugins_alpaka_JetOutput_h

#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include <cassert>
#include <type_traits>
#include "L1TriggerScouting/Phase2/plugins/alpaka/radixSort64.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {
  using namespace cms::alpakatools;
  class JetMinConstituentsFilterKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  BxLookupSoA::ConstView bxLookup,        // particle ranges per BX
                                  ClusterObjDeviceCollection::View jets,  // proto-jets (compacted in place)
                                  ClustersDeviceCollection::View clusters,  // per-particle assignment (remapped)
                                  BxLookupSoA::View jetBxLookup,          // per-BX jet counts (updated)
                                  unsigned int minConstituents,
                                  int32_t* remap,  // scratch: old local id -> new local id (per BX slot)
                                  unsigned int* nJetsTotal) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
        if (end <= begin)
          continue;

        // number of proto-jets in this BX (counts, before the prefix scan in finalize)
        uint32_t njets = jetBxLookup.offset()[block_idx + 1].offset();

        // serial pass over the proto-jets of this BX:
        // build the old -> new local id map and compact the kept jets in place
        if (once_per_block(acc)) {
          uint32_t kept = 0;
          for (uint32_t j = 0; j < njets; ++j) {
            uint32_t slot = begin + j;
            if (jets.numberOfDaughters()[slot] >= minConstituents) {
              remap[slot] = static_cast<int32_t>(kept);
              if (kept != j) {
                jets.pt()[begin + kept] = jets.pt()[slot];
                jets.eta()[begin + kept] = jets.eta()[slot];
                jets.phi()[begin + kept] = jets.phi()[slot];
                jets.cluster()[begin + kept] = kept;  // new local jet id
                jets.numberOfDaughters()[begin + kept] = jets.numberOfDaughters()[slot];
              }
              ++kept;
            } else {
              remap[slot] = -1;
            }
          }
          jetBxLookup.offset()[block_idx + 1].offset() = kept;
          alpaka::atomicAdd(acc, nJetsTotal, kept, alpaka::hierarchy::Blocks{});
        }

        alpaka::syncBlockThreads(acc);

        // remap the per-particle assignments to the compacted local jet ids
        for (uint32_t tid : independent_group_elements(acc, end - begin)) {
          uint32_t i = tid + begin;
          int32_t jcluster = clusters.cluster()[i];
          if (jcluster < 0)
            continue;

          // -1 for dropped jets, new local id for kept jets
          clusters.cluster()[i] = remap[begin + jcluster];
          if (clusters.cluster()[i] < 0)
            clusters.is_seed()[i] = 0;
        }

        alpaka::syncBlockThreads(acc);
      }
    }
  };

  // Common jet compaction and ordered constituent-map construction.
  class JetZSKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  BxLookupSoA::ConstView bxLookup,
                                  ClusterObjDeviceCollection::ConstView jetsNonZS,
                                  BxLookupSoA::ConstView jetBxLookup,
                                  ClusterObjDeviceCollection::View jets,
                                  unsigned int* nClustered) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t beginSrc = bxLookup.offset().offset()[block_idx];
        uint32_t beginDst = jetBxLookup.offset().offset()[block_idx];
        uint32_t endDst = jetBxLookup.offset().offset()[block_idx + 1];
        if (endDst <= beginDst)
          continue;
        for (uint32_t tid : independent_group_elements(acc, endDst - beginDst)) {
          uint32_t ijetDst = tid + beginDst;
          uint32_t ijetSrc = tid + beginSrc;
          jets.pt()[ijetDst] = jetsNonZS.pt()[ijetSrc];
          jets.eta()[ijetDst] = jetsNonZS.eta()[ijetSrc];
          jets.phi()[ijetDst] = jetsNonZS.phi()[ijetSrc];
          jets.cluster()[ijetDst] = jetsNonZS.cluster()[ijetSrc];
          jets.numberOfDaughters()[ijetDst] = jetsNonZS.numberOfDaughters()[ijetSrc];
        }
        alpaka::syncBlockThreads(acc);
        if (once_per_block(acc)) {
          unsigned int nClusteredBlock = 0;
          for (uint32_t iJet = beginDst; iJet < endDst; ++iJet) {
            nClusteredBlock += jets.numberOfDaughters()[iJet];
          }
          alpaka::atomicAdd(acc, nClustered, nClusteredBlock, alpaka::hierarchy::Blocks{});
        }
      }
    }
  };

  // Build constituent maps in the original candidate order within each jet.
  class JetToAssociationMapKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView jetBxLookup,
                                  AssociationMapSoA::ConstView clusterdParticleOffsets,
                                  ClustersDeviceCollection::ConstView clusters,
                                  const unsigned int njets,
                                  const unsigned int nclustered,
                                  uint32_t* key,
                                  uint16_t* idx,
                                  AssociationMapSoA::View map) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      unsigned int nparticles = clusters.metadata().size();
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset().offset()[block_idx];
        uint32_t end = bxLookup.offset().offset()[block_idx + 1];
        if (end <= begin)
          continue;
        for (uint32_t tid : independent_group_elements(acc, end - begin)) {
          uint32_t i = tid + begin;
          assert(i < nparticles);
          assert(clusters.cluster()[i] == -1 || clusters.cluster()[i] < int(end - begin));
          if (clusters.cluster()[i] == -1) {
            key[i] = 0xFFFFFFFF;
          } else {
            uint16_t ptcode = static_cast<uint16_t>(std::max(0xFFFF - puppi.pt()[i] * 32.f, 0.0f));
            key[i] = (static_cast<uint32_t>(clusters.cluster()[i]) << 16) | ptcode;
          }
          idx[i] = tid;
#ifdef L1TSC_VERBOSE_DEBUG
          if (block_idx <= 2)
            printf("In BX %u cand %3u/%3u of pt %6.2f assigned to cluster %3d -> key %08x\n",
                   block_idx + 1,
                   tid,
                   end - begin,
                   puppi.pt()[i],
                   clusters.cluster()[i],
                   unsigned(key[i]));
#endif
        }
      }
      alpaka::syncBlockThreads(acc);

      // sort
      radixSortMulti(acc, key, idx, bxLookup.offset().offset().data(), nullptr);

      // rearrange jet constituents
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset().offset()[block_idx];
        uint32_t end = bxLookup.offset().offset()[block_idx + 1];
        if (end <= begin)
          continue;

        uint32_t jetsBegin = jetBxLookup.offset().offset()[block_idx];
        uint32_t jetsEnd = jetBxLookup.offset().offset()[block_idx + 1];
        assert(jetsEnd >= jetsBegin);
        assert(jetsEnd <= njets);
        assert(jetBxLookup.offset().offset()[grid_dim] == njets);
        assert(clusterdParticleOffsets.offset().metadata().size() == int(njets + 1));
        assert(clusterdParticleOffsets.offset().offset()[njets] == nclustered);

        uint32_t clusteredBegin = clusterdParticleOffsets.offset().offset()[jetsBegin];
        uint32_t clusteredEnd = clusterdParticleOffsets.offset().offset()[jetsEnd];

#ifdef L1TSC_VERBOSE_DEBUG
        if (once_per_block(acc) && block_idx <= 2000)
          printf("In BX %u, particle offsets %u - %u (%u), jets %u - %u (%u), clustered particles %u - %u (%u) \n",
                 block_idx + 1,
                 begin,
                 end,
                 end - begin,
                 jetsBegin,
                 jetsEnd,
                 jetsEnd - jetsBegin,
                 clusteredBegin,
                 clusteredEnd,
                 clusteredEnd - clusteredBegin);
#endif
        assert(clusteredEnd >= clusteredBegin);
        assert((clusteredEnd - clusteredBegin) <= (end - begin));

        for (uint32_t tid : independent_group_elements(acc, clusteredEnd - clusteredBegin)) {
          uint32_t i = tid + begin;
          uint32_t icand = idx[i] + begin;
          map.index().index()[clusteredBegin + tid] = icand;
          assert(clusteredBegin + tid < nclustered);
          assert(map.index().metadata().size() == int(nclustered));
#ifdef L1TSC_VERBOSE_DEBUG
          uint32_t ijet = key[icand] >> 16;
          if (block_idx <= 2)
            printf("In BX %u index %3u/%3u is cand %3u/%3u of index %6u assigned to jet %3u (local)\n",
                   block_idx + 1,
                   tid,
                   clusteredEnd - clusteredBegin,
                   idx[i],
                   end - begin,
                   icand,
                   ijet);
#endif
        }
      }
    }
  };

  template <unsigned int Unroll>
  class JetAssociationRankKernel {
  public:
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView jetBxLookup,
                                  AssociationMapSoA::ConstView particleOffsets,
                                  ClustersDeviceCollection::ConstView clusters,
                                  uint32_t* scratch,
                                  AssociationMapSoA::View map) const {
      auto& keys = alpaka::declareSharedVar<uint32_t[128], __COUNTER__>(acc);
      uint32_t grid = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t bx : independent_groups(acc, grid)) {
        uint32_t begin = bxLookup.offset().offset()[bx], size = bxLookup.offset().offset()[bx + 1] - begin;
        uint32_t* key = size <= 128 ? keys : scratch + begin;
        for (uint32_t tid : independent_group_elements(acc, size)) {
          uint32_t i = begin + tid;
          int cid = clusters.cluster()[i];
          uint16_t code = static_cast<uint16_t>(std::max(0xFFFF - puppi.pt()[i] * 32.f, 0.f));
          key[tid] = cid < 0 ? 0xFFFFFFFFu : (static_cast<uint32_t>(cid) << 16) | code;
        }
        alpaka::syncBlockThreads(acc);
        uint32_t base = particleOffsets.offset().offset()[jetBxLookup.offset().offset()[bx]];
        for (uint32_t tid : independent_group_elements(acc, size)) {
          if (clusters.cluster()[begin + tid] < 0)
            continue;
          uint32_t own = key[tid], rank = 0;
          for (uint32_t j = 0; j < size; j += Unroll) {
            for (uint32_t k = 0; k < Unroll; ++k) {
              if (j + k >= size)
                break;
              uint32_t other = key[j + k];
              rank += other < own || (other == own && j + k < tid);
            }
          }
          map.index().index()[base + rank] = begin + tid;
        }
        alpaka::syncBlockThreads(acc);
      }
    }
  };

  // Compact jets and construct their constituent offsets in one kernel.
  class JetCompactOffsetsKernel {
  public:
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  BxLookupSoA::ConstView bxLookup,
                                  ClusterObjDeviceCollection::ConstView source,
                                  BxLookupSoA::ConstView jetBxLookup,
                                  uint32_t const* particleBxOffsets,
                                  ClusterObjDeviceCollection::View jets,
                                  AssociationMapSoA::View mapOffsets) const {
      uint32_t grid = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      if (once_per_grid(acc))
        mapOffsets.offset().offset()[mapOffsets.offset().metadata().size() - 1] = particleBxOffsets[grid];
      for (uint32_t bx : independent_groups(acc, grid)) {
        uint32_t from = bxLookup.offset().offset()[bx], to = jetBxLookup.offset().offset()[bx];
        uint32_t n = jetBxLookup.offset().offset()[bx + 1] - to;
        for (uint32_t j : independent_group_elements(acc, n)) {
          jets.pt()[to + j] = source.pt()[from + j];
          jets.eta()[to + j] = source.eta()[from + j];
          jets.phi()[to + j] = source.phi()[from + j];
          jets.cluster()[to + j] = source.cluster()[from + j];
          jets.numberOfDaughters()[to + j] = source.numberOfDaughters()[from + j];
          uint32_t offset = particleBxOffsets[bx];
          for (uint32_t k = 0; k < j; ++k)
            offset += source.numberOfDaughters()[from + k];
          mapOffsets.offset().offset()[to + j] = offset;
        }
      }
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif
