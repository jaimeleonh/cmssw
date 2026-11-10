#ifndef L1TriggerScouting_Phase2_plugins_alpaka_SeededCone_h
#define L1TriggerScouting_Phase2_plugins_alpaka_SeededCone_h

#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {
  using namespace cms::alpakatools;
  // Greedy seeded cone: fixed pT priority, assigned-input flags, scalar centroids.
  class SeededConeKernel {
  public:
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView bxIndex,
                                  float R2,
                                  unsigned int nIters,
                                  ClustersDeviceCollection::View clusters,
                                  uint32_t* order,
                                  ClusterObjDeviceCollection::View jets,
                                  BxLookupSoA::View jetBxLookup,
                                  BxLookupSoA::View jetBxIndex,
                                  unsigned int* nJetsTotal) const {
      if (once_per_grid(acc))
        jetBxLookup.offset().offset()[0] = 0;
      auto& cursor = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      auto& seedIndex = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      auto& sumPt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
      auto& sumEta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
      auto& sumPhi = alpaka::declareSharedVar<float, __COUNTER__>(acc);
      auto& daughters = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      uint32_t gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t bx : independent_groups(acc, gridDim)) {
        uint32_t begin = bxLookup.offset().offset()[bx], end = bxLookup.offset().offset()[bx + 1];
        uint32_t size = end - begin;
        for (uint32_t tid : independent_group_elements(acc, size)) {
          uint32_t i = begin + tid, rank = 0;
          float pt = puppi.pt()[i];
          for (uint32_t j = begin; j < end; ++j)
            rank += puppi.pt()[j] > pt || (puppi.pt()[j] == pt && j < i);
          order[begin + rank] = i;
          clusters.cluster()[i] = -1;
          clusters.is_seed()[i] = 0;
        }
        if (once_per_block(acc))
          cursor = 0;
        alpaka::syncBlockThreads(acc);
        uint32_t iter = 0;
        for (; iter < nIters; ++iter) {
          if (once_per_block(acc)) {
            seedIndex = end;
            while (cursor < size) {
              uint32_t i = order[begin + cursor++];
              if (puppi.pt()[i] <= 0.f)
                break;
              if (clusters.cluster()[i] < 0) {
                seedIndex = i;
                break;
              }
            }
            sumPt = sumEta = sumPhi = 0.f;
            daughters = 0;
          }
          alpaka::syncBlockThreads(acc);
          // Copy the shared seed before the membership barrier; only the
          // selecting thread reads the shared sums after that barrier.
          uint32_t seed = seedIndex;
          if (seed == end)
            break;
          float seedEta = puppi.eta()[seed], seedPhi = puppi.phi()[seed];
          for (uint32_t tid : independent_group_elements(acc, size)) {
            uint32_t i = begin + tid;
            if (clusters.cluster()[i] >= 0 || puppi.pt()[i] <= 0.f)
              continue;
            float deta = puppi.eta()[i] - seedEta;
            float dphi = cms::alpakatools::deltaPhi(acc, puppi.phi()[i], seedPhi);
            if (deta * deta + dphi * dphi < R2) {
              float pt = puppi.pt()[i];
              clusters.cluster()[i] = iter;
              clusters.is_seed()[i] = i == seed;
              alpaka::atomicAdd(acc, &sumPt, pt, alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sumEta, pt * deta, alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sumPhi, pt * dphi, alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &daughters, 1u, alpaka::hierarchy::Threads{});
            }
          }
          alpaka::syncBlockThreads(acc);
          if (once_per_block(acc)) {
            uint32_t slot = begin + iter;
            jets.pt()[slot] = sumPt;
            jets.eta()[slot] = seedEta + sumEta / sumPt;
            jets.phi()[slot] = cms::alpakatools::reducePhiRange(acc, seedPhi + sumPhi / sumPt);
            jets.cluster()[slot] = iter;
            jets.numberOfDaughters()[slot] = daughters;
          }
        }
        if (once_per_block(acc)) {
          jetBxIndex.bx().bx()[bx] = bxIndex.bx().bx()[bx];
          jetBxLookup.offset().offset()[bx + 1] = iter;
          alpaka::atomicAdd(acc, nJetsTotal, iter, alpaka::hierarchy::Blocks{});
        }
        alpaka::syncBlockThreads(acc);
      }
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif
