#ifndef L1TriggerScouting_Phase2_plugins_alpaka_LinkTree_h
#define L1TriggerScouting_Phase2_plugins_alpaka_LinkTree_h

#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"
#include <limits>

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {
  using namespace cms::alpakatools;
  // Link tree: nearest higher-pT parent, root chasing and scalar centroids.
  class LinkTreeKernel {
  public:
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView bxIndex,
                                  float RLink2,
                                  float ptMin,
                                  int32_t* parent,
                                  ClustersDeviceCollection::View clusters,
                                  ClusterObjDeviceCollection::View jets,
                                  BxLookupSoA::View jetBxLookup,
                                  BxLookupSoA::View jetBxIndex,
                                  unsigned int* nJetsTotal) const {
      if (once_per_grid(acc))
        jetBxLookup.offset().offset()[0] = 0;
      auto& nroots = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      uint32_t gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t bx : independent_groups(acc, gridDim)) {
        uint32_t begin = bxLookup.offset().offset()[bx], end = bxLookup.offset().offset()[bx + 1];
        uint32_t size = end - begin;
        if (once_per_block(acc))
          nroots = 0;
        // Scan original indices in order to resolve nearest-distance ties.
        for (uint32_t tid : independent_group_elements(acc, size)) {
          uint32_t i = begin + tid;
          float pt = puppi.pt()[i], eta = puppi.eta()[i], phi = puppi.phi()[i];
          int32_t bestParent = -1;
          float bestDr2 = RLink2;
          if (pt >= ptMin) {
            bestParent = i;

            for (uint32_t j = begin; j < end; ++j) {
              if (!(puppi.pt()[j] > pt || (puppi.pt()[j] == pt && j < i)))
                continue;
              float deta = puppi.eta()[j] - eta;
              float dphi = cms::alpakatools::deltaPhi(acc, puppi.phi()[j], phi);
              float dr2 = deta * deta + dphi * dphi;
              if (dr2 < bestDr2) {
                bestDr2 = dr2;
                bestParent = j;
              }
            }
          }
          parent[i] = bestParent;
          clusters.is_seed()[i] = bestParent == static_cast<int32_t>(i);
        }
        alpaka::syncBlockThreads(acc);
        // The cluster buffer temporarily stores original root indices.
        for (uint32_t tid : independent_group_elements(acc, size)) {
          uint32_t i = begin + tid;
          int32_t root = parent[i] < 0 ? -1 : static_cast<int32_t>(i);
          if (root >= 0)
            while (parent[root] != root)
              root = parent[root];
          clusters.cluster()[i] = root;
        }
        alpaka::syncBlockThreads(acc);
        // Chasing is complete: reuse each root's parent slot for its jet ID.
        for (uint32_t tid : independent_group_elements(acc, size)) {
          uint32_t i = begin + tid;
          if (!clusters.is_seed()[i])
            continue;
          uint32_t cid = alpaka::atomicAdd(acc, &nroots, 1u, alpaka::hierarchy::Threads{});
          parent[i] = cid;
          jets.eta()[begin + cid] = puppi.eta()[i];
          jets.phi()[begin + cid] = puppi.phi()[i];
        }
        alpaka::syncBlockThreads(acc);
        for (uint32_t tid : independent_group_elements(acc, size)) {
          uint32_t i = begin + tid;
          int32_t root = clusters.cluster()[i];
          clusters.cluster()[i] = root >= 0 ? parent[root] : -1;
        }
        alpaka::syncBlockThreads(acc);
        // One writer per jet, deterministic original-order scalar sums.
        for (uint32_t cid : independent_group_elements(acc, nroots)) {
          uint32_t slot = begin + cid;
          float axisEta = jets.eta()[slot], axisPhi = jets.phi()[slot];
          float sumPt = 0.f, sumEta = 0.f, sumPhi = 0.f;
          uint32_t daughters = 0;
          for (uint32_t i = begin; i < end; ++i) {
            if (clusters.cluster()[i] != static_cast<int32_t>(cid))
              continue;
            float pt = puppi.pt()[i];
            sumPt += pt;
            sumEta += pt * (puppi.eta()[i] - axisEta);
            sumPhi += pt * cms::alpakatools::deltaPhi(acc, puppi.phi()[i], axisPhi);
            ++daughters;
          }
          jets.pt()[slot] = sumPt;
          jets.eta()[slot] = sumPt > 0.f ? axisEta + sumEta / sumPt : axisEta;
          jets.phi()[slot] = sumPt > 0.f ? cms::alpakatools::reducePhiRange(acc, axisPhi + sumPhi / sumPt) : axisPhi;
          jets.cluster()[slot] = cid;
          jets.numberOfDaughters()[slot] = daughters;
        }
        if (once_per_block(acc)) {
          jetBxIndex.bx().bx()[bx] = bxIndex.bx().bx()[bx];
          jetBxLookup.offset().offset()[bx + 1] = nroots;
          alpaka::atomicAdd(acc, nJetsTotal, nroots, alpaka::hierarchy::Blocks{});
        }
        alpaka::syncBlockThreads(acc);
      }
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif
