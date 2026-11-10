#ifndef L1TriggerScouting_Phase2_plugins_alpaka_SCNMSW_h
#define L1TriggerScouting_Phase2_plugins_alpaka_SCNMSW_h

#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"
#include <limits>
#include <algorithm>
#include <type_traits>

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {
  using namespace cms::alpakatools;
  // Weighted NMS with scalar final-constituent centroids. Alpha=0 selects nearest axes.
  // One-round weighted seeded cone. The common path stores particles in shared
  // memory; the large-BX path uses only a sorted index into the original input.
  // Recombination keeps the reference eta order, without floating-point atomics.
  template <bool AlphaZero>
  class SCNMSWOneRoundKernel {
  public:
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView bxIndex,
                                  float RSeed2,
                                  float RCen2,
                                  float RClu2,
                                  float RScan,
                                  float alpha,
                                  uint32_t* scratch,
                                  ClustersDeviceCollection::View clusters,
                                  ClusterObjDeviceCollection::View jets,
                                  BxLookupSoA::View jetBxLookup,
                                  BxLookupSoA::View jetBxIndex,
                                  unsigned int* nJetsTotal,
                                  unsigned int* nClusteredTotal,
                                  uint32_t* particleBxOffsets) const {
      if (once_per_grid(acc)) {
        jetBxLookup.offset().offset()[0] = 0;
        particleBxOffsets[0] = 0;
      }
      constexpr uint32_t capacity = 128;
      auto& cachePt = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
      auto& cacheEta = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
      auto& cachePhi = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
      using LocalIndex = uint16_t;
      using LocalLabel = int16_t;
      auto& cacheOrig = alpaka::declareSharedVar<LocalIndex[capacity], __COUNTER__>(acc);
      auto& cacheLabel = alpaka::declareSharedVar<LocalLabel[capacity], __COUNTER__>(acc);
      auto& cacheAxisEta = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
      auto& cacheAxisPhi = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      auto gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t bx : independent_groups(acc, gridDim)) {
        uint32_t begin = bxLookup.offset().offset()[bx], size = bxLookup.offset().offset()[bx + 1] - begin;
        auto process = [&](auto cacheTag) {
          constexpr bool cached = decltype(cacheTag)::value;
          auto origAt = [&](uint32_t k) -> uint32_t {
            return cached ? static_cast<uint32_t>(begin + cacheOrig[k]) : scratch[begin + k];
          };
          auto ptAt = [&](uint32_t k) { return cached ? cachePt[k] : puppi.pt()[origAt(k)]; };
          auto etaAt = [&](uint32_t k) { return cached ? cacheEta[k] : puppi.eta()[origAt(k)]; };
          auto phiAt = [&](uint32_t k) { return cached ? cachePhi[k] : puppi.phi()[origAt(k)]; };
          auto labelAt = [&](uint32_t k) { return cached ? cacheLabel[k] : clusters.cluster()[origAt(k)]; };
          float* axisEta = cached ? cacheAxisEta : jets.eta().data() + begin;
          float* axisPhi = cached ? cacheAxisPhi : jets.phi().data() + begin;
          if (once_per_block(acc))
            nseeds = 0;
          if constexpr (cached) {
            for (uint32_t k : independent_group_elements(acc, size))
              cacheEta[k] = puppi.eta()[begin + k];
          }
          alpaka::syncBlockThreads(acc);
          // Stable exact eta ranking; no quantization or input range assumption.
          for (uint32_t k : independent_group_elements(acc, size)) {
            float eta = cached ? cacheEta[k] : puppi.eta()[begin + k];
            uint32_t rank = 0;
            for (uint32_t j = 0; j < size; ++j) {
              float other = cached ? cacheEta[j] : puppi.eta()[begin + j];
              rank += other < eta || (other == eta && j < k);
            }
            if constexpr (cached)
              cacheOrig[rank] = k;
            else
              scratch[begin + rank] = begin + k;
          }
          alpaka::syncBlockThreads(acc);
          for (uint32_t k : independent_group_elements(acc, size)) {
            uint32_t i = origAt(k);
            if constexpr (cached) {
              cachePt[k] = puppi.pt()[i];
              cacheEta[k] = puppi.eta()[i];
              cachePhi[k] = puppi.phi()[i];
              cacheLabel[k] = -1;
            }
            clusters.cluster()[i] = -1;
            clusters.is_seed()[i] = 0;
          }
          alpaka::syncBlockThreads(acc);
          for (uint32_t k : independent_group_elements(acc, size)) {
            float pt = ptAt(k), eta = etaAt(k), phi = phiAt(k);
            float sumPt = pt, sumEta = 0.f, sumPhi = 0.f;
            bool seed = true;
            for (uint32_t j = k; j > 0;) {
              --j;
              float deta = etaAt(j) - eta;
              if (deta < -RScan)
                break;
              float dphi = cms::alpakatools::deltaPhi(acc, phiAt(j), phi);
              float dr2 = deta * deta + dphi * dphi;
              float other = ptAt(j);
              if (dr2 < RSeed2 && (other > pt || (other == pt && origAt(j) < origAt(k)))) {
                seed = false;
                break;
              }
              if (dr2 < RCen2) {
                sumPt += other;
                sumEta += other * deta;
                sumPhi += other * dphi;
              }
            }
            if (seed) {
              for (uint32_t j = k + 1; j < size; ++j) {
                float deta = etaAt(j) - eta;
                if (deta > RScan)
                  break;
                float dphi = cms::alpakatools::deltaPhi(acc, phiAt(j), phi);
                float dr2 = deta * deta + dphi * dphi;
                float other = ptAt(j);
                if (dr2 < RSeed2 && (other > pt || (other == pt && origAt(j) < origAt(k)))) {
                  seed = false;
                  break;
                }
                if (dr2 < RCen2) {
                  sumPt += other;
                  sumEta += other * deta;
                  sumPhi += other * dphi;
                }
              }
            }
            if (!seed)
              continue;
            uint32_t cid = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{});
            axisEta[cid] = sumPt > 0.f ? eta + sumEta / sumPt : eta;
            axisPhi[cid] = sumPt > 0.f ? cms::alpakatools::reducePhiRange(acc, phi + sumPhi / sumPt) : phi;
            // Alpha=0 uses only distance: no strength calculation or storage.
            if constexpr (!AlphaZero) {
              float safePt = sumPt > 1.e-6f ? sumPt : 1.e-6f;
              jets.pt()[begin + cid] = alpha == 2.f ? safePt * safePt : alpaka::math::pow(acc, safePt, alpha);
            }
            // These fields hold protojet metadata until assignment finishes.
            jets.cluster()[begin + cid] = origAt(k);
            clusters.cluster()[origAt(k)] = cid;
            clusters.is_seed()[origAt(k)] = 1;
            if constexpr (cached)
              cacheLabel[k] = cid;
          }
          alpaka::syncBlockThreads(acc);
          for (uint32_t k : independent_group_elements(acc, size)) {
            int cid = labelAt(k);
            if (cid < 0) {
              float eta = etaAt(k), phi = phiAt(k);
              float best = std::numeric_limits<float>::infinity(), bestDr2 = best;
              uint32_t bestSeed = begin + size;
              for (uint32_t j = 0; j < nseeds; ++j) {
                float deta = eta - axisEta[j];
                if (deta * deta >= RClu2)
                  continue;
                // Preserve equal-distance candidates for the seed-index tie.
                if constexpr (AlphaZero) {
                  if (deta * deta > bestDr2)
                    continue;
                }
                float dphi = cms::alpakatools::deltaPhi(acc, phi, axisPhi[j]);
                float dr2 = deta * deta + dphi * dphi;
                if (dr2 >= RClu2)
                  continue;
                float score = dr2;
                if constexpr (!AlphaZero)
                  score /= jets.pt()[begin + j];
                uint32_t seed = jets.cluster()[begin + j];
                if (score < best || (score == best && (dr2 < bestDr2 || (dr2 == bestDr2 && seed < bestSeed)))) {
                  best = score;
                  bestDr2 = dr2;
                  bestSeed = seed;
                  cid = j;
                }
              }
              clusters.cluster()[origAt(k)] = cid;
              if constexpr (cached)
                cacheLabel[k] = cid;
            }
          }
          alpaka::syncBlockThreads(acc);
          // Ordered gathers avoid atomics and match reference kinematic bits.
          for (uint32_t j : independent_group_elements(acc, nseeds)) {
            float eta = axisEta[j], phi = axisPhi[j];
            float sumPt = 0.f, sumEta = 0.f, sumPhi = 0.f;
            uint32_t daughters = 0;
            for (uint32_t k = 0; k < size; ++k) {
              if (labelAt(k) != static_cast<int>(j))
                continue;
              ++daughters;
              float pt = ptAt(k);
              sumPt += pt;
              sumEta += pt * (etaAt(k) - eta);
              sumPhi += pt * cms::alpakatools::deltaPhi(acc, phiAt(k), phi);
            }
            uint32_t slot = begin + j;
            jets.pt()[slot] = sumPt;
            jets.eta()[slot] = sumPt > 0.f ? eta + sumEta / sumPt : eta;
            jets.phi()[slot] = sumPt > 0.f ? cms::alpakatools::reducePhiRange(acc, phi + sumPhi / sumPt) : phi;
            jets.cluster()[slot] = j;
            jets.numberOfDaughters()[slot] = daughters;
          }
          alpaka::syncBlockThreads(acc);
          if (once_per_block(acc)) {
            uint32_t clustered = 0;
            for (uint32_t j = 0; j < nseeds; ++j)
              clustered += jets.numberOfDaughters()[begin + j];
            jetBxLookup.offset().offset()[bx + 1] = nseeds;
            jetBxIndex.bx().bx()[bx] = bxIndex.bx().bx()[bx];
            particleBxOffsets[bx + 1] = clustered;
            alpaka::atomicAdd(acc, nJetsTotal, nseeds, alpaka::hierarchy::Blocks{});
            alpaka::atomicAdd(acc, nClusteredTotal, clustered, alpaka::hierarchy::Blocks{});
          }
          alpaka::syncBlockThreads(acc);
        };
        if (size <= capacity)
          process(std::true_type{});
        else
          process(std::false_type{});
      }
    }
  };

  // Cached repeated rounds; earlier-round jets stay frozen.
  class SCNMSWMultiRoundKernel {
  public:
    template <typename TAcc>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView bxIndex,
                                  float RSeed2,
                                  float RCen2,
                                  float RClu2,
                                  float RScan,
                                  float alpha,
                                  int maxRounds,
                                  uint32_t* indices,
                                  ClusterObjDeviceCollection::View work,
                                  ClustersDeviceCollection::View clusters,
                                  ClusterObjDeviceCollection::View jets,
                                  BxLookupSoA::View jetBxLookup,
                                  BxLookupSoA::View jetBxIndex,
                                  unsigned int* nJetsTotal,
                                  unsigned int* nClusteredTotal,
                                  uint32_t* particleBxOffsets) const {
      if (once_per_grid(acc)) {
        jetBxLookup.offset().offset()[0] = 0;
        particleBxOffsets[0] = 0;
      }
      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      auto& remaining = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      uint32_t gridDim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t bx : independent_groups(acc, gridDim)) {
        uint32_t begin = bxLookup.offset().offset()[bx], end = bxLookup.offset().offset()[bx + 1];
        uint32_t size = end - begin;
        // Small-BX shared cache; larger BXs retain the global-memory path.
        constexpr uint32_t capacity = 128;
        auto& cachePt = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
        auto& cacheEta = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
        auto& cachePhi = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
        auto& cacheOrig = alpaka::declareSharedVar<int32_t[capacity], __COUNTER__>(acc);
        auto& cacheIndex = alpaka::declareSharedVar<uint32_t[capacity], __COUNTER__>(acc);
        auto& cacheLabels = alpaka::declareSharedVar<int32_t[capacity], __COUNTER__>(acc);
        auto& cacheJetEta = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
        auto& cacheJetPhi = alpaka::declareSharedVar<float[capacity], __COUNTER__>(acc);
        auto process = [&](auto cacheTag) {
          constexpr bool cached = decltype(cacheTag)::value;
          float* ptBuf = cached ? cachePt : work.pt().data() + begin;
          float* etaBuf = cached ? cacheEta : work.eta().data() + begin;
          float* phiBuf = cached ? cachePhi : work.phi().data() + begin;
          int32_t* origBuf = cached ? cacheOrig : work.cluster().data() + begin;
          uint32_t* indexBuf = cached ? cacheIndex : indices + begin;
          uint32_t* ownBuf = indexBuf;
          float* jetEtaBuf = cached ? cacheJetEta : jets.eta().data() + begin;
          float* jetPhiBuf = cached ? cacheJetPhi : jets.phi().data() + begin;
          if (cached) {
            for (uint32_t tid : independent_group_elements(acc, size))
              cacheEta[tid] = puppi.eta()[begin + tid];
          }
          alpaka::syncBlockThreads(acc);
          // Exact stable eta rank: no quantization/range assumption, and no
          // multi-pass radix-sort barriers for the small per-BX input arrays.
          for (uint32_t tid : independent_group_elements(acc, size)) {
            uint32_t i = begin + tid, rank = 0;
            float eta = cached ? cacheEta[tid] : puppi.eta()[i];
            for (uint32_t j = begin; j < end; ++j) {
              float other = cached ? cacheEta[j - begin] : puppi.eta()[j];
              rank += (other < eta || (other == eta && j < i));
            }
            indexBuf[rank] = tid;
          }
          alpaka::syncBlockThreads(acc);
          for (uint32_t tid : independent_group_elements(acc, size)) {
            uint32_t i = begin + tid;
            uint32_t orig = begin + indexBuf[i - begin];
            ptBuf[i - begin] = puppi.pt()[orig];
            etaBuf[i - begin] = puppi.eta()[orig];
            phiBuf[i - begin] = puppi.phi()[orig];
            origBuf[i - begin] = orig;
            clusters.cluster()[orig] = -1;
            if (cached)
              cacheLabels[i - begin] = -1;
            clusters.is_seed()[orig] = 0;
          }
          alpaka::syncBlockThreads(acc);
          uint32_t total = 0;
          uint32_t limit = maxRounds < 0 ? size : std::min(size, uint32_t(maxRounds));
          for (uint32_t round = 0; round < limit && size != 0; ++round) {
            if (once_per_block(acc)) {
              nseeds = 0;
              remaining = 0;
            }
            alpaka::syncBlockThreads(acc);
            for (uint32_t tid : independent_group_elements(acc, size)) {
              uint32_t i = begin + tid;
              // Reuse the completed eta-sort index buffer for seed ownership.
              ownBuf[i - begin] = std::numeric_limits<uint32_t>::max();
              if ((cached ? cacheLabels[i - begin] : clusters.cluster()[origBuf[i - begin]]) >= 0)
                continue;
              float pt = ptBuf[i - begin], eta = etaBuf[i - begin], phi = phiBuf[i - begin];
              float sumPt = pt, sumEta = 0.f, sumPhi = 0.f;
              bool seed = true;
              // Fused scans support independent seed and centroid radii.
              // Summation order remains seed, downward eta, upward eta.
              for (uint32_t j = i; j > begin;) {
                --j;
                float deta = etaBuf[j - begin] - eta;
                if (deta < -RScan)
                  break;
                if ((cached ? cacheLabels[j - begin] : clusters.cluster()[origBuf[j - begin]]) >= 0)
                  continue;
                float dphi = cms::alpakatools::deltaPhi(acc, phiBuf[j - begin], phi);
                float dr2 = deta * deta + dphi * dphi;
                if (dr2 < RSeed2 &&
                    (ptBuf[j - begin] > pt || (ptBuf[j - begin] == pt && origBuf[j - begin] < origBuf[i - begin]))) {
                  seed = false;
                  break;
                }
                if (dr2 < RCen2) {
                  sumPt += ptBuf[j - begin];
                  sumEta += ptBuf[j - begin] * deta;
                  sumPhi += ptBuf[j - begin] * dphi;
                }
              }
              if (seed) {
                for (uint32_t j = i + 1; j < end; ++j) {
                  float deta = etaBuf[j - begin] - eta;
                  if (deta > RScan)
                    break;
                  if ((cached ? cacheLabels[j - begin] : clusters.cluster()[origBuf[j - begin]]) >= 0)
                    continue;
                  float dphi = cms::alpakatools::deltaPhi(acc, phiBuf[j - begin], phi);
                  float dr2 = deta * deta + dphi * dphi;
                  if (dr2 < RSeed2 &&
                      (ptBuf[j - begin] > pt || (ptBuf[j - begin] == pt && origBuf[j - begin] < origBuf[i - begin]))) {
                    seed = false;
                    break;
                  }
                  if (dr2 < RCen2) {
                    sumPt += ptBuf[j - begin];
                    sumEta += ptBuf[j - begin] * deta;
                    sumPhi += ptBuf[j - begin] * dphi;
                  }
                }
              }
              if (!seed)
                continue;
              uint32_t local = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{});
              uint32_t cid = total + local, slot = begin + cid;

              jetEtaBuf[slot - begin] = sumPt > 0.f ? eta + sumEta / sumPt : eta;
              jetPhiBuf[slot - begin] = sumPt > 0.f ? cms::alpakatools::reducePhiRange(acc, phi + sumPhi / sumPt) : phi;
              jets.cluster()[slot] = origBuf[i - begin];
              jets.numberOfDaughters()[slot] = 0;
              float safePt = sumPt > 1.e-6f ? sumPt : 1.e-6f;
              float strength = alpha == 2.f ? safePt * safePt : alpaka::math::pow(acc, safePt, alpha);
              jets.pt()[slot] = strength;
              ownBuf[i - begin] = cid;
              clusters.is_seed()[origBuf[i - begin]] = 1;
            }
            alpaka::syncBlockThreads(acc);
            if (nseeds == 0)
              break;
            for (uint32_t tid : independent_group_elements(acc, size)) {
              uint32_t i = begin + tid, orig = origBuf[i - begin];
              if ((cached ? cacheLabels[i - begin] : clusters.cluster()[orig]) >= 0)
                continue;
              uint32_t own = ownBuf[i - begin];
              int cid = own == std::numeric_limits<uint32_t>::max() ? -1 : static_cast<int>(own);
              if (cid < 0) {
                float best = std::numeric_limits<float>::infinity();
                float bestDr2 = std::numeric_limits<float>::infinity();
                uint32_t bestSeed = end;
                for (uint32_t j = 0; j < nseeds; ++j) {
                  uint32_t slot = begin + total + j;
                  float deta = etaBuf[i - begin] - jetEtaBuf[slot - begin];
                  if (deta * deta >= RClu2)
                    continue;
                  float dphi = cms::alpakatools::deltaPhi(acc, phiBuf[i - begin], jetPhiBuf[slot - begin]);
                  float dr2 = deta * deta + dphi * dphi;
                  if (dr2 >= RClu2)
                    continue;
                  float score = dr2 / jets.pt()[slot];
                  uint32_t seedOrig = static_cast<uint32_t>(jets.cluster()[slot]);
                  if (score < best || (score == best && (dr2 < bestDr2 || (dr2 == bestDr2 && seedOrig < bestSeed)))) {
                    best = score;
                    bestDr2 = dr2;
                    bestSeed = seedOrig;
                    cid = total + j;
                  }
                }
              }
              clusters.cluster()[orig] = cid;
              if (cached)
                cacheLabels[i - begin] = cid;
              if (cid < 0)
                alpaka::atomicAdd(acc, &remaining, 1u, alpaka::hierarchy::Threads{});
              // Own-seed metadata is dead now; reuse it for sorted assignments.
              ownBuf[i - begin] = static_cast<uint32_t>(cid);
            }
            alpaka::syncBlockThreads(acc);
            for (uint32_t j : independent_group_elements(acc, nseeds)) {
              uint32_t cid = total + j, slot = begin + cid;
              jets.cluster()[slot] = cid;
              float axisEta = jetEtaBuf[slot - begin], axisPhi = jetPhiBuf[slot - begin];
              float sum0 = 0.f, sum1 = 0.f, sum2 = 0.f;
              uint32_t daughters = 0;
              for (uint32_t i = begin; i < end; ++i) {
                if (ownBuf[i - begin] != cid)
                  continue;
                ++daughters;
                float pt = ptBuf[i - begin];
                sum0 += pt;
                sum1 += pt * (etaBuf[i - begin] - axisEta);
                sum2 += pt * cms::alpakatools::deltaPhi(acc, phiBuf[i - begin], axisPhi);
              }
              jets.numberOfDaughters()[slot] = daughters;
              jets.pt()[slot] = sum0;
              if (sum0 > 0.f) {
                jetEtaBuf[slot - begin] = axisEta + sum1 / sum0;
                jetPhiBuf[slot - begin] = cms::alpakatools::reducePhiRange(acc, axisPhi + sum2 / sum0);
              }
              if (cached) {
                jets.eta()[slot] = jetEtaBuf[slot - begin];
                jets.phi()[slot] = jetPhiBuf[slot - begin];
              }
            }
            alpaka::syncBlockThreads(acc);
            uint32_t roundSeeds = nseeds;
            uint32_t left = remaining;
            // All warps must finish reading the shared counters before
            // the next round resets them.
            alpaka::syncBlockThreads(acc);
            total += roundSeeds;
            if (left == 0)
              break;
          }
          if (once_per_block(acc)) {
            jetBxIndex.bx().bx()[bx] = bxIndex.bx().bx()[bx];
            jetBxLookup.offset().offset()[bx + 1] = total;
            alpaka::atomicAdd(acc, nJetsTotal, total, alpaka::hierarchy::Blocks{});
            uint32_t clustered = 0;
            for (uint32_t j = 0; j < total; ++j)
              clustered += jets.numberOfDaughters()[begin + j];
            particleBxOffsets[bx + 1] = clustered;
            alpaka::atomicAdd(acc, nClusteredTotal, clustered, alpaka::hierarchy::Blocks{});
          }
          alpaka::syncBlockThreads(acc);
        };
        if (size <= capacity)
          process(std::true_type{});
        else
          process(std::false_type{});
      }
    }
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif
