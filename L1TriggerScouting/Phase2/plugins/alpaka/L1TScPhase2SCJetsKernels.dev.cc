#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
// #include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
// #include "L1TriggerScouting/Phase2/plugins/alpaka/radixSort128.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/radixSort64.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"

//#define L1TSC_VERBOSE_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  constexpr uint32_t kThreadsPerBlock = 64;

  // ----------------------------------------------------------------------
  // SCNMS = old non-iterative seeded cone with split radii
  //
  // This is intentionally the old kernel logic, with only one conceptual
  // change:
  //
  //   - RSeed is used in the old "pre-cluster" step:
  //       * local-max seed finding
  //       * accumulation of the old pT-weighted centroid around the seed
  //
  //   - RClu is used only in the final assignment step:
  //       * each particle is assigned to the nearest accepted seed axis
  //         if that axis is within RClu
  //
  // So compared to the old single-radius kernel:
  //   old R  -> split into RSeed for step-2 and RClu for step-3
  //
  // Everything else is kept as close as possible to the old behavior.
  // ----------------------------------------------------------------------
  class JetKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi, // original input p.
                                  OffsetsSoA::ConstView bxLookup, // p. ranges per BX
                                  BxIndexSoA::ConstView bxIndex, // BX label per BX
                                  float RSeed,
                                  float RSeed2,
                                  float RClu2,
                                  uint16_t* uieta, // tmp eta sort key
                                  uint16_t* idx, // tmp sort index buffer
                                  ClusterObjDeviceCollection::View work, // eta sorted working p. cp
                                  ClustersDeviceCollection::View clusters, // output p. assignment info
                                  ClusterObjDeviceCollection::View jets, // tmp jet collection
                                  OffsetsSoA::View jetBxLookup, // number jets per BX
                                  BxIndexSoA::View jetBxIndex, // BX labels for jet collection
                                  unsigned int* nJetsTotal) const {
      // check whether this backend runs 1 thread per block (for optimization with seed veto markers later)
      constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;

      // init first lookup entry to 0
      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      // get number of BX blocks
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      // step-0: sort by eta in blocks
      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get particle range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;

        // for each p. in BX build quantized eta value (map form [-5,5] to [0,65535]) in uieta & store its BX index in idx
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uieta[tid + begin] = (puppi.eta()[tid + begin] + 5.f) * (std::numeric_limits<uint16_t>::max() / 10.0f);
          idx[tid + begin] = tid;  // important for GPU radixSort implementation
        }
      }

      // sort p. (uieta,idx) by eta inside each BX
      radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

      // step-1: rearrange
      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get p. range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;

        // copy eta sorted p. into work
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          auto isrc = idx[ipart] + begin;
          work.pt()[ipart] = puppi.pt()[isrc];
          work.eta()[ipart] = puppi.eta()[isrc];
          work.phi()[ipart] = puppi.phi()[isrc];
          work.cluster()[ipart] = isrc;
        }
      }

      // per BX shared counter nseeds
      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      if (once_per_block(acc))
        nseeds = 0;
      alpaka::syncBlockThreads(acc);

      // step-2: seed finding + old centroid accumulation using RSeed
      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get p. range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;

        // pre-cluster
        // loop over every p. as seed candidate
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          // define cand seed; init jet sums
          uint32_t iseed = tid + begin, icluster = work.cluster()[iseed];
          float seed_pt = work.pt()[iseed], seed_eta = work.eta()[iseed], seed_phi = work.phi()[iseed];
          float sum_pt = seed_pt, sum_eta = 0, sum_phi = 0;
          bool is_seed = true;

          // on single-thread backends we can skip already vetoed seeds
          if constexpr (single_thread) {
            if (clusters.is_seed()[icluster] == -1)
              continue;
          }

          // scan downward in eta
          for (uint32_t j = tid; j > 0; --j) {
            uint32_t ipart = j - 1 + begin;
            float deta = work.eta()[ipart] - seed_eta;
            if (deta < -RSeed)
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < RSeed2) {
              if (work.pt()[ipart] >= seed_pt) {
                is_seed = false;
                break;
              } else {
                if constexpr (single_thread) {
                  clusters.is_seed()[work.cluster()[ipart]] = -1;
                }
                sum_pt += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }

          // scan upward in eta
          for (uint32_t j = tid + 1; j < block_dim; ++j) {
            uint32_t ipart = j + begin;
            float deta = work.eta()[ipart] - seed_eta;
            if (deta > RSeed)
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < RSeed2) {
              if (work.pt()[ipart] > seed_pt) {
                is_seed = false;
                break;
              } else {
                if constexpr (single_thread) {
                  clusters.is_seed()[work.cluster()[ipart]] = -1;
                }
                sum_pt += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }

          // convert centroid sums (offsets) to jet axis
          // important: this is still the old centroid logic, so it is built
          // from the same RSeed neighborhood used during seed competition
          sum_eta = seed_eta + sum_eta / sum_pt;
          sum_phi = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);

          // seed cand survived
          if (is_seed) {
            auto ijet = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{}) + begin;
            jets.pt()[ijet] = sum_pt;
            jets.eta()[ijet] = sum_eta;
            jets.phi()[ijet] = sum_phi;
            jets.cluster()[ijet] = ijet - begin;
            jets.numberOfDaughters()[ijet] = 0;
            clusters.is_seed()[icluster] = 1;
          }
        }

        alpaka::syncBlockThreads(acc);

        // step 3: re-associate p. to nearest accepted seed using RClu
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          uint32_t icluster = work.cluster()[ipart];

          // keep nearest seed only if inside RClu
          float nearest = RClu2;
          int jcluster = -1;

          for (uint32_t j = 0; j < nseeds; ++j) {
            auto jseed = j + begin;
            float deta = work.eta()[ipart] - jets.eta()[jseed];
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], jets.phi()[jseed]);
            float dr2 = deta * deta + dphi * dphi;
            if (dr2 < nearest) {
              jcluster = jets.cluster()[jseed];
              nearest = dr2;
            }
          }

          if (jcluster != -1) {
            alpaka::atomicAdd(acc, &jets.numberOfDaughters()[jcluster + begin], 1u, alpaka::hierarchy::Threads{});
          }

          clusters.cluster()[icluster] = jcluster;

          // clear veto markers (single-thread backends) -> clusters.is_seed() 1 for seed, 0 else
          if constexpr (single_thread) {
            clusters.is_seed()[icluster] = std::max(clusters.is_seed()[icluster], 0);
          }
        }

        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx];
          jetBxLookup.offsets()[block_idx + 1] = nseeds;
          alpaka::atomicAdd(acc, nJetsTotal, nseeds, alpaka::hierarchy::Blocks{});
        }

        alpaka::syncBlockThreads(acc);
      }
    }
  };

  // ----------------------------------------------------------------------
  // SCNMSWeighted = old weighted non-iterative seeded cone with split radii
  //
  // Same idea as JetKernel above:
  //   - RSeed is used for seed competition + old centroid accumulation
  //   - RClu is used only for the final assignment step
  //
  // The only difference with unweighted SCNMS is the assignment metric:
  //
  //   metric = dr2 / pt_jet^2
  //
  // This preserves the old weighted behavior while allowing RSeed and RClu
  // to differ.
  // ----------------------------------------------------------------------
  class JetKernelWeighted {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi, // original input p.
                                  OffsetsSoA::ConstView bxLookup, // p. ranges per BX
                                  BxIndexSoA::ConstView bxIndex, // BX label per BX
                                  float RSeed,
                                  float RSeed2,
                                  float RClu2,
                                  uint16_t* uieta, // tmp eta sort key
                                  uint16_t* idx, // tmp sort index buffer
                                  ClusterObjDeviceCollection::View work, // eta sorted working p. cp
                                  ClustersDeviceCollection::View clusters, // output p. assignment info
                                  ClusterObjDeviceCollection::View jets, // tmp jet collection
                                  OffsetsSoA::View jetBxLookup, // number jets per BX
                                  BxIndexSoA::View jetBxIndex, // BX labels for jet collection
                                  unsigned int* nJetsTotal) const {
      // check whether this backend runs 1 thread per block (for optimization with seed veto markers later)
      constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;

      // init first lookup entry to 0
      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      // get number of BX blocks
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      // step-0: sort by eta in blocks
      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get particle range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;

        // for each p. in BX build quantized eta value (map form [-5,5] to [0,65535]) in uieta & store its BX index in idx
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uieta[tid + begin] = (puppi.eta()[tid + begin] + 5.f) * (std::numeric_limits<uint16_t>::max() / 10.0f);
          idx[tid + begin] = tid;  // this is important for the GPU implementation of radixSort
        }
      }

      // sort p. (uieta,idx) by eta inside each BX
      radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

      // step-1: rearrange
      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get p. range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;

        // copy eta sorted p. into work
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          auto isrc = idx[ipart] + begin;
          work.pt()[ipart] = puppi.pt()[isrc];
          work.eta()[ipart] = puppi.eta()[isrc];
          work.phi()[ipart] = puppi.phi()[isrc];
          work.cluster()[ipart] = isrc;
        }
      }

      // per BX shared counter nseeds
      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      if (once_per_block(acc))
        nseeds = 0;
      alpaka::syncBlockThreads(acc);

      // step-2: seed finding + old centroid accumulation using RSeed
      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get p. range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;

        // pre-cluster
        // loop over every p. as seed candidate
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          // define cand seed; init jet sums
          uint32_t iseed = tid + begin, icluster = work.cluster()[iseed];
          float seed_pt = work.pt()[iseed], seed_eta = work.eta()[iseed], seed_phi = work.phi()[iseed];
          float sum_pt = seed_pt, sum_eta = 0, sum_phi = 0;
          bool is_seed = true;

          // on single-thread backends we can skip already vetoed seeds
          if constexpr (single_thread) {
            if (clusters.is_seed()[icluster] == -1)
              continue;
          }

          // scan downward in eta
          for (uint32_t j = tid; j > 0; --j) {
            uint32_t ipart = j - 1 + begin;
            float deta = work.eta()[ipart] - seed_eta;
            if (deta < -RSeed)
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < RSeed2) {
              if (work.pt()[ipart] >= seed_pt) {
                is_seed = false;
                break;
              } else {
                if constexpr (single_thread) {
                  clusters.is_seed()[work.cluster()[ipart]] = -1;
                }
                sum_pt += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }

          // scan upward in eta
          for (uint32_t j = tid + 1; j < block_dim; ++j) {
            uint32_t ipart = j + begin;
            float deta = work.eta()[ipart] - seed_eta;
            if (deta > RSeed)
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < RSeed2) {
              if (work.pt()[ipart] > seed_pt) {
                is_seed = false;
                break;
              } else {
                if constexpr (single_thread) {
                  clusters.is_seed()[work.cluster()[ipart]] = -1;
                }
                sum_pt += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }

          // convert centroid sums (offsets) to jet axis
          // important: this is still the old centroid logic, so it is built
          // from the same RSeed neighborhood used during seed competition
          sum_eta = seed_eta + sum_eta / sum_pt;
          sum_phi = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);

          // seed cand survived
          if (is_seed) {
            auto ijet = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{}) + begin;
            jets.pt()[ijet] = sum_pt;
            jets.eta()[ijet] = sum_eta;
            jets.phi()[ijet] = sum_phi;
            jets.cluster()[ijet] = ijet - begin;
            jets.numberOfDaughters()[ijet] = 0;
            clusters.is_seed()[icluster] = 1;
          }
        }

        // sync after seed finding
        alpaka::syncBlockThreads(acc);

        // step 3: reassociate p. to best accepted seed using RClu
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          uint32_t icluster = work.cluster()[ipart];

          // weighted mode:
          // keep the best weighted metric among seeds inside RClu
          float nearest = 1e30f;
          int jcluster = -1;

          for (uint32_t j = 0; j < nseeds; ++j) {
            auto jseed = j + begin;
            float deta = work.eta()[ipart] - jets.eta()[jseed];
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], jets.phi()[jseed]);
            float dr2 = deta * deta + dphi * dphi;

            if (dr2 < RClu2) {
              float jet_pt = jets.pt()[jseed];
              float metric = dr2 / (jet_pt * jet_pt + 1e-12f);
              if (metric < nearest) {
                jcluster = jets.cluster()[jseed];
                nearest = metric;
              }
            }
          }

          if (jcluster != -1) {
            alpaka::atomicAdd(acc, &jets.numberOfDaughters()[jcluster + begin], 1u, alpaka::hierarchy::Threads{});
          }

          clusters.cluster()[icluster] = jcluster;

          // clear veto markers (single-thread backends) -> clusters.is_seed() 1 for seed, 0 else
          if constexpr (single_thread) {
            clusters.is_seed()[icluster] = std::max(clusters.is_seed()[icluster], 0);
          }
        }

        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx];
          jetBxLookup.offsets()[block_idx + 1] = nseeds;
          alpaka::atomicAdd(acc, nJetsTotal, nseeds, alpaka::hierarchy::Blocks{});
        }

        alpaka::syncBlockThreads(acc);
      }
    }
  };

  class JetZSKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  OffsetsSoA::ConstView bxLookup,
                                  ClusterObjDeviceCollection::ConstView jetsNonZS,
                                  OffsetsSoA::ConstView jetBxLookup,
                                  ClusterObjDeviceCollection::View jets,
                                  unsigned int* nClustered) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t beginSrc = bxLookup.offsets()[block_idx];
        uint32_t beginDst = jetBxLookup.offsets()[block_idx];
        uint32_t endDst = jetBxLookup.offsets()[block_idx + 1];
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

  // ----------------------------------------------------------------------
  // generic NMS-style seeded cone kernel
  //
  // This one is used for:
  //   SCNMSWeightedMultiIter
  //
  // Main logic:
  //   Step 0: sort PF by eta
  //   Step 1: rearrange into eta-sorted work array
  //   Step 2: local-max seed finding inside RSeed
  //   Step 3: optional centroid refinement inside RCen
  //   Step 4: assignment to best proto-jet inside RClu
  //
  // unweighted mode:   score = dr2
  // weighted mode:     score = dr2 / pt_seed^alpha
  // ----------------------------------------------------------------------
  class JetKernelNMS {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        PuppiDeviceCollection::ConstView puppi,
        OffsetsSoA::ConstView bxLookup,
        BxIndexSoA::ConstView bxIndex,
        float RSeed,
        float RSeed2,
        float RCen,
        float RCen2,
        float RClu,
        float RClu2,
        float alphaSeed,
        float minSeedPt,
        unsigned int nCentroidIters,
        bool useWeightedAssignment,
        uint16_t* uieta,
        uint16_t* idx,
        float* seedWeight,
        uint16_t* seedPos,
        ClusterObjDeviceCollection::View work,
        ClustersDeviceCollection::View clusters,
        ClusterObjDeviceCollection::View jetsNonZS,
        OffsetsSoA::View jetBxLookup,
        BxIndexSoA::View jetBxIndex,
        unsigned int* nJetsTotal) const {

      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      // ============================================================
      // Step 0: eta sort
      // ============================================================
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uint32_t i = tid + begin;
          uieta[i] = (puppi.eta()[i] + 5.f) * (std::numeric_limits<uint16_t>::max() / 10.0f);
          idx[i] = tid;
        }
      }

      radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

      // ============================================================
      // Step 1: rearrange into eta-sorted work array
      // ============================================================
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uint32_t ip   = tid + begin;
          uint32_t isrc = idx[ip] + begin;

          work.pt()[ip]      = puppi.pt()[isrc];
          work.eta()[ip]     = puppi.eta()[isrc];
          work.phi()[ip]     = puppi.phi()[isrc];
          work.cluster()[ip] = isrc; // original particle index
        }
      }

      // per-BX shared counter
      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

      // ============================================================
      // BX loop
      // ============================================================
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {

        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        if (once_per_block(acc))
          nseeds = 0u;

        alpaka::syncBlockThreads(acc);

        // ============================================================
        // Step 2: seed finding (local pT maxima within RSeed)
        // ============================================================
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {

          uint32_t iseed = tid + begin;

          float seed_pt  = work.pt()[iseed];
          if (seed_pt < minSeedPt)
            continue;

          float seed_eta = work.eta()[iseed];
          float seed_phi = work.phi()[iseed];
          uint32_t seed_orig = static_cast<uint32_t>(work.cluster()[iseed]);

          bool is_seed = true;

          // scan down
          for (uint32_t j = tid; j > 0; --j) {
            uint32_t ip = begin + (j - 1);

            float deta = work.eta()[ip] - seed_eta;
            if (deta < -RSeed) break;

            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ip], seed_phi);
            float dr2 = deta*deta + dphi*dphi;

            if (dr2 < RSeed2) {
              float ptj = work.pt()[ip];
              uint32_t origj = static_cast<uint32_t>(work.cluster()[ip]);

              // higher-pt wins; for exact ties choose lower original index
              if (ptj > seed_pt || (ptj == seed_pt && origj < seed_orig)) {
                is_seed = false;
                break;
              }
            }
          }

          // scan up
          if (is_seed) {
            for (uint32_t j = tid + 1; j < block_dim; ++j) {
              uint32_t ip = begin + j;

              float deta = work.eta()[ip] - seed_eta;
              if (deta > RSeed) break;

              float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ip], seed_phi);
              float dr2 = deta*deta + dphi*dphi;

              if (dr2 < RSeed2) {
                float ptj = work.pt()[ip];
                uint32_t origj = static_cast<uint32_t>(work.cluster()[ip]);

                if (ptj > seed_pt || (ptj == seed_pt && origj < seed_orig)) {
                  is_seed = false;
                  break;
                }
              }
            }
          }

          if (!is_seed)
            continue;

          // accepted seed -> create proto-jet using seed axis as initial axis
          uint32_t jlocal = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{});
          uint32_t jseed  = begin + jlocal;

          jetsNonZS.pt()[jseed]  = seed_pt;
          jetsNonZS.eta()[jseed] = seed_eta;
          jetsNonZS.phi()[jseed] = seed_phi;
          jetsNonZS.cluster()[jseed] = jlocal; // local jet id inside BX
          jetsNonZS.numberOfDaughters()[jseed] = 0u;

          seedPos[jseed] = static_cast<uint16_t>(tid);

          // in unweighted mode this value is ignored
          {
            float denom = alpaka::math::pow(
                acc,
                alpaka::math::max(acc, seed_pt, 1e-6f),
                alphaSeed);
            seedWeight[jseed] = 1.f / denom;
          }

          clusters.is_seed()[seed_orig] = 1;
        }

        alpaka::syncBlockThreads(acc);

        // ============================================================
        // Step 3: optional centroid iterations
        // nCentroidIters = 0 -> keep seed axis
        // ============================================================
        if (nCentroidIters > 0u) {
          for (uint32_t j : independent_group_elements(acc, nseeds)) {

            uint32_t jseed = begin + j;

            float axis_eta = jetsNonZS.eta()[jseed];
            float axis_phi = jetsNonZS.phi()[jseed];
            uint32_t pos   = static_cast<uint32_t>(seedPos[jseed]);

            for (unsigned int iter = 0; iter < nCentroidIters; ++iter) {

              float sum_pt  = 0.f;
              float sum_eta = 0.f;
              float sum_sin = 0.f;
              float sum_cos = 0.f;

              // scan down from seed position
              for (int32_t k = static_cast<int32_t>(pos); k >= 0; --k) {
                uint32_t ip = begin + static_cast<uint32_t>(k);

                float deta = work.eta()[ip] - axis_eta;
                if (deta < -RCen) break;

                float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ip], axis_phi);
                float dr2 = deta*deta + dphi*dphi;

                if (dr2 < RCen2) {
                  float w = work.pt()[ip];
                  sum_pt  += w;
                  sum_eta += w * work.eta()[ip];

                  float ph = work.phi()[ip];
                  sum_sin += w * alpaka::math::sin(acc, ph);
                  sum_cos += w * alpaka::math::cos(acc, ph);
                }
              }

              // scan up from seed position
              for (uint32_t k = pos + 1; k < block_dim; ++k) {
                uint32_t ip = begin + k;

                float deta = work.eta()[ip] - axis_eta;
                if (deta > RCen) break;

                float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ip], axis_phi);
                float dr2 = deta*deta + dphi*dphi;

                if (dr2 < RCen2) {
                  float w = work.pt()[ip];
                  sum_pt  += w;
                  sum_eta += w * work.eta()[ip];

                  float ph = work.phi()[ip];
                  sum_sin += w * alpaka::math::sin(acc, ph);
                  sum_cos += w * alpaka::math::cos(acc, ph);
                }
              }

              if (sum_pt > 0.f) {
                axis_eta = sum_eta / sum_pt;
                axis_phi = alpaka::math::atan2(acc, sum_sin, sum_cos);

                // store refined proto-jet axis
                // for weighted assignment we also use proto-jet pt from the centroid catchment
                jetsNonZS.pt()[jseed]  = sum_pt;
                jetsNonZS.eta()[jseed] = axis_eta;
                jetsNonZS.phi()[jseed] = axis_phi;
              }
            }
          }
        }

        alpaka::syncBlockThreads(acc);

        // ============================================================
        // Step 4: assignment to best proto-jet
        // SCNMSWeightedMultiIter -> weighted distance = dr2 / pt_seed^alpha
        // ============================================================
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {

          uint32_t ip   = begin + tid;
          uint32_t orig = static_cast<uint32_t>(work.cluster()[ip]);

          float eta_i = work.eta()[ip];
          float phi_i = work.phi()[ip];

          float bestScore = std::numeric_limits<float>::infinity();
          int bestJ = -1;

          for (uint32_t j = 0; j < nseeds; ++j) {

            uint32_t jseed = begin + j;

            float deta = jetsNonZS.eta()[jseed] - eta_i;
            if (deta < -RClu || deta > RClu)
              continue;

            float dphi = cms::alpakatools::deltaPhi(acc, jetsNonZS.phi()[jseed], phi_i);
            float dr2 = deta*deta + dphi*dphi;

            if (dr2 >= RClu2)
              continue;

            float score = dr2;
            if (useWeightedAssignment) {
              score = dr2 * seedWeight[jseed]; // = dr2 / pt_seed^alpha
            }

            if (score < bestScore) {
              bestScore = score;
              bestJ = static_cast<int>(j);
            }
          }

          clusters.cluster()[orig] = bestJ;

          // non-seed particles: explicitly set is_seed = 0
          // seed particles remain 1 from seed-finding step
          if (clusters.is_seed()[orig] != 1)
            clusters.is_seed()[orig] = 0;

          if (bestJ >= 0) {
            alpaka::atomicAdd(
              acc,
              &jetsNonZS.numberOfDaughters()[begin + bestJ],
              1u,
              alpaka::hierarchy::Threads{});
          }
        }

        alpaka::syncBlockThreads(acc);

        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx];
          jetBxLookup.offsets()[block_idx + 1] = nseeds;
          alpaka::atomicAdd(acc, nJetsTotal, nseeds, alpaka::hierarchy::Blocks{});
        }

        alpaka::syncBlockThreads(acc);
      }
    }
  };

  class JetKernelLinkTree {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        PuppiDeviceCollection::ConstView puppi,
        OffsetsSoA::ConstView bxLookup,
        BxIndexSoA::ConstView bxIndex,
        float RLink,
        float RLink2,
        float ptMin,
        uint16_t* uieta,
        uint16_t* idx,
        int32_t* parent,
        int32_t* root,
        int32_t* rootCid,
        float* sumPx,
        float* sumPy,
        float* sumPz,
        ClusterObjDeviceCollection::View work,
        ClustersDeviceCollection::View clusters,
        ClusterObjDeviceCollection::View jetsNonZS,
        OffsetsSoA::View jetBxLookup,
        BxIndexSoA::View jetBxIndex,
        unsigned int* nJetsTotal) const {

      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      // Step 0: sort PF by eta (per BX)
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uint32_t i = tid + begin;
          uieta[i] = (puppi.eta()[i] + 5.f) * (std::numeric_limits<uint16_t>::max() / 10.0f);
          idx[i] = tid;
        }
      }

      radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

      // Step 1: rearrange into eta-sorted work array
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uint32_t ip   = tid + begin;
          uint32_t isrc = idx[ip] + begin;

          work.pt()[ip]      = puppi.pt()[isrc];
          work.eta()[ip]     = puppi.eta()[isrc];
          work.phi()[ip]     = puppi.phi()[isrc];
          work.cluster()[ip] = isrc;
        }
      }

      auto& nroots = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {

        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        if (once_per_block(acc))
          nroots = 0u;

        alpaka::syncBlockThreads(acc);

        // Step 2: nearest-higher parent within RLink
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {

          uint32_t ip = tid + begin;

          float pti = work.pt()[ip];
          if (pti < ptMin) {
            parent[ip] = -1;
            continue;
          }

          float etai = work.eta()[ip];
          float phii = work.phi()[ip];
          uint32_t origi = static_cast<uint32_t>(work.cluster()[ip]);

          float bestDr2 = std::numeric_limits<float>::infinity();
          int32_t bestParent = -1;

          // scan down
          for (uint32_t j = tid; j > 0; --j) {
            uint32_t jp = begin + (j - 1);

            float deta = work.eta()[jp] - etai;
            if (deta < -RLink)
              break;

            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[jp], phii);
            float dr2  = deta * deta + dphi * dphi;

            if (dr2 >= RLink2)
              continue;

            float ptj = work.pt()[jp];
            uint32_t origj = static_cast<uint32_t>(work.cluster()[jp]);

            bool higher = (ptj > pti) || ((ptj == pti) && (origj < origi));
            if (!higher)
              continue;

            if (dr2 < bestDr2) {
              bestDr2 = dr2;
              bestParent = static_cast<int32_t>(jp);
            }
          }

          // scan up
          for (uint32_t j = tid + 1; j < block_dim; ++j) {
            uint32_t jp = begin + j;

            float deta = work.eta()[jp] - etai;
            if (deta > RLink)
              break;

            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[jp], phii);
            float dr2  = deta * deta + dphi * dphi;

            if (dr2 >= RLink2)
              continue;

            float ptj = work.pt()[jp];
            uint32_t origj = static_cast<uint32_t>(work.cluster()[jp]);

            bool higher = (ptj > pti) || ((ptj == pti) && (origj < origi));
            if (!higher)
              continue;

            if (dr2 < bestDr2) {
              bestDr2 = dr2;
              bestParent = static_cast<int32_t>(jp);
            }
          }

          if (bestParent < 0)
            parent[ip] = static_cast<int32_t>(ip);
          else
            parent[ip] = bestParent;
        }

        alpaka::syncBlockThreads(acc);

        // Step 3: root chasing
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {

          uint32_t ip = tid + begin;

          if (parent[ip] < 0) {
            root[ip] = -1;
            continue;
          }

          int32_t j = static_cast<int32_t>(ip);

          while (true) {
            int32_t pj = parent[j];

            if (pj < 0) {
              root[ip] = -1;
              break;
            }

            if (pj == j) {
              root[ip] = j;
              break;
            }

            j = pj;
          }
        }

        alpaka::syncBlockThreads(acc);

        // Step 4: compact roots to cluster ids
        if (once_per_block(acc)) {

          nroots = 0u;

          for (uint32_t ip = begin; ip < end; ++ip) {
            if (parent[ip] == static_cast<int32_t>(ip)) {
              rootCid[ip] = static_cast<int32_t>(nroots);
              ++nroots;
            }
          }
        }

        alpaka::syncBlockThreads(acc);

        // Step 5: initialize jet sums
        for (uint32_t tid : independent_group_elements(acc, nroots)) {

          uint32_t slot = begin + tid;

          sumPx[slot] = 0.f;
          sumPy[slot] = 0.f;
          sumPz[slot] = 0.f;

          jetsNonZS.pt()[slot] = 0.f;
          jetsNonZS.eta()[slot] = 0.f;
          jetsNonZS.phi()[slot] = 0.f;
          jetsNonZS.cluster()[slot] = tid;
          jetsNonZS.numberOfDaughters()[slot] = 0u;
        }

        alpaka::syncBlockThreads(acc);

        // Step 6: accumulate jet momentum
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {

          uint32_t ip   = tid + begin;
          uint32_t orig = static_cast<uint32_t>(work.cluster()[ip]);

          if (parent[ip] < 0 || root[ip] < 0) {
            clusters.cluster()[orig] = -1;
            clusters.is_seed()[orig] = 0;
            continue;
          }

          int32_t r   = root[ip];
          int32_t cid = rootCid[r];

          clusters.cluster()[orig] = cid;
          clusters.is_seed()[orig] = (parent[ip] == static_cast<int32_t>(ip)) ? 1 : 0;

          if (cid >= 0) {

            uint32_t slot = begin + static_cast<uint32_t>(cid);

            float pt  = work.pt()[ip];
            float phi = work.phi()[ip];
            float eta = work.eta()[ip];

            float px = pt * alpaka::math::cos(acc, phi);
            float py = pt * alpaka::math::sin(acc, phi);
            float pz = pt * alpaka::math::sinh(acc, eta);

            alpaka::atomicAdd(acc, &sumPx[slot], px, alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &sumPy[slot], py, alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &sumPz[slot], pz, alpaka::hierarchy::Threads{});
            alpaka::atomicAdd(acc, &jetsNonZS.numberOfDaughters()[slot], 1u, alpaka::hierarchy::Threads{});
          }
        }

        alpaka::syncBlockThreads(acc);

        // Step 7: finalize jet kinematics
        if (once_per_block(acc)) {

          for (uint32_t j = 0; j < nroots; ++j) {

            uint32_t slot = begin + j;

            float px = sumPx[slot];
            float py = sumPy[slot];
            float pz = sumPz[slot];

            float pt  = alpaka::math::sqrt(acc, px * px + py * py);
            float phi = alpaka::math::atan2(acc, py, px);

            float eta = 0.f;
            if (pt > 0.f)
              eta = alpaka::math::asinh(acc, pz / pt);

            jetsNonZS.pt()[slot]  = pt;
            jetsNonZS.eta()[slot] = eta;
            jetsNonZS.phi()[slot] = phi;
            jetsNonZS.cluster()[slot] = j;
          }

          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx];
          jetBxLookup.offsets()[block_idx + 1] = nroots;

          alpaka::atomicAdd(acc, nJetsTotal, nroots, alpaka::hierarchy::Blocks{});
        }

        alpaka::syncBlockThreads(acc);
      }
    }
  };

  class JetToAssociationMapKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  OffsetsSoA::ConstView bxLookup,
                                  OffsetsSoA::ConstView jetBxLookup,
                                  OffsetsSoA::ConstView clusterdParticleOffsets,
                                  ClustersDeviceCollection::ConstView clusters,
                                  const unsigned int njets,
                                  const unsigned int nclustered,
                                  uint32_t* key,
                                  uint16_t* idx,
                                  IndexSoA::View map) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      unsigned int nparticles = clusters.metadata().size();
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
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
      radixSortMulti(acc, key, idx, bxLookup.offsets().data(), nullptr);

      // rearrange jet constituents
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;

        uint32_t jetsBegin = jetBxLookup.offsets()[block_idx];
        uint32_t jetsEnd = jetBxLookup.offsets()[block_idx + 1];
        assert(jetsEnd >= jetsBegin);
        assert(jetsEnd <= njets);
        assert(jetBxLookup.offsets()[grid_dim] == njets);
        assert(clusterdParticleOffsets.metadata().size() == int(njets + 1));
        assert(clusterdParticleOffsets.offsets()[njets] == nclustered);

        uint32_t clusteredBegin = clusterdParticleOffsets.offsets()[jetsBegin];
        uint32_t clusteredEnd = clusterdParticleOffsets.offsets()[jetsEnd];

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
          map.indexes()[clusteredBegin + tid] = icand;
          assert(clusteredBegin + tid < nclustered);
          assert(map.metadata().size() == int(nclustered));
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

  ///////////////////////////
  // iterative Seeded Cone Algorithm (greedy, L1T style)
  ///////////////////////////
  // repeat njets times: find highest-pT particle; form jet in cone around it; remove particles from event
  class JetIterKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  ClusterObjDeviceCollection::View puppi, // current particle buffer
                                  OffsetsSoA::ConstView bxLookup, // BX particle ranges
                                  BxIndexSoA::ConstView bxIndex, // BX number for each block
                                  float R2,
                                  unsigned int nIters,
                                  ClustersDeviceCollection::View clusters, // particle -> jet assignment storage
                                  uint32_t* tag, // tmp array for prefix scan compaction
                                  ClusterObjDeviceCollection::View work2, // 2nd particle buffer
                                  ClusterObjDeviceCollection::View jets, // tmp jets
                                  OffsetsSoA::View jetBxLookup, // jet counts per BX
                                  BxIndexSoA::View jetBxIndex, // BX numbers for jets
                                  unsigned int* nJetsTotal) const {

      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      uint32_t* ws = nullptr;
      [[maybe_unused]] constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<uint32_t>(acc);
      }

      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];

        if (end <= begin)
          continue;

        auto& size = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
        size = end - begin;

#ifdef L1TSC_VERBOSE_DEBUG
        if (once_per_block(acc) && (block_idx <= 2))
          printf("In BX %u begin with %u PF candidates: \n", block_idx + 1, end - begin);
#endif
        auto& seed_pt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_eta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_phi = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_i = alpaka::declareSharedVar<unsigned int, __COUNTER__>(acc);

        auto& sum_pt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_eta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_phi = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_dau = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

        unsigned int iter = 0;
        for (iter = 0; iter < nIters; ++iter) {
          bool even = (iter % 2 == 0);
          auto pt = even ? puppi.pt() : work2.pt();
          auto eta = even ? puppi.eta() : work2.eta();
          auto phi = even ? puppi.phi() : work2.phi();
          auto cluster = even ? puppi.cluster() : work2.cluster();
          auto pt2 = !even ? puppi.pt() : work2.pt();
          auto eta2 = !even ? puppi.eta() : work2.eta();
          auto phi2 = !even ? puppi.phi() : work2.phi();
          auto cluster2 = !even ? puppi.cluster() : work2.cluster();

          if (once_per_block(acc)) {
            float spt = 0, seta = 0, sphi = 0;
            unsigned int iseed = end;
            for (unsigned int j = begin, myend = begin + size; j < myend; ++j) {
              if (pt[j] > spt) {
                spt = pt[j];
                seta = eta[j];
                sphi = phi[j];
                iseed = j;
              }
            }
            seed_pt = spt;
            seed_eta = seta;
            seed_phi = sphi;
            seed_i = iseed;
            sum_pt = 0;
            sum_eta = 0;
            sum_phi = 0;
            sum_dau = 0;
          }

          alpaka::syncBlockThreads(acc);

          if (seed_pt == 0)
            break;

          for (uint32_t tid : independent_group_elements(acc, size)) {
            auto ipart = tid + begin;
            float deta = eta[ipart] - seed_eta;
            float dphi = cms::alpakatools::deltaPhi(acc, phi[ipart], seed_phi);
            float dr2 = deta * deta + dphi * dphi;
            if (dr2 < R2) {
              clusters.is_seed()[cluster[ipart]] = (ipart == seed_i ? 1 : 0);
              clusters.cluster()[cluster[ipart]] = iter;
              tag[ipart] = 0;
              alpaka::atomicAdd(acc, &sum_pt, pt[ipart], alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sum_eta, deta * pt[ipart], alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sum_phi, dphi * pt[ipart], alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sum_dau, 1u, alpaka::hierarchy::Threads{});
            } else {
              tag[ipart] = 1;
            }
          }

          alpaka::syncBlockThreads(acc);

          if (once_per_block(acc)) {
            jets.pt()[begin + iter] = sum_pt;
            jets.eta()[begin + iter] = seed_eta + sum_eta / sum_pt;
            jets.phi()[begin + iter] = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);
            jets.numberOfDaughters()[begin + iter] = sum_dau;
            jets.cluster()[begin + iter] = iter;
          }

          blockPrefixScan(acc, tag + begin, size, ws);

          for (uint32_t tid : independent_group_elements(acc, size)) {
            auto ipart = tid + begin;
            if (tag[ipart] > (tid == 0 ? 0 : tag[ipart - 1])) {
              int dest = begin + tag[ipart] - 1;
              pt2[dest] = pt[ipart];
              eta2[dest] = eta[ipart];
              phi2[dest] = phi[ipart];
              cluster2[dest] = cluster[ipart];
            }
          }

          if (once_per_block(acc)) {
            size = tag[begin + size - 1];
          }

          alpaka::syncBlockThreads(acc);

        }  // iter loop

        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx];
          jetBxLookup.offsets()[block_idx + 1] = iter;
          alpaka::atomicAdd(acc, nJetsTotal, iter, alpaka::hierarchy::Blocks{});
        }

      }  // BX loop
    }
  };

  L1TScPhase2SCJetsKernels::L1TScPhase2SCJetsKernels() {}

  // ----------------------------------------------------------------------
  // legacy single-radius non-iterative seeded cone
  // ----------------------------------------------------------------------
  std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::run(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float R2,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;

    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    unsigned int npart = src.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    auto work = ClusterObjDeviceCollection(npart, queue);
    auto jetsNonZS = ClusterObjDeviceCollection(npart, queue);
    clusters.zeroInitialise(queue);
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernel{},
                        src.const_view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        bxLookup.const_view<BxIndexSoA>(),
                        std::sqrt(R2),
                        R2,
                        R2,
                        h_key_device.data(),
                        h_idx_device.data(),
                        work.view(),
                        clusters.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view<OffsetsSoA>(),
                        jetBxLookup.view<BxIndexSoA>(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

  // ----------------------------------------------------------------------
  // finalize()
  // ----------------------------------------------------------------------
  std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice>
  L1TScPhase2SCJetsKernels::finalize(Queue& queue,
                                     const PuppiDeviceCollection& src,
                                     const BxLookupDeviceCollection& bxLookup,
                                     const ClustersDeviceCollection& clusters,
                                     const CounterDevice& nJetsTotalDevice,
                                     const ClusterObjDeviceCollection& jetsNonZS,
                                     BxLookupDeviceCollection& jetBxLookup) const {

    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto nJetsTotalHost = CounterHost(queue);
    alpaka::memcpy(queue, nJetsTotalHost.buffer(), nJetsTotalDevice.buffer());
    alpaka::wait(queue);
    auto njets = nJetsTotalHost.value();

    auto jets = ClusterObjDeviceCollection(njets, queue);

    auto pc = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D{1});
    alpaka::memset(queue, pc, 0x00);

    uint32_t jet_threads_per_block = 1024;
    uint32_t jet_blocks_per_grid =
        cms::alpakatools::divide_up_by(jetBxLookup.view<OffsetsSoA>().metadata().size(), jet_threads_per_block);
    auto jet_grid = cms::alpakatools::make_workdiv<Acc1D>(jet_blocks_per_grid, jet_threads_per_block);

    alpaka::exec<Acc1D>(queue,
                        jet_grid,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        jetBxLookup.view<OffsetsSoA>().offsets().data(),
                        jetBxLookup.view<OffsetsSoA>().offsets().data(),
                        jetBxLookup.view<OffsetsSoA>().metadata().size(),
                        jet_blocks_per_grid,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    auto nClusteredDevice = CounterDevice(queue);
    nClusteredDevice.zeroInitialise(queue);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetZSKernel{},
                        bxLookup.const_view<OffsetsSoA>(),
                        jetsNonZS.const_view(),
                        jetBxLookup.const_view<OffsetsSoA>(),
                        jets.view(),
                        nClusteredDevice.data());

    auto nClusteredHost = CounterHost(queue);
    alpaka::memcpy(queue, nClusteredHost.buffer(), nClusteredDevice.buffer());
    alpaka::wait(queue);
    unsigned int nclustered = nClusteredHost.value();

    auto map = AssociationMapDevice({{int(nclustered), int(njets + 1)}}, queue);
    map.zeroInitialise(queue);

    uint32_t constit_threads_per_block = 1024;
    uint32_t constit_blocks_per_grid = cms::alpakatools::divide_up_by(njets, constit_threads_per_block);
    auto constit_grid = cms::alpakatools::make_workdiv<Acc1D>(constit_blocks_per_grid, constit_threads_per_block);
    uint32_t constit_blocks_per_grid1 = cms::alpakatools::divide_up_by(njets + 1, constit_threads_per_block);
    auto constit_grid1 = cms::alpakatools::make_workdiv<Acc1D>(constit_blocks_per_grid1, constit_threads_per_block);

    alpaka::exec<Acc1D>(
        queue,
        constit_grid,
        [] ALPAKA_FN_ACC(Acc1D const& acc, ClusterObjDeviceCollection::ConstView jets, OffsetsSoA::View offsets) {
          if (cms::alpakatools::once_per_grid(acc))
            offsets.offsets()[0] = 0;
          for (int32_t idx : cms::alpakatools::uniform_elements(acc, offsets.metadata().size() - 1)) {
            offsets.offsets()[idx + 1] = jets.numberOfDaughters()[idx];
          }
        },
        jets.const_view(),
        map.view<OffsetsSoA>());

    alpaka::memset(queue, pc, 0x00);

    alpaka::exec<Acc1D>(queue,
                        constit_grid1,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        map.view<OffsetsSoA>().offsets().data(),
                        map.view<OffsetsSoA>().offsets().data(),
                        njets + 1,
                        constit_blocks_per_grid1,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    unsigned int npart = clusters.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetToAssociationMapKernel{},
                        src.const_view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        jetBxLookup.const_view<OffsetsSoA>(),
                        map.const_view<OffsetsSoA>(),
                        clusters.const_view(),
                        njets,
                        nclustered,
                        h_key_device.data(),
                        h_idx_device.data(),
                        map.view<IndexSoA>());

    return std::make_tuple(std::move(jetBxLookup), std::move(jets), std::move(map));
  }

  // ----------------------------------------------------------------------
  // SCGreedy = iterative seeded cone
  // ----------------------------------------------------------------------
  std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::run(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float R2,
      unsigned int nJets,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    unsigned int npf = src.const_view().metadata().size();

    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto work = ClusterObjDeviceCollection(npf, queue);
    auto work2 = ClusterObjDeviceCollection(npf, queue);

    auto h_tag_device = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npf));
    alpaka::memset(queue, h_tag_device, 0x00);

    uint32_t threads_per_flatblock = 1024;
    uint32_t blocks_per_flatgrid = cms::alpakatools::divide_up_by(npf, threads_per_flatblock);
    auto flatgrid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_flatgrid, threads_per_flatblock);

    auto jetsNonZS = ClusterObjDeviceCollection(npf, queue);
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    alpaka::exec<Acc1D>(
        queue,
        flatgrid,
        [] ALPAKA_FN_ACC(Acc1D const& acc,
                         PuppiDeviceCollection::ConstView puppi,
                         ClusterObjDeviceCollection::View work,
                         ClustersDeviceCollection::View clusters) {
          for (int32_t idx : cms::alpakatools::uniform_elements(acc, clusters.metadata().size())) {
            work.pt()[idx] = puppi.pt()[idx];
            work.eta()[idx] = puppi.eta()[idx];
            work.phi()[idx] = puppi.phi()[idx];
            work.cluster()[idx] = idx;
            clusters.cluster()[idx] = -1;
            clusters.is_seed()[idx] = 0;
          }
        },
        src.const_view(),
        work.view(),
        clusters.view());

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetIterKernel{},
                        work.view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        bxLookup.const_view<BxIndexSoA>(),
                        R2,
                        nJets,
                        clusters.view(),
                        h_tag_device.data(),
                        work2.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view<OffsetsSoA>(),
                        jetBxLookup.view<BxIndexSoA>(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

  // ----------------------------------------------------------------------
  // SCNMS
  //
  // Run the old non-iterative seeded-cone kernel with split radii:
  //   RSeed -> seed competition + old centroid accumulation
  //   RClu  -> final particle assignment
  //
  // RCen / minSeedPt / nCentroidIters are intentionally not part of this
  // implementation, because this function is meant to preserve the old
  // SCNMS logic as closely as possible.
  // ----------------------------------------------------------------------
  L1TScPhase2SCJetsKernels::return_type
  L1TScPhase2SCJetsKernels::runSCNMS(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float RSeed2,
      float RClu2,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;

    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    unsigned int npart = src.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    auto work = ClusterObjDeviceCollection(npart, queue);
    auto jetsNonZS = ClusterObjDeviceCollection(npart, queue);
    clusters.zeroInitialise(queue);
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernel{},
                        src.const_view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        bxLookup.const_view<BxIndexSoA>(),
                        std::sqrt(RSeed2),
                        RSeed2,
                        RClu2,
                        h_key_device.data(),
                        h_idx_device.data(),
                        work.view(),
                        clusters.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view<OffsetsSoA>(),
                        jetBxLookup.view<BxIndexSoA>(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

  // ----------------------------------------------------------------------
  // SCNMSWeighted
  //
  // Run the old weighted non-iterative seeded-cone kernel with split radii:
  //   RSeed -> seed competition + old centroid accumulation
  //   RClu  -> final particle assignment
  //
  // The weighted assignment metric is kept identical to the old weighted
  // kernel:
  //   metric = dr2 / pt_jet^2
  //
  // RCen / alphaSeed / minSeedPt / nCentroidIters are intentionally not
  // part of this implementation, because this function is meant to preserve
  // the old weighted SCNMS logic as closely as possible.
  // ----------------------------------------------------------------------
  L1TScPhase2SCJetsKernels::return_type
  L1TScPhase2SCJetsKernels::runSCNMSWeighted(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float RSeed2,
      float RClu2,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;

    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    unsigned int npart = src.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    auto work = ClusterObjDeviceCollection(npart, queue);
    auto jetsNonZS = ClusterObjDeviceCollection(npart, queue);
    clusters.zeroInitialise(queue);
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernelWeighted{},
                        src.const_view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        bxLookup.const_view<BxIndexSoA>(),
                        std::sqrt(RSeed2),
                        RSeed2,
                        RClu2,
                        h_key_device.data(),
                        h_idx_device.data(),
                        work.view(),
                        clusters.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view<OffsetsSoA>(),
                        jetBxLookup.view<BxIndexSoA>(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

  // ----------------------------------------------------------------------
  // SCNMSWeightedMultiIter = newer weighted multi-radius NMS-style kernel
  // ----------------------------------------------------------------------
  L1TScPhase2SCJetsKernels::return_type
  L1TScPhase2SCJetsKernels::runSCNMSWeightedMultiIter(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float RSeed2,
      float RCen2,
      float RClu2,
      float alphaSeed,
      float minSeedPt,
      unsigned int nCentroidIters,
      ClustersDeviceCollection& clusters) const {

    unsigned int nbx   = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();

    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid   = nbx;
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto uieta_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto idx_device   = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    alpaka::memset(queue, idx_device, 0x00);

    auto seedWeight_device = alpaka::allocAsyncBuf<float, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto seedPos_device    = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    alpaka::memset(queue, seedWeight_device, 0x00);
    alpaka::memset(queue, seedPos_device, 0x00);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    auto work      = ClusterObjDeviceCollection(npart, queue);
    auto jetsNonZS = ClusterObjDeviceCollection(npart, queue);
    jetsNonZS.zeroInitialise(queue);
    clusters.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    float RSeed = std::sqrt(RSeed2);
    float RCen  = std::sqrt(RCen2);
    float RClu  = std::sqrt(RClu2);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernelNMS{},
                        src.const_view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        bxLookup.const_view<BxIndexSoA>(),
                        RSeed,
                        RSeed2,
                        RCen,
                        RCen2,
                        RClu,
                        RClu2,
                        alphaSeed,
                        minSeedPt,
                        nCentroidIters,
                        true,                    // weighted assignment
                        uieta_device.data(),
                        idx_device.data(),
                        seedWeight_device.data(),
                        seedPos_device.data(),
                        work.view(),
                        clusters.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view<OffsetsSoA>(),
                        jetBxLookup.view<BxIndexSoA>(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

  // ----------------------------------------------------------------------
  // LinkTree
  // ----------------------------------------------------------------------
  L1TScPhase2SCJetsKernels::return_type
  L1TScPhase2SCJetsKernels::runLinkTree(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float RLink2,
      float ptMin,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx   = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();

    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid   = nbx;
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto uieta_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto idx_device   = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    alpaka::memset(queue, idx_device, 0x00);

    auto parent_device  = alpaka::allocAsyncBuf<int32_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto root_device    = alpaka::allocAsyncBuf<int32_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto rootCid_device = alpaka::allocAsyncBuf<int32_t, Idx>(queue, cms::alpakatools::Vec1D(npart));

    auto sumPx_device = alpaka::allocAsyncBuf<float, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto sumPy_device = alpaka::allocAsyncBuf<float, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto sumPz_device = alpaka::allocAsyncBuf<float, Idx>(queue, cms::alpakatools::Vec1D(npart));

    alpaka::memset(queue, parent_device,  0xFF);
    alpaka::memset(queue, root_device,    0xFF);
    alpaka::memset(queue, rootCid_device, 0xFF);
    alpaka::memset(queue, sumPx_device,   0x00);
    alpaka::memset(queue, sumPy_device,   0x00);
    alpaka::memset(queue, sumPz_device,   0x00);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    auto work      = ClusterObjDeviceCollection(npart, queue);
    auto jetsNonZS = ClusterObjDeviceCollection(npart, queue);
    jetsNonZS.zeroInitialise(queue);
    clusters.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    float RLink = std::sqrt(RLink2);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernelLinkTree{},
                        src.const_view(),
                        bxLookup.const_view<OffsetsSoA>(),
                        bxLookup.const_view<BxIndexSoA>(),
                        RLink,
                        RLink2,
                        ptMin,
                        uieta_device.data(),
                        idx_device.data(),
                        parent_device.data(),
                        root_device.data(),
                        rootCid_device.data(),
                        sumPx_device.data(),
                        sumPy_device.data(),
                        sumPz_device.data(),
                        work.view(),
                        clusters.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view<OffsetsSoA>(),
                        jetBxLookup.view<BxIndexSoA>(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels