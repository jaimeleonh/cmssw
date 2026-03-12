#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"

#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"

//#define L1TSC_VERBOSE_DEBUG

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  // original non-iter clustering (device kernel functor)

  // class JetKernel {
  // public:
  //   template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
  //   ALPAKA_FN_ACC void operator()(TAcc const& acc,
  //                                 PuppiDeviceCollection::ConstView puppi, // original input p.
  //                                 OffsetsSoA::ConstView bxLookup, // p. ranges per BX
  //                                 BxIndexSoA::ConstView bxIndex, // BX label per BX
  //                                 float R,
  //                                 float R2,
  //                                 uint16_t* uieta, // tmp eta sort key
  //                                 uint16_t* idx, // tmp sort index buffer
  //                                 ClusterObjDeviceCollection::View work, // eta sorted working p. cp
  //                                 ClustersDeviceCollection::View clusters, // output p. assignment info
  //                                 ClusterObjDeviceCollection::View jets, // tmp jet collection
  //                                 OffsetsSoA::View jetBxLookup, // number jets per BX
  //                                 BxIndexSoA::View jetBxIndex, // BX labels for jet collection
  //                                 unsigned int* nJetsTotal) const {
  //     // check whether this backend runs 1 thread per block (for optimization with seed veto markers later)
  //     constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;
  //     // init first lookup entry to 0
  //     if (cms::alpakatools::once_per_grid(acc))
  //       jetBxLookup.offsets()[0] = 0;
  //     // get number of BX blocks
  //     uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
  //     // step-0: sort by eta in blocks
  //     // BX loop
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       // get particle range for this BX
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin)
  //         continue;
  //       uint32_t block_dim = end - begin;
  //       // for each p. in BX build quantized eta value (map form [-5,5] to [0,65535]) in uieta & store its BX index in idx
  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         uieta[tid + begin] = (puppi.eta()[tid + begin] + 5.f) * (std::numeric_limits<uint16_t>::max() / 10.0f);
  //         idx[tid + begin] = tid;  // this is important for the GPU implementation of radixSort
  //       }
  //     }
  //     // sort p. (uieta,idx) by eta inside each BX
  //     radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

  //     // step-1: rearrange
  //     // BX loop
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       // get p. range for this BX
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin)
  //         continue;
  //       uint32_t block_dim = end - begin;
  //       // copy eta sorted p. into work
  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         auto ipart = tid + begin;  // global index
  //         auto isrc = idx[ipart] + begin;
  //         work.pt()[ipart] = puppi.pt()[isrc];
  //         work.eta()[ipart] = puppi.eta()[isrc];
  //         work.phi()[ipart] = puppi.phi()[isrc];
  //         work.cluster()[ipart] = isrc;
  //       }
  //     }

  //     // per BX shared counter nseeds
  //     auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
  //     if (once_per_block(acc)) // one thread set it to 0
  //       nseeds = 0;
  //     alpaka::syncBlockThreads(acc);

  //     // step-2: seed finding
  //     // bx loop
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       // get p. range for this BX
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end = bxLookup.offsets()[block_idx + 1];
  //       // skip if malformed or empty
  //       if (end <= begin)
  //         continue;
  //       uint32_t block_dim = end - begin;
        
  //       // pre-cluster
  //       // loop over every p. as seed candidate
  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         // define cand seed; init jet sums
  //         uint32_t iseed = tid + begin, icluster = work.cluster()[iseed];  // global index
  //         float seed_pt = work.pt()[iseed], seed_eta = work.eta()[iseed], seed_phi = work.phi()[iseed];
  //         float sum_pt = seed_pt, sum_eta = 0, sum_phi = 0;
  //         bool is_seed = true;
  //         // on single-thread backends we can skip already vetoed seeds
  //         if constexpr (single_thread) {
  //           if (clusters.is_seed()[icluster] == -1) 
  //             continue;
  //         }
  //         // scan downward in eta
  //         for (uint32_t j = tid; j > 0; --j) {
  //           uint32_t ipart = j - 1 + begin;  // global index
  //           float deta = work.eta()[ipart] - seed_eta;
  //           if (deta < -R)  // since list sorted in eta, break
  //             break;
  //           float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
  //           if (deta * deta + dphi * dphi < R2) {
  //             if (work.pt()[ipart] >= seed_pt) {  // here we use >=, since we're for indices above ipart
  //               is_seed = false;
  //               break;
  //             } else {
  //               if constexpr (single_thread) {
  //                 clusters.is_seed()[work.cluster()[ipart]] = -1;
  //               }
  //               sum_pt += work.pt()[ipart];
  //               sum_eta += work.pt()[ipart] * deta;
  //               sum_phi += work.pt()[ipart] * dphi;
  //             }
  //           }
  //         }
  //         // scan upward in eta
  //         for (uint32_t j = tid + 1; j < block_dim; ++j) {
  //           uint32_t ipart = j + begin;  // global index
  //           float deta = work.eta()[ipart] - seed_eta;
  //           if (deta > R)  // sorted in eta, so we can stop here
  //             break;
  //           float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
  //           if (deta * deta + dphi * dphi < R2) {
  //             if (work.pt()[ipart] > seed_pt) {  // note: here we use > instead of >=
  //               is_seed = false;
  //               break;
  //             } else {
  //               if constexpr (single_thread) {
  //                 clusters.is_seed()[work.cluster()[ipart]] = -1;
  //               }
  //               sum_pt += work.pt()[ipart];
  //               sum_eta += work.pt()[ipart] * deta;
  //               sum_phi += work.pt()[ipart] * dphi;
  //             }
  //           }
  //         }
  //         // convert centroid sums (offsets) to jet axis
  //         sum_eta = seed_eta + sum_eta / sum_pt;
  //         sum_phi = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);

  //         // seed cand survived
  //         if (is_seed) {
  //           // increment shared seed counter
  //           auto ijet = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{}) + begin;
  //           // store cand jet into jets
  //           jets.pt()[ijet] = sum_pt;
  //           jets.eta()[ijet] = sum_eta;
  //           jets.phi()[ijet] = sum_phi;
  //           jets.cluster()[ijet] = ijet - begin;
  //           jets.numberOfDaughters()[ijet] = 0;
  //           clusters.is_seed()[icluster] = 1;
  //         }
  //       }
  //       // sync after seed finding
  //       alpaka::syncBlockThreads(acc);

  //       // step 3: reassociate p. to neareast accepted seed
  //       // loop over all eta-sorted p. (parallel)
  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         auto ipart = tid + begin;                   // global index
  //         uint32_t icluster = work.cluster()[ipart];  // original index
  //         // init nearest search state
  //         float nearest = R2;
  //         int jcluster = -1;
  //         // loop over accepted seeds
  //         for (uint32_t j = 0; j < nseeds; ++j) {
  //           auto jseed = j + begin;  // global index
  //           float deta = work.eta()[ipart] - jets.eta()[jseed];
  //           float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], jets.phi()[jseed]);
  //           float dr2 = deta * deta + dphi * dphi;
  //           // keep neareast inside R
  //           if (dr2 < nearest) {
  //             jcluster = jets.cluster()[jseed]; // local (in this BX) jet assignment for this p.
  //             nearest = dr2;
  //           }
  //         }

  //         // increment daughter count of assigned jet
  //         if (jcluster != -1) {
  //           alpaka::atomicAdd(acc, &jets.numberOfDaughters()[jcluster + begin], 1u, alpaka::hierarchy::Threads{});
  //         }

  //         // write final p. assignment using original p. index (icluster) 
  //         clusters.cluster()[icluster] = jcluster;

  //         // clear veto markers (single-threads backends) -> clusters.is_seed() 1 for seed, 0 else
  //         if constexpr (single_thread) {
  //           clusters.is_seed()[icluster] = std::max(clusters.is_seed()[icluster], 0); // clear -1 markers
  //         }
  //       }

  //       // store per BX output metadata (1 thread per block)
  //       if (once_per_block(acc)) {
  //         jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx]; // cp BX label from input p. to output jet BX structure
  //         jetBxLookup.offsets()[block_idx + 1] = nseeds; // store number of jets in BX
  //         // add this BX's jet count to total (over all BX)
  //         alpaka::atomicAdd(acc, nJetsTotal, nseeds, alpaka::hierarchy::Blocks{});
  //       }
  //       // per block consistency before finishing BX loop
  //       alpaka::syncBlockThreads(acc);
  //     }  // BX block
  //   }  // operator()
  // };  // class



  // non-iter clustering with weighted (device kernel functor)

  class JetKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi, // original input p.
                                  OffsetsSoA::ConstView bxLookup, // p. ranges per BX
                                  BxIndexSoA::ConstView bxIndex, // BX label per BX
                                  float R,
                                  float R2,
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
          auto ipart = tid + begin;  // global index
          auto isrc = idx[ipart] + begin;
          work.pt()[ipart] = puppi.pt()[isrc];
          work.eta()[ipart] = puppi.eta()[isrc];
          work.phi()[ipart] = puppi.phi()[isrc];
          work.cluster()[ipart] = isrc;
        }
      }

      // per BX shared counter nseeds
      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      if (once_per_block(acc)) // one thread set it to 0
        nseeds = 0;
      alpaka::syncBlockThreads(acc);

      // step-2: seed finding
      // bx loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get p. range for this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        // skip if malformed or empty
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;
        
        // pre-cluster
        // loop over every p. as seed candidate
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          // define cand seed; init jet sums
          uint32_t iseed = tid + begin, icluster = work.cluster()[iseed];  // global index
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
            uint32_t ipart = j - 1 + begin;  // global index
            float deta = work.eta()[ipart] - seed_eta;
            if (deta < -R)  // since list sorted in eta, break
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < R2) {
              if (work.pt()[ipart] >= seed_pt) {  // here we use >=, since we're for indices above ipart
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
            uint32_t ipart = j + begin;  // global index
            float deta = work.eta()[ipart] - seed_eta;
            if (deta > R)  // sorted in eta, so we can stop here
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < R2) {
              if (work.pt()[ipart] > seed_pt) {  // note: here we use > instead of >=
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
          sum_eta = seed_eta + sum_eta / sum_pt;
          sum_phi = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);

          // seed cand survived
          if (is_seed) {
            // increment shared seed counter
            auto ijet = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{}) + begin;
            // store cand jet into jets
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

        // step 3: reassociate p. to neareast accepted seed
        // loop over all eta-sorted p. (parallel)
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;                   // global index
          uint32_t icluster = work.cluster()[ipart];  // original index
          // init nearest search state
          float nearest = 1e30f;
          int jcluster = -1;
          // loop over accepted seeds
          for (uint32_t j = 0; j < nseeds; ++j) {
            auto jseed = j + begin;  // global index
            float deta = work.eta()[ipart] - jets.eta()[jseed];
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], jets.phi()[jseed]);
            float dr2 = deta * deta + dphi * dphi;
            // keep neareast inside R, weighted by jet pt (anti-kt-like)
            if (dr2 < R2) {
              float jet_pt = jets.pt()[jseed];
              float metric = dr2 / (jet_pt * jet_pt + 1e-12f);
              if (metric < nearest) {
                jcluster = jets.cluster()[jseed]; // local (in this BX) jet assignment for this p.
                nearest = metric;
              }
            }
          }

          // increment daughter count of assigned jet
          if (jcluster != -1) {
            alpaka::atomicAdd(acc, &jets.numberOfDaughters()[jcluster + begin], 1u, alpaka::hierarchy::Threads{});
          }

          // write final p. assignment using original p. index (icluster) 
          clusters.cluster()[icluster] = jcluster;

          // clear veto markers (single-threads backends) -> clusters.is_seed() 1 for seed, 0 else
          if constexpr (single_thread) {
            clusters.is_seed()[icluster] = std::max(clusters.is_seed()[icluster], 0); // clear -1 markers
          }
        }

        // store per BX output metadata (1 thread per block)
        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx]; // cp BX label from input p. to output jet BX structure
          jetBxLookup.offsets()[block_idx + 1] = nseeds; // store number of jets in BX
          // add this BX's jet count to total (over all BX)
          alpaka::atomicAdd(acc, nJetsTotal, nseeds, alpaka::hierarchy::Blocks{});
        }
        // per block consistency before finishing BX loop
        alpaka::syncBlockThreads(acc);
      }  // BX block
    }  // operator()
  };  // class





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
        // get event range
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
          // simpler to do it on a single thread
          unsigned int nClusteredBlock = 0;
          for (uint32_t iJet = beginDst; iJet < endDst; ++iJet) {
            nClusteredBlock += jets.numberOfDaughters()[iJet];
          }
          alpaka::atomicAdd(acc, nClustered, nClusteredBlock, alpaka::hierarchy::Blocks{});
        }  // once per block
      }  // block loop
    }  // operator()
  };  // class







  // class JetKernelNMSWeighted {
  // public:
  //   template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
  //   ALPAKA_FN_ACC void operator()(
  //       TAcc const& acc,
  //       PuppiDeviceCollection::ConstView puppi,
  //       OffsetsSoA::ConstView bxLookup,
  //       BxIndexSoA::ConstView bxIndex,
  //       float RSeed,
  //       float RSeed2,
  //       float RCen,
  //       float RCen2,
  //       float RClu,
  //       float RClu2,
  //       float alphaSeed,
  //       float minSeedPt,
  //       unsigned int nCentroidIters,
  //       uint16_t* uieta,
  //       uint16_t* idx,
  //       float* seedWeight,
  //       uint16_t* seedPos,
  //       ClusterObjDeviceCollection::View work,
  //       ClustersDeviceCollection::View clusters,
  //       ClusterObjDeviceCollection::View jetsNonZS,
  //       OffsetsSoA::View jetBxLookup,
  //       BxIndexSoA::View jetBxIndex,
  //       unsigned int* nJetsTotal) const
  //   {
  //     if (cms::alpakatools::once_per_grid(acc))
  //       jetBxLookup.offsets()[0] = 0;

  //     uint32_t grid_dim =
  //       alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

  //     // ============================================================
  //     // Step 0: sort PF by eta (per BX)
  //     // ============================================================
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end   = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin) continue;

  //       uint32_t block_dim = end - begin;

  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         uint32_t i = tid + begin;
  //         uieta[i] =
  //           (puppi.eta()[i] + 5.f) *
  //           (std::numeric_limits<uint16_t>::max() / 10.f);
  //         idx[i] = tid;
  //       }
  //     }

  //     radixSortMulti(acc, uieta, idx,
  //                   bxLookup.offsets().data(), nullptr);

  //     // ============================================================
  //     // Step 1: rearrange into eta-sorted work array
  //     // ============================================================
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end   = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin) continue;

  //       uint32_t block_dim = end - begin;

  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         uint32_t ip   = tid + begin;
  //         uint32_t isrc = idx[ip] + begin;

  //         work.pt()[ip]      = puppi.pt()[isrc];
  //         work.eta()[ip]     = puppi.eta()[isrc];
  //         work.phi()[ip]     = puppi.phi()[isrc];
  //         work.cluster()[ip] = isrc;  // original index
  //       }
  //     }

  //     auto& nseeds =
  //       alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

  //     // ============================================================
  //     // BX loop
  //     // ============================================================
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {

  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end   = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin) continue;

  //       uint32_t block_dim = end - begin;

  //       if (once_per_block(acc)) nseeds = 0u;
  //       alpaka::syncBlockThreads(acc);

  //       // ============================================================
  //       // Step 2: seed finding (local pT maxima within RSeed)
  //       // ============================================================
  //       for (uint32_t tid :
  //           independent_group_elements(acc, block_dim)) {

  //         uint32_t iseed = tid + begin;

  //         float seed_pt = work.pt()[iseed];
  //         // if (seed_pt < minSeedPt) continue;

  //         float seed_eta = work.eta()[iseed];
  //         float seed_phi = work.phi()[iseed];
  //         uint32_t seed_orig =
  //           static_cast<uint32_t>(work.cluster()[iseed]);

  //         bool is_seed = true;

  //         // scan down
  //         for (uint32_t j = tid; j > 0; --j) {
  //           uint32_t ip = begin + (j - 1);
  //           float deta = work.eta()[ip] - seed_eta;
  //           if (deta < -RSeed) break;

  //           float dphi =
  //             cms::alpakatools::deltaPhi(acc,
  //                                       work.phi()[ip],
  //                                       seed_phi);

  //           float dr2 = deta*deta + dphi*dphi;

  //           if (dr2 < RSeed2) {
  //             float ptj = work.pt()[ip];
  //             uint32_t origj =
  //               static_cast<uint32_t>(work.cluster()[ip]);

  //             if (ptj > seed_pt ||
  //               (ptj == seed_pt && origj < seed_orig)) {
  //               is_seed = false;
  //               break;
  //             }
  //           }
  //         }

  //         // scan up
  //         if (is_seed) {
  //           for (uint32_t j = tid + 1; j < block_dim; ++j) {
  //             uint32_t ip = begin + j;
  //             float deta = work.eta()[ip] - seed_eta;
  //             if (deta > RSeed) break;

  //             float dphi =
  //               cms::alpakatools::deltaPhi(acc,
  //                                         work.phi()[ip],
  //                                         seed_phi);

  //             float dr2 = deta*deta + dphi*dphi;

  //             if (dr2 < RSeed2) {
  //               float ptj = work.pt()[ip];
  //               uint32_t origj =
  //                 static_cast<uint32_t>(work.cluster()[ip]);

  //               if (ptj > seed_pt ||
  //                 (ptj == seed_pt && origj < seed_orig)) {
  //                 is_seed = false;
  //                 break;
  //               }
  //             }
  //           }
  //         }

  //         if (!is_seed) continue;

  //         uint32_t jlocal =
  //           alpaka::atomicAdd(acc, &nseeds, 1u,
  //                             alpaka::hierarchy::Threads{});

  //         uint32_t jseed = begin + jlocal;

  //         jetsNonZS.pt()[jseed]  = seed_pt;
  //         jetsNonZS.eta()[jseed] = seed_eta;
  //         jetsNonZS.phi()[jseed] = seed_phi;
  //         jetsNonZS.cluster()[jseed] = jlocal;
  //         jetsNonZS.numberOfDaughters()[jseed] = 0u;

  //         seedPos[jseed] = static_cast<uint16_t>(tid);

  //         // float denom =
  //         //   alpaka::math::pow(acc,
  //         //     alpaka::math::max(acc, seed_pt, 1e-6f),
  //         //     alphaSeed);

  //         // seedWeight[jseed] = 1.f / denom;
  //         seedWeight[jseed] = 1.f;

  //         clusters.is_seed()[seed_orig] = 1;
  //         clusters.cluster()[seed_orig] = jlocal;
  //       }

  //       alpaka::syncBlockThreads(acc);

  //       // ============================================================
  //       // Step 3: optional centroid iterations
  //       // ============================================================
  //       if (nCentroidIters > 0u) {

  //         for (uint32_t j :
  //             independent_group_elements(acc, nseeds)) {

  //           uint32_t jseed = begin + j;

  //           float axis_eta = jetsNonZS.eta()[jseed];
  //           float axis_phi = jetsNonZS.phi()[jseed];
  //           uint32_t pos   =
  //             static_cast<uint32_t>(seedPos[jseed]);

  //           for (unsigned int iter = 0;
  //               iter < nCentroidIters;
  //               ++iter) {

  //             float sum_pt  = 0.f;
  //             float sum_eta = 0.f;
  //             float sum_sin = 0.f;
  //             float sum_cos = 0.f;

  //             // scan down
  //             for (int32_t k =
  //                   static_cast<int32_t>(pos);
  //                 k >= 0; --k) {

  //               uint32_t ip =
  //                 begin + static_cast<uint32_t>(k);

  //               float deta =
  //                 work.eta()[ip] - axis_eta;

  //               if (deta < -RCen) break;

  //               float dphi =
  //                 cms::alpakatools::deltaPhi(acc,
  //                                           work.phi()[ip],
  //                                           axis_phi);

  //               float dr2 = deta*deta + dphi*dphi;

  //               if (dr2 < RCen2) {
  //                 float w = work.pt()[ip];
  //                 sum_pt  += w;
  //                 sum_eta += w * work.eta()[ip];

  //                 float ph = work.phi()[ip];
  //                 sum_sin +=
  //                   w * alpaka::math::sin(acc, ph);
  //                 sum_cos +=
  //                   w * alpaka::math::cos(acc, ph);
  //               }
  //             }

  //             // scan up
  //             for (uint32_t k = pos + 1;
  //                 k < block_dim; ++k) {

  //               uint32_t ip = begin + k;

  //               float deta =
  //                 work.eta()[ip] - axis_eta;

  //               if (deta > RCen) break;

  //               float dphi =
  //                 cms::alpakatools::deltaPhi(acc,
  //                                           work.phi()[ip],
  //                                           axis_phi);

  //               float dr2 = deta*deta + dphi*dphi;

  //               if (dr2 < RCen2) {
  //                 float w = work.pt()[ip];
  //                 sum_pt  += w;
  //                 sum_eta += w * work.eta()[ip];

  //                 float ph = work.phi()[ip];
  //                 sum_sin +=
  //                   w * alpaka::math::sin(acc, ph);
  //                 sum_cos +=
  //                   w * alpaka::math::cos(acc, ph);
  //               }
  //             }

  //             if (sum_pt > 0.f) {
  //               axis_eta = sum_eta / sum_pt;
  //               axis_phi =
  //                 alpaka::math::atan2(acc,
  //                                     sum_sin,
  //                                     sum_cos);

  //               jetsNonZS.pt()[jseed]  = sum_pt;
  //               jetsNonZS.eta()[jseed] = axis_eta;
  //               jetsNonZS.phi()[jseed] = axis_phi;
  //             }
  //           }
  //         }
  //       }

  //       alpaka::syncBlockThreads(acc);

  //       // ============================================================
  //       // Step 4: assignment (plain loop over seeds)
  //       // ============================================================
  //       for (uint32_t tid :
  //           independent_group_elements(acc, block_dim)) {

  //         uint32_t ip   = begin + tid;
  //         uint32_t orig =
  //           static_cast<uint32_t>(work.cluster()[ip]);

  //         if (clusters.is_seed()[orig] == 1) {
  //           int jfix = clusters.cluster()[orig];
  //           if (jfix >= 0) {
  //             alpaka::atomicAdd(
  //               acc,
  //               &jetsNonZS.numberOfDaughters()
  //                 [begin + jfix],
  //               1u,
  //               alpaka::hierarchy::Threads{});
  //           }
  //           continue;
  //         }

  //         float eta_i = work.eta()[ip];
  //         float phi_i = work.phi()[ip];

  //         float bestScore =
  //           std::numeric_limits<float>::infinity();

  //         int bestJ = -1;

  //         for (uint32_t j = 0; j < nseeds; ++j) {

  //           uint32_t jseed = begin + j;

  //           float deta =
  //             jetsNonZS.eta()[jseed] - eta_i;

  //           if (deta < -RClu ||
  //               deta >  RClu)
  //             continue;

  //           float dphi =
  //             cms::alpakatools::deltaPhi(acc,
  //                                       jetsNonZS.phi()[jseed],
  //                                       phi_i);

  //           float dr2 = deta*deta + dphi*dphi;

  //           if (dr2 >= RClu2)
  //             continue;

  //           float score =
  //             dr2 * seedWeight[jseed];

  //           if (score < bestScore) {
  //           // if (dr2 < bestScore)
  //             bestScore = score;
  //             bestJ = j;
  //           }
  //         }

  //         clusters.is_seed()[orig] = 0;
  //         clusters.cluster()[orig] = bestJ;

  //         if (bestJ >= 0) {
  //           alpaka::atomicAdd(
  //             acc,
  //             &jetsNonZS.numberOfDaughters()
  //               [begin + bestJ],
  //             1u,
  //             alpaka::hierarchy::Threads{});
  //         }
  //       }

  //       alpaka::syncBlockThreads(acc);

  //       if (once_per_block(acc)) {
  //         jetBxIndex.bx()[block_idx] =
  //           bxIndex.bx()[block_idx];

  //         jetBxLookup.offsets()
  //           [block_idx + 1] = nseeds;

  //         alpaka::atomicAdd(acc,
  //                           nJetsTotal,
  //                           nseeds,
  //                           alpaka::hierarchy::Blocks{});
  //       }

  //       alpaka::syncBlockThreads(acc);
  //     }
  //   }
  // };






  // class JetKernelNMSWeighted {
  // public:
  //   template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
  //   ALPAKA_FN_ACC void operator()(
  //       TAcc const& acc,
  //       PuppiDeviceCollection::ConstView puppi,
  //       OffsetsSoA::ConstView bxLookup,
  //       BxIndexSoA::ConstView bxIndex,
  //       float RSeed,
  //       float RSeed2,
  //       float RCen,
  //       float RCen2,
  //       float RClu,
  //       float RClu2,
  //       float alphaSeed,
  //       float minSeedPt,
  //       unsigned int nCentroidIters,
  //       uint16_t* uieta,
  //       uint16_t* idx,
  //       float* seedWeight,
  //       uint16_t* seedPos,
  //       ClusterObjDeviceCollection::View work,
  //       ClustersDeviceCollection::View clusters,
  //       ClusterObjDeviceCollection::View jetsNonZS,
  //       OffsetsSoA::View jetBxLookup,
  //       BxIndexSoA::View jetBxIndex,
  //       unsigned int* nJetsTotal) const {

  //     constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;

  //     if (cms::alpakatools::once_per_grid(acc))
  //       jetBxLookup.offsets()[0] = 0;

  //     uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

  //     // Step 0: eta sort
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end   = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin) continue;

  //       uint32_t block_dim = end - begin;

  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         uint32_t i = tid + begin;
  //         uieta[i] =
  //           (puppi.eta()[i] + 5.f) *
  //           (std::numeric_limits<uint16_t>::max() / 10.0f);
  //         idx[i] = tid;
  //       }
  //     }

  //     radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

  //     // Step 1: rearrange
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end   = bxLookup.offsets()[block_idx + 1];
  //       if (end <= begin) continue;

  //       uint32_t block_dim = end - begin;

  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {
  //         uint32_t ip   = tid + begin;
  //         uint32_t isrc = idx[ip] + begin;

  //         work.pt()[ip]      = puppi.pt()[isrc];
  //         work.eta()[ip]     = puppi.eta()[isrc];
  //         work.phi()[ip]     = puppi.phi()[isrc];
  //         work.cluster()[ip] = isrc;
  //       }
  //     }

  //     auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

  //     if (once_per_block(acc))
  //       nseeds = 0;

  //     alpaka::syncBlockThreads(acc);

  //     // BX loop
  //     for (uint32_t block_idx : independent_groups(acc, grid_dim)) {

  //       uint32_t begin = bxLookup.offsets()[block_idx];
  //       uint32_t end   = bxLookup.offsets()[block_idx + 1];

  //       if (end <= begin) continue;

  //       uint32_t block_dim = end - begin;

  //       // Seed finding
  //       for (uint32_t tid : independent_group_elements(acc, block_dim)) {

  //         uint32_t iseed = tid + begin;
  //         uint32_t icluster = work.cluster()[iseed];

  //         float seed_pt  = work.pt()[iseed];
  //         if (seed_pt < minSeedPt) continue;

  //         float seed_eta = work.eta()[iseed];
  //         float seed_phi = work.phi()[iseed];

  //         float sum_pt  = seed_pt;
  //         float sum_eta = 0.f;
  //         float sum_phi = 0.f;

  //         bool is_seed = true;

  //         if constexpr (single_thread) {
  //           if (clusters.is_seed()[icluster] == -1)
  //             continue;
  //         }

  //         // scan down
  //         for (uint32_t j = tid; j > 0; --j) {

  //           uint32_t ipart = begin + (j - 1);

  //           float deta = work.eta()[ipart] - seed_eta;
  //           if (deta < -RSeed) break;

  //           float dphi =
  //             cms::alpakatools::deltaPhi(acc,
  //               work.phi()[ipart], seed_phi);

  //           float dr2 = deta*deta + dphi*dphi;

  //           if (dr2 < RSeed2) {

  //             if (work.pt()[ipart] >= seed_pt) {
  //               is_seed = false;
  //               break;
  //             } else {

  //               if constexpr (single_thread) {
  //                 clusters.is_seed()[work.cluster()[ipart]] = -1;
  //               }

  //               sum_pt  += work.pt()[ipart];
  //               sum_eta += work.pt()[ipart] * deta;
  //               sum_phi += work.pt()[ipart] * dphi;
  //             }
  //           }
  //         }

  //         if (!is_seed) continue;

  //         // scan up
  //         for (uint32_t j = tid + 1; j < block_dim; ++j) {

  //           uint32_t ipart = begin + j;

  //           float deta = work.eta()[ipart] - seed_eta;
  //           if (deta > RSeed) break;

  //           float dphi =
  //             cms::alpakatools::deltaPhi(acc,
  //               work.phi()[ipart], seed_phi);

  //           float dr2 = deta*deta + dphi*dphi;

  //           if (dr2 < RSeed2) {

  //             if (work.pt()[ipart] > seed_pt) {
  //               is_seed = false;
  //               break;
  //             } else {

  //               if constexpr (single_thread) {
  //                 clusters.is_seed()[work.cluster()[ipart]] = -1;
  //               }

  //               sum_pt  += work.pt()[ipart];
  //               sum_eta += work.pt()[ipart] * deta;
  //               sum_phi += work.pt()[ipart] * dphi;
  //             }
  //           }
  //         }

  //         if (!is_seed) continue;

  //         float jet_eta = seed_eta + sum_eta / sum_pt;
  //         float jet_phi =
  //           cms::alpakatools::reducePhiRange(
  //             acc,
  //             seed_phi + sum_phi / sum_pt
  //           );

  //         uint32_t jlocal =
  //           alpaka::atomicAdd(acc,
  //             &nseeds,
  //             1u,
  //             alpaka::hierarchy::Threads{});

  //         uint32_t jseed = begin + jlocal;

  //         jetsNonZS.pt()[jseed]  = sum_pt;
  //         jetsNonZS.eta()[jseed] = jet_eta;
  //         jetsNonZS.phi()[jseed] = jet_phi;
  //         jetsNonZS.cluster()[jseed] = jlocal;
  //         jetsNonZS.numberOfDaughters()[jseed] = 0;

  //         float denom =
  //           alpaka::math::pow(
  //             acc,
  //             alpaka::math::max(acc, seed_pt, 1e-6f),
  //             alphaSeed
  //           );

  //         seedWeight[jseed] = 1.f / denom;

  //         clusters.is_seed()[icluster] = 1;
  //       }

  //       alpaka::syncBlockThreads(acc);

  //       // Assignment
  //       for (uint32_t tid :
  //           independent_group_elements(acc, block_dim)) {

  //         uint32_t ipart = tid + begin;
  //         uint32_t icluster = work.cluster()[ipart];

  //         float bestScore =
  //           std::numeric_limits<float>::infinity();

  //         int jcluster = -1;

  //         for (uint32_t j = 0; j < nseeds; ++j) {

  //           uint32_t jseed = begin + j;

  //           float deta =
  //             work.eta()[ipart] -
  //             jetsNonZS.eta()[jseed];

  //           if (deta < -RClu || deta > RClu)
  //             continue;

  //           float dphi =
  //             cms::alpakatools::deltaPhi(
  //               acc,
  //               work.phi()[ipart],
  //               jetsNonZS.phi()[jseed]
  //             );

  //           float dr2 = deta*deta + dphi*dphi;

  //           if (dr2 >= RClu2)
  //             continue;

  //           float score =
  //             dr2 * seedWeight[jseed];

  //           if (score < bestScore) {
  //             bestScore = score;
  //             jcluster  = jetsNonZS.cluster()[jseed];
  //           }
  //         }

  //         if (jcluster != -1) {

  //           alpaka::atomicAdd(
  //             acc,
  //             &jetsNonZS.numberOfDaughters()
  //               [begin + jcluster],
  //             1u,
  //             alpaka::hierarchy::Threads{}
  //           );
  //         }

  //         clusters.cluster()[icluster] = jcluster;

  //         if constexpr (single_thread) {
  //           clusters.is_seed()[icluster] =
  //             std::max(clusters.is_seed()[icluster], 0);
  //         }
  //       }

  //       if (once_per_block(acc)) {

  //         jetBxIndex.bx()[block_idx] =
  //           bxIndex.bx()[block_idx];

  //         jetBxLookup.offsets()[block_idx + 1] =
  //           nseeds;

  //         alpaka::atomicAdd(
  //           acc,
  //           nJetsTotal,
  //           nseeds,
  //           alpaka::hierarchy::Blocks{}
  //         );
  //       }

  //       alpaka::syncBlockThreads(acc);
  //     }
  //   }
  // };






  class JetKernelNMSWeighted {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        PuppiDeviceCollection::ConstView puppi,
        OffsetsSoA::ConstView bxLookup,
        BxIndexSoA::ConstView bxIndex,
        float RSeed,
        float RSeed2,
        float /*RCen*/,
        float /*RCen2*/,
        float /*RClu*/,
        float /*RClu2*/,
        float alphaSeed,
        float /*minSeedPt*/,
        unsigned int /*nCentroidIters*/,
        uint16_t* uieta,
        uint16_t* idx,
        float* seedWeight,
        uint16_t* /*seedPos*/,
        ClusterObjDeviceCollection::View work,
        ClustersDeviceCollection::View clusters,
        ClusterObjDeviceCollection::View jetsNonZS,
        OffsetsSoA::View jetBxLookup,
        BxIndexSoA::View jetBxIndex,
        unsigned int* nJetsTotal) const {

      constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;

      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      // Step 0: eta sort
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uieta[tid + begin] =
            (puppi.eta()[tid + begin] + 5.f) *
            (std::numeric_limits<uint16_t>::max() / 10.0f);
          idx[tid + begin] = tid;
        }
      }

      radixSortMulti(acc, uieta, idx, bxLookup.offsets().data(), nullptr);

      // Step 1: rearrange
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          auto isrc  = idx[ipart] + begin;

          work.pt()[ipart]      = puppi.pt()[isrc];
          work.eta()[ipart]     = puppi.eta()[isrc];
          work.phi()[ipart]     = puppi.phi()[isrc];
          work.cluster()[ipart] = isrc;
        }
      }

      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

      if (once_per_block(acc))
        nseeds = 0;

      alpaka::syncBlockThreads(acc);

      // BX loop
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {

        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end   = bxLookup.offsets()[block_idx + 1];
        if (end <= begin) continue;

        uint32_t block_dim = end - begin;

        // Seed finding
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {

          uint32_t iseed    = tid + begin;
          uint32_t icluster = work.cluster()[iseed];

          float seed_pt  = work.pt()[iseed];
          float seed_eta = work.eta()[iseed];
          float seed_phi = work.phi()[iseed];

          float sum_pt  = seed_pt;
          float sum_eta = 0.f;
          float sum_phi = 0.f;

          bool is_seed = true;

          if constexpr (single_thread) {
            if (clusters.is_seed()[icluster] == -1)
              continue;
          }

          // scan down
          for (uint32_t j = tid; j > 0; --j) {

            uint32_t ipart = j - 1 + begin;

            float deta = work.eta()[ipart] - seed_eta;
            if (deta < -RSeed)
              break;

            float dphi =
              cms::alpakatools::deltaPhi(acc,
                                        work.phi()[ipart],
                                        seed_phi);

            float dr2 = deta*deta + dphi*dphi;

            if (dr2 < RSeed2) {

              if (work.pt()[ipart] >= seed_pt) {
                is_seed = false;
                break;
              } else {

                if constexpr (single_thread)
                  clusters.is_seed()[work.cluster()[ipart]] = -1;

                sum_pt  += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }

          // if (!is_seed) continue;

          // scan up
          for (uint32_t j = tid + 1; j < block_dim; ++j) {

            uint32_t ipart = j + begin;

            float deta = work.eta()[ipart] - seed_eta;
            if (deta > RSeed)
              break;

            float dphi =
              cms::alpakatools::deltaPhi(acc,
                                        work.phi()[ipart],
                                        seed_phi);

            float dr2 = deta*deta + dphi*dphi;

            if (dr2 < RSeed2) {

              if (work.pt()[ipart] > seed_pt) {
                is_seed = false;
                break;
              } else {

                if constexpr (single_thread)
                  clusters.is_seed()[work.cluster()[ipart]] = -1;

                sum_pt  += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }

          if (!is_seed) continue;

          float jet_eta = seed_eta + sum_eta / sum_pt;
          float jet_phi =
            cms::alpakatools::reducePhiRange(
              acc,
              seed_phi + sum_phi / sum_pt
            );

          uint32_t jlocal =
            alpaka::atomicAdd(acc,
                              &nseeds,
                              1u,
                              alpaka::hierarchy::Threads{});

          uint32_t jseed = begin + jlocal;

          jetsNonZS.pt()[jseed]  = sum_pt;
          jetsNonZS.eta()[jseed] = jet_eta;
          jetsNonZS.phi()[jseed] = jet_phi;
          jetsNonZS.cluster()[jseed] = jlocal;
          jetsNonZS.numberOfDaughters()[jseed] = 0;

          // float denom =
          //   alpaka::math::pow(
          //     acc,
          //     alpaka::math::max(acc, seed_pt, 1e-6f),
          //     alphaSeed
          //   );

          // seedWeight[jseed] = 1.f / denom;
          seedWeight[jseed] = 1.f;

          clusters.is_seed()[icluster] = 1;
        }

        alpaka::syncBlockThreads(acc);

        // Assignment (weighted score instead of pure distance)
        for (uint32_t tid :
            independent_group_elements(acc, block_dim)) {

          uint32_t ipart = tid + begin;
          uint32_t icluster = work.cluster()[ipart];

          float bestScore =
            std::numeric_limits<float>::infinity();

          int jcluster = -1;

          for (uint32_t j = 0; j < nseeds; ++j) {

            uint32_t jseed = begin + j;

            float deta =
              work.eta()[ipart] -
              jetsNonZS.eta()[jseed];

            float dphi =
              cms::alpakatools::deltaPhi(
                acc,
                work.phi()[ipart],
                jetsNonZS.phi()[jseed]
              );

            float dr2 = deta*deta + dphi*dphi;

            if (dr2 < RSeed2) {

              float score =
                // dr2 * seedWeight[jseed];
                dr2;

              if (score < bestScore) {
                bestScore = score;
                jcluster  = jetsNonZS.cluster()[jseed];
              }
            }
          }

          if (jcluster != -1) {

            alpaka::atomicAdd(
              acc,
              &jetsNonZS.numberOfDaughters()
                [begin + jcluster],
              1u,
              alpaka::hierarchy::Threads{}
            );
          }

          clusters.cluster()[icluster] = jcluster;

          if constexpr (single_thread) {
            clusters.is_seed()[icluster] =
              std::max(clusters.is_seed()[icluster], 0);
          }
        }

        if (once_per_block(acc)) {

          jetBxIndex.bx()[block_idx] =
            bxIndex.bx()[block_idx];

          jetBxLookup.offsets()[block_idx + 1] =
            nseeds;

          alpaka::atomicAdd(
            acc,
            nJetsTotal,
            nseeds,
            alpaka::hierarchy::Blocks{}
          );
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
        // get event range
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
        // get range of particles in this BX
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];
        if (end <= begin)
          continue;
        // get the offset, i.e. global index, of the first and one-past-the-last jet in this BX
        uint32_t jetsBegin = jetBxLookup.offsets()[block_idx];
        uint32_t jetsEnd = jetBxLookup.offsets()[block_idx + 1];
        assert(jetsEnd >= jetsBegin);
        assert(jetsEnd <= njets);
        assert(jetBxLookup.offsets()[grid_dim] == njets);
        assert(clusterdParticleOffsets.metadata().size() == int(njets + 1));
        assert(clusterdParticleOffsets.offsets()[njets] == nclustered);
        // get the offset of the consitutents of the first jet in this BX
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
          // tid is the index of the constituent within this BX, which is also the order of the constituent across the whole list of clustered particles in this BX
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
        }  // elements
      }  // blocks
    }  // operator()
  };  // class





  
  ///////////////////////////
  // iterative Seeded Cone Algorithm (greedy, L1T style)
  ///////////////////////////
  // repeat njets times: find highest-pT particle; form jet in cone around it; remove particles from event
  class JetIterKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,  // acc: Alpaka accelerator context (provides thread id, block id, shared memory, atomics, synchronization)
                                  ClusterObjDeviceCollection::View puppi, // current particle buffer (puppi.pt(),...,puppi.cluster() -> store original index)
                                  OffsetsSoA::ConstView bxLookup, // BX particle ranges
                                  BxIndexSoA::ConstView bxIndex, // BX number for each block
                                  float R2,
                                  unsigned int nIters,
                                  ClustersDeviceCollection::View clusters, // particle -> jet assignment storage
                                  uint32_t* tag, //tmp array for prefix scan compaction
                                  ClusterObjDeviceCollection::View work2, // 2nd particle buffer
                                  ClusterObjDeviceCollection::View jets, // tmp jets (each BX write jets starting at 'begin')
                                  OffsetsSoA::View jetBxLookup, // jet counts per BX
                                  BxIndexSoA::View jetBxIndex, // BX numbers for jets
                                  unsigned int* nJetsTotal) const { // global jet counter

      // one thread in grid: make first entry 0 (for prefix scans) -> jetBxLookup.offsets = [0,...]
      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offsets()[0] = 0;

      // for prefix scan (only on GPU)
      // workspace for block prefix scan
      uint32_t* ws = nullptr;
      [[maybe_unused]] constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        // allocate shared memory buffer per block (used by blockPrefixScan)
        ws = alpaka::getDynSharedMem<uint32_t>(acc);
      }
      // determine number of BX
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      // loop over BX (each block processes one BX)
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        // get particle range
        uint32_t begin = bxLookup.offsets()[block_idx];
        uint32_t end = bxLookup.offsets()[block_idx + 1];

        // skip empty BX
        if (end <= begin)
          continue;

        // shared variable size (shrinks each iteration)
        auto& size = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
        size = end - begin;

#ifdef L1TSC_VERBOSE_DEBUG
        if (once_per_block(acc) && (block_idx <= 2))
          printf("In BX %u begin with %u PF candidates: \n", block_idx + 1, end - begin);
#endif
        // store seed particles
        auto& seed_pt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_eta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_phi = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& seed_i = alpaka::declareSharedVar<unsigned int, __COUNTER__>(acc);

        // accumulate jet properties (running sum on multiple threads)
        auto& sum_pt = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_eta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_phi = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& sum_dau = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);

        unsigned int iter = 0;
        // iteration loop (1 jet per iter)
        for (iter = 0; iter < nIters; ++iter) {
          // ping-pong buffer selection (read puppi, write work2 or read work2, write puppi)
          bool even = (iter % 2 == 0);
          auto pt = even ? puppi.pt() : work2.pt();
          auto eta = even ? puppi.eta() : work2.eta();
          auto phi = even ? puppi.phi() : work2.phi();
          auto cluster = even ? puppi.cluster() : work2.cluster();
          auto pt2 = !even ? puppi.pt() : work2.pt();
          auto eta2 = !even ? puppi.eta() : work2.eta();
          auto phi2 = !even ? puppi.phi() : work2.phi();
          auto cluster2 = !even ? puppi.cluster() : work2.cluster();

          // seed search (run by one thread per block)
          if (once_per_block(acc)) {
            // tmp vars (best particle found so far)
            float spt = 0, seta = 0, sphi = 0;
            unsigned int iseed = end;
            // scan all particles in BX from begin to begin+size (only unclustered p.! (at front of current working buffer) not until end)
            for (unsigned int j = begin, myend = begin + size; j < myend; ++j) {
              // highest-pT selection
              if (pt[j] > spt) { //tie-break: earlier best particle
                spt = pt[j];
                seta = eta[j];
                sphi = phi[j];
                iseed = j;
              }
            }
            // reset jet accumulators after choosing seed
            seed_pt = spt;
            seed_eta = seta;
            seed_phi = sphi;
            seed_i = iseed;
            sum_pt = 0;
            sum_eta = 0;
            sum_phi = 0;
            sum_dau = 0;
          }

          // all threads in block wait for seed search to finish and write seed_pt,...,sum_dau
          alpaka::syncBlockThreads(acc);

          // break if no p. left
          if (seed_pt == 0)
            break;

          // loop over current active p. ( [begin, begin + size) )
          for (uint32_t tid : independent_group_elements(acc, size)) {
            auto ipart = tid + begin;  // global index
            // compute distance to seed
            float deta = eta[ipart] - seed_eta;
            float dphi = cms::alpakatools::deltaPhi(acc, phi[ipart], seed_phi);
            float dr2 = deta * deta + dphi * dphi;
            if (dr2 < R2) {
              clusters.is_seed()[cluster[ipart]] = (ipart == seed_i ? 1 : 0); // save whether p. is seed
              clusters.cluster()[cluster[ipart]] = iter; // assign current p. to jet iter
              tag[ipart] = 0; // 0 means clustered now -> remove from active list
              // accumulate jet pt, pt weighted centroid around seed (via eta & phi displacements)...
              alpaka::atomicAdd(acc, &sum_pt, pt[ipart], alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sum_eta, deta * pt[ipart], alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sum_phi, dphi * pt[ipart], alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &sum_dau, 1u, alpaka::hierarchy::Threads{});

            } else {
              tag[ipart] = 1;
            }
          }  // elements

          // ensure all per-particle updates finished (so sum_pt,... done)
          alpaka::syncBlockThreads(acc);

          // one thread per block writes new jet into tmp jet collection
          if (once_per_block(acc)) {
            jets.pt()[begin + iter] = sum_pt;
            jets.eta()[begin + iter] = seed_eta + sum_eta / sum_pt;
            jets.phi()[begin + iter] = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);
            jets.numberOfDaughters()[begin + iter] = sum_dau;
          }

          // prefix scan tag to transform keep flags into destination positions ([0,1,1,0,1] -> [0,1,2,2,3])
          blockPrefixScan(acc, tag + begin, size, ws);

          // loop over p. (1 per thread)
          for (uint32_t tid : independent_group_elements(acc, size)) {
            auto ipart = tid + begin;  // global index
            // for unclustered p.
            if (tag[ipart] > (tid == 0 ? 0 : tag[ipart - 1])) {
              // write unclustered p. to indices begin, begin+1,... in the next working buffer (ping-pong)
              int dest = begin + tag[ipart] - 1;
              pt2[dest] = pt[ipart];
              eta2[dest] = eta[ipart];
              phi2[dest] = phi[ipart];
              cluster2[dest] = cluster[ipart]; // preserve original p. index
            }
          }
          // update number of unclustered p. (one thread per block)
          if (once_per_block(acc)) {
            size = tag[begin + size - 1];
          }

          // ensure all updates finished
          alpaka::syncBlockThreads(acc);

        }  // iter (jet in BX) loop continues

        // one thread per block 
        if (once_per_block(acc)) {
          // jet produced from input BX block_idx belong to same BX number as particles in that block
          jetBxIndex.bx()[block_idx] = bxIndex.bx()[block_idx];
          // store jet count for this BX
          jetBxLookup.offsets()[block_idx + 1] = iter;
          // add this BX's jets to global total
          alpaka::atomicAdd(acc, nJetsTotal, iter, alpaka::hierarchy::Blocks{});
        }

      }  // BX block loop continues
    }  // operator()
  };  // class




  L1TScPhase2SCJetsKernels::L1TScPhase2SCJetsKernels() {}

  std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::run(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float R2,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;

    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // index buffer for reordering
    unsigned int npart = src.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    // buffer for the number of jets per BX
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




  //////////////////////////////////
  // finalize()
  // -> take intermediate results produced by clustering kernel and convert to: jetBxLookup (BX offsets of jets), jets (compact jet collection) & map (jet -> const association map)
  //////////////////////////////////
  std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice>
  L1TScPhase2SCJetsKernels::finalize(Queue& queue,
                                     const PuppiDeviceCollection& src,
                                     const BxLookupDeviceCollection& bxLookup,
                                     const ClustersDeviceCollection& clusters,
                                     const CounterDevice& nJetsTotalDevice,
                                     const ClusterObjDeviceCollection& jetsNonZS,
                                     BxLookupDeviceCollection& jetBxLookup) const {

    // computes number of BX
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    // create launch grid with one block per BX
    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // allocate host-side counter
    auto nJetsTotalHost = CounterHost(queue);
    // copy total jet count from device to host
    alpaka::memcpy(queue, nJetsTotalHost.buffer(), nJetsTotalDevice.buffer());
    // wait for completion
    alpaka::wait(queue);
    // read number of jets into njets
    auto njets = nJetsTotalHost.value();

    // create final jet collection (compact with size njets)
    auto jets = ClusterObjDeviceCollection(njets, queue);

    //// prefix sum to build jet offsets
    // allocate small scratch buffer and zero it
    auto pc = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D{1});
    alpaka::memset(queue, pc, 0x00);
    // create separate launch geometry for prefix scan over jet-offset array (jetBxLookup.offsets())
    uint32_t jet_threads_per_block = 1024;
    uint32_t jet_blocks_per_grid =
        cms::alpakatools::divide_up_by(jetBxLookup.view<OffsetsSoA>().metadata().size(), jet_threads_per_block);
    auto jet_grid = cms::alpakatools::make_workdiv<Acc1D>(jet_blocks_per_grid, jet_threads_per_block);
    // prefix-scan jetBxLookup.offsets(); before counts per BX (e.g. [0,3,2,5]) -> later: cumulative offsets (e.g. [0,3,5,10])
    alpaka::exec<Acc1D>(queue,
                        jet_grid,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        jetBxLookup.view<OffsetsSoA>().offsets().data(),
                        jetBxLookup.view<OffsetsSoA>().offsets().data(),
                        jetBxLookup.view<OffsetsSoA>().metadata().size(),
                        jet_blocks_per_grid,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    //// copy jets to final collection; compute running sum of total clustered particles
    // allocate device-side counter for total number of clustered particles & init to 0
    auto nClusteredDevice = CounterDevice(queue);
    nClusteredDevice.zeroInitialise(queue);
    // launch JetZSKernel (copy jets from jetsNonZS into compact jets; use prefix-scanned jetBxLookup to know where each BX's jets belong;
    // sum numberOfDaughters() over all jets into nClusteredDevice)
    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetZSKernel{},
                        bxLookup.const_view<OffsetsSoA>(),
                        jetsNonZS.const_view(),
                        jetBxLookup.const_view<OffsetsSoA>(),
                        jets.view(),
                        nClusteredDevice.data());

    // copy total clustered-particle count from device to host; wait; read value into nclustered
    auto nClusteredHost = CounterHost(queue);
    alpaka::memcpy(queue, nClusteredHost.buffer(), nClusteredDevice.buffer());
    alpaka::wait(queue);
    unsigned int nclustered = nClusteredHost.value();

    // allocate association map; init to 0 (dim: flat const index list x offset array)
    auto map = AssociationMapDevice({{int(nclustered), int(njets + 1)}}, queue);
    map.zeroInitialise(queue);

    // let's first do the stupid thing and just copy the offsets
    // build launch grids for operations over jet arrays & jet-offset arrays
    uint32_t constit_threads_per_block = 1024;
    uint32_t constit_blocks_per_grid = cms::alpakatools::divide_up_by(njets, constit_threads_per_block);
    auto constit_grid = cms::alpakatools::make_workdiv<Acc1D>(constit_blocks_per_grid, constit_threads_per_block);
    uint32_t constit_blocks_per_grid1 = cms::alpakatools::divide_up_by(njets + 1, constit_threads_per_block);
    auto constit_grid1 = cms::alpakatools::make_workdiv<Acc1D>(constit_blocks_per_grid1, constit_threads_per_block);

    // write jet daughter counts into map.offsets(); prepare offset array for prefix scan ([4,2,5] -> [0,4,2,5])
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

    //// prefix sum for jet constituent offsets
    // reset scan scratch buffer
    alpaka::memset(queue, pc, 0x00);
    // prefix-scan map.offsets() ([0,4,2,5] -> [0,4,6,11])
    alpaka::exec<Acc1D>(queue,
                        constit_grid1,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        map.view<OffsetsSoA>().offsets().data(),
                        map.view<OffsetsSoA>().offsets().data(),
                        njets + 1,
                        constit_blocks_per_grid1,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    // indices of the clustered particles within each jet
    // needs temporary index buffer for reordering, and a key buffer for sorting (lengths npart)
    unsigned int npart = clusters.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    // launch JetToAssociationMapKernel (use clusters to know which jet each part belongs to; use jetBxLookup and map.offsets() to know final structure;
    // fill map.indexes() with actual particle indices for each jet)
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

    // return [jetBxLookup, jets, map]
    // jetBxLookup -> which jet belongs to which BX ([0,2,3,6])
    // jets -> SoA jets.pt(), jets.eta(), jets.phi(), jets.numberOfDaughters(), jets.cluster() (gives BX, redundant with jetBxLookup)
    // map.indexes -> which particles are in jets ([0,1,2,4,5,7,8,9,10])
    // map.offsets -> which particles belong to which jet ([0,2,3,5,6,8,9])
    return std::make_tuple(std::move(jetBxLookup), std::move(jets), std::move(map));
  }






  std::tuple<BxLookupDeviceCollection, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::run(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float R2,
      unsigned int nJets,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    unsigned int npf = src.const_view().metadata().size();

    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // space for output and for reordering inputs
    auto work = ClusterObjDeviceCollection(npf, queue);
    auto work2 = ClusterObjDeviceCollection(npf, queue);

    // a buffer space for counting items
    auto h_tag_device = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npf));
    alpaka::memset(queue, h_tag_device, 0x00);

    // one flat grid per particle, with arbitrary block size
    uint32_t threads_per_flatblock = 1024;
    uint32_t blocks_per_flatgrid = cms::alpakatools::divide_up_by(npf, threads_per_flatblock);
    auto flatgrid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_flatgrid, threads_per_flatblock);

    auto jetsNonZS = ClusterObjDeviceCollection(npf, queue);
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    // buffer for the number of jets per BX
    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    // copy input to work buffer
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
          }
        },
        src.const_view(),
        work.view(),
        clusters.view());

    // iterative clustering
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







  L1TScPhase2SCJetsKernels::return_type
  L1TScPhase2SCJetsKernels::runSeededConeNMSWeighted(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float RSeed2,
      float RCen2,
      float RClu2,
      float alphaSeed,
      float minSeedPt,
      unsigned int nCentroidIters,
      ClustersDeviceCollection& clusters) const
  {
    unsigned int nbx   = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();

    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid   = nbx;
    auto grid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // eta sort buffers
    auto uieta_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto idx_device   = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    alpaka::memset(queue, idx_device, 0x00);

    // seed auxiliary buffers (max size npart; used at indices begin+jseed)
    auto seedWeight_device = alpaka::allocAsyncBuf<float, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto seedPos_device    = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    auto seedOrder_device  = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, cms::alpakatools::Vec1D(npart));
    alpaka::memset(queue, seedWeight_device, 0x00);
    alpaka::memset(queue, seedPos_device, 0x00);
    alpaka::memset(queue, seedOrder_device, 0x00);

    // total jet counter
    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    // work + jets + clusters
    auto work      = ClusterObjDeviceCollection(npart, queue);
    auto jetsNonZS = ClusterObjDeviceCollection(npart, queue);
    jetsNonZS.zeroInitialise(queue);

    clusters.zeroInitialise(queue);

    // per-BX jet lookup (offsets later prefix-scanned in finalize)
    auto jetBxLookup = BxLookupDeviceCollection({{int(nbx), int(nbx + 1)}}, queue);
    jetBxLookup.zeroInitialise(queue);

    float RSeed = std::sqrt(RSeed2);
    float RCen  = std::sqrt(RCen2);
    float RClu  = std::sqrt(RClu2);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernelNMSWeighted{},
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

  L1TScPhase2SCJetsKernels::return_type
  L1TScPhase2SCJetsKernels::runLinkTree(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDeviceCollection& bxLookup,
      float RLink2,
      float ptMin,
      ClustersDeviceCollection& clusters) const
  {
    unsigned int nbx   = bxLookup.const_view<OffsetsSoA>().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();

    uint32_t threads_per_block = 256;
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