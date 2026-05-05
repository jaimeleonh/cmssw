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
  using namespace ::l1sc;

  class JetKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  PuppiDeviceCollection::ConstView puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView bxIndex,
                                  float R,
                                  float R2,
                                  uint16_t* uieta,
                                  uint16_t* idx,
                                  ClusterObjDeviceCollection::View work,
                                  ClustersDeviceCollection::View clusters,
                                  ClusterObjDeviceCollection::View jets,
                                  BxLookupSoA::View jetBxLookup,
                                  BxLookupSoA::View jetBxIndex,
                                  unsigned int* nJetsTotal) const {
      constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;
      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offset()[0].offset() = 0;
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uieta[tid + begin] = (puppi.eta()[tid + begin] + 5.f) * (std::numeric_limits<uint16_t>::max() / 10.0f);
          idx[tid + begin] = tid;
        }
      }
      radixSortMulti(acc, uieta, idx, bxLookup.offset().offset().data(), nullptr);

      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
        if (end <= begin)
          continue;
        uint32_t block_dim = end - begin;
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          auto isrc = idx[ipart] + begin;
          work.pt()[ipart] = puppi.pt()[isrc];
          work.eta()[ipart] = puppi.eta()[isrc];
          work.phi()[ipart] = puppi.phi()[isrc];
          work.cluster()[ipart] = isrc;
        }
      }

      auto& nseeds = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
      if (once_per_block(acc))
        nseeds = 0;
      alpaka::syncBlockThreads(acc);

      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
        if (end <= begin)
          continue;

        uint32_t block_dim = end - begin;
        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          uint32_t iseed = tid + begin, icluster = work.cluster()[iseed];
          float seed_pt = work.pt()[iseed], seed_eta = work.eta()[iseed], seed_phi = work.phi()[iseed];
          float sum_pt = seed_pt, sum_eta = 0, sum_phi = 0;
          bool is_seed = true;
          if constexpr (single_thread) {
            if (false)
              continue;
          }
          for (uint32_t j = tid; j > 0; --j) {
            uint32_t ipart = j - 1 + begin;
            float deta = work.eta()[ipart] - seed_eta;
            if (deta < -R)
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < R2) {
              if (work.pt()[ipart] >= seed_pt) {
                is_seed = false;
                break;
              } else {
                if constexpr (single_thread) {
                  /* no is_seed column in 16_1_0_pre3 */
                }
                sum_pt += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }
          for (uint32_t j = tid + 1; j < block_dim; ++j) {
            uint32_t ipart = j + begin;
            float deta = work.eta()[ipart] - seed_eta;
            if (deta > R)
              break;
            float dphi = cms::alpakatools::deltaPhi(acc, work.phi()[ipart], seed_phi);
            if (deta * deta + dphi * dphi < R2) {
              if (work.pt()[ipart] > seed_pt) {
                is_seed = false;
                break;
              } else {
                if constexpr (single_thread) {
                  /* no is_seed column in 16_1_0_pre3 */
                }
                sum_pt += work.pt()[ipart];
                sum_eta += work.pt()[ipart] * deta;
                sum_phi += work.pt()[ipart] * dphi;
              }
            }
          }
          sum_eta = seed_eta + sum_eta / sum_pt;
          sum_phi = cms::alpakatools::reducePhiRange(acc, seed_phi + sum_phi / sum_pt);
          if (is_seed) {
            auto ijet = alpaka::atomicAdd(acc, &nseeds, 1u, alpaka::hierarchy::Threads{}) + begin;
            jets.pt()[ijet] = sum_pt;
            jets.eta()[ijet] = sum_eta;
            jets.phi()[ijet] = sum_phi;
            jets.cluster()[ijet] = ijet - begin;
            jets.numberOfDaughters()[ijet] = 0;
            /* no is_seed column in 16_1_0_pre3 */
          }
        }
        alpaka::syncBlockThreads(acc);

        for (uint32_t tid : independent_group_elements(acc, block_dim)) {
          auto ipart = tid + begin;
          uint32_t icluster = work.cluster()[ipart];
          float nearest = R2;
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
          if constexpr (single_thread) {
          }
        }

        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx].bx() = bxIndex.bx()[block_idx].bx();
          jetBxLookup.offset()[block_idx + 1].offset() = nseeds;
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
                                  BxLookupSoA::ConstView bxLookup,
                                  ClusterObjDeviceCollection::ConstView jetsNonZS,
                                  BxLookupSoA::ConstView jetBxLookup,
                                  ClusterObjDeviceCollection::View jets,
                                  unsigned int* nClustered) const {
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t beginSrc = bxLookup.offset()[block_idx].offset();
        uint32_t beginDst = jetBxLookup.offset()[block_idx].offset();
        uint32_t endDst = jetBxLookup.offset()[block_idx + 1].offset();
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
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
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
        }
      }
      alpaka::syncBlockThreads(acc);
      radixSortMulti(acc, key, idx, bxLookup.offset().offset().data(), nullptr);
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
        if (end <= begin)
          continue;
        uint32_t jetsBegin = jetBxLookup.offset()[block_idx].offset();
        uint32_t jetsEnd = jetBxLookup.offset()[block_idx + 1].offset();
        assert(jetsEnd >= jetsBegin);
        assert(jetsEnd <= njets);
        assert(jetBxLookup.offset()[grid_dim].offset() == njets);
        assert(clusterdParticleOffsets.offset().metadata().size() == int(njets + 1));
        assert(clusterdParticleOffsets.offset()[njets].offset() == nclustered);
        uint32_t clusteredBegin = clusterdParticleOffsets.offset()[jetsBegin].offset();
        uint32_t clusteredEnd = clusterdParticleOffsets.offset()[jetsEnd].offset();
        assert(clusteredEnd >= clusteredBegin);
        assert((clusteredEnd - clusteredBegin) <= (end - begin));
        for (uint32_t tid : independent_group_elements(acc, clusteredEnd - clusteredBegin)) {
          uint32_t i = tid + begin;
          uint32_t icand = idx[i] + begin;
          map.index()[clusteredBegin + tid].index() = icand;
          assert(clusteredBegin + tid < nclustered);
          assert(map.index().metadata().size() == int(nclustered));
        }
      }
    }
  };

  class JetIterKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                  ClusterObjDeviceCollection::View puppi,
                                  BxLookupSoA::ConstView bxLookup,
                                  BxLookupSoA::ConstView bxIndex,
                                  float R2,
                                  unsigned int nIters,
                                  ClustersDeviceCollection::View clusters,
                                  uint32_t* tag,
                                  ClusterObjDeviceCollection::View work2,
                                  ClusterObjDeviceCollection::View jets,
                                  BxLookupSoA::View jetBxLookup,
                                  BxLookupSoA::View jetBxIndex,
                                  unsigned int* nJetsTotal) const {
      if (cms::alpakatools::once_per_grid(acc))
        jetBxLookup.offset()[0].offset() = 0;
      uint32_t* ws = nullptr;
      [[maybe_unused]] constexpr bool single_thread = requires_single_thread_per_block<TAcc>::value;
      if constexpr (!requires_single_thread_per_block_v<TAcc>) {
        ws = alpaka::getDynSharedMem<uint32_t>(acc);
      }
      uint32_t grid_dim = alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];
      for (uint32_t block_idx : independent_groups(acc, grid_dim)) {
        uint32_t begin = bxLookup.offset()[block_idx].offset();
        uint32_t end = bxLookup.offset()[block_idx + 1].offset();
        if (end <= begin)
          continue;

        auto& size = alpaka::declareSharedVar<uint32_t, __COUNTER__>(acc);
        size = end - begin;
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
              /* no is_seed column in 16_1_0_pre3 */
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
        }
        if (once_per_block(acc)) {
          jetBxIndex.bx()[block_idx].bx() = bxIndex.bx()[block_idx].bx();
          jetBxLookup.offset()[block_idx + 1].offset() = iter;
          alpaka::atomicAdd(acc, nJetsTotal, iter, alpaka::hierarchy::Blocks{});
        }
      }
    }
  };

  L1TScPhase2SCJetsKernels::L1TScPhase2SCJetsKernels() {}

  std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::run(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDevice& bxLookup,
      float R2,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;

    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    unsigned int npart = src.const_view().metadata().size();
    auto h_key_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    auto h_idx_device = alpaka::allocAsyncBuf<uint16_t, Idx>(queue, Vec1D(npart));
    alpaka::memset(queue, h_idx_device, 0x00);

    auto nJetsTotalDevice = CounterDevice(queue);
    nJetsTotalDevice.zeroInitialise(queue);

    auto work = ClusterObjDeviceCollection(queue, int(npart));
    auto jetsNonZS = ClusterObjDeviceCollection(queue, int(npart));
    clusters.zeroInitialise(queue);
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDevice(queue, int(nbx), int(nbx + 1));
    jetBxLookup.zeroInitialise(queue);

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetKernel{},
                        src.const_view(),
                        bxLookup.const_view(),
                        bxLookup.const_view(),
                        std::sqrt(R2),
                        R2,
                        h_key_device.data(),
                        h_idx_device.data(),
                        work.view(),
                        clusters.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view(),
                        jetBxLookup.view(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

  std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::finalize(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDevice& bxLookup,
      const ClustersDeviceCollection& clusters,
      const CounterDevice& nJetsTotalDevice,
      const ClusterObjDeviceCollection& jetsNonZS,
      BxLookupDevice& jetBxLookup) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto nJetsTotalHost = CounterHost(queue);
    alpaka::memcpy(queue, nJetsTotalHost.buffer(), nJetsTotalDevice.buffer());
    alpaka::wait(queue);
    auto njets = nJetsTotalHost.value();

    auto jets = ClusterObjDeviceCollection(queue, int(njets));

    auto pc = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D{1});
    alpaka::memset(queue, pc, 0x00);
    uint32_t jet_threads_per_block = 1024;
    uint32_t jet_blocks_per_grid =
        cms::alpakatools::divide_up_by(jetBxLookup.view().offset().metadata().size(), jet_threads_per_block);
    auto jet_grid = cms::alpakatools::make_workdiv<Acc1D>(jet_blocks_per_grid, jet_threads_per_block);
    alpaka::exec<Acc1D>(queue,
                        jet_grid,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        jetBxLookup.view().offset().offset().data(),
                        jetBxLookup.view().offset().offset().data(),
                        jetBxLookup.view().offset().metadata().size(),
                        jet_blocks_per_grid,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));

    auto nClusteredDevice = CounterDevice(queue);
    nClusteredDevice.zeroInitialise(queue);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetZSKernel{},
                        bxLookup.const_view(),
                        jetsNonZS.const_view(),
                        jetBxLookup.const_view(),
                        jets.view(),
                        nClusteredDevice.data());

    auto nClusteredHost = CounterHost(queue);
    alpaka::memcpy(queue, nClusteredHost.buffer(), nClusteredDevice.buffer());
    alpaka::wait(queue);
    unsigned int nclustered = nClusteredHost.value();

    auto map = AssociationMapDevice(queue, int(nclustered), int(njets + 1));
    map.zeroInitialise(queue);

    uint32_t constit_threads_per_block = 1024;
    uint32_t constit_blocks_per_grid = cms::alpakatools::divide_up_by(njets, constit_threads_per_block);
    auto constit_grid = cms::alpakatools::make_workdiv<Acc1D>(constit_blocks_per_grid, constit_threads_per_block);
    uint32_t constit_blocks_per_grid1 = cms::alpakatools::divide_up_by(njets + 1, constit_threads_per_block);
    auto constit_grid1 = cms::alpakatools::make_workdiv<Acc1D>(constit_blocks_per_grid1, constit_threads_per_block);
    alpaka::exec<Acc1D>(
        queue,
        constit_grid,
        [] ALPAKA_FN_ACC(Acc1D const& acc, ClusterObjDeviceCollection::ConstView jets, AssociationMapSoA::View offsets) {
          if (cms::alpakatools::once_per_grid(acc))
            offsets.offset()[0].offset() = 0;
          for (int32_t idx : cms::alpakatools::uniform_elements(acc, offsets.offset().metadata().size() - 1)) {
            offsets.offset()[idx + 1].offset() = jets.numberOfDaughters()[idx];
          }
        },
        jets.const_view(),
        map.view());

    alpaka::memset(queue, pc, 0x00);
    alpaka::exec<Acc1D>(queue,
                        constit_grid1,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        map.view().offset().offset().data(),
                        map.view().offset().offset().data(),
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
                        bxLookup.const_view(),
                        jetBxLookup.const_view(),
                        map.const_view(),
                        clusters.const_view(),
                        njets,
                        nclustered,
                        h_key_device.data(),
                        h_idx_device.data(),
                        map.view());

    return std::make_tuple(std::move(jetBxLookup), std::move(jets), std::move(map));
  }

  std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice> L1TScPhase2SCJetsKernels::run(
      Queue& queue,
      const PuppiDeviceCollection& src,
      const BxLookupDevice& bxLookup,
      float R2,
      unsigned int nJets,
      ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    unsigned int npf = src.const_view().metadata().size();

    uint32_t threads_per_block = 256;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto work = ClusterObjDeviceCollection(queue, int(npf));
    auto work2 = ClusterObjDeviceCollection(queue, int(npf));

    auto h_tag_device = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npf));
    alpaka::memset(queue, h_tag_device, 0x00);

    uint32_t threads_per_flatblock = 1024;
    uint32_t blocks_per_flatgrid = cms::alpakatools::divide_up_by(npf, threads_per_flatblock);
    auto flatgrid = cms::alpakatools::make_workdiv<Acc1D>(blocks_per_flatgrid, threads_per_flatblock);

    auto jetsNonZS = ClusterObjDeviceCollection(queue, int(npf));
    jetsNonZS.zeroInitialise(queue);

    auto jetBxLookup = BxLookupDevice(queue, int(nbx), int(nbx + 1));
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
          }
        },
        src.const_view(),
        work.view(),
        clusters.view());

    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetIterKernel{},
                        work.view(),
                        bxLookup.const_view(),
                        bxLookup.const_view(),
                        R2,
                        nJets,
                        clusters.view(),
                        h_tag_device.data(),
                        work2.view(),
                        jetsNonZS.view(),
                        jetBxLookup.view(),
                        jetBxLookup.view(),
                        nJetsTotalDevice.data());

    return finalize(queue, src, bxLookup, clusters, nJetsTotalDevice, jetsNonZS, jetBxLookup);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
