#include "L1TriggerScouting/Phase2/plugins/alpaka/L1TScPhase2SCJetsKernels.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <limits>
#include <optional>
#include <type_traits>
#include "HeterogeneousCore/AlpakaInterface/interface/memory.h"
#include "HeterogeneousCore/AlpakaInterface/interface/host.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "HeterogeneousCore/AlpakaInterface/interface/prefixScan.h"
#include "HeterogeneousCore/AlpakaMath/interface/deltaPhi.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/radixSort64.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/SeededCone.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/SCNMSW.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/LinkTree.h"
#include "L1TriggerScouting/Phase2/plugins/alpaka/JetOutput.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {
  using namespace cms::alpakatools;
  constexpr uint32_t kThreadsPerBlock = 64;
  std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice>
  L1TScPhase2SCJetsKernels::finalize(Queue& queue,
                                     const PuppiDeviceCollection& src,
                                     const BxLookupDevice& bxLookup,
                                     ClustersDeviceCollection& clusters,
                                     CounterDevice& nJetsTotalDevice,
                                     ClusterObjDeviceCollection& jetsNonZS,
                                     BxLookupDevice& jetBxLookup) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    // Preserve the target branch's optional constituent cut before compaction.
    if (minConstituents_ > 1) {
      nJetsTotalDevice.zeroInitialise(queue);
      auto remap = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D(src.const_view().metadata().size()));
      alpaka::exec<Acc1D>(queue,
                          grid,
                          JetMinConstituentsFilterKernel{},
                          bxLookup.const_view(),
                          jetsNonZS.view(),
                          clusters.view(),
                          jetBxLookup.view(),
                          minConstituents_,
                          remap.data(),
                          nJetsTotalDevice.data());
    }

    auto nJetsTotalHost = CounterHost(queue);
    alpaka::memcpy(queue, nJetsTotalHost.buffer(), nJetsTotalDevice.buffer());
    alpaka::wait(queue);
    auto njets = nJetsTotalHost.value();

    auto jets = ClusterObjDeviceCollection(queue, njets);

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

    // CUDA rejects a zero-block launch for the constituent-count kernel.
    // Keep the final offset sentinel even when the orbit contains no jets.
    if (njets == 0) {
      auto map = AssociationMapDevice(queue, 0, 1);
      map.zeroInitialise(queue);
      return std::make_tuple(std::move(jetBxLookup), std::move(jets), std::move(map));
    }

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
            offsets.offset().offset()[0] = 0;
          for (int32_t idx : cms::alpakatools::uniform_elements(acc, offsets.offset().metadata().size() - 1)) {
            offsets.offset().offset()[idx + 1] = jets.numberOfDaughters()[idx];
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
  std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice>
  L1TScPhase2SCJetsKernels::runSeededCone(Queue& queue,
                                          const PuppiDeviceCollection& src,
                                          const BxLookupDevice& bxLookup,
                                          float R2,
                                          unsigned int nJets,
                                          ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();
    auto grid = make_workdiv<Acc1D>(nbx, kThreadsPerBlock);
    auto order = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npart));
    auto jets = ClusterObjDeviceCollection(queue, npart);
    auto jetBxLookup = BxLookupDevice(queue, int(nbx), int(nbx + 1));
    jetBxLookup.zeroInitialise(queue);
    auto count = CounterDevice(queue);
    count.zeroInitialise(queue);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        SeededConeKernel{},
                        src.const_view(),
                        bxLookup.const_view(),
                        bxLookup.const_view(),
                        R2,
                        nJets,
                        clusters.view(),
                        order.data(),
                        jets.view(),
                        jetBxLookup.view(),
                        jetBxLookup.view(),
                        count.data());
    return finalize(queue, src, bxLookup, clusters, count, jets, jetBxLookup);
  }
  std::tuple<BxLookupDevice, ClusterObjDeviceCollection, AssociationMapDevice>
  L1TScPhase2SCJetsKernels::finalizeSCNMSW(Queue& queue,
                                           const PuppiDeviceCollection& src,
                                           const BxLookupDevice& bxLookup,
                                           const ClustersDeviceCollection& clusters,
                                           const PortableObject<SCNMSWCounts>& countsDevice,
                                           const ClusterObjDeviceCollection& jetsNonZS,
                                           BxLookupDevice& jetBxLookup,
                                           uint32_t* scratch,
                                           uint32_t* particleBxOffsets) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    uint32_t threads_per_block = kThreadsPerBlock;
    uint32_t blocks_per_grid = nbx;
    auto grid = make_workdiv<Acc1D>(blocks_per_grid, threads_per_block);

    auto countsHost = PortableHostObject<SCNMSWCounts>(queue);
    alpaka::memcpy(queue, countsHost.buffer(), countsDevice.buffer());
    alpaka::wait(queue);
    auto njets = countsHost.value().jets;
    unsigned int nclustered = countsHost.value().clustered;

    auto jets = ClusterObjDeviceCollection(queue, njets);

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

    auto map = AssociationMapDevice(queue, int(nclustered), int(njets + 1));

    alpaka::memset(queue, pc, 0x00);
    alpaka::exec<Acc1D>(queue,
                        jet_grid,
                        cms::alpakatools::multiBlockPrefixScan<uint32_t>{},
                        particleBxOffsets,
                        particleBxOffsets,
                        nbx + 1,
                        jet_blocks_per_grid,
                        pc.data(),
                        alpaka::getPreferredWarpSize(alpaka::getDev(queue)));
    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetCompactOffsetsKernel{},
                        bxLookup.const_view(),
                        jetsNonZS.const_view(),
                        jetBxLookup.const_view(),
                        particleBxOffsets,
                        jets.view(),
                        map.view());
    alpaka::exec<Acc1D>(queue,
                        grid,
                        JetAssociationRankKernel<4>{},
                        src.const_view(),
                        bxLookup.const_view(),
                        jetBxLookup.const_view(),
                        map.const_view(),
                        clusters.const_view(),
                        scratch,
                        map.view());

    return std::make_tuple(std::move(jetBxLookup), std::move(jets), std::move(map));
  }

  L1TScPhase2SCJetsKernels::return_type L1TScPhase2SCJetsKernels::runSCNMSW(Queue& queue,
                                                                            const PuppiDeviceCollection& src,
                                                                            const BxLookupDevice& bxLookup,
                                                                            float RSeed2,
                                                                            float RCen2,
                                                                            float RClu2,
                                                                            float alpha,
                                                                            int maxRounds,
                                                                            ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();
    auto grid = make_workdiv<Acc1D>(nbx, kThreadsPerBlock);
    using IndexBuffer = decltype(alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npart)));
    std::optional<IndexBuffer> indices;
    std::optional<ClusterObjDeviceCollection> work;
    if (maxRounds == 1)
      indices.emplace(alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(npart)));
    else
      work.emplace(queue, npart);
    uint32_t* scratch = maxRounds == 1 ? indices->data() : work->view().numberOfDaughters().data();
    auto jets = ClusterObjDeviceCollection(queue, npart);
    auto jetBxLookup = BxLookupDevice(queue, int(nbx), int(nbx + 1));
    auto particleBxOffsets = alpaka::allocAsyncBuf<uint32_t, Idx>(queue, Vec1D(nbx + 1));
    auto counts = PortableObject<SCNMSWCounts>(queue);
    counts.zeroInitialise(queue);
    if (maxRounds == 1) {
      auto launch = [&](auto kernel) {
        alpaka::exec<Acc1D>(queue,
                            grid,
                            kernel,
                            src.const_view(),
                            bxLookup.const_view(),
                            bxLookup.const_view(),
                            RSeed2,
                            RCen2,
                            RClu2,
                            std::sqrt(std::max(RSeed2, RCen2)),
                            alpha,
                            scratch,
                            clusters.view(),
                            jets.view(),
                            jetBxLookup.view(),
                            jetBxLookup.view(),
                            &counts.data()->jets,
                            &counts.data()->clustered,
                            particleBxOffsets.data());
      };
      if (alpha == 0.f)
        launch(SCNMSWOneRoundKernel<true>{});
      else
        launch(SCNMSWOneRoundKernel<false>{});
    } else {
      alpaka::exec<Acc1D>(queue,
                          grid,
                          SCNMSWMultiRoundKernel{},
                          src.const_view(),
                          bxLookup.const_view(),
                          bxLookup.const_view(),
                          RSeed2,
                          RCen2,
                          RClu2,
                          std::sqrt(std::max(RSeed2, RCen2)),
                          alpha,
                          maxRounds,
                          scratch,
                          work->view(),
                          clusters.view(),
                          jets.view(),
                          jetBxLookup.view(),
                          jetBxLookup.view(),
                          &counts.data()->jets,
                          &counts.data()->clustered,
                          particleBxOffsets.data());
    }
    if (minConstituents_ > 1) {
      // Filtering changes both jet and constituent counts; recompute them through
      // the shared finalizer rather than using the one-round precomputed counts.
      auto filteredCount = CounterDevice(queue);
      return finalize(queue, src, bxLookup, clusters, filteredCount, jets, jetBxLookup);
    }
    return finalizeSCNMSW(queue, src, bxLookup, clusters, counts, jets, jetBxLookup, scratch, particleBxOffsets.data());
  }
  L1TScPhase2SCJetsKernels::return_type L1TScPhase2SCJetsKernels::runLinkTree(Queue& queue,
                                                                              const PuppiDeviceCollection& src,
                                                                              const BxLookupDevice& bxLookup,
                                                                              float RLink2,
                                                                              float ptMin,
                                                                              ClustersDeviceCollection& clusters) const {
    unsigned int nbx = bxLookup.const_view().offset().metadata().size() - 1;
    unsigned int npart = src.const_view().metadata().size();
    auto grid = make_workdiv<Acc1D>(nbx, kThreadsPerBlock);
    auto parent = alpaka::allocAsyncBuf<int32_t, Idx>(queue, Vec1D(npart));
    auto jets = ClusterObjDeviceCollection(queue, npart);
    auto jetBxLookup = BxLookupDevice(queue, int(nbx), int(nbx + 1));
    jetBxLookup.zeroInitialise(queue);
    auto count = CounterDevice(queue);
    count.zeroInitialise(queue);
    alpaka::exec<Acc1D>(queue,
                        grid,
                        LinkTreeKernel{},
                        src.const_view(),
                        bxLookup.const_view(),
                        bxLookup.const_view(),
                        RLink2,
                        ptMin,
                        parent.data(),
                        clusters.view(),
                        jets.view(),
                        jetBxLookup.view(),
                        jetBxLookup.view(),
                        count.data());
    return finalize(queue, src, bxLookup, clusters, count, jets, jetBxLookup);
  }
}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
