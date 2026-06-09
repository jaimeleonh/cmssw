#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernel.h"

#include "HeterogeneousCore/AlpakaInterface/interface/HistoContainer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/radixSort.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T clip(const TAcc& acc, T x, T lo, T hi) {
    return alpaka::math::min(acc, alpaka::math::max(acc, x, lo), hi);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T safe_log(const TAcc& acc, T x) {
    return alpaka::math::log(acc, alpaka::math::max(acc, x, 1.0e-8f));
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T preprocess(
      const TAcc& acc,
      T val,
      T subtract_by,
      T multiply_by,
      T clip_min = -5.0f,
      T clip_max = 5.0f) {
    return clip(acc, (val - subtract_by) * multiply_by, clip_min, clip_max);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T delta_phi(const TAcc& acc, T phi1, T phi2) {
    auto kPi = alpaka::math::constants::pi;
    auto dphi = phi1 - phi2;
    dphi = (dphi > kPi) ? dphi - 2.0f * kPi : dphi;
    dphi = (dphi < -kPi) ? dphi + 2.0f * kPi : dphi;
    return dphi;
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T px_from_pt_phi(const TAcc& acc, T pt, T phi) {
    return pt * alpaka::math::cos(acc, phi);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T py_from_pt_phi(const TAcc& acc, T pt, T phi) {
    return pt * alpaka::math::sin(acc, phi);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T pz_from_pt_eta(const TAcc& acc, T pt, T eta) {
    return pt * alpaka::math::sinh(acc, eta);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T e_from_pt_eta_massless(const TAcc& acc, T pt, T eta) {
    return pt * alpaka::math::cosh(acc, eta);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T pt_from_px_py(const TAcc& acc, T px, T py) {
    return alpaka::math::sqrt(acc, px * px + py * py);
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T eta_from_pt_pz(const TAcc& acc, T pt, T pz) {
    return alpaka::math::asinh(acc, pz / alpaka::math::max(acc, pt, 1.0e-8f));
  }

  template <typename TAcc, typename T>
  ALPAKA_FN_ACC T phi_from_px_py(const TAcc& acc, T px, T py) {
    return alpaka::math::atan2(acc, py, px);
  }

  class ComputeParticleTransformerFeaturesKernel {
  public:
    ALPAKA_FN_ACC void operator()(
        Acc1D const& acc,
        PFCandidateDeviceCollection::ConstView pf,
        AssociationMapDevice::ConstView clusters_cands,
        SoftJetInputDeviceTensor::View input_tensors) const {
      constexpr uint32_t maxParticles = PFFeatures::RowsAtCompileTime;

      for (auto block_idx : independent_groups(acc, clusters_cands.offset().metadata().size() - 1)) {
        auto begin = clusters_cands.offset()[block_idx].offset();
        auto end = clusters_cands.offset()[block_idx + 1].offset();
        auto block_dim = end - begin;
  
        if (false)
          printf("block %i %i\n", block_idx, block_dim);

        if (block_dim == 0)
          continue;

        auto input_tensor = input_tensors[block_idx];
        auto input_size = block_dim > maxParticles ? maxParticles : block_dim;

        auto& shPx = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& shPy = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& shPz = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& shE  = alpaka::declareSharedVar<float, __COUNTER__>(acc);

        auto& shPt  = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& shEta = alpaka::declareSharedVar<float, __COUNTER__>(acc);
        auto& shPhi = alpaka::declareSharedVar<float, __COUNTER__>(acc);

        if (once_per_block(acc)) {
          shPx = 0.0f;
          shPy = 0.0f;
          shPz = 0.0f;
          shE  = 0.0f;

          for (uint32_t tid = 0; tid < input_size; ++tid) {
            auto p = clusters_cands.index()[begin + tid].index(); // global candidate index

            auto pt  = pf.pt()[p];
            auto eta = pf.eta()[p];
            auto phi = pf.phi()[p];

            if (pt <= 0.0f)
              continue;

            auto px = px_from_pt_phi(acc, pt, phi);
            auto py = py_from_pt_phi(acc, pt, phi);
            auto pz = pz_from_pt_eta(acc, pt, eta);
            auto en = e_from_pt_eta_massless(acc, pt, eta);

            shPx += px;
            shPy += py;
            shPz += pz;
            shE  += en;
          }

          shPt  = pt_from_px_py(acc, shPx, shPy);
          shEta = eta_from_pt_pz(acc, shPt, shPz);
          shPhi = phi_from_px_py(acc, shPx, shPy);
        }

        alpaka::syncBlockThreads(acc);

        float etaSign = (shEta >= 0.0f) ? 1.0f : -1.0f;

        for (auto tid : independent_group_elements(acc, input_size)) {
          auto thread_idx = tid + begin; 
          auto p = clusters_cands.index()[thread_idx].index();

          auto pt  = pf.pt()[p];
          auto eta = pf.eta()[p];
          auto phi = pf.phi()[p];

          if (pt <= 0.0f)
            continue;

          auto px = px_from_pt_phi(acc, pt, phi);
          auto py = py_from_pt_phi(acc, pt, phi);
          auto pz = pz_from_pt_eta(acc, pt, eta);
          auto en = e_from_pt_eta_massless(acc, pt, eta);

          auto deta = (eta - shEta) * etaSign;
          auto dphi = delta_phi(acc, phi, shPhi);
          auto dr   = alpaka::math::sqrt(acc, deta * deta + dphi * dphi);

          auto pdgid_abs = alpaka::math::abs(acc, static_cast<int>(pf.pdgid()[p]));

          float isChargedHadron =
              (pdgid_abs == 211 || pdgid_abs == 321 || pdgid_abs == 2212) ? 1.0f : 0.0f;
          float isNeutralHadron =
              (pdgid_abs == 130 || pdgid_abs == 2112) ? 1.0f : 0.0f;
          float isPhoton = (pdgid_abs == 22) ? 1.0f : 0.0f;
          float isElectron = (pdgid_abs == 11) ? 1.0f : 0.0f;
          float isMuon = (pdgid_abs == 13) ? 1.0f : 0.0f;

          if (false)
            printf(
              "P%u %f %f %f %f %f %f %f %f %f %f %d %f %f %f %f %f %f %f %f\n",
              tid,
              pt,
              eta,
              phi,
              px,
              py,
              pz,
              en,
              deta,
              dphi,
              dr,
              pdgid_abs,
              isChargedHadron,
              isNeutralHadron,
              isPhoton,
              isElectron,
              isMuon,
              pf.z0()[p],
              pf.dxy()[p],
              pf.puppiw()[p]
            );

          // pf_points: same order as b_kinadd.yaml / preprocess.json
          input_tensor.points()(tid, 0) = deta;
          input_tensor.points()(tid, 1) = dphi;

          // pf_features, 15 channels:
          // 0 part_pt_log
          // 1 part_e_log
          // 2 part_logptrel
          // 3 part_logerel
          // 4 part_deltaR
          // 5 part_deta
          // 6 part_dphi
          // 7 part_z0
          // 8 part_dxy
          // 9 part_pweight
          // 10 part_isChargedHadron
          // 11 part_isNeutralHadron
          // 12 part_isPhoton
          // 13 part_isElectron
          // 14 part_isMuon
          input_tensor.features()(tid, 0) =
              preprocess(acc, safe_log(acc, pt), 1.7f, 0.7f);
          input_tensor.features()(tid, 1) =
              preprocess(acc, safe_log(acc, en), 2.0f, 0.7f);
          input_tensor.features()(tid, 2) =
              preprocess(acc, safe_log(acc, pt / alpaka::math::max(acc, shPt, 1.0e-8f)), -4.7f, 0.7f);
          input_tensor.features()(tid, 3) =
              preprocess(acc, safe_log(acc, en / alpaka::math::max(acc, shE, 1.0e-8f)), -4.7f, 0.7f);
          input_tensor.features()(tid, 4) =
              preprocess(acc, dr, 0.2f, 4.0f);

          input_tensor.features()(tid, 5) = deta;
          input_tensor.features()(tid, 6) = dphi;

          input_tensor.features()(tid, 7) = pf.z0()[p];

          // If your PFCandidate SoA has dxy/puppiw accessors, use them.
          // Otherwise set these to 0 or update the SoA first.
          input_tensor.features()(tid, 8) =
              pf.dxy()[p];

          input_tensor.features()(tid, 9) =
              pf.puppiw()[p];

          input_tensor.features()(tid, 10) = isChargedHadron;
          input_tensor.features()(tid, 11) = isNeutralHadron;
          input_tensor.features()(tid, 12) = isPhoton;
          input_tensor.features()(tid, 13) = isElectron;
          input_tensor.features()(tid, 14) = isMuon;

          // pf_vectors
          input_tensor.vectors()(tid, 0) = px;
          input_tensor.vectors()(tid, 1) = py;
          input_tensor.vectors()(tid, 2) = pz;
          input_tensor.vectors()(tid, 3) = en;

          // pf_mask
          input_tensor.mask()(tid) = 1.0f;
        }
      }
    }
  };

  SoftJetInputDeviceTensor transform(Queue& queue, 
                 const PFCandidateDeviceCollection& pf,
                 const AssociationMapDevice& association_map) {
    const auto num_clusters = association_map.const_view().offset().metadata().size() - 1;
    // auto input_tensor = SoftJetInputDeviceTensor(32, queue);  // FIXME
    auto input_tensor = SoftJetInputDeviceTensor(queue, num_clusters);
    input_tensor.zeroInitialise(queue);

    // work division
    auto threadsPerBlock = 256;
    auto workDiv = make_workdiv<Acc1D>(num_clusters, threadsPerBlock);

    alpaka::exec<Acc1D>(queue, 
      workDiv, 
      ComputeParticleTransformerFeaturesKernel{}, 
      pf.const_view(),
      association_map.const_view(),
      input_tensor.view());

    // alpaka::exec<Acc1D>(queue,
    //     make_workdiv<Acc1D>(1, 1),
    //     [] ALPAKA_FN_ACC(Acc1D const& acc, SoftJetInputDeviceTensor::View input_tensor) {
    //       if (once_per_grid(acc)) {
    //         for (int c = 0; c < input_tensor.metadata().size(); c++) {
    //           auto jet_cluster = input_tensor[c];
    //           printf("Cluster %d:\n", c);
    //           for (int i = 0; i < JetFeatures::RowsAtCompileTime; i++) {
    //             printf("  PF %d: ", i);
    //             for (int f = 0; f < JetFeatures::ColsAtCompileTime; f++) {
    //               printf("%.2f ", jet_cluster.features()(i, f));
    //             }
    //             printf("\n");
    //           }
    //           printf("Mask: ");
    //           for (int i = 0; i < PaddingMask::RowsAtCompileTime; i++) {
    //             printf("%.1f ", jet_cluster.pad_mask()(i));
    //           }
    //           printf("\n\n");
    //         }
    //       }
    //     },
    //     input_tensor.view());

    return input_tensor;
  }

  // SoftJetInputDeviceTensor transform(Queue& queue, 
  //                const PFCandidateDeviceCollection& pf, 
  //                const BxLookupDeviceCollection& bx_lookup, 
  //                const ClustersDeviceCollection& clusters) {
  //   auto input_tensor = SoftJetInputDeviceTensor(1, queue);
  //   input_tensor.zeroInitialise(queue);
  //   return input_tensor; 
  // }

  // SoftJetInputDeviceTensor transform(Queue& queue, 
  //                const PFCandidateDeviceCollection& pf, 
  //                const ClustersDeviceCollection& clusters) {
  //   auto input_tensor = SoftJetInputDeviceTensor(1, queue);
  //   input_tensor.zeroInitialise(queue);
  //   return input_tensor;     
  // }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
