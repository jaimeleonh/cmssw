#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernelTorchScript.h"

#include <cstdint>
#include <cstdio>

#include "DataFormats/L1ScoutingSoA/interface/ClustersSoA.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetTensorSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels::torchscript {
  namespace {
    using namespace cms::alpakatools;

    template <typename TAcc>
    ALPAKA_FN_ACC float px(TAcc const& acc, float pt, float phi) {
      return pt * alpaka::math::cos(acc, phi);
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float py(TAcc const& acc, float pt, float phi) {
      return pt * alpaka::math::sin(acc, phi);
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float pz(TAcc const& acc, float pt, float eta) {
      return pt * alpaka::math::sinh(acc, eta);
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float energy(TAcc const& acc, float pt, float eta) {
      return pt * alpaka::math::cosh(acc, eta);
    }

    // log(clamp(x, min=eps)), matching the reference preprocessing (log_clamped).
    template <typename TAcc>
    ALPAKA_FN_ACC float logClamped(TAcc const& acc, float x, float eps = 1.e-6f) {
      return alpaka::math::log(acc, alpaka::math::max(acc, x, eps));
    }

    // Wrapped phi difference, matching the reference preprocessing (atan2(sin, cos)).
    template <typename TAcc>
    ALPAKA_FN_ACC float deltaPhi(TAcc const& acc, float phi1, float phi2) {
      float const d = phi1 - phi2;
      return alpaka::math::atan2(acc, alpaka::math::sin(acc, d), alpaka::math::cos(acc, d));
    }

    inline constexpr uint32_t kInputThreads = 32;
    inline constexpr uint32_t kMaxParticlesPerJet = ::l1sc::SoftJetNewPFFeatures::RowsAtCompileTime;

    ALPAKA_FN_ACC inline uint32_t bunchCrossingOf(::l1sc::BxLookupSoA::ConstView jetBxLookup,
                                                  uint32_t numberOfBxs,
                                                  uint32_t jet) {
      uint32_t low = 0;
      uint32_t high = numberOfBxs;
      while (high - low > 1) {
        uint32_t const middle = low + (high - low) / 2;
        if (static_cast<uint32_t>(jetBxLookup.offset()[middle].offset()) <= jet)
          low = middle;
        else
          high = middle;
      }
      return low;
    }

    class ComputeTorchScriptInputsKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    l1sc::PFCandidateDeviceCollection::ConstView pf,
                                    l1sc::ClusterObjDeviceCollection::ConstView jets,
                                    l1sc::AssociationMapDevice::ConstView associationMap,
                                    ::l1sc::BxLookupSoA::ConstView jetBxLookup,
                                    l1sc::VertexDeviceCollection::ConstView vertices,
                                    ::l1sc::BxLookupSoA::ConstView vertexBxLookup,
                                    l1sc::SoftJetNewInputDeviceTensor::View inputs,
                                    uint32_t numberOfJets,
                                    uint32_t numberOfBxs) const {
        auto& cachePt = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePx = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePy = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePz = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cacheEnergy = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& sortIdx = alpaka::declareSharedVar<uint32_t[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& jet4 = alpaka::declareSharedVar<float[7], __COUNTER__>(acc);

        auto const lane = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const lanes = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const block = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        for (uint32_t jet = block; jet < numberOfJets; jet += blocks) {
          auto const bx = bunchCrossingOf(jetBxLookup, numberOfBxs, jet);
          auto const vertexBegin = vertexBxLookup.offset()[bx].offset();
          auto const vertexEnd = vertexBxLookup.offset()[bx + 1].offset();
          if (vertexEnd <= vertexBegin)
            continue;
          float const primaryVertexZ = vertices.z0()[vertexBegin];

          auto const begin = associationMap.offset()[jet].offset();
          auto const end = associationMap.offset()[jet + 1].offset();
          auto const available = static_cast<uint32_t>(end - begin);
          auto const count = available < kMaxParticlesPerJet ? available : kMaxParticlesPerJet;
          if (count == 0)
            continue;

          for (uint32_t i = lane; i < count; i += lanes) {
            auto const candidate = associationMap.index()[begin + i].index();
            auto const candidatePt = pf.pt()[candidate];
            auto const candidateEta = pf.eta()[candidate];
            auto const candidatePhi = pf.phi()[candidate];
            bool const usable = candidatePt > 0.f;
            cachePt[i] = candidatePt;
            cachePx[i] = usable ? px(acc, candidatePt, candidatePhi) : 0.f;
            cachePy[i] = usable ? py(acc, candidatePt, candidatePhi) : 0.f;
            cachePz[i] = usable ? pz(acc, candidatePt, candidateEta) : 0.f;
            cacheEnergy[i] = usable ? energy(acc, candidatePt, candidateEta) : 0.f;
          }
          alpaka::syncBlockThreads(acc);

          if (lane == 0) {
            float sumX = 0.f, sumY = 0.f, sumZ = 0.f, sumE = 0.f;
            for (uint32_t i = 0; i < count; ++i) {
              sumX += cachePx[i];
              sumY += cachePy[i];
              sumZ += cachePz[i];
              sumE += cacheEnergy[i];
            }
            float const jetPt = alpaka::math::sqrt(acc, sumX * sumX + sumY * sumY);
            jet4[0] = sumX;
            jet4[1] = sumY;
            jet4[2] = sumZ;
            jet4[3] = sumE;
            jet4[4] = jetPt;
            jet4[5] = alpaka::math::asinh(acc, sumZ / alpaka::math::max(acc, jetPt, 1.e-8f));
            jet4[6] = alpaka::math::atan2(acc, sumY, sumX);
          }
          alpaka::syncBlockThreads(acc);

          // The model was trained on the clustered jet's own pt/eta/phi, so use those
          // directly from the ClusterObj collection. The jet SoA has no mass column, so
          // the mass is recomputed from the constituent 4-vector sum accumulated above.
          float const jetEnergy = jet4[3];
          float const jetMass2 = jetEnergy * jetEnergy - (jet4[0] * jet4[0] + jet4[1] * jet4[1] + jet4[2] * jet4[2]);
          float const jetMass = alpaka::math::sqrt(acc, alpaka::math::max(acc, jetMass2, 0.f));

          float const jetPt = jets.pt()[jet];
          float const jetEta = jets.eta()[jet];
          float const jetPhi = jets.phi()[jet];

          auto output = inputs[jet];

          // Per-jet features consumed by the model:
          //   (log(1+jet_pt), jet_eta, log(1+jet_mass), log1p(retained_nconst)).
          // The retained-constituent count is filled after the sort below.
          if (lane == 0) {
            output.jet_features()(0) = logClamped(acc, jetPt + 1.f);
            output.jet_features()(1) = jetEta;
            // output.jet_features()(2) = logClamped(acc, jetMass + 1.f);
            output.jet_features()(2) = 0.f; // FIXME after comparison
          }

          // Sort constituent indices by decreasing pt (bitonic sort on shared memory).
          // Padded entries (pt <= 0) sink to the end and stay zero-initialised.
          for (uint32_t i = lane; i < kMaxParticlesPerJet; i += lanes)
            sortIdx[i] = i;
          alpaka::syncBlockThreads(acc);

          for (uint32_t k = 2; k <= kMaxParticlesPerJet; k <<= 1) {
            for (uint32_t j = k >> 1; j > 0; j >>= 1) {
              for (uint32_t i = lane; i < kMaxParticlesPerJet; i += lanes) {
                uint32_t const partner = i ^ j;
                if (partner > i) {
                  bool const up = (i & k) == 0;
                  // Guard on the permuted index, not the position: after swaps, a
                  // position i < count can hold a padded slot and a position i >= count
                  // can hold a real constituent. Only indices < count have a cached pt.
                  uint32_t const idxI = sortIdx[i];
                  uint32_t const idxP = sortIdx[partner];
                  float const ptI = idxI < count ? cachePt[idxI] : 0.f;
                  float const ptP = idxP < count ? cachePt[idxP] : 0.f;
                  // "up" sorts descending (larger pt first).
                  if ((up && ptI < ptP) || (!up && ptI > ptP)) {
                    uint32_t const tmp = sortIdx[i];
                    sortIdx[i] = sortIdx[partner];
                    sortIdx[partner] = tmp;
                  }
                }
              }
              alpaka::syncBlockThreads(acc);
            }
          }
          alpaka::syncBlockThreads(acc);

          // log1p of the number of retained (pt-sorted, capped) constituents.
          if (lane == 0)
            output.jet_features()(3) = alpaka::math::log(acc, 1.f + static_cast<float>(count));

          for (uint32_t i = lane; i < count; i += lanes) {
            auto const src = sortIdx[i];
            if (cachePt[src] <= 0.f)
              continue;
            auto const candidate = associationMap.index()[begin + src].index();

            float const candPt = cachePt[src];
            float const candEta = pf.eta()[candidate];
            float const candPhi = pf.phi()[candidate];
            float const deta = candEta - jetEta;
            float const dphi = deltaPhi(acc, candPhi, jetPhi);

            // Prepared points: (pt, deta, dphi) in local jet coordinates.
            output.points()(i, 0) = candPt;
            output.points()(i, 1) = deta;
            output.points()(i, 2) = dphi;

            // Prepared features (14), matching the reference preprocessing column order:
            //   log(1+pt), log(pt/jet_pt), deta, sin(dphi), cos(dphi),
            //   quality, dxy, z0, puppiw,
            //   pid_charged_hadron, pid_neutral_hadron, pid_photon, pid_electron, pid_muon
            output.features()(i, 0) = logClamped(acc, candPt + 1.f);
            output.features()(i, 1) = logClamped(acc, candPt / alpaka::math::max(acc, jetPt, 1.e-6f));
            output.features()(i, 2) = deta;
            output.features()(i, 3) = alpaka::math::sin(acc, dphi);
            output.features()(i, 4) = alpaka::math::cos(acc, dphi);
            output.features()(i, 5) = 0.f;
            // output.features()(i, 5) = static_cast<float>(pf.quality()[candidate]);
            output.features()(i, 6) = pf.dxy()[candidate];
            output.features()(i, 7) = pf.z0()[candidate];
            output.features()(i, 8) = pf.puppiw()[candidate];

            int const absid = pf.pdgid()[candidate] < 0 ? -pf.pdgid()[candidate] : pf.pdgid()[candidate];
            output.features()(i, 9) = (absid == 211 || absid == 321 || absid == 2212) ? 1.f : 0.f;
            output.features()(i, 10) =
                (absid == 111 || absid == 130 || absid == 310 || absid == 311 || absid == 2112) ? 1.f : 0.f;
            output.features()(i, 11) = (absid == 22) ? 1.f : 0.f;
            output.features()(i, 12) = (absid == 11) ? 1.f : 0.f;
            output.features()(i, 13) = (absid == 13) ? 1.f : 0.f;

            output.mask()(i) = 1.f;
          }
          alpaka::syncBlockThreads(acc);

          // Debug: dump the first jet's preprocessed inputs, mirroring the structure of
          // inference.preprocess: points (pt, deta, dphi), features (14), mask, jet_features (4).
          // if (jet == 0 && lane == 0) {
          //   printf("=== preprocess inputs for jet 0 (count=%u) ===\n", count);
          //   printf("jet_features [log(1+pt), eta, log(1+mass), log1p(nconst)]: [%f, %f, %f, %f]\n",
          //          output.jet_features()(0),
          //          output.jet_features()(1),
          //          output.jet_features()(2),
          //          output.jet_features()(3));
          //   for (uint32_t i = 0; i < kMaxParticlesPerJet; ++i) {
          //     printf("cand %2u mask=%.0f points(pt,deta,dphi)=[%.6f, %.6f, %.6f] features=[",
          //            i,
          //            output.mask()(i),
          //            output.points()(i, 0),
          //            output.points()(i, 1),
          //            output.points()(i, 2));
          //     for (uint32_t j = 0; j < 14; ++j)
          //       printf("%s%.6f", (j ? ", " : ""), output.features()(i, j));
          //     printf("]\n");
          //   }
          // }
        }
      }
    };
  }  // namespace

  SoftJetNewInputDeviceTensor transform(
      Queue& queue,
      const PFCandidateDeviceCollection& pf,
      const ClusterObjDeviceCollection& jets,
      const AssociationMapDevice& association_map,
      const BxLookupDevice& jetBxLookup,
      const VertexDeviceCollection& vertices,
      const BxLookupDevice& vertexBxLookup) {
    auto const numberOfOffsets = association_map.const_view().offset().metadata().size();
    auto const numberOfJets = numberOfOffsets == 0 ? 0u : numberOfOffsets - 1u;

    SoftJetNewInputDeviceTensor inputs(queue, numberOfJets);
    inputs.zeroInitialise(queue);

    auto const numberOfJetBxOffsets = jetBxLookup.const_view().offset().metadata().size();
    auto const numberOfVertexBxOffsets = vertexBxLookup.const_view().offset().metadata().size();
    auto const numberOfBxs = numberOfJetBxOffsets == 0 ? 0u : numberOfJetBxOffsets - 1u;
    auto const numberOfVertexBxs = numberOfVertexBxOffsets == 0 ? 0u : numberOfVertexBxOffsets - 1u;
    if (numberOfJets == 0 || numberOfBxs == 0 || numberOfBxs != numberOfVertexBxs)
      return inputs;

    auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(numberOfJets, kInputThreads);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        ComputeTorchScriptInputsKernel{},
                        pf.const_view(),
                        jets.const_view(),
                        association_map.const_view(),
                        jetBxLookup.const_view(),
                        vertices.const_view(),
                        vertexBxLookup.const_view(),
                        inputs.view(),
                        static_cast<uint32_t>(numberOfJets),
                        static_cast<uint32_t>(numberOfBxs));
    return inputs;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels::torchscript
