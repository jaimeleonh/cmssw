#include "PhysicsTools/ParticleTransformerAlpakaFP32/plugins/alpaka/TransformKernel.h"

#include <cstdint>
#include <type_traits>

#include "DataFormats/L1ScoutingSoA/interface/SoftJetTensorSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"
#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelConstants.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::partfp32 {
  namespace {
    using namespace cms::alpakatools;

    template <typename TAcc>
    ALPAKA_FN_ACC float clip(TAcc const& acc, float x, float low, float high) {
      return alpaka::math::min(acc, alpaka::math::max(acc, x, low), high);
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float safeLog(TAcc const& acc, float x) {
      return alpaka::math::log(acc, alpaka::math::max(acc, x, 1.e-8f));
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float preprocess(
        TAcc const& acc, float value, float subtractBy, float multiplyBy, float clipMin = -5.f, float clipMax = 5.f) {
      return clip(acc, (value - subtractBy) * multiplyBy, clipMin, clipMax);
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float deltaPhi(TAcc const&, float phi1, float phi2) {
      constexpr float pi = 3.14159265358979323846f;
      constexpr float twoPi = 2.f * pi;
      float dphi = phi1 - phi2;
      if (dphi > pi)
        dphi -= twoPi;
      if (dphi < -pi)
        dphi += twoPi;
      return dphi;
    }

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
      // This deliberately reproduces the massless four-vector convention used
      // by the reference L1 scouting transform.
      return pt * alpaka::math::cosh(acc, eta);
    }

    inline constexpr uint32_t kInputThreads = 32;
    inline constexpr uint32_t kMaxParticlesPerJet = ::partfp32::kMaxParticles;

    // jetBxLookup.offset() is a monotonic prefix sum, so the bunch crossing a
    // jet belongs to is found in log2(numberOfBxs) steps instead of by
    // scanning.  Every thread of the block evaluates it redundantly, which is
    // cheaper than a broadcast through shared memory.
    ALPAKA_FN_ACC inline uint32_t bunchCrossingOf(::l1sc::BxLookupSoA::ConstView jetBxLookup,
                                                  uint32_t numberOfBxs,
                                                  uint32_t jet) {
      uint32_t low = 0;
      uint32_t high = numberOfBxs;  // answer lies in [low, high)
      while (high - low > 1) {
        uint32_t const middle = low + (high - low) / 2;
        if (static_cast<uint32_t>(jetBxLookup.offset()[middle].offset()) <= jet)
          low = middle;
        else
          high = middle;
      }
      return low;
    }

    // One block of kInputThreads threads owns one jet.  The previous version
    // looped over bunch crossings and put a barrier inside a divergent
    // `independent_group_elements` loop, which is undefined behaviour on the
    // GPU back-ends whenever the number of jets in a bunch crossing is not a
    // multiple of the block size; it also left one thread summing every jet
    // four-vector while the rest of the block waited.  Here every block-wide
    // decision is uniform, so all barriers are reached by the whole block.
    class ComputeInputsKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    l1sc::PFCandidateDeviceCollection::ConstView pf,
                                    l1sc::AssociationMapDevice::ConstView associationMap,
                                    ::l1sc::BxLookupSoA::ConstView jetBxLookup,
                                    l1sc::VertexDeviceCollection::ConstView vertices,
                                    ::l1sc::BxLookupSoA::ConstView vertexBxLookup,
                                    l1sc::SoftJetInputDeviceTensor::View inputs,
                                    uint32_t numberOfJets,
                                    uint32_t numberOfBxs) const {
        auto& cachePt = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cacheEta = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePhi = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePx = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePy = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cachePz = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& cacheEnergy = alpaka::declareSharedVar<float[kMaxParticlesPerJet], __COUNTER__>(acc);
        auto& jet4 = alpaka::declareSharedVar<float[7], __COUNTER__>(acc);  // px py pz e pt eta phi

        auto const lane = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const lanes = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const block = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        for (uint32_t jet = block; jet < numberOfJets; jet += blocks) {
          auto const bx = bunchCrossingOf(jetBxLookup, numberOfBxs, jet);
          auto const vertexBegin = vertexBxLookup.offset()[bx].offset();
          auto const vertexEnd = vertexBxLookup.offset()[bx + 1].offset();
          // A bunch crossing without a reconstructed vertex has no usable z0
          // reference; its jets keep the zero-initialised mask, exactly as in
          // the reference transform.
          if (vertexEnd <= vertexBegin)
            continue;
          float const primaryVertexZ = vertices.z0()[vertexBegin];

          auto const begin = associationMap.offset()[jet].offset();
          auto const end = associationMap.offset()[jet + 1].offset();
          auto const available = static_cast<uint32_t>(end - begin);
          auto const count = available < kMaxParticlesPerJet ? available : kMaxParticlesPerJet;
          if (count == 0)
            continue;

          // Phase 1: every particle's four-vector in parallel.  Candidates
          // with a non-positive pt contribute exact zeros, which leaves the
          // jet sum below bit-identical to the serial reference.
          for (uint32_t i = lane; i < count; i += lanes) {
            auto const candidate = associationMap.index()[begin + i].index();
            auto const candidatePt = pf.pt()[candidate];
            auto const candidateEta = pf.eta()[candidate];
            auto const candidatePhi = pf.phi()[candidate];
            bool const usable = candidatePt > 0.f;
            cachePt[i] = candidatePt;
            cacheEta[i] = candidateEta;
            cachePhi[i] = candidatePhi;
            cachePx[i] = usable ? px(acc, candidatePt, candidatePhi) : 0.f;
            cachePy[i] = usable ? py(acc, candidatePt, candidatePhi) : 0.f;
            cachePz[i] = usable ? pz(acc, candidatePt, candidateEta) : 0.f;
            cacheEnergy[i] = usable ? energy(acc, candidatePt, candidateEta) : 0.f;
          }
          alpaka::syncBlockThreads(acc);

          // Phase 2: the jet four-vector.  The accumulation stays serial and
          // in index order; a tree reduction would reassociate the additions
          // and change the last bits of every downstream feature.
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

          float const jetEnergy = jet4[3];
          float const jetPt = jet4[4];
          float const jetEta = jet4[5];
          float const jetPhi = jet4[6];
          // Preserve the reference transform's eta-reflection convention.
          float const etaSign = jetEta >= 0.f ? 1.f : -1.f;

          // Phase 3: the per-particle record, one thread per particle.
          auto output = inputs[jet];
          for (uint32_t i = lane; i < count; i += lanes) {
            if (cachePt[i] <= 0.f)
              continue;
            auto const candidate = associationMap.index()[begin + i].index();
            auto const candidatePt = cachePt[i];
            auto const candidateEta = cacheEta[i];
            auto const candidatePhi = cachePhi[i];
            auto const candidateEnergy = cacheEnergy[i];
            auto const deta = (candidateEta - jetEta) * etaSign;
            auto const dphi = deltaPhi(acc, candidatePhi, jetPhi);
            auto const dr = alpaka::math::sqrt(acc, deta * deta + dphi * dphi);
            auto const absPdgId = alpaka::math::abs(acc, static_cast<int>(pf.pdgid()[candidate]));

            output.points()(i, 0) = deta;
            output.points()(i, 1) = dphi;

            // Exact pf_features order and preprocessing from b_kinadd.yaml.
            output.features()(i, 0) = preprocess(acc, safeLog(acc, candidatePt), 1.7f, 0.7f);
            output.features()(i, 1) = preprocess(acc, safeLog(acc, candidateEnergy), 2.0f, 0.7f);
            output.features()(i, 2) =
                preprocess(acc, safeLog(acc, candidatePt / alpaka::math::max(acc, jetPt, 1.e-8f)), -4.7f, 0.7f);
            output.features()(i, 3) = preprocess(
                acc, safeLog(acc, candidateEnergy / alpaka::math::max(acc, jetEnergy, 1.e-8f)), -4.7f, 0.7f);
            output.features()(i, 4) = preprocess(acc, dr, 0.2f, 4.f);
            output.features()(i, 5) = deta;
            output.features()(i, 6) = dphi;
            output.features()(i, 7) = pf.z0()[candidate] - primaryVertexZ;
            output.features()(i, 8) = pf.dxy()[candidate];
            output.features()(i, 9) = pf.puppiw()[candidate];
            output.features()(i, 10) = (absPdgId == 211 || absPdgId == 321 || absPdgId == 2212) ? 1.f : 0.f;
            output.features()(i, 11) = (absPdgId == 130 || absPdgId == 2112) ? 1.f : 0.f;
            output.features()(i, 12) = absPdgId == 22 ? 1.f : 0.f;
            output.features()(i, 13) = absPdgId == 11 ? 1.f : 0.f;
            output.features()(i, 14) = absPdgId == 13 ? 1.f : 0.f;

            output.vectors()(i, 0) = cachePx[i];
            output.vectors()(i, 1) = cachePy[i];
            output.vectors()(i, 2) = cachePz[i];
            output.vectors()(i, 3) = candidateEnergy;
            output.mask()(i) = 1.f;
          }
          // The cache is rewritten by the next jet of this block.
          alpaka::syncBlockThreads(acc);
        }
      }
    };
  }  // namespace

  l1sc::SoftJetInputDeviceTensor makeParticleTransformerFP32Inputs(
      Queue& queue,
      l1sc::PFCandidateDeviceCollection const& pf,
      l1sc::AssociationMapDevice const& associationMap,
      l1sc::BxLookupDevice const& jetBxLookup,
      l1sc::VertexDeviceCollection const& vertices,
      l1sc::BxLookupDevice const& vertexBxLookup) {
    auto const numberOfOffsets = associationMap.const_view().offset().metadata().size();
    auto const numberOfJets = numberOfOffsets == 0 ? 0 : numberOfOffsets - 1;
    l1sc::SoftJetInputDeviceTensor inputs(queue, numberOfJets);
    inputs.zeroInitialise(queue);  // zero padding and mask for particles [count,16)

    auto const numberOfJetBxOffsets = jetBxLookup.const_view().offset().metadata().size();
    auto const numberOfVertexBxOffsets = vertexBxLookup.const_view().offset().metadata().size();
    auto const numberOfBxs = numberOfJetBxOffsets == 0 ? 0 : numberOfJetBxOffsets - 1;
    auto const numberOfVertexBxs = numberOfVertexBxOffsets == 0 ? 0 : numberOfVertexBxOffsets - 1;
    if (numberOfJets == 0 || numberOfBxs == 0 || numberOfBxs != numberOfVertexBxs)
      return inputs;

    // One block per jet: the work is now spread over every jet of the event
    // instead of over bunch crossings, so a single busy bunch crossing can no
    // longer serialise the launch.
    auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(static_cast<uint32_t>(numberOfJets), kInputThreads);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        ComputeInputsKernel{},
                        pf.const_view(),
                        associationMap.const_view(),
                        jetBxLookup.const_view(),
                        vertices.const_view(),
                        vertexBxLookup.const_view(),
                        inputs.view(),
                        static_cast<uint32_t>(numberOfJets),
                        static_cast<uint32_t>(numberOfBxs));
    return inputs;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::partfp32
