#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernelONNX.h"

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

  namespace {

    template <typename TAcc>
    ALPAKA_FN_ACC float clip(TAcc const& acc, float value, float low, float high) {
      return alpaka::math::min(acc, alpaka::math::max(acc, value, low), high);
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float safeLog(TAcc const& acc, float value) {
      return alpaka::math::log(acc, alpaka::math::max(acc, value, 1.e-8f));
    }

    template <typename TAcc>
    ALPAKA_FN_ACC float preprocess(TAcc const& acc,
                                   float value,
                                   float subtractBy,
                                   float multiplyBy,
                                   float clipMin = -5.f,
                                   float clipMax = 5.f) {
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

    // ---------------------------------------------------------------------------
    // Warp-level reduction via shfl_xor butterfly, as in the PyTorch transform.
    // On the single-thread CPU backend warp::getSize returns 1 so the loop
    // body never executes -- zero overhead.
    // ---------------------------------------------------------------------------
    template <typename TAcc>
    ALPAKA_FN_ACC ALPAKA_FN_INLINE float warpReduceSum(const TAcc& acc, float val) {
      for (int32_t offset = alpaka::warp::getSize(acc) / 2; offset > 0; offset >>= 1) {
        val += alpaka::warp::shfl_xor(acc, val, offset);
      }
      return val;
    }

    // ---------------------------------------------------------------------------
    // Main kernel. Same two-pass algorithm as the PyTorch transform
    // (ComputeParticleTransformerFeaturesKernel): the features written are
    // identical, but the per-jet eigen blocks are transposed to the
    // [features, particles] layout expected by the ONNX models.
    //
    // Grid:  one block per BX
    // Block: 256 threads, each handling one jet constituent
    // ---------------------------------------------------------------------------

    class ComputeParticleTransformerFeaturesKernelONNX {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    PFCandidateDeviceCollection::ConstView pf,
                                    AssociationMapDevice::ConstView clusters_cands,
                                    BxLookupDevice::ConstView jetBxLookup,
                                    VertexDeviceCollection::ConstView vertices,
                                    BxLookupDevice::ConstView vertexBxLookup,
                                    SoftJetONNXInputDeviceTensor::View input_tensors) const {
        constexpr uint32_t maxParticles = ParTPFMask::ColsAtCompileTime;
        [[maybe_unused]] constexpr bool single_thread =
            requires_single_thread_per_block<TAcc>::value;

        // Shared accumulator: one struct = one cache line.
        struct JetAxis {
          float px, py, pz, e, pt, eta, phi;
        };
        auto& ax = alpaka::declareSharedVar<JetAxis, __COUNTER__>(acc);
        auto& vertex_z0 = alpaka::declareSharedVar<float, __COUNTER__>(acc);

        const uint32_t grid_dim =
            alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

        for (uint32_t bx_idx : independent_groups(acc, grid_dim)) {
          // ---- vertex z0 for this BX ----------------------------------------
          const uint32_t vbegin = vertexBxLookup.offset().offset()[bx_idx];
          const uint32_t vend = vertexBxLookup.offset().offset()[bx_idx + 1];
          if (vend <= vbegin)
            continue;
          if (once_per_block(acc))
            vertex_z0 = vertices.z0()[vbegin];
          alpaka::syncBlockThreads(acc);

          // ---- jet range for this BX ----------------------------------------
          const uint32_t jet_begin = jetBxLookup.offset().offset()[bx_idx];
          const uint32_t jet_end = jetBxLookup.offset().offset()[bx_idx + 1];
          if (jet_end <= jet_begin)
            continue;

          for (uint32_t jet_idx : independent_group_elements(acc, jet_end - jet_begin)) {
            const uint32_t block_idx = jet_begin + jet_idx;

            const uint32_t begin = clusters_cands.offset().offset()[block_idx];
            const uint32_t end = clusters_cands.offset().offset()[block_idx + 1];
            if (begin >= end)
              continue;

            const uint32_t input_size =
                (end - begin > maxParticles) ? maxParticles : (end - begin);
            auto input_tensor = input_tensors[block_idx];

            // ------------------------------------------------------------------
            // Pass 1: parallel 4-momentum sum, reduced via warp butterfly and one
            // shared-memory atomic per warp.
            // ------------------------------------------------------------------
            if (once_per_block(acc)) {
              ax.px = 0.f;
              ax.py = 0.f;
              ax.pz = 0.f;
              ax.e = 0.f;
            }
            alpaka::syncBlockThreads(acc);

            float tpx = 0.f, tpy = 0.f, tpz = 0.f, te = 0.f;

            for (uint32_t tid : independent_group_elements(acc, input_size)) {
              const uint32_t p = clusters_cands.index().index()[begin + tid];
              const float pt = pf.pt()[p];
              const float eta = pf.eta()[p];
              const float phi = pf.phi()[p];

              if (pt > 0.f) {
                float sinPhi, cosPhi;
                alpaka::math::sincos(acc, phi, sinPhi, cosPhi);
                const float sinhEta = alpaka::math::sinh(acc, eta);
                const float coshEta = alpaka::math::cosh(acc, eta);
                tpx += pt * cosPhi;
                tpy += pt * sinPhi;
                tpz += pt * sinhEta;
                te += pt * coshEta;
              }
            }

            tpx = warpReduceSum(acc, tpx);
            tpy = warpReduceSum(acc, tpy);
            tpz = warpReduceSum(acc, tpz);
            te = warpReduceSum(acc, te);

            // lane 0 of each warp writes the warp sum into shared memory
            {
              const uint32_t thread_idx_in_block =
                  alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0];
              const uint32_t warp_size =
                  static_cast<uint32_t>(alpaka::warp::getSize(acc));
              if (thread_idx_in_block % warp_size == 0) {
                alpaka::atomicAdd(acc, &ax.px, tpx, alpaka::hierarchy::Threads{});
                alpaka::atomicAdd(acc, &ax.py, tpy, alpaka::hierarchy::Threads{});
                alpaka::atomicAdd(acc, &ax.pz, tpz, alpaka::hierarchy::Threads{});
                alpaka::atomicAdd(acc, &ax.e, te, alpaka::hierarchy::Threads{});
              }
            }
            alpaka::syncBlockThreads(acc);

            // Derive the jet axis on thread 0 and broadcast via shared memory.
            if (once_per_block(acc)) {
              const float pt2 = ax.px * ax.px + ax.py * ax.py;
              ax.pt = alpaka::math::sqrt(acc, pt2);
              ax.eta = alpaka::math::asinh(acc, ax.pz /
                                               alpaka::math::max(acc, ax.pt, 1.0e-8f));
              ax.phi = alpaka::math::atan2(acc, ax.py, ax.px);
            }
            alpaka::syncBlockThreads(acc);

            const float jetPt = ax.pt;
            const float jetEta = ax.eta;
            const float jetPhi = ax.phi;

            // ------------------------------------------------------------------
            // Pass 2: fill the ParticleTransformer tensors, transposed with
            // respect to the PyTorch layout: feature/component first, then the
            // constituent index.
            // ------------------------------------------------------------------
            for (uint32_t tid : independent_group_elements(acc, input_size)) {
              const uint32_t p = clusters_cands.index().index()[begin + tid];
              const float pt = pf.pt()[p];
              const float eta = pf.eta()[p];
              const float phi = pf.phi()[p];

              if (pt <= 0.f)
                continue;  // tensor was zero-initialised; nothing to write

              const float energy = pt * alpaka::math::cosh(acc, eta);
              const float px = pt * alpaka::math::cos(acc, phi);
              const float py = pt * alpaka::math::sin(acc, phi);
              const float pz = pt * alpaka::math::sinh(acc, eta);
              const float etaSign = jetEta >= 0.f ? 1.f : -1.f;
              const float deta = (eta - jetEta) * etaSign;
              const float dphi = deltaPhi(acc, phi, jetPhi);
              const float dr = alpaka::math::sqrt(acc, deta * deta + dphi * dphi);
              const int absPdgId = alpaka::math::abs(acc, static_cast<int>(pf.pdgid()[p]));

              input_tensor.features()(0, tid) = preprocess(acc, safeLog(acc, pt), 1.7f, 0.7f);
              input_tensor.features()(1, tid) = preprocess(acc, safeLog(acc, energy), 2.f, 0.7f);
              input_tensor.features()(2, tid) = preprocess(
                  acc, safeLog(acc, pt / alpaka::math::max(acc, jetPt, 1.e-8f)), -4.7f, 0.7f);
              input_tensor.features()(3, tid) = preprocess(
                  acc, safeLog(acc, energy / alpaka::math::max(acc, ax.e, 1.e-8f)), -4.7f, 0.7f);
              input_tensor.features()(4, tid) = preprocess(acc, dr, 0.2f, 4.f);
              input_tensor.features()(5, tid) = deta;
              input_tensor.features()(6, tid) = dphi;
              input_tensor.features()(7, tid) = pf.z0()[p] - vertex_z0;
              input_tensor.features()(8, tid) = pf.dxy()[p];
              input_tensor.features()(9, tid) = pf.puppiw()[p];
              input_tensor.features()(10, tid) =
                  (absPdgId == 211 || absPdgId == 321 || absPdgId == 2212) ? 1.f : 0.f;
              input_tensor.features()(11, tid) =
                  (absPdgId == 130 || absPdgId == 2112) ? 1.f : 0.f;
              input_tensor.features()(12, tid) = absPdgId == 22 ? 1.f : 0.f;
              input_tensor.features()(13, tid) = absPdgId == 11 ? 1.f : 0.f;
              input_tensor.features()(14, tid) = absPdgId == 13 ? 1.f : 0.f;

              input_tensor.vectors()(0, tid) = px;
              input_tensor.vectors()(1, tid) = py;
              input_tensor.vectors()(2, tid) = pz;
              input_tensor.vectors()(3, tid) = energy;

              input_tensor.mask()(0, tid) = 1.f;
            }  // pass 2

          }  // jet loop
        }  // BX loop
      }  // operator()
    };  // class ComputerParticleTransformerFeaturesKernelONNX

    // Copy the signal (b jet) softmax score into the per-jet output column.
    class ExtractSoftJetScoresKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    SoftJetONNXScoresDevice::ConstView scores,
                                    SoftJetOutputDeviceTensor::View output) const {
        const uint32_t size = output.metadata().size();
        for (uint32_t i : uniform_elements(acc, size)) {
          output.output()[i] = scores.b()[i];
        }
      }
    };

  }  // namespace

  // ---------------------------------------------------------------------------
  // Host-side launches
  // ---------------------------------------------------------------------------

  SoftJetONNXInputDeviceTensor transformONNX(
      Queue& queue,
      const PFCandidateDeviceCollection& pf,
      const AssociationMapDevice& association_map,
      const BxLookupDevice& jetBxLookup,
      const VertexDeviceCollection& vertices,
      const BxLookupDevice& vertexBxLookup) {
    const auto num_clusters =
        association_map.const_view().offset().metadata().size() - 1;
    auto input_tensor = SoftJetONNXInputDeviceTensor(queue, num_clusters);
    input_tensor.zeroInitialise(queue);

    const unsigned int nbx =
        jetBxLookup.const_view().offset().metadata().size() - 1;

    constexpr uint32_t threadsPerBlock = 256;
    auto workDiv = make_workdiv<Acc1D>(nbx, threadsPerBlock);

    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        ComputeParticleTransformerFeaturesKernelONNX{},
                        pf.const_view(),
                        association_map.const_view(),
                        jetBxLookup.const_view(),
                        vertices.const_view(),
                        vertexBxLookup.const_view(),
                        input_tensor.view());

    return input_tensor;
  }

  void extractScores(Queue& queue, const SoftJetONNXScoresDevice& scores, SoftJetOutputDeviceTensor& output) {
    const uint32_t size = output.view().metadata().size();
    if (size == 0)
      return;

    constexpr uint32_t threadsPerBlock = 256;
    auto workDiv = make_workdiv<Acc1D>(size, threadsPerBlock);

    alpaka::exec<Acc1D>(queue, workDiv, ExtractSoftJetScoresKernel{}, scores.const_view(), output.view());
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels
