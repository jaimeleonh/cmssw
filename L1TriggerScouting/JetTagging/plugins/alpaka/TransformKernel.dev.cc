#include "L1TriggerScouting/JetTagging/plugins/alpaka/TransformKernel.h"

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace cms::alpakatools;

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
  // Warp-level reduction via shfl_xor butterfly.
  //
  // alpaka::warp::shfl_xor(acc, v, mask) XORs the caller's lane index with
  // mask and exchanges values -- equivalent to __shfl_xor_sync on CUDA.
  // A butterfly reduction over log2(warpSize) steps reduces all lanes to
  // the same value in every lane (i.e. a warp-level sum broadcast).
  //
  // On the single-thread CPU backend warp::getSize returns 1 so the loop
  // body never executes -- zero overhead.
  // ---------------------------------------------------------------------------
  template <typename TAcc>
  ALPAKA_FN_ACC ALPAKA_FN_INLINE float warpReduceSum(const TAcc& acc, float val) {
    // shfl_xor is unconditionally available for all alpaka backends that
    // support warp intrinsics; on WarpSingleThread it is a no-op returning val.
    for (int32_t offset = alpaka::warp::getSize(acc) / 2; offset > 0; offset >>= 1) {
      val += alpaka::warp::shfl_xor(acc, val, offset);
    }
    return val;  // lane 0 (and every other lane) now holds the warp sum
  }

  // ---------------------------------------------------------------------------
  // Main kernel
  //
  // Grid:  one block per BX
  // Block: 256 threads, each handling one jet constituent
  //
  // Two-pass design, single global-memory read per particle:
  //   Pass 1  -- parallel 4-momentum accumulation via warp butterfly + one
  //             shared-memory atomic per warp.  After syncBlockThreads the
  //             jet axis is derived on thread 0 and broadcast via shared mem.
  //   Pass 2  -- feature writing, re-reading the same cache lines (hot in L1).
  // ---------------------------------------------------------------------------

  class ComputeParticleTransformerFeaturesKernel {
  public:
    template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
    ALPAKA_FN_ACC void operator()(
        TAcc const& acc,
        PFCandidateDeviceCollection::ConstView pf,
        AssociationMapDevice::ConstView clusters_cands,
        BxLookupDevice::ConstView jetBxLookup,
        VertexDeviceCollection::ConstView vertices,
        BxLookupDevice::ConstView vertexBxLookup,
        SoftJetInputDeviceTensor::View input_tensors) const {

      constexpr uint32_t maxParticles = PFFeatures::RowsAtCompileTime;
      [[maybe_unused]] constexpr bool single_thread =
          requires_single_thread_per_block<TAcc>::value;

      // Shared accumulator: one struct = one cache line.
      struct JetAxis { float px, py, pz, e, pt, eta, phi; };
      auto& ax = alpaka::declareSharedVar<JetAxis, __COUNTER__>(acc);
      auto& vertex_z0 = alpaka::declareSharedVar<float, __COUNTER__>(acc);

      const uint32_t grid_dim =
          alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0];

      for (uint32_t bx_idx : independent_groups(acc, grid_dim)) {

        // ---- vertex z0 for this BX ----------------------------------------
        const uint32_t vbegin = vertexBxLookup.offset().offset()[bx_idx];
        const uint32_t vend   = vertexBxLookup.offset().offset()[bx_idx + 1];
        if (vend <= vbegin)
          continue;
        if (once_per_block(acc))
          vertex_z0 = vertices.z0()[vbegin];
        alpaka::syncBlockThreads(acc);

        // ---- jet range for this BX ----------------------------------------
        const uint32_t jet_begin = jetBxLookup.offset().offset()[bx_idx];
        const uint32_t jet_end   = jetBxLookup.offset().offset()[bx_idx + 1];
        if (jet_end <= jet_begin)
          continue;

        for (uint32_t jet_idx : independent_group_elements(acc, jet_end - jet_begin)) {
          const uint32_t block_idx = jet_begin + jet_idx;

          const uint32_t begin = clusters_cands.offset().offset()[block_idx];
          const uint32_t end   = clusters_cands.offset().offset()[block_idx + 1];
          if (begin >= end)
            continue;

          const uint32_t input_size =
              (end - begin > maxParticles) ? maxParticles : (end - begin);
          auto input_tensor = input_tensors[block_idx];

          // ------------------------------------------------------------------
          // Pass 1: parallel 4-momentum sum
          //
          // Each thread accumulates a thread-local partial sum, then a
          // warp butterfly reduces it.  One atomic per warp (not per thread)
          // writes into shared memory.  On single-thread backend the warp
          // reduce is a no-op and the atomic degenerates to a single store.
          // ------------------------------------------------------------------

          // Reset shared accumulator exactly once.
          if (once_per_block(acc)) {
            ax.px = 0.f; ax.py = 0.f; ax.pz = 0.f; ax.e = 0.f;
          }
          alpaka::syncBlockThreads(acc);

          float tpx = 0.f, tpy = 0.f, tpz = 0.f, te = 0.f;

          for (uint32_t tid : independent_group_elements(acc, input_size)) {
            const uint32_t p = clusters_cands.index().index()[begin + tid];
            const float pt  = pf.pt()[p];
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
              te  += pt * coshEta;
            }
          }

          // Warp-level butterfly reduction (shfl_xor, no lane-ID needed).
          tpx = warpReduceSum(acc, tpx);
          tpy = warpReduceSum(acc, tpy);
          tpz = warpReduceSum(acc, tpz);
          te  = warpReduceSum(acc, te);

          // Only the thread whose index is a multiple of warpSize holds the
          // warp sum after shfl_xor (actually *all* lanes do in a butterfly),
          // so we use atomicAdd from every thread -- this is safe because
          // shfl_xor broadcasts the full sum to every lane in the warp,
          // which would double-count. Instead we gate on lane index 0:
          // derive lane index portably as threadIdx % warpSize.
          {
            const uint32_t thread_idx_in_block =
                alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0];
            const uint32_t warp_size =
                static_cast<uint32_t>(alpaka::warp::getSize(acc));
            // lane 0 of each warp writes the warp sum into shared memory
            if (thread_idx_in_block % warp_size == 0) {
              alpaka::atomicAdd(acc, &ax.px, tpx, alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &ax.py, tpy, alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &ax.pz, tpz, alpaka::hierarchy::Threads{});
              alpaka::atomicAdd(acc, &ax.e,  te,  alpaka::hierarchy::Threads{});
            }
          }
          alpaka::syncBlockThreads(acc);

          // Derive jet axis on thread 0 and broadcast via shared memory.
          if (once_per_block(acc)) {
            const float pt2 = ax.px * ax.px + ax.py * ax.py;
            ax.pt  = alpaka::math::sqrt(acc, pt2);
            ax.eta = alpaka::math::asinh(acc, ax.pz /
                         alpaka::math::max(acc, ax.pt, 1.0e-8f));
            ax.phi = alpaka::math::atan2(acc, ax.py, ax.px);
          }
          alpaka::syncBlockThreads(acc);

          // Read jet axis into registers once so all threads share the values
          // without repeated shared-memory loads inside the feature loop.
          const float jetPt   = ax.pt;
          const float jetEta  = ax.eta;
          const float jetPhi  = ax.phi;
          // ------------------------------------------------------------------
          // Pass 2: fill the legacy ParticleTransformer tensor layout.
          // ------------------------------------------------------------------
          for (uint32_t tid : independent_group_elements(acc, input_size)) {
            const uint32_t p = clusters_cands.index().index()[begin + tid];
            const float pt  = pf.pt()[p];
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

            input_tensor.points()(tid, 0) = deta;
            input_tensor.points()(tid, 1) = dphi;

            input_tensor.features()(tid, 0) = preprocess(acc, safeLog(acc, pt), 1.7f, 0.7f);
            input_tensor.features()(tid, 1) = preprocess(acc, safeLog(acc, energy), 2.f, 0.7f);
            input_tensor.features()(tid, 2) = preprocess(
                acc, safeLog(acc, pt / alpaka::math::max(acc, jetPt, 1.e-8f)), -4.7f, 0.7f);
            input_tensor.features()(tid, 3) = preprocess(
                acc, safeLog(acc, energy / alpaka::math::max(acc, ax.e, 1.e-8f)), -4.7f, 0.7f);
            input_tensor.features()(tid, 4) = preprocess(acc, dr, 0.2f, 4.f);
            input_tensor.features()(tid, 5) = deta;
            input_tensor.features()(tid, 6) = dphi;
            input_tensor.features()(tid, 7) = pf.z0()[p] - vertex_z0;
            input_tensor.features()(tid, 8) = pf.dxy()[p];
            input_tensor.features()(tid, 9) = pf.puppiw()[p];
            input_tensor.features()(tid, 10) =
                (absPdgId == 211 || absPdgId == 321 || absPdgId == 2212) ? 1.f : 0.f;
            input_tensor.features()(tid, 11) =
                (absPdgId == 130 || absPdgId == 2112) ? 1.f : 0.f;
            input_tensor.features()(tid, 12) = absPdgId == 22 ? 1.f : 0.f;
            input_tensor.features()(tid, 13) = absPdgId == 11 ? 1.f : 0.f;
            input_tensor.features()(tid, 14) = absPdgId == 13 ? 1.f : 0.f;

            input_tensor.vectors()(tid, 0) = px;
            input_tensor.vectors()(tid, 1) = py;
            input_tensor.vectors()(tid, 2) = pz;
            input_tensor.vectors()(tid, 3) = energy;

            input_tensor.mask()(tid) = 1.f;
          }  // pass 2

        }  // jet loop
      }  // BX loop
    }  // operator()
  };  // class

  // ---------------------------------------------------------------------------
  // Host-side launch
  // ---------------------------------------------------------------------------

  SoftJetInputDeviceTensor transform(
      Queue& queue,
      const PFCandidateDeviceCollection& pf,
      const AssociationMapDevice& association_map,
      const BxLookupDevice& jetBxLookup,
      const VertexDeviceCollection& vertices,
      const BxLookupDevice& vertexBxLookup) {

    const auto num_clusters =
        association_map.const_view().offset().metadata().size() - 1;
    auto input_tensor = SoftJetInputDeviceTensor(queue, num_clusters);
    input_tensor.zeroInitialise(queue);

    const unsigned int nbx =
        jetBxLookup.const_view().offset().metadata().size() - 1;

    constexpr uint32_t threadsPerBlock = 256;
    auto workDiv = make_workdiv<Acc1D>(nbx, threadsPerBlock);

    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        ComputeParticleTransformerFeaturesKernel{},
                        pf.const_view(),
                        association_map.const_view(),
                        jetBxLookup.const_view(),
                        vertices.const_view(),
                        vertexBxLookup.const_view(),
                        input_tensor.view());

    return input_tensor;
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels