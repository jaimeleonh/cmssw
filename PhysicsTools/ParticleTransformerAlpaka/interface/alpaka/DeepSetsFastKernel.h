#ifndef PhysicsTools_ParticleTransformerAlpaka_DeepSetsFastKernel_h
#define PhysicsTools_ParticleTransformerAlpaka_DeepSetsFastKernel_h

// Block-cooperative INT8 DeepSets (deepset_int8.py).
//
// One Alpaka block evaluates one jet, exactly as the ParticleTransformer fast
// kernel does, but the graph is the DeepSets one:
//
//   phi (per particle)    FloatLinear(15 -> 32) -> BatchNorm -> ReLU
//                         Linear INT8(32 -> 32) -> BatchNorm -> ReLU
//   aggregation           masked FP32 sum over the particle dimension
//   rho                   Linear INT8(32 -> 32) -> BatchNorm -> ReLU
//                         FloatLinear(32 -> 2)  -> softmax
//
// The two INT8 projections quantize their input activation with the scale
// learned during quantization-aware training, accumulate the raw integer
// products in INT32 and apply s_x * s_w once per output channel -- what the
// training graph simulated.  phi.0 and the classifier were left in FP32 by the
// training, and stay FP32: phi.0 reads the raw features, the classifier is 64
// weights directly under the softmax.
//
// Every activation lives in block shared memory (about 4 kB, versus about
// 47 kB for the transformer), so the global traffic of an inference is the
// tiny weight blob, which all resident blocks share through L1/L2.
//
// Differences from the transformer kernel, and why:
//
//   * the INT8 weights are read in the exported row-major [out][in] layout --
//     no transposition pass.  With one 32x32 projection per stage there is no
//     coalescing to win back, and the host side stays a plain blob copy;
//   * masked particles are selected with a ternary, not a multiply: a stale
//     input row may contain arbitrary garbage, and its phi state is then
//     discarded with exactly 0.f, with no arithmetic that could leak a NaN;
//   * all kMaxParticles rows of phi are evaluated unconditionally: the whole
//     graph is a few thousand multiplies, so the work-premultiply by the
//     valid-token count is not worth a dependence on the mask contents.
//
// GeneratedDeepSetsModel.h and GeneratedModel.h declare the same names in
// namespace generated, so this header and ParticleTransformerFastKernel.h can
// never be included in the same translation unit; the helpers duplicated here
// (and kept in a nested namespace) make this header stand alone.

#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/DeepSetsModelDefinition.h"

namespace part::deepsets {

  // The padded sequence length of the jet input SoA, shared with the
  // transformer kernel so the same input tensors can drive either model.
  inline constexpr int kMaxParticles = 16;
  inline constexpr int kThreadsPerJet = 256;

  // Warp layout.  A DeepSets graph is ~25 kMACs per jet: at block granularity
  // the per-jet cost is barrier latency, not arithmetic, so one jet runs per
  // warp instead -- eight jets share a block, each waiting on its own
  // warp-level sync, and the staging of every jet hides behind the evaluation
  // of the other seven.  This is the same shape the graph interpreter uses
  // for small networks (graph/GraphShared.h warp mode).
  inline constexpr int kWarpLanes = 32;
  inline constexpr int kWarpJetsPerBlock = kThreadsPerJet / kWarpLanes;  // 8
  static_assert(kThreadsPerJet % kWarpLanes == 0, "a block must hold a whole number of warps");

  // ---------------------------------------------------------------- shared state

  struct alignas(16) DeepSetsShared {
    alignas(16) float features[kMaxParticles * generated::input_features];
    alignas(16) float x[kMaxParticles * generated::embed_width];        // phi state of every token
    alignas(16) float state[generated::embed_width];                    // the single rho input row
    alignas(16) std::int8_t quant[kMaxParticles * generated::embed_width];  // INT8 staging, row-major
    float mask[kMaxParticles];
    float logits[generated::num_classes];
  };

  static_assert(sizeof(DeepSetsShared) <= 47u * 1024u, "the per-jet shared state grew past its budget");

  // One warp per jet: kWarpJetsPerBlock states share the block's arena.
  struct DeepSetsWarpBlock {
    DeepSetsShared jet[kWarpJetsPerBlock];
  };

  static_assert(sizeof(DeepSetsWarpBlock) <= 47u * 1024u,
                "the warp-mode block state grew past the shared-memory budget");

  // Sixteen INT8 activations, or four packed weight words, in one 16-byte load.
  struct alignas(16) Int4 {
    int v[4];
  };

  // Four INT8 multiply-accumulates in one instruction.  On the GPU this is
  // dp4a: the FP32 pipe retires one multiply-add per lane per clock, dp4a
  // retires four.
  template <typename TAcc>
  ALPAKA_FN_ACC inline int dot4(TAcc const&, int a, int b, int c) {
#if defined(__CUDA_ARCH__) && __CUDA_ARCH__ >= 610
    return __dp4a(a, b, c);
#else
    auto const* x = reinterpret_cast<std::int8_t const*>(&a);
    auto const* y = reinterpret_cast<std::int8_t const*>(&b);
    return c + x[0] * y[0] + x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
#endif
  }

  // ---------------------------------------------------------------- block helpers

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneOf(TAcc const& acc) {
    return static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneCount(TAcc const& acc) {
    return static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
  }

  // A slice of a block that evaluates one jet on its own -- in warp mode, one
  // warp.  It stands in for the accelerator wherever the helpers below ask
  // which lane it is, how many lanes there are, or wait for the others, so
  // the same code runs at block or at warp granularity.
  template <typename TAcc>
  struct WarpAcc {
    TAcc const& acc;
    int lane;
  };

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneOf(WarpAcc<TAcc> const& group) {
    return group.lane;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneCount(WarpAcc<TAcc> const&) {
    return kWarpLanes;
  }

  // Compile-time lane count, for the GEMM register tiles.  The runtime count
  // comes from the launcher and always matches this.
  template <typename>
  struct GroupLanes {
    static constexpr int value = kThreadsPerJet;
  };
  template <typename TAcc>
  struct GroupLanes<WarpAcc<TAcc>> {
    static constexpr int value = kWarpLanes;
  };

  // Wait for the other lanes of whatever is evaluating this jet.
  //
  // For a WarpAcc this has to be a warp-level sync and not a block barrier:
  // warps of a block evaluate different jets, whose staging arrival times
  // differ, and a block barrier would make them wait for each other.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void groupBarrier(TAcc const& acc) {
    alpaka::syncBlockThreads(acc);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void groupBarrier(WarpAcc<TAcc> const& group) {
#if defined(__CUDA_ARCH__)
    __syncwarp();
#else
    // Only correct when the group is the whole block: the CPU test harness,
    // which gives every real thread its own acc, so this degenerates to a
    // per-warp barrier.  The warp launch is compiled for CUDA only anyway
    // (see DeepSetsAlgo.dev.cc), so a multi-group block never reaches this.
    alpaka::syncBlockThreads(group.acc);
#endif
  }

  // ---------------------------------------------------------------- dense layers

  // Choose the register tile: every lane of the group must get work, so the
  // task space is (token tiles) x (column groups) and the tile is the largest
  // one that still covers the group's lanes, exactly as the transformer
  // kernel tiles its GEMMs against the block size.
  constexpr bool tileCovers(int lanes, int outputs, int tokens, int tb, int cb) {
    return outputs % cb == 0 && ((tokens + tb - 1) / tb) * (outputs / cb) >= lanes;
  }

  constexpr int tokenBlockFor(int lanes, int outputs, int tokens) {
    for (int tb = 8; tb > 1; tb /= 2)
      for (int cb = 4; cb >= 1; cb /= 2)
        if (tb <= tokens && tileCovers(lanes, outputs, tokens, tb, cb))
          return tb;
    return 1;
  }

  constexpr int columnBlockFor(int lanes, int outputs, int tokens) {
    for (int tb = 8; tb > 1; tb /= 2)
      for (int cb = 4; cb >= 1; cb /= 2)
        if (tb <= tokens && tileCovers(lanes, outputs, tokens, tb, cb))
          return cb;
    for (int cb = 4; cb > 1; cb /= 2)
      if (tileCovers(lanes, outputs, tokens, 1, cb))
        return cb;
    return 1;
  }

  // out[t][c] = bias[c] + activationScale * scale[c] * sum_k q(input[t][k]) * q(w[c][k])
  //
  // The activations arrive as INT8 in row-major [token][KDIM] and the weights
  // in the exported [out][in] layout, so one 16-byte load fetches sixteen
  // inputs of a token and each 32-bit weight word feeds one dp4a.  The
  // accumulator is INT32: with KDIM = 32 the largest partial sum is
  // 32 * 127 * 127, about 5.2e5, orders below overflow.  Both scales are
  // applied once per output, which is arithmetically identical to what
  // quantization-aware training simulated.
  //
  // KDIM, NOUT, TB and CB are compile-time so the inner loops unroll and the
  // TB x CB accumulator tile stays in registers.  Threads are laid out
  // column-fastest so each thread's weight fetch is one load per four input
  // channels, and the whole warp covers the output row before the next token
  // tile.
  template <int KDIM, int NOUT, int TB, int CB, typename TAcc>
  ALPAKA_FN_ACC inline void gemmInt8(TAcc const& acc,
                                     std::int8_t const* __restrict__ weight,
                                     float const* __restrict__ scale,
                                     float const* __restrict__ bias,
                                     float activationScale,
                                     std::int8_t const* __restrict__ input,
                                     float* __restrict__ output,
                                     int tokens) {
    static_assert(KDIM % 16 == 0, "the input width must allow 16-byte activation loads");
    static_assert(NOUT % CB == 0, "the column block must divide the output width");
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    constexpr int groups = NOUT / CB;
    constexpr int words = KDIM / 4;  // packed 32-bit words per token row
    int const tiles = (tokens + TB - 1) / TB;
    int const tasks = tiles * groups;

    auto const* __restrict__ weightWords = reinterpret_cast<int const*>(weight);
    auto const* __restrict__ inputWords = reinterpret_cast<int const*>(input);

    for (int task = lane; task < tasks; task += lanes) {
      int const tile = task / groups;
      int const column = (task - tile * groups) * CB;
      int const token = tile * TB;

      int sum[TB][CB];
      for (int t = 0; t < TB; ++t)
        for (int c = 0; c < CB; ++c)
          sum[t][c] = 0;

      for (int group = 0; group < words; group += 4) {
        Int4 activation[TB];
        for (int t = 0; t < TB; ++t)
          activation[t] = *reinterpret_cast<Int4 const*>(inputWords + (token + t) * words + group);

        for (int c = 0; c < CB; ++c) {
          int weightWords4[4];
          for (int r = 0; r < 4; ++r)
            weightWords4[r] = weightWords[(column + c) * words + group + r];
          for (int t = 0; t < TB; ++t)
            for (int r = 0; r < 4; ++r)
              sum[t][c] = dot4(acc, activation[t].v[r], weightWords4[r], sum[t][c]);
        }
      }

      for (int t = 0; t < TB; ++t)
        if (token + t < tokens)
          for (int c = 0; c < CB; ++c) {
            int const channel = column + c;
            output[(token + t) * NOUT + channel] =
                static_cast<float>(sum[t][c]) * (activationScale * scale[channel]) + bias[channel];
          }
    }
    groupBarrier(acc);
  }

  // The FP32 phi projection that reads the raw features.  Its weights are the
  // exported FP32 rows [out][in] in the parameter blob, so a column-owning
  // thread walks memory contiguously down k; the projection is 480 multiplies
  // per token and needs no further sophistication.
  template <int KDIM, int NOUT, typename TAcc>
  ALPAKA_FN_ACC inline void gemmFloat(TAcc const& acc,
                                      float const* __restrict__ weight,
                                      float const* __restrict__ bias,
                                      float const* __restrict__ input,
                                      float* __restrict__ output,
                                      int tokens) {
    for (int index = laneOf(acc); index < tokens * NOUT; index += laneCount(acc)) {
      int const row = index / NOUT;
      int const column = index - row * NOUT;
      float value = bias[column];
      for (int k = 0; k < KDIM; ++k)
        value += input[row * KDIM + k] * weight[column * KDIM + k];
      output[index] = value;
    }
    groupBarrier(acc);
  }

  // ---------------------------------------------------------------- pointwise passes

  // Round a block of activations onto the INT8 grid in place.  Quantizing
  // into the INT8 staging area is what feeds every integer GEMM.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void quantizeRows(TAcc const& acc,
                                         float const* __restrict__ input,
                                         std::int8_t* __restrict__ output,
                                         int width,
                                         int rows,
                                         float inverseScale) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    for (int index = lane; index < rows * width; index += lanes)
      output[index] = quantizeOne(input[index], inverseScale);
    groupBarrier(acc);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void applyRelu(TAcc const& acc, float* __restrict__ data, int count) {
    for (int index = laneOf(acc); index < count; index += laneCount(acc))
      data[index] = alpaka::math::max(acc, data[index], 0.f);
    groupBarrier(acc);
  }

  // Elementwise and safe in place: every index is read once and written once
  // by the same thread, so no aliasing hazard arises when input == output.
  template <int SIZE, typename TAcc>
  ALPAKA_FN_ACC inline void applyBatchNorm(TAcc const& acc,
                                           generated::BatchNorm const& parameters,
                                           float const* input,
                                           float* output,
                                           int rows) {
    for (int index = laneOf(acc); index < rows * SIZE; index += laneCount(acc)) {
      int const c = index % SIZE;
      output[index] = (input[index] - parameters.mean[c]) /
                          alpaka::math::sqrt(acc, parameters.var[c] + parameters.eps) * parameters.gamma[c] +
                      parameters.beta[c];
    }
    groupBarrier(acc);
  }

  // ---------------------------------------------------------------- full graph

  // Expects shared.features[t * generated::input_features + f] and shared.mask
  // to be filled and a barrier to have been executed.  Leaves the two class
  // probabilities in shared.logits.
  //
  // `acc` is either a bare accelerator (block layout: one jet per block) or a
  // WarpAcc covering one warp (warp layout: one jet per warp); every pass
  // sizes itself from laneCount(acc), and the register tile of the phi
  // projection from GroupLanes.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void runDeepSets(TAcc const& acc,
                                        generated::ModelView model,
                                        DeepSetsShared& shared) {
    constexpr int kFeatures = generated::input_features;
    constexpr int kWidth = generated::embed_width;
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);

    // ---- phi: FloatLinear(15->32) -> BatchNorm -> ReLU -> INT8 Linear(32->32) -> BatchNorm
    {
      auto const input = generated::inputLinear(model);
      gemmFloat<kFeatures, kWidth>(acc, input.weight, input.bias, shared.features, shared.x, kMaxParticles);
    }
    applyBatchNorm<kWidth>(acc, generated::phiNorm(model, 0), shared.x, shared.x, kMaxParticles);
    applyRelu(acc, shared.x, kMaxParticles * kWidth);
    {
      auto const phi = generated::phiProject(model);
      quantizeRows(acc, shared.x, shared.quant, kWidth, kMaxParticles, 1.f / phi.activation[0]);
      gemmInt8<kWidth, kWidth, tokenBlockFor(GroupLanes<TAcc>::value, kWidth, kMaxParticles), columnBlockFor(GroupLanes<TAcc>::value, kWidth, kMaxParticles)>(
          acc, phi.weight, phi.scale, phi.bias, phi.activation[0], shared.quant, shared.x, kMaxParticles);
    }
    applyBatchNorm<kWidth>(acc, generated::phiNorm(model, 1), shared.x, shared.x, kMaxParticles);
    for (int index = lane; index < kMaxParticles * kWidth; index += lanes) {
      int const token = index / kWidth;
      shared.x[index] = shared.mask[token] > .5f ? alpaka::math::max(acc, shared.x[index], 0.f) : 0.f;
    }
    groupBarrier(acc);

    // ---- aggregation: masked FP32 sum over the particle dimension.  The sum
    // stays in FP32 because its dynamic range grows with the multiplicity;
    // rho.0's activation scale absorbs that at the quantization point.
    for (int c = lane; c < kWidth; c += lanes) {
      float value = 0.f;
      for (int t = 0; t < kMaxParticles; ++t)
        value += shared.x[t * kWidth + c];
      shared.state[c] = value;
    }
    groupBarrier(acc);

    // ---- rho: INT8 Linear(32->32) -> BatchNorm -> ReLU -> FP32 classifier -> softmax
    {
      auto const rho = generated::rhoProject(model);
      quantizeRows(acc, shared.state, shared.quant, kWidth, 1, 1.f / rho.activation[0]);
      gemmInt8<kWidth, kWidth, 1, 1>(acc, rho.weight, rho.scale, rho.bias, rho.activation[0], shared.quant, shared.state, 1);
    }
    applyBatchNorm<kWidth>(acc, generated::rhoNorm(model), shared.state, shared.state, 1);
    applyRelu(acc, shared.state, kWidth);

    {
      auto const head = generated::classifier(model);
      // FP32, row-major [out][in], two outputs under the softmax.
      for (int column = lane; column < generated::num_classes; column += lanes) {
        float value = head.bias[column];
        for (int k = 0; k < kWidth; ++k)
          value += shared.state[k] * head.weight[column * kWidth + k];
        shared.logits[column] = value;
      }
      groupBarrier(acc);
    }
    if (lane == 0) {
      float const maximum = shared.logits[0] > shared.logits[1] ? shared.logits[0] : shared.logits[1];
      float const e0 = alpaka::math::exp(acc, shared.logits[0] - maximum);
      float const e1 = alpaka::math::exp(acc, shared.logits[1] - maximum);
      shared.logits[0] = e0 / (e0 + e1);
      shared.logits[1] = e1 / (e0 + e1);
    }
    groupBarrier(acc);
  }

}  // namespace part::deepsets

#endif