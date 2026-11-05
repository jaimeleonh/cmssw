#ifndef PhysicsTools_ParticleTransformerAlpakaFP32_ParticleTransformerFastKernel_h
#define PhysicsTools_ParticleTransformerAlpakaFP32_ParticleTransformerFastKernel_h

// Block-cooperative FP32 ParticleTransformer.
//
// One Alpaka block evaluates one jet.  Every activation of the graph lives in
// block shared memory (about 47 kB), so the only global traffic during
// inference is the streaming read of the model weights, which is shared by all
// resident blocks through L1/L2.
//
// Three properties drive the layout:
//
//   * dense weights are stored transposed as [in][out] (see ModelTranspose.h),
//     so consecutive threads read consecutive output columns;
//   * every GEMM is register-blocked, with the tile chosen per layer so
//     that the (token tile x column group) task space fills the whole block
//     while reusing each weight load across as many tokens as possible;
//   * only the tokens up to the last valid particle are evaluated, and only
//     the i <= j half of the pairwise interaction matrix, which is exactly
//     symmetric.
//
// Masked keys are forced to -3.4e38 before the softmax, so their weight
// underflows to exactly zero and they cannot influence any surviving token.
// Skipping them is therefore an identity, not an approximation.

#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelDefinition.h"

namespace partfp32 {

  // ---------------------------------------------------------------- shared state

  struct alignas(16) JetShared {
    float x[kMaxParticles * kEmbedDim];         // token states, live for the whole graph
    alignas(16) float na[kNormScratch];         // one tile of normalised tokens, token-contiguous
    alignas(16) float big[kBigScratch];         // QKV | FFN tile | pair tile | class scratch
    float pair[kHeads * kMaxPairs];             // upper-triangular attention bias
    float cls[kEmbedDim];
    float clsNorm[kEmbedDim];
    float reduceA[kReduceThreads];
    float reduceB[kReduceThreads];
    float mean[kMaxParticles + 4];
    float rstd[kMaxParticles + 4];
    float p4[kMaxParticles * kFourVector];
    float mask[kMaxParticles];
    float logits[4];
    int active;  // last unmasked particle plus one
    int tokens;  // `active` rounded up to a multiple of kTokenBlock
    int pairs;   // active * (active + 1) / 2
  };

  // CUDA allows 48 kB of statically declared shared memory per block, and the
  // Alpaka CPU back-ends default to a 47 kB block-shared arena.  Stay clearly
  // below both so that neither limit is a build- or run-time surprise.
  static_assert(sizeof(JetShared) <= 44u * 1024u, "the per-jet shared state grew past its budget");

  // A 16-byte aggregate load.  Reading four adjacent activations in one
  // instruction is what makes the token-contiguous layout worth having: it
  // turns four shared-memory instructions into one, and instruction issue -
  // not shared-memory bandwidth - is what limits the inner loop.
  struct alignas(16) Float4 {
    float v[4];
  };

  // ---------------------------------------------------------------- block helpers

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneOf(TAcc const& acc) {
    return static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneCount(TAcc const& acc) {
    return static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline float geluOf(TAcc const& acc, float value) {
    return .5f * value * (1.f + alpaka::math::erf(acc, value * .7071067811865475f));
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void fill(TAcc const& acc, float* data, int size, float value) {
    for (int index = laneOf(acc); index < size; index += laneCount(acc))
      data[index] = value;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void applyGelu(TAcc const& acc, float* data, int size) {
    for (int index = laneOf(acc); index < size; index += laneCount(acc))
      data[index] = geluOf(acc, data[index]);
    alpaka::syncBlockThreads(acc);
  }

  // ---------------------------------------------------------------- dense layers

  // Choose the register tile.  Every lane of the block must get work, so the
  // task space is (token tiles) x (column groups) and the tile is the largest
  // one that still covers kThreadsPerJet lanes.  A larger TB is preferred
  // first because it divides the weight traffic: each thread loads a weight
  // element once and reuses it for TB tokens.  CB then trades registers for
  // broadcast reads of the activation, which come from shared memory and are
  // nearly free.
  constexpr bool tileCovers(int outputs, int tokens, int tb, int cb) {
    return outputs % cb == 0 && ((tokens + tb - 1) / tb) * (outputs / cb) >= kThreadsPerJet;
  }

  constexpr int tokenBlockFor(int outputs, int tokens) {
    for (int tb = kTokenBlockMax; tb > 1; tb /= 2)
      for (int cb = kColumnBlockMax; cb >= 1; cb /= 2)
        if (tb <= tokens && tileCovers(outputs, tokens, tb, cb))
          return tb;
    return 1;
  }

  constexpr int columnBlockFor(int outputs, int tokens) {
    for (int tb = kTokenBlockMax; tb > 1; tb /= 2)
      for (int cb = kColumnBlockMax; cb >= 1; cb /= 2)
        if (tb <= tokens && tileCovers(outputs, tokens, tb, cb))
          return cb;
    for (int cb = kColumnBlockMax; cb > 1; cb /= 2)
      if (tileCovers(outputs, tokens, 1, cb))
        return cb;
    return 1;
  }

  // out[t][columnBegin + n] = bias[columnBegin + n]
  //                        + sum_k input[t][k] * weightT[k][columnBegin + n]
  //
  // KDIM, WSTRIDE, NOUT, TB and CB are compile-time so the inner loops unroll
  // and the TB x CB accumulator tile stays in registers.  Occupancy here is
  // limited by shared memory rather than by registers, so a large tile is
  // free.  Threads are laid out column-fastest, which keeps each warp's weight
  // fetch a single coalesced transaction.
  // INPUT_T / OUTPUT_T select the token-contiguous layout for the activation
  // operands: element (t, k) then lives at [k * stride + t] instead of
  // [t * stride + k].  With the tokens adjacent, a thread's TB activations for
  // one k are fetched by TB/4 vector loads rather than TB scalar ones.  That
  // matters because an SM retires four instructions per clock and reaches peak
  // FP32 only when all four are fused multiply-adds: at TB=8, CB=1 the scalar
  // form issues 8 LDS + 1 LDG per 8 FFMA, capping the kernel at 47% of the
  // pipe before occupancy is even considered.  The vector form issues 2 + 1
  // per 8, lifting that ceiling to 73%.
  template <int KDIM, int WSTRIDE, int NOUT, int TB, int CB, bool INPUT_T = false, bool OUTPUT_T = false,
            typename TAcc>
  ALPAKA_FN_ACC inline void gemmTiled(TAcc const& acc,
                                      float const* __restrict__ weightT,
                                      float const* __restrict__ bias,
                                      int columnBegin,
                                      float const* __restrict__ input,
                                      int inputStride,
                                      float* __restrict__ output,
                                      int outputStride,
                                      int tokens) {
    static_assert(NOUT % CB == 0, "the column block must divide the output width");
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    constexpr int groups = NOUT / CB;
    int const tiles = (tokens + TB - 1) / TB;
    int const tasks = tiles * groups;

    for (int task = lane; task < tasks; task += lanes) {
      int const tile = task / groups;
      int const column = (task - tile * groups) * CB;
      int const token = tile * TB;

      float const* __restrict__ weight = weightT + columnBegin + column;
      float sum[TB][CB];
      for (int t = 0; t < TB; ++t)
        for (int c = 0; c < CB; ++c)
          sum[t][c] = bias[columnBegin + column + c];

      float const* __restrict__ row = input + (INPUT_T ? token : token * inputStride);
      for (int k = 0; k < KDIM; ++k) {
        float value[CB];
        for (int c = 0; c < CB; ++c)
          value[c] = weight[k * WSTRIDE + c];

        float element[TB];
        if constexpr (INPUT_T) {
          float const* __restrict__ base = row + k * inputStride;
          if constexpr (TB % 4 == 0) {
            // `inputStride` and `token` are both multiples of four, so `base`
            // keeps the 16-byte alignment of the shared buffer.
            for (int q = 0; q < TB / 4; ++q) {
              Float4 const chunk = *reinterpret_cast<Float4 const*>(base + 4 * q);
              for (int r = 0; r < 4; ++r)
                element[4 * q + r] = chunk.v[r];
            }
          } else {
            for (int t = 0; t < TB; ++t)
              element[t] = base[t];
          }
        } else {
          for (int t = 0; t < TB; ++t)
            element[t] = row[t * inputStride + k];
        }

        for (int t = 0; t < TB; ++t)
          for (int c = 0; c < CB; ++c)
            sum[t][c] += element[t] * value[c];
      }

      // A partial trailing tile computes garbage rows from padding that is
      // finite but meaningless; they are simply not written back.
      for (int t = 0; t < TB; ++t)
        if (token + t < tokens)
          for (int c = 0; c < CB; ++c) {
            float const result = sum[t][c];
            if constexpr (OUTPUT_T)
              output[(column + c) * outputStride + token + t] = result;
            else
              output[(token + t) * outputStride + column + c] = result;
          }
    }
    alpaka::syncBlockThreads(acc);
  }

  // Narrow projections (the two-class head) where register blocking is pointless.
  template <int KDIM, int NOUT, typename TAcc>
  ALPAKA_FN_ACC inline void gemmNarrow(TAcc const& acc,
                                       float const* __restrict__ weightT,
                                       float const* __restrict__ bias,
                                       float const* __restrict__ input,
                                       float* __restrict__ output) {
    for (int column = laneOf(acc); column < NOUT; column += laneCount(acc)) {
      float value = bias[column];
      for (int k = 0; k < KDIM; ++k)
        value += input[k] * weightT[k * NOUT + column];
      output[column] = value;
    }
    alpaka::syncBlockThreads(acc);
  }

  // ---------------------------------------------------------------- normalisation

  // Two-pass LayerNorm.  A group of `group` threads reduces one token; the two
  // partial arrays are distinct so no barrier is needed between reading the
  // mean partials and writing the variance partials.
  template <int SIZE, bool INPUT_T = false, bool OUTPUT_T = false, typename TAcc>
  ALPAKA_FN_ACC inline void layerNorm(TAcc const& acc,
                                      generated::Norm const& parameters,
                                      float const* __restrict__ input,
                                      int inputStride,
                                      float* __restrict__ output,
                                      int outputStride,
                                      int tokens,
                                      JetShared& shared) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const reducing = lanes < kReduceThreads ? lanes : kReduceThreads;

    int group = 1;
    while (group < 32 && group * 2 * tokens <= reducing)
      group *= 2;
    int const groups = reducing / group;

    for (int base = 0; base < tokens; base += groups) {
      int const slot = lane % group;
      int const token = base + lane / group;
      bool const contributes = lane < reducing && token < tokens;

      // No execution-space annotation here: nvcc derives it from the enclosing
      // device function and rejects an explicit one on a lambda.
      auto const element = [&](int c) {
        return INPUT_T ? input[c * inputStride + token] : input[token * inputStride + c];
      };

      float partial = 0.f;
      if (contributes) {
        for (int c = slot; c < SIZE; c += group)
          partial += element(c);
      }
      if (lane < reducing)
        shared.reduceA[lane] = partial;
      alpaka::syncBlockThreads(acc);

      float mean = 0.f;
      float variance = 0.f;
      if (contributes) {
        int const first = (lane / group) * group;
        for (int k = 0; k < group; ++k)
          mean += shared.reduceA[first + k];
        mean /= static_cast<float>(SIZE);
        for (int c = slot; c < SIZE; c += group) {
          float const difference = element(c) - mean;
          variance += difference * difference;
        }
      }
      if (lane < reducing)
        shared.reduceB[lane] = variance;
      alpaka::syncBlockThreads(acc);

      if (contributes && slot == 0) {
        int const first = (lane / group) * group;
        float total = 0.f;
        for (int k = 0; k < group; ++k)
          total += shared.reduceB[first + k];
        total /= static_cast<float>(SIZE);
        shared.mean[token] = mean;
        shared.rstd[token] = 1.f / alpaka::math::sqrt(acc, total + parameters.eps);
      }
      alpaka::syncBlockThreads(acc);

      int const covered = (tokens - base < groups ? tokens - base : groups) * SIZE;
      for (int index = lane; index < covered; index += lanes) {
        int const local = index / SIZE;
        int const c = index - local * SIZE;
        int const t = base + local;
        float const source = INPUT_T ? input[c * inputStride + t] : input[t * inputStride + c];
        float const value = (source - shared.mean[t]) * shared.rstd[t] * parameters.gamma[c] +
                            parameters.beta[c];
        if constexpr (OUTPUT_T)
          output[c * outputStride + t] = value;
        else
          output[t * outputStride + c] = value;
      }
      alpaka::syncBlockThreads(acc);
    }
  }

  template <int SIZE, typename TAcc>
  ALPAKA_FN_ACC inline void applyBatchNorm(TAcc const& acc,
                                      generated::BatchNorm const& parameters,
                                      float const* __restrict__ input,
                                      float* __restrict__ output,
                                      int tokens) {
    for (int index = laneOf(acc); index < tokens * SIZE; index += laneCount(acc)) {
      int const c = index % SIZE;
      output[index] = (input[index] - parameters.mean[c]) /
                          alpaka::math::sqrt(acc, parameters.var[c] + parameters.eps) * parameters.gamma[c] +
                      parameters.beta[c];
    }
    alpaka::syncBlockThreads(acc);
  }

  // ---------------------------------------------------------------- pair helpers

  ALPAKA_FN_ACC inline int triangularBase(int row, int size) { return row * size - (row * (row - 1)) / 2; }

  ALPAKA_FN_ACC inline int pairSlot(int i, int j, int size) {
    int const low = i < j ? i : j;
    int const high = i < j ? j : i;
    return triangularBase(low, size) + (high - low);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void pairFromSlot(TAcc const& acc, int slot, int size, int& i, int& j) {
    float const b = 2.f * static_cast<float>(size) + 1.f;
    int row = static_cast<int>((b - alpaka::math::sqrt(acc, b * b - 8.f * static_cast<float>(slot))) * .5f);
    if (row < 0)
      row = 0;
    if (row > size - 1)
      row = size - 1;
    while (row > 0 && triangularBase(row, size) > slot)
      --row;
    while (row < size - 1 && triangularBase(row + 1, size) <= slot)
      ++row;
    i = row;
    j = row + (slot - triangularBase(row, size));
  }

  // ---------------------------------------------------------------- pair embedding

  template <typename TAcc>
  ALPAKA_FN_ACC inline void encodePairs(TAcc const& acc, generated::ModelView model, JetShared& shared) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const size = shared.active;
    int const pairs = shared.pairs;

    float* const first = shared.big;
    float* const second = shared.big + kPairTile * kPairHidden;
    auto const inputNorm = generated::pairInputBN(model);

    for (int begin = 0; begin < pairs; begin += kPairTile) {
      int const remaining = pairs - begin;
      int const raw = remaining < kPairTile ? remaining : kPairTile;
      int const rows = ((raw + kTokenBlock - 1) / kTokenBlock) * kTokenBlock;

      for (int local = lane; local < rows; local += lanes) {
        int const slot = begin + local < pairs ? begin + local : pairs - 1;
        int i = 0;
        int j = 0;
        pairFromSlot(acc, slot, size, i, j);

        float const pxi = shared.p4[i * 4], pyi = shared.p4[i * 4 + 1];
        float const pzi = shared.p4[i * 4 + 2], ei = shared.p4[i * 4 + 3];
        float const pxj = shared.p4[j * 4], pyj = shared.p4[j * 4 + 1];
        float const pzj = shared.p4[j * 4 + 2], ej = shared.p4[j * 4 + 3];

        float const pti = alpaka::math::sqrt(acc, pxi * pxi + pyi * pyi);
        float const ptj = alpaka::math::sqrt(acc, pxj * pxj + pyj * pyj);
        float const rapidityI =
            .5f * alpaka::math::log(acc, 1.f + 2.f * pzi / alpaka::math::max(acc, ei - pzi, 1e-20f));
        float const rapidityJ =
            .5f * alpaka::math::log(acc, 1.f + 2.f * pzj / alpaka::math::max(acc, ej - pzj, 1e-20f));
        float deltaPhi = alpaka::math::atan2(acc, pyi, pxi) - alpaka::math::atan2(acc, pyj, pxj);
        constexpr float pi = 3.14159265358979323846f;
        constexpr float twoPi = 2.f * pi;
        while (deltaPhi > pi)
          deltaPhi -= twoPi;
        while (deltaPhi <= -pi)
          deltaPhi += twoPi;
        float const distance = alpaka::math::sqrt(
            acc, (rapidityI - rapidityJ) * (rapidityI - rapidityJ) + deltaPhi * deltaPhi);
        float const minimumPt = alpaka::math::min(acc, pti, ptj);
        float const sumX = pxi + pxj, sumY = pyi + pyj, sumZ = pzi + pzj, sumE = ei + ej;

        float value[kPairInputs];
        value[0] = alpaka::math::log(acc, alpaka::math::max(acc, minimumPt * distance, 1e-8f));
        value[1] = alpaka::math::log(
            acc, alpaka::math::max(acc, minimumPt / alpaka::math::max(acc, pti + ptj, 1e-8f), 1e-8f));
        value[2] = alpaka::math::log(acc, alpaka::math::max(acc, distance, 1e-8f));
        value[3] = alpaka::math::log(
            acc, alpaka::math::max(acc, sumE * sumE - sumX * sumX - sumY * sumY - sumZ * sumZ, 1e-8f));

        for (int c = 0; c < kPairInputs; ++c)
          first[local * kPairHidden + c] =
              (value[c] - inputNorm.mean[c]) / alpaka::math::sqrt(acc, inputNorm.var[c] + inputNorm.eps) *
                  inputNorm.gamma[c] +
              inputNorm.beta[c];
      }
      alpaka::syncBlockThreads(acc);

      float* source = first;
      float* destination = second;
      for (int layer = 0; layer < 4; ++layer) {
        auto const dense = generated::pairLinear(model, layer);
        auto const norm = generated::pairBN(model, layer);
        if (layer == 0)
          gemmTiled<kPairInputs, kPairHidden, kPairHidden, tokenBlockFor(kPairHidden, kPairTile), columnBlockFor(kPairHidden, kPairTile)>(
              acc, dense.weight, dense.bias, 0, source, kPairHidden, destination, kPairHidden, rows);
        else if (layer < 3)
          gemmTiled<kPairHidden, kPairHidden, kPairHidden, tokenBlockFor(kPairHidden, kPairTile), columnBlockFor(kPairHidden, kPairTile)>(
              acc, dense.weight, dense.bias, 0, source, kPairHidden, destination, kPairHidden, rows);
        else
          gemmTiled<kPairHidden, kHeads, kHeads, tokenBlockFor(kHeads, kPairTile), columnBlockFor(kHeads, kPairTile)>(
              acc, dense.weight, dense.bias, 0, source, kPairHidden, destination, kPairHidden, rows);

        int const width = layer < 3 ? kPairHidden : kHeads;
        for (int index = lane; index < rows * width; index += lanes) {
          int const local = index / width;
          int const c = index - local * width;
          float value = destination[local * kPairHidden + c];
          value = (value - norm.mean[c]) / alpaka::math::sqrt(acc, norm.var[c] + norm.eps) * norm.gamma[c] +
                  norm.beta[c];
          if (layer < 3 || generated::pair_final_gelu)
            value = geluOf(acc, value);
          destination[local * kPairHidden + c] = value;
        }
        alpaka::syncBlockThreads(acc);

        float* const swap = source;
        source = destination;
        destination = swap;
      }

      for (int index = lane; index < raw * kHeads; index += lanes) {
        int const local = index / kHeads;
        int const head = index - local * kHeads;
        int const slot = begin + local;
        int i = 0;
        int j = 0;
        pairFromSlot(acc, slot, size, i, j);
        bool const valid = shared.mask[i] > .5f && shared.mask[j] > .5f;
        shared.pair[head * kMaxPairs + slot] = valid ? source[local * kPairHidden + head] : 0.f;
      }
      alpaka::syncBlockThreads(acc);
    }
  }

  // ---------------------------------------------------------------- particle block

  // Token feed-forward.  Only kFfnTokenTile hidden vectors of width kFfnDim
  // are live at a time, which is what keeps the whole graph inside 44 kB of
  // shared memory.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void feedForward(TAcc const& acc,
                                        generated::Norm const& preFc,
                                        generated::Linear const& fc1,
                                        generated::Norm const& postFc,
                                        generated::Linear const& fc2,
                                        float const* residualScale,
                                        JetShared& shared,
                                        int tokens) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    float* const hidden = shared.big;
    float* const projected = shared.big + kFfnTokenTile * kFfnDim;

    for (int begin = 0; begin < tokens; begin += kFfnTokenTile) {
      int const rows = tokens - begin < kFfnTokenTile ? tokens - begin : kFfnTokenTile;
      // `na` and `hidden` are token-contiguous: [channel][kFfnTokenTile].  The
      // padding lanes of a partial tile are stale but finite, and GELU is
      // applied to the whole tile because the valid entries are no longer a
      // contiguous prefix.
      layerNorm<kEmbedDim, false, true>(
          acc, preFc, shared.x + begin * kEmbedDim, kEmbedDim, shared.na, kFfnTokenTile, rows, shared);
      gemmTiled<kEmbedDim, kFfnDim, kFfnDim, tokenBlockFor(kFfnDim, kFfnTokenTile), columnBlockFor(kFfnDim, kFfnTokenTile), true, true>(
          acc, fc1.weight, fc1.bias, 0, shared.na, kFfnTokenTile, hidden, kFfnTokenTile, rows);
      applyGelu(acc, hidden, kFfnTokenTile * kFfnDim);
      layerNorm<kFfnDim, true, true>(acc, postFc, hidden, kFfnTokenTile, hidden, kFfnTokenTile, rows, shared);
      gemmTiled<kFfnDim, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kFfnTokenTile), columnBlockFor(kEmbedDim, kFfnTokenTile), true>(
          acc, fc2.weight, fc2.bias, 0, hidden, kFfnTokenTile, projected, kEmbedDim, rows);

      for (int index = lane; index < rows * kEmbedDim; index += lanes) {
        int const channel = index % kEmbedDim;
        int const target = begin * kEmbedDim + index;
        shared.x[target] = residualScale[channel] * shared.x[target] + projected[index];
      }
      alpaka::syncBlockThreads(acc);
    }
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void particleBlock(TAcc const& acc,
                                          generated::BlockData const& parameters,
                                          JetShared& shared) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const tokens = shared.tokens;
    int const active = shared.active;

    // Q, K and V of every token live in `big` with a stride of kQkvDim.  The
    // pre-attention LayerNorm is tiled so that only kFfnTokenTile normalised
    // tokens need a buffer of their own.
    for (int begin = 0; begin < tokens; begin += kFfnTokenTile) {
      int const rows = tokens - begin < kFfnTokenTile ? tokens - begin : kFfnTokenTile;
      layerNorm<kEmbedDim, false, true>(
          acc, parameters.preAttn, shared.x + begin * kEmbedDim, kEmbedDim, shared.na, kFfnTokenTile, rows, shared);
      gemmTiled<kEmbedDim, kQkvDim, kQkvDim, tokenBlockFor(kQkvDim, kFfnTokenTile), columnBlockFor(kQkvDim, kFfnTokenTile), true>(acc,
                                                                       parameters.qkv.weight,
                                                                       parameters.qkv.bias,
                                                                       0,
                                                                       shared.na,
                                                                       kFfnTokenTile,
                                                                       shared.big + begin * kQkvDim,
                                                                       kQkvDim,
                                                                       rows);
    }

    // One thread owns one (query, head): it reads its own slice of Q and then
    // overwrites that same slice with the attention output, so no extra buffer
    // is needed.  Scores stay in registers because the key loop has the
    // compile-time bound kMaxParticles.
    for (int task = lane; task < tokens * kHeads; task += lanes) {
      int const query = task / kHeads;
      int const head = task - query * kHeads;
      float* const slice = shared.big + query * kQkvDim + head * kHeadDim;

      float score[kMaxParticles];
      float maximum = -3.4e38f;
      for (int key = 0; key < kMaxParticles; ++key) {
        float const* const k = shared.big + key * kQkvDim + kEmbedDim + head * kHeadDim;
        float value = 0.f;
        for (int d = 0; d < kHeadDim; ++d)
          value += slice[d] * k[d];
        value *= .25f;
        if (query < active && key < active)
          value += shared.pair[head * kMaxPairs + pairSlot(query, key, active)];
        if (shared.mask[key] <= .5f)
          value = -3.4e38f;
        score[key] = value;
        if (value > maximum)
          maximum = value;
      }
      float denominator = 0.f;
      for (int key = 0; key < kMaxParticles; ++key) {
        score[key] = alpaka::math::exp(acc, score[key] - maximum);
        denominator += score[key];
      }
      float const inverse = 1.f / denominator;
      for (int key = 0; key < kMaxParticles; ++key)
        score[key] *= inverse;

      float merged[kHeadDim];
      for (int d = 0; d < kHeadDim; ++d) {
        float value = 0.f;
        for (int key = 0; key < kMaxParticles; ++key)
          value += score[key] * shared.big[key * kQkvDim + 2 * kEmbedDim + head * kHeadDim + d];
        merged[d] = value;
      }
      for (int d = 0; d < kHeadDim; ++d)
        slice[d] = merged[d];
    }
    alpaka::syncBlockThreads(acc);

    // The rest of the attention branch reuses the dead K and V columns of the
    // same token rows: [0,128) holds the heads, [128,256) the projection and
    // [256,384) the scaled and normalised result.
    gemmTiled<kEmbedDim, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kMaxParticles), columnBlockFor(kEmbedDim, kMaxParticles)>(acc,
                                                                         parameters.out.weight,
                                                                         parameters.out.bias,
                                                                         0,
                                                                         shared.big,
                                                                         kQkvDim,
                                                                         shared.big + kEmbedDim,
                                                                         kQkvDim,
                                                                         tokens);

    for (int index = lane; index < tokens * kEmbedDim; index += lanes) {
      int const token = index / kEmbedDim;
      int const channel = index - token * kEmbedDim;
      int const head = channel / kHeadDim;
      int const d = channel - head * kHeadDim;
      int const target = generated::transpose_scaled_heads ? d * kHeads + head : channel;
      shared.big[token * kQkvDim + 2 * kEmbedDim + target] =
          shared.big[token * kQkvDim + kEmbedDim + channel] * parameters.headScale[head];
    }
    alpaka::syncBlockThreads(acc);

    layerNorm<kEmbedDim>(acc,
                         parameters.postAttn,
                         shared.big + 2 * kEmbedDim,
                         kQkvDim,
                         shared.big,
                         kQkvDim,
                         tokens,
                         shared);
    for (int index = lane; index < tokens * kEmbedDim; index += lanes) {
      int const token = index / kEmbedDim;
      int const channel = index - token * kEmbedDim;
      shared.x[index] += shared.big[token * kQkvDim + channel];
    }
    alpaka::syncBlockThreads(acc);

    feedForward(acc, parameters.preFc, parameters.fc1, parameters.postFc, parameters.fc2,
                parameters.residualScale, shared, tokens);
  }

  // ---------------------------------------------------------------- class block

  template <typename TAcc>
  ALPAKA_FN_ACC inline void classBlock(TAcc const& acc,
                                       generated::BlockData const& parameters,
                                       JetShared& shared) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const tokens = shared.tokens;
    int const keys = tokens + 1;

    constexpr int kvStride = 2 * kEmbedDim;
    float* const kv = shared.big;
    float* const query = shared.big + (kMaxParticles + 1) * kvStride;
    float* const scores = query + kEmbedDim;
    float* const merged = scores + kHeads * (kMaxParticles + 1);
    float* const temporary = merged + kEmbedDim;

    layerNorm<kEmbedDim>(acc, parameters.preAttn, shared.cls, kEmbedDim, shared.clsNorm, kEmbedDim, 1, shared);

    // Keys and values come from the normalised [class, particles] sequence.
    gemmTiled<kEmbedDim, kQkvDim, kvStride, 1, columnBlockFor(kvStride, 1)>(
        acc, parameters.qkv.weight, parameters.qkv.bias, kEmbedDim, shared.clsNorm, kEmbedDim, kv, kvStride, 1);
    for (int begin = 0; begin < tokens; begin += kFfnTokenTile) {
      int const rows = tokens - begin < kFfnTokenTile ? tokens - begin : kFfnTokenTile;
      layerNorm<kEmbedDim, false, true>(
          acc, parameters.preAttn, shared.x + begin * kEmbedDim, kEmbedDim, shared.na, kFfnTokenTile, rows, shared);
      gemmTiled<kEmbedDim, kQkvDim, kvStride, tokenBlockFor(kvStride, kFfnTokenTile), columnBlockFor(kvStride, kFfnTokenTile), true>(acc,
                                                                        parameters.qkv.weight,
                                                                        parameters.qkv.bias,
                                                                        kEmbedDim,
                                                                        shared.na,
                                                                        kFfnTokenTile,
                                                                        kv + (begin + 1) * kvStride,
                                                                        kvStride,
                                                                        rows);
    }
    // Weaver projects the class query from the unnormalised class token.
    gemmTiled<kEmbedDim, kQkvDim, kEmbedDim, 1, columnBlockFor(kEmbedDim, 1)>(
        acc, parameters.qkv.weight, parameters.qkv.bias, 0, shared.cls, kEmbedDim, query, kEmbedDim, 1);

    for (int task = lane; task < kHeads * keys; task += lanes) {
      int const head = task / keys;
      int const key = task - head * keys;
      float const* const k = kv + key * kvStride + head * kHeadDim;
      float value = 0.f;
      for (int d = 0; d < kHeadDim; ++d)
        value += query[head * kHeadDim + d] * k[d];
      value *= .25f;
      if (key > 0 && shared.mask[key - 1] <= .5f)
        value = -3.4e38f;
      scores[head * keys + key] = value;
    }
    alpaka::syncBlockThreads(acc);

    for (int head = lane; head < kHeads; head += lanes) {
      float maximum = -3.4e38f;
      for (int key = 0; key < keys; ++key)
        if (scores[head * keys + key] > maximum)
          maximum = scores[head * keys + key];
      float denominator = 0.f;
      for (int key = 0; key < keys; ++key) {
        float const value = alpaka::math::exp(acc, scores[head * keys + key] - maximum);
        scores[head * keys + key] = value;
        denominator += value;
      }
      float const inverse = 1.f / denominator;
      for (int key = 0; key < keys; ++key)
        scores[head * keys + key] *= inverse;
    }
    alpaka::syncBlockThreads(acc);

    for (int channel = lane; channel < kEmbedDim; channel += lanes) {
      int const head = channel / kHeadDim;
      int const d = channel - head * kHeadDim;
      float value = 0.f;
      for (int key = 0; key < keys; ++key)
        value += scores[head * keys + key] * kv[key * kvStride + kEmbedDim + head * kHeadDim + d];
      merged[channel] = value;
    }
    alpaka::syncBlockThreads(acc);

    gemmTiled<kEmbedDim, kEmbedDim, kEmbedDim, 1, columnBlockFor(kEmbedDim, 1)>(
        acc, parameters.out.weight, parameters.out.bias, 0, merged, kEmbedDim, temporary, kEmbedDim, 1);
    for (int channel = lane; channel < kEmbedDim; channel += lanes) {
      int const head = channel / kHeadDim;
      int const d = channel - head * kHeadDim;
      int const target = generated::transpose_scaled_heads ? d * kHeads + head : channel;
      merged[target] = temporary[channel] * parameters.headScale[head];
    }
    alpaka::syncBlockThreads(acc);

    layerNorm<kEmbedDim>(acc, parameters.postAttn, merged, kEmbedDim, temporary, kEmbedDim, 1, shared);
    for (int channel = lane; channel < kEmbedDim; channel += lanes)
      shared.cls[channel] += temporary[channel];
    alpaka::syncBlockThreads(acc);

    layerNorm<kEmbedDim>(acc, parameters.preFc, shared.cls, kEmbedDim, temporary, kEmbedDim, 1, shared);
    gemmTiled<kEmbedDim, kFfnDim, kFfnDim, 1, columnBlockFor(kFfnDim, 1)>(
        acc, parameters.fc1.weight, parameters.fc1.bias, 0, temporary, kEmbedDim, shared.big, kFfnDim, 1);
    applyGelu(acc, shared.big, kFfnDim);
    layerNorm<kFfnDim>(acc, parameters.postFc, shared.big, kFfnDim, shared.big, kFfnDim, 1, shared);
    gemmTiled<kFfnDim, kEmbedDim, kEmbedDim, 1, columnBlockFor(kEmbedDim, 1)>(
        acc, parameters.fc2.weight, parameters.fc2.bias, 0, shared.big, kFfnDim, merged, kEmbedDim, 1);
    for (int channel = lane; channel < kEmbedDim; channel += lanes)
      shared.cls[channel] = parameters.residualScale[channel] * shared.cls[channel] + merged[channel];
    alpaka::syncBlockThreads(acc);
  }

  // ---------------------------------------------------------------- full graph

  // Expects shared.na[t * kInputFeatures + f], shared.p4 and shared.mask to be
  // filled and a barrier to have been executed.  Leaves the two class
  // probabilities in shared.logits.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void runParticleTransformer(TAcc const& acc,
                                                   generated::ModelView model,
                                                   JetShared& shared) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);

    if (lane == 0) {
      int last = 0;
      for (int t = 0; t < kMaxParticles; ++t)
        if (shared.mask[t] > .5f)
          last = t + 1;
      shared.active = last;
      int padded = ((last + kTokenBlock - 1) / kTokenBlock) * kTokenBlock;
      shared.tokens = padded < kMaxParticles ? padded : kMaxParticles;
      shared.pairs = last * (last + 1) / 2;
    }
    // Stale rows are never read by an unmasked token, but they must stay
    // finite so that a zero softmax weight cannot produce a NaN, and so that a
    // trailing partial GEMM tile cannot multiply an uninitialised value.  The
    // head of `na` already holds the input features written by the caller.
    for (int index = lane + kMaxParticles * kInputFeatures; index < kNormScratch; index += lanes)
      shared.na[index] = 0.f;
    fill(acc, shared.x, kMaxParticles * kEmbedDim, 0.f);
    fill(acc, shared.big, kBigScratch, 0.f);
    fill(acc, shared.pair, kHeads * kMaxPairs, 0.f);
    alpaka::syncBlockThreads(acc);

    int const tokens = shared.tokens;
    if (tokens == 0) {
      if (lane == 0) {
        shared.logits[0] = 0.f;
        shared.logits[1] = 0.f;
      }
      alpaka::syncBlockThreads(acc);
      return;
    }

    // ---- input embedding: BN -> LN(15) -> 15x128 -> LN -> 128x512 -> LN -> 512x128
    applyBatchNorm<kInputFeatures>(acc, generated::inputBN(model), shared.na, shared.big, tokens);
    // The whole sequence is normalised at once here, so `na` is
    // [feature][kMaxParticles] rather than [channel][kFfnTokenTile].
    layerNorm<kInputFeatures, false, true>(
        acc, generated::embedNorm(model, 0), shared.big, kInputFeatures, shared.na, kMaxParticles, tokens, shared);
    {
      auto const dense = generated::embedLinear(model, 0);
      gemmTiled<kInputFeatures, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kMaxParticles), columnBlockFor(kEmbedDim, kMaxParticles), true>(
          acc, dense.weight, dense.bias, 0, shared.na, kMaxParticles, shared.big, kEmbedDim, tokens);
    }
    applyGelu(acc, shared.big, tokens * kEmbedDim);
    {
      auto const first = generated::embedLinear(model, 1);
      auto const second = generated::embedLinear(model, 2);
      auto const firstNorm = generated::embedNorm(model, 1);
      auto const secondNorm = generated::embedNorm(model, 2);
      // big[0, 16*128) keeps the 128-wide embedding of every token while the
      // 512-wide hidden state of one tile occupies the rest of the buffer.
      float* const hidden = shared.big + kMaxParticles * kEmbedDim;
      for (int begin = 0; begin < tokens; begin += kFfnTokenTile) {
        int const rows = tokens - begin < kFfnTokenTile ? tokens - begin : kFfnTokenTile;
        layerNorm<kEmbedDim, false, true>(acc,
                             firstNorm,
                             shared.big + begin * kEmbedDim,
                             kEmbedDim,
                             shared.na,
                             kFfnTokenTile,
                             rows,
                             shared);
        gemmTiled<kEmbedDim, kFfnDim, kFfnDim, tokenBlockFor(kFfnDim, kFfnTokenTile), columnBlockFor(kFfnDim, kFfnTokenTile), true, true>(
            acc, first.weight, first.bias, 0, shared.na, kFfnTokenTile, hidden, kFfnTokenTile, rows);
        applyGelu(acc, hidden, kFfnTokenTile * kFfnDim);
        layerNorm<kFfnDim, true, true>(acc, secondNorm, hidden, kFfnTokenTile, hidden, kFfnTokenTile, rows, shared);
        gemmTiled<kFfnDim, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kFfnTokenTile), columnBlockFor(kEmbedDim, kFfnTokenTile), true>(
            acc, second.weight, second.bias, 0, hidden, kFfnTokenTile, shared.na, kEmbedDim, rows);
        for (int index = lane; index < rows * kEmbedDim; index += lanes) {
          int const token = begin + index / kEmbedDim;
          shared.x[begin * kEmbedDim + index] = shared.mask[token] > .5f ? geluOf(acc, shared.na[index]) : 0.f;
        }
        alpaka::syncBlockThreads(acc);
      }
    }

    encodePairs(acc, model, shared);

    for (int layer = 0; layer < kParticleBlocks; ++layer)
      particleBlock(acc, generated::block(model, layer), shared);

    auto const classToken = generated::clsToken(model);
    for (int channel = lane; channel < kEmbedDim; channel += lanes)
      shared.cls[channel] = classToken[channel];
    alpaka::syncBlockThreads(acc);

    for (int layer = 0; layer < kClassBlocks; ++layer)
      classBlock(acc, generated::clsBlock(model, layer), shared);

    layerNorm<kEmbedDim>(
        acc, generated::finalNorm(model), shared.cls, kEmbedDim, shared.clsNorm, kEmbedDim, 1, shared);
    {
      auto const head = generated::classifier(model);
      gemmNarrow<kEmbedDim, kClasses>(acc, head.weight, head.bias, shared.clsNorm, shared.logits);
    }
    if (lane == 0) {
      float const maximum = shared.logits[0] > shared.logits[1] ? shared.logits[0] : shared.logits[1];
      float const e0 = alpaka::math::exp(acc, shared.logits[0] - maximum);
      float const e1 = alpaka::math::exp(acc, shared.logits[1] - maximum);
      shared.logits[0] = e0 / (e0 + e1);
      shared.logits[1] = e1 / (e0 + e1);
    }
    alpaka::syncBlockThreads(acc);
  }

}  // namespace partfp32

#endif
