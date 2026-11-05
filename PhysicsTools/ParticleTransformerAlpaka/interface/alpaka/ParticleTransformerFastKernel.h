#ifndef PhysicsTools_ParticleTransformerAlpaka_ParticleTransformerFastKernel_h
#define PhysicsTools_ParticleTransformerAlpaka_ParticleTransformerFastKernel_h

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

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelDefinition.h"

namespace part {

  // ---------------------------------------------------------------- shared state

  struct alignas(16) JetShared {
    float x[kMaxParticles * kEmbedDim];         // token states, live for the whole graph
    alignas(16) float na[kNormScratch];         // one tile of normalised tokens
    alignas(16) std::int8_t quant[kQuantScratch];  // one tile of INT8 activations, row-major
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
  // 47872 B: the INT8 staging area costs 4 kB on top of the FP32 layout.  Still
  // inside the 48 kB CUDA static limit and the 47 KiB Alpaka CPU arena, and it
  // does not change residency on Ada (100 kB/SM still gives two blocks).
  static_assert(sizeof(JetShared) <= 47u * 1024u, "the per-jet shared state grew past its budget");

  // A 16-byte aggregate load.  Reading four adjacent activations in one
  // instruction is what makes the token-contiguous layout worth having: it
  // turns four shared-memory instructions into one, and instruction issue -
  // not shared-memory bandwidth - is what limits the inner loop.
  struct alignas(16) Float4 {
    float v[4];
  };

  // Sixteen INT8 activations, or four packed weight words, in one 16-byte load.
  struct alignas(16) Int4 {
    int v[4];
  };

  // Four INT8 multiply-accumulates in one instruction.  This is the whole point
  // of quantizing activations: the FP32 pipe retires one multiply-add per lane
  // per clock, dp4a retires four.
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
  // element once, converts it to float once, and reuses it for TB tokens.  CB
  // then trades registers for broadcast reads of the activation, which come
  // from shared memory and are nearly free.
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
  //         + scale[columnBegin + n] * sum_k input[t][k] * q[k][columnBegin + n]
  //
  // The weights are symmetric per-output-channel INT8, so the scale is a
  // property of the column and comes out of the k sum entirely: the inner loop
  // accumulates against the raw integer values and the scale is applied once
  // per output, instead of once per multiply-add as the scalar reference does.
  // That is the same quantity in exact arithmetic and slightly more accurate in
  // FP32, because the reference rounds `scale * q` before every product.
  //
  // KDIM, WSTRIDE, NOUT, TB and CB are compile-time so the inner loops unroll
  // and the TB x CB accumulator tile stays in registers.  Occupancy here is
  // limited by shared memory rather than by registers, so a large tile is
  // free.  Threads are laid out column-fastest, which keeps each warp's weight
  // fetch a single coalesced transaction.  A large TB also amortises the
  // integer-to-float conversion of each weight over TB tokens, which matters
  // more here than in the FP32 build because the conversion competes with the
  // multiply-add pipeline.
  // out[t][columnBegin + n] = bias[.] + s_x * s_w[.] * sum_k q(x[t][k]) * q(w[.][k])
  //
  // The activations arrive as INT8 in row-major [token][KDIM] and the weights
  // in the [k/4][out][4] image, so one 16-byte load fetches sixteen inputs of a
  // token and each 32-bit weight word feeds one dp4a.  The accumulator is
  // INT32: with KDIM at most 512 the largest partial sum is 512 * 127 * 127,
  // about 8.3e6, so it cannot overflow.
  //
  // Both scales are applied once per output, exactly as the weight-only path
  // applied the weight scale, which makes this arithmetically identical to what
  // quantization-aware training simulated.
  template <int KDIM, int WSTRIDE, int NOUT, int TB, int CB, typename TAcc>
  ALPAKA_FN_ACC inline void gemmInt8(TAcc const& acc,
                                     std::int8_t const* __restrict__ packed,
                                     float const* __restrict__ scale,
                                     float const* __restrict__ bias,
                                     float activationScale,
                                     int columnBegin,
                                     std::int8_t const* __restrict__ input,
                                     float* __restrict__ output,
                                     int outputStride,
                                     int tokens) {
    static_assert(KDIM % 16 == 0, "the input width must allow 16-byte activation loads");
    static_assert(NOUT % CB == 0, "the column block must divide the output width");
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    constexpr int groups = NOUT / CB;
    constexpr int words = KDIM / 4;  // packed 32-bit words per token
    int const tiles = (tokens + TB - 1) / TB;
    int const tasks = tiles * groups;

    auto const* __restrict__ weightWords = reinterpret_cast<int const*>(packed);
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
          int weight[4];
          for (int r = 0; r < 4; ++r)
            weight[r] = weightWords[(group + r) * WSTRIDE + columnBegin + column + c];
          for (int t = 0; t < TB; ++t)
            for (int r = 0; r < 4; ++r)
              sum[t][c] = dot4(acc, activation[t].v[r], weight[r], sum[t][c]);
        }
      }

      for (int t = 0; t < TB; ++t)
        if (token + t < tokens)
          for (int c = 0; c < CB; ++c) {
            int const channel = columnBegin + column + c;
            output[(token + t) * outputStride + column + c] =
                static_cast<float>(sum[t][c]) * (activationScale * scale[channel]) + bias[channel];
          }
    }
    alpaka::syncBlockThreads(acc);
  }

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
                                      std::int8_t const* __restrict__ weightT,
                                      float const* __restrict__ scale,
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

      std::int8_t const* __restrict__ weight = weightT + columnBegin + column;
      float sum[TB][CB];
      for (int t = 0; t < TB; ++t)
        for (int c = 0; c < CB; ++c)
          sum[t][c] = 0.f;

      float const* __restrict__ row = input + (INPUT_T ? token : token * inputStride);
      for (int k = 0; k < KDIM; ++k) {
        float value[CB];
        for (int c = 0; c < CB; ++c)
          value[c] = static_cast<float>(weight[k * WSTRIDE + c]);

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
            int const channel = columnBegin + column + c;
            float const result = sum[t][c] * scale[channel] + bias[channel];
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
                                       float const* __restrict__ weight,
                                       float const* __restrict__ bias,
                                       float const* __restrict__ input,
                                       float* __restrict__ output) {
    // The classifier is 256 weights directly under the softmax and was left in
    // FP32 by quantization-aware training, so it is evaluated in FP32 here
    // too.  Its weights are row-major [out][in]; with two outputs there is
    // nothing to gain from a transposed layout.
    for (int column = laneOf(acc); column < NOUT; column += laneCount(acc)) {
      float value = bias[column];
      for (int k = 0; k < KDIM; ++k)
        value += input[k] * weight[column * KDIM + k];
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

  // Round a block of activations onto the INT8 grid without changing its type.
  // Used for the one layer whose input width is not a multiple of 16: it keeps
  // the FP32 GEMM, but its activation still has to be quantized so the result
  // matches what quantization-aware training simulated.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void roundToGrid(TAcc const& acc, float* __restrict__ data, int count, float scale) {
    if (scale <= 0.f)
      return;
    float const inverse = 1.f / scale;
    for (int index = laneOf(acc); index < count; index += laneCount(acc))
      data[index] = static_cast<float>(quantizeOne(data[index], inverse)) * scale;
    alpaka::syncBlockThreads(acc);
  }

  // Quantize a block of activations into the INT8 staging area.  Every integer
  // GEMM is fed by this: a LayerNorm writes FP32 as before, then one pass turns
  // the tile into INT8.  Fusing the two would save a pass over at most 4 kB and
  // would mean duplicating the group reduction, so it is kept separate.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void quantizeRows(TAcc const& acc,
                                         float const* __restrict__ input,
                                         int inputStride,
                                         std::int8_t* __restrict__ output,
                                         int width,
                                         int rows,
                                         float inverseScale) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    for (int index = lane; index < rows * width; index += lanes) {
      int const row = index / width;
      int const c = index - row * width;
      output[row * width + c] = quantizeOne(input[row * inputStride + c], inverseScale);
    }
    alpaka::syncBlockThreads(acc);
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
        // The first layer takes four raw log-scale kinematics and keeps the
        // FP32 path; the other three are 64 wide and go through dp4a.  The
        // source rows are strided by kPairHidden even when only kHeads columns
        // are live, so the quantized copy is packed down to its true width.
        if (layer == 0) {
          gemmTiled<kPairInputs, kPairHidden, kPairHidden, tokenBlockFor(kPairHidden, kPairTile), columnBlockFor(kPairHidden, kPairTile)>(
              acc, dense.weight, dense.scale, dense.bias, 0, source, kPairHidden, destination, kPairHidden, rows);
        } else {
          quantizeRows(
              acc, source, kPairHidden, shared.quant, kPairHidden, rows, 1.f / dense.activation[0]);
          if (layer < 3)
            gemmInt8<kPairHidden, kPairHidden, kPairHidden, tokenBlockFor(kPairHidden, kPairTile), columnBlockFor(kPairHidden, kPairTile)>(
                acc, dense.weight, dense.scale, dense.bias, dense.activation[0], 0, shared.quant, destination, kPairHidden, rows);
          else
            gemmInt8<kPairHidden, kHeads, kHeads, tokenBlockFor(kHeads, kPairTile), columnBlockFor(kHeads, kPairTile)>(
                acc, dense.weight, dense.scale, dense.bias, dense.activation[0], 0, shared.quant, destination, kPairHidden, rows);
        }

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
      layerNorm<kEmbedDim>(
          acc, preFc, shared.x + begin * kEmbedDim, kEmbedDim, shared.na, kEmbedDim, rows, shared);
      quantizeRows(acc, shared.na, kEmbedDim, shared.quant, kEmbedDim, rows, 1.f / fc1.activation[0]);
      gemmInt8<kEmbedDim, kFfnDim, kFfnDim, tokenBlockFor(kFfnDim, kFfnTokenTile), columnBlockFor(kFfnDim, kFfnTokenTile)>(
          acc, fc1.weight, fc1.scale, fc1.bias, fc1.activation[0], 0, shared.quant, hidden, kFfnDim, rows);
      applyGelu(acc, hidden, rows * kFfnDim);
      layerNorm<kFfnDim>(acc, postFc, hidden, kFfnDim, hidden, kFfnDim, rows, shared);
      quantizeRows(acc, hidden, kFfnDim, shared.quant, kFfnDim, rows, 1.f / fc2.activation[0]);
      gemmInt8<kFfnDim, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kFfnTokenTile), columnBlockFor(kEmbedDim, kFfnTokenTile)>(
          acc, fc2.weight, fc2.scale, fc2.bias, fc2.activation[0], 0, shared.quant, projected, kEmbedDim, rows);

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
      layerNorm<kEmbedDim>(
          acc, parameters.preAttn, shared.x + begin * kEmbedDim, kEmbedDim, shared.na, kEmbedDim, rows, shared);
      quantizeRows(
          acc, shared.na, kEmbedDim, shared.quant, kEmbedDim, rows, 1.f / parameters.qkv.activation[0]);
      gemmInt8<kEmbedDim, kQkvDim, kQkvDim, tokenBlockFor(kQkvDim, kFfnTokenTile), columnBlockFor(kQkvDim, kFfnTokenTile)>(acc,
                                                                       parameters.qkv.weight,
                                                                       parameters.qkv.scale,
                                                                       parameters.qkv.bias,
                                                                       parameters.qkv.activation[0],
                                                                       0,
                                                                       shared.quant,
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
    // The only integer GEMM not fed by a normalization: its input is the
    // attention output sitting in the Q columns of `big`.
    quantizeRows(
        acc, shared.big, kQkvDim, shared.quant, kEmbedDim, tokens, 1.f / parameters.out.activation[0]);
    gemmInt8<kEmbedDim, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kMaxParticles), columnBlockFor(kEmbedDim, kMaxParticles)>(acc,
                                                                         parameters.out.weight,
                                                                         parameters.out.scale,
                                                                         parameters.out.bias,
                                                                         parameters.out.activation[0],
                                                                         0,
                                                                         shared.quant,
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
    quantizeRows(acc, shared.clsNorm, kEmbedDim, shared.quant, kEmbedDim, 1, 1.f / parameters.qkv.activation[0]);
    gemmInt8<kEmbedDim, kQkvDim, kvStride, 1, columnBlockFor(kvStride, 1)>(
        acc, parameters.qkv.weight, parameters.qkv.scale, parameters.qkv.bias, parameters.qkv.activation[0],
        kEmbedDim, shared.quant, kv, kvStride, 1);
    for (int begin = 0; begin < tokens; begin += kFfnTokenTile) {
      int const rows = tokens - begin < kFfnTokenTile ? tokens - begin : kFfnTokenTile;
      layerNorm<kEmbedDim>(
          acc, parameters.preAttn, shared.x + begin * kEmbedDim, kEmbedDim, shared.na, kEmbedDim, rows, shared);
      quantizeRows(
          acc, shared.na, kEmbedDim, shared.quant, kEmbedDim, rows, 1.f / parameters.qkv.activation[0]);
      gemmInt8<kEmbedDim, kQkvDim, kvStride, tokenBlockFor(kvStride, kFfnTokenTile), columnBlockFor(kvStride, kFfnTokenTile)>(acc,
                                                                        parameters.qkv.weight,
                                                                        parameters.qkv.scale,
                                                                        parameters.qkv.bias,
                                                                        parameters.qkv.activation[0],
                                                                        kEmbedDim,
                                                                        shared.quant,
                                                                        kv + (begin + 1) * kvStride,
                                                                        kvStride,
                                                                        rows);
    }
    // Weaver projects the class query from the unnormalised class token.
    // The query has its own activation scale: it observes the raw class token
    // rather than the normalised sequence.
    quantizeRows(acc, shared.cls, kEmbedDim, shared.quant, kEmbedDim, 1, 1.f / parameters.qkv.activation[1]);
    gemmInt8<kEmbedDim, kQkvDim, kEmbedDim, 1, columnBlockFor(kEmbedDim, 1)>(
        acc, parameters.qkv.weight, parameters.qkv.scale, parameters.qkv.bias, parameters.qkv.activation[1],
        0, shared.quant, query, kEmbedDim, 1);

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

    quantizeRows(acc, merged, kEmbedDim, shared.quant, kEmbedDim, 1, 1.f / parameters.out.activation[0]);
    gemmInt8<kEmbedDim, kEmbedDim, kEmbedDim, 1, columnBlockFor(kEmbedDim, 1)>(
        acc, parameters.out.weight, parameters.out.scale, parameters.out.bias, parameters.out.activation[0],
        0, shared.quant, temporary, kEmbedDim, 1);
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
    quantizeRows(acc, temporary, kEmbedDim, shared.quant, kEmbedDim, 1, 1.f / parameters.fc1.activation[0]);
    gemmInt8<kEmbedDim, kFfnDim, kFfnDim, 1, columnBlockFor(kFfnDim, 1)>(
        acc, parameters.fc1.weight, parameters.fc1.scale, parameters.fc1.bias, parameters.fc1.activation[0],
        0, shared.quant, shared.big, kFfnDim, 1);
    applyGelu(acc, shared.big, kFfnDim);
    layerNorm<kFfnDim>(acc, parameters.postFc, shared.big, kFfnDim, shared.big, kFfnDim, 1, shared);
    quantizeRows(acc, shared.big, kFfnDim, shared.quant, kFfnDim, 1, 1.f / parameters.fc2.activation[0]);
    gemmInt8<kFfnDim, kEmbedDim, kEmbedDim, 1, columnBlockFor(kEmbedDim, 1)>(
        acc, parameters.fc2.weight, parameters.fc2.scale, parameters.fc2.bias, parameters.fc2.activation[0],
        0, shared.quant, merged, kEmbedDim, 1);
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
    // 15 inputs is not a multiple of 16, so this one projection keeps the FP32
    // path; it is 0.1% of the arithmetic.
    layerNorm<kInputFeatures>(
        acc, generated::embedNorm(model, 0), shared.big, kInputFeatures, shared.na, kInputFeatures, tokens, shared);
    {
      auto const dense = generated::embedLinear(model, 0);
      roundToGrid(acc, shared.na, tokens * kInputFeatures, dense.activation[0]);
      gemmTiled<kInputFeatures, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kMaxParticles), columnBlockFor(kEmbedDim, kMaxParticles)>(
          acc, dense.weight, dense.scale, dense.bias, 0, shared.na, kInputFeatures, shared.big, kEmbedDim, tokens);
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
        layerNorm<kEmbedDim>(acc,
                             firstNorm,
                             shared.big + begin * kEmbedDim,
                             kEmbedDim,
                             shared.na,
                             kEmbedDim,
                             rows,
                             shared);
        quantizeRows(acc, shared.na, kEmbedDim, shared.quant, kEmbedDim, rows, 1.f / first.activation[0]);
        gemmInt8<kEmbedDim, kFfnDim, kFfnDim, tokenBlockFor(kFfnDim, kFfnTokenTile), columnBlockFor(kFfnDim, kFfnTokenTile)>(
            acc, first.weight, first.scale, first.bias, first.activation[0], 0, shared.quant, hidden, kFfnDim, rows);
        applyGelu(acc, hidden, rows * kFfnDim);
        layerNorm<kFfnDim>(acc, secondNorm, hidden, kFfnDim, hidden, kFfnDim, rows, shared);
        quantizeRows(acc, hidden, kFfnDim, shared.quant, kFfnDim, rows, 1.f / second.activation[0]);
        gemmInt8<kFfnDim, kEmbedDim, kEmbedDim, tokenBlockFor(kEmbedDim, kFfnTokenTile), columnBlockFor(kEmbedDim, kFfnTokenTile)>(
            acc, second.weight, second.scale, second.bias, second.activation[0], 0, shared.quant, shared.na, kEmbedDim, rows);
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

}  // namespace part

#endif
