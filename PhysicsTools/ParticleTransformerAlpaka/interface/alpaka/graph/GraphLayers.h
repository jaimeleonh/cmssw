#ifndef PhysicsTools_ParticleTransformerAlpaka_graph_GraphLayers_h
#define PhysicsTools_ParticleTransformerAlpaka_graph_GraphLayers_h

// Block-cooperative implementations of the layers the interpreter dispatches.
//
// These are the same algorithms as the hand-written ParticleTransformer kernel,
// with the shapes moved from template parameters to arguments.  What stays at
// compile time is the part that has to: the register tile of a GEMM.  Each GEMM
// picks its tile from the runtime shapes and jumps into one of twelve
// instantiations, so the accumulator still lives in registers and the innermost
// loops still unroll, while the layer widths come from the model file.
//
// Every routine here is called by all lanes of the block and ends with a
// barrier, so the interpreter can treat an op as a single step.

#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/GraphShared.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelFormat.h"

namespace nn {

  // ---------------------------------------------------------------- lanes

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneOf(TAcc const& acc) {
    return static_cast<int>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneCount(TAcc const& acc) {
    return static_cast<int>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
  }

  // A slice of a block that evaluates one jet on its own -- in warp mode, one
  // warp.  It stands in for the accelerator wherever the interpreter asks
  // which lane it is, how many lanes there are, or waits for the others, so
  // the same op code runs at block or at warp granularity.
  //
  // `lanes` must be a whole warp (kWarpLanes) and the group warp-aligned: the
  // barrier below is a full-mask __syncwarp.
  template <typename TAcc>
  struct SubGroup {
    TAcc const& acc;
    int lane;
    int lanes;
  };

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneOf(SubGroup<TAcc> const& group) {
    return group.lane;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline int laneCount(SubGroup<TAcc> const& group) {
    return group.lanes;
  }

  // The accelerator the Alpaka math functions need, which a SubGroup wraps.
  template <typename TAcc>
  ALPAKA_FN_ACC inline TAcc const& inner(TAcc const& acc) {
    return acc;
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline TAcc const& inner(SubGroup<TAcc> const& group) {
    return group.acc;
  }

  // Wait for the other lanes of whatever is evaluating this jet.
  //
  // For a SubGroup this has to be a warp-level sync and not a block barrier,
  // and not only for speed.  Each warp is a different jet, so their control
  // flow differs -- an op with no rows returns before its barrier, a tiled
  // region runs as many times as that jet has tiles -- and a block barrier
  // would deadlock the moment two jets in a block disagreed.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void barrier(TAcc const& acc) {
    alpaka::syncBlockThreads(acc);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void barrier(SubGroup<TAcc> const& group) {
#if defined(__CUDA_ARCH__)
    __syncwarp();
#else
    // Only correct when the group is the whole block: the CPU test harness,
    // which runs one group at a time.  Warp mode is launched on CUDA only (see
    // GraphInferenceAlgo), so a multi-group block never reaches this line.
    alpaka::syncBlockThreads(group.acc);
#endif
  }

  // A 16-byte aggregate load: four packed weight words, or sixteen INT8
  // activations, in one instruction.
  struct alignas(16) Int4 {
    int v[4];
  };

  // Four INT8 multiply-accumulates in one instruction.  This is the whole
  // point of quantizing activations: the FP32 pipe retires one multiply-add
  // per lane per clock, dp4a retires four.
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

  // ---------------------------------------------------------------- lattice

  // Walking a [rows][columns] rectangle with one lane per element, without an
  // integer division per step.
  //
  // `index / columns` looks free and is not.  An NVIDIA SM has no integer
  // divide unit, so a division by a value the compiler cannot see expands to
  // roughly twenty instructions -- a reciprocal estimate, two multiplies and a
  // fix-up.  The hand-written kernel never pays that because every width is a
  // template parameter; the interpreter reads its widths from the file, so
  // every elementwise op was paying it once per element, and the GEMMs once
  // per register tile.  Measured on an L4, that alone was most of the gap
  // between the two.
  //
  // Stepping the row and column forward by a precomputed stride costs an add
  // and a compare instead.  The two divisions in the constructor are amortised
  // over the whole loop.
  struct Lattice {
    // `at` rather than `index`: a template of that name is visible from the
    // CUDA headers, and `walk.index < n` inside a template function then
    // parses as the start of a template argument list.
    int at;
    int row;
    int column;
    int columns;
    int indexStep;
    int rowStep;
    int columnStep;

    ALPAKA_FN_ACC Lattice(int lane, int lanes, int width) : at(lane), indexStep(lanes) {
      // The width comes from the model file.  ModelLoader rejects a buffer
      // narrower than one column, but an op's own width is not covered by
      // that, and a zero here would be a division by zero on the device.
      columns = width > 0 ? width : 1;
      row = lane / columns;
      column = lane - row * columns;
      rowStep = lanes / columns;
      columnStep = lanes - rowStep * columns;
    }

    ALPAKA_FN_ACC void advance() {
      at += indexStep;
      row += rowStep;
      column += columnStep;
      if (column >= columns) {
        column -= columns;
        ++row;
      }
    }
  };

  // ---------------------------------------------------------------- activations

  template <typename TAcc>
  ALPAKA_FN_ACC inline float applyActivation(TAcc const& acc, int kind, float value) {
    switch (static_cast<Activation>(kind)) {
      case Activation::Gelu:
        // The exact erf form, which is what training used.  The tanh
        // approximation differs by up to 1e-3 and would show up directly in
        // the output.
        return .5f * value * (1.f + alpaka::math::erf(inner(acc), value * .7071067811865475f));
      case Activation::Relu:
        return value > 0.f ? value : 0.f;
      case Activation::Silu:
        return value / (1.f + alpaka::math::exp(inner(acc), -value));
      case Activation::Tanh:
        return alpaka::math::tanh(inner(acc), value);
      case Activation::Sigmoid:
        return 1.f / (1.f + alpaka::math::exp(inner(acc), -value));
      case Activation::Identity:
      default:
        return value;
    }
  }

  // Symmetric INT8 quantization with round-half-to-even, matching torch.round
  // so the kernel reproduces what quantization-aware training simulated.
  ALPAKA_FN_HOST_ACC inline std::int8_t quantizeOne(float value, float inverseScale) {
    float const scaled = value * inverseScale;
    float rounded = static_cast<float>(static_cast<int>(scaled >= 0.f ? scaled + .5f : scaled - .5f));
    float const difference = scaled - rounded;
    if (difference == .5f || difference == -.5f) {
      if (static_cast<int>(rounded) % 2 != 0)
        rounded -= (difference > 0.f ? -1.f : 1.f);
    }
    if (rounded > 127.f)
      rounded = 127.f;
    if (rounded < -127.f)
      rounded = -127.f;
    return static_cast<std::int8_t>(rounded);
  }

  // ---------------------------------------------------------------- tiling

  // The task space of a GEMM is (row tiles) x (column groups) and the tile is
  // the largest one that still gives every lane work.  A larger row block is
  // preferred because each weight element is loaded once, converted once, and
  // then reused for that many rows.
  ALPAKA_FN_ACC inline bool tileCovers(int outputs, int rows, int tb, int cb, int lanes) {
    return outputs % cb == 0 && ((rows + tb - 1) / tb) * (outputs / cb) >= lanes;
  }

  ALPAKA_FN_ACC inline void chooseTile(int outputs, int rows, int lanes, int& tb, int& cb) {
    for (int candidateTb = 8; candidateTb > 1; candidateTb /= 2)
      for (int candidateCb = 4; candidateCb >= 1; candidateCb /= 2)
        if (candidateTb <= rows && tileCovers(outputs, rows, candidateTb, candidateCb, lanes)) {
          tb = candidateTb;
          cb = candidateCb;
          return;
        }
    tb = 1;
    cb = 1;
    for (int candidateCb = 4; candidateCb > 1; candidateCb /= 2)
      if (tileCovers(outputs, rows, 1, candidateCb, lanes)) {
        cb = candidateCb;
        return;
      }
  }

  // ---------------------------------------------------------------- int8 GEMM

  // out[r][column + n] = bias[.] + activationScale * scale[.]
  //                        * sum_k q(input[r][k]) * q(weight[.][k])
  //
  // Weights are already in the [k/4][out][4] image the exporter wrote, so one
  // 32-bit word feeds one dp4a and a warp's fetch of consecutive output columns
  // is one coalesced transaction.  The accumulator is INT32: with an input
  // width of at most 512 the largest partial sum is 512 * 127 * 127, about
  // 8.3e6, so it cannot overflow.
  //
  // Both scales are applied once per output rather than once per multiply-add,
  // which is the same quantity in exact arithmetic and slightly more accurate
  // in FP32.
  template <int TB, int CB, typename TAcc>
  ALPAKA_FN_ACC inline void gemmInt8Tile(TAcc const& acc,
                                         std::int8_t const* __restrict__ packed,
                                         float const* __restrict__ scale,
                                         float const* __restrict__ bias,
                                         float activationScale,
                                         int weightOut,
                                         int columnBegin,
                                         int outputs,
                                         int kdim,
                                         std::int8_t const* __restrict__ input,
                                         int inputStride,
                                         float* __restrict__ output,
                                         int outputStride,
                                         int rows) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const groups = outputs / CB;
    int const words = kdim / 4;  // packed 32-bit words per row
    int const inputWordsPerRow = inputStride / 4;
    int const tiles = (rows + TB - 1) / TB;
    int const tasks = tiles * groups;

    auto const* __restrict__ weightWords = reinterpret_cast<int const*>(packed);
    auto const* __restrict__ inputWords = reinterpret_cast<int const*>(input);
    // Sixteen activations per load when the widths allow it, which turns four
    // shared-memory instructions into one.  Instruction issue, not bandwidth,
    // is what limits this loop.
    bool const vectorised = (words % 4 == 0) && (inputWordsPerRow % 4 == 0);

    for (Lattice walk(lane, lanes, groups); walk.at < tasks; walk.advance()) {
      int const column = walk.column * CB;
      int const row = walk.row * TB;

      int sum[TB][CB];
      for (int t = 0; t < TB; ++t)
        for (int c = 0; c < CB; ++c)
          sum[t][c] = 0;

      // Both strides are runtime values here, where the hand-written kernel
      // has them as template parameters.  Recomputing `k * stride + column`
      // per access is a multiply-add per weight word -- as many instructions
      // as the dp4a they feed.  Walking pointers forward instead turns each
      // one into an add, which is what makes the interpreter's inner loop cost
      // the same as the compiled kernel's.
      int const* __restrict__ weightRow = weightWords + columnBegin + column;
      int const* __restrict__ activationRow = inputWords + row * inputWordsPerRow;

      if (vectorised) {
        for (int group = 0; group < words; group += 4) {
          Int4 activation[TB];
          int const* source = activationRow + group;
          for (int t = 0; t < TB; ++t) {
            activation[t] = *reinterpret_cast<Int4 const*>(source);
            source += inputWordsPerRow;
          }
          for (int q = 0; q < 4; ++q) {
            for (int c = 0; c < CB; ++c) {
              int const w = weightRow[c];
              for (int t = 0; t < TB; ++t)
                sum[t][c] = dot4(acc, activation[t].v[q], w, sum[t][c]);
            }
            weightRow += weightOut;
          }
        }
      } else {
        for (int group = 0; group < words; ++group) {
          int activation[TB];
          int const* source = activationRow + group;
          for (int t = 0; t < TB; ++t) {
            activation[t] = *source;
            source += inputWordsPerRow;
          }
          for (int c = 0; c < CB; ++c) {
            int const w = weightRow[c];
            for (int t = 0; t < TB; ++t)
              sum[t][c] = dot4(acc, activation[t], w, sum[t][c]);
          }
          weightRow += weightOut;
        }
      }

      // A partial trailing tile computes rows that do not exist; they are
      // simply not written back.
      for (int t = 0; t < TB; ++t)
        if (row + t < rows)
          for (int c = 0; c < CB; ++c) {
            int const o = columnBegin + column + c;
            output[(row + t) * outputStride + column + c] =
                static_cast<float>(sum[t][c]) * (activationScale * scale[o]) + bias[o];
          }
    }
    barrier(acc);
  }

  // ---------------------------------------------------------------- float GEMM

  // The same shape, for a layer whose weights are INT8 but whose input width
  // is not a multiple of 16.  Weights are in the [in][out] image; the
  // per-channel scale still comes out of the k sum.
  template <int TB, int CB, typename TAcc>
  ALPAKA_FN_ACC inline void gemmFloatTile(TAcc const& acc,
                                          std::int8_t const* __restrict__ weight,
                                          float const* __restrict__ scale,
                                          float const* __restrict__ bias,
                                          int weightOut,
                                          int columnBegin,
                                          int outputs,
                                          int kdim,
                                          float const* __restrict__ input,
                                          int inputStride,
                                          float* __restrict__ output,
                                          int outputStride,
                                          int rows) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const groups = outputs / CB;
    int const tiles = (rows + TB - 1) / TB;
    int const tasks = tiles * groups;

    for (Lattice walk(lane, lanes, groups); walk.at < tasks; walk.advance()) {
      int const column = walk.column * CB;
      int const row = walk.row * TB;

      float sum[TB][CB];
      for (int t = 0; t < TB; ++t)
        for (int c = 0; c < CB; ++c)
          sum[t][c] = 0.f;

      // Pointers walked forward rather than `k * stride` recomputed, for the
      // same reason as the integer GEMM above: both strides come from the file.
      std::int8_t const* __restrict__ weightRow = weight + columnBegin + column;
      float const* __restrict__ activationColumn = input + row * inputStride;
      for (int k = 0; k < kdim; ++k) {
        float value[CB];
        for (int c = 0; c < CB; ++c)
          value[c] = static_cast<float>(weightRow[c]);
        float element[TB];
        float const* source = activationColumn + k;
        for (int t = 0; t < TB; ++t) {
          element[t] = *source;
          source += inputStride;
        }
        for (int t = 0; t < TB; ++t)
          for (int c = 0; c < CB; ++c)
            sum[t][c] += element[t] * value[c];
        weightRow += weightOut;
      }

      for (int t = 0; t < TB; ++t)
        if (row + t < rows)
          for (int c = 0; c < CB; ++c) {
            int const o = columnBegin + column + c;
            output[(row + t) * outputStride + column + c] = sum[t][c] * scale[o] + bias[o];
          }
    }
    barrier(acc);
  }

  // ---------------------------------------------------------------- dispatch

  // Jump from the runtime tile choice into the matching instantiation.  Twelve
  // of them covers every shape a model has produced so far, and each one keeps
  // its accumulator in registers.
#define NN_TILE_DISPATCH(CALL)                         \
  do {                                                 \
    int tb = 1, cb = 1;                                \
    chooseTile(outputs, rows, laneCount(acc), tb, cb); \
    switch (tb * 8 + cb) {                             \
      case 8 * 8 + 4:                                  \
        CALL(8, 4);                                    \
        break;                                         \
      case 8 * 8 + 2:                                  \
        CALL(8, 2);                                    \
        break;                                         \
      case 8 * 8 + 1:                                  \
        CALL(8, 1);                                    \
        break;                                         \
      case 4 * 8 + 4:                                  \
        CALL(4, 4);                                    \
        break;                                         \
      case 4 * 8 + 2:                                  \
        CALL(4, 2);                                    \
        break;                                         \
      case 4 * 8 + 1:                                  \
        CALL(4, 1);                                    \
        break;                                         \
      case 2 * 8 + 4:                                  \
        CALL(2, 4);                                    \
        break;                                         \
      case 2 * 8 + 2:                                  \
        CALL(2, 2);                                    \
        break;                                         \
      case 2 * 8 + 1:                                  \
        CALL(2, 1);                                    \
        break;                                         \
      case 1 * 8 + 4:                                  \
        CALL(1, 4);                                    \
        break;                                         \
      case 1 * 8 + 2:                                  \
        CALL(1, 2);                                    \
        break;                                         \
      default:                                         \
        CALL(1, 1);                                    \
        break;                                         \
    }                                                  \
  } while (false)

  template <typename TAcc>
  ALPAKA_FN_ACC inline void gemmInt8(TAcc const& acc,
                                     std::int8_t const* packed,
                                     float const* scale,
                                     float const* bias,
                                     float activationScale,
                                     int weightOut,
                                     int columnBegin,
                                     int outputs,
                                     int kdim,
                                     std::int8_t const* input,
                                     int inputStride,
                                     float* output,
                                     int outputStride,
                                     int rows) {
#define NN_CALL(TB, CB)                                                                   \
  gemmInt8Tile<TB, CB>(acc, packed, scale, bias, activationScale, weightOut, columnBegin, \
                       outputs, kdim, input, inputStride, output, outputStride, rows)
    NN_TILE_DISPATCH(NN_CALL);
#undef NN_CALL
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void gemmFloat(TAcc const& acc,
                                      std::int8_t const* weight,
                                      float const* scale,
                                      float const* bias,
                                      int weightOut,
                                      int columnBegin,
                                      int outputs,
                                      int kdim,
                                      float const* input,
                                      int inputStride,
                                      float* output,
                                      int outputStride,
                                      int rows) {
#define NN_CALL(TB, CB)                                                                  \
  gemmFloatTile<TB, CB>(acc, weight, scale, bias, weightOut, columnBegin, outputs, kdim, \
                        input, inputStride, output, outputStride, rows)
    NN_TILE_DISPATCH(NN_CALL);
#undef NN_CALL
  }

  // A layer whose weights stay FP32, row-major [out][in] -- every layer of a
  // `--float` TorchScript export, and ParT's classifier.
  //
  // This used to compute one output per lane with a global weight load and a
  // shared input load for every multiply-add.  For ParT that was a handful of
  // tiny layers and did not matter; for a small FP32 network it is nearly all
  // the arithmetic, and in warp mode it was what the SM stalled on (memory
  // instruction throttling, once the barrier stall was gone).  The register
  // tile is the same idea as the INT8 GEMMs above: a weight is loaded once and
  // used for TB rows, an input element once and used for CB columns.
  //
  // Each accumulator starts at its bias and walks k in order, which is the
  // arithmetic the one-output-per-lane version did, so only the compiler's
  // FMA contraction can move the result.
  //
  // Out of line on purpose.  The interpreter is one very large kernel and
  // nvcc's inliner works to a budget for it: inlining these twelve tiles
  // pushed it over, and it responded by outlining ~57 *other* functions --
  // putting real calls and register spills into ParT's hot path, which never
  // runs this code, and costing it ~3%.  Kept out of line, each tile is one
  // call per Linear op, which is nothing next to the GEMM it performs, and
  // everything else inlines as before.  Measured, not assumed: see the
  // warp-mode notes in GRAPH_FORMAT.md.
  //
  // `__attribute__((noinline))` rather than `__noinline__`: it is what CUDA's
  // `__noinline__` expands to, and gcc and clang accept it as well, so the CPU
  // back-end and the test harness build too.
  template <int TB, int CB, typename TAcc>
  ALPAKA_FN_ACC __attribute__((noinline)) void gemmPlainTile(TAcc const& acc,
                                          float const* __restrict__ weight,
                                          float const* __restrict__ bias,
                                          int columnBegin,
                                          int outputs,
                                          int kdim,
                                          float const* __restrict__ input,
                                          int inputStride,
                                          float* __restrict__ output,
                                          int outputStride,
                                          int rows) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const groups = outputs / CB;
    int const tiles = (rows + TB - 1) / TB;
    int const tasks = tiles * groups;

    for (Lattice walk(lane, lanes, groups); walk.at < tasks; walk.advance()) {
      int const column = walk.column * CB;
      int const row = walk.row * TB;

      float sum[TB][CB];
      for (int c = 0; c < CB; ++c)
        for (int t = 0; t < TB; ++t)
          sum[t][c] = bias[columnBegin + column + c];

      // One base pointer each and small integer offsets, not a pointer per
      // row and per column: twelve 64-bit pointers were enough to make this
      // spill, and the spill cost the ParT models ~3% although they barely
      // use this routine -- the stack frame is the whole kernel's.
      float const* __restrict__ weightBase = weight + (columnBegin + column) * kdim;
      float const* __restrict__ inputBase = input + row * inputStride;

      // A trailing partial tile reads its missing rows from the last real
      // one rather than past the buffer; they are computed and never stored.
      int const last = rows - 1 - row;
      int inputOffset[TB];
      for (int t = 0; t < TB; ++t)
        inputOffset[t] = (t < last ? t : last) * inputStride;

      for (int k = 0; k < kdim; ++k) {
        float w[CB];
        for (int c = 0; c < CB; ++c)
          w[c] = weightBase[c * kdim + k];
        float x[TB];
        for (int t = 0; t < TB; ++t)
          x[t] = inputBase[inputOffset[t] + k];
        for (int t = 0; t < TB; ++t)
          for (int c = 0; c < CB; ++c)
            sum[t][c] += x[t] * w[c];
      }

      for (int t = 0; t < TB; ++t)
        if (row + t < rows)
          for (int c = 0; c < CB; ++c)
            output[(row + t) * outputStride + column + c] = sum[t][c];
    }
    barrier(acc);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void gemmPlain(TAcc const& acc,
                                      float const* weight,
                                      float const* bias,
                                      int columnBegin,
                                      int outputs,
                                      int kdim,
                                      float const* input,
                                      int inputStride,
                                      float* output,
                                      int outputStride,
                                      int rows) {
#define NN_CALL(TB, CB)                                                                        \
  gemmPlainTile<TB, CB>(acc, weight, bias, columnBegin, outputs, kdim, input, inputStride, \
                        output, outputStride, rows)
    NN_TILE_DISPATCH(NN_CALL);
#undef NN_CALL
  }

  // ---------------------------------------------------------------- normalization

  // Two-pass LayerNorm.  A group of `group` lanes reduces one row; the two
  // partial arrays are distinct so no barrier is needed between reading the
  // mean partials and writing the variance partials.
  template <typename TAcc, typename TShared>
  ALPAKA_FN_ACC inline void layerNorm(TAcc const& acc,
                                      float const* __restrict__ gamma,
                                      float const* __restrict__ beta,
                                      float eps,
                                      float const* __restrict__ input,
                                      int inputStride,
                                      float* __restrict__ output,
                                      int outputStride,
                                      int size,
                                      int rows,
                                      TShared& shared,
                                      std::uint16_t flags = 0) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    int const reducing = lanes < TShared::kReduce ? lanes : TShared::kReduce;

    // A power of two, so the lane's position within its group is a mask and a
    // shift rather than a division by a value the compiler cannot see.
    int group = 1;
    int groupShift = 0;
    while (group < 32 && group * 2 * rows <= reducing) {
      group *= 2;
      ++groupShift;
    }
    int const groups = reducing >> groupShift;

    for (int base = 0; base < rows; base += groups) {
      int const slot = lane & (group - 1);
      int const row = base + (lane >> groupShift);
      bool const contributes = lane < reducing && row < rows;

      float partial = 0.f;
      if (contributes)
        for (int c = slot; c < size; c += group)
          partial += input[row * inputStride + c];
      if (lane < reducing)
        shared.reduceA[lane] = partial;
      barrier(acc);

      float mean = 0.f;
      float variance = 0.f;
      if (contributes) {
        int const first = (lane >> groupShift) << groupShift;
        for (int k = 0; k < group; ++k)
          mean += shared.reduceA[first + k];
        mean /= static_cast<float>(size);
        for (int c = slot; c < size; c += group) {
          float const difference = input[row * inputStride + c] - mean;
          variance += difference * difference;
        }
      }
      if (lane < reducing)
        shared.reduceB[lane] = variance;
      barrier(acc);

      if (contributes && slot == 0) {
        int const first = (lane >> groupShift) << groupShift;
        float total = 0.f;
        for (int k = 0; k < group; ++k)
          total += shared.reduceB[first + k];
        // torch.std() is unbiased by default, and a model written with
        // mean()/std() rather than nn.LayerNorm inherits that; over 128
        // channels the two denominators differ by 0.4%, which is far too much
        // to paper over.
        int const denominator = (flags & kFlagNormUnbiased) && size > 1 ? size - 1 : size;
        total /= static_cast<float>(denominator);
        shared.mean[row - base] = mean;
        // Likewise, `(x - mean) / (std + eps)` is not `(x - mean) /
        // sqrt(var + eps)`; which one a model meant is read from its graph.
        shared.rstd[row - base] = (flags & kFlagNormEpsOutside)
                                      ? 1.f / (alpaka::math::sqrt(inner(acc), total) + eps)
                                      : 1.f / alpaka::math::sqrt(inner(acc), total + eps);
      }
      barrier(acc);

      int const covered = (rows - base < groups ? rows - base : groups) * size;
      for (Lattice walk(lane, lanes, size); walk.at < covered; walk.advance()) {
        int const local = walk.row;
        int const c = walk.column;
        int const r = base + local;
        output[r * outputStride + c] =
            (input[r * inputStride + c] - shared.mean[local]) * shared.rstd[local] * gamma[c] +
            beta[c];
      }
      barrier(acc);
    }
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void batchNorm(TAcc const& acc,
                                      float const* __restrict__ gamma,
                                      float const* __restrict__ beta,
                                      float const* __restrict__ mean,
                                      float const* __restrict__ variance,
                                      float eps,
                                      float const* __restrict__ input,
                                      int inputStride,
                                      float* __restrict__ output,
                                      int outputStride,
                                      int size,
                                      int rows) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    for (Lattice walk(lane, lanes, size); walk.at < rows * size; walk.advance()) {
      int const row = walk.row;
      int const c = walk.column;
      output[row * outputStride + c] =
          (input[row * inputStride + c] - mean[c]) / alpaka::math::sqrt(inner(acc), variance[c] + eps) *
              gamma[c] +
          beta[c];
    }
    barrier(acc);
  }

  // ---------------------------------------------------------------- pair helpers

  ALPAKA_FN_ACC inline int triangularBase(int row, int size) {
    return row * size - (row * (row - 1)) / 2;
  }

  ALPAKA_FN_ACC inline int pairSlot(int i, int j, int size) {
    int const low = i < j ? i : j;
    int const high = i < j ? j : i;
    return triangularBase(low, size) + (high - low);
  }

  template <typename TAcc>
  ALPAKA_FN_ACC inline void pairFromSlot(TAcc const& acc, int slot, int size, int& i, int& j) {
    float const b = 2.f * static_cast<float>(size) + 1.f;
    int row =
        static_cast<int>((b - alpaka::math::sqrt(inner(acc), b * b - 8.f * static_cast<float>(slot))) * .5f);
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

  // The four pairwise Lorentz features of ParticleTransformer's interaction
  // embedding.  Every step is FP32 and in this order on purpose: for the
  // diagonal pair the invariant mass is a cancellation of two nearly equal
  // squares, so what survives is rounding, and computing it any other way puts
  // the 1e-8 clamp on the other side and changes the attention bias visibly.
  template <typename TAcc>
  ALPAKA_FN_ACC inline void lorentzFeatures(
      TAcc const& acc, float const* a, float const* b, float* out, int count) {
    float const pxi = a[0], pyi = a[1], pzi = a[2], ei = a[3];
    float const pxj = b[0], pyj = b[1], pzj = b[2], ej = b[3];
    float const pti = alpaka::math::sqrt(inner(acc), pxi * pxi + pyi * pyi);
    float const ptj = alpaka::math::sqrt(inner(acc), pxj * pxj + pyj * pyj);
    float const yi =
        .5f * alpaka::math::log(inner(acc), 1.f + 2.f * pzi / alpaka::math::max(inner(acc), ei - pzi, 1e-20f));
    float const yj =
        .5f * alpaka::math::log(inner(acc), 1.f + 2.f * pzj / alpaka::math::max(inner(acc), ej - pzj, 1e-20f));
    float dphi = alpaka::math::atan2(inner(acc), pyi, pxi) - alpaka::math::atan2(inner(acc), pyj, pxj);
    constexpr float pi = 3.14159265358979323846f;
    constexpr float twoPi = 2.f * pi;
    while (dphi > pi)
      dphi -= twoPi;
    while (dphi <= -pi)
      dphi += twoPi;
    float const distance = alpaka::math::sqrt(inner(acc), (yi - yj) * (yi - yj) + dphi * dphi);
    float const minimum = alpaka::math::min(inner(acc), pti, ptj);
    float const sx = pxi + pxj, sy = pyi + pyj, sz = pzi + pzj, se = ei + ej;

    float value[4];
    value[0] = alpaka::math::log(inner(acc), alpaka::math::max(inner(acc), minimum * distance, 1e-8f));
    value[1] = alpaka::math::log(inner(acc), alpaka::math::max(inner(acc), minimum / alpaka::math::max(inner(acc), pti + ptj, 1e-8f), 1e-8f));
    value[2] = alpaka::math::log(inner(acc), alpaka::math::max(inner(acc), distance, 1e-8f));
    value[3] =
        alpaka::math::log(inner(acc), alpaka::math::max(inner(acc), se * se - sx * sx - sy * sy - sz * sz, 1e-8f));
    for (int c = 0; c < count && c < 4; ++c)
      out[c] = value[c];
  }

}  // namespace nn

#endif
