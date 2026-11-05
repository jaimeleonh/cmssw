#ifndef PhysicsTools_ParticleTransformerAlpaka_graph_GraphKernel_h
#define PhysicsTools_ParticleTransformerAlpaka_graph_GraphKernel_h

// The interpreter.  A group of lanes evaluates one jet by walking the op table
// of a model file; every activation lives in shared memory, so the only global
// traffic is the streaming read of the weights, which all resident blocks
// share through L1 and L2.
//
// The group is either a whole 256-lane block (block layout, the full arena --
// ParT, tau) or a single warp (warp layout, eight jets per block with a small
// arena each -- a small network such as a DeepSets).  The op code below does
// not know which: it asks laneOf / laneCount / barrier / inner, and an
// nn::SubGroup answers those for a warp.  It is templated on the state type
// (GraphShared or GraphWarpState) for the same reason.
//
// The whole point is that this file contains no model.  It knows how to
// evaluate a Linear, a LayerNorm, an Attention and a dozen other ops, and it
// knows how to walk a tiled region; which of them to run, over what shapes,
// and where to put the results are all read from the file.  A new network is a
// new file.
//
// Control flow is uniform across the block: every lane walks the same op list
// and the same tile loop, so the tile base is a register rather than shared
// state and no barrier is needed to agree on it.  Each layer routine ends with
// a barrier of its own, so an op is a single step as far as this file cares.

#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/GraphLayers.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/GraphShared.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelFormat.h"

namespace nn {

  // Where one of the model's inputs lives in the arena.  A caller that would
  // rather write its inputs straight into shared memory than stage a copy --
  // which is what the producer does, because a staging buffer would cost the
  // arena memory it needs -- uses this to find the destination.
  struct InputSlot {
    int offset = -1;
    int stride = 0;
    int rows = 0;
    int cols = 0;
    ALPAKA_FN_HOST_ACC bool valid() const { return offset >= 0; }
  };

  // Everything the kernel needs about a loaded model.  All pointers are into
  // one device upload; nothing here is owned.
  struct ModelView {
    OpRecord const* ops;
    BufferRecord const* buffers;
    std::int8_t const* weights;
    float const* params;
    int opCount;
    int maxTokens;
    int outputSize;
    InputSlot features;
    InputSlot vectors;
  };

  // Per-jet inputs.  A null pointer means the caller has already written that
  // input straight into its arena buffer -- which is what the producer does,
  // because staging it a second time would cost shared memory the arena needs.
  // ModelLoader reports where each input buffer lives.
  struct GraphInputs {
    float const* features;
    float const* vectors;
  };

  // ---------------------------------------------------------------- addressing

  template <typename TShared>
  ALPAKA_FN_ACC inline int resolveRows(int spec, TShared const& shared) {
    if (spec >= 0)
      return spec;
    switch (spec) {
      case kRowsTokens:
        return shared.tokens;
      case kRowsTokensPlusOne:
        return shared.tokens + 1;
      case kRowsPairs:
        return shared.pairs;
      case kRowsActive:
        return shared.active;
      case kRowsOne:
      default:
        return 1;
    }
  }

  // Element (0, 0) of a buffer in the frame of the enclosing tiled region.  A
  // tiled buffer holds only the tile and is indexed from zero; anything else
  // spans the whole sequence and is indexed from the tile base.
  ALPAKA_FN_ACC inline int bufferOrigin(BufferRecord const buffer, int tileBase) {
    int const base = buffer.tiled ? 0 : tileBase;
    return buffer.layout == kRowMajor ? buffer.offset + base * buffer.stride
                                      : buffer.offset + base;
  }

  // ---------------------------------------------------------------- attention

  // Fused multi-head attention.  One lane owns one (query, head): it reads its
  // own slice of Q, keeps the whole score row in registers, and then writes the
  // result back over that same slice.  Writing in place is what lets a
  // transformer block reuse its Q columns for the attention output instead of
  // finding another 8 kB of shared memory -- which means the output buffer may
  // alias the queries, but must never alias the keys or the values, since
  // other lanes are still reading those.
  template <int MAXKEYS, typename TAcc, typename TShared>
  ALPAKA_FN_ACC inline void attentionImpl(TAcc const& acc,
                                          float const* __restrict__ queries,
                                          int queryStride,
                                          float const* __restrict__ keys,
                                          int keyStride,
                                          float const* __restrict__ values,
                                          int valueStride,
                                          float* __restrict__ output,
                                          int outputStride,
                                          float const* __restrict__ bias,
                                          int biasStride,
                                          int rows,
                                          int keyRows,
                                          int heads,
                                          int headDim,
                                          float scale,
                                          bool maskKeys,
                                          int keyOffset,
                                          TShared const& shared) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);

    for (Lattice walk(lane, lanes, heads); walk.at < rows * heads; walk.advance()) {
      int const query = walk.row;
      int const head = walk.column;
      int const channel = head * headDim;

      float score[MAXKEYS];
      float maximum = -3.4e38f;
      for (int key = 0; key < keyRows; ++key) {
        float value = 0.f;
        for (int d = 0; d < headDim; ++d)
          value += queries[query * queryStride + channel + d] * keys[key * keyStride + channel + d];
        value *= scale;
        if (bias != nullptr && query < shared.active && key < shared.active)
          value += bias[head * biasStride + pairSlot(query, key, shared.active)];
        if (maskKeys) {
          int const particle = key - keyOffset;
          // Masked keys are forced to -3.4e38, so their weight underflows to
          // exactly zero: skipping them is an identity, not an approximation.
          //
          // A key past the last unmasked particle is padding by construction,
          // and has to be masked too.  A model run at its own multiplicity
          // never produces one, but a fixed-length model does: it evaluates
          // the full padded sequence, and the rows beyond the jet are exactly
          // the ones the original module's key_padding_mask removed.
          if (particle >= 0 && (particle >= shared.tokens || shared.mask[particle] <= .5f))
            value = -3.4e38f;
        }
        score[key] = value;
        if (value > maximum)
          maximum = value;
      }

      float denominator = 0.f;
      for (int key = 0; key < keyRows; ++key) {
        score[key] = alpaka::math::exp(inner(acc), score[key] - maximum);
        denominator += score[key];
      }
      float const inverse = 1.f / denominator;

      // Q is finished with at this point, so each output element can be
      // written as it is computed even when `output` aliases `queries`.
      for (int d = 0; d < headDim; ++d) {
        float value = 0.f;
        for (int key = 0; key < keyRows; ++key)
          value += score[key] * values[key * valueStride + channel + d];
        output[query * outputStride + channel + d] = value * inverse;
      }
    }
    barrier(acc);
  }

  // ---------------------------------------------------------------- one op

  // The op is taken by value, not by reference.  A reference leaves every
  // field access as its own load from the op table, and the table lives in
  // global memory: the interpreter was issuing about a dozen scalar loads per
  // op per warp just to find out what to do.  A 64-byte `alignas(16)` copy is
  // four vector loads, and every lane of the block wants the same bytes, so
  // they broadcast.  The same reasoning applies to the BufferRecords below.
  template <typename TAcc, typename TShared>
  ALPAKA_FN_ACC inline void executeOp(TAcc const& acc,
                                      ModelView const& model,
                                      OpRecord const op,
                                      TShared& shared,
                                      GraphInputs const& inputs,
                                      int tileBase,
                                      int tileRows) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);
    float* const arena = shared.arena;
    std::int8_t* const staging = shared.quant;
    float const* const params = model.params;

    // Inside a region a literal row count is the tile size, and the last tile
    // may be shorter.
    int const rows = (tileRows >= 0 && op.rows >= 0)
                         ? (op.rows < tileRows ? op.rows : tileRows)
                         : resolveRows(op.rows, shared);
    int const cols = op.cols;
    if (rows <= 0)
      return;

    BufferRecord const* const buffers = model.buffers;
    auto origin = [&](int id) { return bufferOrigin(buffers[id], tileBase); };

    switch (static_cast<OpCode>(op.opcode)) {
      case OpCode::Nop:
      case OpCode::TileBegin:
      case OpCode::TileEnd:
        return;

      case OpCode::LoadInput: {
        float const* source = op.i[0] == 0 ? inputs.features : inputs.vectors;
        if (source == nullptr)
          return;  // already staged in the arena by the caller
        BufferRecord const out = buffers[op.out];
        int const base = origin(op.out);
        for (Lattice walk(lane, lanes, cols); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row;
          int const c = walk.column;
          arena[base + row * out.stride + c] = source[row * cols + c];
        }
        barrier(acc);
        return;
      }

      case OpCode::LoadConst: {
        int const base = origin(op.out);
        for (int c = lane; c < op.i[0]; c += lanes)
          arena[base + c] = params[op.t[0] + c];
        barrier(acc);
        return;
      }

      case OpCode::Linear: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        int const kdim = src.cols;
        int const columnBegin = op.i[0];
        int const weightOut = op.i[2];
        float const* const bias = params + op.t[2];
        float* const destination = arena + origin(op.out);

        if (op.flags & kFlagInt8Input) {
          gemmInt8(acc, model.weights + op.t[0], params + op.t[1], bias, params[op.t[3]],
                   weightOut, columnBegin, cols, kdim, staging + origin(op.in[0]), src.stride,
                   destination, dst.stride, rows);
        } else if (op.flags & kFlagInt8Weights) {
          gemmFloat(acc, model.weights + op.t[0], params + op.t[1], bias, weightOut, columnBegin,
                    cols, kdim, arena + origin(op.in[0]), src.stride, destination, dst.stride,
                    rows);
        } else {
          gemmPlain(acc, params + op.t[0], bias, columnBegin, cols, kdim,
                    arena + origin(op.in[0]), src.stride, destination, dst.stride, rows);
        }
        return;
      }

      case OpCode::LayerNorm: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        layerNorm(acc, params + op.t[0], params + op.t[1], op.f[0], arena + origin(op.in[0]),
                  src.stride, arena + origin(op.out), dst.stride, cols, rows, shared, op.flags);
        return;
      }

      case OpCode::BatchNorm: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        batchNorm(acc, params + op.t[0], params + op.t[1], params + op.t[2], params + op.t[3],
                  op.f[0], arena + origin(op.in[0]), src.stride, arena + origin(op.out),
                  dst.stride, cols, rows);
        return;
      }

      case OpCode::Activate: {
        BufferRecord const dst = buffers[op.out];
        int const base = origin(op.out);
        for (Lattice walk(lane, lanes, cols); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row;
          int const c = walk.column;
          int const at = base + row * dst.stride + c;
          arena[at] = applyActivation(acc, op.i[0], arena[at]);
        }
        barrier(acc);
        return;
      }

      case OpCode::Residual: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        int const sourceBase = origin(op.in[0]);
        int const targetBase = origin(op.out);
        float const* const alpha = op.t[0] >= 0 ? params + op.t[0] : nullptr;
        for (Lattice walk(lane, lanes, cols); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row;
          int const c = walk.column;
          int const at = targetBase + row * dst.stride + c;
          float const addend = arena[sourceBase + row * src.stride + c];
          arena[at] = alpha != nullptr ? alpha[c] * arena[at] + addend : arena[at] + addend;
        }
        barrier(acc);
        return;
      }

      case OpCode::Quantize: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        int const sourceBase = origin(op.in[0]);
        int const targetBase = origin(op.out);
        float const inverse = 1.f / params[op.t[0]];
        for (Lattice walk(lane, lanes, cols); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row;
          int const c = walk.column;
          staging[targetBase + row * dst.stride + c] =
              quantizeOne(arena[sourceBase + row * src.stride + c], inverse);
        }
        barrier(acc);
        return;
      }

      case OpCode::RoundToGrid: {
        BufferRecord const dst = buffers[op.out];
        int const base = origin(op.out);
        float const scale = params[op.t[0]];
        if (scale <= 0.f)
          return;
        float const inverse = 1.f / scale;
        for (Lattice walk(lane, lanes, cols); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row;
          int const c = walk.column;
          int const at = base + row * dst.stride + c;
          arena[at] = static_cast<float>(quantizeOne(arena[at], inverse)) * scale;
        }
        barrier(acc);
        return;
      }

      case OpCode::Attention: {
        BufferRecord const keys = buffers[op.in[0]];
        BufferRecord const values = buffers[op.in[1]];
        BufferRecord const queries = buffers[op.in[2]];
        BufferRecord const out = buffers[op.out];
        int const keyRows = resolveRows(keys.rows, shared);
        int const heads = op.i[0];
        int const headDim = op.i[1];

        float const* bias = nullptr;
        int biasStride = 0;
        if ((op.flags & kFlagPairBias) && op.i[2] >= 0) {
          BufferRecord const table = buffers[op.i[2]];
          bias = arena + table.offset;  // column-major [head][pair slot]
          biasStride = table.stride;
        }

        bool const maskKeys = (op.flags & kFlagMaskKeys) != 0;
        int const keyOffset = (op.flags & kFlagKeyOffsetOne) ? 1 : 0;

        // Two instantiations keep the score row in registers without paying
        // for the largest sequence the interpreter allows on every model.
        // Two instantiations keep the score row in registers without paying
        // for the largest sequence the interpreter allows on every model.
        if (keyRows <= 17) {
          attentionImpl<17>(acc, arena + origin(op.in[2]), queries.stride,
                            arena + origin(op.in[0]), keys.stride, arena + origin(op.in[1]),
                            values.stride, arena + origin(op.out), out.stride, bias, biasStride,
                            rows, keyRows, heads, headDim, op.f[0], maskKeys, keyOffset, shared);
        } else {
          attentionImpl<kMaxTokens + 1>(
              acc, arena + origin(op.in[2]), queries.stride, arena + origin(op.in[0]), keys.stride,
              arena + origin(op.in[1]), values.stride, arena + origin(op.out), out.stride, bias,
              biasStride, rows, keyRows, heads, headDim, op.f[0], maskKeys, keyOffset, shared);
        }
        return;
      }

      case OpCode::HeadScale: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        int const sourceBase = origin(op.in[0]);
        int const targetBase = origin(op.out);
        int const heads = op.i[0];
        int const headDim = op.i[1];
        float const* const scale = params + op.t[0];
        bool const transpose = (op.flags & kFlagTransposeHeads) != 0;
        // Two nested walks: the outer over (row, head), the inner over the
        // channels of one head, so neither index needs a division.
        for (Lattice walk(lane, lanes, headDim); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row / heads;
          int const head = walk.row - row * heads;
          int const d = walk.column;
          int const channel = head * headDim + d;
          int const target = transpose ? d * heads + head : channel;
          arena[targetBase + row * dst.stride + target] =
              arena[sourceBase + row * src.stride + channel] * scale[head];
        }
        barrier(acc);
        return;
      }

      case OpCode::PairFeatures: {
        // The four-vectors are indexed by particle while this region tiles
        // over pairs, so the input is addressed absolutely.
        BufferRecord const source = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        float const* const p4 = arena + source.offset;
        int const targetBase = origin(op.out);
        for (int local = lane; local < rows; local += lanes) {
          int const slot = tileBase + local < shared.pairs ? tileBase + local : shared.pairs - 1;
          int i = 0;
          int j = 0;
          pairFromSlot(acc, slot, shared.active, i, j);
          lorentzFeatures(acc, p4 + i * source.stride, p4 + j * source.stride,
                          arena + targetBase + local * dst.stride, op.i[0]);
        }
        barrier(acc);
        return;
      }

      case OpCode::PairScatter: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const table = buffers[op.out];
        int const sourceBase = origin(op.in[0]);
        int const heads = op.i[0];
        for (Lattice walk(lane, lanes, heads); walk.at < rows * heads; walk.advance()) {
          int const local = walk.row;
          int const head = walk.column;
          int const slot = tileBase + local;
          if (slot >= shared.pairs)
            continue;
          int i = 0;
          int j = 0;
          pairFromSlot(acc, slot, shared.active, i, j);
          bool const valid = shared.mask[i] > .5f && shared.mask[j] > .5f;
          arena[table.offset + head * table.stride + slot] =
              valid ? arena[sourceBase + local * src.stride + head] : 0.f;
        }
        barrier(acc);
        return;
      }

      case OpCode::Copy: {
        BufferRecord const src = buffers[op.in[0]];
        BufferRecord const dst = buffers[op.out];
        int const sourceBase = origin(op.in[0]);
        int const targetBase = origin(op.out);
        bool const maskRows = (op.flags & kFlagMaskRows) != 0;
        for (Lattice walk(lane, lanes, cols); walk.at < rows * cols; walk.advance()) {
          int const row = walk.row;
          int const c = walk.column;
          int const token = tileBase + row;
          float value = arena[sourceBase + row * src.stride + c];
          if (maskRows && (token >= shared.tokens || shared.mask[token] <= .5f))
            value = 0.f;
          arena[targetBase + row * dst.stride + c] = value;
        }
        barrier(acc);
        return;
      }

      case OpCode::Softmax: {
        BufferRecord const dst = buffers[op.out];
        int const base = origin(op.out);

        // Down the columns instead of along the rows: the softmax an
        // attention-pooling layer takes over the particles.  One lane owns one
        // column and walks the rows, so nothing crosses lanes here either.
        if (op.flags & kFlagSoftmaxRows) {
          for (int c = lane; c < cols; c += lanes) {
            float* const column = arena + base + c;
            float maximum = -3.4e38f;
            for (int row = 0; row < rows; ++row) {
              float const value = column[row * dst.stride];
              if (value > maximum)
                maximum = value;
            }
            float total = 0.f;
            for (int row = 0; row < rows; ++row) {
              float const value = alpaka::math::exp(inner(acc), column[row * dst.stride] - maximum);
              column[row * dst.stride] = value;
              total += value;
            }
            float const inverse = 1.f / total;
            for (int row = 0; row < rows; ++row)
              column[row * dst.stride] *= inverse;
          }
          barrier(acc);
          return;
        }

        for (int row = lane; row < rows; row += lanes) {
          float* const line = arena + base + row * dst.stride;
          float maximum = line[0];
          for (int c = 1; c < cols; ++c)
            if (line[c] > maximum)
              maximum = line[c];
          float total = 0.f;
          for (int c = 0; c < cols; ++c) {
            line[c] = alpaka::math::exp(inner(acc), line[c] - maximum);
            total += line[c];
          }
          float const inverse = 1.f / total;
          for (int c = 0; c < cols; ++c)
            line[c] *= inverse;
        }
        barrier(acc);
        return;
      }

      case OpCode::Reduce: {
        // One lane owns one channel and walks the rows, so nothing has to be
        // reduced across lanes.  `rows` is the number of input rows consumed;
        // the output is a single row.
        BufferRecord const src = buffers[op.in[0]];
        int const sourceBase = origin(op.in[0]);
        int const targetBase = origin(op.out);
        bool const maskRows = (op.flags & kFlagMaskRows) != 0;
        auto const kind = static_cast<Reduction>(op.i[0]);

        // An optional per-row weight, which is what turns a softmax over the
        // particles plus a weighted sum into one op.
        float const* weights = nullptr;
        int weightStride = 0;
        if ((op.flags & kFlagWeightedRows) && op.in[1] >= 0) {
          BufferRecord const column = buffers[op.in[1]];
          weights = arena + origin(op.in[1]);
          weightStride = column.stride;
        }

        for (int c = lane; c < cols; c += lanes) {
          float total = kind == Reduction::Max ? -3.4e38f : 0.f;
          int counted = 0;
          for (int row = 0; row < rows; ++row) {
            int const token = tileBase + row;
            if (maskRows && (token >= shared.tokens || shared.mask[token] <= .5f))
              continue;
            float value = arena[sourceBase + row * src.stride + c];
            if (weights != nullptr)
              value *= weights[row * weightStride];
            if (kind == Reduction::Max) {
              if (value > total)
                total = value;
            } else {
              total += value;
            }
            ++counted;
          }
          if (kind == Reduction::Mean)
            total = counted > 0 ? total / static_cast<float>(counted) : 0.f;
          else if (kind == Reduction::Max && counted == 0)
            total = 0.f;
          arena[targetBase + c] = total;
        }
        barrier(acc);
        return;
      }

      case OpCode::Output: {
        int const base = origin(op.in[0]);
        for (int c = lane; c < op.i[0] && c < kMaxOutputs; c += lanes)
          shared.output[c] = arena[base + c];
        barrier(acc);
        return;
      }

      default:
        // ModelLoader rejects an unknown opcode before the model ever reaches
        // the device, so reaching here would mean the file changed underneath
        // us.  Do nothing rather than read arbitrary memory.
        return;
    }
  }

  // ---------------------------------------------------------------- the graph

  // Expects shared.mask to be filled and a barrier to have been executed.
  // Leaves the model's outputs in shared.output.
  template <typename TAcc, typename TShared>
  ALPAKA_FN_ACC inline void runGraph(TAcc const& acc,
                                     ModelView const& model,
                                     TShared& shared,
                                     GraphInputs const& inputs) {
    int const lane = laneOf(acc);
    int const lanes = laneCount(acc);

    if (lane == 0) {
      int last = 0;
      for (int t = 0; t < model.maxTokens; ++t)
        if (shared.mask[t] > .5f)
          last = t + 1;
      shared.active = last;
      shared.tokens = last;
      shared.pairs = last * (last + 1) / 2;
    }
    barrier(acc);

    // A null pointer in GraphInputs means the caller staged that input straight
    // into its arena buffer before calling us, so this clear must leave those
    // buffers alone.  It did not, and the wipe reduced the entire network to a
    // function of the token count -- shared.mask lives outside the arena and
    // was the only input left, so every jet with the same particle count got
    // an identical score.  The standalone benchmarks never caught it because
    // they pass real pointers and are loaded by LoadInput *after* this point.
    int featuresBegin = -1, featuresEnd = -1, vectorsBegin = -1, vectorsEnd = -1;
    if (inputs.features == nullptr && model.features.valid()) {
      featuresBegin = model.features.offset;
      featuresEnd = featuresBegin + model.features.rows * model.features.stride;
    }
    if (inputs.vectors == nullptr && model.vectors.valid()) {
      vectorsBegin = model.vectors.offset;
      vectorsEnd = vectorsBegin + model.vectors.rows * model.vectors.stride;
    }

    // Stale values are never read by an unmasked token, but they must stay
    // finite: a trailing partial GEMM tile multiplies rows that do not exist,
    // and a NaN there would not stay in its own accumulator.
    for (int index = lane; index < TShared::kFloats; index += lanes) {
      if ((index >= featuresBegin && index < featuresEnd) ||
          (index >= vectorsBegin && index < vectorsEnd))
        continue;
      shared.arena[index] = 0.f;
    }
    for (int index = lane; index < TShared::kInt8; index += lanes)
      shared.quant[index] = 0;
    for (int index = lane; index < model.outputSize && index < kMaxOutputs; index += lanes)
      shared.output[index] = 0.f;
    barrier(acc);

    if (shared.tokens == 0)
      return;

    int index = 0;
    while (index < model.opCount) {
      OpRecord const& op = model.ops[index];
      if (static_cast<OpCode>(op.opcode) == OpCode::TileBegin) {
        int end = index + 1;
        while (end < model.opCount &&
               static_cast<OpCode>(model.ops[end].opcode) != OpCode::TileEnd)
          ++end;
        int const total = resolveRows(op.rows, shared);
        int const step = op.i[0] > 0 ? op.i[0] : total;
        for (int base = 0; base < total; base += step) {
          int const remaining = total - base;
          int const rows = remaining < step ? remaining : step;
          for (int inner = index + 1; inner < end; ++inner)
            executeOp(acc, model, model.ops[inner], shared, inputs, base, rows);
        }
        index = end + 1;
        continue;
      }
      executeOp(acc, model, op, shared, inputs, 0, -1);
      ++index;
    }
  }

}  // namespace nn

#endif
