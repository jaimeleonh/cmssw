#ifndef PhysicsTools_ParticleTransformerAlpaka_graph_ModelLoader_h
#define PhysicsTools_ParticleTransformerAlpaka_graph_ModelLoader_h

// Host-side validation of a model file.
//
// The interpreter does no bounds checking: it trusts that a buffer id names a
// buffer, that a tensor offset is inside its blob, and that the arena the
// exporter packed into is the one this build provides.  All of that is checked
// exactly once, here, before anything is uploaded -- so a corrupt or
// mismatched file is a clear exception at configure time rather than a wrong
// number or a silent out-of-bounds read per jet.
//
// The file is used in place: the sections are already aligned and in device
// order, so the loader hands back pointers into the same bytes the producer
// uploads. There is no repacking step -- the exporter wrote each dense tensor
// in the layout its arithmetic path reads.

#include <cstddef>
#include <cstdint>
#include <string>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/GraphKernel.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/GraphShared.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/graph/ModelFormat.h"

namespace nn {

  struct LoadedModel {
    FileHeader const* header = nullptr;
    OpRecord const* ops = nullptr;
    BufferRecord const* buffers = nullptr;
    std::int8_t const* weights = nullptr;
    float const* params = nullptr;
    InputSlot features;
    InputSlot vectors;

    // How much of each arena the model actually touches -- the highest
    // offset + span over its buffers -- as opposed to the capacity it was
    // packed against.  The exporter packs from zero, so a small network
    // exported with the default --arena-floats still only reaches its own
    // high-water mark, and that is what decides whether it fits warp mode.
    int floatsUsed = 0;
    int int8Used = 0;

    // Whether this model can run one jet per warp (see GraphShared.h).
    bool fitsWarpMode() const {
      return floatsUsed <= kWarpArenaFloats && int8Used <= kWarpArenaInt8;
    }

    ModelView view(std::int8_t const* deviceWeights, float const* deviceParams,
                   OpRecord const* deviceOps, BufferRecord const* deviceBuffers) const {
      return ModelView{deviceOps,
                       deviceBuffers,
                       deviceWeights,
                       deviceParams,
                       static_cast<int>(header->opCount),
                       static_cast<int>(header->maxTokens),
                       static_cast<int>(header->outputSize),
                       features,
                       vectors};
    }
  };

  namespace detail {

    inline bool fits(std::uint64_t offset, std::uint64_t count, std::uint64_t elementSize,
                     std::size_t total) {
      return offset <= total && count <= (total - offset) / (elementSize ? elementSize : 1);
    }

  }  // namespace detail

  // Returns true and fills `model` when the file is usable; otherwise returns
  // false and puts a one-line reason in `error`.
  inline bool loadModel(void const* image, std::size_t bytes, LoadedModel& model,
                        std::string& error) {
    auto fail = [&error](std::string message) {
      error = std::move(message);
      return false;
    };

    if (bytes < sizeof(FileHeader))
      return fail("the model file is shorter than its header");
    auto const* base = static_cast<char const*>(image);
    auto const* header = reinterpret_cast<FileHeader const*>(base);

    if (header->magic != kFormatMagic)
      return fail("the model file does not start with the expected magic number");
    if (header->version != kFormatVersion)
      return fail("the model file is format version " + std::to_string(header->version) +
                  ", but this build reads version " + std::to_string(kFormatVersion));

    if (!detail::fits(header->opOffset, header->opCount, sizeof(OpRecord), bytes))
      return fail("the op table runs past the end of the file");
    if (!detail::fits(header->bufferOffset, header->bufferCount, sizeof(BufferRecord), bytes))
      return fail("the buffer table runs past the end of the file");
    if (!detail::fits(header->weightsOffset, header->weightsCount, 1, bytes))
      return fail("the weight blob runs past the end of the file");
    if (!detail::fits(header->paramsOffset, header->paramsCount, sizeof(float), bytes))
      return fail("the parameter blob runs past the end of the file");

    auto const* ops = reinterpret_cast<OpRecord const*>(base + header->opOffset);
    auto const* buffers = reinterpret_cast<BufferRecord const*>(base + header->bufferOffset);
    auto const* weights = reinterpret_cast<std::int8_t const*>(base + header->weightsOffset);
    auto const* params = reinterpret_cast<float const*>(base + header->paramsOffset);

    if (fnv1a(weights, header->weightsCount) != header->weightsHash)
      return fail("the weight blob does not match the hash in the header");
    if (fnv1a(params, header->paramsCount * sizeof(float)) != header->paramsHash)
      return fail("the parameter blob does not match the hash in the header");

    // What this build can actually provide.
    if (header->arenaFloats > static_cast<std::uint32_t>(kArenaFloats))
      return fail("the model was packed for a " + std::to_string(header->arenaFloats) +
                  "-float arena but this build provides " + std::to_string(kArenaFloats) +
                  "; re-export with --arena-floats " + std::to_string(kArenaFloats));
    if (header->arenaInt8 > static_cast<std::uint32_t>(kArenaInt8))
      return fail("the model was packed for a " + std::to_string(header->arenaInt8) +
                  "-byte int8 arena but this build provides " + std::to_string(kArenaInt8));
    if (header->maxTokens > static_cast<std::uint32_t>(kMaxTokens))
      return fail("the model has a sequence length of " + std::to_string(header->maxTokens) +
                  ", above the interpreter's limit of " + std::to_string(kMaxTokens));
    if (header->outputSize > static_cast<std::uint32_t>(kMaxOutputs))
      return fail("the model publishes " + std::to_string(header->outputSize) +
                  " outputs, above the interpreter's limit of " + std::to_string(kMaxOutputs));

    // Buffers: every one must lie inside its arena.
    auto const maxRows = [&](BufferRecord const& buffer) -> std::int64_t {
      if (buffer.rows >= 0)
        return buffer.rows;
      switch (buffer.rows) {
        case kRowsTokensPlusOne:
          return header->maxTokens + 1;
        case kRowsPairs:
          return static_cast<std::int64_t>(header->maxTokens) * (header->maxTokens + 1) / 2;
        case kRowsOne:
          return 1;
        default:
          return header->maxTokens;
      }
    };

    int floatsUsed = 0;
    int int8Used = 0;
    for (std::uint32_t index = 0; index < header->bufferCount; ++index) {
      BufferRecord const& buffer = buffers[index];
      if (buffer.offset < 0 || buffer.stride <= 0 || buffer.cols <= 0)
        return fail("buffer " + std::to_string(index) + " has a nonsensical shape");
      std::int64_t const span = buffer.layout == kColumnMajor
                                    ? static_cast<std::int64_t>(buffer.cols) * buffer.stride
                                    : maxRows(buffer) * buffer.stride;
      std::int64_t const capacity = buffer.arena == kInt8Arena ? kArenaInt8 : kArenaFloats;
      if (buffer.offset + span > capacity)
        return fail("buffer " + std::to_string(index) + " runs past the end of its arena");
      std::int64_t const reach = buffer.offset + span;
      if (buffer.arena == kInt8Arena)
        int8Used = reach > int8Used ? static_cast<int>(reach) : int8Used;
      else
        floatsUsed = reach > floatsUsed ? static_cast<int>(reach) : floatsUsed;
    }

    // Ops: opcodes this build knows, buffer ids in range, tensors inside their
    // blobs, and tiles no deeper than the normalization scratch.
    for (std::uint32_t index = 0; index < header->opCount; ++index) {
      OpRecord const& op = ops[index];
      if (op.opcode >= static_cast<std::uint16_t>(OpCode::kOpCodeCount))
        return fail("op " + std::to_string(index) + " uses opcode " +
                    std::to_string(op.opcode) + ", which this build does not implement");

      for (int slot = 0; slot < kMaxOpInputs; ++slot)
        if (op.in[slot] >= 0 && static_cast<std::uint32_t>(op.in[slot]) >= header->bufferCount)
          return fail("op " + std::to_string(index) + " names an input buffer that does not exist");
      if (op.out >= 0 && static_cast<std::uint32_t>(op.out) >= header->bufferCount)
        return fail("op " + std::to_string(index) + " names an output buffer that does not exist");

      auto const code = static_cast<OpCode>(op.opcode);
      if (code == OpCode::TileBegin) {
        if (op.i[0] <= 0 || op.i[0] > kMaxReduceRows)
          return fail("op " + std::to_string(index) + " opens a tile of " +
                      std::to_string(op.i[0]) + " rows; the limit is " +
                      std::to_string(kMaxReduceRows));
      }

      if (code == OpCode::Attention && (op.flags & kFlagPairBias) && op.i[2] >= 0 &&
          static_cast<std::uint32_t>(op.i[2]) >= header->bufferCount)
        return fail("op " + std::to_string(index) + " names a pair-bias buffer that does not exist");

      if (code == OpCode::Linear) {
        if (op.in[0] < 0 || op.out < 0)
          return fail("op " + std::to_string(index) + " is a Linear without both buffers");
        std::int64_t const kdim = buffers[op.in[0]].cols;
        std::int64_t const total = op.i[2];
        if (total <= 0 || op.i[0] < 0 || op.i[0] + op.cols > total)
          return fail("op " + std::to_string(index) + " projects columns outside its weight");
        std::int64_t const elements = kdim * total;
        bool const int8 = (op.flags & kFlagInt8Weights) != 0;
        std::uint64_t const limit = int8 ? header->weightsCount : header->paramsCount;
        if (op.t[0] < 0 || static_cast<std::uint64_t>(op.t[0]) + elements > limit)
          return fail("op " + std::to_string(index) + " has a weight tensor outside its blob");
        if ((op.flags & kFlagInt8Input) && op.t[3] < 0)
          return fail("op " + std::to_string(index) +
                      " reads INT8 activations but carries no activation scale");
      }

      // Every parameter vector this op names has to be inside the blob.  The
      // widths differ per op, so check the offset and one element; a tensor
      // that starts inside the blob and is sized by `cols` is covered by the
      // buffer check above.
      for (int slot = (code == OpCode::Linear ? 1 : 0); slot < kMaxOpTensors; ++slot)
        if (op.t[slot] >= 0 && static_cast<std::uint64_t>(op.t[slot]) >= header->paramsCount)
          return fail("op " + std::to_string(index) + " names a parameter outside the blob");

      if (code == OpCode::LoadInput && op.out >= 0) {
        InputSlot slot;
        slot.offset = buffers[op.out].offset;
        slot.stride = buffers[op.out].stride;
        slot.rows = static_cast<int>(maxRows(buffers[op.out]));
        slot.cols = op.cols;
        if (op.i[0] == 0)
          model.features = slot;
        else if (op.i[0] == 1)
          model.vectors = slot;
      }
    }

    model.header = header;
    model.floatsUsed = floatsUsed;
    model.int8Used = int8Used;
    model.ops = ops;
    model.buffers = buffers;
    model.weights = weights;
    model.params = params;
    return true;
  }

}  // namespace nn

#endif
