#ifndef PhysicsTools_ParticleTransformerAlpaka_graph_ModelFormat_h
#define PhysicsTools_ParticleTransformerAlpaka_graph_ModelFormat_h

// The on-disk contract between the exporter (scripts/nnexport) and the
// interpreter (GraphKernel.h).  It is deliberately free of Alpaka, CMSSW and
// C++ library dependencies so that the host loader, the device kernel and any
// standalone tool can all include it.
//
// A model file is one blob:
//
//   [FileHeader]         magic, version, section extents, payload hashes
//   [OpRecord    x nOps] the program, in execution order
//   [BufferRecord x nBufs] the activation arena layout, precomputed by the
//                        exporter's allocator
//   [int8  weights]      quantized dense tensors, row-major [out][in]
//   [float params]       scales, biases, normalisation parameters, FP32 weights
//   [string table]       NUL-separated names, diagnostics only
//
// Every section is 16-byte aligned so that the loader can hand the device a
// pointer into a single upload rather than four separate buffers.
//
// Nothing in the file is a pointer.  A tensor is named by its element offset
// into the weight or parameter blob, and a buffer by its index in the buffer
// table, so the same bytes are valid on the host and on every device.

#include <cstdint>

namespace nn {

  // Bumped whenever the meaning of an existing field changes.  Adding a new
  // opcode at the end of the enum, or a new flag bit, does not require a bump:
  // the loader rejects an opcode it does not know, which is the behaviour we
  // want for a model that uses a layer this release cannot evaluate.
  inline constexpr std::uint32_t kFormatMagic = 0x4e4e4731u;  // "NNG1"
  inline constexpr std::uint32_t kFormatVersion = 1u;

  inline constexpr int kMaxOpInputs = 3;
  inline constexpr int kMaxOpTensors = 6;
  inline constexpr int kMaxOpInts = 3;
  inline constexpr int kMaxOpFloats = 2;

  // ---------------------------------------------------------------- row specs
  //
  // Most buffers have a row count that is only known once the jet is read --
  // the number of unmasked particles.  A non-negative `rows` field is a
  // literal; a negative one names one of these runtime quantities.  The
  // interpreter resolves them once per jet in resolveRows().

  enum RowSpec : std::int32_t {
    kRowsTokens = -1,      // padded active token count
    kRowsTokensPlusOne = -2,  // tokens + 1, the class-attention key sequence
    kRowsPairs = -3,       // padded count of i <= j pairs among active tokens
    kRowsActive = -4,      // unpadded active token count
    kRowsOne = -5,         // a single row; equivalent to a literal 1
  };

  // ---------------------------------------------------------------- opcodes

  enum class OpCode : std::uint16_t {
    Nop = 0,

    // Region markers.  TileBegin opens a loop over the rows named by `rows`,
    // in strides of i[0].  Buffers whose `tiled` flag is set are addressed
    // relative to the current tile base; the others are addressed absolutely.
    // Regions do not nest.
    TileBegin = 1,
    TileEnd = 2,

    // out[r][c] = bias[c] + scale[c] * sum_k in[r][k] * weight[k][c]
    //
    // in[0] input buffer, out output buffer.
    //   t[0] weight     -- int8 blob when kFlagInt8Weights, else float blob
    //   t[1] scale      -- per-output-channel weight scale, -1 for an FP32 layer
    //   t[2] bias
    //   t[3] activation -- the input's quantization scale, required when
    //        kFlagInt8Input is set so the INT32 accumulator can be scaled back
    //   i[0] column offset into the weight tensor's output axis, so that one
    //        fused projection can be evaluated in slices (QKV -> KV and Q)
    //   i[1] number of output columns actually produced
    //   i[2] the full extent of that output axis, which is the tensor's row
    //        stride; a sliced projection still indexes the whole tensor
    Linear = 3,

    // Row-wise LayerNorm over `cols` channels.  t[0] gamma, t[1] beta,
    // f[0] epsilon.
    //   flags: kFlagNormUnbiased divides the sum of squares by cols - 1
    //          rather than cols, which is what torch.std() does by default;
    //          kFlagNormEpsOutside makes the scale 1 / (sqrt(var) + eps)
    //          instead of 1 / sqrt(var + eps).  Both exist because a model
    //          may spell its normalization with mean() and std() rather than
    //          nn.LayerNorm, and the difference is not small: over 128
    //          channels the two variance conventions differ by 0.4%.
    LayerNorm = 4,

    // Per-channel BatchNorm in inference form.  t[0] gamma, t[1] beta,
    // t[2] running mean, t[3] running variance, f[0] epsilon.
    BatchNorm = 5,

    // Elementwise activation; i[0] is an Activation kind.
    Activate = 6,

    // out[r][c] = alpha[c] * out[r][c] + in[0][r][c].  t[0] is a per-channel
    // alpha, or -1 for a plain sum.
    Residual = 7,

    // Quantize a float buffer into an int8 buffer using a symmetric scale.
    // t[0] points at the scale (one float), f[0] is a literal fallback used
    // when t[0] is -1.
    Quantize = 8,

    // Snap a float buffer onto the int8 grid in place, leaving it float.  The
    // projections whose input width is not a multiple of 16 keep FP32
    // arithmetic but must still see the activation the trained model saw.
    RoundToGrid = 9,

    // Fused multi-head attention.  in[0] keys, in[1] values, in[2] queries,
    // out the merged per-token result.  The key count is the `keys` buffer's
    // own row specification, so self-attention over 16 particles and class
    // attention over 17 tokens are the same op.
    //   i[0] heads, i[1] head dimension, i[2] pair-bias buffer id (-1 if none)
    //   f[0] score scale
    //   flags: kFlagPairBias adds the pair bias, kFlagMaskKeys forces masked
    //   keys to -inf, kFlagKeyOffsetOne makes key r correspond to particle
    //   r-1 (the class-attention sequence has the class token first).
    Attention = 10,

    // out[perm(h,d)] = in[0][h * headDim + d] * headScale[h].
    // t[0] head scale, i[0] heads, i[1] head dimension.  kFlagTransposeHeads
    // selects the (h,d) -> (d,h) destination permutation.
    HeadScale = 11,

    // Build the pairwise Lorentz features of the i <= j pairs into `out`,
    // one pair per row.  i[0] is the number of features produced.
    //
    // in[0] is the four-vector buffer, and it is addressed ABSOLUTELY: the
    // region around this op tiles over pairs, while the four-vectors are
    // indexed by particle, so the current tile base does not apply to it.
    PairFeatures = 12,

    // Scatter a pair-embedding tile into the attention-bias buffer, zeroing
    // pairs that involve a masked particle.  i[0] heads.
    PairScatter = 13,

    // Copy a parameter vector into a buffer row.  t[0] source, i[0] length.
    LoadConst = 14,

    // out[r][c] = in[0][r][c], with kFlagMaskRows zeroing masked rows.  The
    // row index is taken in the source buffer's own frame, so this is also how
    // a tile is written back into a persistent buffer.
    Copy = 15,

    // Row-wise softmax over `cols`.  kFlagSoftmaxRows transposes the sense of
    // it: the normalization then runs down each column over `rows` rows, which
    // is the softmax an attention-pooling layer takes over the particles.
    Softmax = 16,

    // Copy the jet's input features, four-vectors or mask into a buffer.
    // i[0] selects which: 0 features, 1 four-vectors, 2 mask.
    LoadInput = 17,

    // Publish a buffer row as the model output.  i[0] is the number of values.
    Output = 18,

    // Pool over the row axis: out[0][c] = reduce over r of in[0][r][c].  This
    // is what turns a per-particle stack into a jet-level one, and it is the
    // only op whose output has fewer rows than its input.
    //   i[0] a Reduction kind
    //   flags: kFlagMaskRows skips masked rows, so padding cannot contribute
    //   and a mean divides by the number of real particles;
    //   kFlagWeightedRows multiplies row r by in[1][r][0] first, which is how
    //   attention pooling -- softmax over the particles, then a weighted sum
    //   -- becomes a single op.  A weighted Mean divides by the number of
    //   rows counted, not by the sum of the weights.
    Reduce = 19,

    kOpCodeCount
  };

  enum class Reduction : std::int32_t {
    Sum = 0,
    Mean = 1,
    Max = 2,
    kReductionCount
  };

  enum class Activation : std::int32_t {
    Identity = 0,
    Gelu = 1,   // exact erf form, as trained
    Relu = 2,
    Silu = 3,
    Tanh = 4,
    Sigmoid = 5,
    kActivationCount
  };

  // ---------------------------------------------------------------- flags

  inline constexpr std::uint16_t kFlagInt8Weights = 1u << 0;   // Linear
  inline constexpr std::uint16_t kFlagInt8Input = 1u << 1;     // Linear reads an int8 buffer
  inline constexpr std::uint16_t kFlagPairBias = 1u << 2;      // Attention
  inline constexpr std::uint16_t kFlagMaskKeys = 1u << 3;      // Attention
  inline constexpr std::uint16_t kFlagKeyOffsetOne = 1u << 4;  // Attention
  inline constexpr std::uint16_t kFlagTransposeHeads = 1u << 5;  // HeadScale
  inline constexpr std::uint16_t kFlagMaskRows = 1u << 6;      // Copy, Reduce
  inline constexpr std::uint16_t kFlagInPlace = 1u << 7;       // output aliases input
  inline constexpr std::uint16_t kFlagNormUnbiased = 1u << 8;    // LayerNorm
  inline constexpr std::uint16_t kFlagNormEpsOutside = 1u << 9;  // LayerNorm
  inline constexpr std::uint16_t kFlagSoftmaxRows = 1u << 10;    // Softmax
  inline constexpr std::uint16_t kFlagWeightedRows = 1u << 11;   // Reduce

  // ---------------------------------------------------------------- records

  // 64 bytes, so a warp's read of the next op is one transaction and the op
  // table of a large model still fits comfortably in L1.
  struct alignas(16) OpRecord {
    std::uint16_t opcode;
    std::uint16_t flags;
    std::int16_t in[kMaxOpInputs];  // input buffer ids, -1 when unused
    std::int16_t out;               // output buffer id, -1 when unused
    std::int32_t rows;              // RowSpec or a literal row count
    std::int32_t cols;              // channels per row
    std::int32_t i[kMaxOpInts];
    std::int32_t t[kMaxOpTensors];  // blob element offsets, -1 when absent
    float f[kMaxOpFloats];
  };
  static_assert(sizeof(OpRecord) == 64, "the op record must stay 64 bytes");

  // Which arena a buffer lives in.  Named apart from the arena *sizes* in
  // GraphShared.h, which are a different kind of thing entirely.
  enum Arena : std::uint8_t {
    kFloatArena = 0,
    kInt8Arena = 1,
  };

  // How a dense tensor is laid out in the weight blob.  The exporter writes
  // each one in the layout its arithmetic path reads, so the host uploads the
  // file verbatim and there is no repacking step on the way to the device.
  //
  // With `out` the full output extent (i[2]) and `in` the input width, element
  // (o, k) of a tensor lives at
  //
  //   kWeightOutIn   o * in + k              the narrow FP32 head
  //   kWeightInOut   k * out + o             consecutive threads read
  //                                          consecutive output columns, so a
  //                                          warp's fetch is one transaction
  //   kWeightDp4a    (k/4) * out * 4         four inputs packed per 32-bit
  //                    + o * 4 + (k%4)       word, each word feeding one dp4a
  //
  // The layout is implied by the Linear flags: kFlagInt8Input means dp4a,
  // kFlagInt8Weights alone means [in][out], neither means [out][in].
  enum WeightLayout : std::uint8_t {
    kWeightOutIn = 0,
    kWeightInOut = 1,
    kWeightDp4a = 2,
  };

  enum Layout : std::uint8_t {
    kRowMajor = 0,     // element (r, c) at offset + r * stride + c
    kColumnMajor = 1,  // element (r, c) at offset + c * stride + r
  };

  struct alignas(8) BufferRecord {
    std::int32_t offset;  // floats into the float arena, or bytes into the int8 arena
    std::int32_t rows;    // RowSpec or literal; the capacity the allocator reserved
    std::int32_t cols;
    std::int32_t stride;  // distance between consecutive rows (or columns, if column-major)
    std::uint8_t arena;
    std::uint8_t tiled;   // rows are relative to the enclosing TileBegin
    std::uint8_t layout;
    std::uint8_t pad_;
  };
  static_assert(sizeof(BufferRecord) == 24, "the buffer record must stay 24 bytes");

  struct alignas(16) FileHeader {
    std::uint32_t magic;
    std::uint32_t version;

    std::uint32_t opCount;
    std::uint32_t bufferCount;

    std::uint64_t opOffset;      // byte offsets from the start of the file
    std::uint64_t bufferOffset;
    std::uint64_t weightsOffset;
    std::uint64_t paramsOffset;
    std::uint64_t stringsOffset;

    std::uint64_t weightsCount;  // int8 elements
    std::uint64_t paramsCount;   // float elements
    std::uint64_t stringsBytes;

    // What the interpreter must be able to provide for this model to run.
    std::uint32_t arenaFloats;   // float arena the allocator packed into
    std::uint32_t arenaInt8;     // int8 staging arena, in bytes
    std::uint32_t maxTokens;     // padded sequence length
    std::uint32_t inputFeatures;
    std::uint32_t fourVectorSize;
    std::uint32_t outputSize;
    std::uint32_t heads;         // largest head count of any Attention op
    std::uint32_t pairFeatures;  // 0 when the model has no pair embedding

    std::uint64_t weightsHash;   // FNV-1a over the payload bytes
    std::uint64_t paramsHash;
  };
  static_assert(sizeof(FileHeader) == 128, "the file header must stay 128 bytes");

  // The loader and the exporter both compute this over the payload sections, so
  // a truncated or mismatched file is rejected before it reaches the device.
  inline std::uint64_t fnv1a(void const* data, std::uint64_t bytes) {
    auto const* byte = static_cast<unsigned char const*>(data);
    std::uint64_t hash = 14695981039346656037ull;  // 0xCBF29CE484222325
    for (std::uint64_t index = 0; index < bytes; ++index) {
      hash ^= byte[index];
      hash *= 1099511628211ull;
    }
    return hash;
  }

}  // namespace nn

#endif
