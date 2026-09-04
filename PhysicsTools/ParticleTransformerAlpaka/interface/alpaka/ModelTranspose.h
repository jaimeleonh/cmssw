#ifndef PhysicsTools_ParticleTransformerAlpaka_ModelTranspose_h
#define PhysicsTools_ParticleTransformerAlpaka_ModelTranspose_h

// Host-only helper.  The exporter stores every quantized weight row-major as
// [out][in], which makes a GPU thread that owns one output column walk memory
// with a stride of `in` bytes.  The kernel instead wants [in][out], so that
// the threads of a warp read consecutive columns of the same input row and a
// warp's weight fetch is a single coalesced transaction.
//
// The transposition is a pure relabelling: no value is modified, rounded or
// requantized, and the per-output-channel scales are untouched, so the model
// is bit-for-bit the one the exporter produced.

#include <cstddef>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelDefinition.h"

namespace part {

  // Visit every dense tensor of the graph exactly once.
  template <typename TFunction>
  inline void forEachLinear(generated::ModelView view, TFunction&& function) {
    for (int i = 0; i < 3; ++i)
      function(generated::embedLinear(view, i));
    for (int i = 0; i < 4; ++i)
      function(generated::pairLinear(view, i));
    for (int i = 0; i < kParticleBlocks; ++i) {
      auto const block = generated::block(view, i);
      function(block.qkv);
      function(block.out);
      function(block.fc1);
      function(block.fc2);
    }
    for (int i = 0; i < kClassBlocks; ++i) {
      auto const block = generated::clsBlock(view, i);
      function(block.qkv);
      function(block.out);
      function(block.fc1);
      function(block.fc2);
    }
    // The classifier keeps FP32 weights in the parameter blob and so has no
    // entry in the INT8 weight image.
  }

  // Two layouts, selected per layer:
  //
  //   integer layers   [k/4][out][4] -- four consecutive inputs packed into one
  //                    32-bit word per output, so a warp's weight fetch is one
  //                    coalesced transaction and each word feeds one dp4a;
  //   float layers     [in][out]     -- the transposed layout the FP32 path
  //                    already used.
  //
  // Neither changes a value; both are pure relabellings of the exported bytes.
  // Write the device image of `source` into `destination`.  Both buffers hold
  // generated::weights_size int8 values and every tensor keeps the offset it
  // has in the generated header.
  inline std::size_t transposeWeights(std::int8_t const* source, float const* params, std::int8_t* destination) {
    generated::ModelView const view{source, params};
    std::size_t covered = 0;
    forEachLinear(view, [&](generated::Linear const& layer) {
      std::size_t const offset = static_cast<std::size_t>(layer.weight - source);
      std::int8_t* out = destination + offset;
      if (usesIntegerArithmetic(layer)) {
        for (int o = 0; o < layer.out; ++o)
          for (int i = 0; i < layer.in; ++i)
            out[(i / 4) * layer.out * 4 + o * 4 + (i % 4)] =
                layer.weight[static_cast<std::size_t>(o) * layer.in + i];
        covered += static_cast<std::size_t>(layer.in) * layer.out;
        return;
      }
      for (int o = 0; o < layer.out; ++o)
        for (int i = 0; i < layer.in; ++i)
          out[static_cast<std::size_t>(i) * layer.out + o] = layer.weight[static_cast<std::size_t>(o) * layer.in + i];
      covered += static_cast<std::size_t>(layer.in) * layer.out;
    });
    return covered;
  }

}  // namespace part

#endif
