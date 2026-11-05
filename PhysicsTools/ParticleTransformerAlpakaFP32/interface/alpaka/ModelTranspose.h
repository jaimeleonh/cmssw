#ifndef PhysicsTools_ParticleTransformerAlpakaFP32_ModelTranspose_h
#define PhysicsTools_ParticleTransformerAlpakaFP32_ModelTranspose_h

// Host-only helper.  The checkpoint stores every dense weight row-major as
// [out][in], which makes a GPU thread that owns one output column walk memory
// with a stride of `in` floats.  The kernel instead wants [in][out], so that
// the threads of a warp read consecutive columns of the same input row and
// every weight fetch is a fully coalesced 128-byte transaction.
//
// The transposition is a pure relabelling: no value is modified, rounded or
// reordered within a column, so the FP32 bytes of the model are preserved.

#include <cstddef>

#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelDefinition.h"

namespace partfp32 {

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
    function(generated::classifier(view));
  }

  // Write the column-major (i.e. [in][out]) image of `source` into
  // `destination`.  Both buffers hold generated::weights_size floats and every
  // tensor keeps the offset it has in the generated header.
  inline std::size_t transposeWeights(float const* source, float* destination) {
    generated::ModelView const view{source, nullptr};
    std::size_t covered = 0;
    forEachLinear(view, [&](generated::Linear const& layer) {
      std::size_t const offset = static_cast<std::size_t>(layer.weight - source);
      float* out = destination + offset;
      for (int o = 0; o < layer.out; ++o)
        for (int i = 0; i < layer.in; ++i)
          out[static_cast<std::size_t>(i) * layer.out + o] = layer.weight[static_cast<std::size_t>(o) * layer.in + i];
      covered += static_cast<std::size_t>(layer.in) * layer.out;
    });
    return covered;
  }

}  // namespace partfp32

#endif
