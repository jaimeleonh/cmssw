#ifndef PhysicsTools_ParticleTransformerAlpaka_ModelDefinition_h
#define PhysicsTools_ParticleTransformerAlpaka_ModelDefinition_h

// Every translation unit that needs the exported tensor descriptors includes
// this header instead of GeneratedModel.h, so the generated file is opened in
// exactly one namespace and included exactly once.

#include <alpaka/alpaka.hpp>
#include <cstdint>

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelConstants.h"
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ModelTypes.h"

namespace part {
#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/GeneratedModel.h"

// Shared by the host packer, the scalar reference and the block kernel, so all
// three agree on which layers are evaluated with integer arithmetic and on how
// an activation is quantized.

  // Integer arithmetic needs an activation scale from quantization-aware
  // training and an input width that is a multiple of 16, which is what lets
  // the kernel fetch sixteen INT8 activations in one 16-byte load.
  ALPAKA_FN_HOST_ACC inline bool usesIntegerArithmetic(generated::Linear const& layer) {
    return layer.activation != nullptr && layer.activation[0] > 0.f && layer.in % 16 == 0;
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

}  // namespace part

#endif
