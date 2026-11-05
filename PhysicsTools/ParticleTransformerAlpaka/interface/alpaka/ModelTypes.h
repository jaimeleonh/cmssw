#ifndef PhysicsTools_ParticleTransformerAlpaka_ModelTypes_h
#define PhysicsTools_ParticleTransformerAlpaka_ModelTypes_h
#include <cstdint>
namespace part::generated {
// `scale` is the per-output-channel weight scale.  `activation` points at two
// floats: the symmetric activation scale of this layer's input, and the scale
// of the class-block query input for the fused QKV projection.  A zero means
// the layer was left in FP32 during quantization-aware training, so the kernel
// dequantizes its weights instead of quantizing the activation.
struct Linear { int in, out; int8_t const* weight; float const* scale; float const* bias; float const* activation; };
// A layer that quantization-aware training left in FP32; its weights live in
// the parameter blob, row-major [out][in].
struct FloatLinear { int in, out; float const* weight; float const* bias; };
struct LayerNorm { int size; float eps; float const* gamma; float const* beta; };
struct NamedLinear { char const* name; Linear const* value; };
struct NamedLayerNorm { char const* name; LayerNorm const* value; };
using Norm=LayerNorm;
struct BatchNorm { int size; float eps; float const* gamma; float const* beta; float const* mean; float const* var; };
struct BlockData {
  Norm preAttn; Linear qkv; Linear out; Norm postAttn;
  Norm preFc; Linear fc1; Norm postFc; Linear fc2;
  float const* headScale; float const* residualScale;
};
}
#endif
