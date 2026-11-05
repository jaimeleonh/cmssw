#ifndef PhysicsTools_ParticleTransformerAlpakaFP32_ModelTypes_h
#define PhysicsTools_ParticleTransformerAlpakaFP32_ModelTypes_h
namespace partfp32::generated {
struct Linear { int in, out; float const* weight; float const* bias; };
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
