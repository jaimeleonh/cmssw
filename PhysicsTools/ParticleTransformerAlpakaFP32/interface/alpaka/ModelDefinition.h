#ifndef PhysicsTools_ParticleTransformerAlpakaFP32_ModelDefinition_h
#define PhysicsTools_ParticleTransformerAlpakaFP32_ModelDefinition_h

// Every translation unit that needs the exported tensor descriptors includes
// this header instead of GeneratedModel.h, so the generated file is opened in
// exactly one namespace and included exactly once.

#include <alpaka/alpaka.hpp>

#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelConstants.h"
#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/ModelTypes.h"

namespace partfp32 {
#include "PhysicsTools/ParticleTransformerAlpakaFP32/interface/alpaka/GeneratedModel.h"
}  // namespace partfp32

#endif
