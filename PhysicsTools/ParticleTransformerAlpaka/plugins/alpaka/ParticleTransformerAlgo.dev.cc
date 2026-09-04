#include "PhysicsTools/ParticleTransformerAlpaka/plugins/alpaka/ParticleTransformerAlgo.h"

#if !defined(PARTICLE_TRANSFORMER_ALGO_API_REVISION) || PARTICLE_TRANSFORMER_ALGO_API_REVISION != 5
#error "ParticleTransformerAlgo.h is stale: replace the complete ParticleTransformerAlpaka package before rebuilding"
#endif

#include "PhysicsTools/ParticleTransformerAlpaka/interface/alpaka/ParticleTransformerFastKernel.h"

#include <algorithm>
#include <type_traits>

#include "HeterogeneousCore/AlpakaInterface/interface/workdivision.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {
  namespace {

    // One block evaluates one jet.  The grid-stride loop only matters when the
    // caller caps the grid; with the default configuration each block runs a
    // single iteration.
    class JetInferenceKernel {
    public:
      template <typename TAcc, typename = std::enable_if_t<alpaka::isAccelerator<TAcc>>>
      ALPAKA_FN_ACC void operator()(TAcc const& acc,
                                    ::part::generated::ModelView model,
                                    l1sc::SoftJetInputDeviceTensor::ConstView inputs,
                                    l1sc::SoftJetOutputDeviceTensor::View outputs,
                                    uint32_t numberOfJets) const {
        auto& shared = alpaka::declareSharedVar<::part::JetShared, __COUNTER__>(acc);
        auto const lane = static_cast<uint32_t>(alpaka::getIdx<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const lanes = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Block, alpaka::Threads>(acc)[0u]);
        auto const block = static_cast<uint32_t>(alpaka::getIdx<alpaka::Grid, alpaka::Blocks>(acc)[0u]);
        auto const blocks = static_cast<uint32_t>(alpaka::getWorkDiv<alpaka::Grid, alpaka::Blocks>(acc)[0u]);

        constexpr uint32_t particles = ::part::kMaxParticles;
        constexpr uint32_t features = ::part::kInputFeatures;
        constexpr uint32_t components = ::part::kFourVector;

        for (uint32_t jet = block; jet < numberOfJets; jet += blocks) {
          auto const input = inputs[jet];

          // The Eigen-backed SoA must be addressed by matrix coordinates; it
          // cannot be reinterpreted as the kernel's row-major arrays.
          for (uint32_t index = lane; index < particles * features; index += lanes) {
            auto const particle = index / features;
            shared.na[index] = input.features()(particle, index - particle * features);
          }
          for (uint32_t index = lane; index < particles * components; index += lanes) {
            auto const particle = index / components;
            shared.p4[index] = input.vectors()(particle, index - particle * components);
          }
          for (uint32_t particle = lane; particle < particles; particle += lanes)
            shared.mask[particle] = input.mask()(particle);
          alpaka::syncBlockThreads(acc);

          ::part::runParticleTransformer(acc, model, shared);

          if (lane == 0)
            outputs.output()[jet] = shared.logits[0];  // YAML class jet_isB
          alpaka::syncBlockThreads(acc);
        }
      }
    };

  }  // namespace

  void ParticleTransformerAlgo::run(Queue& queue,
                                        std::int8_t const* weights,
                                        float const* params,
                                        l1sc::SoftJetInputDeviceTensor const& inputs,
                                        l1sc::SoftJetOutputDeviceTensor& outputs,
                                        uint32_t numberOfJets,
                                        uint32_t maxBlocks) const {
    if (numberOfJets == 0)
      return;
    auto const blocks = maxBlocks == 0 ? numberOfJets : std::min(maxBlocks, numberOfJets);
    auto const workDiv = cms::alpakatools::make_workdiv<Acc1D>(blocks, ::part::kThreadsPerJet);
    alpaka::exec<Acc1D>(queue,
                        workDiv,
                        JetInferenceKernel{},
                        ::part::generated::ModelView{weights, params},
                        inputs.const_view(),
                        outputs.view(),
                        numberOfJets);
  }

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE
