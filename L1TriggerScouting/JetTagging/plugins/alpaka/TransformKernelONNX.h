#ifndef L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernelONNX_h
#define L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernelONNX_h

#include <cstdint>

#include <alpaka/alpaka.hpp>

#include <Eigen/Core>

#include "DataFormats/SoATemplate/interface/SoALayout.h"
#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/SoftJetONNXDeviceTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

// Raw [batch, 2] softmax scores of the ONNX tagger, before the signal-class
// score is extracted into the SoftJetOutputTensorSoA presentation. Local to
// the JetTagging plugins, not meant to be put into the event.
namespace l1sc {

  GENERATE_SOA_LAYOUT(SoftJetONNXScoresLayout,
                      SOA_COLUMN(float, b),     // score of the signal (b jet) class
                      SOA_COLUMN(float, notB))  // score of the background class
  using SoftJetONNXScoresSoA = SoftJetONNXScoresLayout<>;

}  // namespace l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  using namespace ::l1sc;
  using SoftJetONNXScoresDevice = PortableCollection<SoftJetONNXScoresSoA>;

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels {

  using namespace ::l1sc;

  // Build the ONNX input tensors: same physics content as the PyTorch transform
  // (TransformKernel), but with the per-jet eigen blocks transposed to the
  // [features, particles] layout expected by the ONNX models.
  SoftJetONNXInputDeviceTensor transformONNX(Queue& queue,
                                             const PFCandidateDeviceCollection& pf,
                                             const AssociationMapDevice& association_map,
                                             const BxLookupDevice& jetBxLookup,
                                             const VertexDeviceCollection& vertices,
                                             const BxLookupDevice& vertexBxLookup);

  // Extract the b-jet score into the final output tensor.
  void extractScores(Queue& queue, const SoftJetONNXScoresDevice& scores, SoftJetOutputDeviceTensor& output);

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc::kernels

#endif  // L1TriggerScouting_JetTagging_plugins_alpaka_TransformKernelONNX_h
