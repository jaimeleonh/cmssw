#ifndef DataFormats_L1ScoutingSoA_interface_SoftJetONNXTensorSoA_h
#define DataFormats_L1ScoutingSoA_interface_SoftJetONNXTensorSoA_h

#include <Eigen/Core>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

constexpr int kParTMaxParticles = 16;
constexpr int kParTFeatures = 15;
constexpr int kParTVectorSize = 4;

namespace l1sc {

  // Input tensors for the ONNX Runtime soft jet taggers (e.g. SoftJetIdONNX).
  //
  // Same content as SoftJetInputTensorSoA, but with the constituent and feature
  // indices transposed: each per-jet eigen block is [features, particles], so
  // that once gathered by ONNXRuntimeAlpaka it matches the layout expected by
  // the exported models (e.g. pf_features [batch, 15, n_particles]) with no
  // transpose pass. The "points" (pairwise interaction coordinates) column is
  // not present: the ONNX ParticleTransformer derives them from the vectors.
  using ParTPFFeatures = Eigen::Matrix<float, kParTFeatures, kParTMaxParticles>;
  using ParTPFVectors = Eigen::Matrix<float, kParTVectorSize, kParTMaxParticles>;
  // row vector, presented to the model as [batch, 1, n_particles]
  using ParTPFMask = Eigen::Matrix<float, 1, kParTMaxParticles>;
  GENERATE_SOA_LAYOUT(SoftJetONNXInputTensorLayout,
                      SOA_EIGEN_COLUMN(ParTPFFeatures, features),
                      SOA_EIGEN_COLUMN(ParTPFVectors, vectors),
                      SOA_EIGEN_COLUMN(ParTPFMask, mask))
  using SoftJetONNXInputTensorSoA = SoftJetONNXInputTensorLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_SoftJetONNXTensorSoA_h
