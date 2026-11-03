#ifndef DataFormats_L1ScoutingSoA_interface_SoftJetTensorSoA_h
#define DataFormats_L1ScoutingSoA_interface_SoftJetTensorSoA_h

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

#define NUM_CANDIDATES_PER_JET 16

namespace l1sc {

  using SoftJetPFPoints = Eigen::Matrix<float, 16, 2>;
  using SoftJetPFFeatures = Eigen::Matrix<float, 16, 15>;
  using SoftJetPFVectors = Eigen::Matrix<float, 16, 4>;
  using SoftJetPaddingMask = Eigen::Vector<float, 16>;
  GENERATE_SOA_LAYOUT(SoftJetInputTensorLayout,
                      SOA_EIGEN_COLUMN(SoftJetPFPoints, points),
                      SOA_EIGEN_COLUMN(SoftJetPFFeatures, features),
                      SOA_EIGEN_COLUMN(SoftJetPFVectors, vectors),
                      SOA_EIGEN_COLUMN(SoftJetPaddingMask, mask))
  using SoftJetInputTensorSoA = SoftJetInputTensorLayout<>;

  using SoftJetNewPFPoints = Eigen::Matrix<float, NUM_CANDIDATES_PER_JET, 3>;
  using SoftJetNewPFFeatures = Eigen::Matrix<float, NUM_CANDIDATES_PER_JET, 14>;
  using SoftJetNewPaddingMask = Eigen::Vector<float, NUM_CANDIDATES_PER_JET>;
  using SoftJetNewJetFeatures = Eigen::Vector<float, 4>;
  GENERATE_SOA_LAYOUT(SoftJetNewInputTensorLayout,
                      SOA_EIGEN_COLUMN(SoftJetNewPFPoints, points),
                      SOA_EIGEN_COLUMN(SoftJetNewPFFeatures, features),
                      SOA_EIGEN_COLUMN(SoftJetNewPaddingMask, mask),
                      SOA_EIGEN_COLUMN(SoftJetNewJetFeatures, jet_features))
  using SoftJetNewInputTensorSoA = SoftJetNewInputTensorLayout<>;

  GENERATE_SOA_LAYOUT(SoftJetOutputTensorLayout,
                      SOA_COLUMN(float, output)
                    )
  using SoftJetOutputTensorSoA = SoftJetOutputTensorLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_SoftJetTensorSoA_h
