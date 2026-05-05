#ifndef DataFormats_L1ScoutingSoA_interface_SoftJetTensorSoA_h
#define DataFormats_L1ScoutingSoA_interface_SoftJetTensorSoA_h

#include <Eigen/Core>
#include <Eigen/Dense>

#include "DataFormats/SoATemplate/interface/SoACommon.h"
#include "DataFormats/SoATemplate/interface/SoALayout.h"

namespace l1sc {

  using PFPoints = Eigen::Matrix<float, 16, 2>;
  using PFFeatures = Eigen::Matrix<float, 16, 15>;
  using PFVectors = Eigen::Matrix<float, 16, 4>;
  using PaddingMask = Eigen::Vector<float, 16>;
  GENERATE_SOA_LAYOUT(SoftJetInputTensorLayout,
                      SOA_EIGEN_COLUMN(PFPoints, points),
                      SOA_EIGEN_COLUMN(PFFeatures, features),
                      SOA_EIGEN_COLUMN(PFVectors, vectors),
                      SOA_EIGEN_COLUMN(PaddingMask, mask))
  using SoftJetInputTensorSoA = SoftJetInputTensorLayout<>;

  GENERATE_SOA_LAYOUT(SoftJetOutputTensorLayout,
                      SOA_COLUMN(float, output)
                    )
  using SoftJetOutputTensorSoA = SoftJetOutputTensorLayout<>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_SoftJetTensorSoA_h