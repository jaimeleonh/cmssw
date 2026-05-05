#ifndef DataFormats_L1ScoutingSoA_interface_SoftJetHostTensor_h
#define DataFormats_L1ScoutingSoA_interface_SoftJetHostTensor_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetTensorSoA.h"

namespace l1sc {

  using SoftJetInputHostTensor = PortableHostCollection<SoftJetInputTensorSoA>;
  using SoftJetOutputHostTensor = PortableHostCollection<SoftJetOutputTensorSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_SoftJetHostTensor_h