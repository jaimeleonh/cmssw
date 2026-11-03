#ifndef DataFormats_L1ScoutingSoA_interface_SoftJetONNXHostTensor_h
#define DataFormats_L1ScoutingSoA_interface_SoftJetONNXHostTensor_h

#include "DataFormats/Portable/interface/PortableHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetONNXTensorSoA.h"

namespace l1sc {

  using SoftJetONNXInputHostTensor = PortableHostCollection<SoftJetONNXInputTensorSoA>;

}  // namespace l1sc

#endif  // DataFormats_L1ScoutingSoA_interface_SoftJetONNXHostTensor_h
