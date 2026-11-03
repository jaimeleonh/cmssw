#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_SoftJetONNXDeviceTensor_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_SoftJetONNXDeviceTensor_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetONNXHostTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetONNXTensorSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace l1sc {

    // make the names from the top-level l1sc namespace visible for unqualified lookup
    // inside the ALPAKA_ACCELERATOR_NAMESPACE::l1sc namespace
    using namespace ::l1sc;

    using SoftJetONNXInputDeviceTensor = PortableCollection<SoftJetONNXInputTensorSoA>;

  }  // namespace l1sc

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// heterogeneous ml data checks
ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::SoftJetONNXInputDeviceTensor, l1sc::SoftJetONNXInputHostTensor);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_SoftJetONNXDeviceTensor_h
