#ifndef DataFormats_L1ScoutingSoA_interface_alpaka_SoftJetDeviceTensor_h
#define DataFormats_L1ScoutingSoA_interface_alpaka_SoftJetDeviceTensor_h

#include "DataFormats/Portable/interface/alpaka/PortableCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetHostTensor.h"
#include "DataFormats/L1ScoutingSoA/interface/SoftJetTensorSoA.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE {

  namespace l1sc {

    // make the names from the top-level portabletest namespace visible for unqualified lookup
    // inside the ALPAKA_ACCELERATOR_NAMESPACE::portabletest namespace
    using namespace ::l1sc;

    using SoftJetInputDeviceTensor = PortableCollection<SoftJetInputTensorSoA>;
    using SoftJetOutputDeviceTensor = PortableCollection<SoftJetOutputTensorSoA>;

  }  // namespace l1sc

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE

// heterogeneous ml data checks
ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::SoftJetInputDeviceTensor,
                                      l1sc::SoftJetInputHostTensor);
ASSERT_DEVICE_MATCHES_HOST_COLLECTION(l1sc::SoftJetOutputDeviceTensor,
                                      l1sc::SoftJetOutputHostTensor);

#endif  // DataFormats_L1ScoutingSoA_interface_alpaka_SoftJetDeviceTensor_h