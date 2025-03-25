///
/// Description: Firmware headers
///
/// Implementation:
///    Concrete firmware implementations
///
/// \author: Jaime Leon Holgado - Imperial College

//
//

#ifndef Stage2Layer2DemuxNetMETAlgoFirmware_H
#define Stage2Layer2DemuxNetMETAlgoFirmware_H

#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxNetMETAlgo.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"

namespace l1t {

  class Stage2Layer2DemuxNetMETAlgoFirmwareImp1 : public Stage2Layer2DemuxNetMETAlgo {
  public:
    Stage2Layer2DemuxNetMETAlgoFirmwareImp1(CaloParamsHelper const* params);
    ~Stage2Layer2DemuxNetMETAlgoFirmwareImp1() override = default;
    void processEvent(const std::vector<l1t::Jet>& inputJets,
        const std::vector<l1t::EtSum>& inputSums, std::vector<l1t::EtSum>& outputSums) override;

  private:
    conifer::BDT <input_t, score_t, true> *bdt_;
  };

}  // namespace l1t

#endif
