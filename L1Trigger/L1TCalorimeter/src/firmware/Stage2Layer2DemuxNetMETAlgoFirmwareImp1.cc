///
/// \class l1t::Stage2Layer2NetMETAlgorithmFirmwareImp1
///
/// \author:
///
/// Description:

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "L1Trigger/L1TCalorimeter/interface/Stage2Layer2DemuxNetMETAlgoFirmware.h"

#include "L1Trigger/L1TCalorimeter/interface/CaloParamsHelper.h"
#include "L1Trigger/L1TCalorimeter/interface/CaloTools.h"

#ifndef CMSSW_GIT_HASH
#include "hls_math.h"
#endif

#include <vector>
#include <numeric>
#include <algorithm>

#include "L1Trigger/L1TCalorimeter/interface/BitonicSort.h"

l1t::Stage2Layer2DemuxNetMETAlgoFirmwareImp1::Stage2Layer2DemuxNetMETAlgoFirmwareImp1(CaloParamsHelper const* params) {
  bdt_ = new conifer::BDT<input_t, score_t, true>(params->netMetModelFile());
}

void l1t::Stage2Layer2DemuxNetMETAlgoFirmwareImp1::processEvent(const std::vector<l1t::Jet>& inputJets,
    const std::vector<l1t::EtSum>& inputSums, std::vector<l1t::EtSum>& outputSums) {
  outputSums = inputSums;

  std::vector<input_t> inputs = 
  {
    inputJets.size() > 0 ? inputJets[0].hwEta() : 0,
    inputJets.size() > 0 ? inputJets[0].hwPhi() : 0,
    inputJets.size() > 0 ? inputJets[0].hwPt() : 0,
    inputJets.size() > 1 ? inputJets[1].hwEta() : 0,
    inputJets.size() > 1 ? inputJets[1].hwPhi() : 0,
    inputJets.size() > 1 ? inputJets[1].hwPt() : 0,
    inputJets.size() > 2 ? inputJets[2].hwEta() : 0,
    inputJets.size() > 2 ? inputJets[2].hwPhi() : 0,
    inputJets.size() > 2 ? inputJets[2].hwPt() : 0,
    inputJets.size() > 3 ? inputJets[3].hwEta() : 0,
    inputJets.size() > 3 ? inputJets[3].hwPhi() : 0,
    inputJets.size() > 3 ? inputJets[3].hwPt() : 0,
    inputSums[11].hwPt(),
    inputSums[12].hwPt(),
    inputSums[13].hwPt()
  };
  std::vector<score_t> netmet = bdt_->decision_function(inputs);

  // Replace met noHF with the NetMET
  outputSums[7].setHwPt(netmet[0].to_int());
}

