///
/// \class l1t::Stage2Layer2DemuxJetAlgo
///
/// Description: interface for demux Net MET algorithm
///
/// Implementation:
///
/// \author: Jaime Leon Holgado - Imperial College
///

//

#ifndef Stage2Layer2DemuxNetMETAlgo_h
#define Stage2Layer2DemuxNetMETAlgo_h

#include "conifer.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include <vector>

#include "ap_fixed.h"

typedef ap_fixed<24, 16> input_t;
typedef ap_fixed<24, 16, AP_RND_CONV, AP_SAT> score_t;

namespace l1t {

  class Stage2Layer2DemuxNetMETAlgo {
  public:
    virtual ~Stage2Layer2DemuxNetMETAlgo() {}
    virtual void processEvent(const std::vector<l1t::Jet>& inputJets,
        const std::vector<l1t::EtSum>& inputSums, std::vector<l1t::EtSum>& outputSums) = 0;
  private:
  };

}  // namespace l1t

#endif
