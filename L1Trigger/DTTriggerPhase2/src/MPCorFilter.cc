#include "L1Trigger/DTTriggerPhase2/interface/MPCorFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;
using namespace std;
using namespace cmsdt;


// ============================================================================
// Constructors and destructor
// ============================================================================
MPCorFilter::MPCorFilter(const ParameterSet &pset)
    : MPFilter(pset), debug_(pset.getUntrackedParameter<bool>("debug")) {}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MPCorFilter::initialise(const edm::EventSetup &iEventSetup) {}

void MPCorFilter::run(edm::Event &iEvent,
                                  const edm::EventSetup &iEventSetup,
                                  std::vector<metaPrimitive> &inSLMPaths,
                                  std::vector<metaPrimitive> &inCorMPaths,
                                  std::vector<metaPrimitive> &outMPaths) {

  if (debug_)
    LogDebug("MPCorFilter") << "MPCorFilter: run";

  std::map<int, std::vector<metaPrimitive>> SL1metaPrimitivesPerBX;
  std::map<int, std::vector<metaPrimitive>> SL3metaPrimitivesPerBX;
  std::map<int, std::vector<metaPrimitive>> CormetaPrimitivesPerBX;
  uint32_t sl1Id_rawid = -1, sl3Id_rawid = -1;
  if (inSLMPaths.size() > 0) {
    int dum_sl_rawid = inSLMPaths[0].rawId;
    DTSuperLayerId dumSlId(dum_sl_rawid);
    DTChamberId ChId(dumSlId.wheel(), dumSlId.station(), dumSlId.sector());
    DTSuperLayerId sl1Id(ChId.rawId(), 1);
    sl1Id_rawid = sl1Id.rawId();
    DTSuperLayerId sl3Id(ChId.rawId(), 3);
    sl3Id_rawid = sl3Id.rawId();

    for (const auto &metaprimitiveIt : inSLMPaths) {
      int BX = metaprimitiveIt.t0 / 25;
      if (metaprimitiveIt.rawId == sl1Id_rawid)
        SL1metaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
      else if (metaprimitiveIt.rawId == sl3Id_rawid)
        SL3metaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
    }
  }
  for (const auto &metaprimitiveIt : inCorMPaths) {
    int BX = metaprimitiveIt.t0 / 25;
    CormetaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
  }

  auto filteredMPs = filter(SL1metaPrimitivesPerBX, SL3metaPrimitivesPerBX, CormetaPrimitivesPerBX);
  for (auto & mp: filteredMPs)
    outMPaths.push_back(mp);
}

void MPCorFilter::finish(){};

///////////////////////////
///  OTHER METHODS

std::vector<metaPrimitive> MPCorFilter::filter(
    std::map<int, std::vector<metaPrimitive>> SL1mpsPerBX,
    std::map<int, std::vector<metaPrimitive>> SL3mpsPerBX,
    std::map<int, std::vector<metaPrimitive>> CormpsPerBX)
{
  std::map<int, valid_cor_tp_arr_t> mp_valid_per_bx;
  for (auto &elem: SL1mpsPerBX) {
    mp_valid_per_bx[elem.first] = valid_cor_tp_arr_t(12);
    int imp = 0;
    for (auto &mp : elem.second) {
      auto coarsed = coarsify(mp, 1);
      mp_valid_per_bx[elem.first][imp] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
      imp += 2;
    }
  }
  for (auto &elem: SL3mpsPerBX) {
    if (mp_valid_per_bx.find(elem.first) == mp_valid_per_bx.end())
      mp_valid_per_bx[elem.first] = valid_cor_tp_arr_t(12);
    int imp = 1;
    for (auto &mp : elem.second) {
      auto coarsed = coarsify(mp, 3);
      mp_valid_per_bx[elem.first][imp] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
      imp += 2;
    }
  }
  for (auto &elem: CormpsPerBX) {
    if (mp_valid_per_bx.find(elem.first) == mp_valid_per_bx.end()) {
      mp_valid_per_bx[elem.first] = valid_cor_tp_arr_t(12);
    }
    for (auto &mp : elem.second) {
      auto coarsed = coarsify(mp, 0);
      if (isDead(mp, coarsed, mp_valid_per_bx)) continue;
      auto index = killTps(mp, coarsed, elem.first, mp_valid_per_bx);
      mp_valid_per_bx[elem.first][index] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    }
  }

  std::vector<metaPrimitive> outTPs;
  for (auto &elem: mp_valid_per_bx) {
    for (auto &mp_valid : elem.second) {
      if (mp_valid.valid) {
        outTPs.push_back(mp_valid.mp);
      }
    }
  }

  return outTPs;
}

std::vector<int> MPCorFilter::coarsify(cmsdt::metaPrimitive mp, int sl) {
  float sign = 0;
  if (sl == 1) sign = -1;
  else if (sl == 3) sign = 1;
  float pos_ch_f = mp.x + sign * mp.tanPhi * VERT_PHI1_PHI3 / 2;

  // translating into tdc counts
  int pos_ch = int(round(pos_ch_f / (((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC) / 10)));
  int slope = (int) (-mp.tanPhi / (SLOPE_LSB * INCREASED_RES_SLOPE_POW));

  std::cout << sl << " " << pos_ch_f << " " << pos_ch << std::endl; 

  std::vector<int> t0_slv, t0_coarse, pos_slv, pos_coarse, slope_slv, slope_coarse;
  vhdl_int_to_unsigned(mp.t0, t0_slv);
  vhdl_int_to_signed(pos_ch, pos_slv);
  vhdl_int_to_signed(slope, slope_slv);

  for (size_t i = 0; i < pos_slv.size(); i++)
    std::cout << pos_slv[i];
  std::cout << std::endl;

  vhdl_resize_unsigned(t0_slv, WIDTH_FULL_TIME);
  vhdl_resize_signed(pos_slv, WIDTH_FULL_POS);
  vhdl_resize_signed(slope_slv, WIDTH_FULL_SLOPE);

  for (size_t i = 0; i < pos_slv.size(); i++)
    std::cout << pos_slv[i];
  std::cout << std::endl;

  t0_coarse = vhdl_slice(t0_slv, FSEG_T0_BX_LSB + 4, FSEG_T0_DISCARD_LSB - 1);
  pos_coarse = vhdl_slice(pos_slv, WIDTH_FULL_POS - 1, FSEG_POS_DISCARD_LSB - 1);
  slope_coarse = vhdl_slice(slope_slv, WIDTH_FULL_SLOPE - 1, FSEG_SLOPE_DISCARD_LSB - 1);

  for (size_t i = 0; i < pos_coarse.size(); i++)
    std::cout << pos_coarse[i];
  std::cout << std::endl;

  std::vector <int> results;
  int t0_coarse_int = vhdl_unsigned_to_int(t0_coarse);
  int pos_coarse_int = vhdl_signed_to_int(pos_coarse);
  int slope_coarse_int = vhdl_signed_to_int(slope_coarse);

  for (int index = 0; index <= 2; index++) {
    auto aux_t0_coarse_int = t0_coarse_int + (2 * index - 1);
    auto aux_pos_coarse_int = pos_coarse_int + (2 * index - 1);
    auto aux_slope_coarse_int = slope_coarse_int + (2 * index - 1);
    results.push_back(aux_t0_coarse_int >> 1);
    results.push_back(aux_pos_coarse_int >> 1);
    results.push_back(aux_slope_coarse_int >> 1);
  }
  return results;  
}


int MPCorFilter::match(cmsdt::metaPrimitive mp, std::vector<int> coarsed, valid_cor_tp_t valid_cor_tp2) {
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
    std::cout << coarsed[3 * i + j] << " ";
  std::cout << std::endl;
  std::cout << valid_cor_tp2.coarsed_t0 << " " << valid_cor_tp2.coarsed_pos << " " << valid_cor_tp2.coarsed_slope << std::endl;

  bool matched = (
    (coarsed[0] == valid_cor_tp2.coarsed_t0    || coarsed[3] == valid_cor_tp2.coarsed_t0    || coarsed[6] == valid_cor_tp2.coarsed_t0)  &&
    (coarsed[1] == valid_cor_tp2.coarsed_pos   || coarsed[4] == valid_cor_tp2.coarsed_pos   || coarsed[7] == valid_cor_tp2.coarsed_pos) &&
    (coarsed[2] == valid_cor_tp2.coarsed_slope || coarsed[5] == valid_cor_tp2.coarsed_slope || coarsed[8] == valid_cor_tp2.coarsed_slope)
  );
  return ((int) matched) * 2 + (int) (mp.quality > valid_cor_tp2.mp.quality);
}

bool MPCorFilter::isDead(cmsdt::metaPrimitive mp, std::vector<int> coarsed, std::map<int, valid_cor_tp_arr_t> tps_per_bx) {
  for (auto &elem: tps_per_bx) {
    for (auto &mp_valid : elem.second) {
      if (!mp_valid.valid)
        continue;
      int isMatched = match(mp, coarsed, mp_valid);
      if (isMatched == 2) return true; // matched and quality <= store tp
    }
  }
  return false;
}

int MPCorFilter::killTps(cmsdt::metaPrimitive mp, std::vector<int> coarsed,
    int bx, std::map<int, valid_cor_tp_arr_t> &tps_per_bx) {
  int index_to_occupy = -1;
  int index_to_kill = -1;
  for (auto &elem: tps_per_bx) {
    if (abs(bx - elem.first) > 2) continue;
    for (size_t i = 0; i < elem.second.size(); i++) {
      if (elem.second[i].valid == 1) {
        int isMatched = match(mp, coarsed, elem.second[i]);
        if (isMatched == 3) {
          elem.second[i].valid = false;
          if (elem.first == bx && index_to_kill == -1) index_to_kill = i;
        }
      } else if (elem.first == bx && index_to_occupy == -1) index_to_occupy = i;
    }
  }
  // My first option is to replace the one from my BX that I killed first
  if (index_to_kill != -1) return index_to_kill;
  // If I wasn't able to kill anyone from my BX, I fill the first empty space
  return index_to_occupy;
}

void MPCorFilter::printmP(metaPrimitive mP) {
  DTSuperLayerId slId(mP.rawId);
  LogDebug("MPCorFilter") << slId << "\t"
                                      << " " << setw(2) << left << mP.wi1 << " " << setw(2) << left << mP.wi2 << " "
                                      << setw(2) << left << mP.wi3 << " " << setw(2) << left << mP.wi4 << " " << setw(5)
                                      << left << mP.tdc1 << " " << setw(5) << left << mP.tdc2 << " " << setw(5) << left
                                      << mP.tdc3 << " " << setw(5) << left << mP.tdc4 << " " << setw(10) << right
                                      << mP.x << " " << setw(9) << left << mP.tanPhi << " " << setw(5) << left << mP.t0
                                      << " " << setw(13) << left << mP.chi2;
}
