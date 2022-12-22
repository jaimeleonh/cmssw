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

  // std::map<int, std::vector<metaPrimitive>> SL1metaPrimitivesPerBX;
  std::vector<metaPrimitive> SL1metaPrimitives;
  // std::map<int, std::vector<metaPrimitive>> SL3metaPrimitivesPerBX;
  std::vector<metaPrimitive> SL3metaPrimitives;
  // std::map<int, std::vector<metaPrimitive>> CormetaPrimitivesPerBX;
  std::vector<metaPrimitive> CormetaPrimitives;
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
      // int BX = metaprimitiveIt.t0 / 25;
      if (metaprimitiveIt.rawId == sl1Id_rawid) {
        SL1metaPrimitives.push_back(metaprimitiveIt);
        // SL1metaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
      }
      else if (metaprimitiveIt.rawId == sl3Id_rawid)
        SL3metaPrimitives.push_back(metaprimitiveIt);
        // SL3metaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
    }
  }
  // for (const auto &metaprimitiveIt : inCorMPaths) {
    // int BX = metaprimitiveIt.t0 / 25;
    // std::cout << metaprimitiveIt.quality << " " << metaprimitiveIt.t0 << std::endl;
    // CormetaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
    // CormetaPrimitives.push_back(metaprimitiveIt);
  // }

  // auto filteredMPs = filter(SL1metaPrimitivesPerBX, SL3metaPrimitivesPerBX, CormetaPrimitivesPerBX);
  auto filteredMPs = filter(SL1metaPrimitives, SL3metaPrimitives, inCorMPaths);
  for (auto & mp: filteredMPs)
    outMPaths.push_back(mp);
}

void MPCorFilter::finish(){};

///////////////////////////
///  OTHER METHODS

std::vector<metaPrimitive> MPCorFilter::filter(
    std::vector<metaPrimitive> SL1mps,
    // std::map<int, std::vector<metaPrimitive>> SL1mpsPerBX,
    std::vector<metaPrimitive> SL3mps,
    // std::map<int, std::vector<metaPrimitive>> SL3mpsPerBX,
    std::vector<metaPrimitive> Cormps)
    // std::map<int, std::vector<metaPrimitive>> CormpsPerBX)
{
  // std::cout << "SL1 " << SL1mps.size() << " SL3 " << SL3mps.size() << std::endl;
  std::map<int, valid_cor_tp_arr_t> mp_valid_per_bx;
  // for (auto &elem: SL1mpsPerBX) {
  std::map<int, int> imp_per_bx_sl1;
  for (auto &mp: SL1mps) {
    int BX = mp.t0 / 25;
    // mp_valid_per_bx[elem.first] = valid_cor_tp_arr_t(12);
    if (mp_valid_per_bx.find(BX) == mp_valid_per_bx.end()) {
        mp_valid_per_bx[BX] = valid_cor_tp_arr_t(12);
    }
    
    if (imp_per_bx_sl1.find(BX) == imp_per_bx_sl1.end()) {
        imp_per_bx_sl1[BX] = 0;
    }
    
    // std::cout << "SL1 " << mp.t0 << std::endl;
    // int imp = 0;
    // for (auto &mp : elem.second) {
    auto coarsed = coarsify(mp, 1);
    // mp_valid_per_bx[elem.first][imp] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    mp_valid_per_bx[BX][imp_per_bx_sl1[BX]] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    // imp += 2;
    imp_per_bx_sl1[BX] += 2;
    // }
  }
  // for (auto &elem: SL3mpsPerBX) {
  std::map<int, int> imp_per_bx_sl3;
  for (auto &mp: SL3mps) {
    int BX = mp.t0 / 25;
    // if (mp_valid_per_bx.find(elem.first) == mp_valid_per_bx.end())
      // mp_valid_per_bx[elem.first] = valid_cor_tp_arr_t(12);
    if (mp_valid_per_bx.find(BX) == mp_valid_per_bx.end()) {
        mp_valid_per_bx[BX] = valid_cor_tp_arr_t(12);
    }

    if (imp_per_bx_sl3.find(BX) == imp_per_bx_sl3.end()) {
        imp_per_bx_sl3[BX] = 1;
    }

    // int imp = 1;
    // for (auto &mp : elem.second) {
    auto coarsed = coarsify(mp, 3);
    // mp_valid_per_bx[elem.first][imp] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    mp_valid_per_bx[BX][imp_per_bx_sl3[BX]] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    imp_per_bx_sl3[BX] += 2;
    // }
  }
  // for (auto &elem: CormpsPerBX) {
  for (auto &mp: Cormps) {
    // if (mp_valid_per_bx.find(elem.first) == mp_valid_per_bx.end()) {
      // mp_valid_per_bx[elem.first] = valid_cor_tp_arr_t(12);
    // }
    int BX = mp.t0 / 25;
    if (mp_valid_per_bx.find(BX) == mp_valid_per_bx.end()) {
      mp_valid_per_bx[BX] = valid_cor_tp_arr_t(12);
    }
    // for (auto &mp : elem.second) {
    auto coarsed = coarsify(mp, 0);
    auto dead = isDead(mp, coarsed, mp_valid_per_bx);
    // std::cout << mp.quality << " " << mp.t0 << " " << dead << std::endl;
    if (isDead(mp, coarsed, mp_valid_per_bx)) continue;
    // auto index = killTps(mp, coarsed, elem.first, mp_valid_per_bx);
    auto index = killTps(mp, coarsed, BX, mp_valid_per_bx);
    // mp_valid_per_bx[elem.first][index] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    mp_valid_per_bx[BX][index] = valid_cor_tp_t({true, mp, coarsed[3], coarsed[4], coarsed[5]});
    // }
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
  // float sign = 0;
  // if (sl == 1) sign = -1;
  // else if (sl == 3) sign = 1;

  // DTSuperLayerId MuonPathSLId(mp.rawId);  // SL1
  // DTChamberId ChId(MuonPathSLId.wheel(), MuonPathSLId.station(), MuonPathSLId.sector());  
  // DTSuperLayerId MuonPathSL1Id(ChId.wheel(), ChId.station(), ChId.sector(), 1);
  // DTWireId wireIdSL1(MuonPathSL1Id, 2, 1);

  // float pos_ch_f = mp.x + sign * mp.tanPhi * VERT_PHI1_PHI3 / 2;
  float pos_ch_f = mp.x;

  // translating into tdc counts
  // int pos_ch = int(round(pos_ch_f / (((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC) / 10)));
  int pos_ch = int(round(pos_ch_f));
  // int slope = (int) (-mp.tanPhi / (SLOPE_LSB * INCREASED_RES_SLOPE_POW));
  int slope = (int) (mp.tanPhi);

  // std::cout << sl << " " << mp.t0 << " " << pos_ch_f << " " << pos_ch << std::endl; 

  std::vector<int> t0_slv, t0_coarse, pos_slv, pos_coarse, slope_slv, slope_coarse;
  vhdl_int_to_unsigned(mp.t0, t0_slv);
  vhdl_int_to_signed(pos_ch, pos_slv);
  vhdl_int_to_signed(slope, slope_slv);

  // for (size_t i = 0; i < pos_slv.size(); i++)
    // std::cout << pos_slv[i];
  // std::cout << std::endl;

  vhdl_resize_unsigned(t0_slv, WIDTH_FULL_TIME);
  vhdl_resize_signed(pos_slv, WIDTH_FULL_POS);
  vhdl_resize_signed(slope_slv, WIDTH_FULL_SLOPE);

  // for (size_t i = 0; i < pos_slv.size(); i++)
    // std::cout << pos_slv[i];
  // std::cout << std::endl;

  t0_coarse = vhdl_slice(t0_slv, FSEG_T0_BX_LSB + 4, FSEG_T0_DISCARD_LSB);
  pos_coarse = vhdl_slice(pos_slv, WIDTH_FULL_POS - 1, FSEG_POS_DISCARD_LSB);
  slope_coarse = vhdl_slice(slope_slv, WIDTH_FULL_SLOPE - 1, FSEG_SLOPE_DISCARD_LSB);

  // for (size_t i = 0; i < pos_coarse.size(); i++)
    // std::cout << pos_coarse[i];
  // std::cout << std::endl;

  std::vector <int> results;
  int t0_coarse_int = vhdl_unsigned_to_int(t0_coarse);
  int pos_coarse_int = vhdl_signed_to_int(pos_coarse);
  int slope_coarse_int = vhdl_signed_to_int(slope_coarse);

  for (int index = 0; index <= 2; index++) {
    auto aux_t0_coarse_int = (t0_coarse_int + (index - 1))
      % (int) std::pow(2, FSEG_T0_BX_LSB + 4 - (FSEG_T0_DISCARD_LSB));
    auto aux_pos_coarse_int = pos_coarse_int + (index - 1);
    auto aux_slope_coarse_int = slope_coarse_int + (index - 1);
    // std::cout << "aux_pos_coarse_int " << aux_pos_coarse_int << " " << (aux_pos_coarse_int >> 1) << std::endl;
    results.push_back(aux_t0_coarse_int);
    results.push_back(aux_pos_coarse_int);
    results.push_back(aux_slope_coarse_int);
  }
  return results;  
}


int MPCorFilter::match(cmsdt::metaPrimitive mp, std::vector<int> coarsed, valid_cor_tp_t valid_cor_tp2) {
  // for (int i = 0; i < 3; i++)
    // for (int j = 0; j < 3; j++)
    // std::cout << coarsed[3 * i + j] << " ";
  // std::cout << std::endl;
  // std::cout << valid_cor_tp2.coarsed_t0 << " " << valid_cor_tp2.coarsed_pos << " " << valid_cor_tp2.coarsed_slope << std::endl;

  bool matched = (
    (coarsed[0] == valid_cor_tp2.coarsed_t0    || coarsed[3] == valid_cor_tp2.coarsed_t0    || coarsed[6] == valid_cor_tp2.coarsed_t0)  &&
    (coarsed[1] == valid_cor_tp2.coarsed_pos   || coarsed[4] == valid_cor_tp2.coarsed_pos   || coarsed[7] == valid_cor_tp2.coarsed_pos) &&
    (coarsed[2] == valid_cor_tp2.coarsed_slope || coarsed[5] == valid_cor_tp2.coarsed_slope || coarsed[8] == valid_cor_tp2.coarsed_slope)
  ) && (abs(mp.t0 / 25 - valid_cor_tp2.mp.t0 / 25) <= 1);
  // std::cout << mp.quality << " " << valid_cor_tp2.mp.quality << " " << mp.t0 << " " << valid_cor_tp2.mp.t0 << " " << matched << std::endl;
  return ((int) matched) * 2 +
    (int) (mp.quality > valid_cor_tp2.mp.quality) +
    (int) (mp.quality == valid_cor_tp2.mp.quality) * (int) (get_chi2(mp) < get_chi2(valid_cor_tp2.mp));
}

bool MPCorFilter::isDead(cmsdt::metaPrimitive mp, std::vector<int> coarsed, std::map<int, valid_cor_tp_arr_t> tps_per_bx) {
  for (auto &elem: tps_per_bx) {
    for (auto &mp_valid : elem.second) {
      if (!mp_valid.valid)
        continue;
      int isMatched = match(mp, coarsed, mp_valid);
      if (isMatched == 2) return true; // matched and quality <= stored tp
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
      // std::cout << elem.first << " " << i << " " << elem.second[i].mp.t0 << std::endl;
      if (elem.second[i].valid == 1) {
        int isMatched = match(mp, coarsed, elem.second[i]);
        // std::cout << i << " " << mp.t0 << " " << mp.quality << " " << elem.second[i].mp.t0 << " " << elem.second[i].mp.quality << " " << isMatched << std::endl;
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

int MPCorFilter::get_chi2(cmsdt::metaPrimitive mp) {
  // chi2 is coarsified to the index of the chi2's highest bit set to 1

  int chi2 = (int) round(mp.chi2 / (std::pow(((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC), 2) / 100));

  std::vector<int> chi2_unsigned, chi2_unsigned_msb;
  vhdl_int_to_unsigned(chi2, chi2_unsigned);
  // std::cout << mp.t0 << " " << chi2 << " ";
  // for (auto &elem: chi2_unsigned)
    // std::cout << elem;
  // std::cout << std::endl;

  for (int i = (int) chi2_unsigned.size() - 1; i >= 0; i--) {
    if (chi2_unsigned[i] == 1) {
      return i;
    }
  }
  return -1;
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
