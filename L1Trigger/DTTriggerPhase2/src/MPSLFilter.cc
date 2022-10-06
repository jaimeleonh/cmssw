#include "L1Trigger/DTTriggerPhase2/interface/MPSLFilter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

using namespace edm;
using namespace std;
using namespace cmsdt;


// ============================================================================
// Constructors and destructor
// ============================================================================
MPSLFilter::MPSLFilter(const ParameterSet &pset)
    : MPFilter(pset), debug_(pset.getUntrackedParameter<bool>("debug")) {}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MPSLFilter::initialise(const edm::EventSetup &iEventSetup) {}

void MPSLFilter::run(edm::Event &iEvent,
                                  const edm::EventSetup &iEventSetup,
                                  std::vector<metaPrimitive> &inMPaths,
                                  std::vector<metaPrimitive> &outMPaths) {

  if (debug_)
    LogDebug("MPSLFilter") << "MPSLFilter: run";
  if (inMPaths.size() > 0) {
    int dum_sl_rawid = inMPaths[0].rawId;
    DTSuperLayerId dumSlId(dum_sl_rawid);
    DTChamberId ChId(dumSlId.wheel(), dumSlId.station(), dumSlId.sector());
    DTSuperLayerId sl1Id(ChId.rawId(), 1);
    DTSuperLayerId sl3Id(ChId.rawId(), 3);

    std::map<int, std::vector<metaPrimitive>> SL1metaPrimitivesPerBX;
    std::map<int, std::vector<metaPrimitive>> SL3metaPrimitivesPerBX;
    for (const auto &metaprimitiveIt : inMPaths) {
      int BX = metaprimitiveIt.t0 / 25;
      if (metaprimitiveIt.rawId == sl1Id.rawId())
        SL1metaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
      else if (metaprimitiveIt.rawId == sl3Id.rawId())
        SL3metaPrimitivesPerBX[BX].push_back(metaprimitiveIt);
    }

    auto filteredSL1MPs = filter(SL1metaPrimitivesPerBX);
    auto filteredSL3MPs = filter(SL3metaPrimitivesPerBX);

    for (auto & mp: filteredSL1MPs)
      outMPaths.push_back(mp);
    for (auto & mp: filteredSL3MPs)
      outMPaths.push_back(mp);
  }
}

void MPSLFilter::finish(){};

///////////////////////////
///  OTHER METHODS

std::vector<metaPrimitive> MPSLFilter::filter(std::map<int, std::vector<metaPrimitive>> mpsPerBX) {
  std::map<int, valid_tp_arr_t> mp_valid_per_bx;
  for (auto &elem: mpsPerBX) {
    mp_valid_per_bx[elem.first] = valid_tp_arr_t(6);
    for (auto &mp : elem.second) {
      // is this mp getting killed?
      if (isDead(mp, mp_valid_per_bx)) continue;
      // if not, let's kill other mps
      auto index = killTps(mp, elem.first, mp_valid_per_bx);
      if (index == -1) continue;
      mp_valid_per_bx[elem.first][index] = valid_tp_t({true, mp});
    }
  }

  std::vector<metaPrimitive> outTPs;
  for (auto &elem: mp_valid_per_bx) {
    for (auto &mp_valid : elem.second) {
      if (mp_valid.valid) outTPs.push_back(mp_valid.mp);
    }
  }
    
  return outTPs;
}


int MPSLFilter::match(cmsdt::metaPrimitive mp, cmsdt::metaPrimitive mp2) {
  if (mp.quality == LOWQ && mp2.quality == LOWQ)
    return 1;
  if (mp.quality != LOWQ && mp2.quality == LOWQ) {
    if (share_hit(mp, mp2)) return 2;
    return 3;
  }
  if (mp.quality == LOWQ && mp2.quality != LOWQ) {
    if (share_hit(mp, mp2)) return 4;
    return 5;
  }
  if (share_hit(mp, mp2)) {
    if (smaller_chi2(mp, mp2) == 0) return 6;
    return 7;
  }
  if (smaller_chi2(mp, mp2) == 0) return 8;
  return 9;
}

bool MPSLFilter::isDead(cmsdt::metaPrimitive mp, std::map<int, valid_tp_arr_t> tps_per_bx) {
  for (auto &elem: tps_per_bx) {
    for (auto &mp_valid : elem.second) {
      if (!mp_valid.valid)
        continue;
      int isMatched = match(mp, mp_valid.mp);
      if (isMatched == 4 || isMatched == 7) return true;
    }
  }
  return false;
}

int MPSLFilter::smaller_chi2(cmsdt::metaPrimitive mp, cmsdt::metaPrimitive mp2) {
  auto chi2_1 = get_chi2(mp);
  auto chi2_2 = get_chi2(mp2);
  if (chi2_1 < chi2_2) return 0;
  return 1;
}

int MPSLFilter::get_chi2(cmsdt::metaPrimitive mp) {
  // CHI2 is converted to an unsigned in which 4 msb are the exponent
  // of a float-like value and the rest of the bits are the mantissa 
  // (without the first 1). So comparing these reduced-width unsigned
  // values is equivalent to comparing rounded versions of the chi2

  int chi2 = (int) (mp.chi2 / (std::pow(((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC), 2) / 100));
  std::vector<int> chi2_unsigned, chi2_unsigned_msb;
  vhdl_int_to_unsigned(chi2, chi2_unsigned);

  if (chi2_unsigned.size() > 2) {
    for (int i = (int) chi2_unsigned.size() - 1; i >= 2; i--) {
      if (chi2_unsigned[i] == 1) {
        vhdl_int_to_unsigned(i - 1, chi2_unsigned_msb);
        for (int j = i - 1; j > i - 3; j--) {
          chi2_unsigned_msb.push_back(chi2_unsigned[j]);
        }
        return vhdl_unsigned_to_int(chi2_unsigned_msb);
      }
    }
  }
  vhdl_resize_unsigned(chi2_unsigned, 2);
  return vhdl_unsigned_to_int(vhdl_slice(chi2_unsigned, 1, 0));
}


int MPSLFilter::killTps(cmsdt::metaPrimitive mp, int bx, std::map<int, valid_tp_arr_t> &tps_per_bx) {
  int index_to_occupy = -1;
  int index_to_kill = -1;
  for (auto &elem: tps_per_bx) {
    if (abs(bx - elem.first) > 16) continue;
    for (size_t i = 0; i < elem.second.size(); i++) {
      if (elem.second[i].valid == 1) {
        int isMatched = match(mp, elem.second[i].mp);
        if (isMatched == 2 || isMatched == 6) {
          elem.second[i].valid = false;
          if (elem.first == bx && index_to_kill == -1) index_to_kill = i;
        }
      } else if (elem.first == bx && index_to_occupy == -1) index_to_occupy = i;
    }
  }
  // My first option is to replace the one from my BX that I killed first
  if (index_to_kill != -1) return index_to_kill;
  // If I wasn't able to kill anyone from my BX, I fill the first empty space
  if (index_to_occupy != -1) return index_to_occupy;
  // If I'm a 3h and there were no empty spaces, I don't replace any tp
  if (mp.quality == LOWQ) return -1;
  // If I'm a 4h, I replace the first 3h or the 4h with the biggest chi2.
  // Let's try to find both
  int biggest_chi2 = 0;
  for (size_t i = 0; i < tps_per_bx[bx].size(); i++) {
    if (tps_per_bx[bx][i].mp.quality == LOWQ) return i;
    auto chi2 = get_chi2(tps_per_bx[bx][i].mp);
    if (chi2 > biggest_chi2) {
      index_to_kill = i;
      biggest_chi2 = chi2;
    }
  }
  // If all stored tps are 4h and their chi2 is smaller than mine, I don't replace any
  if (biggest_chi2 < get_chi2(mp)) return -1;
  // If at least one chi2 is bigger than mine, I replace the corresponding tp
  return index_to_kill;
}

int MPSLFilter::share_hit(cmsdt::metaPrimitive mp, cmsdt::metaPrimitive mp2) {
  if (mp.rawId != mp2.rawId)
    return 0;
  if ((mp.wi1 == mp2.wi1 and mp.tdc1 == mp2.tdc1 and mp.wi1 != -1 and mp.tdc1 != -1) ||
      (mp.wi5 == mp2.wi5 and mp.tdc5 == mp2.tdc5 and mp.wi5 != -1 and mp.tdc5 != -1))
    return 1;
  if ((mp.wi2 == mp2.wi2 and mp.tdc2 == mp2.tdc2 and mp.wi2 != -1 and mp.tdc2 != -1) ||
  (mp.wi6 == mp2.wi6 and mp.tdc6 == mp2.tdc6 and mp.wi6 != -1 and mp.tdc6 != -1))
    return 2;
  if ((mp.wi3 == mp2.wi3 and mp.tdc3 == mp2.tdc3 and mp.wi3 != -1 and mp.tdc3 != -1) ||
      (mp.wi7 == mp2.wi7 and mp.tdc7 == mp2.tdc7 and mp.wi7 != -1 and mp.tdc7 != -1))
    return 3;
  if ((mp.wi4 == mp2.wi4 and mp.tdc4 == mp2.tdc4 and mp.wi4 != -1 and mp.tdc4 != -1) ||
      (mp.wi8 == mp2.wi8 and mp.tdc8 == mp2.tdc8 and mp.wi8 != -1 and mp.tdc8 != -1))
    return 4;
  return 0;
}

// int get_quality


void MPSLFilter::printmP(metaPrimitive mP) {
  DTSuperLayerId slId(mP.rawId);
  LogDebug("MPSLFilter") << slId << "\t"
                                      << " " << setw(2) << left << mP.wi1 << " " << setw(2) << left << mP.wi2 << " "
                                      << setw(2) << left << mP.wi3 << " " << setw(2) << left << mP.wi4 << " " << setw(5)
                                      << left << mP.tdc1 << " " << setw(5) << left << mP.tdc2 << " " << setw(5) << left
                                      << mP.tdc3 << " " << setw(5) << left << mP.tdc4 << " " << setw(10) << right
                                      << mP.x << " " << setw(9) << left << mP.tanPhi << " " << setw(5) << left << mP.t0
                                      << " " << setw(13) << left << mP.chi2;
}
