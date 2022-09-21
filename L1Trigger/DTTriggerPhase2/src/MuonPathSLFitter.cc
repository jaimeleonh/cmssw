#include "L1Trigger/DTTriggerPhase2/interface/MuonPathSLFitter.h"
#include <cmath>
#include <memory>

using namespace edm;
using namespace std;
using namespace cmsdt;
// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathSLFitter::MuonPathSLFitter(const ParameterSet &pset,
                                 edm::ConsumesCollector &iC,
                                 std::shared_ptr<GlobalCoordsObtainer> &globalcoordsobtainer)
    : MuonPathFitter(pset, iC, globalcoordsobtainer) {
  if (debug_)
    LogDebug("MuonPathSLFitter") << "MuonPathSLFitter: constructor";

  //shift theta
  int rawId;
  double shift;
  shift_theta_filename_ = pset.getParameter<edm::FileInPath>("shift_theta_filename");
  std::ifstream ifin4(shift_theta_filename_.fullPath());
  if (ifin4.fail()) {
    throw cms::Exception("Missing Input File")
        << "MuonPathSLFitter::MuonPathSLFitter() -  Cannot find " << shift_theta_filename_.fullPath();
  }

  while (ifin4.good()) {
    ifin4 >> rawId >> shift;
    shiftthetainfo_[rawId] = shift;
  }

  // LUTs
  sl1_filename_ = pset.getParameter<edm::FileInPath>("lut_sl1");
  sl3_filename_ = pset.getParameter<edm::FileInPath>("lut_sl3");

  fillLuts();
  
  setChi2Th(pset.getParameter<double>("chi2Th"));
  setTanPhiTh(pset.getParameter<double>("tanPhiTh"));
}

MuonPathSLFitter::~MuonPathSLFitter() {
  if (debug_)
    LogDebug("MuonPathSLFitter") << "MuonPathSLFitter: destructor";
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathSLFitter::initialise(const edm::EventSetup &iEventSetup) {
  if (debug_)
    LogDebug("MuonPathSLFitter") << "MuonPathSLFitter::initialiase";

  auto geom = iEventSetup.getHandle(dtGeomH);
  dtGeo_ = &(*geom);

}

void MuonPathSLFitter::run(edm::Event &iEvent,
                                   const edm::EventSetup &iEventSetup,
                                   MuonPathPtrs &muonpaths,
                                   std::vector<lat_vector>& lateralities,
                                   std::vector<metaPrimitive> &metaPrimitives) {
  if (debug_)
    LogDebug("MuonPathSLFitter") << "MuonPathSLFitter: run";

  // fit per SL (need to allow for multiple outputs for a single mpath)
  // for (auto &muonpath : muonpaths) {
  for (size_t i = 0; i < muonpaths.size(); i++) {
    // std::cout << "Starting with muon path " << i << std::endl;
    auto muonpath = muonpaths[i];
    auto lats = lateralities[i];
    analyze(muonpath, lats, metaPrimitives);
  }
}

void MuonPathSLFitter::finish() {
  if (debug_)
    LogDebug("MuonPathSLFitter") << "MuonPathSLFitter: finish";
};

//------------------------------------------------------------------
//--- Metodos privados
//------------------------------------------------------------------

void MuonPathSLFitter::analyze(MuonPathPtr &inMPath, lat_vector lat_combs, std::vector<cmsdt::metaPrimitive> &metaPrimitives) {
  auto sl = inMPath->primitive(0)->superLayerId(); // 0, 1, 2

  int selected_lay = 1;
  if (inMPath->primitive(0)->tdcTimeStamp() != -1)
    selected_lay = 0;

  int dumLayId = inMPath->primitive(selected_lay)->cameraId();
  auto dtDumlayerId = DTLayerId(dumLayId);
  DTSuperLayerId MuonPathSLId(dtDumlayerId.wheel(), dtDumlayerId.station(), dtDumlayerId.sector(), sl + 1);

  // if (MuonPathSLId.rawId() != 580788224)
    // return;

  DTChamberId ChId(MuonPathSLId.wheel(), MuonPathSLId.station(), MuonPathSLId.sector());
  // std::cout << "SL" << sl << std::endl;
  if (sl == 1)
    return;
  fit_common_in_t fit_common_in;

  // 8-element vectors, for the 8 layers. As here we are fitting one SL only, we leave the other SL values as dummy ones
  fit_common_in.hits = {};
  fit_common_in.hits_valid = {};

  int quality = 3;
  if (inMPath->missingLayer() != -1)
    quality = 1;

  for (int isl = 0; isl < 2; isl++) {
    for (int i = 0; i < NUM_LAYERS; i++) {
      if (isl * 2 == sl && inMPath->missingLayer() != i) {
        // Include both valid and non-valid hits. Non-valid values can be whatever, leaving all as -1 to make debugging easier.
        auto ti = inMPath->primitive(i)->tdcTimeStamp();
        if (ti != -1) ti = (int) round(((float) TIME_TO_TDC_COUNTS/ (float) LHC_CLK_FREQ) * ti);
        // std::cout << "ti: " << ti << std::endl;
        auto wi = inMPath->primitive(i)->channelId();
        auto ly = inMPath->primitive(i)->layerId();
        int layId = inMPath->primitive(i)->cameraId();
        auto dtlayerId = DTLayerId(layId);
        auto wireId = DTWireId(dtlayerId, wi + 1); // wire start from 1, mixer groups them starting from 0
        int rawId = wireId.rawId();
        // wp in tdc counts (still in floating point)
        float wp_f = ((10. * shiftinfo_[rawId] / CELL_SEMILENGTH) * MAXDRIFTTDC);
        // std::cout << "WPF: " << wp_f << " " <<  MAXDRIFTTDC << " " << shiftinfo_[rawId] << " " << (10. * shiftinfo_[rawId] / CELL_SEMILENGTH) << " " << rawId << std::endl;
        int wp = (int) ((long int)(round(wp_f * std::pow(2, WIREPOS_WIDTH))) / (int) std::pow(2, WIREPOS_WIDTH));
        // std::cout << "WP: " << wp << std::endl;
        fit_common_in.hits.push_back({ti, wi, ly, wp});
        // fill valids as well
        if (inMPath->missingLayer() == i) fit_common_in.hits_valid.push_back(0);
        else fit_common_in.hits_valid.push_back(1);
      } else {
        fit_common_in.hits.push_back({-1, -1, -1, -1});
        fit_common_in.hits_valid.push_back(0);
      }
    }
  }

  int smallest_time = 999999, tmp_coarse_wirepos_1 = -1, tmp_coarse_wirepos_3 = -1;
  // coarse_bctr is the 12 MSB of the smallest tdc
  for (int isl = 0; isl < 2; isl++) {
    if (isl * 2 != sl) continue;
    for (size_t i = 0; i < NUM_LAYERS; i++) {
      if (fit_common_in.hits_valid[NUM_LAYERS * isl + i] == 0) continue;
      else if (fit_common_in.hits[NUM_LAYERS * isl + i].ti < smallest_time)
        smallest_time = fit_common_in.hits[NUM_LAYERS * isl + i].ti;
    }
    if (fit_common_in.hits_valid[NUM_LAYERS * isl + 0] == 1)
      tmp_coarse_wirepos_1 = fit_common_in.hits[NUM_LAYERS * isl + 0].wp;
    else                                                     
      tmp_coarse_wirepos_1 = fit_common_in.hits[NUM_LAYERS * isl + 1].wp;
    if (fit_common_in.hits_valid[NUM_LAYERS * isl + 3] == 1)
      tmp_coarse_wirepos_3 = fit_common_in.hits[NUM_LAYERS * isl + 3].wp;
    else
      tmp_coarse_wirepos_3 = fit_common_in.hits[NUM_LAYERS * isl + 2].wp;

    tmp_coarse_wirepos_1 = tmp_coarse_wirepos_1 >> WIREPOS_NORM_LSB_IGNORED;
    tmp_coarse_wirepos_3 = tmp_coarse_wirepos_3 >> WIREPOS_NORM_LSB_IGNORED;
  }
  fit_common_in.coarse_bctr = smallest_time >> (WIDTH_FULL_TIME - WIDTH_COARSED_TIME);
  fit_common_in.coarse_wirepos = (tmp_coarse_wirepos_1 + tmp_coarse_wirepos_3) >> 1;

  for (auto &lat_comb : lat_combs) {
    if (lat_comb[0] == 0 && lat_comb[1] == 0 && lat_comb[2] == 0 && lat_comb[3] == 0)
      continue;
    fit_common_in.lateralities.clear();

    auto rom_addr = get_rom_addr(inMPath, lat_comb);
    coeffs_t coeffs;
    if (sl == 0) {
      coeffs = RomDataConvert(lut_sl1[rom_addr], COEFF_WIDTH_SL_T0, COEFF_WIDTH_SL_POSITION, COEFF_WIDTH_SL_SLOPE, 2 * sl, 2 * sl + 3);
    } else {
      coeffs = RomDataConvert(lut_sl3[rom_addr], COEFF_WIDTH_SL_T0, COEFF_WIDTH_SL_POSITION, COEFF_WIDTH_SL_SLOPE, 2 * sl, 2 * sl + 3);
    }

    // Filling lateralities
    for (int isl = 0; isl < 2; isl++) {
      for (size_t i = 0; i < NUM_LAYERS; i++) {
        if (isl * 2 == sl) {
          fit_common_in.lateralities.push_back(lat_comb[i]);
        }
        else fit_common_in.lateralities.push_back(-1);
      }
    }
    fit_common_in.coeffs = coeffs;
    // std::cout << "Starting to fit" << std::endl;
    // std::cout << inMPath->primitive(0)->channelId() << " ";
    // std::cout << inMPath->primitive(1)->channelId() << " ";
    // std::cout << inMPath->primitive(2)->channelId() << " ";
    // std::cout << inMPath->primitive(3)->channelId() << " ";
    // std::cout << inMPath->primitive(0)->tdcTimeStamp() << " ";
    // std::cout << inMPath->primitive(1)->tdcTimeStamp() << " ";
    // std::cout << inMPath->primitive(2)->tdcTimeStamp() << " ";
    // std::cout << inMPath->primitive(3)->tdcTimeStamp() << std::endl;

    auto fit_common_out = fit(fit_common_in,
                              XI_SL_WIDTH,
                              COEFF_WIDTH_SL_T0,
                              COEFF_WIDTH_SL_POSITION,
                              COEFF_WIDTH_SL_SLOPE,
                              PRECISSION_SL_T0,
                              PRECISSION_SL_POSITION,
                              PRECISSION_SL_SLOPE,
                              PROD_RESIZE_SL_T0,
                              PROD_RESIZE_SL_POSITION,
                              PROD_RESIZE_SL_SLOPE);
                              
    // std::cout << "Valid fit: " << fit_common_out.valid_fit << std::endl;
    if (fit_common_out.valid_fit == 1) {
      float t0_f = ((float) fit_common_out.t0) * (float) LHC_CLK_FREQ / (float) TIME_TO_TDC_COUNTS;
      float slope_f = -fit_common_out.slope * SLOPE_LSB;
      // std::cout << std::abs(slope_f) << " " << tanPhiTh_ << std::endl;
      if (std::abs(slope_f) > tanPhiTh_)
        return;

      // std::cout << "SLOPE: " << fit_common_out.slope  << " " << SLOPE_LSB << " " << slope_f << std::endl;
      // float pos_sl_f = ((float) (fit_common_out.position) + (sl - 1) * (fit_common_out.slope / 16.))
        // * ((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC);
      // std::cout << "POSITION: " << fit_common_out.position << " " << ((float) (fit_common_out.position) + (sl - 1) * (fit_common_out.slope / 16.)) << " " << ((float) (fit_common_out.position) + (sl - 1) * (fit_common_out.slope / 16.))
         // ((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC) << std::endl;
      // pos_sl_f /= 10.;
      float pos_ch_f = (float) (fit_common_out.position) * ((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC) / 10;
      float pos_sl_f = pos_ch_f - (sl - 1) * slope_f * VERT_PHI1_PHI3 / 2;

      float chi2_f = fit_common_out.chi2 * std::pow(((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC), 2) / 100;

      // obtention of global coordinates using luts
      DTWireId wireId(MuonPathSLId, 2, 1);      
      int pos = (int) (10 * (pos_sl_f - shiftinfo_[wireId.rawId()]) * INCREASED_RES_POS_POW);
      int slope = (int) (-slope_f * INCREASED_RES_SLOPE_POW);
      auto global_coords =
        globalcoordsobtainer_->get_global_coordinates(ChId.rawId(), sl + 1, pos, slope);
      float phi = global_coords[0];
      float phiB = global_coords[1];

      // obtention of global coordinates using cmssw geometry
      double z = 0;
      double z1 = Z_POS_SL;
      double z3 = -1. * z1;
      if (ChId.station() == 3 or ChId.station() == 4) {
        z1 = z1 + Z_SHIFT_MB4;
        z3 = z3 + Z_SHIFT_MB4;
      }
      if (MuonPathSLId.superLayer() == 1)
        z = z1;
      else if (MuonPathSLId.superLayer() == 3)
        z = z3;
      GlobalPoint jm_x_cmssw_global = dtGeo_->chamber(ChId)->toGlobal(LocalPoint(pos_sl_f, 0., z));
      int thisec = ChId.sector();
      if (thisec == 13)
        thisec = 4;
      if (thisec == 14)
        thisec = 10;
      float phi_cmssw = jm_x_cmssw_global.phi() - PHI_CONV * (thisec - 1);
      float psi = atan(slope_f);
      float phiB_cmssw = hasPosRF(ChId.wheel(), ChId.sector()) ? psi - phi_cmssw : -psi - phi_cmssw;

      metaPrimitives.emplace_back(metaPrimitive({MuonPathSLId.rawId(),
                                               t0_f,
                                               pos_sl_f,
                                               slope_f,
                                               phi,
                                               phiB,
                                               phi_cmssw,
                                               phiB_cmssw,
                                               chi2_f,
                                               quality,
                                               inMPath->primitive(0)->channelId(),
                                               inMPath->primitive(0)->tdcTimeStamp(),
                                               lat_comb[0],
                                               inMPath->primitive(1)->channelId(),
                                               inMPath->primitive(1)->tdcTimeStamp(),
                                               lat_comb[1],
                                               inMPath->primitive(2)->channelId(),
                                               inMPath->primitive(2)->tdcTimeStamp(),
                                               lat_comb[2],
                                               inMPath->primitive(3)->channelId(),
                                               inMPath->primitive(3)->tdcTimeStamp(),
                                               lat_comb[3],
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1,
                                               -1}));
    }
  }
  return;
}

void MuonPathSLFitter::fillLuts() {
  std::ifstream ifinsl1(sl1_filename_.fullPath());
  std::string line;
  while (ifinsl1.good()) {
    ifinsl1 >> line;

    std::vector<int> myNumbers;
    for (size_t i = 0; i < line.size(); i++) {
      // This converts the char into an int and pushes it into vec
      myNumbers.push_back(line[i] - '0');  // The digits will be in the same order as before
    }
    std::reverse(myNumbers.begin(), myNumbers.end());
    lut_sl1.push_back(myNumbers);
  }

  std::ifstream ifinsl3(sl3_filename_.fullPath());
  while (ifinsl3.good()) {
    ifinsl3 >> line;

    std::vector<int> myNumbers;
    for (size_t i = 0; i < line.size(); i++) {
      // This converts the char into an int and pushes it into vec
      myNumbers.push_back(line[i] - '0');  // The digits will be in the same order as before
    }
    std::reverse(myNumbers.begin(), myNumbers.end());
    lut_sl3.push_back(myNumbers);
  }

  return;  
}


int MuonPathSLFitter::get_rom_addr(MuonPathPtr &inMPath, latcomb lats) {
  /*
    vhdl code:
    rom_addr(5) <= reg.c1_input.segment.is4hit;
    if reg.c1_input.segment.is4hit = '1' then -- 4 layers fit
      rom_addr(4) <= '0';
      rom_addr(3 downto 0) <= reg.c1_input.segment.lateralities;
    else -- 3 layers fit
      rom_addr(4 downto 3) <= std_logic_vector(reg.c1_input.segment.missing_layer);
      rom_addr(2 downto 0) <= reg.c1_zeroSupprLats;
    end if;
  */
  std::vector<int> rom_addr;
  auto missing_layer = inMPath->missingLayer();
  if (missing_layer == -1) {
    rom_addr.push_back(1);
    rom_addr.push_back(0);
  } else {
    if (missing_layer == 0) {
      rom_addr.push_back(0); rom_addr.push_back(0);
    } else if (missing_layer == 1) {
      rom_addr.push_back(0); rom_addr.push_back(1);
    } else if (missing_layer == 2) {
      rom_addr.push_back(1); rom_addr.push_back(0);
    } else { // missing_layer == 3
      rom_addr.push_back(1); rom_addr.push_back(1);
    }
  }
  for (size_t ilat = 0; ilat < lats.size(); ilat++) {
    if ((int) ilat == missing_layer) // only applies to 3-hit, as in 4-hit missL=-1
      continue;
    auto lat = lats[ilat];
    if (lat == -1)
      lat = 0;
    rom_addr.push_back(lat);
  }
  std::reverse(rom_addr.begin(), rom_addr.end());
  return vhdl_unsigned_to_int(rom_addr);
}
