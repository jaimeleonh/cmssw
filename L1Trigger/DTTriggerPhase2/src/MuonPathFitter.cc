#include "L1Trigger/DTTriggerPhase2/interface/MuonPathFitter.h"
#include <cmath>
#include <memory>

using namespace edm;
using namespace std;
using namespace cmsdt;
// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathFitter::MuonPathFitter(const ParameterSet &pset,
                                 edm::ConsumesCollector &iC,
                                 std::shared_ptr<GlobalCoordsObtainer> &globalcoordsobtainer)
    : MuonPathAnalyzer(pset, iC),
      // chi2Th_(pset.getParameter<double>("chi2Th")),
      tanPhiTh_(pset.getParameter<double>("tanPhiTh")),
      debug_(pset.getUntrackedParameter<bool>("debug")) {
  if (debug_)
    LogDebug("MuonPathFitter") << "MuonPathAnalyzer: constructor";

  //shift phi
  int rawId;
  shift_filename_ = pset.getParameter<edm::FileInPath>("shift_filename");
  std::ifstream ifin3(shift_filename_.fullPath());
  double shift;
  if (ifin3.fail()) {
    throw cms::Exception("Missing Input File")
        << "MuonPathFitter::MuonPathFitter() -  Cannot find " << shift_filename_.fullPath();
  }
  while (ifin3.good()) {
    ifin3 >> rawId >> shift;
    shiftinfo_[rawId] = shift;
  }

  //shift theta

  shift_theta_filename_ = pset.getParameter<edm::FileInPath>("shift_theta_filename");
  std::ifstream ifin4(shift_theta_filename_.fullPath());
  if (ifin4.fail()) {
    throw cms::Exception("Missing Input File")
        << "MuonPathAnalyzerPerSL::MuonPathAnalyzerPerSL() -  Cannot find " << shift_theta_filename_.fullPath();
  }

  while (ifin4.good()) {
    ifin4 >> rawId >> shift;
    shiftthetainfo_[rawId] = shift;
  }

  chosen_sl_ = pset.getParameter<int>("trigger_with_sl");

  if (chosen_sl_ != 1 && chosen_sl_ != 3 && chosen_sl_ != 4) {
    LogDebug("MuonPathFitter") << "chosen sl must be 1,3 or 4(both superlayers)";
    assert(chosen_sl_ != 1 && chosen_sl_ != 3 && chosen_sl_ != 4);  //4 means run using the two superlayers
  }

  dtGeomH = iC.esConsumes<DTGeometry, MuonGeometryRecord, edm::Transition::BeginRun>();
  globalcoordsobtainer_ = globalcoordsobtainer;

  // LUTs
  sl1_filename_ = pset.getParameter<edm::FileInPath>("lut_sl1");
  sl3_filename_ = pset.getParameter<edm::FileInPath>("lut_sl3");

  fillLuts();
}

MuonPathFitter::~MuonPathFitter() {
  if (debug_)
    LogDebug("MuonPathFitter") << "MuonPathAnalyzer: destructor";
}

// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathFitter::initialise(const edm::EventSetup &iEventSetup) {
  if (debug_)
    LogDebug("MuonPathFitter") << "MuonPathFitter::initialiase";

  auto geom = iEventSetup.getHandle(dtGeomH);
  dtGeo_ = &(*geom);

}

void MuonPathFitter::run(edm::Event &iEvent,
                                   const edm::EventSetup &iEventSetup,
                                   MuonPathPtrs &muonpaths,
                                   std::vector<lat_vector>& lateralities,
                                   std::vector<metaPrimitive> &metaPrimitives) {
  if (debug_)
    LogDebug("MuonPathFitter") << "MuonPathFitter: run";

  // fit per SL (need to allow for multiple outputs for a single mpath)
  // for (auto &muonpath : muonpaths) {
  for (size_t i = 0; i < muonpaths.size(); i++) {
    // std::cout << "Starting with muon path " << i << std::endl;
    auto muonpath = muonpaths[i];
    auto lats = lateralities[i];
    analyze(muonpath, lats, metaPrimitives);
  }
}

void MuonPathFitter::finish() {
  if (debug_)
    LogDebug("MuonPathFitter") << "MuonPathAnalyzer: finish";
};

//------------------------------------------------------------------
//--- Metodos privados
//------------------------------------------------------------------

void MuonPathFitter::analyze(MuonPathPtr &inMPath, lat_vector lat_combs, std::vector<cmsdt::metaPrimitive> &metaPrimitives) {
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

    auto fit_common_out = fit(inMPath,
                              fit_common_in,
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

fit_common_out_t MuonPathFitter::fit(MuonPathPtr &inMPath,
                                     fit_common_in_t fit_common_in,
                                     int XI_WIDTH,
                                     int COEFF_WIDTH_T0,
                                     int COEFF_WIDTH_POSITION,
                                     int COEFF_WIDTH_SLOPE,
                                     int PRECISSION_T0,
                                     int PRECISSION_POSITION,
                                     int PRECISSION_SLOPE,
                                     int PROD_RESIZE_T0,
                                     int PROD_RESIZE_POSITION,
                                     int PROD_RESIZE_SLOPE) {

  const int PARTIALS_PRECISSION = 4;
  // const int NORM_TIME_WIDTH = 10;
  const int PARTIALS_SHR_T0 = PRECISSION_T0 - PARTIALS_PRECISSION;
  const int PARTIALS_SHR_POSITION = PRECISSION_POSITION - PARTIALS_PRECISSION;
  const int PARTIALS_SHR_SLOPE = PRECISSION_SLOPE - PARTIALS_PRECISSION;
  const int PARTIALS_WIDTH_T0 = PROD_RESIZE_T0 - PARTIALS_SHR_T0;
  const int PARTIALS_WIDTH_POSITION = PROD_RESIZE_POSITION - PARTIALS_SHR_POSITION;
  const int PARTIALS_WIDTH_SLOPE = PROD_RESIZE_SLOPE - PARTIALS_SHR_SLOPE;

  const int WIDTH_TO_PREC = 11 + PARTIALS_PRECISSION;
  const int WIDTH_SLOPE_PREC = 14 + PARTIALS_PRECISSION;
  const int WIDTH_POSITION_PREC = WIDTH_SLOPE_PREC + 1;

  const int SEMICHAMBER_H_PRECISSION = 13 + PARTIALS_PRECISSION;
  // const int SEMICHAMBER_H_WIDTH = 2 + SEMICHAMBER_H_PRECISSION;
  const float SEMICHAMBER_H_REAL = ((235. / 2.) / (16. * 6.5)) * std::pow(2, SEMICHAMBER_H_PRECISSION);
  const int SEMICHAMBER_H = (int) SEMICHAMBER_H_REAL; // signed(SEMICHAMBER_H_WIDTH-1 downto 0)

  const int SEMICHAMBER_RES_SHR = SEMICHAMBER_H_PRECISSION;
  // const int SEMICHAMBER_RES_WIDTH = WIDTH_POSITION_PREC;
  // const int SEMICHAMBER_RES_RESIZE = SEMICHAMBER_RES_WIDTH + SEMICHAMBER_RES_SHR;

  const int LYRANDAHALF_RES_SHR = 4;
  // const int LYRANDAHALF_RES_WIDTH = WIDTH_SLOPE_PREC - 2;
  // const int LYRANDAHALF_RES_RESIZE = LYRANDAHALF_RES_WIDTH + LYRANDAHALF_RES_SHR;

  const int CHI2_CALC_RES_BITS = 7;

  /*******************************
            clock cycle 1
  *******************************/
  std::vector<int> normalized_times;
  std::vector<int> normalized_wirepos;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    // std::cout << i << std::endl;
    // normalized times
    // this should be resized to an unsigned of 10 bits (max drift time ~508 TDC counts, using 9+1 to include tolerance)
    // leaving it as an integer for now
    // we are obtaining the difference as the difference in BX + the LS bits from the hit time

    if (fit_common_in.hits_valid[i] == 1) {
      // std::cout << fit_common_in.hits[i].ti << " " << fit_common_in.coarse_bctr << std::endl;
      int dif_bx = (fit_common_in.hits[i].ti >> (WIDTH_FULL_TIME - WIDTH_COARSED_TIME))
        - fit_common_in.coarse_bctr;

      int tmp_norm_time = (dif_bx << (WIDTH_FULL_TIME - WIDTH_COARSED_TIME)) +
        (fit_common_in.hits[i].ti % (int) std::pow(2, WIDTH_FULL_TIME - WIDTH_COARSED_TIME));
      // std::cout << dif_bx << " " << tmp_norm_time << std::endl;
      // resize test
      // this has implications in the FW (reducing number of bits).
      // we keep here the int as it is, but we do the same check done in the fw
      std::vector<int> tmp_dif_bx_vector;
      vhdl_int_to_unsigned(dif_bx, tmp_dif_bx_vector);
      // for (auto & elem: tmp_dif_bx_vector) {
        // std::cout << elem;
      // }
      // std::cout << std::endl;
      vhdl_resize_unsigned(tmp_dif_bx_vector, 12);
      // for (auto & elem: tmp_dif_bx_vector) {
        // std::cout << elem;
      // }
      // std::cout << std::endl;
      if (!vhdl_resize_unsigned_ok(tmp_dif_bx_vector, WIDTH_DIFBX))
        return fit_common_out_t();
      // std::cout << "After resizing tmp_dif_bx_vector" << std::endl;

      normalized_times.push_back(tmp_norm_time);
      // std::cout << "normalized_times[" << i << "]=" << normalized_times[i] << std::endl;
      int tmp_wirepos = fit_common_in.hits[i].wp - 
        (fit_common_in.coarse_wirepos << WIREPOS_NORM_LSB_IGNORED);
      // std::cout << fit_common_in.hits[i].wp << " " << fit_common_in.coarse_wirepos << " " << tmp_wirepos << std::endl;
      // resize test
      std::vector<int> tmp_wirepos_vector;
      vhdl_int_to_signed(tmp_wirepos, tmp_wirepos_vector);
      // for (auto & elem: tmp_wirepos_vector) {
        // std::cout << elem;
      // }
      // std::cout << std::endl;
      vhdl_resize_signed(tmp_wirepos_vector, WIREPOS_WIDTH);
      // for (auto & elem: tmp_wirepos_vector) {
        // std::cout << elem;
      // }
      // std::cout << std::endl;
      if (!vhdl_resize_signed_ok(tmp_wirepos_vector, XI_WIDTH))
        return fit_common_out_t();
      // std::cout << "After resizing tmp_wirepos_vector" << std::endl;

      normalized_wirepos.push_back(tmp_wirepos);
    } else { // dummy hit
      normalized_times.push_back(-1);
      normalized_wirepos.push_back(-1);
    }
  }

  // std::cout << "Clock cycle 1 finished" << std::endl;

  /*******************************
            clock cycle 2
  *******************************/

  std::vector<int> xi_arr;
  // min and max times are computed throught several clk cycles in the fw, 
  // here we compute it at once
  int min_hit_time = 999999, max_hit_time = 0;

  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 1) {
      // calculate xi array
      auto tmp_xi_incr = normalized_wirepos[i];
      tmp_xi_incr += (-1 + 2 * fit_common_in.lateralities[i]) * normalized_times[i];

      // resize test
      std::vector<int> tmp_xi_incr_vector;
      vhdl_int_to_signed(tmp_xi_incr, tmp_xi_incr_vector);
      vhdl_resize_signed(tmp_xi_incr_vector, XI_WIDTH + 1);
      if (!vhdl_resize_signed_ok(tmp_xi_incr_vector, XI_WIDTH))
        return fit_common_out_t();
      xi_arr.push_back(tmp_xi_incr);

      // std::cout << "xi_arr[" << i << "]=" << xi_arr[i] << std::endl;

      // calculate min and max times
      if (normalized_times[i] < min_hit_time) {
        min_hit_time = normalized_times[i];
      }
      if (normalized_times[i] > max_hit_time) {
        max_hit_time = normalized_times[i];
      }
    } else {
      xi_arr.push_back(-1);
    }
  }

  // std::cout << "Clock cycle 2 finished" << std::endl;
  /*******************************
            clock cycle 3
  *******************************/

  std::vector <int> products_t0;
  std::vector <int> products_position;
  std::vector <int> products_slope;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 0) {
      products_t0.push_back(       -1);
      products_position.push_back( -1);
      products_slope.push_back(    -1);
    } else {
      // std::cout << "time coeff ";
      // for (auto & elem: fit_common_in.coeffs.t0       [i])
        // std::cout << elem;

      // std::cout << " " << vhdl_signed_to_int(fit_common_in.coeffs.t0       [i]) << std::endl;
      products_t0.push_back(       xi_arr[i] * vhdl_signed_to_int(fit_common_in.coeffs.t0       [i]));
      products_position.push_back( xi_arr[i] * vhdl_signed_to_int(fit_common_in.coeffs.position [i]));
      products_slope.push_back(    xi_arr[i] * vhdl_signed_to_int(fit_common_in.coeffs.slope    [i]));
    }
  }

  // std::cout << "Clock cycle 3 finished" << std::endl;
  /*******************************
            clock cycle 4
  *******************************/
  // Do the 8 element sums
  int t0_prec = 0, position_prec = 0, slope_prec = 0;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 0) {
      continue;
    } else {
      t0_prec       += products_t0[i] >> PARTIALS_SHR_T0;
      position_prec += products_position[i] >> PARTIALS_SHR_POSITION;
      slope_prec    += products_slope[i] >> PARTIALS_SHR_SLOPE;
    }
  }

  // std::cout << "Clock cycle 4 finished" << std::endl;
  /*******************************
            clock cycle 5
  *******************************/
  // Do resize tests for the computed sums with full precision
  std::vector<int> t0_prec_vector, position_prec_vector, slope_prec_vector;
  // std::cout << "T0" << std::endl;
  vhdl_int_to_signed(t0_prec, t0_prec_vector);
  // std::cout << t0_prec << " ";
  // for (auto & elem: t0_prec_vector) {
    // std::cout << elem;
  // }
  // std::cout << " ";

  vhdl_resize_signed(t0_prec_vector, PARTIALS_WIDTH_T0);
  // for (auto & elem: t0_prec_vector) {
    // std::cout << elem;
  // }
  // std::cout << std::endl;
  if (!vhdl_resize_signed_ok(t0_prec_vector, WIDTH_TO_PREC))
    return fit_common_out_t();

  vhdl_int_to_signed(position_prec, position_prec_vector);
  vhdl_resize_signed(position_prec_vector, PARTIALS_WIDTH_POSITION);
  // std::cout << "Position" << std::endl;
  // std::cout << position_prec << " ";
  if (!vhdl_resize_signed_ok(position_prec_vector, WIDTH_POSITION_PREC))
    return fit_common_out_t();

  vhdl_int_to_signed(slope_prec, slope_prec_vector);
  vhdl_resize_signed(slope_prec_vector, PARTIALS_WIDTH_SLOPE);
  // std::cout << "Slope" << std::endl;
  // std::cout << slope_prec << " ";
  if (!vhdl_resize_signed_ok(slope_prec_vector, WIDTH_SLOPE_PREC))
    return fit_common_out_t();

  // std::cout << "Clock cycle 5 finished" << std::endl;
  /*******************************
            clock cycle 6
  *******************************/
  // Round the fitting parameters to the final resolution;
  // in vhdl something more sofisticated is done, here we do a float division, round
  // and cast again to integer
  int norm_t0 = (int)(round(t0_prec / std::pow(2, PARTIALS_PRECISSION)));
  int norm_position = (int)(round((float) position_prec / std::pow(2, PARTIALS_PRECISSION)));
  int norm_slope = (int)(round((float) slope_prec / std::pow(2, PARTIALS_PRECISSION)));

  // std::cout << "normt0 " << norm_t0 << " norm_position " << norm_position << " norm_slope " << norm_slope << std::endl; 

  // Calculate the (-xi) + pos (+/-) t0, which only is lacking the slope term to become the residuals
  std::vector<int> res_partials_arr;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 0) {
      res_partials_arr.push_back(-1);
    } else {
      int tmp_position_prec = position_prec - (xi_arr[i] << PARTIALS_PRECISSION);
      // rounding
      tmp_position_prec += std::pow(2, PARTIALS_PRECISSION - 1);

      // std::cout << "position_prec=" << position_prec;
      // std::cout << " xi_arr[i]=" << xi_arr[i] * (int) std::pow(2, PARTIALS_PRECISSION);
      // std::cout << " tmp_position_prec[" << i << "]=" << tmp_position_prec << std::endl; 

      tmp_position_prec += (-1 + 2 * fit_common_in.lateralities[i]) * t0_prec;
      res_partials_arr.push_back(tmp_position_prec);
      // std::cout << "c6.res_partials_arr[" << i << "]=" << res_partials_arr[i] << std::endl;
    }
  }

  // calculate the { slope x semichamber, slope x 1.5 layers, slope x 0.5 layers }
  // these 3 values are later combined with different signs to get the slope part
  // of the residual for each of the layers.
  // std::cout << slope_prec << " " << SEMICHAMBER_H << " " << (slope_prec * SEMICHAMBER_H) << " " << ((slope_prec * SEMICHAMBER_H)>> SEMICHAMBER_RES_SHR) << std::endl;
  int slope_x_halfchamb = (((long int)slope_prec * (long int) SEMICHAMBER_H)) >> SEMICHAMBER_RES_SHR;
  int slope_x_3semicells = (slope_prec * 3) >> LYRANDAHALF_RES_SHR;
  int slope_x_1semicell = (slope_prec * 1) >> LYRANDAHALF_RES_SHR;

  // std::cout << "slope times stuff: " <<  SEMICHAMBER_H << " " << slope_x_halfchamb << " " << slope_x_3semicells << " " << slope_x_1semicell << std::endl;

  // std::cout << "Clock cycle 6 finished" << std::endl;
  /*******************************
            clock cycle 7
  *******************************/
  // Complete the residuals calculation by constructing the slope term (1/2)
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 1) {
      if      (i % 4 == 0) res_partials_arr[i] -= slope_x_3semicells;
      else if (i % 4 == 1) res_partials_arr[i] -= slope_x_1semicell;
      else if (i % 4 == 2) res_partials_arr[i] += slope_x_1semicell;
      else                 res_partials_arr[i] += slope_x_3semicells;
      // std::cout << "c7.res_partials_arr[" << i << "]=" << res_partials_arr[i] << std::endl;
    }
  }

  // std::cout << "Clock cycle 7 finished" << std::endl;
  /*******************************
            clock cycle 8
  *******************************/
  // Complete the residuals calculation by constructing the slope term (2/2)
  std::vector<int> residuals, position_prec_arr;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 0) {
      residuals.push_back(-1);
      position_prec_arr.push_back(-1);
    } else {
      int tmp_position_prec = res_partials_arr[i];
      // std::cout << "tmp_position_prec[" << i << "]=" << tmp_position_prec << std::endl;
      tmp_position_prec += (-1 + 2 * (int)(i >= NUM_LAYERS)) * slope_x_halfchamb;
      // std::cout << "tmp_position_prec[" << i << "]=" << tmp_position_prec << std::endl;
      position_prec_arr.push_back(tmp_position_prec);
      residuals.push_back(abs(tmp_position_prec >> PARTIALS_PRECISSION));
      // std::cout << "residuals[" << i << "]=" << residuals[i] << std::endl;
    }
  }

  // minimum and maximum fit t0
  int min_t0 = max_hit_time - MAXDRIFTTDC - T0_CUT_TOLERANCE;
  int max_t0 = min_hit_time + T0_CUT_TOLERANCE;

  // std::cout << "Clock cycle 8 finished" << std::endl;
  /*******************************
            clock cycle 9
  *******************************/
  // Prepare addition of coarse_offset to T0 (T0 de-normalization)
  int t0_fine = norm_t0 & (int) (std::pow(2, 5) - 1);
  int t0_bx_sign = ((int) (norm_t0 < 0)) * 1;
  int t0_bx_abs = abs(norm_t0 >> 5);

  // De-normalize Position and slope
  int position = (fit_common_in.coarse_wirepos << WIREPOS_NORM_LSB_IGNORED) + norm_position;
  int slope = norm_slope;

  // std::cout << "t0_fine=" << t0_fine << " " ;
  // std::cout << "t0_bx_sign=" << t0_bx_sign << " ";
  // std::cout << "t0_bx_abs=" << t0_bx_abs << " ";
  // std::cout << "position=" << position << " ";
  // std::cout << "slope=" << slope << " " << std::endl;

  // Apply T0 cuts
  // std::cout << norm_t0 << " " << min_t0 << std::endl;
  if (norm_t0 < min_t0) return fit_common_out_t();
  // std::cout << norm_t0 << " " << max_t0 << std::endl;
  if (norm_t0 > max_t0) return fit_common_out_t();

  // double slope_f = -(double(slope) / INCREASED_RES_SLOPE_POW);
  // if (std::abs(slope_f) > tanPhiTh_)
    // return;

  // std::cout << slope << " " << (506 / 2) * 16 << std::endl;
  // if (abs(slope) > (506 / 2) * 16)
    // return fit_common_out_t();

  // square the residuals
  std::vector<int> squared_residuals;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 0) {
      squared_residuals.push_back(-1);
    } else {
      squared_residuals.push_back(residuals[i] * residuals[i]);
    }
  }

  // check for residuals overflow
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 1) {
      std::vector<int> tmp_vector;
      int tmp_position_prec = (position_prec_arr[i] >> PARTIALS_PRECISSION);
      vhdl_int_to_signed(tmp_position_prec, tmp_vector);
      vhdl_resize_signed(tmp_vector, WIDTH_POSITION_PREC);
      // for (auto & elem: tmp_vector)
        // std::cout << elem;
      // std::cout << " " << CHI2_CALC_RES_BITS + 1 << std::endl;
      if (!vhdl_resize_signed_ok(tmp_vector, CHI2_CALC_RES_BITS + 1))
        return fit_common_out_t();
      // Commented for now, maybe later we need to do something here
      // if ((tmp_position_prec / (int) std::pow(2, CHI2_CALC_RES_BITS)) > 0)
        // return fit_common_out_t();
    }
  }

  // std::cout << "Clock cycle 9 finished" << std::endl;
  /*******************************
        clock cycle 10, 11, 12
  *******************************/
  int t0 = t0_fine;
  t0 += (fit_common_in.coarse_bctr - (- 1 + 2 * t0_bx_sign) * t0_bx_abs) * (int) std::pow(2, 5);

  int chi2 = 0;
  for (int i = 0; i < 2 * NUM_LAYERS; i++) {
    if (fit_common_in.hits_valid[i] == 1) {
      chi2 += squared_residuals[i];
    }
  }

  // Impose the thresholds
  // if (chi2 > 16 * 16)
  // std::cout << "chi2 " << chi2 << std::endl;
  if (chi2 > (0.01 / (std::pow(((float) CELL_SEMILENGTH / (float) MAXDRIFTTDC), 2) / 100))) // FIXME
    return fit_common_out_t();

  // double chi2_f = double(chi2) / (16. * 64. * 100.);

  // std::cout << "Final position: " << position << std::endl;
  // std::cout << "Final slope: " << slope << std::endl;
  // std::cout << "Final t0: " << t0 << std::endl;
  // std::cout << "Final chi2: " << chi2 << std::endl;

  fit_common_out_t fit_common_out;
  fit_common_out.position = position;
  fit_common_out.slope = slope;
  fit_common_out.t0 = t0;
  fit_common_out.chi2 = chi2;
  fit_common_out.valid_fit = 1;

  // std::cout << "Clock cycle 10,11,12 finished" << std::endl;

  return fit_common_out;  
}



void MuonPathFitter::fillLuts() {
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

coeffs_t MuonPathFitter::RomDataConvert(std::vector<int> slv, short COEFF_WIDTH_T0, short COEFF_WIDTH_POSITION, short COEFF_WIDTH_SLOPE, short LOLY, short HILY) {
  coeffs_t res;
  int ctr = 0;
  for (int i = LOLY; i <= HILY; i++) {
    res.t0[i] = vhdl_slice(slv, COEFF_WIDTH_T0 + ctr - 1, ctr);
    vhdl_resize_unsigned(res.t0[i], GENERIC_COEFF_WIDTH);
    res.t0[i] = vhdl_slice(res.t0[i], COEFF_WIDTH_T0 - 1, 0);
    ctr += COEFF_WIDTH_T0;
  }
  for (int i = LOLY; i <= HILY; i++) {
    res.position[i] = vhdl_slice(slv, COEFF_WIDTH_POSITION + ctr - 1, ctr);
    vhdl_resize_unsigned(res.position[i], GENERIC_COEFF_WIDTH);
    res.position[i] = vhdl_slice(res.position[i], COEFF_WIDTH_POSITION - 1, 0);
    ctr += COEFF_WIDTH_POSITION;
  }
  for (int i = LOLY; i <= HILY; i++) {
    res.slope[i] = vhdl_slice(slv, COEFF_WIDTH_SLOPE + ctr - 1, ctr);
    vhdl_resize_unsigned(res.slope[i], GENERIC_COEFF_WIDTH);
    res.slope[i] = vhdl_slice(res.slope[i], COEFF_WIDTH_SLOPE - 1, 0);
    ctr += COEFF_WIDTH_SLOPE;
  }
  return res;
}


int MuonPathFitter::get_rom_addr(MuonPathPtr &inMPath, latcomb lats) {
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
