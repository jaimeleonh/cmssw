#ifndef L1Trigger_DTTriggerPhase2_MuonPathFitter_h
#define L1Trigger_DTTriggerPhase2_MuonPathFitter_h

#include "L1Trigger/DTTriggerPhase2/interface/MuonPathAnalyzer.h"

// ===============================================================================
// Previous definitions and declarations
// ===============================================================================

using coeff_arr_t = std::vector<std::vector<int>>;
struct coeffs_t {
  coeff_arr_t t0;
  coeff_arr_t position;
  coeff_arr_t slope;
  coeffs_t():
    t0(cmsdt::N_COEFFS, std::vector<int>(cmsdt::GENERIC_COEFF_WIDTH, 0)),
    position(cmsdt::N_COEFFS, std::vector<int>(cmsdt::GENERIC_COEFF_WIDTH, 0)),
    slope(cmsdt::N_COEFFS, std::vector<int> (cmsdt::GENERIC_COEFF_WIDTH, 0)) {}
};

struct SLhitP {
  int ti; // unsigned(16 downto 0); -- 12 msb = bunch_ctr, 5 lsb = tdc counts, resolution 25/32 ns
  int wi; // unsigned(6 downto 0); -- ~ 96 channels per layer
  int ly; // unsigned(1 downto 0); -- 4 layers
  int wp; // signed(WIREPOS_WIDTH-1 downto 0);
};

struct fit_common_in_t {
  // int valid; not needed, we will not propagate the mpath to the fitter
  std::vector<SLhitP> hits;
  std::vector<int> hits_valid; // slv(0 to 7)
  std::vector<int> lateralities; // slv(0 to 7)
  coeffs_t coeffs;
  int coarse_bctr; // unsigned(11 downto 0)
  int coarse_wirepos; // signed(WIDTH_FULL_POS-1 downto WIREPOS_NORM_LSB_IGNORED);
};

struct fit_common_out_t {
  int t0;
  int slope;
  int position;
  int chi2;
  int valid_fit;
  fit_common_out_t(): t0(0), slope(0), position(0), chi2(0), valid_fit(0) {}
};


// "à la vhdl" functions

std::vector<int> vhdl_slice(std::vector<int> v, int upper, int lower)
{
    int final_value = lower;
    if (final_value < 0) final_value = 0;

    std::vector<int> v1;
    for (int i = final_value; i <= upper; i++) {
        v1.push_back(v[i]);
    }
    return v1;
}

int vhdl_unsigned_to_int(std::vector<int> v) {
  int res = 0;
  
  for (size_t i = 0; i < v.size(); i++) {
    res = res + v[i] * std::pow(2, i);
  }
  return res;
}

int vhdl_signed_to_int(std::vector<int> v) {
  if (v[v.size() - 1] == 0) return vhdl_unsigned_to_int(v);
  else return -(std::pow(2, v.size()) - vhdl_unsigned_to_int(v));
}

void vhdl_int_to_unsigned(int value, std::vector<int> &v) {
  if (value == 0) {
    v.push_back(0);
  } else if (value != 1) {
    v.push_back(value % 2);
    vhdl_int_to_unsigned(value / 2, v);
  } else {
    v.push_back(1);
  }
  return;
}

void vhdl_int_to_signed(int value, std::vector<int> &v) {
  if (value < 0) {
    int val = 1;
    while (val < value) {
      val *= 2;
    }
    vhdl_int_to_unsigned(val - value, v);
    v.push_back(0);
    v.push_back(1);
  } else {
    vhdl_int_to_unsigned(value, v);
  }
  return;
}

void vhdl_resize_unsigned(std::vector<int> &v, int new_size) {
  for (int i = v.size(); i < new_size; i++) {
    v.push_back(0);
  }
}

void vhdl_resize_signed(std::vector<int> &v, int new_size) {
  int elem = 0;
  if (v[v.size() - 1] == 1) elem = 1;
  for (int i = v.size(); i < new_size; i++) {
    v.push_back(elem);
  }
}


bool vhdl_resize_signed_ok(std::vector<int> v, int new_size) {
  for (size_t i = v.size() - 1 - 1; i >= v.size() - 1 - (v.size() - new_size); i--) {
    if (v[i] != v[v.size() - 1]) return false;
  }
  return true;
};


bool vhdl_resize_unsigned_ok(std::vector<int> v, int new_size) {
  for (size_t i = v.size() - 1; i >= v.size() - 1 + 1 - (v.size() - new_size); i--) {
    if (v[i] != 0) return false;
  }
  return true;
};


// ===============================================================================
// Class declarations
// ===============================================================================


class MuonPathFitter : public MuonPathAnalyzer {
public:
  // Constructors and destructor
  MuonPathFitter(const edm::ParameterSet &pset,
                           edm::ConsumesCollector &iC,
                           std::shared_ptr<GlobalCoordsObtainer> &globalcoordsobtainer);
  ~MuonPathFitter() override;

  // Main methods
  void initialise(const edm::EventSetup &iEventSetup) override;
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           MuonPathPtrs &inMpath,
           std::vector<cmsdt::metaPrimitive> &metaPrimitives) override{};
  virtual void run(edm::Event& iEvent,
           const edm::EventSetup& iEventSetup,
           MuonPathPtrs& inMpath,
           std::vector<lat_vector>& lateralities,
           std::vector<cmsdt::metaPrimitive>& metaPrimitives) override;
  void run(edm::Event &iEvent,
           const edm::EventSetup &iEventSetup,
           MuonPathPtrs &inMpath,
           MuonPathPtrs &outMPath) override{};

  void finish() override;

  // Other public methods

  bool hasPosRF(int wh, int sec) { return wh > 0 || (wh == 0 && sec % 4 > 1); };

  // Public attributes
  DTGeometry const *dtGeo_;
  edm::ESGetToken<DTGeometry, MuonGeometryRecord> dtGeomH;

  //shift
  edm::FileInPath shift_filename_;
  std::map<int, float> shiftinfo_;

  //shift theta
  edm::FileInPath shift_theta_filename_;
  std::map<int, float> shiftthetainfo_;

  // luts
  edm::FileInPath sl1_filename_;
  edm::FileInPath sl3_filename_;

  int chosen_sl_;

private:
  // Private methods
  void analyze(MuonPathPtr &inMPath, lat_vector lat_combs, std::vector<cmsdt::metaPrimitive> &metaPrimitives);
  void fillLuts();
  coeffs_t RomDataConvert(std::vector<int> slv, short COEFF_WIDTH_T0, short COEFF_WIDTH_POSITION, short COEFF_WIDTH_SLOPE, short LOLY, short HILY);
  int get_rom_addr(MuonPathPtr &inMPath, latcomb lats);
  fit_common_out_t fit(MuonPathPtr &inMPath,
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
                       int PROD_RESIZE_SLOPE);

  // Private attributes
  double tanPhiTh_;
  const bool debug_;
  // double chi2Th_;
  std::vector<std::vector<int>> lut_sl1;
  std::vector<std::vector<int>> lut_sl3;

  // global coordinates
  std::shared_ptr<GlobalCoordsObtainer> globalcoordsobtainer_;
};

#endif
