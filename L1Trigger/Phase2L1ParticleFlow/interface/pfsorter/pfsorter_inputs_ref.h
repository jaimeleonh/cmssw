#ifndef pfsorter_inputs_ref_h
#define pfsorter_inputs_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"

#include <vector>

namespace l1ct {
  namespace pfsorter {

    // One particle coming out of the ParticleFlow algorithm, as seen on one input link of
    // the sorter. PF particles come in two flavours, charged (l1ct::PFChargedObjEmu, also
    // used for muons) and neutral (l1ct::PFNeutralObjEmu, also used for photons), and their
    // coordinates are local to the region they were reconstructed in; the sorter tree
    // instead stores l1ct::PuppiObjEmu, in global coordinates, so every particle has to be
    // converted with toPuppi() and the region it comes from before entering the tree.
    struct PFParticleEmu {
      enum Kind { None = 0, Charged = 1, Neutral = 2 };

      Kind kind;
      l1ct::PFChargedObjEmu charged;
      l1ct::PFNeutralObjEmu neutral;
      l1ct::puppiWgt_t puppiWgt;  // weight of a neutral particle (1 = no puppi correction)
      l1ct::pt_t puppiPt;         // pt of a neutral particle after the weight

      PFParticleEmu() { clear(); }
      explicit PFParticleEmu(const l1ct::PFChargedObjEmu &src) { set(src); }
      explicit PFParticleEmu(const l1ct::PFNeutralObjEmu &src, l1ct::puppiWgt_t wgt = 1.0) { set(src, wgt); }

      void clear() {
        kind = None;
        charged.clear();
        neutral.clear();
        puppiWgt = 0;
        puppiPt = 0;
      }

      void set(const l1ct::PFChargedObjEmu &src) {
        clear();
        if (src.hwPt == 0)
          return;
        kind = Charged;
        charged = src;
      }

      void set(const l1ct::PFNeutralObjEmu &src, l1ct::puppiWgt_t wgt = 1.0) {
        clear();
        if (src.hwPt == 0)
          return;
        kind = Neutral;
        neutral = src;
        puppiWgt = wgt;
        puppiPt = l1ct::pt_t(src.hwPt * wgt);
      }

      // pt the particle will have once converted (0 for an empty link)
      l1ct::pt_t hwPt() const {
        if (kind == Charged)
          return charged.hwPt;
        else if (kind == Neutral)
          return puppiPt;
        return l1ct::pt_t(0);
      }
      bool valid() const { return kind != None && hwPt() != 0; }
      int intPt() const { return l1ct::Scales::intPt(hwPt()); }
    };

    // convert one PF particle into a PuppiObj with global coordinates; an invalid particle
    // (empty link, or one whose pt was zeroed by the puppi weight) gives an empty object
    inline void toPuppi(const l1ct::PFRegionEmu &region, const PFParticleEmu &in, l1ct::PuppiObjEmu &out) {
      if (!in.valid()) {
        out.clear();
      } else if (in.kind == PFParticleEmu::Charged) {
        out.fill(region, in.charged);
      } else {
        out.fill(region, in.neutral, in.puppiPt, in.puppiWgt);
      }
    }
    inline l1ct::PuppiObjEmu toPuppi(const l1ct::PFRegionEmu &region, const PFParticleEmu &in) {
      l1ct::PuppiObjEmu ret;
      toPuppi(region, in, ret);
      return ret;
    }

    // convert a whole PF region: charged particles, muons, photons and neutral hadrons, in
    // this order, skipping the empty slots. Returns the number of particles appended.
    inline unsigned int toPFParticles(const l1ct::OutputRegion &pf,
                                      std::vector<PFParticleEmu> &out,
                                      l1ct::puppiWgt_t wgt = 1.0) {
      unsigned int n = 0;
      for (const auto &c : pf.pfcharged) {
        if (c.hwPt != 0) {
          out.emplace_back(c);
          n++;
        }
      }
      for (const auto &m : pf.pfmuon) {
        if (m.hwPt != 0) {
          out.emplace_back(m);
          n++;
        }
      }
      for (const auto &p : pf.pfphoton) {
        if (p.hwPt != 0) {
          out.emplace_back(p, wgt);
          n++;
        }
      }
      for (const auto &h : pf.pfneutral) {
        if (h.hwPt != 0) {
          out.emplace_back(h, wgt);
          n++;
        }
      }
      return n;
    }

    // same, converting the particles into PuppiObj at the same time
    inline unsigned int toPuppi(const l1ct::PFRegionEmu &region,
                                const l1ct::OutputRegion &pf,
                                std::vector<l1ct::PuppiObjEmu> &out,
                                l1ct::puppiWgt_t wgt = 1.0) {
      std::vector<PFParticleEmu> particles;
      unsigned int n = toPFParticles(pf, particles, wgt);
      for (const auto &p : particles)
        out.push_back(toPuppi(region, p));
      return n;
    }

  }  // namespace pfsorter
}  // namespace l1ct

#endif
