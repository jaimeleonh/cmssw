#ifndef pfsorter_ref_h
#define pfsorter_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_tree_elements_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_inputs_ref.h"

#ifdef CMSSW_GIT_HASH
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#endif

namespace l1ct {

  // Clock-cycle emulation of the pfsorter merge tree, built on the same fifo / staging
  // area / queue structure as the multififo regionizer (see pfsorter_tree_elements_ref.h).
  //
  //   link 0 --> fifo 0 --.
  //                        >-- o --.
  //   link 1 --> fifo 1 --'         \
  //   link 2 --> fifo 2 --.          >-- o --> one object per clock cycle
  //                        >-- o --'
  //   link 3 --> fifo 3 --'
  //
  // One buffer per link (1:1:1 link / buffer / fifo of the first layer), any number of
  // links, at most 3 output nodes, one object per node.
  //
  // Usage:
  //    l1ct::PFSorterEmulator sorter(/*nlinks=*/12, /*noutputs=*/3);
  //    l1ct::PuppiObjEmu out;
  //    for (int iclock = 0; iclock < nclocks; ++iclock) {
  //      bool newEvent = (iclock % nclocksPerEvent == 0);
  //      sorter.step(newEvent, links[iclock], out);  // out is empty if nothing came out
  //    }
  template <typename T>
  class PFSorterEmulatorT {
  public:
    PFSorterEmulatorT() : nevt_(0) {}
    PFSorterEmulatorT(unsigned int nlinks, unsigned int noutputs = 3, unsigned int nclocks = 162) : nevt_(0), nclocks_(nclocks) {
      buffer_.initFifos(nlinks, noutputs);
    }

    void initFifos(unsigned int nlinks, unsigned int noutputs = 3) {
      buffer_.initFifos(nlinks, noutputs);
    }

    void reset() {
      buffer_.reset();
      nevt_ = 0;
    }

    // single clock emulation; on a new event the tree is emptied first.
    // 'inputs' must have one (possibly empty, i.e. hwPt == 0) object per link, and they go
    // into the fifo with the same index. 'out' is the object leaving the tree in this clock
    // cycle; the return value says whether it is a valid one.
    bool step(bool newEvent, const std::vector<T>& inputs, T& out) {
      if (newEvent) {
        buffer_.flush();
        nevt_++;
      }
      assert(inputs.size() == buffer_.nfifos());
      for (unsigned int i = 0, n = inputs.size(); i < n; ++i) {
        buffer_.maybe_push(i, inputs[i]);
      }
      out = buffer_.pop();
      return out.hwPt != 0;
    }

    // keep clocking with no new inputs until the tree is empty, appending to 'out'
    // everything that comes out. Returns the number of clock cycles used.
    unsigned int drain(std::vector<T>& out, unsigned int maxclocks = 10000) {
      unsigned int nclocks = 0;
      while (buffer_.inFlight() != 0 && nclocks < maxclocks) {
        T obj = buffer_.pop();
        nclocks++;
        if (obj.hwPt != 0)
          out.push_back(obj);
      }
      return nclocks;
    }

    const l1ct::pfsorter::RegionBuffer<T>& buffer() const { return buffer_; }
    l1ct::pfsorter::RegionBuffer<T>& buffer() { return buffer_; }
    unsigned int nEvents() const { return nevt_; }
    unsigned int nClocks() const { return nclocks_; }

  private:
    l1ct::pfsorter::RegionBuffer<T> buffer_;
    unsigned int nevt_;
    unsigned int nclocks_;
  };

  // The sorter as it is used in practice: it stores PuppiObj, and it can be fed either
  // directly with PuppiObj or with the particles coming out of the ParticleFlow algorithm,
  // which are converted into PuppiObj on the way in (see pfsorter_inputs_ref.h).
  class PFSorterEmulator : public PFSorterEmulatorT<l1ct::PuppiObjEmu> {
  public:
    typedef l1ct::pfsorter::PFParticleEmu PFParticle;
    using PFSorterEmulatorT<l1ct::PuppiObjEmu>::PFSorterEmulatorT;
    using PFSorterEmulatorT<l1ct::PuppiObjEmu>::step;

#ifdef CMSSW_GIT_HASH
    PFSorterEmulator(const edm::ParameterSet& iConfig)
        : PFSorterEmulator(iConfig.getParameter<uint32_t>("nLinks"),
                           iConfig.getParameter<uint32_t>("nOutputs"),
                           iConfig.getParameter<uint32_t>("nClocks")) {}

    static edm::ParameterSetDescription getParameterSetDescription() {
      edm::ParameterSetDescription description;
      description.add<uint32_t>("nLinks", 18);
      description.add<uint32_t>("nOutputs", 3);
      description.add<uint32_t>("nClocks", 162);
      return description;
    }
#endif

    // single clock cycle taking one PF particle per link, each with the region it was
    // reconstructed in (PF coordinates are local to the region): the particles are
    // converted into PuppiObj and then pushed into the fifos.
    bool step(bool newEvent,
              const std::vector<l1ct::PFRegionEmu>& regions,
              const std::vector<PFParticle>& inputs,
              l1ct::PuppiObjEmu& out) {
      assert(regions.size() == inputs.size());
      std::vector<l1ct::PuppiObjEmu> converted(inputs.size());
      for (unsigned int i = 0, n = inputs.size(); i < n; ++i) {
        l1ct::pfsorter::toPuppi(regions[i], inputs[i], converted[i]);
      }
      return step(newEvent, converted, out);
    }

    // same, when all the links carry particles of the same region
    bool step(bool newEvent,
              const l1ct::PFRegionEmu& region,
              const std::vector<PFParticle>& inputs,
              l1ct::PuppiObjEmu& out) {
      return step(newEvent, std::vector<l1ct::PFRegionEmu>(inputs.size(), region), inputs, out);
    }

    struct PFLinkItem {
      l1ct::pfsorter::PFParticleEmu particle;
      l1ct::PFRegionEmu region;
    };

#ifdef CMSSW_GIT_HASH
    void makePFLinks(const l1ct::Event &event,
                                unsigned int nlinks,
                                std::vector<std::vector<PFLinkItem>> &links,
                                l1ct::puppiWgt_t wgt = 1.0) {
      links.clear();
      links.resize(nlinks);
      // every output region must come with its input region, or the fiducial cut below has
      // no geometry to cut on and would silently drop the whole region
      assert(event.out.size() == event.pfinputs.size());
      for (unsigned int ireg = 0, nreg = event.out.size(); ireg < nreg; ++ireg) {
        const l1ct::PFRegionEmu &region = event.pfinputs[ireg].region;
        std::vector<l1ct::pfsorter::PFParticleEmu> particles;
        l1ct::pfsorter::toPFParticles(event.out[ireg], particles, wgt);
        for (auto &p : particles) {
          // Skip an empty slot, a neutral whose pt was zeroed by the puppi weight, and
          // anything sitting in the overlap pad of the region: regions are padded by
          // etaExtra/phiExtra, so a particle in the pad is reconstructed a second time by
          // the neighbouring region that holds it as fiducial, and keeping both copies
          // would duplicate it in the sorter output. This is the same cut that fetchPF()
          // and linpuppi_ref() apply to these very same PF collections.
          if (!p.valid() || !region.isFiducial(p))
            continue;
          links[ireg % nlinks].push_back(PFLinkItem{p, region});
        }
      }
      return;
    }

    void run(const l1ct::Event &event,
            std::vector<l1ct::PuppiObjEmu>& out) {
      std::vector<std::vector<PFLinkItem>> links;
      makePFLinks(event, event.pfinputs.size(), links);
      l1ct::PuppiObjEmu outObj;
      for (unsigned int iclock = 0; iclock < nClocks(); ++iclock) {
        bool newevt = (iclock == 0);
        std::vector<l1ct::pfsorter::PFParticleEmu> pfin( event.pfinputs.size());
        std::vector<l1ct::PFRegionEmu> pfregions(event.pfinputs.size());
        for (unsigned int l = 0; l < event.pfinputs.size(); ++l) {
          if (iclock < links[l].size()) {
            pfin[l] = links[l][iclock].particle;
            pfregions[l] = links[l][iclock].region;
          }
        }
        step(newevt, pfregions, pfin, outObj);
        if (outObj.hwPt != 0)
          out.push_back(outObj);
      }
    }
#else
    // standalone version: run from raw links (one vector of particles per link)
    void run(const std::vector<l1ct::PFRegionEmu>& regions,
             const std::vector<std::vector<PFParticle>>& links,
             std::vector<l1ct::PuppiObjEmu>& out) {
      l1ct::PuppiObjEmu outObj;
      for (unsigned int iclock = 0; iclock < nClocks(); ++iclock) {
        bool newevt = (iclock == 0);
        std::vector<l1ct::pfsorter::PFParticleEmu> pfin(regions.size());
        std::vector<l1ct::PFRegionEmu> pfregions(regions.size());
        for (unsigned int l = 0; l < regions.size(); ++l) {
          if (iclock < links[l].size()) {
            pfin[l] = links[l][iclock];
            pfregions[l] = regions[l];
          }
        }
        step(newevt, pfregions, pfin, outObj);
        if (outObj.hwPt != 0)
          out.push_back(outObj);
      }
    }
#endif
  };
}  // namespace l1ct

#endif
