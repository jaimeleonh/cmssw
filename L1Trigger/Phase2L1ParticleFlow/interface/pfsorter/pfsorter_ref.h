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
    unsigned int nclocks_ = 162;  // also covers the default constructor
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

    // flush the tree and the conversion register
    void reset() {
      PFSorterEmulatorT<l1ct::PuppiObjEmu>::reset();
      clearConverted_();
    }

    // single clock cycle taking one PF particle per link, each with the region it was
    // reconstructed in (PF coordinates are local to the region): the particles are
    // converted into PuppiObj and then pushed into the fifos.
    //
    // The conversion costs one clock cycle, as it does in firmware (pf_to_puppi registers
    // its output), so what is pushed into the fifos in this clock cycle is what was
    // converted in the previous one, and newEvent is delayed by the same cycle.
    bool step(bool newEvent,
              const std::vector<l1ct::PFRegionEmu>& regions,
              const std::vector<PFParticle>& inputs,
              l1ct::PuppiObjEmu& out) {
      assert(regions.size() == inputs.size());
      std::vector<l1ct::PuppiObjEmu> converted(inputs.size());
      for (unsigned int i = 0, n = inputs.size(); i < n; ++i) {
        l1ct::pfsorter::toPuppi(regions[i], inputs[i], converted[i]);
      }
      if (converted_.size() != inputs.size()) {
        converted_.resize(inputs.size());
        clearConverted_();
      }
      bool ret = step(convertedNewEvent_, converted_, out);
      converted_.swap(converted);
      convertedNewEvent_ = newEvent;
      return ret;
    }

    // same, when all the links carry particles of the same region
    bool step(bool newEvent,
              const l1ct::PFRegionEmu& region,
              const std::vector<PFParticle>& inputs,
              l1ct::PuppiObjEmu& out) {
      return step(newEvent, std::vector<l1ct::PFRegionEmu>(inputs.size(), region), inputs, out);
    }

    // keep clocking with no new inputs until the tree is empty. The first of those clock
    // cycles is the one that lets the last converted particles into the fifos.
    unsigned int drain(std::vector<l1ct::PuppiObjEmu>& out, unsigned int maxclocks = 10000) {
      unsigned int nclocks = 0;
      if (!converted_.empty()) {
        l1ct::PuppiObjEmu obj;
        if (PFSorterEmulatorT<l1ct::PuppiObjEmu>::step(convertedNewEvent_, converted_, obj))
          out.push_back(obj);
        clearConverted_();
        nclocks++;
      }
      return nclocks + PFSorterEmulatorT<l1ct::PuppiObjEmu>::drain(out, maxclocks);
    }

    struct PFLinkItem {
      l1ct::pfsorter::PFParticleEmu particle;
      l1ct::PFRegionEmu region;
    };

    // Collect the puppi candidates of an event and distribute them over 'nlinks' input
    // links, one link per region, round robin if there are more regions than links.
    // Returns the number of objects. If 'regionIndex' is given, it is filled in parallel
    // to 'links' with the index of the region each object came from.
    static unsigned int makeLinks(const l1ct::Event &event,
                                  unsigned int nlinks,
                                  std::vector<std::vector<l1ct::PuppiObjEmu>> &links,
                                  std::vector<std::vector<unsigned int>> *regionIndex = nullptr) {
      links.clear();
      links.resize(nlinks);
      if (regionIndex) {
        regionIndex->clear();
        regionIndex->resize(nlinks);
      }
      unsigned int nobj = 0;
      for (unsigned int ireg = 0, nreg = event.out.size(); ireg < nreg; ++ireg) {
        for (const auto &p : event.out[ireg].puppi) {
          if (p.hwPt == 0)
            continue;
          links[ireg % nlinks].push_back(p);
          if (regionIndex)
            (*regionIndex)[ireg % nlinks].push_back(ireg);
          nobj++;
        }
      }
      return nobj;
    }

    // Same, taking the particles produced by the PF algorithm (charged, muons, photons and
    // neutral hadrons of each region) instead of the puppi candidates. They are handed to
    // the sorter as PF particles and converted into PuppiObj on the way in.
    static unsigned int makePFLinks(const std::vector<l1ct::PFInputRegion> &pfinputs,
                                const std::vector<l1ct::OutputRegion> &pfouts,
                                unsigned int nlinks,
                                std::vector<std::vector<PFLinkItem>> &links,
                                l1ct::puppiWgt_t wgt = 1.0,
                                std::vector<std::vector<unsigned int>> *regionIndex = nullptr) {
      links.clear();
      links.resize(nlinks);
      if (regionIndex) {
        regionIndex->clear();
        regionIndex->resize(nlinks);
      }
      unsigned int nobj = 0;
      // every output region must come with its input region, or the fiducial cut below has
      // no geometry to cut on and would silently drop the whole region
      assert(pfouts.size() == pfinputs.size());
      for (unsigned int ireg = 0, nreg = pfouts.size(); ireg < nreg; ++ireg) {
        const l1ct::PFRegionEmu &region = pfinputs[ireg].region;
        std::vector<l1ct::pfsorter::PFParticleEmu> particles;
        l1ct::pfsorter::toPFParticles(pfouts[ireg], particles, wgt);
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
          if (regionIndex)
            (*regionIndex)[ireg % nlinks].push_back(ireg);
          nobj++;
        }
      }
      return nobj;
    }

    // same, from a full event (the two vectors above are event.pfinputs and event.out)
    static unsigned int makePFLinks(const l1ct::Event &event,
                                unsigned int nlinks,
                                std::vector<std::vector<PFLinkItem>> &links,
                                l1ct::puppiWgt_t wgt = 1.0,
                                std::vector<std::vector<unsigned int>> *regionIndex = nullptr) {
      return makePFLinks(event.pfinputs, event.out, nlinks, links, wgt, regionIndex);
    }

    // Run the sorter over one event's worth of PF output, one link per region, and collect
    // the sorted puppi candidates. This is the whole of the algorithm; callers that do not
    // have an l1ct::Event at hand (a standalone testbench working region by region) can
    // pass the two vectors directly.
    //
    // The tree is sized on the number of regions if it has not been sized already, or if it
    // was sized for a different number of links, so that a default constructed emulator can
    // be used as is.
    void run(const std::vector<l1ct::PFInputRegion> &pfinputs,
             const std::vector<l1ct::OutputRegion> &pfouts,
             std::vector<l1ct::PuppiObjEmu> &out) {
      const unsigned int nlinks = pfinputs.size();
      if (buffer().nfifos() != nlinks)
        initFifos(nlinks, buffer().noutputs() ? buffer().noutputs() : 3);
      std::vector<std::vector<PFLinkItem>> links;
      makePFLinks(pfinputs, pfouts, nlinks, links);
      l1ct::PuppiObjEmu outObj;
      for (unsigned int iclock = 0; iclock < nClocks(); ++iclock) {
        bool newevt = (iclock == 0);
        std::vector<l1ct::pfsorter::PFParticleEmu> pfin(nlinks);
        std::vector<l1ct::PFRegionEmu> pfregions(nlinks);
        for (unsigned int l = 0; l < nlinks; ++l) {
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

#ifdef CMSSW_GIT_HASH
    // writes the sorted candidates into event.sortedpf, so 'event' can't be const here
    void run(l1ct::Event &event) {
      run(event.pfinputs, event.out, event.sortedpf);
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

  private:
    // the output register of the conversion: what it holds enters the fifos on the next
    // clock cycle, together with the newEvent that came with it
    std::vector<l1ct::PuppiObjEmu> converted_;
    bool convertedNewEvent_ = false;

    void clearConverted_() {
      for (auto& o : converted_)
        o.clear();
      convertedNewEvent_ = false;
    }
  };
}  // namespace l1ct

#endif
