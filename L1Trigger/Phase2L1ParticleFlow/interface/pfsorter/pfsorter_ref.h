#ifndef pfsorter_ref_h
#define pfsorter_ref_h

#include "DataFormats/L1TParticleFlow/interface/layer1_emulator.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_tree_elements_ref.h"
#include "L1Trigger/Phase2L1ParticleFlow/interface/pfsorter/pfsorter_inputs_ref.h"

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
    PFSorterEmulatorT(unsigned int nlinks, unsigned int noutputs = 3) : nevt_(0) {
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

  private:
    l1ct::pfsorter::RegionBuffer<T> buffer_;
    unsigned int nevt_;
  };

  // The sorter as it is used in practice: it stores PuppiObj, and it can be fed either
  // directly with PuppiObj or with the particles coming out of the ParticleFlow algorithm,
  // which are converted into PuppiObj on the way in (see pfsorter_inputs_ref.h).
  class PFSorterEmulator : public PFSorterEmulatorT<l1ct::PuppiObjEmu> {
  public:
    typedef l1ct::pfsorter::PFParticleEmu PFParticle;
    using PFSorterEmulatorT<l1ct::PuppiObjEmu>::PFSorterEmulatorT;
    using PFSorterEmulatorT<l1ct::PuppiObjEmu>::step;

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
  };

}  // namespace l1ct

#endif
