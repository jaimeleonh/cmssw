#include <cstdint>
#include <vector>

#include "DataFormats/L1ScoutingSoA/interface/AssociationMapHost.h"
#include "DataFormats/L1ScoutingSoA/interface/BxLookupHost.h"
#include "DataFormats/L1ScoutingSoA/interface/PFCandidateHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/VertexHostCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/AssociationMapDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/BxLookupDevice.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/PFCandidateDeviceCollection.h"
#include "DataFormats/L1ScoutingSoA/interface/alpaka/VertexDeviceCollection.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EDPutToken.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/Event.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/EventSetup.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/MakerMacros.h"
#include "HeterogeneousCore/AlpakaCore/interface/alpaka/stream/FixedQueueEDProducer.h"
#include "HeterogeneousCore/AlpakaInterface/interface/config.h"

namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc {

  // Produce the device collections consumed by SoftJetIdONNX from explicit
  // per-particle parameters, for testing. The jets are given as ranges of a
  // flat PF candidate list (the association map index is the identity); each
  // BX entry of jetBxLookup/vertexBxLookup selects jets/vertices per BX.
  // The same data is uploaded for every event.
  class SoftJetONNXTestSource : public stream::FixedQueueEDProducer<> {
  public:
    SoftJetONNXTestSource(const edm::ParameterSet& params)
        : FixedQueueEDProducer<>(params),
          pf_token_{produces()},
          association_map_token_{produces()},
          jet_bx_lookup_token_{produces("jetBxLookup")},
          vertex_token_{produces()},
          vertex_bx_lookup_token_{produces("vertexBxLookup")},
          pt_{params.getParameter<std::vector<double>>("pt")},
          eta_{params.getParameter<std::vector<double>>("eta")},
          phi_{params.getParameter<std::vector<double>>("phi")},
          z0_{params.getParameter<std::vector<double>>("z0")},
          dxy_{params.getParameter<std::vector<double>>("dxy")},
          puppiw_{params.getParameter<std::vector<double>>("puppiw")},
          pdgid_{params.getParameter<std::vector<int>>("pdgid")},
          jet_offsets_{params.getParameter<std::vector<unsigned int>>("jetOffsets")},
          jet_bx_offsets_{params.getParameter<std::vector<unsigned int>>("jetBxOffsets")},
          vertex_z0_{params.getParameter<std::vector<double>>("vertexZ0")},
          vertex_bx_offsets_{params.getParameter<std::vector<unsigned int>>("vertexBxOffsets")} {
      check_sizes();
    }

    static void fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
      edm::ParameterSetDescription desc;
      // PF candidates, concatenated over all jets
      desc.add<std::vector<double>>("pt");
      desc.add<std::vector<double>>("eta");
      desc.add<std::vector<double>>("phi");
      desc.add<std::vector<double>>("z0");
      desc.add<std::vector<double>>("dxy");
      desc.add<std::vector<double>>("puppiw");
      desc.add<std::vector<int>>("pdgid");
      // jets: constituent ranges of the flat candidate list, n_jets + 1 entries
      desc.add<std::vector<unsigned int>>("jetOffsets");
      // jets per BX: begin offset of each BX in the jet list, n_bx + 1 entries
      desc.add<std::vector<unsigned int>>("jetBxOffsets");
      // reconstructed vertices
      desc.add<std::vector<double>>("vertexZ0");
      // vertices per BX: begin offset of each BX in the vertex list, n_bx + 1 entries
      desc.add<std::vector<unsigned int>>("vertexBxOffsets");
      descriptions.addWithDefaultLabel(desc);
    }

    void produce(device::Event& event, const device::EventSetup& event_setup) override {
      const uint32_t n_particles = pt_.size();
      const uint32_t n_jets = jet_offsets_.size() - 1;
      const uint32_t n_vertices = vertex_z0_.size();
      const uint32_t n_bx = jet_bx_offsets_.size() - 1;

      // PF candidates
      PFCandidateHostCollection host_pf(n_particles);
      auto pf_view = host_pf.view();
      for (uint32_t i = 0; i < n_particles; ++i) {
        pf_view.pt()[i] = pt_[i];
        pf_view.eta()[i] = eta_[i];
        pf_view.phi()[i] = phi_[i];
        pf_view.z0()[i] = z0_[i];
        pf_view.dxy()[i] = dxy_[i];
        pf_view.puppiw()[i] = puppiw_[i];
        pf_view.quality()[i] = 0;
        pf_view.pdgid()[i] = static_cast<int16_t>(pdgid_[i]);
      }
      PFCandidateDeviceCollection pf(event.queue(), n_particles);
      alpaka::memcpy(event.queue(), pf.buffer(), host_pf.buffer());

      // jet -> constituents association map (identity index)
      AssociationMapHost host_map(n_particles, n_jets + 1);
      auto map_view = host_map.view();
      for (uint32_t i = 0; i < n_particles; ++i)
        map_view.index().index()[i] = i;
      for (uint32_t j = 0; j <= n_jets; ++j)
        map_view.offset().offset()[j] = jet_offsets_[j];
      AssociationMapDevice association_map(event.queue(), n_particles, n_jets + 1);
      alpaka::memcpy(event.queue(), association_map.buffer(), host_map.buffer());

      // jets per BX
      BxLookupHost host_jet_bx(n_bx, n_bx + 1);
      auto jet_bx_view = host_jet_bx.view();
      for (uint32_t bx = 0; bx < n_bx; ++bx)
        jet_bx_view.bx().bx()[bx] = static_cast<uint16_t>(bx);
      for (uint32_t bx = 0; bx <= n_bx; ++bx)
        jet_bx_view.offset().offset()[bx] = jet_bx_offsets_[bx];
      BxLookupDevice jet_bx_lookup(event.queue(), n_bx, n_bx + 1);
      alpaka::memcpy(event.queue(), jet_bx_lookup.buffer(), host_jet_bx.buffer());

      // vertices
      VertexHostCollection host_vertices(n_vertices);
      auto vertex_view = host_vertices.view();
      for (uint32_t v = 0; v < n_vertices; ++v) {
        vertex_view.valid()[v] = true;
        vertex_view.z0()[v] = vertex_z0_[v];
        vertex_view.multIn()[v] = 0;
        vertex_view.sumPt()[v] = 0.f;
        vertex_view.quality()[v] = 0;
        vertex_view.multOut()[v] = 0;
        vertex_view.unassigned()[v] = 0;
      }
      VertexDeviceCollection vertices(event.queue(), n_vertices);
      alpaka::memcpy(event.queue(), vertices.buffer(), host_vertices.buffer());

      // vertices per BX
      BxLookupHost host_vertex_bx(n_bx, n_bx + 1);
      auto vertex_bx_view = host_vertex_bx.view();
      for (uint32_t bx = 0; bx < n_bx; ++bx)
        vertex_bx_view.bx().bx()[bx] = static_cast<uint16_t>(bx);
      for (uint32_t bx = 0; bx <= n_bx; ++bx)
        vertex_bx_view.offset().offset()[bx] = vertex_bx_offsets_[bx];
      BxLookupDevice vertex_bx_lookup(event.queue(), n_bx, n_bx + 1);
      alpaka::memcpy(event.queue(), vertex_bx_lookup.buffer(), host_vertex_bx.buffer());

      // put device-side products into the event
      event.emplace(pf_token_, std::move(pf));
      event.emplace(association_map_token_, std::move(association_map));
      event.emplace(jet_bx_lookup_token_, std::move(jet_bx_lookup));
      event.emplace(vertex_token_, std::move(vertices));
      event.emplace(vertex_bx_lookup_token_, std::move(vertex_bx_lookup));
    }

  private:
    void check_sizes() const {
      const size_t n_particles = pt_.size();
      const std::vector<std::pair<const char*, size_t>> fields{
          {"eta", eta_.size()},        {"phi", phi_.size()}, {"z0", z0_.size()},
          {"dxy", dxy_.size()},        {"puppiw", puppiw_.size()},
          {"pdgid", pdgid_.size()},
      };
      for (const auto& [name, size] : fields) {
        if (size != n_particles)
          throw cms::Exception("Configuration")
              << "SoftJetONNXTestSource: parameter \"" << name << "\" has " << size
              << " entries, expected " << n_particles << " (the size of \"pt\").";
      }
      if (jet_offsets_.size() < 2 || jet_offsets_.front() != 0 ||
          jet_offsets_.back() != n_particles)
        throw cms::Exception("Configuration")
            << "SoftJetONNXTestSource: \"jetOffsets\" must start at 0 and end at the number of PF candidates ("
            << n_particles << ").";
      if (jet_bx_offsets_.size() < 2 || jet_bx_offsets_.front() != 0 ||
          jet_bx_offsets_.back() != jet_offsets_.size() - 1)
        throw cms::Exception("Configuration")
            << "SoftJetONNXTestSource: \"jetBxOffsets\" must start at 0 and end at the number of jets ("
            << jet_offsets_.size() - 1 << ").";
      if (vertex_bx_offsets_.size() != jet_bx_offsets_.size())
        throw cms::Exception("Configuration")
            << "SoftJetONNXTestSource: \"vertexBxOffsets\" and \"jetBxOffsets\" must have the same number of BX "
               "entries.";
      if (vertex_bx_offsets_.front() != 0 || vertex_bx_offsets_.back() != vertex_z0_.size())
        throw cms::Exception("Configuration")
            << "SoftJetONNXTestSource: \"vertexBxOffsets\" must start at 0 and end at the number of vertices ("
            << vertex_z0_.size() << ").";
    }

    const device::EDPutToken<PFCandidateDeviceCollection> pf_token_;
    const device::EDPutToken<AssociationMapDevice> association_map_token_;
    const device::EDPutToken<BxLookupDevice> jet_bx_lookup_token_;
    const device::EDPutToken<VertexDeviceCollection> vertex_token_;
    const device::EDPutToken<BxLookupDevice> vertex_bx_lookup_token_;

    const std::vector<double> pt_, eta_, phi_, z0_, dxy_, puppiw_;
    const std::vector<int> pdgid_;
    const std::vector<unsigned int> jet_offsets_, jet_bx_offsets_;
    const std::vector<double> vertex_z0_;
    const std::vector<unsigned int> vertex_bx_offsets_;
  };

}  // namespace ALPAKA_ACCELERATOR_NAMESPACE::l1sc

DEFINE_FWK_ALPAKA_MODULE(l1sc::SoftJetONNXTestSource);
