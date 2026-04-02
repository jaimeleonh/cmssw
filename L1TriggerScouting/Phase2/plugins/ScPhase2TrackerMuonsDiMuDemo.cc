#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/PFCandidate.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1TMuonPhase2/interface/L1ScoutingTrackerMuon.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include <DataFormats/Math/interface/deltaR.h>
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2TrackerMuonDiMuDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TrackerMuonDiMuDemo(const edm::ParameterSet &);
  ~ScPhase2TrackerMuonDiMuDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &src,
              const OrbitCollection<l1Scouting::Puppi> &srcPuppi,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  edm::EDGetTokenT<OrbitCollection<l1Scouting::TrackerMuon>> structToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> puppiToken_;

  float isolation(
      float eta, float phi, float z0, const l1Scouting::Puppi *cands, unsigned int size, float &cache) const {
    if (cache == -1)
      cache = isolation(eta, phi, z0, cands, size);
    return cache;
  }

  float isolation(float eta, float phi, float z0, const l1Scouting::Puppi *cands, unsigned int size) const;

  struct Cuts {
    float minptOverMass = 0.25;
    unsigned qualityMin = 8;
    float massMin = 0.5;
    bool doOppositeCharge = true;
    float etaMax = 2.0;
    std::array<float, 2> ptMin{{2.0, 2.0}};
    float dzMax = 1.0;
    float minRelIso = 0.2;
    float isoMaxDeltaR2 = 0.02 * 0.02;
    float isoMinDeltaR2 = 0.40 * 0.40;
  } cuts;

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2TrackerMuonDiMuDemo::ScPhase2TrackerMuonDiMuDemo(const edm::ParameterSet &iConfig)
    : structToken_(consumes<OrbitCollection<l1Scouting::TrackerMuon>>(iConfig.getParameter<edm::InputTag>("src"))),
      puppiToken_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("srcPuppi"))) {
  produces<std::vector<unsigned>>("selectedBx");
  produces<l1ScoutingRun3::OrbitFlatTable>("dimu");
  auto ptMin = iConfig.getParameter<std::vector<double>>("ptMin");
  if (ptMin.size() != 2) {
    throw cms::Exception("Configuration")
        << "ptMin should be a vector of size 2, for the leading and subleading electron";
  }
  std::copy_n(ptMin.begin(), 2, cuts.ptMin.begin());
  cuts.etaMax = iConfig.getParameter<double>("etaMax");
  cuts.qualityMin = iConfig.getParameter<unsigned>("quality");
  cuts.massMin = iConfig.getParameter<double>("massMin");
  cuts.doOppositeCharge = iConfig.getParameter<bool>("oppositeCharge");
  cuts.dzMax = iConfig.getParameter<double>("dzMax");
  cuts.minptOverMass = iConfig.getParameter<double>("minptOverMass");
  cuts.minRelIso = iConfig.getParameter<double>("relIsoMax");
  cuts.isoMinDeltaR2 = std::pow(iConfig.getParameter<double>("isolationMinDeltaR"), 2);
  cuts.isoMaxDeltaR2 = std::pow(iConfig.getParameter<double>("isolationMaxDeltaR"), 2);
}

ScPhase2TrackerMuonDiMuDemo::~ScPhase2TrackerMuonDiMuDemo() {};

void ScPhase2TrackerMuonDiMuDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2TrackerMuonDiMuDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::TrackerMuon>> src;
  edm::Handle<OrbitCollection<l1Scouting::Puppi>> srcPuppi;
  iEvent.getByToken(structToken_, src);
  iEvent.getByToken(puppiToken_, srcPuppi);
  runObj(*src, *srcPuppi, iEvent, countStruct_, passStruct_, "");
}

void ScPhase2TrackerMuonDiMuDemo::endStream() {
  edm::LogImportant("ScPhase2AnalysisSummary")
      << "DiTrackerMuon Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2TrackerMuonDiMuDemo::runObj(const OrbitCollection<T> &src,
                                         const OrbitCollection<l1Scouting::Puppi> &srcPuppi,
                                         edm::Event &iEvent,
                                         unsigned long &nTry,
                                         unsigned long &nPass,
                                         const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto selectedBx_idx = std::make_unique<std::vector<unsigned>>();
  float mass;
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s;
  std::vector<float> iso0s, iso1s;
  std::vector<float> isoCache;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = src.bxIterator(bx);
    auto size = range.size();
    const T *cands = (size > 0) ? &range.front() : nullptr;

    auto puppiRange = srcPuppi.bxIterator(bx);
    auto puppiSize = puppiRange.size();
    const l1Scouting::Puppi *puppiCands = (puppiSize > 0) ? &puppiRange.front() : nullptr;

    unsigned int selPairs = 0;

    isoCache.resize(size);
    std::fill(isoCache.begin(), isoCache.end(), -1);

    for (unsigned int i = 0; i < size; ++i) {
      if (cands[i].pt() < cuts.ptMin[0])
        break;  // assumes pt ordering
      if (std::abs(cands[i].eta()) > cuts.etaMax)
        continue;
      if (cands[i].quality() < cuts.qualityMin)
        continue;

      if (isolation(cands[i].eta(), cands[i].phi(), cands[i].z0(), puppiCands, puppiSize, isoCache[i]) >
          cuts.minRelIso * cands[i].pt())
        continue;

      for (unsigned int j = i + 1; j < size; ++j) {
        if (cands[j].pt() < cuts.ptMin[1])
          break;  // assumes pt ordering
        if (cuts.doOppositeCharge and (cands[i].charge() * cands[j].charge() > 0))
          continue;
        if (std::abs(cands[i].z0() - cands[j].z0()) > cuts.dzMax)
          continue;
        if (std::abs(cands[j].eta()) > cuts.etaMax)
          continue;
        if (cands[j].quality() < cuts.qualityMin)
          continue;

        if (isolation(cands[j].eta(), cands[j].phi(), cands[j].z0(), puppiCands, puppiSize, isoCache[j]) >
            cuts.minRelIso * cands[j].pt())
          continue;

        mass = (cands[i].p4() + cands[j].p4()).mass();
        if (mass < cuts.massMin)
          continue;
        if (cands[j].pt() < cuts.minptOverMass * mass)
          continue;
        if (cands[i].pt() < cuts.minptOverMass * mass)
          continue;
        if (selPairs == 0)
          selectedBx_idx->emplace_back(bx);
        masses.push_back(mass);
        i0s.push_back(i);
        i1s.push_back(j);
        iso0s.push_back(isoCache[i]);
        iso1s.push_back(isoCache[j]);
        selPairs++;
      }
    }

    if (selPairs > 0)
      nPass++;

    bxOffsetsFiller.addBx(bx, selPairs);
  }  // loop on BXs

  iEvent.put(std::move(selectedBx_idx), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "DiMu" + label, false);
  tab->addColumn<float>("mass", masses, "Dimuon invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading muon");
  tab->addColumn<uint8_t>("i1", i1s, "subleading muon");
  tab->addColumn<float>("iso0", iso0s, "isolation of leading muon");
  tab->addColumn<float>("iso1", iso1s, "isolation of subleading muon");
  iEvent.put(std::move(tab), "dimu" + label);
}

float ScPhase2TrackerMuonDiMuDemo::isolation(
    float eta, float phi, float z0, const l1Scouting::Puppi *cands, unsigned int size) const {
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (std::abs(cands[j].pdgId()) == 13)
      continue;
    if (cands[j].charge() != 0 && std::abs(cands[j].z0() - z0) > cuts.dzMax)
      continue;
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.isoMinDeltaR2 && dr2 <= cuts.isoMaxDeltaR2)
      psum += cands[j].pt();
  }
  return psum;
}

void ScPhase2TrackerMuonDiMuDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<std::vector<double>>("ptMin", {2.0, 2.0});
  desc.add<double>("etaMax", 2.4);
  desc.add<unsigned>("quality", 8);
  desc.add<double>("massMin", 0.5);
  desc.add<bool>("oppositeCharge", true);
  desc.add<double>("dzMax", 0.5);
  desc.add<double>("minptOverMass", 0.25);
  desc.add<edm::InputTag>("srcPuppi");
  desc.add<double>("relIsoMax", 0.4);
  desc.add<double>("isolationMinDeltaR", 0.02);
  desc.add<double>("isolationMaxDeltaR", 0.4);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TrackerMuonDiMuDemo);
