#include <memory>
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/L1Scouting/interface/OrbitCollection.h"
#include "DataFormats/NanoAOD/interface/OrbitFlatTable.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingPuppi.h"
#include "DataFormats/L1TParticleFlow/interface/L1ScoutingTkEm.h"
#include "L1TriggerScouting/Utilities/interface/BxOffsetsFiller.h"

#include "DataFormats/Math/interface/deltaR.h"
#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2PuppiWDsGammaDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiWDsGammaDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiWDsGammaDemo() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T, typename U>
  void runObj(const OrbitCollection<T> &srcPuppi,
              const OrbitCollection<U> &srcTkEm,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  edm::EDGetTokenT<OrbitCollection<l1Scouting::Puppi>> structPuppiToken_;
  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEm>> structTkEmToken_;

  struct Cuts {
    float minpt1 = 3;
    float minpt2 = 5;
    float minpt3 = 10;
    float minpt4 = 25;
    float maxdeltar2 = 0.15 * 0.15;
    float maxdeltarDsGamma2 = 3.5 * 3.5;
    float mindeltaphiDsGamma = 2.5;
    float minmass = 60;
    float maxmass = 100;
    float minmass3 = 1.75;
    float maxmass3 = 2.30;
    float mindr2 = 0.00 * 0.00;
    float maxdr2 = 0.50 * 0.50;
    float maxiso = 0.45;
    float mindr2tkem = 0.02 * 0.02;
    float maxdr2tkem = 0.50 * 0.50;
    float maxisotkem = 0.25;
  } cuts;

  template <typename T>
  bool isolationDs(
      unsigned int pidex1, unsigned int pidex2, unsigned int pidex3, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  bool deltar(float eta1, float eta2, float phi1, float phi2) const;
  bool deltarmin(float eta1, float eta2, float phi1, float phi2) const;
  bool deltaphi(float phi1, float phi2) const;

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiWDsGammaDemo::ScPhase2PuppiWDsGammaDemo(const edm::ParameterSet &iConfig)
    : structPuppiToken_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("srcPuppi"))),
      structTkEmToken_(consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcTkEm"))) {
  produces<std::vector<unsigned>>("selectedBx");
  produces<l1ScoutingRun3::OrbitFlatTable>("wdsgamma");
  auto ptPis = iConfig.getParameter<std::vector<double>>("ptHad");
  if (ptPis.size() != 3 || ptPis[1] > ptPis[0] || ptPis[2] > ptPis[1]) {
    throw cms::Exception("InvalidConfiguration") << "ptPi must have exactly 3 elements, in descending order";
  }
  cuts.minpt1 = ptPis[2];
  cuts.minpt2 = ptPis[1];
  cuts.minpt3 = ptPis[0];
  cuts.minpt4 = iConfig.getParameter<double>("ptTkEm");
  cuts.maxdeltar2 = std::pow(iConfig.getParameter<double>("maxDeltaRHad"), 2);
  cuts.maxdeltarDsGamma2 = std::pow(iConfig.getParameter<double>("maxDeltaRDsTkEm"), 2);
  cuts.mindeltaphiDsGamma = iConfig.getParameter<double>("minDeltaPhiDsTkEm");
  cuts.minmass = iConfig.getParameter<double>("minMass");
  cuts.maxmass = iConfig.getParameter<double>("maxMass");
  cuts.minmass3 = iConfig.getParameter<double>("minMassDs");
  cuts.maxmass3 = iConfig.getParameter<double>("maxMassDs");
  cuts.maxiso = iConfig.getParameter<double>("relIsoDs");
  cuts.maxisotkem = iConfig.getParameter<double>("relIsoTkEm");
  cuts.mindr2 = std::pow(iConfig.getParameter<double>("isolationMinDeltaRDs"), 2);
  cuts.maxdr2 = std::pow(iConfig.getParameter<double>("isolationMaxDeltaRDs"), 2);
  cuts.mindr2tkem = std::pow(iConfig.getParameter<double>("isolationMinDeltaRTkEm"), 2);
  cuts.maxdr2tkem = std::pow(iConfig.getParameter<double>("isolationMaxDeltaRTkEm"), 2);
}

ScPhase2PuppiWDsGammaDemo::~ScPhase2PuppiWDsGammaDemo() {};

void ScPhase2PuppiWDsGammaDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiWDsGammaDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::Puppi>> srcPuppi;
  iEvent.getByToken(structPuppiToken_, srcPuppi);

  edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
  iEvent.getByToken(structTkEmToken_, srcTkEm);

  runObj(*srcPuppi, *srcTkEm, iEvent, countStruct_, passStruct_, "");
}

void ScPhase2PuppiWDsGammaDemo::endStream() {
  edm::LogImportant("ScPhase2AnalysisSummary")
      << "WDsGammma Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T, typename U>
void ScPhase2PuppiWDsGammaDemo::runObj(const OrbitCollection<T> &srcPuppi,
                                       const OrbitCollection<U> &srcTkEm,
                                       edm::Event &iEvent,
                                       unsigned long &nTry,
                                       unsigned long &nPass,
                                       const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s, i2s, i3s;  //i3s is the photon
  ROOT::RVec<unsigned int> ix;              // pions, kaons
  ROOT::RVec<unsigned int> ig;              // photons
  std::array<unsigned int, 4> bestQuadruplet{{0, 0, 0, 0}};
  float bestQuadrupletScore, bestQuadrupletMass;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = srcPuppi.bxIterator(bx);
    auto size = range.size();
    const T *cands = (size > 0) ? &range.front() : nullptr;

    auto range2 = srcTkEm.bxIterator(bx);
    auto size2 = range2.size();
    const U *cands2 = (size2 > 0) ? &range2.front() : nullptr;

    ix.clear();
    int intermediatecut = 0;
    int highcut = 0;
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minpt1) {
          ix.push_back(i);
          if (cands[i].pt() >= cuts.minpt2)
            intermediatecut++;
          if (cands[i].pt() >= cuts.minpt3)
            highcut++;
        }
      }
    }
    unsigned int npions = ix.size();
    if (highcut < 1 || intermediatecut < 2 || npions < 3)
      continue;

    ig.clear();
    for (unsigned int i = 0; i < size2; ++i) {  //make list of isolated photons
      if (cands2[i].pt() >= cuts.minpt4) {
        if (isolationTkEm(cands2[i].pt(), cands2[i].eta(), cands2[i].phi(), cands, size)) {
          ig.push_back(i);
        }
      }
    }

    unsigned int ngammas = ig.size();
    if (ngammas < 1)
      continue;

    bestQuadrupletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (cands[ix[i1]].pt() < cuts.minpt3)
        continue;  //high pt cut
      for (unsigned int i2 = 0; i2 < npions; ++i2) {
        if (i2 == i1 || cands[ix[i2]].pt() < cuts.minpt2)
          continue;
        if (cands[ix[i2]].pt() > cands[ix[i1]].pt() || (cands[ix[i2]].pt() == cands[ix[i1]].pt() and i2 < i1))
          continue;  //intermediate pt cut
        if (!deltar(cands[ix[i1]].eta(), cands[ix[i2]].eta(), cands[ix[i1]].phi(), cands[ix[i2]].phi()))
          continue;  //angular sep of top 2 tracks
        for (unsigned int i3 = 0; i3 < npions; ++i3) {
          if (i3 == i1 or i3 == i2)
            continue;
          if (cands[ix[i3]].pt() > cands[ix[i1]].pt() || (cands[ix[i3]].pt() == cands[ix[i1]].pt() and i3 < i1))
            continue;
          if (cands[ix[i3]].pt() > cands[ix[i2]].pt() || (cands[ix[i3]].pt() == cands[ix[i2]].pt() and i3 < i2))
            continue;

          if (!deltar(cands[ix[i1]].eta(), cands[ix[i3]].eta(), cands[ix[i1]].phi(), cands[ix[i3]].phi()) ||
              !deltar(cands[ix[i2]].eta(), cands[ix[i3]].eta(), cands[ix[i2]].phi(), cands[ix[i3]].phi()))
            continue;  //angular sep of 3rd track from the other two

          if (std::abs(cands[ix[i1]].charge() + cands[ix[i2]].charge() + cands[ix[i3]].charge()) != 1)
            continue;  //Ds charge requirement

          auto dsP4 = cands[ix[i1]].p4() + cands[ix[i2]].p4() + cands[ix[i3]].p4();
          auto mass3 = dsP4.mass();
          if (!(cuts.minmass3 <= mass3 && mass3 <= cuts.maxmass3))
            continue;  // Ds mass cut

          if (!isolationDs(ix[i1], ix[i2], ix[i3], cands, size))
            continue;

          for (unsigned int i4 = 0; i4 < ngammas; ++i4) {
            auto mass = (dsP4 + cands2[ig[i4]].p4()).mass();
            if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
              bool pass_deltaphi = deltaphi(dsP4.phi(), cands2[ig[i4]].phi());
              bool pass_deltar = deltarmin(dsP4.eta(), cands2[ig[i4]].eta(), dsP4.phi(), cands2[ig[i4]].phi());
              if (pass_deltaphi && pass_deltar) {
                float ptsum = cands[ix[i1]].pt() + cands[ix[i2]].pt() + cands[ix[i3]].pt() + cands2[ig[i4]].pt();
                if (ptsum > bestQuadrupletScore) {
                  bestQuadruplet[0] = ix[i1];
                  bestQuadruplet[1] = ix[i2];
                  bestQuadruplet[2] = ix[i3];
                  bestQuadruplet[3] = ig[i4];
                  bestQuadrupletScore = ptsum;
                  bestQuadrupletMass = mass;
                }  // best
              }  // delta R
            }  // mass
          }  // photon loop
        }  //low pt cut
      }  //intermediate pt cut
    }  //high pt cut

    if (bestQuadrupletScore > 0) {
      ret->emplace_back(bx);
      nPass++;
      masses.push_back(bestQuadrupletMass);
      i0s.push_back(bestQuadruplet[0]);
      i1s.push_back(bestQuadruplet[1]);
      i2s.push_back(bestQuadruplet[2]);
      i3s.push_back(bestQuadruplet[3]);
      bxOffsetsFiller.addBx(bx, 1);
    }
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "WDsGamma" + label, true);
  tab->addColumn<float>("mass", masses, "W invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "leading pion");
  tab->addColumn<uint8_t>("i1", i1s, "subleading pion");
  tab->addColumn<uint8_t>("i2", i2s, "trailing pion");
  tab->addColumn<uint8_t>("i3", i3s, "photon");
  iEvent.put(std::move(tab), "wdsgamma" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiWDsGammaDemo::isolationDs(
    unsigned int pidex1, unsigned int pidex2, unsigned int pidex3, const T *cands, unsigned int size) const {
  float psum = 0;
  float eta = cands[pidex1].eta();  //center cone around leading track
  float phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j or pidex2 == j or pidex3 == j)
      continue;
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  return (psum <= cuts.maxiso * (cands[pidex1].pt() + cands[pidex2].pt() + cands[pidex3].pt()));
}

template <typename T>
bool ScPhase2PuppiWDsGammaDemo::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.mindr2tkem && dr2 <= cuts.maxdr2tkem)
      psum += cands[j].pt();
  }
  return (psum <= cuts.maxisotkem * pt);
}

bool ScPhase2PuppiWDsGammaDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  float dr2 = reco::deltaR2(eta1, phi1, eta2, phi2);
  return (dr2 < cuts.maxdeltar2);
}

bool ScPhase2PuppiWDsGammaDemo::deltarmin(float eta1, float eta2, float phi1, float phi2) const {
  float dr2 = reco::deltaR2(eta1, phi1, eta2, phi2);
  return (dr2 < cuts.maxdeltarDsGamma2);
}

bool ScPhase2PuppiWDsGammaDemo::deltaphi(float phi1, float phi2) const {
  return std::abs(reco::deltaPhi(phi1, phi2)) >= cuts.mindeltaphiDsGamma;
}

void ScPhase2PuppiWDsGammaDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPuppi");
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<std::vector<double>>("ptHad", {10., 5., 3.});
  desc.add<double>("ptTkEm", 20.);
  desc.add<double>("minMassDs", 1.75);
  desc.add<double>("maxMassDs", 2.30);
  desc.add<double>("minMass", 60.);
  desc.add<double>("maxMass", 100.);
  desc.add<double>("maxDeltaRHad", 0.15);
  desc.add<double>("maxDeltaRDsTkEm", 3.5);
  desc.add<double>("minDeltaPhiDsTkEm", 2.5);
  desc.add<double>("relIsoDs", 0.45);
  desc.add<double>("relIsoTkEm", 0.25);
  desc.add<double>("isolationMinDeltaRDs", 0.00);
  desc.add<double>("isolationMaxDeltaRDs", 0.50);
  desc.add<double>("isolationMinDeltaRTkEm", 0.02);
  desc.add<double>("isolationMaxDeltaRTkEm", 0.50);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiWDsGammaDemo);
