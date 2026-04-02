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

class ScPhase2PuppiWPiGammaDemo : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2PuppiWPiGammaDemo(const edm::ParameterSet &);
  ~ScPhase2PuppiWPiGammaDemo() override;
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
    float minpt_pi = 25;
    float minpt_tkem = 20;
    float minmass = 60;
    float maxmass = 100;
    float maxiso_pi = 0.30;
    float maxiso_tkem = 0.30;
    float mindeltar2 = 0.50 * 0.50;
    float mindr2 = 0.00 * 0.00;
    float maxdr2 = 0.50 * 0.50;
    float mindr2tkem = 0.02 * 0.02;
    float maxdr2tkem = 0.50 * 0.50;
  } cuts;

  template <typename T>
  bool isolationPi(unsigned int pidex1, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationPi(unsigned int pidex1, const T *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolationPi(pidex1, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const;

  template <typename T>
  bool isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size, unsigned int &cache) const {
    if (cache == 0)
      cache = isolationTkEm(pt, eta, phi, cands, size) ? 1 : 2;
    return (cache == 1);
  }

  bool deltar(float eta1, float eta2, float phi1, float phi2) const;
  bool deltarmin(float eta1, float eta2, float phi1, float phi2) const;
  bool deltaphi(float phi1, float phi2) const;
  static float doubletmass(const std::array<unsigned int, 2> &t, const float *pts, const float *etas, const float *phis);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2PuppiWPiGammaDemo::ScPhase2PuppiWPiGammaDemo(const edm::ParameterSet &iConfig)
    : structPuppiToken_(consumes<OrbitCollection<l1Scouting::Puppi>>(iConfig.getParameter<edm::InputTag>("srcPuppi"))),
      structTkEmToken_(consumes<OrbitCollection<l1Scouting::TkEm>>(iConfig.getParameter<edm::InputTag>("srcTkEm"))) {
  produces<std::vector<unsigned>>("selectedBx");
  produces<l1ScoutingRun3::OrbitFlatTable>("wpigamma");
  cuts.minpt_pi = iConfig.getParameter<double>("ptPi");
  cuts.minpt_tkem = iConfig.getParameter<double>("ptTkEm");
  cuts.minmass = iConfig.getParameter<double>("minMass");
  cuts.maxmass = iConfig.getParameter<double>("maxMass");
  cuts.maxiso_pi = iConfig.getParameter<double>("relIsoPi");
  cuts.maxiso_tkem = iConfig.getParameter<double>("relIsoTkEm");
  cuts.mindeltar2 = std::pow(iConfig.getParameter<double>("minDeltaRPiTkEm"), 2);
  cuts.mindr2 = std::pow(iConfig.getParameter<double>("isolationMinDeltaRPi"), 2);
  cuts.maxdr2 = std::pow(iConfig.getParameter<double>("isolationMaxDeltaRPi"), 2);
  cuts.mindr2tkem = std::pow(iConfig.getParameter<double>("isolationMinDeltaRTkEm"), 2);
  cuts.maxdr2tkem = std::pow(iConfig.getParameter<double>("isolationMaxDeltaRTkEm"), 2);
}

ScPhase2PuppiWPiGammaDemo::~ScPhase2PuppiWPiGammaDemo() {};

void ScPhase2PuppiWPiGammaDemo::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2PuppiWPiGammaDemo::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::Puppi>> srcPuppi;
  iEvent.getByToken(structPuppiToken_, srcPuppi);

  edm::Handle<OrbitCollection<l1Scouting::TkEm>> srcTkEm;
  iEvent.getByToken(structTkEmToken_, srcTkEm);

  runObj(*srcPuppi, *srcTkEm, iEvent, countStruct_, passStruct_, "");
}

void ScPhase2PuppiWPiGammaDemo::endStream() {
  edm::LogImportant("ScPhase2AnalysisSummary") << "WPiGamma Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T, typename U>
void ScPhase2PuppiWPiGammaDemo::runObj(const OrbitCollection<T> &srcPuppi,
                                       const OrbitCollection<U> &srcTkEm,
                                       edm::Event &iEvent,
                                       unsigned long &nTry,
                                       unsigned long &nPass,
                                       const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();
  std::vector<float> masses;
  std::vector<uint8_t> i0s, i1s;  // i1 is the photon
  ROOT::RVec<unsigned int> ix;    // pions
  ROOT::RVec<unsigned int> ig;    // photons
  ROOT::RVec<unsigned int>
      iso;  //stores whether the Pi or photon passes isolation test so we don't calculate reliso twice
  std::array<unsigned int, 2> bestDoublet{{0, 0}};
  float bestDoubletScore, bestDoubletMass;
  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = srcPuppi.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    auto range2 = srcTkEm.bxIterator(bx);
    const U *cands2 = &range2.front();
    auto size2 = range2.size();

    ix.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all hadrons
      if ((std::abs(cands[i].pdgId()) == 211 or std::abs(cands[i].pdgId()) == 11)) {
        if (cands[i].pt() >= cuts.minpt_pi) {
          ix.push_back(i);
        }
      }
    }
    unsigned int npions = ix.size();
    if (npions < 1)
      continue;

    ig.clear();
    for (unsigned int i = 0; i < size2; ++i) {  //make list of all photons
      if (cands2[i].pt() >= cuts.minpt_tkem) {
        ig.push_back(i);
      }
    }
    unsigned int ngammas = ig.size();
    if (ngammas < 1)
      continue;

    iso.resize(npions + ngammas);  //gamma and Pi isolations
    std::fill(iso.begin(), iso.end(), 0);
    bestDoubletScore = 0;

    for (unsigned int i1 = 0; i1 < npions; ++i1) {
      if (!isolationPi(ix[i1], cands, size, iso[i1]))
        continue;  //ISOLATION test for pion

      for (unsigned int i2 = 0; i2 < ngammas; ++i2) {
        auto mass = (cands[ix[i1]].p4() + cands2[ig[i2]].p4()).mass();
        if (mass >= cuts.minmass and mass <= cuts.maxmass) {  //MASS test
          bool isop = isolationTkEm(
              cands2[ig[i2]].pt(), cands2[ig[i2]].eta(), cands2[ig[i2]].phi(), cands, size, iso[npions + i2]);
          bool pass_deltar = deltar(cands[ix[i1]].eta(), cands[ig[i2]].eta(), cands[ix[i1]].phi(), cands[ig[i2]].phi());
          if (isop == true and pass_deltar == true) {  //ISOLATION and DR tests
            float ptsum = cands[ix[i1]].pt() + cands2[ig[i2]].pt();
            if (ptsum > bestDoubletScore) {
              bestDoublet[0] = ix[i1];
              bestDoublet[1] = ig[i2];
              bestDoubletScore = ptsum;
              bestDoubletMass = mass;
            }
          }  //iso and dr tests
        }  //mass test
      }  //photon loop
    }  //pion loop

    if (bestDoubletScore > 0) {
      ret->emplace_back(bx);
      nPass++;
      masses.push_back(bestDoubletMass);
      i0s.push_back(bestDoublet[0]);
      i1s.push_back(bestDoublet[1]);
      bxOffsetsFiller.addBx(bx, 1);
    }
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "WPiGamma" + label, true);
  tab->addColumn<float>("mass", masses, "pi-gamma invariant mass");
  tab->addColumn<uint8_t>("i0", i0s, "Pion");
  tab->addColumn<uint8_t>("i1", i1s, "Photon");
  iEvent.put(std::move(tab), "wpigamma" + label);
}

//TEST functions
template <typename T>
bool ScPhase2PuppiWPiGammaDemo::isolationPi(unsigned int pidex1, const T *cands, unsigned int size) const {
  float psum = 0;
  float eta = cands[pidex1].eta(), phi = cands[pidex1].phi();
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    if (pidex1 == j)
      continue;
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.mindr2 && dr2 <= cuts.maxdr2)
      psum += cands[j].pt();
  }
  return (psum <= cuts.maxiso_pi * cands[pidex1].pt());
}

template <typename T>
bool ScPhase2PuppiWPiGammaDemo::isolationTkEm(float pt, float eta, float phi, const T *cands, unsigned int size) const {
  float psum = 0;
  for (unsigned int j = 0u; j < size; ++j) {  //loop over other particles
    float dr2 = reco::deltaR2(eta, phi, cands[j].eta(), cands[j].phi());
    if (dr2 >= cuts.mindr2tkem && dr2 <= cuts.maxdr2tkem)
      psum += cands[j].pt();
  }
  return (psum <= cuts.maxiso_tkem * pt);
}

bool ScPhase2PuppiWPiGammaDemo::deltar(float eta1, float eta2, float phi1, float phi2) const {
  float dr2 = reco::deltaR2(eta1, phi1, eta2, phi2);
  return (dr2 >= cuts.mindeltar2);
}

float ScPhase2PuppiWPiGammaDemo::doubletmass(const std::array<unsigned int, 2> &t,
                                             const float *pts,
                                             const float *etas,
                                             const float *phis) {
  ROOT::Math::PtEtaPhiMVector p1(pts[t[0]], etas[t[0]], phis[t[0]], 0.2);
  ROOT::Math::PtEtaPhiMVector p2(pts[t[1]], etas[t[1]], phis[t[1]], 0.00);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2PuppiWPiGammaDemo::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("srcPuppi");
  desc.add<edm::InputTag>("srcTkEm");
  desc.add<double>("ptPi", 25.);
  desc.add<double>("ptTkEm", 20.);
  desc.add<double>("minMass", 60.);
  desc.add<double>("maxMass", 100.);
  desc.add<double>("minDeltaRPiTkEm", 0.50);
  desc.add<double>("relIsoPi", 0.30);
  desc.add<double>("relIsoTkEm", 0.30);
  desc.add<double>("isolationMinDeltaRPi", 0.00);
  desc.add<double>("isolationMaxDeltaRPi", 0.50);
  desc.add<double>("isolationMinDeltaRTkEm", 0.02);
  desc.add<double>("isolationMaxDeltaRTkEm", 0.50);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2PuppiWPiGammaDemo);
