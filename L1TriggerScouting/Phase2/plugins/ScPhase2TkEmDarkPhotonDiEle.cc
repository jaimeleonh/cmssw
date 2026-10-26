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

#include <ROOT/RVec.hxx>
#include <Math/Vector4D.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <algorithm>
#include <array>
#include <iostream>

class ScPhase2TkEmDarkPhotonDiEle : public edm::stream::EDProducer<> {
public:
  explicit ScPhase2TkEmDarkPhotonDiEle(const edm::ParameterSet &);
  ~ScPhase2TkEmDarkPhotonDiEle() override;
  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  void beginStream(edm::StreamID) override;
  void produce(edm::Event &, const edm::EventSetup &) override;
  void endStream() override;
  template <typename T>
  void runObj(const OrbitCollection<T> &src,
              edm::Event &out,
              unsigned long &nTry,
              unsigned long &nPass,
              const std::string &bxLabel);

  edm::EDGetTokenT<OrbitCollection<l1Scouting::TkEle>> structTkEleToken_;

  struct Cuts {
    std::array<float, 2> minpt;
    std::array<float, 2> minid;
    float maxeta = 1.479;
    float maxdz = 1;
    float maxRelIso = 0.4;
  } cuts;

  template <typename T>
  static float pairmass(const std::array<unsigned int, 2> &t, const T *cands);

  unsigned long countStruct_;
  unsigned long passStruct_;
};

ScPhase2TkEmDarkPhotonDiEle::ScPhase2TkEmDarkPhotonDiEle(const edm::ParameterSet &iConfig) {
  structTkEleToken_ = consumes<OrbitCollection<l1Scouting::TkEle>>(iConfig.getParameter<edm::InputTag>("src"));
  produces<std::vector<unsigned>>("selectedBx");
  produces<l1ScoutingRun3::OrbitFlatTable>("zdee");
  auto minPts = iConfig.getParameter<std::vector<double>>("ptMin");
  if (minPts.size() != 2) {
    throw cms::Exception("Configuration")
        << "ptMin should be a vector of size 2, for the leading and subleading electron";
  }
  std::copy_n(minPts.begin(), 2, cuts.minpt.begin());
  auto minIds = iConfig.getParameter<std::vector<double>>("idScore");
  if (minIds.size() != 2) {
    throw cms::Exception("Configuration")
        << "idScore should be a vector of size 2, for the leading and subleading electron";
  }
  std::copy_n(minIds.begin(), 2, cuts.minid.begin());
  cuts.maxRelIso = iConfig.getParameter<double>("relIso");
  cuts.maxeta = iConfig.getParameter<double>("etaMax");
  cuts.maxdz = iConfig.getParameter<double>("dzMax");
}

ScPhase2TkEmDarkPhotonDiEle::~ScPhase2TkEmDarkPhotonDiEle() {};

void ScPhase2TkEmDarkPhotonDiEle::beginStream(edm::StreamID) {
  countStruct_ = 0;
  passStruct_ = 0;
}

void ScPhase2TkEmDarkPhotonDiEle::produce(edm::Event &iEvent, const edm::EventSetup &iSetup) {
  edm::Handle<OrbitCollection<l1Scouting::TkEle>> srcTkEle;
  iEvent.getByToken(structTkEleToken_, srcTkEle);

  runObj(*srcTkEle, iEvent, countStruct_, passStruct_, "");
}

void ScPhase2TkEmDarkPhotonDiEle::endStream() {
  edm::LogImportant("ScPhase2AnalysisSummary") << "zdee Struct analysis: " << countStruct_ << " -> " << passStruct_;
}

template <typename T>
void ScPhase2TkEmDarkPhotonDiEle::runObj(const OrbitCollection<T> &srcTkEle,
                                         edm::Event &iEvent,
                                         unsigned long &nTry,
                                         unsigned long &nPass,
                                         const std::string &label) {
  l1ScoutingRun3::BxOffsetsFillter bxOffsetsFiller;
  bxOffsetsFiller.start();
  auto ret = std::make_unique<std::vector<unsigned>>();

  std::vector<float> masses;
  ROOT::RVec<unsigned int> iEle;
  std::array<unsigned int, 2> bestPair{{0, 0}};

  bool bestPairFound;
  float minDZ;

  for (unsigned int bx = 1; bx <= OrbitCollection<T>::NBX; ++bx) {
    nTry++;
    auto range = srcTkEle.bxIterator(bx);
    const T *cands = &range.front();
    auto size = range.size();

    // Select events with two or more electrons with pT > 5 GeV and in barrel
    iEle.clear();
    for (unsigned int i = 0; i < size; ++i) {  //make list of all electrons
      if ((cands[i].pt() >= cuts.minpt[1]) && (std::abs(cands[i].eta()) <= cuts.maxeta) &&
          (cands[i].idScore() >= std::min(cuts.minid[0], cuts.minid[1]))) {
        iEle.push_back(i);
      }
    }

    unsigned int nEle = iEle.size();
    if (nEle < 2)
      continue;

    // Loop over possible ee pairs; get the best pair
    bestPairFound = false;
    minDZ = cuts.maxdz;
    for (unsigned int i1 = 0; i1 < nEle; ++i1) {
      if (cands[iEle[i1]].pt() < cuts.minpt[0])
        continue;

      if (cands[iEle[i1]].idScore() < cuts.minid[0])
        continue;

      if (cands[iEle[i1]].isolation() > cands[iEle[i1]].pt() * cuts.maxRelIso)
        continue;

      for (unsigned int i2 = i1 + 1; i2 < nEle; ++i2) {
        // pt cut was already applied when filling iEle,
        // but we need to apply id cuts
        if ((cands[iEle[i2]].idScore() < cuts.minid[1]))
          continue;

        // isolation cut
        if (cands[iEle[i2]].isolation() > cands[iEle[i2]].pt() * cuts.maxRelIso)
          continue;

        // OS requirement
        if (!(cands[iEle[i1]].charge() * cands[iEle[i2]].charge() < 0))
          continue;

        // dz requirement
        float pairDZ = std::abs(cands[iEle[i1]].z0() - cands[iEle[i2]].z0());
        if (pairDZ >= cuts.maxdz)
          continue;

        if (pairDZ < minDZ) {
          bestPair[0] = iEle[i1];
          bestPair[1] = iEle[i2];
          minDZ = pairDZ;
          bestPairFound = true;
        }
        // Find the one with the max dPhi
      }
    }
    if (!bestPairFound)
      continue;

    // Best ee pair mass
    auto mass = pairmass(bestPair, cands);

    ret->emplace_back(bx);
    nPass++;

    masses.push_back(mass);
    bxOffsetsFiller.addBx(bx, 1);
  }  // loop on BXs

  iEvent.put(std::move(ret), "selectedBx" + label);
  // now we make the table
  auto bxOffsets = bxOffsetsFiller.done();
  auto tab = std::make_unique<l1ScoutingRun3::OrbitFlatTable>(bxOffsets, "Zdee" + label, true);

  tab->addColumn<float>("mass", masses, "di-electron invariant mass");

  iEvent.put(std::move(tab), "zdee" + label);
}

template <typename T>
float ScPhase2TkEmDarkPhotonDiEle::pairmass(const std::array<unsigned int, 2> &t, const T *cands) {
  const float eleMass = 0.51e-3;
  ROOT::Math::PtEtaPhiMVector p1(cands[t[0]].pt(), cands[t[0]].eta(), cands[t[0]].phi(), eleMass);
  ROOT::Math::PtEtaPhiMVector p2(cands[t[1]].pt(), cands[t[1]].eta(), cands[t[1]].phi(), eleMass);
  float mass = (p1 + p2).M();
  return mass;
}

void ScPhase2TkEmDarkPhotonDiEle::fillDescriptions(edm::ConfigurationDescriptions &descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("src");
  desc.add<std::vector<double>>("ptMin", {5.0, 4.0});
  desc.add<std::vector<double>>("idScore", {0.0, 0.0});
  desc.add<double>("etaMax", 1.479);
  desc.add<double>("dzMax", 1.0);
  desc.add<double>("relIso", 0.4);
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(ScPhase2TkEmDarkPhotonDiEle);
