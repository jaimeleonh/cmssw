#include "NtupleMaker/NtupleMaker/plugins/NtupleMaker.h"

//using namespace std;

NtupleMaker::NtupleMaker(const edm::ParameterSet& iConfig) :
    fillingTriggers(iConfig.getUntrackedParameter<bool>("fillingTriggers")),
    fillingL1(iConfig.getUntrackedParameter<bool>("fillingL1")),
    fillingEventInfo(iConfig.getUntrackedParameter<bool>("fillingEventInfo")),
    fillingTaus(iConfig.getUntrackedParameter<bool>("fillingTaus")),
    fillingJets(iConfig.getUntrackedParameter<bool>("fillingJets")),

    development_(iConfig.getUntrackedParameter<bool>("development")),
    doGenParticles_(iConfig.getUntrackedParameter<bool>("doGenParticles")),
    genParticlesCollection_(consumes<vector<reco::GenParticle>>(iConfig.getUntrackedParameter<edm::InputTag>("genParticleSrc"))),
    tauCollection_(consumes<vector<pat::Tau>>(iConfig.getUntrackedParameter<edm::InputTag>("tauSrc"))),
    //PFTauCollection_(consumes<vector<reco::PFTau>>(iConfig.getUntrackedParameter<edm::InputTag>("PFTauCollection"))),
    // L1 Trigger Primitives
    jetTriggerPrimitives_(consumes<BXVector<l1t::Jet>>(iConfig.getUntrackedParameter<edm::InputTag>("jetTriggerPrimitives"))),
    tauTriggerPrimitives_(consumes<BXVector<l1t::Tau>>(iConfig.getUntrackedParameter<edm::InputTag>("tauTriggerPrimitives"))),

    vtxLabel_(consumes<reco::VertexCollection>(iConfig.getUntrackedParameter<edm::InputTag>("VtxLabel"))),
    rhoLabel_(consumes<double>(iConfig.getUntrackedParameter<edm::InputTag>("rhoLabel"))),
    jetsAK4Label_(consumes<edm::View<pat::Jet>>(iConfig.getUntrackedParameter<edm::InputTag>("ak4JetSrc"))),

    triggerResultToken_(consumes<edm::TriggerResults>(iConfig.getUntrackedParameter<edm::InputTag>("triggerResults"))),
    triggerEventToken_(consumes<trigger::TriggerEvent>(iConfig.getUntrackedParameter<edm::InputTag>("triggerEvent"))),

    triggerEventWithRefsToken_(consumes<trigger::TriggerEventWithRefs>(iConfig.getUntrackedParameter<edm::InputTag>("triggerEventWithRefs")))

{

    edm::Service<TFileService> fs;
    tree_ = fs->make<TTree>("vbf", "vbf");

    if(fillingTriggers) branchesTriggers(tree_);
    if(fillingEventInfo) branchesEventInfo(tree_);
    if(fillingL1) branchesL1Taus(tree_);
    if(fillingL1) branchesL1Jets(tree_);
    if(fillingTaus) branchesTaus(tree_);
    if(fillingJets) branchesJets(tree_);

}

NtupleMaker::~NtupleMaker(){
//destructor
}

void NtupleMaker::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    //using namespace edm;
    //if(doGenParticles_){
        //jetResolution_   = JME::JetResolution::get(iSetup, "AK4PFchs_pt");
        //jetResolutionSF_ = JME::JetResolutionScaleFactor::get(iSetup, "AK4PFchs");
        //AK8jetResolution_   = JME::JetResolution::get(es, "AK8PFchs_pt");
        //AK8jetResolutionSF_ = JME::JetResolutionScaleFactor::get(es, "AK8PFchs");
    //}

    if(fillingTriggers) fillTriggers(iEvent);
    if(fillingEventInfo) fillEventInfo(iEvent);
    if(fillingL1) fillL1Taus(iEvent);
    if(fillingL1) fillL1Jets(iEvent);
    if(fillingTaus) fillTaus(iEvent);
    if(fillingJets) fillJets(iEvent, iSetup);

    tree_->Fill();
}

void NtupleMaker::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(NtupleMaker);
