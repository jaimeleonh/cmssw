#ifndef Ntuplizer_h
#define Ntuplizer_h

// system include files
#include <memory>
#include <vector>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"

#include "TTree.h"
/***
// trigger include files
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

// jet/tau include files
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

//#include "DataFormats/PatCandidates/interface/MET.h"
//#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
//#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
***/
#include "TTree.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/CaloRecHit/interface/CaloCluster.h"
#include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "CondFormats/JetMETObjects/interface/FactorizedJetCorrector.h"
//#include "EgammaAnalysis/ElectronTools/interface/EnergyScaleCorrection_class.h"
//#include "HiggsAnalysis/HiggsTo2photons/interface/CiCPhotonID.h"
#include "JetMETCorrections/Modules/interface/JetResolution.h"
////#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"

using namespace std;

class NtupleMaker : public edm::one::EDAnalyzer<edm::one::SharedResources> {
    public:
        explicit NtupleMaker(const edm::ParameterSet&);
	~NtupleMaker();

	static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

    private:
	//virtual void beginJob() override;
	virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
	//virtual void endJob() override;
	
	void branchesTriggers(TTree*);
	void branchesEventInfo(TTree*);
	void branchesL1Taus(TTree*);
	void branchesL1Jets(TTree*);
        void branchesTaus(TTree*);
        void branchesJets(TTree*);

	void fillTriggers(const edm::Event&);
	void fillEventInfo(const edm::Event&);
	void fillL1Taus(const edm::Event&);
	void fillL1Jets(const edm::Event&);
        void fillTaus(const edm::Event&);
	void fillJets(const edm::Event&, const edm::EventSetup&);

	//-------------member data----------------//
	TTree* tree_; 


	bool fillingTriggers;
	bool fillingL1;
	bool fillingEventInfo;
	bool fillingTaus;
	bool fillingJets;

	bool development_; //had to add these in so the copied tau/jet files would play nice with my config file
	bool doGenParticles_;
	edm::EDGetTokenT<vector<reco::GenParticle> >    genParticlesCollection_;
	edm::EDGetTokenT<vector<pat::Tau> >             tauCollection_; 
        //edm::EDGetTokenT<vector<reco::PFTau>>		PFTauCollection_;
        // trigger primitives
        edm::EDGetTokenT<BXVector<l1t::Jet> >		jetTriggerPrimitives_;
	edm::EDGetTokenT<BXVector<l1t::Tau> >		tauTriggerPrimitives_;

	edm::EDGetTokenT<reco::VertexCollection>        vtxLabel_;
	edm::EDGetTokenT<double>                        rhoLabel_;
	edm::EDGetTokenT<edm::View<pat::Jet> >          jetsAK4Label_;
	JME::JetResolution            jetResolution_;
        JME::JetResolutionScaleFactor jetResolutionSF_;	////
	
	edm::EDGetTokenT<edm::TriggerResults> triggerResultToken_;
	edm::EDGetTokenT<trigger::TriggerEvent> triggerEventToken_;
	edm::EDGetTokenT<trigger::TriggerEventWithRefs> triggerEventWithRefsToken_;

};

#endif
