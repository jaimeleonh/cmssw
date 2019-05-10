#include "L1Trigger/DTPhase2Trigger/interface/MuonPathAssociator.h"
#include "L1Trigger/DTPhase2Trigger/interface/MuonPathAnalyzerPerSL.h"

using namespace edm;
using namespace std;



// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAssociator::MuonPathAssociator(const ParameterSet& pset) {
  // Obtention of parameters
  debug            = pset.getUntrackedParameter<Bool_t>("debug");
  dT0_correlate_TP = pset.getUntrackedParameter<double>("dT0_correlate_TP");
  if (debug) cout <<"MuonPathAssociator: constructor" << endl;
}


MuonPathAssociator::~MuonPathAssociator() {
  if (debug) cout <<"MuonPathAssociator: destructor" << endl;
}



// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAssociator::initialise(const edm::EventSetup& iEventSetup) {
  if(debug) cout << "MuonPathAssociator::initialiase" << endl;

  iEventSetup.get<MuonGeometryRecord>().get(dtGeo);//1103
}


void MuonPathAssociator::run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &muonpaths) {
  if (debug) cout <<"MuonPathAssociator: run" << endl;  
  
  for(auto mpath = muonpaths.begin();mpath!=muonpaths.end();++mpath) {
    if (debug)     
      cout << "MuonPathAssociator::run -> mpath (before association): "
	   << (*mpath)->getBxTimeValue() << " (" <<  (*mpath)->getBxTimeValue(0) << ","<< (*mpath)->getBxTimeValue(2) << ") " 
	   << (*mpath)->getHorizPos() << " (" <<  (*mpath)->getHorizPos(0) << ","<< (*mpath)->getBxTimeValue(2) << ") " 
	   << (*mpath)->getTanPhi() << " (" <<  (*mpath)->getTanPhi(0) << ","<< (*mpath)->getBxTimeValue(2) << ") " 
	   << (*mpath)->getPhi() << " (" <<  (*mpath)->getPhi(0) << ","<< (*mpath)->getPhi(2) << ") " 
	   << (*mpath)->getPhiB() << " (" <<  (*mpath)->getPhiB(0) << ","<< (*mpath)->getPhiB(2) << ") " 
	   << (*mpath)->getQuality() << " (" <<  (*mpath)->getQuality(0) << ","<< (*mpath)->getQuality(2) << ") " 
	   << endl;
    associate(*mpath);
    if (debug)     
      cout << "MuonPathAssociator::run -> SFG mpath (after association): "
	   << (*mpath)->getBxTimeValue() << " " 
	   << (*mpath)->getHorizPos() << " " 
	   << (*mpath)->getTanPhi() << " " 
	   << (*mpath)->getPhi() << " " 
	   << (*mpath)->getPhiB() << " " 
	   << (*mpath)->getQuality() << " " 
	   << endl;

  }

}

void MuonPathAssociator::finish() {
  if (debug) cout <<"MuonPathAssociator: finish" << endl;
};

void MuonPathAssociator::associate(MuonPath *mpath) {
  
  // First try to match 
  if (mpath->getNPrimitivesUp()>=3 && mpath->getNPrimitivesDown()>=3) {
    if(fabs(mpath->getBxTimeValue(0)-mpath->getBxTimeValue(2)) < dT0_correlate_TP) { //time match
      float PosSL1=mpath->getHorizPos(0);
      float PosSL3=mpath->getHorizPos(2);
      float NewSlope=(PosSL1-PosSL3)/23.5;     
      float MeanT0=(mpath->getBxTimeValue(0)+mpath->getBxTimeValue(2))/2;
      float MeanPos=(PosSL3+PosSL1)/2;
      float newChi2=(mpath->getChiSq(0)+mpath->getChiSq(2))*0.5;//to be recalculated
      MP_QUALITY quality=NOPATH;
      
      if (mpath->getQuality(0) <=LOWQ and mpath->getQuality(2) <=LOWQ)  quality=LOWLOWQ;
      if ((mpath->getQuality(0) >=HIGHQ and mpath->getQuality(2) <=LOWQ) or 
	  (mpath->getQuality(0) <=LOWQ and mpath->getQuality(2) >=HIGHQ))
	quality=HIGHLOWQ;
      if (mpath->getQuality(0) >=3 and mpath->getQuality(2) >=3)  quality=HIGHHIGHQ;
      
      DTChamberId ChId(mpath->getRawId());
      GlobalPoint jm_x_cmssw_global = dtGeo->chamber(ChId)->toGlobal(LocalPoint(MeanPos,0.,0.));//jm_x is already extrapolated to the middle of the SL
      int thisec = ChId.sector();
      float phi= jm_x_cmssw_global.phi()-0.5235988*(thisec-1);
      float psi=atan(NewSlope);
      float phiB=(hasPosRF(ChId.wheel(),ChId.sector())) ? psi-phi :-psi-phi ;
			
      mpath->setBxTimeValue(MeanT0);
      mpath->setTanPhi(NewSlope);
      mpath->setHorizPos(MeanPos);
      mpath->setPhi(phi);
      mpath->setPhiB(phiB);
      mpath->setChiSq(newChi2);
      mpath->setQuality(quality);
    }
  }
  else if (mpath->getNPrimitivesUp()>=3 && mpath->getNPrimitivesDown()<3 && mpath->getNPrimitivesDown()>0 ) {
    // IF this is not the case try to confirm with other SL: 
    mpath->setBxTimeValue(mpath->getBxTimeValue(2));
    mpath->setTanPhi(mpath->getTanPhi(2));
    mpath->setHorizPos(mpath->getHorizPos(2));
    mpath->setPhi(mpath->getPhi(2));
    mpath->setPhiB(mpath->getPhiB(2));
    mpath->setChiSq(mpath->getChiSq(2));

    if (mpath->getQuality(2) == HIGHQ) 
      mpath->setQuality(CHIGHQ);
    else if (mpath->getQuality(2) == LOWQ) 
      mpath->setQuality(CLOWQ);
    
  }
  else if (mpath->getNPrimitivesDown()>=3 && mpath->getNPrimitivesDown()<3 && mpath->getNPrimitivesDown()>0 ) {
    // IF this is not the case try to confirm with other SL: 
    mpath->setBxTimeValue(mpath->getBxTimeValue(2));
    mpath->setTanPhi(mpath->getTanPhi(2));
    mpath->setHorizPos(mpath->getHorizPos(2));
    mpath->setPhi(mpath->getPhi(2));
    mpath->setPhiB(mpath->getPhiB(2));
    mpath->setChiSq(mpath->getChiSq(2));
    mpath->setQuality(CHIGHQ);

    if (mpath->getQuality(0) == HIGHQ) 
      mpath->setQuality(CHIGHQ);
    else if (mpath->getQuality(0) == LOWQ) 
      mpath->setQuality(CLOWQ);
    
  }
  
}


