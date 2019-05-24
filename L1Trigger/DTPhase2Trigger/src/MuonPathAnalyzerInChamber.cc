#include "L1Trigger/DTPhase2Trigger/interface/MuonPathAnalyzerInChamber.h"
#include <cmath> 

using namespace edm;
using namespace std;



// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathAnalyzerInChamber::MuonPathAnalyzerInChamber(const ParameterSet& pset) :
  MuonPathAnalyzer(pset),
  bxTolerance(30),
  minQuality(LOWQGHOST),
  chiSquareThreshold(50)
{
  // Obtention of parameters
  debug         = pset.getUntrackedParameter<Bool_t>("debug");
  if (debug) cout <<"MuonPathAnalyzer: constructor" << endl;

  chi2Th = pset.getUntrackedParameter<double>("chi2Th");  
  setChiSquareThreshold(chi2Th*100.); 

  //z
  int rawId;
  z_filename = pset.getUntrackedParameter<std::string>("z_filename");
  std::ifstream ifin2(z_filename.c_str());
  double z;
  while (ifin2.good()){
    ifin2 >> rawId >> z;
    zinfo[rawId]=z;
  }

  //shift
  shift_filename = pset.getUntrackedParameter<std::string>("shift_filename");
  std::ifstream ifin3(shift_filename.c_str());
  double shift;
  while (ifin3.good()){
    ifin3 >> rawId >> shift;
    shiftinfo[rawId]=shift;
  }

  chosen_sl = pset.getUntrackedParameter<int>("trigger_with_sl");

  if(chosen_sl!=1 && chosen_sl!=3 && chosen_sl!=4){
    std::cout<<"chosen sl must be 1,3 or 4(both superlayers)"<<std::endl;
    assert(chosen_sl!=1 && chosen_sl!=3 && chosen_sl!=4); //4 means run using the two superlayers
  }
}

MuonPathAnalyzerInChamber::~MuonPathAnalyzerInChamber() {
  if (debug) cout <<"MuonPathAnalyzer: destructor" << endl;
}



// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathAnalyzerInChamber::initialise(const edm::EventSetup& iEventSetup) {
  if(debug) cout << "MuonPathAnalyzerInChamber::initialiase" << endl;
  iEventSetup.get<MuonGeometryRecord>().get(dtGeo);//1103
}

void MuonPathAnalyzerInChamber::run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &muonpaths, std::vector<MuonPath*> &outmuonpaths) {
  if (debug) cout <<"MuonPathAnalyzerInChamber: run" << endl;

  // fit per SL (need to allow for multiple outputs for a single mpath)
  for(auto muonpath = muonpaths.begin();muonpath!=muonpaths.end();++muonpath) {
    analyze(*muonpath, outmuonpaths);
  }

}

void MuonPathAnalyzerInChamber::finish() {
  if (debug) cout <<"MuonPathAnalyzer: finish" << endl;
};

const int MuonPathAnalyzerInChamber::LAYER_ARRANGEMENTS[4][3] = {
    {0, 1, 2}, {1, 2, 3},                       // Grupos consecutivos
    {0, 1, 3}, {0, 2, 3}                        // Grupos salteados
};


//------------------------------------------------------------------
//--- Métodos get / set
//------------------------------------------------------------------
void MuonPathAnalyzerInChamber::setBXTolerance(int t) { bxTolerance = t; }
int  MuonPathAnalyzerInChamber::getBXTolerance(void)  { return bxTolerance; }

void MuonPathAnalyzerInChamber::setChiSquareThreshold(float ch2Thr) {
    chiSquareThreshold = ch2Thr;
}

void MuonPathAnalyzerInChamber::setMinimumQuality(MP_QUALITY q) {
    if (minQuality >= LOWQGHOST) minQuality = q;
}
MP_QUALITY MuonPathAnalyzerInChamber::getMinimumQuality(void) { return minQuality; }


//------------------------------------------------------------------
//--- Métodos privados
//------------------------------------------------------------------
void MuonPathAnalyzerInChamber::analyze(MuonPath *inMPath,std::vector<MuonPath*>& outMPath) {
  debug=true;
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t starts"<<std::endl;
  
  // Clonamos el objeto analizado.
  if (debug) cout << inMPath->getNPrimitives() << endl;
  MuonPath *mPath = new MuonPath(*inMPath);
  
  if (debug) {
    for (int i=0; i<mPath->getNPrimitives(); i++) 
      std::cout << "DTp2::analyze, looking at mPath: " 
		<< mPath->getPrimitive(i)->getLayerId() << " , " 
		<< mPath->getPrimitive(i)->getSuperLayerId() << " , " 
		<< mPath->getPrimitive(i)->getChannelId() << " , " 
		<< mPath->getPrimitive(i)->getLaterality() << std::endl;
  }

  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t is Analyzable? "<<std::endl;
  if (!mPath->isAnalyzable())  return;
  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t yes it is analyzable "<<mPath->isAnalyzable()<<std::endl;
  
  // first of all, get info from primitives, so we can reduce the number of latereralities:
  buildLateralities(mPath);
  setCellLayout(mPath);  
  setWirePosAndTimeInMP(mPath);
  
  for (int i = 0; i < totalNumValLateralities; i++) {// LOOP for all lateralities:   
    calculateFitParameters(mPath,lateralities[i]);
  }
  return; 

  if ( mPath->getQuality() < minQuality ) return;
  
  /*
  // LOCATE MPATH (do we need this?) 
  int selected_Id=0;
  if     (inMPath->getPrimitive(0)->getTDCTime()!=-1) selected_Id= inMPath->getPrimitive(0)->getCameraId();
  else if(inMPath->getPrimitive(1)->getTDCTime()!=-1) selected_Id= inMPath->getPrimitive(1)->getCameraId(); 
  else if(inMPath->getPrimitive(2)->getTDCTime()!=-1) selected_Id= inMPath->getPrimitive(2)->getCameraId(); 
  else if(inMPath->getPrimitive(3)->getTDCTime()!=-1) selected_Id= inMPath->getPrimitive(3)->getCameraId(); 
  
  DTLayerId thisLId(selected_Id);
  if(debug) std::cout<<"Building up MuonPathSLId from rawId in the Primitive"<<std::endl;
  DTSuperLayerId MuonPathSLId(thisLId.wheel(),thisLId.station(),thisLId.sector(),thisLId.superLayer());
  if(debug) std::cout<<"The MuonPathSLId is"<<MuonPathSLId<<std::endl;
  DTWireId wireId(MuonPathSLId,2,1);
  /// END LOCATE MPATH

  if ( mPath->getQuality() < minQuality ) return;
  
  if(debug) {
    std::cout<<"DTp2:analyze \t\t\t\t min quality achievedCalidad: "<<mPath->getQuality()<<std::endl;
    for (int i = 0; i <= 3; i++)
      std::cout<<"DTp2:analyze \t\t\t\t  Capa: "<<mPath->getPrimitive(i)->getLayerId() <<" Canal: "<<mPath->getPrimitive(i)->getChannelId() <<" TDCTime: "<<mPath->getPrimitive(i)->getTDCTime()   <<std::endl;
  
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t Starting lateralities loop, totalNumValLateralities: "<<totalNumValLateralities<<std::endl;
  }
    
  
  for (int i = 0; i < totalNumValLateralities; i++) {//here
    if(debug) {
      std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<std::endl;
      std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" checking quality:"<<std::endl;
      std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" checking mPath Quality="<<mPath->getQuality()<<std::endl;   
      std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" latQuality[i].val="<<latQuality[i].valid<<std::endl;   
      std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" before if:"<<std::endl;
    }
    
    if (latQuality[i].valid and (((mPath->getQuality()==HIGHQ or mPath->getQuality()==HIGHQGHOST) and latQuality[i].quality==HIGHQ)
				 or
				 ((mPath->getQuality() == LOWQ or mPath->getQuality()==LOWQGHOST) and latQuality[i].quality==LOWQ))){
      
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" inside if"<<std::endl;
      mPath->setBxTimeValue(latQuality[i].bxValue);
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" settingLateralCombination"<<std::endl;
      mPath->setLateralComb(lateralities[i]);
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t laterality #- "<<i<<" done settingLateralCombination"<<std::endl;
      
      // Clonamos el objeto analizado.
      MuonPath *mpAux = new MuonPath(mPath);
      
      int idxHitNotValid = latQuality[i].invalidateHitIdx;
      if (idxHitNotValid >= 0) {
	delete mpAux->getPrimitive(idxHitNotValid);
	mpAux->setPrimitive(std::move(new DTPrimitive()), idxHitNotValid);
      }
      
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  calculating parameters "<<std::endl;
      calculatePathParameters(mpAux);
      if ((mpAux->getQuality() == HIGHQ or mpAux->getQuality() == HIGHQGHOST) && mpAux->getChiSq() > chiSquareThreshold) {//check this if!!!
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  HIGHQ or HIGHQGHOST but min chi2 or Q test not satisfied "<<std::endl;
      }				
      else{
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  inside else, returning values: "<<std::endl;
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  BX Time = "<<mpAux->getBxTimeValue()<<std::endl;
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  BX Id   = "<<mpAux->getBxNumId()<<std::endl;
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  XCoor   = "<<mpAux->getHorizPos()<<std::endl;
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  tan(Phi)= "<<mpAux->getTanPhi()<<std::endl;
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  chi2= "<<mpAux->getChiSq()<<std::endl;
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t  lateralities = "
			   <<" "<<mpAux->getLateralComb()[0]
			   <<" "<<mpAux->getLateralComb()[1]
			   <<" "<<mpAux->getLateralComb()[2]
			   <<" "<<mpAux->getLateralComb()[3]
			   <<std::endl;
	
	DTChamberId ChId(MuonPathSLId.wheel(),MuonPathSLId.station(),MuonPathSLId.sector());
	
	double jm_tanPhi=-1.*mpAux->getTanPhi(); //testing with this line
	double jm_x=(mpAux->getHorizPos()/10.)+shiftinfo[wireId.rawId()]; 
	//changing to chamber frame or reference:
	double jm_t0=mpAux->getBxTimeValue();		      
	//	  int quality= mpAux->getQuality();
	
	//computing phi and phiB
	double z=0;
	double z1=11.75;
	double z3=-1.*z1;
	if (ChId.station() == 3 or ChId.station() == 4){
	  z1=9.95;
	  z3=-13.55;
	}
	
	if(MuonPathSLId.superLayer()==1) z=z1;
	if(MuonPathSLId.superLayer()==3) z=z3;
	
	GlobalPoint jm_x_cmssw_global = dtGeo->chamber(MuonPathSLId)->toGlobal(LocalPoint(jm_x,0.,z));//jm_x is already extrapolated to the middle of the SL
	int thisec = MuonPathSLId.sector();
	if(thisec==13) thisec = 4;
	if(thisec==14) thisec = 10;
	double phi= jm_x_cmssw_global.phi()-0.5235988*(thisec-1);
	double psi=atan(jm_tanPhi);
	double phiB=hasPosRF(MuonPathSLId.wheel(),MuonPathSLId.sector()) ? psi-phi : -psi-phi ;
	double chi2= mpAux->getChiSq()*0.01;//in cmssw we want cm, 1 cm^2 = 100 mm^2
	
	if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaPrimitive at x="<<jm_x<<" tanPhi:"<<jm_tanPhi<<" t0:"<<jm_t0<<std::endl;
	  
	if(mpAux->getQuality() == HIGHQ or mpAux->getQuality() == HIGHQGHOST){//keep only the values with the best chi2 among lateralities
	  if(chi2<best_chi2){
	    chi2_jm_tanPhi=jm_tanPhi;
	      chi2_jm_x=(mpAux->getHorizPos()/10.)+shiftinfo[wireId.rawId()]; 
	      //chi2_jm_x=chi2_jm_x-(zinfo[wireId.rawId()]-0.65)*chi2_jm_tanPhi; //from SL to CH no needed for co
	      chi2_jm_t0=mpAux->getBxTimeValue();		      
	      chi2_phi=phi;
	      chi2_phiB=phiB;
	      chi2_chi2=chi2;
	      chi2_quality= mpAux->getQuality();
	  }
	}else{//write the metaprimitive in case no HIGHQ or HIGHQGHOST
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back metaprimitive no HIGHQ or HIGHQGHOST"<<std::endl;
	  metaPrimitives.push_back(metaPrimitive({MuonPathSLId.rawId(),jm_t0,jm_x,jm_tanPhi,phi,phiB,chi2,quality,
		    wi[0],tdc[0],
		    wi[1],tdc[1],
		    wi[2],tdc[2],
		    wi[3],tdc[3],
		    wi[4],tdc[4],
		    wi[5],tdc[5],
		    wi[6],tdc[6],
		    wi[7],tdc[7]
		    }));
	  if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  done pushing back metaprimitive no HIGHQ or HIGHQGHOST"<<std::endl;	
	}				
      }
      delete mpAux;
    }
    else{
      if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  latQuality[i].valid and (((mPath->getQuality()==HIGHQ or mPath->getQuality()==HIGHQGHOST) and latQuality[i].quality==HIGHQ) or  ((mPath->getQuality() == LOWQ or mPath->getQuality()==LOWQGHOST) and latQuality[i].quality==LOWQ)) not passed"<<std::endl;
    }
  }
  if(chi2_jm_tanPhi!=999){//
    if(debug) std::cout<<"DTp2:analyze \t\t\t\t\t\t\t\t  pushing back best chi2 metaPrimitive"<<std::endl;
     metaPrimitives.push_back(metaPrimitive({MuonPathSLId.rawId(),chi2_jm_t0,chi2_jm_x,chi2_jm_tanPhi,chi2_phi,chi2_phiB,chi2_chi2,chi2_quality,
       wi[0],tdc[0],
       wi[1],tdc[1],
       wi[2],tdc[2],
       wi[3],tdc[3],
       wi[4],tdc[4],
       wi[5],tdc[5],
       wi[6],tdc[6],
       wi[7],tdc[7]
       }));
    
       }
       }
       delete mPath;
       if(debug) std::cout<<"DTp2:analyze \t\t\t\t finishes"<<std::endl;
  */
}


void MuonPathAnalyzerInChamber::setCellLayout(MuonPath *mpath) {
  
  for (int i=0; i<=mpath->getNPrimitives(); i++) {
    if (mpath->getPrimitive(i)->isValidTime())   
      cellLayout[i] = mpath->getPrimitive(i)->getChannelId();
    else 
      cellLayout[i] = -99; 
  }
  
  // putting info back into the mpath:
  mpath->setCellHorizontalLayout(cellLayout);
  for (int i=0; i<=mpath->getNPrimitives(); i++){
    if (cellLayout[i]>=0) {
      mpath->setBaseChannelId(cellLayout[i]);
      break;
    }
  }
}

/**
 * For a combination of up to 8 cells, build all the lateralities to be tested,
 * and discards all  construye de forma automática todas las posibles
 * combinaciones de lateralidad (LLLL, LLRL,...) que son compatibles con una
 * trayectoria recta. Es decir, la partícula no hace un zig-zag entre los hilos
 * de diferentes celdas, al pasar de una a otra.
 */
void MuonPathAnalyzerInChamber::buildLateralities(MuonPath *mpath) {
  
  if (debug) cout << "MPAnalyzer::buildLateralities << setLateralitiesFromPrims " << endl;
  mpath->setLateralCombFromPrimitives();
  
  totalNumValLateralities = 0;
  lateralities.clear();
  latQuality.clear();
  
  /* We generate all the possible laterality combinations compatible with the built 
     group in the previous step*/
  lateralities.push_back(TLateralities());
  for (int ilat = 0; ilat <  NLayers; ilat++) {
    // Get value from input
    LATERAL_CASES lr = (mpath->getLateralComb())[ilat];
    if (debug) std::cout << "[DEBUG] Input[" << ilat << "]: " << lr << std::endl;
    
    // If left/right fill number
    if (lr != NONE ) {
      if (debug) std::cout << "[DEBUG]   - Adding it to " << lateralities.size() << " lists..." << std::endl;
      for (unsigned int iall = 0; iall < lateralities.size(); iall++) {
        lateralities[iall][ilat]  = lr;
	
      }
    }
    // both possibilites
    else {
      // Get the number of possible options now
      auto ncurrentoptions = lateralities.size();
      
      // Duplicate them
      std::cout << "[DEBUG]   - Duplicating " << ncurrentoptions << " lists..." << std::endl;
      copy(lateralities.begin(), lateralities.end(), back_inserter(lateralities));     
      std::cout << "[DEBUG]   - Now we have " << lateralities.size() << " lists..." << std::endl;

      // Asign LEFT to first ncurrentoptions and RIGHT to the last
      for (unsigned int iall = 0; iall < ncurrentoptions; iall++) {
        lateralities[iall][ilat]  = LEFT;
        lateralities[iall+ncurrentoptions][ilat]  = RIGHT;
      }
    } // else
  } // Iterate over input array
  
  totalNumValLateralities = (int) lateralities.size(); 
  for (unsigned int iall=0; iall<lateralities.size(); iall++) {
    latQuality.push_back(LATQ_TYPE());

    latQuality[iall].valid            = false;
    latQuality[iall].bxValue          = 0;
    latQuality[iall].quality          = NOPATH;
    latQuality[iall].invalidateHitIdx = -1;
  }

  if (totalNumValLateralities>33) {
    // ADD PROTECTION!
    cout << "[WARNING]: TOO MANY LATERALITIES TO CHECK !!" << endl;
    cout << "[WARNING]: skipping this muon" << endl;
    lateralities.clear();
    latQuality.clear();
    totalNumValLateralities = 0;
  }
  
  // Dump values
  if (debug) {
    for (unsigned int iall = 0; iall < lateralities.size(); iall++) {
      std::cout << iall << " -> [";
      for (int ilat = 0; ilat < NLayers; ilat++) {
	if (ilat != 0)
	  std::cout << ",";
	std::cout << lateralities[iall][ilat];
      }
      std::cout << "]" << std::endl;
    } 
  }
}
void MuonPathAnalyzerInChamber::setWirePosAndTimeInMP(MuonPath *mpath){
  
  float delta = 42000; //um
  float zwire[8]={-13.7, -12.4, -11.1, -9.8002, 9.79999, 11.1, 12.4, 13.7};
  for (int i=0; i<mpath->getNPrimitives(); i++){ 
    if (mpath->getPrimitive(i)->isValidTime())  {
      mpath->setXWirePos(mpath->getPrimitive(i)->getChannelId() +0.5*(double)(i%2) * delta,i);
      mpath->setZWirePos(zwire[i],i);
      mpath->settWireTDC(mpath->getPrimitive(i)->getTDCTime()*DRIFT_SPEED,i);
    }
    else {
      mpath->setXWirePos(0.,i);
      mpath->setZWirePos(0.,i);
      mpath->settWireTDC(-1*DRIFT_SPEED,i);
    }
  }
}
void MuonPathAnalyzerInChamber::calculateFitParameters(MuonPath *mpath, TLateralities laterality) {
  
  // First prepare mpath for fit: 
  int NMissingHits=0;
  float xwire[8],zwire[8],tTDCvdrift[8];
  int present_layer[8];
  double b[8];
  for (int i=0; i<8; i++){ 
    xwire[i]      = mpath->getXWirePos(i); 
    zwire[i]      = mpath->getZWirePos(i);
    tTDCvdrift[i] = mpath->gettWireTDC(i);
    b[i]          = 1;
    if (xwire[i]==0) { 
      present_layer[i]=0;
      NMissingHits++;      
    }
    else             
      present_layer[i]=1;
  }
  
  //// NOW Start FITTING:  
  
  // fill hit position
  float xhit[8];
  for  (int lay=0; lay<8; lay++){
    if (debug) cout << "In fitPerLat " << lay << " xwire " << xwire[lay] << " zwire "<< zwire[lay]<<" tTDCvdrift "<< tTDCvdrift[lay]<< endl;
    xhit[lay]=xwire[lay]-1*laterality[lay]*tTDCvdrift[lay];
    if (debug) cout << "In fitPerLat " << lay << " xhit "<< xhit[lay]<< endl;
  }  
      
  //Proceed with calculation of fit parameters
  double cbscal={0.000000d};
  double zbscal={0.000000d};
  double czscal={0.000000d};
  double bbscal={0.000000d};
  double zzscal={0.000000d};
  double ccscal={0.000000d};
  
  for  (int lay=0; lay<8; lay++){
    if (debug) cout<< " For layer " << lay+1 << " xwire[lay] " << xwire[lay] << " zwire " << zwire[lay] << " b " << b[lay] << endl;
    if (debug) cout<< " xhit[lat][lay] " << xhit[lay] << endl;
    cbscal=laterality[lay]*b[lay]+cbscal;
    zbscal=zwire[lay]*b[lay]+zbscal; //it actually does not depend on laterality
    czscal=laterality[lay]*zwire[lay]+czscal;
    
    bbscal=b[lay]*b[lay]+bbscal; //it actually does not depend on laterality
    zzscal=zwire[lay]*zwire[lay]+zzscal; //it actually does not depend on laterality
    ccscal=laterality[lay]*laterality[lay]+ccscal;
  }
  
  
  double cz= {0.000000d};
  double cb= {0.000000d};
  double zb= {0.000000d};
  double zc= {0.000000d};
  double bc= {0.000000d};
  double bz= {0.000000d};
  
  cz=(cbscal*zbscal-czscal*bbscal)/(zzscal*bbscal-zbscal*zbscal);
  cb=(czscal*zbscal-cbscal*zzscal)/(zzscal*bbscal-zbscal*zbscal);
  
  zb=(czscal*cbscal-zbscal*ccscal)/(bbscal*ccscal-cbscal*cbscal);
  zc=(zbscal*cbscal-czscal*bbscal)/(bbscal*ccscal-cbscal*cbscal);
  
  bc=(zbscal*czscal-cbscal*zzscal)/(ccscal*zzscal-czscal*czscal);
  bz=(cbscal*czscal-zbscal*ccscal)/(ccscal*zzscal-czscal*czscal);
  
  
  double c_tilde[8]; 
  double z_tilde[8];
  double b_tilde[8];
  
  for  (int lay=0; lay<8; lay++){
    c_tilde[lay]=laterality[lay]+cz*zwire[lay]+cb*b[lay];       	
    z_tilde[lay]=zwire[lay]+zb*b[lay]+zc*laterality[lay];
    b_tilde[lay]=b[lay]+bc*laterality[lay]+bz*zwire[lay];
    
  }
  
  //Calculate results per lat
  double xctilde={0.000000d};
  double xztilde={0.000000d};
  double xbtilde={0.000000d};
  double ctildectilde={0.000000d};
  double ztildeztilde={0.000000d};
  double btildebtilde={0.000000d};
  
  double rect0vdrift={0.000000d};
  double recslope={0.000000d};
  double recpos={0.000000d};
  
  for  (int lay=0; lay<8; lay++){
    xctilde=xhit[lay]*c_tilde[lay]+xctilde;
    ctildectilde=c_tilde[lay]*c_tilde[lay]+ctildectilde;
    xztilde=xhit[lay]*z_tilde[lay]+xztilde;
    ztildeztilde=z_tilde[lay]*z_tilde[lay]+ztildeztilde;
    xbtilde=xhit[lay]*b_tilde[lay]+xbtilde;
    btildebtilde=b_tilde[lay]*b_tilde[lay]+btildebtilde;
  }
  
  //Results for t0vdrift (BX), slope and position per lat
  rect0vdrift=xctilde/ctildectilde;
  recslope=xztilde/ztildeztilde;
  recpos=xbtilde/btildebtilde;
  if(debug) {
    cout<< " In fitPerLat Reconstructed values per lat " << " rect0vdrift "<< rect0vdrift;
    cout <<"rect0 "<< rect0vdrift/DRIFT_SPEED <<" recBX " << rect0vdrift/DRIFT_SPEED/25 << " recslope " << recslope << " recpos " << recpos  << endl;
  }
  
  //Get t*v and residuals per layer, and chi2 per laterality
  double rectdriftvdrift[8]={0.000000d};
  double recres[8]={0.000000d};
  double recchi2={0.000000d};
  int sign_tdriftvdrift={0};    
  int incell_tdriftvdrift={0};    
  int physical_slope={0}; 
  
  for  (int lay=0; lay<8; lay++){
    rectdriftvdrift[lay]= tTDCvdrift[lay]- rect0vdrift;
    recres[lay]=xhit[lay]-zwire[lay]*recslope-b[lay]*recpos-laterality[lay]*rect0vdrift;
    
    if ((present_layer[lay]==1)&&(rectdriftvdrift[lay] <-0.1)) sign_tdriftvdrift=-1;		  
    if ((present_layer[lay]==1)&&(abs(rectdriftvdrift[lay]) >2.15)) incell_tdriftvdrift=-1; //Changed to 2.11 to account for resolution effects		  
  }
  if (abs(recslope)>1)  physical_slope=-1;
  
  if (physical_slope==-1 && debug)  cout << "Combination with UNPHYSICAL slope " <<endl;
  if (sign_tdriftvdrift==-1 && debug) cout << "Combination with negative tdrift-vdrift " <<endl;
  if (incell_tdriftvdrift==-1 && debug) cout << "Combination with tdrift-vdrift larger than half cell " <<endl;
  
  for  (int lay=0; lay<8; lay++){
    recchi2=recres[lay]*recres[lay] + recchi2;
  }
  if(debug) cout << "In fitPerLat Chi2 " << recchi2 << " with sign " << sign_tdriftvdrift << " within cell " << incell_tdriftvdrift << " physical_slope "<< physical_slope << endl;
  
  // LATERALITY IS VALID... 
  if(!(sign_tdriftvdrift==-1) && !(incell_tdriftvdrift==-1) && !(physical_slope==-1)){
    mpath->setBxTimeValue(rect0vdrift/DRIFT_SPEED);
    mpath->setTanPhi(recslope);
    mpath->setHorizPos(recpos);
    mpath->setChiSq(recchi2);
    evaluateQuality(mpath);
    
    if(debug) cout << "In fitPerLat " << "t0 " <<  mpath->getBxTimeValue() <<" slope " << mpath->getTanPhi() <<" pos "<< mpath->getHorizPos() <<" chi2 "<< mpath->getChiSq() << endl;
  }
    
}
/**
 * Recorre las calidades calculadas para todas las combinaciones de lateralidad
 * válidas, para determinar la calidad final asignada al "MuonPath" con el que
 * se está trabajando.
 */
void MuonPathAnalyzerInChamber::evaluateQuality(MuonPath *mPath) {
  
  // here
  int totalHighQ = 0, totalLowQ = 0;
    
  if(debug) {
    std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t En evaluatePathQuality Evaluando PathQ. Celda base: "<<mPath->getBaseChannelId()<<std::endl;
    std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t Total lateralidades: "<<totalNumValLateralities<<std::endl;
  }
  
  // Por defecto.
  mPath->setQuality(NOPATH);
  
  for (int latIdx = 0; latIdx < totalNumValLateralities; latIdx++) {
    if(debug) {
      std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t Analizando combinacion de lateralidad: " ;
      for (int i=0; i<NLayers; i++)  std::cout <<lateralities[latIdx][i]<<" ";
      std::cout << std::endl; 
    }
    
    if (latQuality[latIdx].quality == HIGHQ) {
      totalHighQ++;
      if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t\t Lateralidad HIGHQ"<<std::endl;
    }
    if (latQuality[latIdx].quality == LOWQ) {
      totalLowQ++;
      if(debug) std::cout<<"DTp2:evaluatePathQuality \t\t\t\t\t\t Lateralidad LOWQ"<<std::endl;
    }
  }
  /*
   * Establecimiento de la calidad.
   */
  if (totalHighQ == 1) {
    mPath->setQuality(HIGHQ);
  }
  else if (totalHighQ > 1) {
    mPath->setQuality(HIGHQGHOST);
  }
  else if (totalLowQ == 1) {
    mPath->setQuality(LOWQ);
  }
  else if (totalLowQ > 1) {
    mPath->setQuality(LOWQGHOST);
  }
}

