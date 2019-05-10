#include "L1Trigger/DTPhase2Trigger/interface/MuonPathFilter.h"

using namespace edm;
using namespace std;



// ============================================================================
// Constructors and destructor
// ============================================================================
MuonPathFilter::MuonPathFilter(const ParameterSet& pset) {
  // Obtention of parameters
  debug         = pset.getUntrackedParameter<Bool_t>("debug");
  filter_primos = pset.getUntrackedParameter<bool>("filter_primos");
  tanPhiTh      = pset.getUntrackedParameter<double>("tanPhiTh");
  if (debug) cout <<"MuonPathFilter: constructor" << endl;
}


MuonPathFilter::~MuonPathFilter() {
  if (debug) cout <<"MuonPathFilter: destructor" << endl;
}



// ============================================================================
// Main methods (initialise, run, finish)
// ============================================================================
void MuonPathFilter::initialise(const edm::EventSetup& iEventSetup) {
  if(debug) cout << "MuonPathFilter::initialiase" << endl;
}

int MuonPathFilter::arePrimos(metaPrimitive primera, metaPrimitive segunda) {
  if(primera.rawId!=segunda.rawId) return 0;
  if(primera.wi1==segunda.wi1 and primera.tdc1==segunda.tdc1 and primera.wi1!=-1 and primera.tdc1!=-1) return 1;
  if(primera.wi2==segunda.wi2 and primera.tdc2==segunda.tdc2 and primera.wi2!=-1 and primera.tdc2!=-1) return 2;
  if(primera.wi3==segunda.wi3 and primera.tdc3==segunda.tdc3 and primera.wi3!=-1 and primera.tdc3!=-1) return 3;
  if(primera.wi4==segunda.wi4 and primera.tdc4==segunda.tdc4 and primera.wi4!=-1 and primera.tdc4!=-1) return 4;
  return 0;
}
int MuonPathFilter::rango(metaPrimitive primera) {
  int rango=0;
  if(primera.wi1!=-1)rango++;
  if(primera.wi2!=-1)rango++;
  if(primera.wi3!=-1)rango++;
  if(primera.wi4!=-1)rango++;
  return rango;
}

void MuonPathFilter::run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<metaPrimitive> &inMPath, std::vector<metaPrimitive> &outMpath) {
  
  if (debug) cout <<"MuonPathFilter: run" << endl;  
  if(filter_primos){
    if(debug) std::cout<<"filtering: starting primos filtering"<<std::endl;    
    
    int primo_index=0;
    bool oneof4=false;
    if(inMPath.size()==1){
      if(debug){
	std::cout<<"filtering:";
	std::cout<<" \t is:"<<0<<" "<<primo_index<<" "<<" "<<oneof4<<std::endl;
      }
      if(fabs(inMPath[0].tanPhi)<tanPhiTh){
	outMpath.push_back(inMPath[0]);
	if(debug)std::cout<<"filtering: kept1 i="<<0<<std::endl;
      }
    }
    else {
      for(int i=1; i<int(inMPath.size()); i++){ 
	if(fabs(inMPath[i].tanPhi)>tanPhiTh) continue;
	if(rango(inMPath[i])==4)oneof4=true;
	if(debug){
	  std::cout<<"filtering:";
	  std::cout<<" \t is:"<<i<<" "<<primo_index<<" "<<" "<<oneof4<<std::endl;
	}
	if(arePrimos(inMPath[i],inMPath[i-1])!=0  and arePrimos(inMPath[i],inMPath[i-primo_index-1])!=0){ 
	  primo_index++;
	}
	else{
	  if(primo_index==0){
	    outMpath.push_back(inMPath[i]);
	    if(debug)std::cout<<"filtering: kept2 i="<<i<<std::endl;
	  }
	  else{
	    if(oneof4){
	      double minchi2=99999;
	      int selected_i=0;
	      for(int j=i-1;j>=i-primo_index-1;j--){
		if(rango(inMPath[j])!=4) continue;
		if(minchi2>inMPath[j].chi2){
		  minchi2=inMPath[j].chi2;
		  selected_i=j;
		}
	      }
	      outMpath.push_back(inMPath[selected_i]);
	      if(debug)std::cout<<"filtering: kept4 i="<<selected_i<<std::endl;
	    }
	    else{
	      for(int j=i-1;j>=i-primo_index-1;j--){
		outMpath.push_back(inMPath[j]);
		if(debug)std::cout<<"filtering: kept3 i="<<j<<std::endl;
	      }
	    }
	  }
	  primo_index=0;
	  oneof4=false;
	}
      }
    }
  }
  else{
    for (size_t i=0; i<inMPath.size(); i++){ 
      if(fabs(inMPath[i].tanPhi)>tanPhiTh) continue;
      outMpath.push_back(inMPath[i]); 
    }
  }
}

void MuonPathFilter::finish() {
  if (debug) cout <<"MuonPathFilter: finish" << endl;
};


