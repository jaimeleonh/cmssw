#include "L1Trigger/DTPhase2Trigger/interface/DTTrigPhase2Prod.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Run.h"
#include "L1TriggerConfig/DTTPGConfig/interface/DTConfigManager.h"

#include "L1TriggerConfig/DTTPGConfig/interface/DTConfigManagerRcd.h"
#include "L1Trigger/DTSectorCollector/interface/DTSectCollPhSegm.h"
#include "L1Trigger/DTSectorCollector/interface/DTSectCollThSegm.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1Phase2MuDTPhDigi.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThContainer.h"
#include "DataFormats/L1DTTrackFinder/interface/L1MuDTChambThDigi.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/DTGeometry/interface/DTGeometry.h>
#include "Geometry/DTGeometry/interface/DTLayer.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"


// DT trigger GeomUtils
#include "DQM/DTMonitorModule/interface/DTTrigGeomUtils.h"

#include "CalibMuon/DTDigiSync/interface/DTTTrigBaseSync.h"
#include "CalibMuon/DTDigiSync/interface/DTTTrigSyncFactory.h"

#include <iostream>
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"


#include "L1Trigger/DTPhase2Trigger/interface/muonpath.h"

using namespace edm;
using namespace std;

typedef vector<DTSectCollPhSegm>  SectCollPhiColl;
typedef SectCollPhiColl::const_iterator SectCollPhiColl_iterator;
typedef vector<DTSectCollThSegm>  SectCollThetaColl;
typedef SectCollThetaColl::const_iterator SectCollThetaColl_iterator;

/*
  Channels are labeled following next schema:
    ---------------------------------
    |   6   |   7   |   8   |   9   |
    ---------------------------------
        |   3   |   4   |   5   |
        -------------------------
            |   1   |   2   |
            -----------------
                |   0   |
                ---------
*/

/* Cell's combination, following previous labeling, to obtain every possible  muon's path. Others cells combinations imply non straight paths */
// const int CHANNELS_PATH_ARRANGEMENTS[8][4] = {
//     {0, 1, 3, 6}, {0, 1, 3, 7}, {0, 1, 4, 7}, {0, 1, 4, 8},
//     {0, 2, 4, 7}, {0, 2, 4, 8}, {0, 2, 5, 8}, {0, 2, 5, 9}
// };

/* For each of the previous cell's combinations, this array stores the associated cell's displacement, relative to lower layer cell, measured in semi-cell length units */

// const int CELL_HORIZONTAL_LAYOUTS[8][4] = {
//     {0, -1, -2, -3}, {0, -1, -2, -1}, {0, -1, 0, -1}, {0, -1, 0, 1},
//     {0,  1,  0, -1}, {0,  1,  0,  1}, {0,  1, 2,  1}, {0,  1, 2, 3}
// };


DTTrigPhase2Prod::DTTrigPhase2Prod(const ParameterSet& pset):
  chiSquareThreshold(50),
  bxTolerance(30),
  minQuality(LOWQGHOST)
{
    
    produces<L1MuDTChambContainer>();
    produces<L1MuDTChambPhContainer>();
    produces<L1MuDTChambThContainer>();
    produces<L1Phase2MuDTPhContainer>();
    
    debug = pset.getUntrackedParameter<bool>("debug");
    pinta = pset.getUntrackedParameter<bool>("pinta");
    tanPhiTh = pset.getUntrackedParameter<double>("tanPhiTh");
    dT0_correlate_TP = pset.getUntrackedParameter<double>("dT0_correlate_TP");
    min_dT0_match_segment = pset.getUntrackedParameter<double>("min_dT0_match_segment");
    min_phinhits_match_segment = pset.getUntrackedParameter<int>("min_phinhits_match_segment");
    minx_match_2digis = pset.getUntrackedParameter<double>("minx_match_2digis");

    do_correlation = pset.getUntrackedParameter<bool>("do_correlation");
    p2_df = pset.getUntrackedParameter<int>("p2_df");
    filter_primos = pset.getUntrackedParameter<bool>("filter_primos");
    
    txt_ttrig_bc0 = pset.getUntrackedParameter<bool>("apply_txt_ttrig_bc0");
    
    dt4DSegmentsToken = consumes<DTRecSegment4DCollection>(pset.getParameter < edm::InputTag > ("dt4DSegments"));
    dtDigisToken = consumes< DTDigiCollection >(pset.getParameter<edm::InputTag>("digiTag"));

    rpcRecHitsLabel = consumes<RPCRecHitCollection>(pset.getUntrackedParameter < edm::InputTag > ("rpcRecHits"));
  
    
  // Choosing grouping scheme:
  grcode = pset.getUntrackedParameter<Int_t>("grouping_code");
  
  if (grcode == 0) grouping_obj = new InitialGrouping(pset);
  if (grcode == 2) grouping_obj = new PseudoBayesGrouping(pset);
  else {
    if (debug) cout << "DTp2::constructor: Non-valid grouping code. Choosing InitialGrouping by default." << endl;
    grouping_obj = new InitialGrouping(pset);
  }
  
  mpathanalyzer   = new MuonPathAnalyzerPerSL(pset);
  mpathfilter     = new MuonPathFilter(pset);
  mpathassociator = new MuonPathAssociator(pset);
  
    if(pinta){
	std::cout<<"BOOKING HISTOS"<<std::endl;

	theFileOut = new TFile("dt_phase2.root", "RECREATE");

	Nsegments = new TH1F("Nsegments","Nsegments",21,-0.5,20.5);
	NmetaPrimitives = new TH1F("NmetaPrimitives","NmetaPrimitives",201,-0.5,200.5);
	NfilteredMetaPrimitives = new TH1F("NfilteredMetaPrimitives","NfilteredMetaPrimitives",201,-0.5,200.5);
	NcorrelatedMetaPrimitives = new TH1F("NcorrelatedMetaPrimitives","NcorrelatedMetaPrimitives",201,-0.5,200.5);
	Ngroups = new TH1F("Ngroups","Ngroups",201,-0.5,200.5);
	Nquality = new TH1F("Nquality","Nquality",9,0.5,9.5);
	Nquality_matched = new TH1F("Nquality_matched","Nquality_matched",9,0.5,9.5);
	Nsegosl = new TH1F("Nsegosl","Nsegosl",100,-10,10);
	Nsegosl31 = new TH1F("Nsegosl31","Nsegosl31",100,-10,10);
	Nmd = new TH1F("Nmd","Nmd",11,-0.5,10.5);
	Nmd31 = new TH1F("Nmd31","Nmd31",11,-0.5,10.5);
	Nhits_segment_tp = new TH2F("Nhits_segment_tp","Nhits_segment_tp",10,-0.5,9.5,10,-0.5,9.5);
	
	char name [128];

	for(int wh=-2;wh<=2;wh++){
	    int iwh=wh+2;
	    auto swh = std::to_string(wh);
	    for(int st=1;st<=4;st++){
		int ist=st-1;
		auto sst = std::to_string(st);
		for(int se=1;se<=14;se++){
		    if(se>=13&&st!=4)continue;
		    int ise=se-1;
		    auto sse = std::to_string(se);
		    for(int qu=1;qu<=9;qu++){
			int iqu=qu-1;
			auto squ = std::to_string(qu);
			  
			std::string nameSL = "Wh"+swh+"_St"+sst+"_Se"+sse+"_Qu"+squ;
			  
			//TIME
			sprintf(name,"TIMEPhase2_%s",nameSL.c_str());
			TIMEPhase2[iwh][ist][ise][iqu] = new TH1F(name,name,100,-0.5,89075.5);
			      
			//T0
			sprintf(name,"TOPhase2_%s",nameSL.c_str());
			T0Phase2[iwh][ist][ise][iqu] = new TH1F(name,name,100,-0.5,89075.5);
    
			//2D
			sprintf(name,"segment_vs_jm_x_%s",nameSL.c_str());
			segment_vs_jm_x[iwh][ist][ise][iqu] = new TH2F(name,name,250,-250,250,250,-250,250);
			  
			sprintf(name,"segment_vs_jm_x_gauss_%s",nameSL.c_str());
			segment_vs_jm_x_gauss[iwh][ist][ise][iqu] = new TH1F(name,name,300,-0.04,0.04);

			sprintf(name,"segment_vs_jm_tanPhi_%s",nameSL.c_str());
			segment_vs_jm_tanPhi[iwh][ist][ise][iqu] = new TH2F(name,name,100,-1.5,1.5,100,-1.5,1.5);

			sprintf(name,"segment_vs_jm_tanPhi_gauss_%s",nameSL.c_str());
			segment_vs_jm_tanPhi_gauss[iwh][ist][ise][iqu] = new TH1F(name,name,300,-0.5,0.5); //for single ones resolution
			//segment_vs_jm_tanPhi_gauss[iwh][ist][ise][iqu] = new TH1F(name,name,300,-0.01,0.01); //for correlated
			  
			sprintf(name,"segment_vs_jm_T0_%s",nameSL.c_str());
			segment_vs_jm_T0[iwh][ist][ise][iqu] = new TH2F(name,name,100,0,90000,100,0,90000);

			sprintf(name,"segment_vs_jm_T0_gauss_%s",nameSL.c_str());
			segment_vs_jm_T0_gauss[iwh][ist][ise][iqu] = new TH1F(name,name,300,-100,100);

			sprintf(name,"segment_vs_jm_T0_gauss_all_%s",nameSL.c_str());
			segment_vs_jm_T0_gauss_all[iwh][ist][ise][iqu] = new TH1F(name,name,300,-100,100);
			
			sprintf(name,"observed_tanPsi_%s",nameSL.c_str());
			observed_tanPsi[iwh][ist][ise][iqu] = new TH1F(name,name,100,-1.5,1.5);

			sprintf(name,"all_observed_tanPsi_%s",nameSL.c_str());
			all_observed_tanPsi[iwh][ist][ise][iqu] = new TH1F(name,name,100,-1.5,1.5);
			
			sprintf(name,"observed_x_%s",nameSL.c_str());
			observed_x[iwh][ist][ise][iqu] = new TH1F(name,name,250,-250,250);

			sprintf(name,"all_observed_x_%s",nameSL.c_str());
			all_observed_x[iwh][ist][ise][iqu] = new TH1F(name,name,250,-250,250);

			sprintf(name,"observed_t0_%s",nameSL.c_str());
			observed_t0[iwh][ist][ise][iqu] = new TH1F(name,name,100,0,90000);

			sprintf(name,"all_observed_t0_%s",nameSL.c_str());
			all_observed_t0[iwh][ist][ise][iqu] = new TH1F(name,name,100,-100,100);
			
			sprintf(name,"chi2_%s",nameSL.c_str());
			chi2[iwh][ist][ise][iqu] = new TH1F(name,name,100,0.,0.02);

			sprintf(name,"TPphi_%s",nameSL.c_str());
			TPphi[iwh][ist][ise][iqu] = new TH1F(name,name,250,-1.5,1.5);

			sprintf(name,"TPphiB_%s",nameSL.c_str());
			TPphiB[iwh][ist][ise][iqu] = new TH1F(name,name,250,-1.5,1.5);

			sprintf(name,"MP_x_back_%s",nameSL.c_str());
			MP_x_back[iwh][ist][ise][iqu] = new TH2F(name,name,100,-250,250,100,-250,250);

			sprintf(name,"MP_psi_back_%s",nameSL.c_str());
			MP_psi_back[iwh][ist][ise][iqu] = new TH2F(name,name,100,-1.5,1.5,100,-1.5,1.5);

			
		    }
		    std::string nameSL = "Wh"+swh+"_St"+sst+"_Se"+sse;
		    sprintf(name,"expected_tanPsi_%s",nameSL.c_str());
		    expected_tanPsi[iwh][ist][ise] = new TH1F(name,name,100,-1.5,1.5);

		    sprintf(name,"expected_x_%s",nameSL.c_str());
		    expected_x[iwh][ist][ise] = new TH1F(name,name,250,-250,250);

		    sprintf(name,"expected_t0_%s",nameSL.c_str());
		    expected_t0[iwh][ist][ise] = new TH1F(name,name,100,0,90000);
		    
		}
	    }
	}
	
    }	  

    int rawId;
    
    //ttrig
    ttrig_filename = pset.getUntrackedParameter<std::string>("ttrig_filename");
    if(txt_ttrig_bc0){
	std::ifstream ifin(ttrig_filename.c_str());
	double ttrig;
	while (ifin.good()){
	    ifin >> rawId >> ttrig;
	    ttriginfo[rawId]=ttrig;
	}
    }
    
    //z
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

DTTrigPhase2Prod::~DTTrigPhase2Prod(){

    //delete inMuonPath;
    //delete outValidMuonPath;
    
    if(debug) std::cout<<"DTp2: calling destructor"<<std::endl;

    if(pinta){
	if(debug) std::cout<<"DTp2: writing histograms and files"<<std::endl;
	
	theFileOut->cd();

	for(int wh=-2;wh<=2;wh++){
	    int iwh=wh+2;
	    for(int st=1;st<=4;st++){
		int ist=st-1;
		for(int se=1;se<=14;se++){
		    int ise=se-1;
		    if(se>=13&&st!=4)continue;
		    for(int qu=1;qu<=9;qu++){
			int iqu=qu-1;

			//digi TIME
			TIMEPhase2[iwh][ist][ise][iqu]->Write();
			      
			//digiT0
			T0Phase2[iwh][ist][ise][iqu]->Write();
    
			//2D
			segment_vs_jm_x[iwh][ist][ise][iqu]->Write();
			  
			segment_vs_jm_x_gauss[iwh][ist][ise][iqu]->Write();

			segment_vs_jm_tanPhi[iwh][ist][ise][iqu]->Write();

			segment_vs_jm_tanPhi_gauss[iwh][ist][ise][iqu]->Write();
			  
			segment_vs_jm_T0[iwh][ist][ise][iqu]->Write();

			segment_vs_jm_T0_gauss[iwh][ist][ise][iqu]->Write();
			segment_vs_jm_T0_gauss_all[iwh][ist][ise][iqu]->Write();

			observed_tanPsi[iwh][ist][ise][iqu]->Write();
			
			all_observed_tanPsi[iwh][ist][ise][iqu]->Write();

			observed_x[iwh][ist][ise][iqu]->Write();
			
			all_observed_x[iwh][ist][ise][iqu]->Write();

			observed_t0[iwh][ist][ise][iqu]->Write();
			
			all_observed_t0[iwh][ist][ise][iqu]->Write();

			chi2[iwh][ist][ise][iqu]->Write();
			
			TPphi[iwh][ist][ise][iqu]->Write();
			TPphiB[iwh][ist][ise][iqu]->Write();

			MP_x_back[iwh][ist][ise][iqu]->Write();
			MP_psi_back[iwh][ist][ise][iqu]->Write();

		    }
		    expected_tanPsi[iwh][ist][ise]->Write();
		    expected_x[iwh][ist][ise]->Write();
		    expected_t0[iwh][ist][ise]->Write();
		}
	    }
	}

	Nsegments->Write();
	NmetaPrimitives->Write();
	NfilteredMetaPrimitives->Write();
	NcorrelatedMetaPrimitives->Write();
	Ngroups->Write();
	Nquality->Write();
	Nquality_matched->Write();
	Nsegosl->Write();
	Nsegosl31->Write();
	Nmd->Write();
	Nmd31->Write();
	Nhits_segment_tp->Write();
	
	theFileOut->Write();
	theFileOut->Close();

	delete Nsegments;
	delete NmetaPrimitives;
	delete NfilteredMetaPrimitives;
	delete NcorrelatedMetaPrimitives;
	delete Ngroups;
	delete Nquality;
	delete Nquality_matched;
	delete Nsegosl;
	delete Nsegosl31;
	delete Nmd;
	delete Nmd31;
	delete Nhits_segment_tp;
    
	for(int wh=-2;wh<=2;wh++){
	    int iwh=wh+2;
	    for(int st=1;st<=4;st++){
		int ist=st-1;
		for(int se=1;se<=14;se++){
		    int ise=se-1;
		    if(se>=13&&st!=4)continue;
		    for(int qu=1;qu<=9;qu++){
			int iqu=qu-1;
	    
			//digi TIME
			delete TIMEPhase2[iwh][ist][ise][iqu];
	    
			//digiT0
			delete T0Phase2[iwh][ist][ise][iqu];
	    
			//2D
			delete segment_vs_jm_x[iwh][ist][ise][iqu];    
			delete segment_vs_jm_x_gauss[iwh][ist][ise][iqu];
			delete segment_vs_jm_tanPhi[iwh][ist][ise][iqu];
			delete segment_vs_jm_tanPhi_gauss[iwh][ist][ise][iqu];
			delete segment_vs_jm_T0[iwh][ist][ise][iqu];
			delete segment_vs_jm_T0_gauss[iwh][ist][ise][iqu];
			delete segment_vs_jm_T0_gauss_all[iwh][ist][ise][iqu];
			delete observed_tanPsi[iwh][ist][ise][iqu];
			delete all_observed_tanPsi[iwh][ist][ise][iqu];
			delete observed_x[iwh][ist][ise][iqu];
			delete all_observed_x[iwh][ist][ise][iqu];
			delete observed_t0[iwh][ist][ise][iqu];
			delete all_observed_t0[iwh][ist][ise][iqu];
			
			delete chi2[iwh][ist][ise][iqu];
			delete TPphi[iwh][ist][ise][iqu];
			delete TPphiB[iwh][ist][ise][iqu];

			delete MP_x_back[iwh][ist][ise][iqu];
			delete MP_psi_back[iwh][ist][ise][iqu];
			
		    }
		    delete expected_tanPsi[iwh][ist][ise];
		    delete expected_x[iwh][ist][ise];
		    delete expected_t0[iwh][ist][ise];
		}
	    }
	}
    
	delete theFileOut; 
    }
  
  delete grouping_obj; // Grouping destructor
  delete mpathanalyzer;
  delete mpathfilter;
  delete mpathassociator;

}


void DTTrigPhase2Prod::beginRun(edm::Run const& iRun, const edm::EventSetup& iEventSetup) {
  if(debug) cout << "DTTrigPhase2Prod::beginRun " << iRun.id().run() << endl;
  if(debug) cout << "DTTrigPhase2Prod::beginRun: getting DT geometry" << endl;
    
  if(debug) std::cout<<"getting DT geometry"<<std::endl;
  iEventSetup.get<MuonGeometryRecord>().get(dtGeo);//1103
  
  ESHandle< DTConfigManager > dtConfig ;
  iEventSetup.get< DTConfigManagerRcd >().get( dtConfig );

  grouping_obj->initialise(iEventSetup); // Grouping object initialisation
  mpathanalyzer->initialise(iEventSetup); // pathanalyzer object initialisation
  mpathfilter->initialise(iEventSetup);
  mpathassociator->initialise(iEventSetup);
  
  //trigGeomUtils = new DTTrigGeomUtils(dtGeo);

  //filling up zcn
  for (int ist=1; ist<=4; ++ist) {
    const DTChamberId chId(-2,ist,4);
    const DTChamber *chamb = dtGeo->chamber(chId);
    const DTSuperLayer *sl1 = chamb->superLayer(DTSuperLayerId(chId,1));
    const DTSuperLayer *sl3 = chamb->superLayer(DTSuperLayerId(chId,3));
    zcn[ist-1] = .5*(chamb->surface().toLocal(sl1->position()).z() + chamb->surface().toLocal(sl3->position()).z());
  }
  
  const DTChamber* chamb   = dtGeo->chamber(DTChamberId(-2,4,13));
  const DTChamber* scchamb = dtGeo->chamber(DTChamberId(-2,4,4));
    xCenter[0] = scchamb->toLocal(chamb->position()).x()*.5;
    chamb   = dtGeo->chamber(DTChamberId(-2,4,14));
    scchamb = dtGeo->chamber(DTChamberId(-2,4,10));
    xCenter[1] = scchamb->toLocal(chamb->position()).x()*.5;
}


void DTTrigPhase2Prod::produce(Event & iEvent, const EventSetup& iEventSetup){
  edm::Handle<DTDigiCollection> dtdigis;
  iEvent.getByToken(dtDigisToken, dtdigis);
    
  edm::Handle<DTRecSegment4DCollection> all4DSegments;
  if(pinta){
    iEvent.getByToken(dt4DSegmentsToken, all4DSegments);
    if(debug) std::cout<<"DTp2: I got the segments"<<std::endl;
  }
  //int bx32 = iEvent.eventAuxiliary().bunchCrossing()*32;
  int bx25 = iEvent.eventAuxiliary().bunchCrossing()*25;
  
  if(debug) std::cout <<"\t Getting the RPC RecHits"<<std::endl;
  Handle<RPCRecHitCollection> rpcHits;
  iEvent.getByToken(rpcRecHitsLabel,rpcHits);
  
  
  //Santi's code
  // GROUPING BEGIN
  DTDigiMap digiMap;
  DTDigiCollection::DigiRangeIterator detUnitIt;
  for (detUnitIt=dtdigis->begin(); detUnitIt!=dtdigis->end(); ++detUnitIt) {
    const DTLayerId& layId               = (*detUnitIt).first;
    const DTChamberId chambId            = layId.superlayerId().chamberId();
    const DTDigiCollection::Range& range = (*detUnitIt).second; 
    digiMap[chambId].put(range,layId);
  }
  
  // generate a list muon paths for each event!!!
  std::vector<MuonPath*> muonpaths;
  for (std::vector<const DTChamber*>::const_iterator ich = dtGeo->chambers().begin(); ich != dtGeo->chambers().end(); ich++) {
    const DTChamber* chamb  = (*ich);
    DTChamberId chid        = chamb->id();
    DTDigiMap_iterator dmit = digiMap.find(chid);
    
    //     if (dmit !=digiMap.end()) buildMuonPathCandidates((*dmit).second, &muonpaths);              // Old grouping implementation
    if (dmit !=digiMap.end()) grouping_obj->run(iEvent, iEventSetup, (*dmit).second, &muonpaths);  // New grouping implementation
  }
  
  if (debug && grcode == 2){
    for (std::vector<MuonPath*>::iterator itmPaths = muonpaths.begin(); itmPaths != muonpaths.end(); itmPaths++){
      std::cout << "---------------------------------------------------------------" << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(0)->getLayerId() << " , " << (*itmPaths)->getPrimitive(0)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(0)->getChannelId() << " , " << (*itmPaths)->getPrimitive(0)->getLaterality() << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(1)->getLayerId() << " , " << (*itmPaths)->getPrimitive(1)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(1)->getChannelId() << " , " << (*itmPaths)->getPrimitive(1)->getLaterality()<< std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(2)->getLayerId() << " , " << (*itmPaths)->getPrimitive(2)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(2)->getChannelId() << " , " << (*itmPaths)->getPrimitive(2)->getLaterality() << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(3)->getLayerId() << " , " << (*itmPaths)->getPrimitive(3)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(3)->getChannelId() << " , " << (*itmPaths)->getPrimitive(3)->getLaterality() << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(4)->getLayerId() << " , " << (*itmPaths)->getPrimitive(4)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(4)->getChannelId() << " , " << (*itmPaths)->getPrimitive(4)->getLaterality() << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(5)->getLayerId() << " , " << (*itmPaths)->getPrimitive(5)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(5)->getChannelId() << " , " << (*itmPaths)->getPrimitive(5)->getLaterality() << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(6)->getLayerId() << " , " << (*itmPaths)->getPrimitive(6)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(6)->getChannelId() << " , " << (*itmPaths)->getPrimitive(6)->getLaterality()  << std::endl;
      std::cout << "DTTrigPhase2Prod::beginRun DTTPath " << (*itmPaths)->getPrimitive(7)->getLayerId() << " , " << (*itmPaths)->getPrimitive(7)->getSuperLayerId() << " , " << (*itmPaths)->getPrimitive(7)->getChannelId() << " , " << (*itmPaths)->getPrimitive(7)->getLaterality() << std::endl;
      std::cout << "---------------------------------------------------------------" << std::endl;
    }
  }
    
  digiMap.clear();
  // GROUPING ENDS

  if (debug) cout << "MUON PATHS found: " << muonpaths.size() <<" in event"<<iEvent.id().event()<<endl;
  if(pinta) Ngroups->Fill(muonpaths.size());
  
//  std::vector<MuonPath*> outmpath;
//  mpathanalyzer->run(iEvent, iEventSetup, muonpaths, outmpath);
  std::vector<metaPrimitive> metaPrimitives;
  mpathanalyzer->run(iEvent, iEventSetup, muonpaths, metaPrimitives);

  /*
    //filtro por groupos de TDC times en las mismas celdas... corrobarar si sucede... esta implementacion no existe en software pero existe en firmware
    // loop over vector of muonpahts produced by grouping
    int iGroup=1;
    std::vector<metaPrimitive> metaPrimitives;
    for(auto muonpath = muonpaths.begin();muonpath!=muonpaths.end();++muonpath){
	DTPrimitive testPrim0((*muonpath)->getPrimitive(0));
	DTPrimitive testPrim1((*muonpath)->getPrimitive(1));
	DTPrimitive testPrim2((*muonpath)->getPrimitive(2));
	DTPrimitive testPrim3((*muonpath)->getPrimitive(3));

	if(debug){
	    std::cout<<"test_grouping:group # "<<iGroup<<std::endl;
	    std::cout<<"test_grouping:\t 0 "<<testPrim0.getChannelId()<<" "<<testPrim0.getTDCTime()<<std::endl;
	    std::cout<<"test_grouping:\t 1 "<<testPrim1.getChannelId()<<" "<<testPrim1.getTDCTime()<<std::endl;
	    std::cout<<"test_grouping:\t 2 "<<testPrim2.getChannelId()<<" "<<testPrim2.getTDCTime()<<std::endl;
	    std::cout<<"test_grouping:\t 3 "<<testPrim3.getChannelId()<<" "<<testPrim3.getTDCTime()<<std::endl;
	}

	int selected_Id=0;
	if(testPrim0.getTDCTime()!=-1) selected_Id= testPrim0.getCameraId();
	else if(testPrim1.getTDCTime()!=-1) selected_Id= testPrim1.getCameraId(); 
	else if(testPrim2.getTDCTime()!=-1) selected_Id= testPrim2.getCameraId(); 
	else if(testPrim3.getTDCTime()!=-1) selected_Id= testPrim3.getCameraId(); 
 

	DTLayerId thisLId(selected_Id);
	if(debug) std::cout<<"Building up MuonPathSLId from rawId in the Primitive"<<std::endl;
	DTSuperLayerId MuonPathSLId(thisLId.wheel(),thisLId.station(),thisLId.sector(),thisLId.superLayer());
	if(debug) std::cout<<"The MuonPathSLId is"<<MuonPathSLId<<std::endl;

	if(chosen_sl==4){
	    if(debug) std::cout<<"Generating Primitives with Both Superlayers 1 and 3"<<std::endl;
	    analyze(*muonpath,metaPrimitives,MuonPathSLId);
	}else if(thisLId.superLayer()==chosen_sl){
	    if(debug) std::cout<<"Generating Primitives only with SL "<<chosen_sl<<std::endl;
	    analyze(*muonpath,metaPrimitives,MuonPathSLId);
	}
	iGroup++;
    }
  */
    if(debug) std::cout<<"filling NmetaPrimtives"<<std::endl;
    
    if(pinta) NmetaPrimitives->Fill(metaPrimitives.size());

    if(debug) std::cout<<"deleting muonpaths"<<std::endl;    
    for (unsigned int i=0; i<muonpaths.size(); i++){
      delete muonpaths[i];
    }
    muonpaths.clear();
    
    //FILTER SECTIONS:
    //filtro de duplicados puro popdr'ia ir ac'a mpredundantfilter.cpp primos?
    //filtro en |tanPhi|<~1.?
    //filtro de calidad por chi2 qualityenhancefilter.cpp mirar el metodo filter

    if(debug) std::cout<<"declaring new vector for filtered"<<std::endl;    

    std::vector<metaPrimitive> filteredMetaPrimitives;

    if(filter_primos){

	if(debug) std::cout<<"filtering: starting primos filtering"<<std::endl;    
    
	int primo_index=0;
	bool oneof4=false;
	//    for (auto metaPrimitiveIt = metaPrimitives.begin(); metaPrimitiveIt != metaPrimitives.end(); ++metaPrimitiveIt){
	if(metaPrimitives.size()==1){
	    if(debug){
		std::cout<<"filtering:";
		printmP(metaPrimitives[0]);
		std::cout<<" \t is:"<<0<<" "<<primo_index<<" "<<" "<<oneof4<<std::endl;
	    }
	    if(fabs(metaPrimitives[0].tanPhi)<tanPhiTh){
		filteredMetaPrimitives.push_back(metaPrimitives[0]);
		if(debug)std::cout<<"filtering: kept1 i="<<0<<std::endl;
	    }
	}
	else for(int i=1; i<int(metaPrimitives.size()); i++){ 
		if(fabs(metaPrimitives[i].tanPhi)>tanPhiTh) continue;
		if(rango(metaPrimitives[i])==4)oneof4=true;
		if(debug){
		    std::cout<<"filtering:";
		    printmP(metaPrimitives[i]);
		    std::cout<<" \t is:"<<i<<" "<<primo_index<<" "<<" "<<oneof4<<std::endl;
		}
		if(arePrimos(metaPrimitives[i],metaPrimitives[i-1])!=0  and arePrimos(metaPrimitives[i],metaPrimitives[i-primo_index-1])!=0){
		    primo_index++;
		}else{
		    if(primo_index==0){
			filteredMetaPrimitives.push_back(metaPrimitives[i]);
			if(debug)std::cout<<"filtering: kept2 i="<<i<<std::endl;
		    }else{
			if(oneof4){
			    float minchi2=99999;
			    int selected_i=0;
			    for(int j=i-1;j>=i-primo_index-1;j--){
				if(rango(metaPrimitives[j])!=4) continue;
				if(minchi2>metaPrimitives[j].chi2){
				    minchi2=metaPrimitives[j].chi2;
				    selected_i=j;
				}
			    }
			    filteredMetaPrimitives.push_back(metaPrimitives[selected_i]);
			    if(debug)std::cout<<"filtering: kept4 i="<<selected_i<<std::endl;
			}else{
			    for(int j=i-1;j>=i-primo_index-1;j--){
				filteredMetaPrimitives.push_back(metaPrimitives[j]);
				if(debug)std::cout<<"filtering: kept3 i="<<j<<std::endl;
			    }
			}
		    }
		    primo_index=0;
		    oneof4=false;
		}
	    }
    }else{
	for (size_t i=0; i<metaPrimitives.size(); i++){ 
	    if(fabs(metaPrimitives[i].tanPhi)>tanPhiTh) continue;
	    filteredMetaPrimitives.push_back(metaPrimitives[i]); 
	}
    }
    
    metaPrimitives.clear();
    metaPrimitives.erase(metaPrimitives.begin(),metaPrimitives.end());

    if(pinta) NfilteredMetaPrimitives->Fill(filteredMetaPrimitives.size());
    
    
    if(debug) std::cout<<"filteredMetaPrimitives: starting correlations"<<std::endl;    
    if(!do_correlation){
	if(debug) std::cout<<"DTp2 in event:"<<iEvent.id().event()<<" we found "<<filteredMetaPrimitives.size()<<" filteredMetaPrimitives (superlayer)"<<std::endl;
	vector<L1MuDTChambPhDigi> outPhi;
	vector<L1MuDTChambDigi> outP2;
	vector<L1Phase2MuDTPhDigi> outP2Ph;

	for (auto metaPrimitiveIt = filteredMetaPrimitives.begin(); metaPrimitiveIt != filteredMetaPrimitives.end(); ++metaPrimitiveIt){
	    DTSuperLayerId slId((*metaPrimitiveIt).rawId);
	    if(debug) std::cout<<"looping in final vector: SuperLayerId"<<slId<<" x="<<(*metaPrimitiveIt).x<<" quality="<<(*metaPrimitiveIt).quality<<std::endl;
	    
	    int sectorTP=slId.sector();
	    if(sectorTP==13) sectorTP=4;
	    if(sectorTP==14) sectorTP=10;
	    sectorTP=sectorTP-1;
	  
	    if(p2_df==0){
		outPhi.push_back(L1MuDTChambPhDigi((*metaPrimitiveIt).t0,
						   slId.wheel(),
						   sectorTP,
						   slId.station(),
						   (int)round((*metaPrimitiveIt).phi*65536./0.8),
						   (int)round((*metaPrimitiveIt).phiB*2048./1.4),
						   (*metaPrimitiveIt).quality,
						   1,
						   0
						   ));
	    }else if(p2_df==1){
		if(debug)std::cout<<"pushing back phase-2 dataformat agreement with Oscar for comparison with slice test"<<std::endl;
		outP2.push_back(L1MuDTChambDigi((int)round((*metaPrimitiveIt).t0/25.),   // ubx (m_bx) //bx en la orbita
						slId.wheel(),   // uwh (m_wheel)     
						sectorTP,   // usc (m_sector)    
						slId.station(),   // ust (m_station)
						(int)round((*metaPrimitiveIt).x*1000),   // uphi (_phiAngle)
						(int)round((*metaPrimitiveIt).tanPhi*4096),   // uphib (m_phiBending)
						0,   // uz (m_zCoordinate)
						0,   // uzsl (m_zSlope)
						(*metaPrimitiveIt).quality,  // uqua (m_qualityCode)
						0,  // uind (m_segmentIndex)
						(int)round((*metaPrimitiveIt).t0),  // ut0 (m_t0Segment)
						(int)round((*metaPrimitiveIt).chi2),  // uchi2 (m_chi2Segment)
						-10    // urpc (m_rpcFlag)
						));
	    }else if(p2_df==3){
		if(debug)std::cout<<"pushing back phase-2 dataformat carlo-federica dataformat"<<std::endl;
		outP2Ph.push_back(L1Phase2MuDTPhDigi((int)round((*metaPrimitiveIt).t0/25.),   // ubx (m_bx) //bx en la orbita
						     slId.wheel(),   // uwh (m_wheel)     // FIXME: It is not clear who provides this?
						     sectorTP,   // usc (m_sector)    // FIXME: It is not clear who provides this?
						     slId.station(),   // ust (m_station)
						     (int)round((*metaPrimitiveIt).phi*65536./0.8), // uphi (_phiAngle)
						     (int)round((*metaPrimitiveIt).phiB*2048./1.4), // uphib (m_phiBending)
						     (*metaPrimitiveIt).quality,  // uqua (m_qualityCode)
						     0,  // uind (m_segmentIndex)
						     (int)round((*metaPrimitiveIt).t0),  // ut0 (m_t0Segment)
						     (int)round((*metaPrimitiveIt).chi2),  // uchi2 (m_chi2Segment)
						     -10    // urpc (m_rpcFlag)
						     ));
		
	    }  
	}		

	if(p2_df==0){
	    std::unique_ptr<L1MuDTChambPhContainer> resultPhi (new L1MuDTChambPhContainer);
	    resultPhi->setContainer(outPhi); iEvent.put(std::move(resultPhi));
	    outPhi.clear();
	    outPhi.erase(outPhi.begin(),outPhi.end());
	}else if(p2_df==1){
	    std::unique_ptr<L1MuDTChambContainer> resultP2 (new L1MuDTChambContainer);
	    resultP2->setContainer(outP2); iEvent.put(std::move(resultP2));
	    outP2.clear();
	    outP2.erase(outP2.begin(),outP2.end());
	}else if(p2_df==2){
	    std::unique_ptr<L1Phase2MuDTPhContainer> resultP2Ph (new L1Phase2MuDTPhContainer);
	    resultP2Ph->setContainer(outP2Ph); iEvent.put(std::move(resultP2Ph));
	    outP2Ph.clear();
	    outP2Ph.erase(outP2Ph.begin(),outP2Ph.end());
	}
	    
    }else{
	//Silvia's code for correlationg filteredMetaPrimitives
	
	if(debug) std::cout<<"starting correlation"<<std::endl;
	
	std::vector<metaPrimitive> correlatedMetaPrimitives;
	
	for(int wh=-2;wh<=2;wh++){
	    for(int st=1;st<=4;st++){
		for(int se=1;se<=14;se++){
		    if(se>=13&&st!=4)continue;
		    
		    DTChamberId ChId(wh,st,se);
		    DTSuperLayerId sl1Id(wh,st,se,1);
		    DTSuperLayerId sl3Id(wh,st,se,3);
	      
		    //filterSL1
		    std::vector<metaPrimitive> SL1metaPrimitives;
		    for(auto metaprimitiveIt = filteredMetaPrimitives.begin();metaprimitiveIt!=filteredMetaPrimitives.end();++metaprimitiveIt)
			if(metaprimitiveIt->rawId==sl1Id.rawId())
			    SL1metaPrimitives.push_back(*metaprimitiveIt);
	      
		    //filterSL3
		    std::vector<metaPrimitive> SL3metaPrimitives;
		    for(auto metaprimitiveIt = filteredMetaPrimitives.begin();metaprimitiveIt!=filteredMetaPrimitives.end();++metaprimitiveIt)
			if(metaprimitiveIt->rawId==sl3Id.rawId())
			    SL3metaPrimitives.push_back(*metaprimitiveIt);
		    
		    if(SL1metaPrimitives.size()==0 and SL3metaPrimitives.size()==0) continue;
		    
		    if(debug) std::cout<<"correlating "<<SL1metaPrimitives.size()<<" metaPrim in SL1 and "<<SL3metaPrimitives.size()<<" in SL3 for "<<sl3Id<<std::endl;

		    bool at_least_one_correlation=false;

		    //SL1-SL3

		    for (auto SL1metaPrimitive = SL1metaPrimitives.begin(); SL1metaPrimitive != SL1metaPrimitives.end(); ++SL1metaPrimitive){
			for (auto SL3metaPrimitive = SL3metaPrimitives.begin(); SL3metaPrimitive != SL3metaPrimitives.end(); ++SL3metaPrimitive){
			    if(fabs(SL1metaPrimitive->t0-SL3metaPrimitive->t0) < dT0_correlate_TP){//time match
				float PosSL1=SL1metaPrimitive->x;
				float PosSL3=SL3metaPrimitive->x;
				float NewSlope=(PosSL1-PosSL3)/23.5;     
				float MeanT0=(SL1metaPrimitive->t0+SL3metaPrimitive->t0)/2;
				float MeanPos=(PosSL3+PosSL1)/2;
				float newChi2=(SL1metaPrimitive->chi2+SL3metaPrimitive->chi2)*0.5;//to be recalculated
				int quality = 0;
				if(SL3metaPrimitive->quality <= 2 and SL1metaPrimitive->quality <=2) quality=6;

				if((SL3metaPrimitive->quality >= 3 && SL1metaPrimitive->quality <=2)
				   or (SL1metaPrimitive->quality >= 3 && SL3metaPrimitive->quality <=2) ) quality=8;

				if(SL3metaPrimitive->quality >= 3 && SL1metaPrimitive->quality >=3) quality=9;
			  
				GlobalPoint jm_x_cmssw_global = dtGeo->chamber(ChId)->toGlobal(LocalPoint(MeanPos,0.,0.));//jm_x is already extrapolated to the middle of the SL
				int thisec = ChId.sector();
				if(se==13) thisec = 4;
				if(se==14) thisec = 10;
				float phi= jm_x_cmssw_global.phi()-0.5235988*(thisec-1);
				float psi=TMath::ATan(NewSlope);
				float phiB=hasPosRF(ChId.wheel(),ChId.sector()) ? psi-phi :-psi-phi ;
			
				correlatedMetaPrimitives.push_back(metaPrimitive({ChId.rawId(),MeanT0,MeanPos,NewSlope,phi,phiB,newChi2,quality,
						SL1metaPrimitive->wi1,SL1metaPrimitive->tdc1,
						SL1metaPrimitive->wi2,SL1metaPrimitive->tdc2,
						SL1metaPrimitive->wi3,SL1metaPrimitive->tdc3,
						SL1metaPrimitive->wi4,SL1metaPrimitive->tdc4,
						SL3metaPrimitive->wi1,SL3metaPrimitive->tdc1,
						SL3metaPrimitive->wi2,SL3metaPrimitive->tdc2,
						SL3metaPrimitive->wi3,SL3metaPrimitive->tdc3,
						SL3metaPrimitive->wi4,SL3metaPrimitive->tdc4
						}));
				at_least_one_correlation=true;
			    }
			}
			
			if(at_least_one_correlation==false){//no correlation was found, trying with pairs of two digis in the other SL
			    
			    int matched_digis=0;
			    float minx=minx_match_2digis;
			    int best_tdc=-1;
			    int next_tdc=-1;
			    int best_wire=-1;
			    int next_wire=-1;
			    int best_layer=-1;
			    int next_layer=-1;

			    for (auto dtLayerId_It=dtdigis->begin(); dtLayerId_It!=dtdigis->end(); ++dtLayerId_It){
				const DTLayerId dtLId = (*dtLayerId_It).first;
				DTSuperLayerId dtSLId(dtLId);
				if(dtSLId.rawId()!=sl3Id.rawId()) continue;
				float l_shift=0;
				if(dtLId.layer()==4)l_shift=1.95;
				if(dtLId.layer()==3)l_shift=0.65;
				if(dtLId.layer()==2)l_shift=-0.65;
				if(dtLId.layer()==1)l_shift=-1.95;
				float x_inSL3=SL1metaPrimitive->x-SL1metaPrimitive->tanPhi*(23.5+l_shift);
				for (auto digiIt = ((*dtLayerId_It).second).first;digiIt!=((*dtLayerId_It).second).second; ++digiIt){
				    DTWireId wireId(dtLId,(*digiIt).wire());
				    int x_wire = shiftinfo[wireId.rawId()]+((*digiIt).time()-SL1metaPrimitive->t0)*0.00543; 
				    int x_wire_left = shiftinfo[wireId.rawId()]-((*digiIt).time()-SL1metaPrimitive->t0)*0.00543; 
				    if(fabs(x_inSL3-x_wire)>fabs(x_inSL3-x_wire_left)) x_wire=x_wire_left; //choose the closest laterality
				    if(pinta) Nsegosl->Fill(x_inSL3-x_wire);
				    if(fabs(x_inSL3-x_wire)<minx){
					minx=fabs(x_inSL3-x_wire);
					next_wire=best_wire;
					next_tdc=best_tdc;
					next_layer=best_layer;
					
					best_wire=(*digiIt).wire();
					best_tdc=(*digiIt).time();
					best_layer=dtLId.layer();
					matched_digis++;
				    }
				}
				
			    }
			    if(pinta) Nmd->Fill(matched_digis);    
			    if(matched_digis>=2 and best_layer!=-1 and next_layer!=-1){
				int new_quality=7;
				if(SL1metaPrimitive->quality<=2) new_quality=5;

				int wi1=-1;int tdc1=-1;
				int wi2=-1;int tdc2=-1;
				int wi3=-1;int tdc3=-1;
				int wi4=-1;int tdc4=-1;
				
				if(next_layer==1) {wi1=next_wire; tdc1=next_tdc; }
				if(next_layer==2) {wi2=next_wire; tdc2=next_tdc; }
				if(next_layer==3) {wi3=next_wire; tdc3=next_tdc; }
				if(next_layer==4) {wi4=next_wire; tdc4=next_tdc; }

				if(best_layer==1) {wi1=best_wire; tdc1=best_tdc; }
				if(best_layer==2) {wi2=best_wire; tdc2=best_tdc; }
				if(best_layer==3) {wi3=best_wire; tdc3=best_tdc; }
				if(best_layer==4) {wi4=best_wire; tdc4=best_tdc; } 
				
				

				correlatedMetaPrimitives.push_back(metaPrimitive({ChId.rawId(),SL1metaPrimitive->t0,SL1metaPrimitive->x,SL1metaPrimitive->tanPhi,SL1metaPrimitive->phi,SL1metaPrimitive->phiB,SL1metaPrimitive->chi2,
						new_quality,
						SL1metaPrimitive->wi1,SL1metaPrimitive->tdc1,
						SL1metaPrimitive->wi2,SL1metaPrimitive->tdc2,
						SL1metaPrimitive->wi3,SL1metaPrimitive->tdc3,
						SL1metaPrimitive->wi4,SL1metaPrimitive->tdc4,
						wi1,tdc1,
						wi2,tdc2,
						wi3,tdc3,
						wi4,tdc4
						}));
				at_least_one_correlation=true;
			    }
			}
		    }

		    //finish SL1-SL3

		    //SL3-SL1
		    for (auto SL3metaPrimitive = SL3metaPrimitives.begin(); SL3metaPrimitive != SL3metaPrimitives.end(); ++SL3metaPrimitive){
			for (auto SL1metaPrimitive = SL1metaPrimitives.begin(); SL1metaPrimitive != SL1metaPrimitives.end(); ++SL1metaPrimitive){
			    if(fabs(SL1metaPrimitive->t0-SL3metaPrimitive->t0) < dT0_correlate_TP){//time match
				//this comb was already filled up in the previous loop now we just want to know if there was at least one match
				at_least_one_correlation=true;
			    }
			}
			
			if(at_least_one_correlation==false){//no correlation was found, trying with pairs of two digis in the other SL
				
			    int matched_digis=0;
			    float minx=minx_match_2digis;
			    int best_tdc=-1;
			    int next_tdc=-1;
			    int best_wire=-1;
			    int next_wire=-1;
			    int best_layer=-1;
			    int next_layer=-1;
				
			    for (auto dtLayerId_It=dtdigis->begin(); dtLayerId_It!=dtdigis->end(); ++dtLayerId_It){
				const DTLayerId dtLId = (*dtLayerId_It).first;
				DTSuperLayerId dtSLId(dtLId);
				if(dtSLId.rawId()!=sl1Id.rawId()) continue;
				float l_shift=0;
				if(dtLId.layer()==4)l_shift=1.95;
				if(dtLId.layer()==3)l_shift=0.65;
				if(dtLId.layer()==2)l_shift=-0.65;
				if(dtLId.layer()==1)l_shift=-1.95;
				float x_inSL1=SL3metaPrimitive->x+SL3metaPrimitive->tanPhi*(23.5-l_shift);
				for (auto digiIt = ((*dtLayerId_It).second).first;digiIt!=((*dtLayerId_It).second).second; ++digiIt){
				    DTWireId wireId(dtLId,(*digiIt).wire());
				    int x_wire = shiftinfo[wireId.rawId()]+((*digiIt).time()-SL3metaPrimitive->t0)*0.00543; 
				    int x_wire_left = shiftinfo[wireId.rawId()]-((*digiIt).time()-SL3metaPrimitive->t0)*0.00543; 
				    if(fabs(x_inSL1-x_wire)>fabs(x_inSL1-x_wire_left)) x_wire=x_wire_left; //choose the closest laterality
				    if(pinta) Nsegosl31->Fill(x_inSL1-x_wire);
				    if(fabs(x_inSL1-x_wire)<minx){
					minx=fabs(x_inSL1-x_wire);
					next_wire=best_wire;
					next_tdc=best_tdc;
					next_layer=best_layer;
					    
					best_wire=(*digiIt).wire();
					best_tdc=(*digiIt).time();
					best_layer=dtLId.layer();
					matched_digis++;
				    }
				}
				    
			    }
			    if(pinta) Nmd31->Fill(matched_digis);    
			    if(matched_digis>=2 and best_layer!=-1 and next_layer!=-1){
				int new_quality=7;
				if(SL3metaPrimitive->quality<=2) new_quality=5;
				    
				int wi1=-1;int tdc1=-1;
				int wi2=-1;int tdc2=-1;
				int wi3=-1;int tdc3=-1;
				int wi4=-1;int tdc4=-1;
				    
				if(next_layer==1) {wi1=next_wire; tdc1=next_tdc; }
				if(next_layer==2) {wi2=next_wire; tdc2=next_tdc; }
				if(next_layer==3) {wi3=next_wire; tdc3=next_tdc; }
				if(next_layer==4) {wi4=next_wire; tdc4=next_tdc; }
				    
				if(best_layer==1) {wi1=best_wire; tdc1=best_tdc; }
				if(best_layer==2) {wi2=best_wire; tdc2=best_tdc; }
				if(best_layer==3) {wi3=best_wire; tdc3=best_tdc; }
				if(best_layer==4) {wi4=best_wire; tdc4=best_tdc; } 
				    
				    
				    
				correlatedMetaPrimitives.push_back(metaPrimitive({ChId.rawId(),SL3metaPrimitive->t0,SL3metaPrimitive->x,SL3metaPrimitive->tanPhi,SL3metaPrimitive->phi,SL3metaPrimitive->phiB,SL3metaPrimitive->chi2,
						new_quality,
						wi1,tdc1,
						wi2,tdc2,
						wi3,tdc3,
						wi4,tdc4,
						SL3metaPrimitive->wi1,SL3metaPrimitive->tdc1,
						SL3metaPrimitive->wi2,SL3metaPrimitive->tdc2,
						SL3metaPrimitive->wi3,SL3metaPrimitive->tdc3,
						SL3metaPrimitive->wi4,SL3metaPrimitive->tdc4
						}));
				at_least_one_correlation=true;
			    }
			}
		    }
		
		    //finish SL3-SL1

		    if(at_least_one_correlation==false){
			if(debug) std::cout<<"correlation we found zero correlations, adding both collections as they are to the correlatedMetaPrimitives"<<std::endl;
			if(debug) std::cout<<"correlation sizes:"<<SL1metaPrimitives.size()<<" "<<SL3metaPrimitives.size()<<std::endl;
			for (auto SL1metaPrimitive = SL1metaPrimitives.begin(); SL1metaPrimitive != SL1metaPrimitives.end(); ++SL1metaPrimitive){
			    DTSuperLayerId SLId(SL1metaPrimitive->rawId);
			    DTChamberId(SLId.wheel(),SLId.station(),SLId.sector());
			    correlatedMetaPrimitives.push_back(metaPrimitive({ChId.rawId(),SL1metaPrimitive->t0,SL1metaPrimitive->x,SL1metaPrimitive->tanPhi,SL1metaPrimitive->phi,SL1metaPrimitive->phiB,SL1metaPrimitive->chi2,SL1metaPrimitive->quality,
					    SL1metaPrimitive->wi1,SL1metaPrimitive->tdc1,
					    SL1metaPrimitive->wi2,SL1metaPrimitive->tdc2,
					    SL1metaPrimitive->wi3,SL1metaPrimitive->tdc3,
					    SL1metaPrimitive->wi4,SL1metaPrimitive->tdc4,
					    -1,-1,
					    -1,-1,
					    -1,-1,
					    -1,-1
					    }));
			}
			for (auto SL3metaPrimitive = SL3metaPrimitives.begin(); SL3metaPrimitive != SL3metaPrimitives.end(); ++SL3metaPrimitive){
			    DTSuperLayerId SLId(SL3metaPrimitive->rawId);
			    DTChamberId(SLId.wheel(),SLId.station(),SLId.sector());
			    correlatedMetaPrimitives.push_back(metaPrimitive({ChId.rawId(),SL3metaPrimitive->t0,SL3metaPrimitive->x,SL3metaPrimitive->tanPhi,SL3metaPrimitive->phi,SL3metaPrimitive->phiB,SL3metaPrimitive->chi2,SL3metaPrimitive->quality,
					    -1,-1,
					    -1,-1,
					    -1,-1,
					    -1,-1,
					    SL3metaPrimitive->wi1,SL3metaPrimitive->tdc1,
					    SL3metaPrimitive->wi2,SL3metaPrimitive->tdc2,
					    SL3metaPrimitive->wi3,SL3metaPrimitive->tdc3,
					    SL3metaPrimitive->wi4,SL3metaPrimitive->tdc4
					    }));
			}
		    }

		    SL1metaPrimitives.clear();
		    SL1metaPrimitives.erase(SL1metaPrimitives.begin(),SL1metaPrimitives.end());
		    SL3metaPrimitives.clear();
		    SL3metaPrimitives.erase(SL3metaPrimitives.begin(),SL3metaPrimitives.end());
		}
	    }
	}
	if(pinta) NcorrelatedMetaPrimitives->Fill(correlatedMetaPrimitives.size());

	filteredMetaPrimitives.clear();
	filteredMetaPrimitives.erase(filteredMetaPrimitives.begin(),filteredMetaPrimitives.end());

	if(debug) std::cout<<"DTp2 in event:"<<iEvent.id().event()<<" we found "<<correlatedMetaPrimitives.size()<<" correlatedMetPrimitives (chamber)"<<std::endl;
	

	//Trying confirmation by RPCs
	
	//for (auto metaPrimitiveIt = correlatedMetaPrimitives.begin(); metaPrimitiveIt != correlatedMetaPrimitives.end(); ++metaPrimitiveIt){
	//}
	
	vector<L1MuDTChambPhDigi> outPhiCH;
	vector<L1MuDTChambDigi> outP2CH;
	vector<L1Phase2MuDTPhDigi> outP2PhCH;
	
	for (auto metaPrimitiveIt = correlatedMetaPrimitives.begin(); metaPrimitiveIt != correlatedMetaPrimitives.end(); ++metaPrimitiveIt){
	    DTChamberId chId((*metaPrimitiveIt).rawId);
	    if(debug) std::cout<<"looping in final vector: SuperLayerId"<<chId<<" x="<<(*metaPrimitiveIt).x<<" quality="<<(*metaPrimitiveIt).quality<<std::endl;

	    int sectorTP=chId.sector();
	    if(sectorTP==13) sectorTP=4;
	    if(sectorTP==14) sectorTP=10;
	    sectorTP=sectorTP-1;
	    
	    L1MuDTChambPhDigi thisTP((*metaPrimitiveIt).t0,
				      chId.wheel(),
				      sectorTP,
				      chId.station(),
				      (int)round((*metaPrimitiveIt).phi*65536./0.8),
				      (int)round((*metaPrimitiveIt).phiB*2048./1.4),
				      (*metaPrimitiveIt).quality,
				      1,
				      0
				      );
	    
	    if(pinta){
		all_observed_tanPsi[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).tanPhi);
		all_observed_x[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).x);
		all_observed_t0[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).t0);
 		chi2[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).chi2);
 		TPphi[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).phi);
 		TPphiB[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).phiB);
		
		float x_back = trigPos((*metaPrimitiveIt));
		float psi_back = trigDir((*metaPrimitiveIt));
		
		MP_x_back[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill((*metaPrimitiveIt).x,x_back);
 		MP_psi_back[chId.wheel()+2][chId.station()-1][chId.sector()-1][(*metaPrimitiveIt).quality-1]->Fill( TMath::ATan((*metaPrimitiveIt).tanPhi) ,psi_back);

		if(debug)std::cout<<"back:(x,x_back)= "<<(*metaPrimitiveIt).x<<","<<x_back<<std::endl;
		if(debug)std::cout<<"back:(psi,psi_back)= "<<TMath::ATan((*metaPrimitiveIt).tanPhi)<<","<<psi_back<<std::endl;
	    }

	    
	    if(p2_df==0){
		outPhiCH.push_back(thisTP);
	    }else if(p2_df==1){
		if(debug)std::cout<<"pushing back slice-test dataformat"<<std::endl;
		
		outP2CH.push_back(L1MuDTChambDigi((int)round((*metaPrimitiveIt).t0/25.),
						  chId.wheel(),
						  sectorTP,
						  chId.station(),
						  (int)round((*metaPrimitiveIt).phi*65536./0.8),
						  (int)round((*metaPrimitiveIt).phiB*2048./1.4),
						  0,
						  0,
						  (*metaPrimitiveIt).quality,
						  0,
						  (int)round((*metaPrimitiveIt).t0),
						  (int)round((*metaPrimitiveIt).chi2),
						  -10
						  ));
	    }else if(p2_df==2){
		if(debug)std::cout<<"pushing back carlo-federica dataformat"<<std::endl;
		
		outP2PhCH.push_back(L1Phase2MuDTPhDigi((int)round((*metaPrimitiveIt).t0/25.),
						       chId.wheel(),
						       sectorTP,
						       chId.station(),
						       (int)round((*metaPrimitiveIt).phi*65536./0.8),
						       (int)round((*metaPrimitiveIt).phiB*2048./1.4),
						       (*metaPrimitiveIt).quality,
						       0,
						       (int)round((*metaPrimitiveIt).t0),
						       (int)round((*metaPrimitiveIt).chi2),
						       -10
						       ));
		
	    }
	}

	if(p2_df==0){ 
	    std::unique_ptr<L1MuDTChambPhContainer> resultPhiCH (new L1MuDTChambPhContainer);
	    resultPhiCH->setContainer(outPhiCH); iEvent.put(std::move(resultPhiCH));
	    outPhiCH.clear();
	    outPhiCH.erase(outPhiCH.begin(),outPhiCH.end());
	}else if(p2_df==1){
	    std::unique_ptr<L1MuDTChambContainer> resultP2CH (new L1MuDTChambContainer);
	    resultP2CH->setContainer(outP2CH); iEvent.put(std::move(resultP2CH));
	    outP2CH.clear();
	    outP2CH.erase(outP2CH.begin(),outP2CH.end());
	}else if(p2_df==2){
	    std::unique_ptr<L1Phase2MuDTPhContainer> resultP2PhCH (new L1Phase2MuDTPhContainer);
	    resultP2PhCH->setContainer(outP2PhCH); iEvent.put(std::move(resultP2PhCH));
	    outP2PhCH.clear();
	    outP2PhCH.erase(outP2PhCH.begin(),outP2PhCH.end());
	}

	if(pinta){

	    //ploting all qualities before correlation
	    for (auto metaPrimitiveIt = correlatedMetaPrimitives.begin(); metaPrimitiveIt != correlatedMetaPrimitives.end(); ++metaPrimitiveIt){
		Nquality->Fill(metaPrimitiveIt->quality);
	    }

	    DTRecSegment4DCollection::const_iterator segment;

	    if(pinta) Nsegments->Fill(all4DSegments->size());
	    
	    //if(debug) std::cout<<"min_phinhits_match_segment="<<min_phinhits_match_segment<<std::endl;
	    for (segment = all4DSegments->begin();segment!=all4DSegments->end(); ++segment){
		if(!segment->hasPhi()) continue;
		if(int(segment->phiSegment()->recHits().size())<min_phinhits_match_segment) continue;
		DTChamberId chId(segment->chamberId());
		
	  
		//filter CH correlated MP
		std::vector<metaPrimitive> CHmetaPrimitives;
		for(auto metaprimitiveIt = correlatedMetaPrimitives.begin();metaprimitiveIt!=correlatedMetaPrimitives.end();++metaprimitiveIt)
		    if(metaprimitiveIt->rawId==chId.rawId())
			CHmetaPrimitives.push_back(*metaprimitiveIt);
	  
		if(debug) std::cout<<"plots: In Chamber "<<chId<<" we have a phi segment and "<<CHmetaPrimitives.size()<<" correlatedMetaPrimitives"<<std::endl;
		if(CHmetaPrimitives.size()==0)continue;
	
		//T0
		float segment_t0=segment->phiSegment()->t0();
		float segment_t0Phase2=segment_t0+bx25;
	  
		//tanPhi
		LocalVector segmentDirection=segment->localDirection();
		float dx=segmentDirection.x();          
		float dz=segmentDirection.z();          
		float segment_tanPhi=dx/dz;          
		//cassert(TMath::ATan(segment_tanPhi)==TMath::ACos(dx));

		//x
		LocalPoint segmentPosition= segment->localPosition();
		//if(debug) std::cout<<"building wireId inside sl loop wire="<<1<<std::endl;
		//DTWireId wireId(wh,st,se,sl,2,1);//sl,la,wi          
		float segment_x=segmentPosition.x();          
		
		int i=-1;
		float minT=9999;
	      
		for(auto metaprimitiveIt = CHmetaPrimitives.begin();metaprimitiveIt!=CHmetaPrimitives.end();++metaprimitiveIt){
		    float deltaT0=metaprimitiveIt->t0-segment_t0Phase2;
		    if(fabs(deltaT0)<minT){
			i=std::distance(CHmetaPrimitives.begin(),metaprimitiveIt);
			minT=fabs(deltaT0);
		    }
		}

		int iwh=chId.wheel()+2;
		int ist=chId.station()-1;
		int ise=chId.sector()-1;
		int iqu=CHmetaPrimitives[i].quality-1;
	      
		expected_tanPsi[iwh][ist][ise]->Fill(segment_tanPhi);
		expected_x[iwh][ist][ise]->Fill(segment_x);
		expected_t0[iwh][ist][ise]->Fill(segment_t0Phase2);
		
		float z1=11.75;
		float z3=-1.*z1;
		//if (chId.station == 3 or chId.station == 4){
		//z1=9.95;
		//z3=-13.55;
		//}

		if (chId.station()==3 or chId.station()==4) segment_x = segment_x-segment_tanPhi*1.8; //extrapolating segment position from chamber reference frame to chamber middle SL plane in MB3&MB4
		
		if(!(CHmetaPrimitives[i].quality == 9 or CHmetaPrimitives[i].quality == 8 or CHmetaPrimitives[i].quality == 6)){
		    if(inner(CHmetaPrimitives[i])) segment_x = segment_x+segment_tanPhi*z1;
		    if(outer(CHmetaPrimitives[i])) segment_x = segment_x+segment_tanPhi*z3;
		}
		
		if(minT<min_dT0_match_segment){//the closest segment should be within min_dT0_match_segment 
		    observed_tanPsi[iwh][ist][ise][iqu]->Fill(segment_tanPhi);
		    observed_x[iwh][ist][ise][iqu]->Fill(segment_x);
		    observed_t0[iwh][ist][ise][iqu]->Fill(segment_t0Phase2);
		    
		    if(debug) std::cout<<"seg mpm "<<chId<<" -> "
				       <<segment_x<<" "<<CHmetaPrimitives[i].x<<" "	  
				       <<segment_tanPhi<<" "<<CHmetaPrimitives[i].tanPhi<<" "	    
				       <<segment_t0Phase2<<" "<<CHmetaPrimitives[i].t0<<" "<<std::endl;	  
		    
		    //correlation and matched plots
		    segment_vs_jm_x[iwh][ist][ise][iqu]->Fill(segment_x,CHmetaPrimitives[i].x);	  
		    segment_vs_jm_tanPhi[iwh][ist][ise][iqu]->Fill(segment_tanPhi,CHmetaPrimitives[i].tanPhi);
		    segment_vs_jm_T0[iwh][ist][ise][iqu]->Fill(segment_t0Phase2,CHmetaPrimitives[i].t0);
		    
		    segment_vs_jm_x_gauss[iwh][ist][ise][iqu]->Fill(segment_x-CHmetaPrimitives[i].x);
		    segment_vs_jm_tanPhi_gauss[iwh][ist][ise][iqu]->Fill(segment_tanPhi-CHmetaPrimitives[i].tanPhi);
		    segment_vs_jm_T0_gauss[iwh][ist][ise][iqu]->Fill(segment_t0Phase2-CHmetaPrimitives[i].t0);
		    segment_vs_jm_T0_gauss_all[iwh][ist][ise][iqu]->Fill(segment_t0Phase2-CHmetaPrimitives[i].t0);

		    Nquality_matched->Fill(CHmetaPrimitives[i].quality);
		    Nhits_segment_tp->Fill(segment->phiSegment()->recHits().size(),CHmetaPrimitives[i].quality);
		}else{
		    //segment could not be matched
		    if(debug) std::cout<<segment_x<<" "<<segment_tanPhi<<" "<<segment_t0Phase2<<" "<<std::endl;
		    segment_vs_jm_T0_gauss_all[iwh][ist][ise][iqu]->Fill(segment_t0Phase2-CHmetaPrimitives[i].t0);
		    Nhits_segment_tp->Fill(segment->phiSegment()->recHits().size(),0);
		    if(segment->phiSegment()->recHits().size()==4)
			if(debug)std::cout<<chId<<" ineficient event with 4 hits segments in event"<<iEvent.id().event()<<endl;

		}
	    }
	    
	    correlatedMetaPrimitives.clear();
	    correlatedMetaPrimitives.erase(correlatedMetaPrimitives.begin(),correlatedMetaPrimitives.end());
	}
    }
}

void DTTrigPhase2Prod::endRun(edm::Run const& iRun, const edm::EventSetup& iEventSetup) {
  grouping_obj->finish();
  mpathanalyzer->finish();
};



int DTTrigPhase2Prod::arePrimos(metaPrimitive primera, metaPrimitive segunda) {
    if(primera.rawId!=segunda.rawId) return 0;
    if(primera.wi1==segunda.wi1 and primera.tdc1==segunda.tdc1 and primera.wi1!=-1 and primera.tdc1!=-1) return 1;
    if(primera.wi2==segunda.wi2 and primera.tdc2==segunda.tdc2 and primera.wi2!=-1 and primera.tdc2!=-1) return 2;
    if(primera.wi3==segunda.wi3 and primera.tdc3==segunda.tdc3 and primera.wi3!=-1 and primera.tdc3!=-1) return 3;
    if(primera.wi4==segunda.wi4 and primera.tdc4==segunda.tdc4 and primera.wi4!=-1 and primera.tdc4!=-1) return 4;
    return 0;
}


void DTTrigPhase2Prod::printmP(metaPrimitive mP){
    DTSuperLayerId slId(mP.rawId);
    std::cout<<slId<<"\t"
	     <<" "<<setw(2)<<left<<mP.wi1
	     <<" "<<setw(2)<<left<<mP.wi2
	     <<" "<<setw(2)<<left<<mP.wi3
	     <<" "<<setw(2)<<left<<mP.wi4
	     <<" "<<setw(5)<<left<<mP.tdc1
	     <<" "<<setw(5)<<left<<mP.tdc2
	     <<" "<<setw(5)<<left<<mP.tdc3
	     <<" "<<setw(5)<<left<<mP.tdc4
	     <<" "<<setw(10)<<right<<mP.x
	     <<" "<<setw(9)<<left<<mP.tanPhi
	     <<" "<<setw(5)<<left<<mP.t0
	     <<" "<<setw(13)<<left<<mP.chi2
	     <<" r:"<<rango(mP);
}

int DTTrigPhase2Prod::rango(metaPrimitive primera) {
    int rango=0;
    if(primera.wi1!=-1)rango++;
    if(primera.wi2!=-1)rango++;
    if(primera.wi3!=-1)rango++;
    if(primera.wi4!=-1)rango++;
    return rango;
}


bool DTTrigPhase2Prod::outer(metaPrimitive primera){
    if(primera.wi1==-1 and primera.wi2==-1 and primera.wi3==-1 and primera.wi4==-1)
	return true;
    return false;
}

bool DTTrigPhase2Prod::inner(metaPrimitive primera){
    return !outer(primera);
}


bool DTTrigPhase2Prod::hasPosRF(int wh,int sec){
    return  wh>0 || (wh==0 && sec%4>1);
}

double DTTrigPhase2Prod::trigDir(metaPrimitive mp){
    DTChamberId chId(mp.rawId);
    int wh   = chId.wheel();
    int sec  = chId.sector();
    double phi  = mp.phi;
    double phib  = mp.phiB;    
    //double dir = (phib/512.+phi/4096.);
    double dir = (phib+phi);
    //change sign in case of negative wheels
    if (!hasPosRF(wh,sec)) { dir = -dir; }
    return dir;
}

double DTTrigPhase2Prod::trigPos(metaPrimitive mp){
    DTChamberId chId(mp.rawId);

    if(debug) cout<<"back: chId="<<chId<<endl;
    
    int wh   = chId.wheel();
    int sec  = chId.sector();
    int st   = chId.station();
    double phi  = mp.phi;
    double phin = (sec-1)*Geom::pi()/6;
    double phicenter = 0;
    double r = 0;
    double xcenter = 0;
    
    if (sec==4 && st==4) {
	GlobalPoint gpos = phi>0 ? dtGeo->chamber(DTChamberId(wh,st,13))->position() : dtGeo->chamber(DTChamberId(wh,st,4))->position();
	xcenter = phi>0 ? xCenter[0] : -xCenter[0];
	phicenter =  gpos.phi();
	r = gpos.perp();
    } else if (sec==10 && st==4) {
	GlobalPoint gpos = phi>0 ? dtGeo->chamber(DTChamberId(wh,st,14))->position() : dtGeo->chamber(DTChamberId(wh,st,10))->position();
	xcenter = phi>0 ? xCenter[1] : -xCenter[1];
	phicenter =  gpos.phi();
	r = gpos.perp();
    } else {
	GlobalPoint gpos = dtGeo->chamber(DTChamberId(wh,st,sec))->position();
	phicenter =  gpos.phi();
	r = gpos.perp();
    }

    if(debug)cout<<"back: phicenter="<<phicenter<<" phin="<<phicenter<<endl;
   
    double deltaphi = phicenter-phin;
    if(debug)cout<<"back: deltaphi="<<deltaphi<<endl;
    //double x = (tan(phi/4096.)-tan(deltaphi))*(r*cos(deltaphi) - zcn[st-1]); //zcn is in local coordinates -> z invreases approching to vertex
    double x = (tan(phi)-tan(deltaphi))*(r*cos(deltaphi) - zcn[st-1]); //zcn is in local coordinates -> z invreases approching to vertex
    if(debug)cout<<"back: x="<<x<<endl;
    if (hasPosRF(wh,sec)){ x = -x; } // change sign in case of positive wheels
    if(debug)cout<<"back: hasPosRF="<<hasPosRF(wh,sec)<<endl;
    if(debug)cout<<xcenter<<endl;
    //x+=xcenter; this s the bug found by luigi
    return x;
    
}


//------------------------------------------------------------------
//--- Métodos get / set
//------------------------------------------------------------------
void DTTrigPhase2Prod::setBXTolerance(int t) { bxTolerance = t; }
int  DTTrigPhase2Prod::getBXTolerance(void)  { return bxTolerance; }

void DTTrigPhase2Prod::setChiSquareThreshold(float ch2Thr) {
    chiSquareThreshold = ch2Thr;
}

void DTTrigPhase2Prod::setMinimumQuality(MP_QUALITY q) {
    if (minQuality >= LOWQGHOST) minQuality = q;
}
MP_QUALITY DTTrigPhase2Prod::getMinimumQuality(void) { return minQuality; }




