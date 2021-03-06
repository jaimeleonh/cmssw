#ifndef L1Trigger_DTPhase2Trigger_MuonPathAnalyzerPerSL_cc
#define L1Trigger_DTPhase2Trigger_MuonPathAnalyzerPerSL_cc

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/DTSuperLayerId.h"
#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
#include "DataFormats/DTDigi/interface/DTDigiCollection.h"

#include "L1Trigger/DTPhase2Trigger/interface/muonpath.h"
#include "L1Trigger/DTPhase2Trigger/interface/analtypedefs.h"
#include "L1Trigger/DTPhase2Trigger/interface/constants.h"
#include "L1Trigger/DTPhase2Trigger/interface/MuonPathAnalyzer.h" 

#include "CalibMuon/DTDigiSync/interface/DTTTrigBaseSync.h"
#include "CalibMuon/DTDigiSync/interface/DTTTrigSyncFactory.h"

#include "L1Trigger/DTSectorCollector/interface/DTSectCollPhSegm.h"
#include "L1Trigger/DTSectorCollector/interface/DTSectCollThSegm.h"

#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/DTGeometry/interface/DTLayer.h"

#include <iostream>
#include <fstream>


// ===============================================================================
// Previous definitions and declarations
// ===============================================================================

// ===============================================================================
// Class declarations
// ===============================================================================

class MuonPathAnalyzerPerSL : public MuonPathAnalyzer {
 public:
  // Constructors and destructor
  MuonPathAnalyzerPerSL(const edm::ParameterSet& pset);
  virtual ~MuonPathAnalyzerPerSL();
  
  // Main methods
  void initialise(const edm::EventSetup& iEventSetup);
  void run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &inMpath, std::vector<metaPrimitive> &metaPrimitives);
  void run(edm::Event& iEvent, const edm::EventSetup& iEventSetup, std::vector<MuonPath*> &inMpath, std::vector<MuonPath*> &outMPath) {};

  void finish();
  
  // Other public methods
  void setBXTolerance(int t);
  int getBXTolerance(void);
  
  void setChiSquareThreshold(float ch2Thr);
  
  void setMinimumQuality(MP_QUALITY q);
  MP_QUALITY getMinimumQuality(void);

  bool hasPosRF(int wh,int sec) {    return  wh>0 || (wh==0 && sec%4>1);   };

  // Public attributes
  edm::ESHandle<DTGeometry> dtGeo;

  //ttrig
  edm::FileInPath ttrig_filename;
  std::map<int,float> ttriginfo;
  
  //z 
  edm::FileInPath z_filename;
  std::map<int,float> zinfo;
  
  //shift
  edm::FileInPath shift_filename;
  std::map<int,float> shiftinfo;
  
  int chosen_sl;
  
 private:
  // Private methods
  void analyze(MuonPath *inMPath, std::vector<metaPrimitive> &metaPrimitives);

  void setCellLayout(const int layout[4]);
  void buildLateralities(void);
  bool isStraightPath(LATERAL_CASES sideComb[4]);
  
 
  /* Determina si los valores de 4 primitivas forman una trayectoria
     Los valores tienen que ir dispuestos en el orden de capa:
     0    -> Capa más próxima al centro del detector,
     1, 2 -> Siguientes capas
     3    -> Capa más externa */
  void evaluatePathQuality(MuonPath *mPath);
  void evaluateLateralQuality(int latIdx, MuonPath *mPath,
			      LATQ_TYPE *latQuality);
  /* Función que evalua, mediante el criterio de mean-timer, la bondad
     de una trayectoria. Involucra 3 celdas en 3 capas distintas, ordenadas
     de abajo arriba siguiendo el índice del array.
     Es decir:
     0-> valor temporal de la capa inferior,
     1-> valor temporal de la capa intermedia
     2-> valor temporal de la capa superior
     Internamente implementa diferentes funciones según el paso de la
     partícula dependiendo de la lateralidad por la que atraviesa cada
     celda (p. ej.: LLR => Left (inferior); Left (media); Right (superior))
     
     En FPGA debería aplicarse la combinación adecuada para cada caso,
     haciendo uso de funciones que generen el código en tiempo de síntesis,
     aunque la función software diseñada debería ser exportable directamente
     a VHDL */
  void validate(LATERAL_CASES sideComb[3], int layerIndex[3],
		MuonPath* mPath, PARTIAL_LATQ_TYPE *latq);
  
  int eqMainBXTerm(LATERAL_CASES sideComb[2], int layerIdx[2],
		   MuonPath* mPath);
  
  int eqMainTerm(LATERAL_CASES sideComb[2], int layerIdx[2], MuonPath* mPath,
		 int bxValue);
  
  void getLateralCoeficients(LATERAL_CASES sideComb[2], int *coefs);
  bool sameBXValue(PARTIAL_LATQ_TYPE *latq);
  
  void calculatePathParameters(MuonPath *mPath);
  void calcTanPhiXPosChamber  (MuonPath *mPath);
  void calcCellDriftAndXcoor  (MuonPath *mPath);
  void calcChiSquare          (MuonPath *mPath);
  
  void calcTanPhiXPosChamber3Hits(MuonPath *mPath);
  void calcTanPhiXPosChamber4Hits(MuonPath *mPath);
  
  int getOmittedHit(int idx);
  
  // Private attributes

  /* Combinaciones verticales de 3 celdas sobre las que se va a aplicar el
     mean-timer */
  static const int LAYER_ARRANGEMENTS[4][3];
  
  /* El máximo de combinaciones de lateralidad para 4 celdas es 16 grupos
     Es feo reservar todo el posible bloque de memoria de golpe, puesto que
     algunas combinaciones no serán válidas, desperdiciando parte de la
     memoria de forma innecesaria, pero la alternativa es complicar el
     código con vectores y reserva dinámica de memoria y, ¡bueno! ¡si hay
     que ir se va, pero ir p'a n'á es tontería! */
  LATERAL_CASES lateralities[16][4];
  LATQ_TYPE latQuality[16];
  
  int totalNumValLateralities;
  /* Posiciones horizontales de cada celda (una por capa), en unidades de
     semilongitud de celda, relativas a la celda de la capa inferior
     (capa 0). Pese a que la celda de la capa 0 siempre está en posición
     0 respecto de sí misma, se incluye en el array para que el código que
     hace el procesamiento sea más homogéneo y sencillo */

  int bxTolerance;
  MP_QUALITY minQuality;
  float chiSquareThreshold;
  Bool_t debug;
  double chi2Th;
  double chi2corTh;
  double tanPhiTh;
  int cellLayout[4];
  
};


#endif
