//-------------------------------------------------
//
//   Class L1Phase2MuDTPhDigi	
//
//   Description: trigger primtive data for the
//                muon barrel Phase2 trigger
//
//
//   Author List: Federica Primavera  Bologna INFN
//
//
//--------------------------------------------------
#ifndef L1Phase2MuDTPhDigi_H
#define L1Phase2MuDTPhDigi_H

//------------------------------------
// Collaborating Class Declarations --
//------------------------------------

#include <DataFormats/DTDigi/interface/DTDigiCollection.h>
//#include <DataFormats/MuonData/interface/MuonDigiCollection.h>
/*#include <DataFormats/MuonDetId/interface/DTLayerId.h>
#include <DataFormats/DTDigi/interface/DTDigi.h>
#include <DataFormats/MuonData/interface/MuonDigiCollection.h>
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/Ptr.h"
#include <vector>*/


using namespace edm; 
using namespace std; 

//----------------------
// Base Class Headers --
//----------------------


//---------------
// C++ Headers --
//---------------


//              ---------------------
//              -- Class Interface --
//              ---------------------


class L1Phase2MuDTPhDigi 
{
 
 public:

  //  Constructors
  L1Phase2MuDTPhDigi();
  
  L1Phase2MuDTPhDigi( int bx, int wh, int sc, int st, int sl, int phi, int phib,
		      int qual, int idx, int t0, int chi2, int rpc=-10);

//  L1Phase2MuDTPhDigi( int bx, int wh, int sc, int st, int sl, int phi, int phib,
//		      int qual, int idx, int t0, int chi2, int rpc=-10, std::vector< RefDTDigi_t > digis);
  
  // Operations
  int bxNum()       const;

  int whNum()       const;
  int scNum()       const;
  int stNum()       const;
  int slNum()       const;

  int phi()         const;
  int phiBend()     const;

  int quality()     const;
  int index()       const;
   
  int t0()          const;
  int chi2()        const;

  int rpcFlag()      const;
  
  //RefDTDigis digis() const; 

 private:

  int m_bx;
  int m_wheel;
  int m_sector;
  int m_station;
  int m_superlayer;

  int m_phiAngle;
  int m_phiBending;

  int m_qualityCode;
  int m_index;
  
  int m_t0;
  int m_chi2;
  
  int m_rpcFlag;

//  edm::Ref <DTDigiCollection, DTDigi> m_digis;
  //std::vector <edm::Ref <DTDigiCollection, DTDigi>> m_digis;
   RefDTDigis m_digis;
};

#endif
