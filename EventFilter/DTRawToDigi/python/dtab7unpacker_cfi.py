import FWCore.ParameterSet.Config as cms

dtAB7unpacker = cms.EDProducer("OglezDTAB7RawToDigi",
                               DTAB7_FED_Source = cms.InputTag("rawDataCollector"),
                               feds = cms.untracked.vint32( 1368,),
                               debug = cms.untracked.bool(False),
)
