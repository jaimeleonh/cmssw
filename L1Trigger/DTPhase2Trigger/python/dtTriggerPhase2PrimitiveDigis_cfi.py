import FWCore.ParameterSet.Config as cms

from L1TriggerConfig.DTTPGConfigProducers.L1DTTPGConfigFromDB_cff import *

dtTriggerPhase2PrimitiveDigis = cms.EDProducer("DTTrigPhase2Prod",
                                               digiTag = cms.InputTag("CalibratedDigis"),
                                               trigger_with_sl = cms.untracked.int32(4),
                                               tanPhiTh = cms.untracked.double(10.),
                                               chi2Th = cms.untracked.double(0.01), #in cm^2
                                               do_correlation = cms.untracked.bool(True),
                                               dT0_correlate_TP = cms.untracked.double(25.),
                                               minx_match_2digis = cms.untracked.double(2.1),
                                               p2_df = cms.untracked.int32(2), #1 for phase-1, 2 for slice-test, 3 for phase-2 carlo-federica
                                               filter_cousins = cms.untracked.bool(True),
                                               apply_txt_ttrig_bc0 = cms.untracked.bool(False),
                                               ttrig_filename = cms.untracked.string('data/wire_rawId_ttrig.txt'),
                                               z_filename = cms.untracked.string('data/wire_rawId_z.txt'),
                                               shift_filename = cms.untracked.string('data/wire_rawId_x.txt'),
                                               grouping_code = cms.untracked.int32(2),       # 0 = initial grouping, 1 = Hough transform
                                               min_phinhits_match_segment = cms.untracked.int32(8),
                                               min_dT0_match_segment = cms.untracked.double(12.5),
                                               #debugging
                                               debug = cms.untracked.bool(False),
                                               dump  = cms.untracked.bool(False),
                                               #RPC
                                               rpcRecHits = cms.untracked.InputTag("rpcRecHits"),
                                               ### PseudoBayesPattern parameters ###
                                               pattern_filename = cms.untracked.string("data/patterns.root"),
                                               #Minimum number of layers hit (total). Together with the two parameters under this it means 4+4, 4+3 or 3+3
                                               minNLayerHits   = cms.untracked.int32(6),
                                               #Minimum number of hits in the most hit superlayer
                                               minSingleSLHitsMax = cms.untracked.int32(3),
                                               #Minimum number of hits in the less hit superlayer
                                               minSingleSLHitsMin = cms.untracked.int32(3),
                                               #By default pattern width is 1, 0 can be considered (harder fits but, lower efficiency of high quality), 2 is the absolute limit unless we have extremely bent muons somehow
                                               allowedVariance = cms.untracked.int32(1),
                                               #If true, it will provide all candidate sets with the same hits of the same quality (with lateralities defined). If false only the leading one (with its lateralities).
                                               allowDuplicates = cms.untracked.bool(False),
                                               #Also provide best estimates for the lateralities
                                               setLateralities = cms.untracked.bool(True),
                                               #DTPrimitives are saved in the appropriate element of the muonPath array
                                               saveOnPlace = cms.untracked.bool(True),

                                               )
