import sys
import ROOT
ROOT.gROOT.SetBatch(True)

df = ROOT.RDataFrame("Events", "{OUTPUTFILENAME}")

branches = [
    "nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass",
    "nTau", "Tau_pt", "Tau_eta", "Tau_phi", "Tau_mass", "Tau_dz",
    "Tau_idDeepTau2017v2p1VSmu", "Tau_idDeepTau2017v2p1VSe", "Tau_idDeepTau2017v2p1VSjet",
    "nJet", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId",
    "HLT_IsoMu24_eta2p1",
    "L1_SingleMu22", "L1_Mu18er2p1_Tau26er2p1", 
    "L1_Mu18er2p1_Tau26er2p1_Jet55", "L1_Mu18er2p1_Tau26er2p1_Jet70",
    "nL1Obj", "L1Obj_pt", "L1Obj_eta", "L1Obj_phi", "L1Obj_iso", "L1Obj_type",
]

df.Snapshot("Events", "{OUTPUTPATH}/skimmed_{OUTPUTFILENAME}", tuple(branches))