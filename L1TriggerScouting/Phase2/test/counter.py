import ROOT
df_c = ROOT.RDataFrame("Events", "NanoOutput_cascade_mb.root")
#df_b = ROOT.RDataFrame("Events", "NanoOutput.root")
df_b = ROOT.RDataFrame("Events", "NanoOutput_btag_mb.root")

h_c = df_c.Histo1D("nSC4AlpakaJets")
h_b = df_b.Histo1D("nSC4AlpakaJets")
print("Cascade: ", h_c.Integral() * h_c.GetMean())
print("BTag: ", h_b.Integral() * h_b.GetMean())

print("Cascade:", df_c.Histo1D("SC4AlpakaJets_pt").Integral())
print("BTag:", df_b.Histo1D("SC4AlpakaJets_pt").Integral())
