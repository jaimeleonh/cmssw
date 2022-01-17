a = [
    "hltHpsDoublePFTau30MediumDitauWPDeepTauDz02DoubleTauJet",
    "hltHpsDoublePFTau30MediumDitauWPDeepTauL1HLTMatchedDoubleTauJet",
    "hltHpsPFTauTrack",
    "hltL2DoubleTauTagNNFilterDoubleTauJet",
    "hltPFJets60L1HLTMatched",
    "hltHpsOverlapFilterDeepTauDoublePFTau30PFJet60",
]

for elem in a:
    print 'int pass%s;' % elem
    print 'int %s_nObjs;' % elem
    print 'vector<float> %s_pt;' % elem
    print 'vector<float> %s_eta;' % elem
    print 'vector<float> %s_phi;' % elem
    print 'vector<float> %s_energy;' % elem
    print

print
print
print

for elem in a:
    print 'tree->Branch("pass%s", &pass%s);' % (elem, elem)
    print 'tree->Branch("%s_nObjs", &%s_nObjs);' % (elem, elem)
    print 'tree->Branch("%s_pt", &%s_pt);' % (elem, elem)
    print 'tree->Branch("%s_eta", &%s_eta);' % (elem, elem)
    print 'tree->Branch("%s_phi", &%s_phi);' % (elem, elem)
    print 'tree->Branch("%s_energy", &%s_energy);' % (elem, elem)
    print

print
print
print

for elem in a:
    print 'pass%s = 0;' % elem
    print '%s_nObjs = 0;' % elem
    print '%s_pt.clear();' % elem
    print '%s_eta.clear();' % elem
    print '%s_phi.clear();' % elem
    print '%s_energy.clear();' % elem
    print

print
print
print


for elem in a:
    print 'std::string %s_Tag = "%s::MYHLT";' % (elem, elem)

print 
print
print

for elem in a:
    print ('if (filterTag == %s_Tag && nObjKeys >= 2) {\n'
           '    pass%s = 1;\n'
           '    %s_nObjs = nObjKeys;\n'
           '}' % (elem, elem, elem))

print 
print
print

for elem in a:
    print ('if (filterTag == %s_Tag\n'
           '        && pass%s && pt_>0) {\n'
           '    %s_pt.push_back(pt_);\n'
           '    %s_eta.push_back(eta_);\n'
           '    %s_phi.push_back(phi_);\n'
           '    %s_energy.push_back(energy_);\n'
           '}' % (elem, elem, elem, elem, elem, elem))
   