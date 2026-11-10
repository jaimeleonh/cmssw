import FWCore.ParameterSet.Config as cms
from L1TriggerScouting.Phase2.modules import l1sc_L1TScPhase2SCJets_alpaka

from FWCore.ParameterSet.VarParsing import VarParsing
options = VarParsing('analysis')
options.register('backend', 'serial_sync', VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, 'Alpaka backend to test')
options.parseArguments()
process = cms.Process('TEST')
process.load('Configuration.StandardSequences.Accelerators_cff')
process.source = cms.Source('EmptySource')
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(2))
process.options = cms.untracked.PSet(numberOfThreads=cms.untracked.uint32(1),
                                   numberOfStreams=cms.untracked.uint32(1))
process.input = cms.EDProducer('Phase2JetTestInput')
process.onDevice = cms.EDProducer('l1sc::Phase2JetTestCopy@alpaka', src=cms.InputTag('input'),
    alpaka=cms.untracked.PSet(backend=cms.untracked.string(options.backend)))
sequence = cms.Sequence(process.input + process.onDevice)
cases = [
    ('SeededCone', 2., 1, 16),
    ('SeededCone', 2., 1, 1),
    ('SCNMSW', 0., 1, 16),
    ('SCNMSW', 2., 1, 16),
    ('SCNMSW', 2., 2, 16),
    ('SCNMSW', 2., -1, 16),
    ('LinkTree', 2., 1, 16),
]
for case, (algo, alpha, rounds, njets) in enumerate(cases):
    for cut in (0, 2, 3):
        label = f'jets{case}cut{cut}'
        producer = l1sc_L1TScPhase2SCJets_alpaka(
            alpaka=cms.untracked.PSet(backend=cms.untracked.string(options.backend)),
            src=cms.InputTag('onDevice'), algo=cms.string(algo),
            alphaSeed=cms.double(alpha), maxRounds=cms.int32(rounds),
            nJets=cms.uint32(njets), minConstituents=cms.uint32(cut),
        )
        setattr(process, label, producer)
        checker = cms.EDAnalyzer('Phase2JetTestCheck', src=cms.InputTag(label),
                                minConstituents=cms.uint32(cut),
                                capped=cms.bool(algo == 'SeededCone' and njets == 1))
        setattr(process, label + 'check', checker)
        sequence += producer + checker
process.test = cms.Path(sequence)
