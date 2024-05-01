from PhysicsTools.NanoAOD.common_cff import *
from DPGAnalysis.L1TNanoAOD.l1tnanotables_cff import *
from PhysicsTools.NanoAOD.l1trig_cff import *
from PhysicsTools.NanoAOD.nano_cff import *

l1tnanoMetadata = cms.EDProducer("UniqueStringProducer",
    strings = cms.PSet(
        tag = cms.string("untagged"),
    )
)

l1tNanoTask = cms.Task(nanoMetadata,l1TablesTask)

l1tNanoSequence = cms.Sequence(l1tNanoTask)

def addEmulObjects(process):

    process.l1tNanoTask.add(l1EmulObjTablesTask)
    
    return process


def addUnpackedCaloTPs(process):

    process.l1tNanoTask.add(process.l1CaloTPsNanoTask)
    
    return process

def addEmulCaloTPs(process):

    process.l1tNanoTask.add(process.l1EmulCaloTPsNanoTask)

    return process

def addUnpackedCaloLayer1(process):

    process.l1tNanoTask.add(process.l1CaloLayer1NanoTask)

    return process

def addEmulCaloLayer1(process):

    process.l1tNanoTask.add(process.l1EmulCaloLayer1NanoTask)
         
    return process

def addUnpackedCaloTPsandLayer1(process):

    addUnpackedCaloTPs(process)
    addUnpackedCaloLayer1(process)

    return process

def addEmulCaloTPsandLayer1(process):

    addEmulCaloTPs(process)
    addEmulCaloLayer1(process)

    return process

def addCaloFull(process):

    addEmulCaloTPsandLayer1(process)
    addUnpackedCaloTPsandLayer1(process)
    addEmulObjects(process)

    return process

def addGenParticles(process):
    process.genParticleTable.externalVariables = cms.PSet()

    process.load("PhysicsTools.PatAlgos.slimming.genParticles_cff")
    process.l1tNanoTask.add(process.genParticlesTask)
    process.load("PhysicsTools.NanoAOD.genparticles_cff")
    process.l1tNanoTask.add(process.genParticleTask)
    process.l1tNanoTask.add(process.genParticleTablesTask)
    process.l1tNanoSequence.insert(0, process.finalGenParticles)

    process.load("PhysicsTools.NanoAOD.jetMC_cff")
    process.l1tNanoTask.add(process.genJetTable)
    from PhysicsTools.PatAlgos.slimming.slimmedGenJets_cfi import slimmedGenJets
    process.slimmedGenJets = slimmedGenJets.clone()
    process.l1tNanoSequence.insert(1, process.slimmedGenJets)

    return process

'''
l1tNanoTask = cms.Task(
    #nanoMetadata, 
    l1CaloTPsNanoTask,
    l1CaloLayer1NanoTask,
    l1EmulCaloTPsNanoTask,
    l1EmulCaloLayer1NanoTask,
)
'''
