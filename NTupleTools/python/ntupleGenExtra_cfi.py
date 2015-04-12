import FWCore.ParameterSet.Config as cms

#from SLHCL1TrackTriggerSimulations.NTupleTools.prunedGenParticles_cfi import prunedGenParticles
from SLHCL1TrackTriggerSimulations.NTupleTools.BeamSpotFromSim_cfi import BeamSpotFromSim

ntupleGenParticlesExtra = cms.EDProducer('NTupleGenParticlesExtra',
    #inputTag = cms.InputTag('prunedGenParticles'),
    inputTag = cms.InputTag('genParticles'),
    prefix = cms.string('genParts@'),
    suffix = cms.string(''),
    cut = cms.string(''),
    maxN = cms.uint32(999999)
)

ntupleGenExtra = cms.Sequence(BeamSpotFromSim * ntupleGenParticlesExtra)

