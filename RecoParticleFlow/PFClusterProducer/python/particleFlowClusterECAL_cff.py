import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi import *

particleFlowClusterECALSequence = cms.Sequence(
    particleFlowRecHitECAL+
    particleFlowClusterECALUncorrected +
    particleFlowClusterECAL
    )
                                                   
