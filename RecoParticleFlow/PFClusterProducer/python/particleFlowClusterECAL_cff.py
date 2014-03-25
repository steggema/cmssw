import FWCore.ParameterSet.Config as cms

from RecoParticleFlow.PFClusterProducer.particleFlowClusterECAL_cfi import *
from RecoParticleFlow.PFClusterProducer.particleFlowRecHitECAL_cfi import *

particleFlowClusterECALSequenceOld = cms.Sequence(
    particleFlowRecHitECAL+
    particleFlowClusterECALUncorrected +
    particleFlowClusterECAL
    )
                                                   
