import FWCore.ParameterSet.Config as cms

# load tau-jet specific JEC parameters from SQLlite file
payloads = [
   'AK5tauHPSlooseCombDBcorr',
   'AK5tauHPSlooseCombDBcorrOneProng0Pi0',
   'AK5tauHPSlooseCombDBcorrOneProng1Pi0',
   'AK5tauHPSlooseCombDBcorrOneProng2Pi0',
   'AK5tauHPSlooseCombDBcorrTwoProng0Pi0',
   'AK5tauHPSlooseCombDBcorrTwoProng1Pi0',
   'AK5tauHPSlooseCombDBcorrThreeProng0Pi0',
   'AK5tauHPSlooseCombDBcorrThreeProng1Pi0'
]    

PoolDBESSource_toGet = []
for payload in payloads:
    PoolDBESSource_toGet.append(cms.PSet(
        record = cms.string('JetCorrectionsRecord'),
        tag    = cms.string('JetCorrectorParametersCollection_TauJecSpring15_V1_%s' % payload),
        label  = cms.untracked.string(payload)
    ))

SQLliteInput = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        messageLevel = cms.untracked.int32(2)
    ),
    timetype = cms.string('runnumber'),
    toGet = cms.VPSet(PoolDBESSource_toGet),
    connect = cms.string('sqlite_fip:JetMETCorrections/Modules/test/TauJecSpring15_V1.db')
)
es_prefer_jec = cms.ESPrefer('PoolDBESSource', 'SQLliteInput')

# produce associated tau-jet energy correction factors in a valuemap
from PhysicsTools.PatAlgos.recoLayer0.tauJetCorrFactors_cfi import *
patTauJetCorrections = cms.Sequence(patTauJetCorrFactors)


