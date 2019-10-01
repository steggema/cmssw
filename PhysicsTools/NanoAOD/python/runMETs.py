import FWCore.ParameterSet.Config as cms


year = 2018
dataset_type = 'MC'

globalTag = {
    2016 : {
        'MC' : '102X_mcRun2_asymptotic_v6',
        'data' : '102X_dataRun2_nanoAOD_2016_v1',
        'embedding' : '102X_dataRun2_nanoAOD_2016_v1',
    },
    2017 : {
        'MC' : '102X_mc2017_realistic_v6',
        'data' : '102X_dataRun2_v8',
        'embedding' : '102X_dataRun2_v8',
    },
    2018 : {
        'MC' : '102X_upgrade2018_realistic_v18',
        'data' : '102X_dataRun2_Sep2018ABC_v2',
        'data-prompt' : '102X_dataRun2_Prompt_v13',
        'embedding' : '102X_dataRun2_Sep2018ABC_v2',
        'embedding-prompt' : '102X_dataRun2_Prompt_v13',
    },
}


def addMoreMETs(process):
    # Jets: JEC + DeepJet b-taggers + PU Jet ID
    from PhysicsTools.PatAlgos.tools.jetTools import updateJetCollection
    updateJetCollection(
         process,
         jetSource = cms.InputTag('slimmedJets'),
         labelName = 'UpdatedJEC',
         jetCorrections = ('AK4PFchs', cms.vstring(['L1FastJet', 'L2Relative', 'L3Absolute', 'L2L3Residual']), 'None'),
         btagDiscriminators = ['None'
         ],
    )

    process.load("RecoJets.JetProducers.PileupJetID_cfi")
    process.pileupJetIdUpdated = process.pileupJetId.clone(
      jets=cms.InputTag("slimmedJets"),
      inputIsCorrected=True,
      applyJec=True,
      vertexes=cms.InputTag("offlineSlimmedPrimaryVertices")
    )
    process.looseJetId = cms.EDProducer("PatJetIDValueMapProducer",
                filterParams=cms.PSet(
                  version = cms.string('WINTER17'),
                  quality = cms.string('LOOSE'),
                ),
                            src = cms.InputTag("slimmedJets")
    )
    process.tightJetId = cms.EDProducer("PatJetIDValueMapProducer",
                filterParams=cms.PSet(
                  version = cms.string('WINTER17'),
                  quality = cms.string('TIGHT'),
                ),
                            src = cms.InputTag("slimmedJets")
    )
    process.tightJetIdLepVeto = cms.EDProducer("PatJetIDValueMapProducer",
                filterParams=cms.PSet(
                  version = cms.string('WINTER17'),
                  quality = cms.string('TIGHTLEPVETO'),
                ),
                            src = cms.InputTag("slimmedJets")
    )
    process.updatedPatJetsUpdatedJEC.userData.userInts.src += ['pileupJetIdUpdated:fullId','looseJetId', 'tightJetId', 'tightJetIdLepVeto']
    process.updatedPatJetsUpdatedJEC.userData.userFloats.src += ['pileupJetIdUpdated:fullDiscriminant']

    process.jecSequence = cms.Sequence(
      process.pileupJetIdUpdated
      *process.looseJetId
      *process.tightJetId
      *process.tightJetIdLepVeto
      *process.patJetCorrFactorsUpdatedJEC
      *process.updatedPatJetsUpdatedJEC
      # Following two lines not needed when no b-tagging
      # *process.patJetCorrFactorsTransientCorrectedUpdatedJEC
      # *process.updatedPatJetsTransientCorrectedUpdatedJEC
      *process.selectedUpdatedPatJetsUpdatedJEC
    )

    process.nanoSequenceCommon.insert(0, process.jecSequence)

    from PhysicsTools.NanoAOD.prepareMETDefinitions_cff import prepareMETs # Preparations for MVA MET
    prepareMETs(process, "selectedUpdatedPatJetsUpdatedJEC")
    process.neutralInJets.jetPUDIWP = cms.string("medium")
    process.neutralInJets.jetPUIDMapLabel = cms.string("pileupJetIdUpdated:fullId")
    process.nanoSequenceCommon.insert(1, cms.Sequence(process.pfNeutrals*process.pfChargedPV*process.neutralInJets*process.pfChargedPU))
    process.nanoSequenceCommon.insert(2, cms.Sequence(process.pfTrackMETCands*process.pfTrackMET*process.patpfTrackMET*process.pfNoPUMETCands*process.pfNoPUMET*process.patpfNoPUMET*process.pfPUCorrectedMETCands*process.pfPUCorrectedMET*process.patpfPUCorrectedMET*process.pfPUMETCands*process.pfPUMET*process.patpfPUMET))

