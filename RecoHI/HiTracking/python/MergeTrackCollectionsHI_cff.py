import FWCore.ParameterSet.Config as cms

import RecoTracker.FinalTrackSelectors.trackListMerger_cfi
hiGeneralTracksNoRegitMu = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = (cms.InputTag('hiGlobalPrimTracks'),
                      cms.InputTag('hiDetachedTripletStepTracks'),
                      cms.InputTag('hiLowPtTripletStepTracks'),
                      cms.InputTag('hiPixelPairGlobalPrimTracks'),
                      cms.InputTag('hiJetCoreRegionalStepTracks')
                     ),
    hasSelector=cms.vint32(1,1,1,1,1),
    selectedTrackQuals = cms.VInputTag(
    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
    cms.InputTag("hiDetachedTripletStepSelector","hiDetachedTripletStep"),
    cms.InputTag("hiLowPtTripletStepSelector","hiLowPtTripletStep"),
    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
    ),                    
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3), pQual=cms.bool(True)),  # should this be False?
                             ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
    )
from Configuration.Eras.Modifier_trackingPhase1_cff import trackingPhase1
trackingPhase1.toModify(hiGeneralTracksNoRegitMu,
    TrackProducers = (cms.InputTag('hiGlobalPrimTracks'),
                      cms.InputTag('hiLowPtQuadStepTracks'),
                      cms.InputTag('hiHighPtTripletStepTracks'),
                      cms.InputTag('hiDetachedQuadStepTracks'),
                      cms.InputTag('hiDetachedTripletStepTracks'),
                      cms.InputTag('hiLowPtTripletStepTracks'),
                      cms.InputTag('hiPixelPairGlobalPrimTracks'),
                      cms.InputTag('hiJetCoreRegionalStepTracks')
                     ),
    hasSelector=cms.vint32(1,1,1,1,1,1,1,1),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5,6), pQual=cms.bool(True))),
    selectedTrackQuals = cms.VInputTag(
    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
    cms.InputTag("hiLowPtQuadStepSelector","hiLowPtQuadStep"),
    cms.InputTag("hiHighPtTripletStepSelector","hiHighPtTripletStep"),
    cms.InputTag("hiDetachedQuadStepSelector","hiDetachedQuadStep"),
    cms.InputTag("hiDetachedTripletStepSelector","hiDetachedTripletStep"),
    cms.InputTag("hiLowPtTripletStepSelector","hiLowPtTripletStep"),
    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
    )                    
)    

hiGeneralTracks = RecoTracker.FinalTrackSelectors.trackListMerger_cfi.trackListMerger.clone(
    TrackProducers = (cms.InputTag('hiGlobalPrimTracks'),
                      cms.InputTag('hiDetachedTripletStepTracks'),
                      cms.InputTag('hiLowPtTripletStepTracks'),
                      cms.InputTag('hiPixelPairGlobalPrimTracks'),
                      cms.InputTag('hiJetCoreRegionalStepTracks'),
                      cms.InputTag('hiRegitMuInitialStepTracks'),
                      cms.InputTag('hiRegitMuPixelPairStepTracks'),
                      cms.InputTag('hiRegitMuMixedTripletStepTracks'),
                      cms.InputTag('hiRegitMuPixelLessStepTracks'),
                      cms.InputTag('hiRegitMuDetachedTripletStepTracks'),
                      cms.InputTag('hiRegitMuonSeededTracksOutIn'),
                      cms.InputTag('hiRegitMuonSeededTracksInOut')
                     ),
    hasSelector=cms.vint32(1,1,1,1,1,1,1,1,1,1,1,1),
    selectedTrackQuals = cms.VInputTag(
    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
    cms.InputTag("hiDetachedTripletStepSelector","hiDetachedTripletStep"),
    cms.InputTag("hiLowPtTripletStepSelector","hiLowPtTripletStep"),
    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
    cms.InputTag("hiJetCoreRegionalStepSelector","hiJetCoreRegionalStep"),
    cms.InputTag("hiRegitMuInitialStepSelector","hiRegitMuInitialStepLoose"),
    cms.InputTag("hiRegitMuPixelPairStepSelector","hiRegitMuPixelPairStep"),
    cms.InputTag("hiRegitMuMixedTripletStepSelector","hiRegitMuMixedTripletStep"),
    cms.InputTag("hiRegitMuPixelLessStepSelector","hiRegitMuPixelLessStep"),
    cms.InputTag("hiRegitMuDetachedTripletStepSelector","hiRegitMuDetachedTripletStep"),
    cms.InputTag("hiRegitMuonSeededTracksOutInSelector","hiRegitMuonSeededTracksOutInHighPurity"),
    cms.InputTag("hiRegitMuonSeededTracksInOutSelector","hiRegitMuonSeededTracksInOutHighPurity")
    ),                    
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11), pQual=cms.bool(True)),  # should this be False?
                             ),
    copyExtras = True,
    makeReKeyedSeeds = cms.untracked.bool(False)
)
trackingPhase1.toModify(hiGeneralTracks,
    TrackProducers = (cms.InputTag('hiGlobalPrimTracks'),
                      cms.InputTag('hiLowPtQuadStepTracks'),
                      cms.InputTag('hiHighPtTripletStepTracks'),
                      cms.InputTag('hiDetachedQuadStepTracks'),
                      cms.InputTag('hiDetachedTripletStepTracks'),
                      cms.InputTag('hiLowPtTripletStepTracks'),
                      cms.InputTag('hiPixelPairGlobalPrimTracks'),
                      cms.InputTag('hiMixedTripletStepTracks'),
                      cms.InputTag('hiPixelLessStepTracks'),
                      cms.InputTag('hiTobTecStepTracks'),
                      cms.InputTag('hiJetCoreRegionalStepTracks'),
                      cms.InputTag('hiRegitMuInitialStepTracks'),
                      cms.InputTag('hiRegitMuPixelPairStepTracks'),
                      cms.InputTag('hiRegitMuMixedTripletStepTracks'),
                      cms.InputTag('hiRegitMuPixelLessStepTracks'),
                      cms.InputTag('hiRegitMuDetachedTripletStepTracks'),
                      cms.InputTag('hiRegitMuonSeededTracksOutIn'),
                      cms.InputTag('hiRegitMuonSeededTracksInOut')
                     ),
    hasSelector=cms.vint32(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1),
    setsToMerge = cms.VPSet( cms.PSet( tLists=cms.vint32(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), pQual=cms.bool(True))),  # should this be False?
    selectedTrackQuals = cms.VInputTag(
    cms.InputTag("hiInitialStepSelector","hiInitialStep"),
    cms.InputTag("hiLowPtQuadStepSelector","hiLowPtQuadStep"),
    cms.InputTag("hiHighPtTripletStepSelector","hiHighPtTripletStep"),
    cms.InputTag("hiDetachedQuadStepSelector","hiDetachedQuadStep"),
    cms.InputTag("hiDetachedTripletStepSelector","hiDetachedTripletStep"),
    cms.InputTag("hiLowPtTripletStepSelector","hiLowPtTripletStep"),
    cms.InputTag("hiPixelPairStepSelector","hiPixelPairStep"),
    cms.InputTag("hiMixedTripletStepSelector","hiMixedTripletStep"),
    cms.InputTag("hiPixelLessStepSelector","hiPixelLessStep"),
    cms.InputTag("hiTobTecStepSelector","hiTobTecStep"),
    cms.InputTag("hiJetCoreRegionalStepSelector","hiJetCoreRegionalStep"),
    cms.InputTag("hiRegitMuInitialStepSelector","hiRegitMuInitialStepLoose"),
    cms.InputTag("hiRegitMuPixelPairStepSelector","hiRegitMuPixelPairStep"),
    cms.InputTag("hiRegitMuMixedTripletStepSelector","hiRegitMuMixedTripletStep"),
    cms.InputTag("hiRegitMuPixelLessStepSelector","hiRegitMuPixelLessStep"),
    cms.InputTag("hiRegitMuDetachedTripletStepSelector","hiRegitMuDetachedTripletStep"),
    cms.InputTag("hiRegitMuonSeededTracksOutInSelector","hiRegitMuonSeededTracksOutInHighPurity"),
    cms.InputTag("hiRegitMuonSeededTracksInOutSelector","hiRegitMuonSeededTracksInOutHighPurity")
    )                    
)    
