import math

from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Muon, GenParticle
from CMGTools.RootTools.utils.DeltaR import matchObjectCollection, deltaR, deltaR2

def isFinal(p):
    return not (p.numberOfDaughters() == 1 and p.daughter(0).pdgId() == p.pdgId())

def deltaRObj(obj1, obj2):
    return deltaR(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())

def deltaR2Obj(obj1, obj2):
    return deltaR2(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())

def getWeights(neutral, event):
    sumPUDr2 = 0.
    sumPVDr2 = 0.
    sumPUDr201 = 0.
    sumPVDr201 = 0.
    sumPUDr = 0.
    sumPVDr = 0.
    sumPUPt = 0.
    sumPVPt = 0.
    sumPUPuppi = 0.
    sumPVPuppi = 0.
    sumPUPuppi05 = 0.
    sumPVPuppi05 = 0.
    sumPUPuppi05Dr2 = 0.
    sumPVPuppi05Dr2 = 0.
    sumPUPuppi05Log2 = 0.
    sumPVPuppi05Log2 = 0.

    for charged in event.chargedPfCandidates:
        dr2 = deltaR2Obj(neutral, charged)
        dr = math.sqrt(dr2)
        if dr2 == 0.:
            import pdb; pdb.set_trace()
            continue
        if charged.fromPV:
            sumPVDr2 += 1./dr2
            if dr2 < 1.:
                sumPVDr201 += 1./dr2 if dr > 0.05 else 1./0.05/0.05
                
            sumPVDr += 1./dr
            sumPVPt += charged.pt()/dr2
            if charged.pt() > dr and charged.pt() > 0.05:
                if dr < 0.05:
                    sumPVPuppi += math.log(charged.pt()/0.05)
                    sumPVPuppi05 += math.log(charged.pt()/0.05)
                    sumPVPuppi05Dr2 += math.log(charged.pt()/0.05/0.05)
                    sumPVPuppi05Log2 += math.log(charged.pt()/0.05)**2
                else:
                    sumPVPuppi += math.log(charged.pt()/dr)
                    if dr < 0.5:
                        sumPVPuppi05 += math.log(charged.pt()/dr)
                        sumPVPuppi05Dr2 += math.log(charged.pt()/dr2)
                        sumPVPuppi05Log2 += math.log(charged.pt()/dr)**2
        else:
            sumPUDr2+= 1./dr2
            if dr2 < 1.:
                sumPUDr201 += 1./dr2 if dr > 0.05 else 1./0.05/0.05
            sumPUDr += 1./dr
            sumPUPt += charged.pt()/dr2
            if charged.pt() > dr and charged.pt() > 0.05:
                if dr < 0.05:
                    sumPUPuppi += math.log(charged.pt()/0.05)
                    sumPUPuppi05 += math.log(charged.pt()/dr)
                    sumPUPuppi05Dr2 += math.log(charged.pt()/0.05/0.05)**2
                    sumPUPuppi05Log2 += math.log(charged.pt()/dr)**2
                else:
                    sumPUPuppi += math.log(charged.pt()/dr)
                    if dr < 0.5:
                        sumPUPuppi05 += math.log(charged.pt()/dr)
                        sumPUPuppi05Dr2 += math.log(charged.pt()/dr2)
                        sumPUPuppi05Log2 += math.log(charged.pt()/dr)**2

    neutral.deltaBetaWeightDr2 = sumPVDr2/(sumPVDr2 + sumPUDr2)
    neutral.deltaBetaWeightDr201 = sumPVDr201/(sumPVDr201 + sumPUDr201)
    neutral.deltaBetaWeightDr = sumPVDr/(sumPVDr + sumPUDr)
    neutral.deltaBetaWeightPt = sumPVPt/(sumPVPt + sumPUPt)
    neutral.deltaBetaWeightLogPt = math.log(sumPVPt+1.)/(math.log(sumPVPt+1.) + math.log(sumPUPt+1.)) 
    neutral.deltaBetaWeightPuppi = sumPVPuppi/(sumPVPuppi + sumPUPuppi)
    neutral.deltaBetaWeightPuppi05 = sumPVPuppi05/(sumPVPuppi05 + sumPUPuppi05) if sumPVPuppi05 > 0. else 0.
    neutral.deltaBetaWeightPuppi05Dr2 = sumPVPuppi05Dr2/(sumPVPuppi05Dr2 + sumPUPuppi05Dr2) if sumPVPuppi05Dr2 > 0. else 0.
    neutral.deltaBetaWeightPuppi05Log2 = sumPVPuppi05Log2/(sumPVPuppi05Log2 + sumPUPuppi05Log2) if sumPVPuppi05Log2 > 0. else 0.


class MuonAnalyzer( Analyzer ):
    ''' Basic muon selection and matching to generated muons.
    Saves selected muons as event.muons
    
    genSrc = 'genParticlesPruned',
    absGenIds = [13],
    muonSrc = 'cmgMuonSel',
    minPt = 20.,
    maxEta = 2.4,
    matchDeltaR = 0.3,
    '''

    def declareHandles(self):

        super(MuonAnalyzer, self).declareHandles()
        self.mchandles['genParticles'] =  AutoHandle(
            self.cfg_ana.genSrc,
            'std::vector<reco::GenParticle>'
            )

        self.handles['leptons'] = AutoHandle(
            self.cfg_ana.muonSrc,
            'std::vector<reco::Muon>'
            )

        self.handles['pfCandidates'] = AutoHandle(
            'particleFlow',
            'std::vector<reco::PFCandidate>'
            )

    def process(self, iEvent, event):

        super(MuonAnalyzer, self).process(iEvent, event) # This reads collections

        # import pdb; pdb.set_trace()

        muons = self.handles['leptons'].product()

        selMuons = [Muon(muon) for muon in muons if muon.pt() > self.cfg_ana.minPt and abs(muon.eta()) < self.cfg_ana.maxEta and muon.isTrackerMuon()]

        event.muons = selMuons

        for muon in event.muons:
            muon.associatedVertex = event.goodVertices[0]

        # vertices = event.goodVertices
        event.pfCandidates = self.handles['pfCandidates'].product()
        event.chargedPfCandidates = [p for p in event.pfCandidates if p.particleId() in [1, 2]] # electron or CH
        event.neutralPfCandidates = [p for p in event.pfCandidates if p.particleId() not in [1, 2]] # electron or CH

        vertexZs = [v.z() for v in event.goodVertices]

        # Calculate PV weight for each charged PFCandidate
        for p in event.chargedPfCandidates:
            trackWeights = [v.trackWeight(p.trackRef()) for v in event.goodVertices]
            maxWeight = max(trackWeights)
            sumWeights = sum(trackWeights)
            # print 'Max/sum weights', maxWeight, sumWeights
            p.pvWeight = trackWeights[0]/sumWeights if sumWeights else 0.
            if p.pvWeight and p.pvWeight != 1.:
                print 'Track shared by two vertices; PV weight', p.pvWeight
            if sumWeights:
                p.fromPV = (maxWeight == trackWeights[0])
            else:
                # Check closest z vertex
                p.fromPV = False
                deltaZs = [abs(vz-p.vertex().z()) for vz in vertexZs]

                if p.vertex().z() == 0.:
                    if p.trackRef().get():
                        deltaZs = [abs(vz-p.trackRef().vertex().z()) for vz in vertexZs]
                    elif p.gsfTrackRef().get():
                        deltaZs = [abs(vz-p.gsfTrackRef().vertex().z()) for vz in vertexZs]

                minDeltaZ = min(deltaZs)
                if minDeltaZ == deltaZs[0]:
                    if deltaZs[0] < 2.: # no dramatic outliers if vertex is at the edge
                        p.fromPV = True
                    # else:
                    #     print 'Associate PF candidate by closest vertex, but large deltaZ', deltaZs[0]
                    #     print deltaZs
        
        # for neutral in event.neutralPfCandidates:
        #     getWeights(neutral, event)

        event.sumPtPV = sum(p.pt() for p in event.chargedPfCandidates if p.fromPV)
        event.sumPtPU = sum(p.pt() for p in event.chargedPfCandidates if not p.fromPV)
        event.sumPtNeutral = sum(p.pt() for p in event.neutralPfCandidates)
        event.sumPtPhotons = sum(p.pt() for p in event.neutralPfCandidates if p.particleId() == 4)
        event.sumPtAll= sum(p.pt() for p in event.pfCandidates)

        if self.cfg_comp.isMC:
            # print event.eventId
            
            if not hasattr(event, 'genParticlesTop'):
                genParticles = self.mchandles['genParticles'].product()
                event.genParticles = map( GenParticle, genParticles)
                genParticles = [p for p in event.genParticles if isFinal(p)]
            else:
                genParticles = event.genParticlesTop
            # Allow Higgs/Z/photon/W
            # allowedGenMothers = [6, 15, 21, 23, 24, 25, 35, 36, 37]
            # allowedGenMothers = [6, 15, 24]
            event.generatedMuons = [p for p in genParticles if abs(p.pdgId()) in self.cfg_ana.absGenIds]

            pairs = matchObjectCollection(event.muons, event.generatedMuons, self.cfg_ana.matchDeltaR)

            for muon in event.muons:
                genMuon = pairs[muon]
                muon.parton = genMuon
                self.testMuonIso(muon, 0.1)
                muon.id = self.testMuonID(muon)
                if genMuon:
                    # Can also be done in tree filler:
                    muon.genMatchDeltaR = pairs[muon], deltaRObj(muon, pairs[muon])
                    genMuon.recMuon = muon

                #ISO BLOCK
                muon.phIsoDr = 0.
                muon.phIsoDr2 = 0.
                muon.phIsoDr201 = 0.
                muon.phIsoPt = 0.
                muon.phIsoLogPt = 0.
                muon.phIsoPuppi = 0.
                muon.phIsoPuppi05 = 0.
                muon.phIsoPuppi05Dr2 = 0.
                muon.phIsoPuppi05Log2 = 0.
                muon.nhIsoDr = 0.
                muon.nhIsoDr2 = 0.
                muon.nhIsoDr201 = 0.
                muon.nhIsoPt = 0.
                muon.nhIsoLogPt = 0.
                muon.nhIsoPuppi = 0.
                muon.nhIsoPuppi05 = 0.
                muon.nhIsoPuppi05Dr2 = 0.
                muon.nhIsoPuppi05Log2 = 0.

                for neutral in event.neutralPfCandidates:
                    if neutral.particleId() in [4, 5]:
                        dr2 = deltaR2Obj(neutral, muon)
                        if dr2 < 0.4*0.4 and dr2 > 0.01*0.01: # veto cone
                            getWeights(neutral, event)
                            isoPt = neutral.pt()
                            if neutral.particleId() == 4:
                                muon.phIsoDr += isoPt * neutral.deltaBetaWeightDr
                                muon.phIsoDr2 += isoPt * neutral.deltaBetaWeightDr2
                                muon.phIsoDr201 += isoPt * neutral.deltaBetaWeightDr201
                                muon.phIsoPuppi += isoPt * neutral.deltaBetaWeightPuppi
                                muon.phIsoPt += isoPt * neutral.deltaBetaWeightPt
                                muon.phIsoLogPt += isoPt * neutral.deltaBetaWeightLogPt
                                muon.phIsoPuppi05 += isoPt * neutral.deltaBetaWeightPuppi05
                                muon.phIsoPuppi05Dr2 += isoPt * neutral.deltaBetaWeightPuppi05Dr2
                                muon.phIsoPuppi05Log2 += isoPt * neutral.deltaBetaWeightPuppi05Log2
                            else:
                                muon.nhIsoDr += isoPt * neutral.deltaBetaWeightDr
                                muon.nhIsoDr2 += isoPt * neutral.deltaBetaWeightDr2
                                muon.nhIsoDr201 += isoPt * neutral.deltaBetaWeightDr201
                                muon.nhIsoPuppi += isoPt * neutral.deltaBetaWeightPuppi
                                muon.nhIsoPt += isoPt * neutral.deltaBetaWeightPt
                                muon.nhIsoLogPt += isoPt * neutral.deltaBetaWeightLogPt
                                muon.nhIsoPuppi05 += isoPt * neutral.deltaBetaWeightPuppi05
                                muon.nhIsoPuppi05Dr2 += isoPt * neutral.deltaBetaWeightPuppi05Dr2
                                muon.nhIsoPuppi05Log2 += isoPt * neutral.deltaBetaWeightPuppi05Log2

            if len(event.muons) > 1:
                event.zmass_muon = (event.muons[0].p4() + event.muons[1].p4()).mass()

        return True

    def testMuonIso(self, muon, isocut ):
        '''dbeta corrected pf isolation with all charged particles instead of
        charged hadrons'''
        iso = muon.pfIsolationR04().sumChargedParticlePt + muon.pfIsolationR04().sumNeutralHadronEt +  muon.pfIsolationR04().sumPhotonEt
        iso += max(0., - muon.pfIsolationR04().sumPUPt + 0.5 *(muon.pfIsolationR04().sumNeutralHadronEt +  muon.pfIsolationR04().sumPhotonEt))
        muon.iso = iso
        muon.relIso = iso/muon.pt()
        isoNoNH = muon.pfIsolationR04().sumChargedParticlePt + muon.pfIsolationR04().sumPhotonEt
        isoNoNH += max(0., - muon.pfIsolationR04().sumPUPt + 0.5 *( muon.pfIsolationR04().sumPhotonEt))
        muon.isoNoNH = isoNoNH
        muon.relIsoNoNH = isoNoNH/muon.pt()
        return iso/muon.pt() < isocut

        # return muon.relIsoAllChargedDB05()<isocut
    def testMuonID(self, muon):
        '''Tight muon selection, no isolation requirement'''
        # import pdb; pdb.set_trace()
        return muon.isPFMuon() and \
        muon.isGlobalMuon() and \
               muon.globalTrack().normalizedChi2() < 10 and \
               muon.globalTrack().hitPattern().numberOfValidMuonHits() > 0 and \
               muon.numberOfMatchedStations()>1 and \
               muon.innerTrack().hitPattern().numberOfValidPixelHits()>0 and \
               muon.globalTrack().hitPattern().trackerLayersWithMeasurement() > 5 
