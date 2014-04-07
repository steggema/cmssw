from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Muon, GenParticle
from CMGTools.RootTools.utils.DeltaR import matchObjectCollection, deltaR

def deltaRObj(obj1, obj2):
    return deltaR(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())

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

    def process(self, iEvent, event):

        super(MuonAnalyzer, self).process(iEvent, event) # This reads collections

        # import pdb; pdb.set_trace()

        muons = self.handles['leptons'].product()

        selMuons = [Muon(muon) for muon in muons if muon.pt() > self.cfg_ana.minPt and abs(muon.eta()) < self.cfg_ana.maxEta and muon.isTrackerMuon()]

        event.muons = selMuons

        for muon in event.muons:
            muon.associatedVertex = event.goodVertices[0]

        if self.cfg_comp.isMC:
            # print event.eventId
            
            if not hasattr(event, 'genParticles'):
                genParticles = self.mchandles['genParticles'].product()
                event.genParticles = map( GenParticle, genParticles)
            # Allow Higgs/Z/photon/W
            # allowedGenMothers = [6, 15, 21, 23, 24, 25, 35, 36, 37]
            allowedGenMothers = [6, 15, 24]
            event.generatedMuons = [p for p in event.genParticles if abs(p.pdgId()) in self.cfg_ana.absGenIds and abs(p.mother().pdgId()) in allowedGenMothers and p.status() in range(21, 30)]

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

        return True

    def testMuonIso(self, muon, isocut ):
        '''dbeta corrected pf isolation with all charged particles instead of
        charged hadrons'''
        iso = muon.pfIsolationR04().sumChargedParticlePt + muon.pfIsolationR04().sumNeutralHadronEt +  muon.pfIsolationR04().sumPhotonEt
        iso += max(0., - muon.pfIsolationR04().sumPUPt + 0.5 *(muon.pfIsolationR04().sumNeutralHadronEt +  muon.pfIsolationR04().sumPhotonEt))
        muon.iso = iso
        muon.relIso = iso/muon.pt()
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
