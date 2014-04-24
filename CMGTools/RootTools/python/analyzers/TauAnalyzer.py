from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import PhysicsObject, GenParticle
from CMGTools.RootTools.utils.DeltaR import matchObjectCollection

from CMGTools.RootTools.analyzers.MuonAnalyzer import deltaRObj, isFinal, getWeights


def p4sumvis(particles):
    visparticles = [p for p in particles if abs(p.pdgId()) not in [12, 14, 16]]
    p4 = visparticles[-1].p4() if particles else 0.
    visparticles.pop()
    for p in visparticles:
        p4 += p.p4()
    return p4

def finalDaughters(particle, daughters):
    '''Fills daughters with all the daughters of particle.
    Recursive function.'''
    if particle.numberOfDaughters() == 0:
        daughters.append(particle)
    else:
        foundDaughter = False
        for i in range( particle.numberOfDaughters() ):
            dau = GenParticle(particle.daughter(i))
            if dau.status() >= 2:
                daughters = finalDaughters( dau, daughters )
                foundDaughter = True
        if not foundDaughter:
            daughters.append(particle)

    return daughters


class TauAnalyzer( Analyzer ):
    ''' Basic tau selection and matching to generated taus.
    Saves selected taus as event.taus
    
    genSrc = 'genParticlesPruned',
    absGenIds = [15],
    tauSrc = 'cmgTauSel',
    minPt = 20.,
    maxEta = 2.4,
    matchDeltaR = 0.3,
    '''

    def declareHandles(self):

        super(TauAnalyzer, self).declareHandles()
        self.mchandles['genParticles'] =  AutoHandle(
            self.cfg_ana.genSrc,
            'std::vector<reco::GenParticle>'
            )

        self.handles['taus'] = AutoHandle(
            self.cfg_ana.tauSrc,
            'std::vector<reco::PFTau>'
            )

        self.handles['decayModeFinding'] = AutoHandle(
            'hpsPFTauDiscriminationByDecayModeFinding',
            'reco::PFTauDiscriminator'
            )

        self.handles['chargedIso'] = AutoHandle(
            'hpsPFTauMVA3IsolationChargedIsoPtSum',
            'reco::PFTauDiscriminator'
            )

        self.handles['neutralIso'] = AutoHandle(
            'hpsPFTauMVA3IsolationNeutralIsoPtSum',
            'reco::PFTauDiscriminator'
            )

        self.handles['puIso'] = AutoHandle(
            'hpsPFTauMVA3IsolationPUcorrPtSum',
            'reco::PFTauDiscriminator'
            )



    def process(self, iEvent, event):

        super(TauAnalyzer, self).process(iEvent, event) # This reads collections

        # import pdb; pdb.set_trace()

        taus = self.handles['taus'].product()
        decayModeFinding = self.handles['decayModeFinding'].product()
        chargedIso = self.handles['chargedIso'].product()
        neutralIso = self.handles['neutralIso'].product()
        puIso = self.handles['puIso'].product()

        selTaus = [PhysicsObject(tau) for tau in taus if tau.pt() > self.cfg_ana.minPt and abs(tau.eta()) < self.cfg_ana.maxEta]

        event.taus = selTaus

        for tau in event.taus:
            tau.associatedVertex = event.goodVertices[0]

        for i in range(len(decayModeFinding)):
            dm = decayModeFinding.value(i)
            dmTau = decayModeFinding.key(i).get()

            chIso = chargedIso.value(i)
            chIsoTau = chargedIso.key(i).get()

            nIso = neutralIso.value(i)
            nIsoTau = neutralIso.key(i).get()

            pIso = puIso.value(i)
            pIsoTau = puIso.key(i).get()
            for tau in event.taus:
                if tau.physObj == dmTau:
                    tau.decayModeFinding = dm
                if tau.physObj == chIsoTau:
                    tau.chargedIso = chIso
                if tau.physObj == nIsoTau:
                    tau.neutralIso = nIso
                if tau.physObj == pIsoTau:
                    tau.puIso = pIso

        for tau in event.taus:
            if not hasattr(tau, 'chargedIso'):
                import pdb; pdb.set_trace()

        for tau in event.taus:
            photons = tau.isolationPFGammaCands()
            photonSum03 = 0.
            photonSumAll = 0.

            tau.phIsoDr = 0.
            tau.phIsoDr2 = 0.
            tau.phIsoDr201 = 0.
            tau.phIsoPt = 0.
            tau.phIsoLogPt = 0.
            tau.phIsoPuppi = 0.
            tau.phIsoPuppi05 = 0.
            tau.phIsoPuppi05Dr2 = 0.
            tau.phIsoPuppi05Log2 = 0.


            for photon in photons:
                photonSumAll += photon.get().pt()
                if photon.get().pt() > 0.3:
                    photonSum03 += photon.get().pt()
                isoPt = photon.pt()
                if photon.particleId() != 4:
                    # It's an electron, associate to PV with full weight
                    deltaZ = 100.
                    if photon.trackRef().get():
                        deltaZ = abs(event.goodVertices[0].z()-photon.trackRef().vertex().z())
                    elif photon.gsfTrackRef().get():
                        deltaZ = abs(event.goodVertices[0].z()-photon.gsfTrackRef().vertex().z())
                    if deltaZ < 2.:
                        tau.phIsoDr += isoPt
                        tau.phIsoDr2 += isoPt
                        tau.phIsoDr201 += isoPt
                        tau.phIsoPuppi += isoPt
                        tau.phIsoPt += isoPt
                        tau.phIsoLogPt += isoPt
                        tau.phIsoPuppi05 += isoPt
                        tau.phIsoPuppi05Dr2 += isoPt
                        tau.phIsoPuppi05Log2 += isoPt
                else:
                    getWeights(photon, event)
                    
                    tau.phIsoDr += isoPt * photon.deltaBetaWeightDr
                    tau.phIsoDr2 += isoPt * photon.deltaBetaWeightDr2
                    tau.phIsoDr201 += isoPt * photon.deltaBetaWeightDr201
                    tau.phIsoPuppi += isoPt * photon.deltaBetaWeightPuppi
                    tau.phIsoPt += isoPt * photon.deltaBetaWeightPt
                    tau.phIsoLogPt += isoPt * photon.deltaBetaWeightLogPt
                    tau.phIsoPuppi05 += isoPt * photon.deltaBetaWeightPuppi05
                    tau.phIsoPuppi05Dr2 += isoPt * photon.deltaBetaWeightPuppi05Dr2
                    tau.phIsoPuppi05Log2 += isoPt * photon.deltaBetaWeightPuppi05Log2

            
            tau.photonPtSum03 = photonSum03
            tau.photonPtSumAll = photonSumAll

        if self.cfg_comp.isMC:
            # print event.eventId
            if not hasattr(event, 'genParticlesTop'):
                genParticles = self.mchandles['genParticles'].product()
                event.genParticles = map( GenParticle, genParticles)
                genParticles = [p for p in event.genParticles if isFinal(p)]
            else:
                genParticles = event.genParticlesTop
            # Allow Higgs/Z/photon/W
            # allowedGenMothers = [15, 21, 23, 24, 25, 35, 36, 37]
            event.generatedTaus = [p for p in genParticles if abs(p.pdgId()) in self.cfg_ana.absGenIds]# and abs(p.mother().pdgId()) in allowedGenMothers]

            pairs = matchObjectCollection(event.taus, event.generatedTaus, self.cfg_ana.matchDeltaR)
            
            for tau in event.taus:
                genTau = pairs[tau]
                tau.parton = genTau
                if genTau and abs(genTau.pdgId()==15):
                    finDaughters = finalDaughters(tau.parton, [])
                    tau.genVisP4 = p4sumvis(finDaughters)
                    print genTau.pt(), tau.genVisP4.pt()
                    tau.genMatchDeltaR = pairs[tau], deltaRObj(tau, pairs[tau])
                    genTau.recTau = tau

        return True
        
