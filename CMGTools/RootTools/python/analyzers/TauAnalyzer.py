from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import PhysicsObject, GenParticle
from CMGTools.RootTools.utils.DeltaR import matchObjectCollection, deltaR

def isFinal(p):
    return not (p.numberOfDaughters() == 1 and p.daughter(0).pdgId() == p.pdgId())

def deltaRObj(obj1, obj2):
    return deltaR(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())

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
        
