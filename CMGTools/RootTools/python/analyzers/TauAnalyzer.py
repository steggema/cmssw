from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Tau, GenParticle
from CMGTools.RootTools.utils.DeltaR import matchObjectCollection, deltaR
from CMGTools.RootTools.physicsobjects.genutils import allDaughters

def deltaRObj(obj1, obj2):
    return deltaR(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())

def p4sumvis(particles):
    particles = [p for p in particles if abs(p.pdgId()) not in [12, 14, 16]]
    p4 = particles[-1].p4() if particles else None
    particles.pop()
    for p in particles:
        p4 += p.p4()
    return p4

def finalDaughters(particle, daughters):
    '''Fills daughters with all the daughters of particle.
    Recursive function.'''
    if particle.numberOfDaughters() == 0:
        daughters.append(particle)
    else:
        for i in range( particle.numberOfDaughters() ):
            dau = GenParticle(particle.daughter(i))
            daughters = finalDaughters( dau, daughters )
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

        self.handles['leptons'] = AutoHandle(
            self.cfg_ana.tauSrc,
            'std::vector<reco::Tau>'
            )

    def process(self, iEvent, event):

        result = super(TauAnalyzer, self).process(iEvent, event) # This reads collections

        # import pdb; pdb.set_trace()

        taus = self.handles['leptons'].product()

        selTaus = [Tau(tau) for tau in taus if tau.pt() > self.cfg_ana.minPt and abs(tau.eta()) < self.cfg_ana.maxEta]

        event.taus = selTaus

        for tau in event.taus:
            tau.associatedVertex = event.goodVertices[0]

        if self.cfg_comp.isMC:
            # print event.eventId
            if not hasattr(event, 'genParticles'):
                genParticles = self.mchandles['genParticles'].product()
                event.genParticles = map( GenParticle, genParticles)
            # Allow Higgs/Z/photon/W
            allowedGenMothers = [15, 21, 23, 24, 25, 35, 36, 37]
            event.generatedTaus = [p for p in event.genParticles if abs(p.pdgId()) in self.cfg_ana.absGenIds and abs(p.mother().pdgId()) in allowedGenMothers]

            pairs = matchObjectCollection(event.taus, event.generatedTaus, self.cfg_ana.matchDeltaR)
            
            for tau in event.taus:
                genTau = pairs[tau]
                tau.genTau = genTau
                if genTau:
                    finDaughters = finalDaughters(tau.genTau, [])
                    tau.genVisP4 = p4sumvis(finDaughters)
                    tau.genMatchDeltaR = pairs[tau], deltaRObj(tau, pairs[tau])
                    genTau.recTau = tau

        return True
        
