from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import  GenParticle

from ROOT import Math

def p4suminvis(particles):
    particles = [p for p in particles if abs(p.pdgId()) in [12, 14, 16]]
    if not particles:
        return 0. # no neutrino
    p4 = particles[-1].p4()
    particles.pop()
    for p in particles:
        p4 += p.p4()
    return p4

def p4sumPV(particles):
    particles = [p for p in particles if p.fromPV()]
    p4 = Math.PtEtaPhiEVector()
    # p4 = particles[-1].p4() if particles else None
    # particles.pop()
    for p in particles:
        p4 += Math.PtEtaPhiEVector(p.pt(), p.eta(), p.phi(), p.mass())
    return p4

def p4sumPDG(particles, ids=[]):
    particles = [p for p in particles if abs(p.pdgId()) in ids or not ids]
    # p4 = particles[-1].p4() if particles else None
    # particles.pop()
    p4 = Math.PtEtaPhiEVector()
    for p in particles:
        p4 +=  Math.PtEtaPhiEVector(p.pt(), p.eta(), p.phi(), p.mass())
    return p4

class GenMetAnalyzer( Analyzer ):
    ''' Calculates gen MET based on final
    particles that are no neutrinos
    
    genSrc = 'genParticlesPruned',

    '''

    def declareHandles(self):

        super(GenMetAnalyzer, self).declareHandles()
        self.mchandles['genParticles'] =  AutoHandle(
            self.cfg_ana.genSrc,
            'std::vector<reco::GenParticle>'
            )
        self.handles['met'] =  AutoHandle(
            self.cfg_ana.metSrc,
            'std::vector<reco::PFMET>'
            )
        # self.handles['pfParticles'] =  AutoHandle(
        #     'particleFlow',
        #     'std::vector<reco::PFCandidate>'
        #     )


    def process(self, iEvent, event):

        super(GenMetAnalyzer, self).process(iEvent, event) # This reads collections

        # pfParticles = self.handles['pfParticles'].product()

        # p4CH = p4sumPDG(pfParticles, [211])
        # p4NH = p4sumPDG(pfParticles, [130])
        # p4PH = p4sumPDG(pfParticles, [22])

        # p4PV = p4sumPV(pfParticles)
        # p4all = p4sumPDG(pfParticles)

        # import pdb; pdb.set_trace()

        event.met = self.handles['met'].product()[0]

        if self.cfg_comp.isMC:
            # print event.eventId
            if not event.genParticles:
                genParticles = self.mchandles['genParticles'].product()
                event.genParticles = map( GenParticle, genParticles)
            # finalParticles = [p for p in event.genParticles if p.numberOfDaughters() == 0]
            finalParticles = [p for p in event.genParticles if p.status() in range(21, 30)]
            event.genMet = p4suminvis(finalParticles)


        return True
        
