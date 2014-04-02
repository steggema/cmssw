from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import GenParticle
from CMGTools.RootTools.utils.DeltaR import matchObjectCollection, deltaR
from CMGTools.RootTools.physicsobjects.HTauTauElectron import HTauTauElectron as Electron

def deltaRObj(obj1, obj2):
    return deltaR(obj1.eta(), obj1.phi(), obj2.eta(), obj2.phi())

class ElectronAnalyzer( Analyzer ):
    ''' Basic electron selection and matching to generated electrons.
    Saves selected electrons as event.electrons
    
    genSrc = 'genParticlesPruned',
    absGenIds = [11],
    electronSrc = 'cmgElectronSel',
    minPt = 20.,
    maxEta = 2.4,
    matchDeltaR = 0.3,
    '''

    def declareHandles(self):

        super(ElectronAnalyzer, self).declareHandles()
        self.mchandles['genParticles'] =  AutoHandle(
            self.cfg_ana.genSrc,
            'std::vector<reco::GenParticle>'
            )

        self.handles['leptons'] = AutoHandle(
            self.cfg_ana.electronSrc,
            'std::vector<cmg::Electron>'
            )

    def process(self, iEvent, event):

        result = super(ElectronAnalyzer, self).process(iEvent, event) # This reads collections

        # import pdb; pdb.set_trace()

        electrons = self.handles['leptons'].product()

        selElectrons = [Electron(electron) for electron in electrons if electron.pt() > self.cfg_ana.minPt and abs(electron.eta()) < self.cfg_ana.maxEta]

        event.electrons = selElectrons

        for electron in event.electrons:
            electron.associatedVertex = event.goodVertices[0]

        if self.cfg_comp.isMC:
            # print event.eventId
            if not hasattr(event, 'genParticles'):
                genParticles = self.mchandles['genParticles'].product()
                event.genParticles = map( GenParticle, genParticles)
            # Allow Higgs/Z/photon/W
            allowedGenMothers = [15, 21, 23, 24, 25, 35, 36, 37]
            event.generatedElectrons = [p for p in event.genParticles if abs(p.pdgId()) in self.cfg_ana.absGenIds and abs(p.mother().pdgId()) in allowedGenMothers]

            pairs = matchObjectCollection(event.electrons, event.generatedElectrons, self.cfg_ana.matchDeltaR)
            
            for electron in event.electrons:
                genEle = pairs[electron]
                electron.genElectron = genEle
                if genEle:
                    electron.genMatchDeltaR = genEle, deltaRObj(electron, genEle)
                    genEle.recElectron = electron
                # print electron.pt(), pairs[electron]
        # print 'Electrons', len(event.electrons)
        return True
        
