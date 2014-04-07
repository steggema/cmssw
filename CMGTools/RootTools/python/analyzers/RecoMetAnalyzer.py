import random
from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Jet, GenJet, PhysicsObject
from CMGTools.RootTools.utils.DeltaR import cleanObjectCollection, matchObjectCollection
# from CMGTools.RootTools.physicsobjects.VBF import VBF
from CMGTools.RootTools.statistics.Counter import Counter, Counters
from CMGTools.RootTools.physicsobjects.BTagSF import BTagSF
from CMGTools.RootTools.physicsobjects.PhysicsObjects import GenParticle
from CMGTools.RootTools.utils.DeltaR import deltaR2
from CMGTools.RootTools.utils.cmsswRelease import isNewerThan

class RecoMetAnalyzer( Analyzer ):
    """Analyze jets ;-)

    This analyzer filters the jets that do not correspond to the leptons
    stored in event.selectedLeptons, and puts in the event:
    - jets: all jets passing the pt and eta cuts
    - cleanJets: the collection of clean jets
    - bJets: the bjets passing testBJet (see this method)

    Example configuration:

    jetAna = cfg.Analyzer(
      'RecoJetAnalyzer',
      # cmg jet input collection
      jetCol = 'cmgPFJetSel',
      # pt threshold
      jetPt = 30,
      # eta range definition
      jetEta = 5.0,
      # seed for the btag scale factor
      btagSFseed = 123456,
      # if True, the PF and PU jet ID are not applied, and the jets get flagged
      relaxJetId = False,
    )
    """

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(RecoMetAnalyzer,self).__init__(cfg_ana, cfg_comp, looperName)

    def declareHandles(self):
        super(RecoMetAnalyzer, self).declareHandles()


        self.handles['met'] = AutoHandle( self.cfg_ana.metCol,
                                           'std::vector<reco::PFMET>' )
        if self.cfg_comp.isMC:
            # and ("BB" in self.cfg_comp.name):
            self.mchandles['genParticles'] = AutoHandle( 'genParticles',
                                                         'std::vector<reco::GenParticle>' )
            self.mchandles['genMet'] = AutoHandle('genMetTrue',
                                                   'std::vector<reco::GenMET>')
##             self.mchandles['genJetsP'] = AutoHandle('ak5GenJetsNoNu',
##                                                     'std::vector< reco::GenJet >')

    def beginLoop(self):
        super(RecoMetAnalyzer,self).beginLoop()
        
    def process(self, iEvent, event):
        
        self.readCollections( iEvent )
        mets = self.handles['met'].product()

        event.met = PhysicsObject( mets[0] )

  
        event.genMET = None
        if self.cfg_comp.isMC:
             genMETs = map( PhysicsObject, self.mchandles['genMet'].product() ) 
             if genMETs:
                event.genMET = genMETs[0]
        



        genParticles = self.mchandles['genParticles'].product()
        event.genParticles = map( GenParticle, genParticles)
        event.partonMet = 0
        for gen in genParticles:
            if abs(gen.pdgId()) in [12, 14, 16]:
                if not event.partonMet:
                    event.partonMet = gen.p4()
                else:
                    event.partonMet += gen.p4()


        return True

    

