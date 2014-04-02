import random
from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Jet, GenJet
from CMGTools.RootTools.utils.DeltaR import cleanObjectCollection, matchObjectCollection
from CMGTools.RootTools.statistics.Counter import Counter, Counters
from CMGTools.RootTools.physicsobjects.PhysicsObjects import GenParticle
from CMGTools.RootTools.utils.DeltaR import deltaR2
from CMGTools.RootTools.utils.cmsswRelease import isNewerThan

class PFJetAnalyzer( Analyzer ):
    """Analyze jets ;-)

    This analyzer filters the jets that do not correspond to the leptons
    stored in event.selectedLeptons, and puts in the event:
    - jets: all jets passing the pt and eta cuts
    - cleanJets: the collection of clean jets

    Example configuration:

    jetAna = cfg.Analyzer(
      'PFJetAnalyzer',
      # cmg jet input collection
      jetCol = 'cmgPFJetSel',
      # pt threshold
      jetPt = 30,
      # eta range definition
      jetEta = 5.0,
    )
    """

    def __init__(self, cfg_ana, cfg_comp, looperName):
        super(PFJetAnalyzer,self).__init__(cfg_ana, cfg_comp, looperName)
        self.is2012 = isNewerThan('CMSSW_5_2_0')

    def declareHandles(self):
        super(PFJetAnalyzer, self).declareHandles()

        self.handles['jets'] = AutoHandle( self.cfg_ana.jetCol,
                                           'std::vector<cmg::PFJet>' )
        if self.cfg_comp.isMC:
            self.mchandles['genParticles'] = AutoHandle( 'genParticlesPruned',
                                                         'std::vector<reco::GenParticle>' )
            self.mchandles['genJets'] = AutoHandle('genJetSel',
                                                   'std::vector<cmg::PhysicsObjectWithPtr< edm::Ptr<reco::GenJet> > >')

    def beginLoop(self):
        super(PFJetAnalyzer,self).beginLoop()
        self.counters.addCounter('jets')
        count = self.counters.counter('jets')
        count.register('all events')
        count.register('at least 2 good jets')
        count.register('at least 2 clean jets')
        
    def process(self, iEvent, event):
        
        self.readCollections( iEvent )
        cmgJets = self.handles['jets'].product()

        allJets = []
        event.jets = []
        event.cleanJets = []

        leptons = []
        if hasattr(event, 'muons'):
            leptons += event.muons
        if hasattr(event, 'electrons'):
            leptons += event.electrons


        genJets = None
        

        for cmgJet in cmgJets:
            jet = Jet( cmgJet )
            allJets.append( jet )

            if self.testJet( jet ):
                event.jets.append(jet)
        

        if self.cfg_comp.isMC:
            genJets = map( GenJet, self.mchandles['genJets'].product() ) 
            # Use DeltaR = 0.25 matching like JetMET
            pairs = matchObjectCollection( event.jets, genJets, 0.25*0.25)
            for jet in event.jets:
                jet.genJet = pairs[jet]

        self.counters.counter('jets').inc('all events')

        event.cleanJets, dummy = cleanObjectCollection( event.jets,
                                                        masks = leptons,
                                                        deltaRMin = 0.5 )
        


        # associating a jet to each lepton
        pairs = matchObjectCollection( leptons, allJets, 0.5*0.5)
        for lepton in leptons:
            jet = pairs[lepton]
            lepton.jet = jet
                

        # associating a leg to each clean jet
        invpairs = matchObjectCollection( event.cleanJets, leptons, 99999. )
        for jet in event.cleanJets:
            leg = invpairs[jet]
            jet.leg = leg

        for jet in event.cleanJets:
            jet.matchGenParton = 999.0

        if self.cfg_comp.isMC and "BB" in self.cfg_comp.name:
            genParticles = self.mchandles['genParticles'].product()
            event.genParticles = map( GenParticle, genParticles)
            for gen in genParticles:
                if abs(gen.pdgId())==5 and gen.mother() and abs(gen.mother().pdgId())==21:
                    for jet in event.cleanJets:
                        dR=deltaR2(jet.eta(), jet.phi(), gen.eta(), gen.phi() )
                        if dR<jet.matchGenParton:
                            jet.matchGenParton=dR

        
        if len( event.jets )>=2:
            self.counters.counter('jets').inc('at least 2 good jets')
               
        if len( event.cleanJets )>=2:
            self.counters.counter('jets').inc('at least 2 clean jets')

        return True



    def testJetID(self, jet):
        jet.puJetIdPassed = jet.puJetId(wp53x=True)
        jet.pfJetIdPassed = jet.looseJetId()
        return True # Only save the info
        # return jet.puJetIdPassed and jet.pfJetIdPassed
        
        
    def testJet( self, jet ):
        return jet.pt() > self.cfg_ana.jetPt and \
               abs( jet.eta() ) < self.cfg_ana.jetEta and \
               self.testJetID(jet)

