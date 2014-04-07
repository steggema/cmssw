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
                                           'std::vector<reco::PFJet>' )
        if self.cfg_comp.isMC:
            self.mchandles['genParticles'] = AutoHandle( 'genParticles',
                                                         'std::vector<reco::GenParticle>' )
            self.mchandles['genJets'] = AutoHandle('ak5GenJets',
                                                   'std::vector<reco::GenJet>')

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
        # if hasattr(event, 'muons'):
        #     leptons += event.muons
        # if hasattr(event, 'electrons'):
        #     leptons += event.electrons

        genJets = None

        for cmgJet in cmgJets:
            jet = Jet( cmgJet )
            allJets.append( jet )

            if self.testJet( jet ):
                event.jets.append(jet)
        

        genJets = map( GenJet, self.mchandles['genJets'].product() ) 
        # Use DeltaR = 0.25 matching like JetMET
        pairs = matchObjectCollection( event.jets, genJets, 0.25*0.25)
        for jet in event.jets:
            jet.genJet = pairs[jet]
            if pairs[jet] and not hasattr(pairs[jet], 'jet'):
                pairs[jet].jet = jet

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



        genParticles = self.mchandles['genParticles'].product()
        event.genParticles = map( GenParticle, genParticles)

        # event.genParticlesTop = [p.getFinal() for p in event.genParticles if p.mother() and 
            # (abs(p.mother().pdgId())==6 or abs(p.mother().pdgId())==24) and p.status() in range(21, 30) and abs(p.pdgId()) != 24]

        topsWs = [p.getFinal() for p in event.genParticles if abs(p.pdgId()) in [6, 24] and p.status() in range(21, 30)]

        event.genParticlesTop = []
        for top in topsWs:
            for i in range(top.numberOfDaughters()):
                 d = GenParticle(top.daughter(i))
                 event.genParticlesTop.append(d.getFinal())



        print 'Top and W daughters', [(p.pdgId(), p.eta(), p.status()) for p in event.genParticlesTop]
        # print 'Top stati', [p.getFinal().status() for p in event.genParticles if abs(p.pdgId()) == 6 and p.status() in range(21, 30)]

        for gen in event.genParticlesTop:
            for jet in event.cleanJets:
                dR = deltaR2(jet.eta(), jet.phi(), gen.eta(), gen.phi() )
                if dR<0.25:
                    if hasattr(jet, 'matchGenPartons'):
                        jet.matchGenPartons.append(gen) # for merging
                    else:
                        jet.parton = gen
                        jet.matchGenPartons = [gen]
                    if not hasattr(gen, 'jet'):
                        gen.jet = jet
                        print 'Attaching jet', gen.pdgId(), gen.status()

        for gen in event.genParticles:
            if abs(gen.pdgId()) == 24 and gen.numberOfDaughters() == 2:  # W boson
                jets = []
                for i in range(0, gen.numberOfDaughters()):
                    daughter = GenParticle(gen.daughter(i)).getFinal()
                    print 'Found W daughter', daughter.pdgId(), daughter.status()
                    if abs(daughter.pdgId()) in [1,2,3,4,5] and hasattr(daughter, 'jet'):
                        jets.append(daughter.jet)
                        print 'Has a matched jet'
                if len(jets) == 2:
                    event.wmass = (jets[0].p4() + jets[1].p4()).mass()
                    print 'Reconstruct W', event.wmass


        
        if len( event.jets )>=2:
            self.counters.counter('jets').inc('at least 2 good jets')
               
        if len( event.cleanJets )>=2:
            self.counters.counter('jets').inc('at least 2 clean jets')

        return True




    def testJetID(self, jet):
        try:
            jet.puJetIdPassed = jet.puJetId(wp53x=True)
            jet.pfJetIdPassed = jet.jetID("POG_PFID_Loose")

            if self.cfg_ana.relaxJetId:
                return True
            else:
                return jet.puJetIdPassed and jet.pfJetIdPassed
        except AttributeError:
            if jet.nConstituents() > 1 and jet.photonEnergyFraction() < 0.99 and jet.neutralHadronEnergyFraction() < 0.99 and (abs(jet.eta()) > 2.4 or (jet.chargedHadronMultiplicity() > 0 and jet.chargedHadronEnergyFraction() > 0. and jet.electronEnergyFraction() < 0.99)):
                return True
            return False
        
        
    def testJet( self, jet ):
        # 2 is loose pile-up jet id
        return jet.pt() > self.cfg_ana.jetPt and \
               abs( jet.eta() ) < self.cfg_ana.jetEta and \
               self.testJetID(jet)
               # jet.passPuJetId('full', 2)


    def testBJet(self, jet):
        # medium csv working point
        # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagPerformanceOP#B_tagging_Operating_Points_for_3
        try:
            jet.btagMVA = jet.btag("combinedSecondaryVertexBJetTags")
            jet.btagFlag = self.btagSF.BTagSFcalc.isbtagged(jet.pt(), 
                              jet.eta(),
                              jet.btag("combinedSecondaryVertexBJetTags"),
                              abs(jet.partonFlavour()),
                              not self.cfg_comp.isMC,
                              0,0,
                              self.is2012 )
            return jet.pt()>20 and \
                   abs( jet.eta() ) < 2.4 and \
                   jet.btagFlag and \
                   self.testJetID(jet)
        except AttributeError:
            return False

