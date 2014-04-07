import random
from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle
from CMGTools.RootTools.physicsobjects.PhysicsObjects import Jet, GenJet
from CMGTools.RootTools.utils.DeltaR import cleanObjectCollection, matchObjectCollection
# from CMGTools.RootTools.physicsobjects.VBF import VBF
from CMGTools.RootTools.statistics.Counter import Counter, Counters
from CMGTools.RootTools.physicsobjects.BTagSF import BTagSF
from CMGTools.RootTools.physicsobjects.PhysicsObjects import GenParticle
from CMGTools.RootTools.utils.DeltaR import deltaR2
from CMGTools.RootTools.utils.cmsswRelease import isNewerThan

class RecoJetAnalyzer( Analyzer ):
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
        super(RecoJetAnalyzer,self).__init__(cfg_ana, cfg_comp, looperName)

    def declareHandles(self):
        super(RecoJetAnalyzer, self).declareHandles()

        if hasattr(self.cfg_ana, 'jetColType'):
            self.handles['jets'] = AutoHandle( self.cfg_ana.jetCol,
                                           self.cfg_ana.jetColType )
        else:
            self.handles['jets'] = AutoHandle( self.cfg_ana.jetCol,
                                           'std::vector<reco::PFJet>' )
        if self.cfg_comp.isMC:
            # and ("BB" in self.cfg_comp.name):
            self.mchandles['genParticles'] = AutoHandle( 'genParticles',
                                                         'std::vector<reco::GenParticle>' )
            self.mchandles['genJets'] = AutoHandle('ak5GenJets',
                                                   'std::vector<reco::GenJet>')
##             self.mchandles['genJetsP'] = AutoHandle('ak5GenJetsNoNu',
##                                                     'std::vector< reco::GenJet >')

    def beginLoop(self):
        super(RecoJetAnalyzer,self).beginLoop()
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
        event.bJets = []
        event.cleanJets = []
        event.cleanBJets = []

        leptons = []
        if hasattr(event, 'selectedLeptons'):
            leptons = event.selectedLeptons


        genJets = None
        if self.cfg_comp.isMC:
            genJets = map( GenJet, self.mchandles['genJets'].product() ) 
            
        for cmgJet in cmgJets:
            jet = Jet( cmgJet )
            allJets.append( jet )
            if self.cfg_comp.isMC and hasattr( self.cfg_comp, 'jetScale'):
                scale = random.gauss( self.cfg_comp.jetScale,
                                      self.cfg_comp.jetSmear )
                jet.scaleEnergy( scale )
            if genJets:
                # Use DeltaR = 0.25 matching like JetMET
                pairs = matchObjectCollection( [jet], genJets, 0.25*0.25)
                if pairs[jet] is None:
                    pass
                    #jet.genJet = None
                else:
                    jet.genJet = pairs[jet] 
                    # print 'Matched genjet'
                    # print 'jet pt', jet.pt(), 'parton pt', jet.genJet.pt()
                # print jet, jet.genJet

            #Add JER correction for MC jets. Requires gen-jet matching. 
            if self.cfg_comp.isMC and hasattr(self.cfg_ana, 'jerCorr') and self.cfg_ana.jerCorr:
                self.jerCorrection(jet)

            #Add JES correction for MC jets.
            if self.cfg_comp.isMC and hasattr(self.cfg_ana, 'jesCorr'):
                self.jesCorrection(jet, self.cfg_ana.jesCorr)
            if self.testJet( jet ):
                event.jets.append(jet)

                
        self.counters.counter('jets').inc('all events')

        event.cleanJets, dummy = cleanObjectCollection( event.jets,
                                                        masks = leptons,
                                                        deltaRMin = 0.5 )
        

        event.cleanBJets, dummy = cleanObjectCollection( event.bJets,
                                                         masks = leptons,
                                                         deltaRMin = 0.5 )        

        pairs = matchObjectCollection( leptons, allJets, 0.5*0.5)


        # associating a jet to each lepton
        for lepton in leptons:
            jet = pairs[lepton]
            if jet is None:
                lepton.jet = lepton
            else:
                lepton.jet = jet

        # associating a leg to each clean jet
        invpairs = matchObjectCollection( event.cleanJets, leptons, 99999. )
        for jet in event.cleanJets:
            leg = invpairs[jet]
            jet.leg = leg

        for jet in event.cleanJets:
            jet.matchGenParton=999.0
            jet.leptonDr=999.0

        genParticles = self.mchandles['genParticles'].product()
        event.genParticles = map( GenParticle, genParticles)
        event.cleanJets = []
        for gen in genParticles:
            if abs(gen.pdgId()) in [1,2,3,4,5,21]:
                for jet in event.jets:
                    dR=deltaR2(jet.eta(), jet.phi(), gen.eta(), gen.phi() )
                    if dR<jet.matchGenParton and dR < 0.25:
                        jet.matchGenParton=dR
                        jet.parton=gen
                        # print 'Matched parton'
                        # print 'jet pt', jet.pt(), 'parton pt', gen.pt()
            # if abs(gen.pdgId()) in [11, 12, 13, 14, 15, 16] and gen.status() in [2, 3]:
            if abs(gen.pdgId()) in [12, 13, 14, 15, 16]:
                print gen.pdgId(),
                for jet in event.jets:
                    dR=deltaR2(jet.eta(), jet.phi(), gen.eta(), gen.phi() )
                    if dR<jet.leptonDr and dR < 0.5: # it's in the cone
                        print 'Reject JET'
                        jet.leptonDr = dR
                
        for jet in event.jets:
            if jet.leptonDr > 0.5:
                event.cleanJets.append(jet)


        event.jets30 = [jet for jet in event.jets if jet.pt()>30]
        event.cleanJets30 = [jet for jet in event.cleanJets if jet.pt()>30]
        
        if len( event.jets30 )>=2:
            self.counters.counter('jets').inc('at least 2 good jets')
               
        if len( event.cleanJets30 )>=2:
            self.counters.counter('jets').inc('at least 2 clean jets')

        if len(event.cleanJets)<2:
            return True

        return True

    

    def jesCorrection(self, jet, scale=0.):
        ''' Adds JES correction in number of sigmas (scale)
        '''
        # Do nothing if nothing to change
        if scale == 0.:
            return

        unc = jet.uncOnFourVectorScale()

        totalScale = 1. + scale * unc

        if totalScale < 0.:
            totalScale = 0.
        jet.scaleEnergy(totalScale)

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

