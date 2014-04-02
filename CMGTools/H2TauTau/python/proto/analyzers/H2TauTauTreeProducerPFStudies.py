from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy
from CMGTools.H2TauTau.proto.analyzers.ntuple import *
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle


class H2TauTauTreeProducerPFStudies( TreeAnalyzerNumpy ):
    '''Tree producer for the H->tau tau analysis.

    Some of the functions in this class should be made available to everybody.'''
    
    def declareVariables(self):

        tr = self.tree

        var( tr, 'run', int)
        var( tr, 'lumi', int)
        var( tr, 'evt', int)
        var( tr, 'rho')

        var( tr, 'pfmet')
        var( tr, 'pfmetphi')
        var( tr, 'genmet')
        var( tr, 'genmetphi')

        bookTau(tr, 'tau1')
        bookGenParticle(tr, 'genTau1')
        bookTau(tr, 'tau2')
        bookGenParticle(tr, 'genTau2')

        bookMuon(tr, 'muon1')
        bookGenParticle(tr, 'genMuon1')
        bookMuon(tr, 'muon2')
        bookGenParticle(tr, 'genMuon2')

        bookEle(tr, 'electron1')
        bookGenParticle(tr, 'genElectron1')
        bookEle(tr, 'electron2')
        bookGenParticle(tr, 'genElectron2')

        bookJet(tr, 'jet1')
        bookGenParticle(tr, 'genJet1')
        bookJet(tr, 'jet2')
        bookGenParticle(tr, 'genJet2')

        var( tr, 'nJets')

        var( tr, 'weight')
        var( tr, 'vertexWeight')
        var( tr, 'nVert')
       
       
    def declareHandles(self):
        super(H2TauTauTreeProducerPFStudies, self).declareHandles()
        self.handles['pfmetraw'] = AutoHandle(
            'cmgPFMETRaw',
            'std::vector<cmg::BaseMET>' 
            )
        
    def process(self, iEvent, event):
        self.readCollections( iEvent )
        
        run = iEvent.eventAuxiliary().id().run()
        lumi = iEvent.eventAuxiliary().id().luminosityBlock()
        eventId = iEvent.eventAuxiliary().id().event()

        event.run = run
        event.lumi = lumi
        event.eventId = eventId

        tr = self.tree
        tr.reset()

        fill( tr, 'run', event.run) 
        fill( tr, 'lumi',event.lumi)
        fill( tr, 'evt', event.eventId)
        fill( tr, 'rho', event.rho)

        # import pdb; pdb.set_trace()
        pfmet = self.handles['pfmetraw'].product()[0]
        fill(tr, 'pfmet', pfmet.pt())
        fill(tr, 'pfmetphi', pfmet.phi())

        if hasattr(event, 'genMet'):
            fill(tr, 'genmet', event.genMet.pt())
            fill(tr, 'genmetphi', event.genMet.phi())

        for i, genTau in enumerate(event.generatedTaus):
            if i > 1: break
            fillGenParticle(tr, 'genTau{i}'.format(i=i+1), event.generatedTaus[i])
            if hasattr(event.generatedTaus[i], 'recTau'):
                fillTau(tr, 'tau{i}'.format(i=i+1), event.generatedTaus[i].recTau)

        for i, genMuon in enumerate(event.generatedMuons):
            if i > 1: break
            fillGenParticle(tr, 'genMuon{i}'.format(i=i+1), event.generatedMuons[i])
            if hasattr(event.generatedMuons[i], 'recMuon'):
                fillMuon(tr, 'muon{i}'.format(i=i+1), event.generatedMuons[i].recMuon)

        for i, genElectron in enumerate(event.generatedElectrons):
            if i > 1: break
            fillGenParticle(tr, 'genElectron{i}'.format(i=i+1), event.generatedElectrons[i])
            if hasattr(event.generatedElectrons[i], 'recElectron'):
                fillEle(tr, 'electron{i}'.format(i=i+1), event.generatedElectrons[i].recElectron)




        # if len(event.taus) > 0:
        #     fillTau(tr, 'tau1', event.taus[0])
        #     if event.taus[0].genTau:
        #         fillGenParticle(tr, 'genTau1', event.taus[0].genTau)
        # if len(event.taus) > 1:
        #     fillTau(tr, 'tau2', event.taus[1])
        #     if event.taus[1].genTau:
        #         fillGenParticle(tr, 'genTau2', event.taus[1].genTau)

        # if len(event.muons) > 0:
        #     fillMuon(tr, 'muon1', event.muons[0])
        #     if event.muons[0].genMuon:
        #         fillGenParticle(tr, 'genMuon1', event.muons[0].genMuon)
        # if len(event.muons) > 1:
        #     fillMuon(tr, 'muon2', event.muons[1])
        #     if event.muons[1].genMuon:
        #         fillGenParticle(tr, 'genMuon2', event.muons[1].genMuon)

        # if len(event.electrons) > 0:
        #     fillEle(tr, 'electron1', event.electrons[0])
        #     if event.electrons[0].genElectron:
        #         fillGenParticle(tr, 'genElectron1', event.electrons[0].genElectron)
        # if len(event.electrons) > 1:
        #     fillEle(tr, 'electron2', event.electrons[1])
        #     if event.electrons[1].genElectron:
        #         fillGenParticle(tr, 'genElectron2', event.electrons[1].genElectron)

        nJets = len(event.jets)
        fill(tr, 'nJets', nJets )

        if nJets>=1:
           fillJet(tr, 'jet1', event.jets[0] )
           if event.jets[0].genJet:
                fillGenParticle(tr, 'genJet1', event.jets[0].genJet ) 
        if nJets>=2:
           fillJet(tr, 'jet2', event.jets[1] )
           if event.jets[1].genJet:
                fillGenParticle(tr, 'genJet2', event.jets[1].genJet ) 



        fill(tr, 'weight', event.eventWeight)


        if hasattr( event, 'vertexWeight'): 
            fill(tr, 'vertexWeight', event.vertexWeight)
            fill(tr, 'nVert', len(event.vertices) ) 


        self.tree.tree.Fill()
