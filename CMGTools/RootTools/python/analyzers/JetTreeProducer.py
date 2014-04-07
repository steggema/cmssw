from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class JetTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var( tr, 'njets')
        var( tr, 'jet_pt')
        var( tr, 'jet_clean')
        var( tr, 'jet_eta')
        var( tr, 'jet_mass')
        var( tr, 'genjet_pt')
        var( tr, 'genjet_eta')
        var( tr, 'genjet_mass')
        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_pdgId')
        var( tr, 'jet_phfraction')
        var( tr, 'jet_nhfraction')
        var( tr, 'jet_chfraction')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for jet in event.jets:
            tr.reset()
            fill( tr, 'njets',len(event.cleanJets))
            fill( tr, 'jet_pt', jet.pt())
            fill( tr, 'jet_clean', (jet in event.cleanJets))
            fill( tr, 'jet_eta', jet.eta())
            fill( tr, 'jet_mass', jet.mass())
            if hasattr(jet, 'genJet') and jet.genJet:
                fill( tr, 'genjet_pt', jet.genJet.pt())
                fill( tr, 'genjet_eta', jet.genJet.eta())
                fill( tr, 'genjet_mass', jet.genJet.mass())
            if hasattr(jet, 'parton') and jet.parton:
                fill(tr, 'parton_pt', jet.parton.pt())
                fill(tr, 'parton_eta', jet.parton.eta())
                fill(tr, 'parton_pdgId', jet.parton.pdgId())
            fill( tr, 'jet_phfraction', jet.photonEnergyFraction())
            fill( tr, 'jet_nhfraction', jet.neutralHadronEnergyFraction())
            fill( tr, 'jet_chfraction', jet.chargedHadronEnergyFraction())

            self.tree.tree.Fill()
       
