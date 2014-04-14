from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class PartonTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var( tr, 'njets')
        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_phi')
        var( tr, 'parton_pdgId')
        var( tr, 'jet_pt')
        var( tr, 'jet_phi')
        var( tr, 'jet_eta')
        var( tr, 'jet_mass')
        var( tr, 'muon_pt')
        var( tr, 'muon_phi')
        var( tr, 'muon_eta')
        var( tr, 'muon_id')
        var( tr, 'electron_pt')
        var( tr, 'electron_phi')
        var( tr, 'electron_eta')

        var( tr, 'tau_pt')
        var( tr, 'tau_phi')
        var( tr, 'tau_eta')
        var( tr, 'tau_mass')
        var( tr, 'tau_decayModeFinding')


    def process(self, iEvent, event):
        
        tr = self.tree
        
        for p in event.genParticlesTop:
            tr.reset()
            fill(tr, 'parton_pt', p.pt())
            fill(tr, 'parton_eta', p.eta())
            fill(tr, 'parton_phi', p.phi())
            fill(tr, 'parton_pdgId', p.pdgId())
            if hasattr(p, 'jet'):
                fill( tr, 'jet_pt', p.jet.pt())
                fill( tr, 'jet_eta', p.jet.eta())
                fill( tr, 'jet_phi', p.jet.phi())
                fill( tr, 'jet_mass', p.jet.mass())
            if hasattr(p, 'recMuon'):
                fill( tr, 'muon_pt', p.recMuon.pt())
                fill( tr, 'muon_eta', p.recMuon.eta())
                fill( tr, 'muon_phi', p.recMuon.phi())
                fill( tr, 'muon_id', p.recMuon.id)
            if hasattr(p, 'recElectron'):
                fill( tr, 'electron_pt', p.recElectron.pt())
                fill( tr, 'electron_eta', p.recElectron.eta())
                fill( tr, 'electron_phi', p.recElectron.phi())
            if hasattr(p, 'recTau'):
                fill( tr, 'tau_pt', p.recTau.pt())
                fill( tr, 'tau_eta', p.recTau.eta())
                fill( tr, 'tau_phi', p.recTau.phi())
                fill( tr, 'tau_mass', p.recTau.mass())
                fill( tr, 'tau_decayModeFinding', p.recTau.decayModeFinding)

            self.tree.tree.Fill()
       
