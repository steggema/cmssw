from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class MuonTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var( tr, 'nmuons')
        var( tr, 'muon_pt')
        var( tr, 'muon_eta')

        var( tr, 'muon_iso')
        var( tr, 'muon_reliso')
        var( tr, 'muon_id')

        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_pdgId')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for muon in event.muons:
            tr.reset()
            fill( tr, 'nmuons',len(event.muons))
            fill( tr, 'muon_pt', muon.pt())
            fill( tr, 'muon_eta', muon.eta())

            fill( tr, 'muon_iso', muon.iso)
            fill( tr, 'muon_reliso', muon.relIso)
            fill( tr, 'muon_id', muon.id)

            if hasattr(muon, 'parton') and muon.parton:
                fill(tr, 'parton_pt', muon.parton.pt())
                fill(tr, 'parton_eta', muon.parton.eta())
                fill(tr, 'parton_pdgId', muon.parton.pdgId())


            self.tree.tree.Fill()
       
