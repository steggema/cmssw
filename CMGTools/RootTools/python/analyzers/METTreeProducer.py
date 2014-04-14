from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class METTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var( tr, 'met_pt')
        var( tr, 'met_px')
        var( tr, 'met_py')
        var( tr, 'met_phi')
        var( tr, 'genMet_pt')
        var( tr, 'genMet_px')
        var( tr, 'genMet_py')
        var( tr, 'genMet_phi')
        var( tr, 'partonMet_pt')
        var( tr, 'partonMet_px')
        var( tr, 'partonMet_py')
        var( tr, 'partonMet_phi')

        var( tr, 'wmass')
        var( tr, 'zmass_muon')
        var( tr, 'zmass_ele')
        var(tr, 'njets')

    def process(self, iEvent, event):
        
        tr = self.tree
        tr.reset()

        fill( tr, 'njets', len(event.cleanJets))
        fill( tr, 'met_pt', event.met.pt())
        fill( tr, 'met_px', event.met.px())
        fill( tr, 'met_py', event.met.py())
        fill( tr, 'met_phi', event.met.phi())
        if hasattr(event, 'genMET'):
            fill( tr, 'genMet_pt', event.genMET.pt())
            fill( tr, 'genMet_px', event.genMET.px())
            fill( tr, 'genMet_py', event.genMET.py())
            fill( tr, 'genMet_phi', event.genMET.eta())
        if hasattr(event, 'partonMet'):
            if event.partonMet:
                fill(tr, 'partonMet_pt', event.partonMet.pt())
                fill(tr, 'partonMet_px', event.partonMet.px())
                fill(tr, 'partonMet_py', event.partonMet.py())
                fill(tr, 'partonMet_phi', event.partonMet.phi())

        if hasattr(event, 'wmass'):
            fill(tr, 'wmass', event.wmass)
        if hasattr(event, 'zmass_muon'):
            fill( tr, 'zmass_muon', event.zmass_muon)
        if hasattr(event, 'zmass_ele'):
            fill( tr, 'zmass_ele', event.zmass_ele)
        

        self.tree.tree.Fill()
       
