from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class TauTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var( tr, 'ntaus')
        var( tr, 'tau_pt')
        var( tr, 'tau_eta')
        var( tr, 'tau_mass')
        var( tr, 'tau_charge')

        var( tr, 'tau_ncharged')
        var( tr, 'tau_nchargedpf')
        var( tr, 'tau_npizero')

        var( tr, 'tau_iso_chargedPt')
        var( tr, 'tau_iso_neutralPt')
        var( tr, 'tau_iso_sumPUPt')

        var( tr, 'tau_decayModeFinding')
        var( tr, 'tau_decayMode')

        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_pdgId')
        var( tr, 'parton_ptvis')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for tau in event.taus:
            tr.reset()
            fill( tr, 'ntaus',len(event.taus))
            fill( tr, 'tau_pt', tau.pt())
            fill( tr, 'tau_eta', tau.eta())
            fill( tr, 'tau_mass', tau.mass())
            fill( tr, 'tau_charge', tau.charge())

            fill( tr, 'tau_ncharged', tau.signalTauChargedHadronCandidates().size())
            fill( tr, 'tau_nchargedpf', tau.signalPFChargedHadrCands().size())
            fill( tr, 'tau_npizero', tau.signalPiZeroCandidates().size())

            fill( tr, 'tau_iso_chargedPt', tau.chargedIso)
            fill( tr, 'tau_iso_neutralPt', tau.neutralIso)
            fill( tr, 'tau_iso_sumPUPt', tau.puIso)

            fill( tr, 'tau_decayModeFinding', tau.decayModeFinding)
            fill( tr, 'tau_decayMode', tau.decayMode())

            if hasattr(tau, 'parton') and tau.parton:
                fill(tr, 'parton_pt', tau.parton.pt())
                fill(tr, 'parton_eta', tau.parton.eta())
                fill(tr, 'parton_pdgId', tau.parton.pdgId())
                # fill(tr, 'parton_vispt', tau.genVisP4)
            if hasattr(tau, 'genVisP4'):
                if tau.genVisP4:
                    fill(tr, 'parton_ptvis', tau.genVisP4.pt())


            self.tree.tree.Fill()
       
