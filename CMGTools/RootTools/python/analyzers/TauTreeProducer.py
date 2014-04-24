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
        var( tr, 'tau_iso_photonPtSum03')
        var( tr, 'tau_iso_photonPtSumAll')

        var(tr, 'tau_iso_photonDr')
        var(tr, 'tau_iso_photonDr2')
        var(tr, 'tau_iso_photonDr201')
        var(tr, 'tau_iso_photonPuppi')
        var(tr, 'tau_iso_photonPuppi05')
        var(tr, 'tau_iso_photonPuppi05Dr2')
        var(tr, 'tau_iso_photonPuppi05Log2')
        var(tr, 'tau_iso_photonPtDr2')
        var(tr, 'tau_iso_photonLogPtDr2')

        var( tr, 'tau_decayModeFinding')
        var( tr, 'tau_decayMode')

        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_pdgId')
        var( tr, 'parton_ptvis')

        var(tr, 'ngoodvertices')
        var(tr, 'nvertices')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for tau in event.taus:
            tr.reset()
            fill( tr, 'ntaus',len(event.taus))

            fill(tr, 'ngoodvertices', len(event.goodVertices))
            fill(tr, 'nvertices', len(event.vertices))

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

            fill( tr, 'tau_iso_photonPtSum03', tau.photonPtSum03)
            fill( tr, 'tau_iso_photonPtSumAll', tau.photonPtSumAll)

            fill(tr, 'tau_iso_photonDr', tau.phIsoDr)
            fill(tr, 'tau_iso_photonDr2', tau.phIsoDr2)
            fill(tr, 'tau_iso_photonDr201', tau.phIsoDr201)
            fill(tr, 'tau_iso_photonPuppi', tau.phIsoPuppi)
            fill(tr, 'tau_iso_photonPuppi05', tau.phIsoPuppi05)
            fill(tr, 'tau_iso_photonPuppi05Dr2', tau.phIsoPuppi05Dr2)
            fill(tr, 'tau_iso_photonPuppi05Log2', tau.phIsoPuppi05Log2)
            fill(tr, 'tau_iso_photonPtDr2', tau.phIsoPt)
            fill(tr, 'tau_iso_photonLogPtDr2', tau.phIsoLogPt)

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
       
