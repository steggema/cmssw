from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class ElectronTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var( tr, 'nelectrons')
        var( tr, 'electron_pt')
        var( tr, 'electron_pt1')
        var( tr, 'electron_pt2')
        var( tr, 'electron_eta')
        var( tr, 'electron_phi')

        var( tr, 'electron_iso')
        var( tr, 'electron_reliso')
        var( tr, 'electron_isoNoNH')
        var( tr, 'electron_relisoNoNH')
        var( tr, 'electron_iso_chargedPt')
        var( tr, 'electron_iso_neutralPt')
        var( tr, 'electron_iso_photonEt')
        var( tr, 'electron_iso_sumPUPt')
     
        var( tr, 'electron_eSuperClusterOverP')
        var( tr, 'electron_sigmaEtaEta')
        var( tr, 'electron_sigmaIphiIphi')
        
        var( tr, 'electron_hcalOverEcal')
        var( tr, 'electron_deltaEtaSuperClusterTrackAtVtx')
        var( tr, 'electron_mva')

        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_phi')
        var( tr, 'parton_pdgId')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for electron in event.electrons:
            tr.reset()
            fill( tr, 'nelectrons',len(event.electrons))
            fill( tr, 'electron_pt', electron.p4(0).pt()) #P4_FROM_SUPER_CLUSTER=0, P4_COMBINATION=1, P4_PFLOW_COMBINATION=2
            fill( tr, 'electron_pt1', electron.p4(1).pt())
            fill( tr, 'electron_pt2', electron.p4(2).pt())
            fill( tr, 'electron_eta', electron.p4(0).eta())
            fill( tr, 'electron_eta', electron.p4(0).phi())

            fill( tr, 'electron_iso', electron.iso)
            fill( tr, 'electron_reliso', electron.relIso)

            fill( tr, 'electron_isoNoNH', electron.isoNoNH)
            fill( tr, 'electron_relisoNoNH', electron.relIsoNoNH)

            fill( tr, 'electron_iso_chargedPt', electron.pfIsolationVariables().sumChargedHadronPt)
            fill( tr, 'electron_iso_neutralPt', electron.pfIsolationVariables().sumNeutralHadronEt)
            fill( tr, 'electron_iso_photonEt', electron.pfIsolationVariables().sumPhotonEt)
            fill( tr, 'electron_iso_sumPUPt', electron.pfIsolationVariables().sumPUPt)

            fill( tr, 'electron_eSuperClusterOverP', electron.eSuperClusterOverP())
            fill( tr, 'electron_sigmaEtaEta', electron.sigmaEtaEta())
            fill( tr, 'electron_sigmaIphiIphi', electron.sigmaIphiIphi())
            fill( tr, 'electron_hcalOverEcal', electron.hcalOverEcal())
            fill( tr, 'electron_deltaEtaSuperClusterTrackAtVtx', electron.deltaEtaSuperClusterTrackAtVtx())
            fill( tr, 'electron_mva', electron.mva())
            
            
             

            if hasattr(electron, 'parton') and electron.parton:
                fill(tr, 'parton_pt', electron.parton.pt())
                fill(tr, 'parton_eta', electron.parton.eta())
                fill(tr, 'parton_eta', electron.parton.phi())
                fill(tr, 'parton_pdgId', electron.parton.pdgId())


            self.tree.tree.Fill()
       
