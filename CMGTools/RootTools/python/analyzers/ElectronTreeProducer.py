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

        var(tr, 'electron_iso_photonDr')
        var(tr, 'electron_iso_photonDr2')
        var(tr, 'electron_iso_photonDr201')
        var(tr, 'electron_iso_photonPuppi')
        var(tr, 'electron_iso_photonPuppi05')
        var(tr, 'electron_iso_photonPuppi05Dr2')
        var(tr, 'electron_iso_photonPuppi05Log2')
        var(tr, 'electron_iso_photonPtDr2')
        var(tr, 'electron_iso_photonLogPtDr2')
        var(tr, 'electron_iso_neutralDr')
        var(tr, 'electron_iso_neutralDr2')
        var(tr, 'electron_iso_neutralDr201')
        var(tr, 'electron_iso_neutralPuppi')
        var(tr, 'electron_iso_neutralPuppi05')
        var(tr, 'electron_iso_neutralPuppi05Dr2')
        var(tr, 'electron_iso_neutralPuppi05Log2')
        var(tr, 'electron_iso_neutralPtDr2')
        var(tr, 'electron_iso_neutralLogPtDr2')
     
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

        var(tr, 'ngoodvertices')
        var(tr, 'nvertices')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for electron in event.electrons:
            tr.reset()
            fill( tr, 'nelectrons',len(event.electrons))

            fill(tr, 'ngoodvertices', len(event.goodVertices))
            fill(tr, 'nvertices', len(event.vertices))

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


            fill(tr, 'electron_iso_photonDr', electron.phIsoDr)
            fill(tr, 'electron_iso_photonDr2', electron.phIsoDr2)
            fill(tr, 'electron_iso_photonDr201', electron.phIsoDr201)
            fill(tr, 'electron_iso_photonPuppi', electron.phIsoPuppi)
            fill(tr, 'electron_iso_photonPuppi05', electron.phIsoPuppi05)
            fill(tr, 'electron_iso_photonPuppi05Dr2', electron.phIsoPuppi05Dr2)
            fill(tr, 'electron_iso_photonPuppi05Log2', electron.phIsoPuppi05Log2)
            fill(tr, 'electron_iso_photonPtDr2', electron.phIsoPt)
            fill(tr, 'electron_iso_photonLogPtDr2', electron.phIsoLogPt)
            fill(tr, 'electron_iso_neutralDr', electron.nhIsoDr)
            fill(tr, 'electron_iso_neutralDr2', electron.nhIsoDr2)
            fill(tr, 'electron_iso_neutralDr201', electron.nhIsoDr201)
            fill(tr, 'electron_iso_neutralPuppi', electron.nhIsoPuppi)
            fill(tr, 'electron_iso_neutralPuppi05', electron.nhIsoPuppi05)
            fill(tr, 'electron_iso_neutralPuppi05Dr2', electron.nhIsoPuppi05Dr2)
            fill(tr, 'electron_iso_neutralPuppi05Log2', electron.nhIsoPuppi05Log2)
            fill(tr, 'electron_iso_neutralPtDr2', electron.nhIsoPt)
            fill(tr, 'electron_iso_neutralLogPtDr2', electron.nhIsoLogPt)

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
       
