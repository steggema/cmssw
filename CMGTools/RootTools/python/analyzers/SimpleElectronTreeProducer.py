from CMGTools.RootTools.analyzers.TreeAnalyzerNumpy import TreeAnalyzerNumpy

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )


class SimpleElectronTreeProducer( TreeAnalyzerNumpy ):
    def declareVariables(self):
        tr = self.tree
        var(tr, 'nelectrons')
        var(tr, 'zmass')
        var(tr, 'zpt')
        var(tr, 'zmass2')
        var(tr, 'zpt2')
        var(tr, 'met')
        var(tr, 'njets')
        var(tr, 'jetZPt')
        var(tr, 'jetVectorPt')

        for iEle in range(0, 3):
            var( tr, 'electron{i}_pt'.format(i=iEle))
            var( tr, 'electron{i}_pt1'.format(i=iEle))
            var( tr, 'electron{i}_pt2'.format(i=iEle))
            var( tr, 'electron{i}_eta'.format(i=iEle))
            var( tr, 'electron{i}_phi'.format(i=iEle))
            var( tr, 'electron{i}_iso'.format(i=iEle))
            var( tr, 'electron{i}_reliso'.format(i=iEle))
            var( tr, 'electron{i}_isoNoNH'.format(i=iEle))
            var( tr, 'electron{i}_relisoNoNH'.format(i=iEle))
            var( tr, 'electron{i}_iso_chargedPt'.format(i=iEle))
            var( tr, 'electron{i}_iso_neutralPt'.format(i=iEle))
            var( tr, 'electron{i}_iso_photonEt'.format(i=iEle))
            var( tr, 'electron{i}_iso_sumPUPt'.format(i=iEle))

            var( tr, 'electron{i}_eSuperClusterOverP'.format(i=iEle))
            var( tr, 'electron{i}_sigmaEtaEta'.format(i=iEle))
            var( tr, 'electron{i}_sigmaIphiIphi'.format(i=iEle))
            
            var( tr, 'electron{i}_hcalOverEcal'.format(i=iEle))
            var( tr, 'electron{i}_deltaEtaSuperClusterTrackAtVtx'.format(i=iEle))
            var( tr, 'electron{i}_mva'.format(i=iEle))

    def process(self, iEvent, event):
        
        tr = self.tree
        tr.reset()
        fill( tr, 'nelectrons', len(event.electrons))
        if hasattr(event, 'zmass_ele'):
            fill(tr, 'zmass', event.zmass_ele)
            fill(tr, 'zpt', event.zboson_pt)
            fill(tr, 'zmass2', event.zmass_ele2)
            fill(tr, 'zpt2', event.zboson2_pt)

        fill(tr, 'met', event.met.pt())
        fill(tr, 'njets', len(event.cleanJets))

        if hasattr(event, 'jetZPt'):
            fill(tr, 'jetZPt', event.jetZPt)
        if hasattr(event, 'hadronicP4'):
            fill(tr, 'jetVectorPt', event.hadronicP4.pt())

        for iEle, electron in enumerate(event.electrons):
            if iEle >= 3:
                break
            fill( tr, 'electron{i}_pt'.format(i=iEle), electron.p4(0).pt()) #P4_FROM_SUPER_CLUSTER=0, P4_COMBINATION=1, P4_PFLOW_COMBINATION=2
            fill( tr, 'electron{i}_pt1'.format(i=iEle), electron.p4(1).pt())
            fill( tr, 'electron{i}_pt2'.format(i=iEle), electron.p4(2).pt())
            fill( tr, 'electron{i}_eta'.format(i=iEle), electron.p4(0).eta())
            fill( tr, 'electron{i}_phi'.format(i=iEle), electron.p4(0).phi())

            fill( tr, 'electron{i}_iso'.format(i=iEle), electron.iso)
            fill( tr, 'electron{i}_reliso'.format(i=iEle), electron.relIso)

            fill( tr, 'electron{i}_isoNoNH'.format(i=iEle), electron.isoNoNH)
            fill( tr, 'electron{i}_relisoNoNH'.format(i=iEle), electron.relIsoNoNH)

            fill( tr, 'electron{i}_iso_chargedPt'.format(i=iEle), electron.pfIsolationVariables().sumChargedHadronPt)
            fill( tr, 'electron{i}_iso_neutralPt'.format(i=iEle), electron.pfIsolationVariables().sumNeutralHadronEt)
            fill( tr, 'electron{i}_iso_photonEt'.format(i=iEle), electron.pfIsolationVariables().sumPhotonEt)
            fill( tr, 'electron{i}_iso_sumPUPt'.format(i=iEle), electron.pfIsolationVariables().sumPUPt)

            fill( tr, 'electron{i}_eSuperClusterOverP'.format(i=iEle), electron.eSuperClusterOverP())
            fill( tr, 'electron{i}_sigmaEtaEta'.format(i=iEle), electron.sigmaEtaEta())
            fill( tr, 'electron{i}_sigmaIphiIphi'.format(i=iEle), electron.sigmaIphiIphi())
            fill( tr, 'electron{i}_hcalOverEcal'.format(i=iEle), electron.hcalOverEcal())
            fill( tr, 'electron{i}_deltaEtaSuperClusterTrackAtVtx'.format(i=iEle), electron.deltaEtaSuperClusterTrackAtVtx())
            fill( tr, 'electron{i}_mva'.format(i=iEle), electron.mva())
            
            


            self.tree.tree.Fill()
       
