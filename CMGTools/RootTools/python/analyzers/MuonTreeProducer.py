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
        var( tr, 'muon_charge')

        var( tr, 'muon_iso')
        var( tr, 'muon_reliso')
        var( tr, 'muon_isoNoNH')
        var( tr, 'muon_relisoNoNH')

        var( tr, 'muon_iso_chargedPt')
        var( tr, 'muon_iso_neutralPt')
        var( tr, 'muon_iso_photonEt')
        var( tr, 'muon_iso_sumPUPt')

        var(tr, 'muon_iso_photonDr')
        var(tr, 'muon_iso_photonDr2')
        var(tr, 'muon_iso_photonDr201')
        var(tr, 'muon_iso_photonPuppi')
        var(tr, 'muon_iso_photonPuppi05')
        var(tr, 'muon_iso_photonPuppi05Dr2')
        var(tr, 'muon_iso_photonPuppi05Log2')
        var(tr, 'muon_iso_photonPtDr2')
        var(tr, 'muon_iso_photonLogPtDr2')

        var(tr, 'muon_iso_neutralDr')
        var(tr, 'muon_iso_neutralDr2')
        var(tr, 'muon_iso_neutralDr201')
        var(tr, 'muon_iso_neutralPuppi')
        var(tr, 'muon_iso_neutralPuppi05')
        var(tr, 'muon_iso_neutralPuppi05Dr2')
        var(tr, 'muon_iso_neutralPuppi05Log2')
        var(tr, 'muon_iso_neutralPtDr2')
        var(tr, 'muon_iso_neutralLogPtDr2')

        var( tr, 'muon_id')

        var( tr, 'parton_pt')
        var( tr, 'parton_eta')
        var( tr, 'parton_pdgId')

        var(tr, 'ngoodvertices')
        var(tr, 'nvertices')

    def process(self, iEvent, event):
        
        tr = self.tree
        
        for muon in event.muons:
            tr.reset()
            fill( tr, 'nmuons',len(event.muons))

            fill(tr, 'ngoodvertices', len(event.goodVertices))
            fill(tr, 'nvertices', len(event.vertices))

            fill( tr, 'muon_pt', muon.pt())
            fill( tr, 'muon_eta', muon.eta())
            fill( tr, 'muon_charge', muon.charge())

            fill( tr, 'muon_iso', muon.iso)
            fill( tr, 'muon_reliso', muon.relIso)

            fill( tr, 'muon_isoNoNH', muon.isoNoNH)
            fill( tr, 'muon_relisoNoNH', muon.relIsoNoNH)

            fill( tr, 'muon_iso_chargedPt', muon.pfIsolationR04().sumChargedParticlePt)
            fill( tr, 'muon_iso_neutralPt', muon.pfIsolationR04().sumNeutralHadronEt)
            fill( tr, 'muon_iso_photonEt', muon.pfIsolationR04().sumPhotonEt)
            fill( tr, 'muon_iso_sumPUPt', muon.pfIsolationR04().sumPUPt)

            fill(tr, 'muon_iso_photonDr', muon.phIsoDr)
            fill(tr, 'muon_iso_photonDr2', muon.phIsoDr2)
            fill(tr, 'muon_iso_photonDr201', muon.phIsoDr201)
            fill(tr, 'muon_iso_photonPuppi', muon.phIsoPuppi)
            fill(tr, 'muon_iso_photonPuppi05', muon.phIsoPuppi05)
            fill(tr, 'muon_iso_photonPuppi05Dr2', muon.phIsoPuppi05Dr2)
            fill(tr, 'muon_iso_photonPuppi05Log2', muon.phIsoPuppi05Log2)
            fill(tr, 'muon_iso_photonPtDr2', muon.phIsoPt)
            fill(tr, 'muon_iso_photonLogPtDr2', muon.phIsoLogPt)
            fill(tr, 'muon_iso_neutralDr', muon.nhIsoDr)
            fill(tr, 'muon_iso_neutralDr2', muon.nhIsoDr2)
            fill(tr, 'muon_iso_neutralDr201', muon.nhIsoDr201)
            fill(tr, 'muon_iso_neutralPuppi', muon.nhIsoPuppi)
            fill(tr, 'muon_iso_neutralPuppi05', muon.nhIsoPuppi05)
            fill(tr, 'muon_iso_neutralPuppi05Dr2', muon.nhIsoPuppi05Dr2)
            fill(tr, 'muon_iso_neutralPuppi05Log2', muon.nhIsoPuppi05Log2)
            fill(tr, 'muon_iso_neutralPtDr2', muon.nhIsoPt)
            fill(tr, 'muon_iso_neutralLogPtDr2', muon.nhIsoLogPt)

            fill( tr, 'muon_id', muon.id)

            if hasattr(muon, 'parton') and muon.parton:
                fill(tr, 'parton_pt', muon.parton.pt())
                fill(tr, 'parton_eta', muon.parton.eta())
                fill(tr, 'parton_pdgId', muon.parton.pdgId())


            self.tree.tree.Fill()
       
