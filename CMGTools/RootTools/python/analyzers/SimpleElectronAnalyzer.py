from CMGTools.RootTools.fwlite.Analyzer import Analyzer
from CMGTools.RootTools.fwlite.AutoHandle import AutoHandle


from CMGTools.RootTools.physicsobjects.HTauTauElectron import HTauTauElectron as Electron
from CMGTools.RootTools.utils.DeltaR import deltaR2

class SimpleElectronAnalyzer( Analyzer ):
    ''' Basic electron selection and Z boson reconstruction.
    Saves selected electrons as event.electrons
    
    electronSrc = 'cmgElectronSel',
    minPt = 20.,
    maxEta = 2.4,
    '''

    def declareHandles(self):

        super(SimpleElectronAnalyzer, self).declareHandles()

        self.handles['leptons'] = AutoHandle(
            self.cfg_ana.electronSrc,
            'std::vector<reco::GsfElectron>'
            )

        self.handles['jets'] = AutoHandle(
            ('ak5PFJets', '', 'RERECO'),
            'std::vector<reco::PFJet>'
            )     

        self.handles['MET'] = AutoHandle(
            ('pfMet', '', 'RERECO'),
            'std::vector<reco::PFMET>'
            )      

    def process(self, iEvent, event):

        super(SimpleElectronAnalyzer, self).process(iEvent, event) # This reads collections

        # import pdb; pdb.set_trace()

        electrons = self.handles['leptons'].product()

        selElectrons = [Electron(electron) for electron in electrons if electron.p4(0).pt() > self.cfg_ana.minPt and abs(electron.eta()) < self.cfg_ana.maxEta]

        event.electrons = selElectrons

        event.met = self.handles['MET'].product()[0]
        jets = self.handles['jets'].product()
        event.jets = [jet for jet in jets if jet.pt()>20. and abs(jet.eta()) < 4.7]
        
        for electron in event.electrons:
            electron.iso = electron.pfIsolationVariables().sumChargedHadronPt + max(electron.pfIsolationVariables().sumNeutralHadronEt+electron.pfIsolationVariables().sumPhotonEt-0.5 * electron.pfIsolationVariables().sumPUPt, 0.)
            electron.isoNoNH = electron.pfIsolationVariables().sumChargedHadronPt + max(electron.pfIsolationVariables().sumPhotonEt-0.5 * electron.pfIsolationVariables().sumPUPt, 0.)
                # relIso = muon.pfIsolationR04().sumChargedHadronPt+muon.pfIsolationR04().sumPhotonEt
            electron.relIso = electron.iso/electron.p4(0).pt()
            electron.relIsoNoNH = electron.isoNoNH/electron.p4(0).pt()

        event.cleanJets = [jet for jet in event.jets]
        if len(event.electrons) > 0:
            event.cleanJets = [jet for jet in event.jets if deltaR2(jet.eta(), jet.phi(), event.electrons[0].eta(), event.electrons[0].phi()) < 0.25]
        if len(event.electrons) > 1:
            event.cleanJets = [jet for jet in event.cleanJets if deltaR2(jet.eta(), jet.phi(), event.electrons[1].eta(), event.electrons[1].phi()) < 0.25]
        
        for i, jet in enumerate(event.cleanJets):
            if i == 0:
                event.hadronicP4 = jet.p4()
            else:
                event.hadronicP4 += jet.p4()


        if len(event.electrons) > 1:
            event.zmass_ele = (event.electrons[0].p4() + event.electrons[1].p4()).mass()
            event.zboson_pt = (event.electrons[0].p4() + event.electrons[1].p4()).pt()
            if hasattr(event, 'hadronicP4'):
                event.jetZPt = (event.electrons[0].p4() + event.electrons[1].p4() + event.hadronicP4).pt()
            event.zmass_ele2 = (event.electrons[0].p4(1) + event.electrons[1].p4(1)).mass()
            event.zboson2_pt = (event.electrons[0].p4(1) + event.electrons[1].p4(1)).pt()


        # print 'Electrons', len(event.electrons)
        return True
        
