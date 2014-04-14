import ROOT
import os
import numpy
from officialStyle import officialStyle

ROOT.gROOT.SetBatch(True)
officialStyle(ROOT.gStyle)
ROOT.gStyle.SetOptTitle(0)
ROOT.gStyle.SetPadLeftMargin(0.18)
ROOT.gStyle.SetPadBottomMargin(0.15)

colours = [1, 2, 3, 4, 6, 7, 8]

directory = 'WithZFullProd'

OUTPUT_DIR = 'plots'
mode = 'top'
# mode = 'dy'

if mode == 'top':
    samples = ['TTNoTime', 'TTNoChi2', 'TT3SigmaNeighbour', 'TTChi2',  'TT3SigmaCut', 'TTTimeCutOnly']# 'TTTimeFromSeed',

    sampleDict = {
        'TTNoTime':{'name':'No time'},
        'TTTimeCutOnly':{'name':'Time cut only'},
        'TT3SigmaNeighbour':{'name':'3 sigma neighbour'},
        'TTChi2':{'name':'Chi2'},
        'TTNoChi2':{'name':'Default'},
        'TTTimeFromSeed':{'name':'Time from seed'},
        'TT3SigmaCut':{'name':'3 sigma cut'},

    }
else:
    OUTPUT_DIR = 'dyplots'
    samples = ['DYNoTime', 'DYTimeFromSeed', 'DY3SigmaNeighbour', 'DY3SigmaCut', 'DYTimeCutOnly']# 'TTTimeFromSeed',

    sampleDict = {
        'DYNoTime':{'name':'No time'},
        'DYTimeCutOnly':{'name':'Time cut only'},
        'DY3SigmaNeighbour':{'name':'3 sigma neighbour'},
        # 'TTChi2':{'name':'Chi2'},
        # 'TTNoChi2':{'name':'Default'},
        'DYTimeFromSeed':{'name':'Time from seed'},
        'DY3SigmaCut':{'name':'3 sigma cut'},

    }

setupsROC = {
    'muon':{
        'tree':'MuonTreeProducer/MuonTreeProducer_tree.root',
        'treename':'MuonTreeProducer',
    },
    'electron':{
        'tree':'ElectronTreeProducer/ElectronTreeProducer_tree.root',
        'treename':'ElectronTreeProducer',
    },
    'tau':{
        'tree':'TauTreeProducer/TauTreeProducer_tree.root',
        'treename':'TauTreeProducer',
    }
}

setupsVars = {
    'MET':{
        'tree':'METTreeProducer/METTreeProducer_tree.root',
        'treename':'METTreeProducer',
    },
    'jet':{
        'tree':'JetTreeProducer/JetTreeProducer_tree.root',
        'treename':'JetTreeProducer',
    },
    'tau':{
        'tree':'TauTreeProducer/TauTreeProducer_tree.root',
        'treename':'TauTreeProducer',
    },
    'electron':{
        'tree':'ElectronTreeProducer/ElectronTreeProducer_tree.root',
        'treename':'ElectronTreeProducer',
    }
}

setupsEff = {
    'tau':{
        'tree':'PartonTreeProducer/PartonTreeProducer_tree.root',
        'treename':'PartonTreeProducer',
    },
}

rocVars = {}
rocVars['muon'] = {
    'muon_relIso':{
        'var':'muon_reliso', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel}'
    },
    'muon_relIsoNoNH':{
        'var':'muon_relisoNoNH', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (charged, photon)'
    },
    'muon_relIsoOnlyPH':{
        'var':'muon_iso_photonEt/muon_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (photon)'
    },
    'muon_relIsoOnlyCH':{
        'var':'muon_iso_chargedPt/muon_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (charged)'
    },
    'muon_absIsoOnlyPH':{
        'var':'muon_iso_photonEt', 'nbinsx':200, 'xmin':0., 'xmax':40., 'title':'Iso (photon)'
    }
}
rocVars['electron'] = {
    'electron_relIso':{
        'var':'electron_reliso', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel}'
    },
    'electron_relIsoNoNH':{
        'var':'electron_relisoNoNH', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (charged, photon)'
    },
    'electron_relIsoOnlyPH':{
        'var':'electron_iso_photonEt/electron_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (photon)'
    },
    'electron_relIsoOnlyCH':{
        'var':'electron_iso_chargedPt/electron_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (charged)'
    },
    'electron_absIsoOnlyPH':{
        'var':'electron_iso_photonEt', 'nbinsx':200, 'xmin':0., 'xmax':40., 'title':'Iso (photon)'
    },
    'electron_mva':{
        'var':'-electron_mva', 'nbinsx':200, 'xmin':-1., 'xmax':1., 'title':'Electron MVA'
    }
}

rocVars['tau'] = {
    'tau_relIso':{
        'var':'(tau_iso_chargedPt+tau_iso_neutralPt)/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (no delta-beta)'
    },
    'tau_relIso_dbCorr06':{
        'var':'(tau_iso_chargedPt+TMath::Max(0., tau_iso_neutralPt-0.6*tau_iso_sumPUPt))/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} delta-beta 0.6'
    },
    'tau_relIso_dbCorr05':{
        'var':'(tau_iso_chargedPt+TMath::Max(0., tau_iso_neutralPt-0.5*tau_iso_sumPUPt))/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} delta-beta 0.5'
    },
    'tau_relIso_dbCorr04':{
        'var':'(tau_iso_chargedPt+TMath::Max(0., tau_iso_neutralPt-0.4*tau_iso_sumPUPt))/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} delta-beta 0.4'
    },
    'tau_relIso_dbCorr03':{
        'var':'(tau_iso_chargedPt+TMath::Max(0., tau_iso_neutralPt-0.3*tau_iso_sumPUPt))/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} delta-beta 0.3'
    },
    'tau_relIso_dbCorr02':{
        'var':'(tau_iso_chargedPt+TMath::Max(0., tau_iso_neutralPt-0.2*tau_iso_sumPUPt))/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} delta-beta 0.2'
    },
    'tau_relIsoOnlyPH':{
        'var':'tau_iso_neutralPt/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (photon)'
    },
    'tau_relIsoOnlyCH':{
        'var':'tau_iso_chargedPt/tau_pt', 'nbinsx':200, 'xmin':0., 'xmax':2., 'title':'Iso_{rel} (charged)'
    },
    'tau_absIsoOnlyPH':{
        'var':'tau_iso_neutralPt', 'nbinsx':200, 'xmin':0., 'xmax':40., 'title':'Iso (photon)'
    },
    'tau_absIso_dbCorr03':{
        'var':'(tau_iso_chargedPt+TMath::Max(0., tau_iso_neutralPt-0.3*tau_iso_sumPUPt))', 'nbinsx':200, 'xmin':0., 'xmax':20., 'title':'Iso (photon + charged) delta-beta 0.3'
    },
}

plotVariables = {}
plotVariables['MET'] = {
    'met_pt_res':{
        'var':'met_pt-partonMet_pt', 'nbinsx':50, 'xmin':-150., 'xmax':150., 'title':'#Delta E_{T}^{miss} (GeV)'
    },
    'met_px_res':{
        'var':'met_px-partonMet_px', 'nbinsx':50, 'xmin':-150., 'xmax':150., 'title':'#Delta E_{x}^{miss} (GeV)'
    },
    'met_py_res':{
        'var':'met_py-partonMet_py', 'nbinsx':50, 'xmin':-150., 'xmax':150., 'title':'#Delta E_{y}^{miss} (GeV)'
    },
    'wmass':{
        'var':'wmass', 'nbinsx':50, 'xmin':50., 'xmax':200., 'title':'W boson mass (GeV)'
    },
    'zmass_muon_res':{
        'var':'zmass_muon', 'nbinsx':50, 'xmin':50., 'xmax':125., 'title':'di-muon mass (GeV)'
    },
    'zmass_ele_res':{
        'var':'zmass_ele', 'nbinsx':50, 'xmin':50., 'xmax':125., 'title':'di-electron mass (GeV)'
    }
}

plotVariables['jet'] = {
    'jet_pt_res_parton':{
        'var':'jet_pt - parton_pt', 'nbinsx':50, 'xmin':-50., 'xmax':100., 'title':'#Delta jet p_{T} parton (GeV)'
    },
    'jet_pt_res_genjet':{
        'var':'jet_pt - genjet_pt', 'nbinsx':50, 'xmin':-50., 'xmax':100., 'title':'#Delta jet p_{T} gen. jet (GeV)'
    },
    'jet_mass_res_genjet':{
        'var':'jet_mass - genjet_mass', 'nbinsx':50, 'xmin':-20., 'xmax':50., 'title':'#Delta jet mass} gen. jet (GeV)'
    },
    'jet_phfraction':{
        'var':'jet_phfraction', 'nbinsx':30, 'xmin':0., 'xmax':1., 'title':'jet photon fraction'
    }
}

plotVariables['tau'] = {
    'tau_pt_res_parton':{
        'var':'tau_pt - parton_pt', 'nbinsx':20, 'xmin':-50., 'xmax':50., 'title':'#Delta #tau p_{T} parton (GeV)'
    },
    'tau_pt_res_partonvis':{
        'var':'tau_pt - parton_ptvis', 'nbinsx':20, 'xmin':-50., 'xmax':50., 'title':'#Delta #tau p_{T} parton (vis.) (GeV)'
    },
    'tau_mass':{
        'var':'tau_mass', 'nbinsx':50, 'xmin':0., 'xmax':4., 'title':'#tau mass (GeV)'
    },
    'tau_npizero':{
        'var':'tau_npizero', 'nbinsx':5, 'xmin':-0.5, 'xmax':4.5, 'title':'N #pi_{0}'
    },
    'tau_nchargedpf':{
        'var':'tau_nchargedpf', 'nbinsx':5, 'xmin':-0.5, 'xmax':4.5, 'title':'N charged PF'
    }
}

deltaPhi = 'TVector2::Phi_mpi_pi({phi1} - {phi2})'

plotVariables['electron'] = {
    'electron_pt_res_parton':{
        'var':'electron_pt - parton_pt', 'nbinsx':40, 'xmin':-10., 'xmax':10., 'title':'#Delta electron p_{T} parton (GeV)'
    },
    'electron_sigmaEtaEta':{
        'var':'electron_sigmaEtaEta', 'nbinsx':50, 'xmin':0., 'xmax':0.05, 'title':'#sigma_{#eta#eta}'
    },
    'electron_sigmaIphiIphi':{
        'var':'electron_sigmaIphiIphi', 'nbinsx':50, 'xmin':0., 'xmax':0.05, 'title':'#sigma_{i#phi i#phi}'
    },
    'electron_deltaEtaSuperClusterTrackAtVtx_res_100':{
        'var':'electron_deltaEtaSuperClusterTrackAtVtx*100', 'nbinsx':80, 'xmin':-2., 'xmax':3., 'title':'#Delta#eta supercluster track (*100)'
    },
    'electron_eta_res_parton_1000':{
        'var':'(electron_eta-parton_eta)*1000', 'nbinsx':50, 'xmin':-2., 'xmax':2., 'title':'electron #Delta#eta parton * 1000'
    },
    'electron_eSuperClusterOverP_res':{
        'var':'electron_eSuperClusterOverP', 'nbinsx':50, 'xmin':0.5, 'xmax':2.5, 'title':'electron E_{SC}/p'
    },
    'electron_mva':{
        'var':'electron_mva', 'nbinsx':50, 'xmin':0.5, 'xmax':1.000001, 'title':'electron MVA'
    },
    
    # 'electron_phi_res_parton':{
        # 'var':deltaPhi.format(phi1='electron_phi', phi2='parton_phi'), 'nbinsx':20, 'xmin':-0.15, 'xmax':0.15, 'title':'electron #Delta#phi parton'
    # }
}

effVariables = {}
effVariables['tau'] = {
    'tau_pt':{
        'var':'tau_pt', 'nbinsx':10, 'xmin':20., 'xmax':100., 'title':'#tau p_{T} (GeV)'
    },
    'parton_pt':{
        'var':'parton_pt', 'nbinsx':10, 'xmin':20., 'xmax':100., 'title':'gen. #tau p_{T} (GeV)'
    },
    'tau_eta':{
        'var':'tau_eta', 'nbinsx':10, 'xmin':-2.5, 'xmax':2.5, 'title':'#tau #eta'
    },
}
effVariables['electron'] = {
    'electron_pt':{
        'var':'electron_pt', 'nbinsx':10, 'xmin':20., 'xmax':100., 'title':'#electron p_{T} (GeV)'
    },
    'parton_pt':{
        'var':'parton_pt', 'nbinsx':10, 'xmin':20., 'xmax':100., 'title':'gen. electron p_{T} (GeV)'
    },
    'electron_eta':{
        'var':'electron_eta', 'nbinsx':10, 'xmin':-2.5, 'xmax':2.5, 'title':'electron #eta'
    },
}

selections = {}
selections['muon'] = {
    'signal':'abs(parton_pdgId)==13',
    # 'background':'abs(parton_pdgId)==5',
    'background':'abs(parton_pdgId)==5',
}
selections['electron'] = {
    'signal':'abs(parton_pdgId)==11',
    # 'background':'abs(parton_pdgId)==5',
    'background':'abs(parton_pdgId)==5',
}
selections['tau'] = {
    'signal':'abs(parton_pdgId)==15 && tau_decayModeFinding',
    # 'background':'abs(parton_pdgId)==5',
    'background':'abs(parton_pdgId)<=5 && tau_decayModeFinding',
}
selections['MET'] = {
    'signal':'1.'
}
selections['jet'] = {
    'signal':'jet_pt>30.'
}


effSelections = {}
effSelections['tau'] = {
    'decayModeFinding':{
        'signal':'abs(parton_pdgId)==15 && tau_decayModeFinding',
        'background':'abs(parton_pdgId)==15',
    },    
    'tauReco':{
        'signal':'abs(parton_pdgId)==15 && tau_pt>0.',
        'background':'abs(parton_pdgId)==15',
    },
}
effSelections['electron'] = {
    'electronReco':{
        'signal':'abs(parton_pdgId)==11 && electron_pt>0.',
        'background':'abs(parton_pdgId)==11',
    }
}

def ensureDir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)

def rocCurve(hS, hB):
  ''' Create a ROC TGraph from two input histograms.
  '''
  maxBin = hS.GetNbinsX()

  if hS.Integral() == 0.:
    print 'ROC curve creator, hist', hS.GetName(), 'has zero entries'
    return

  #rocPoints = [(hS.Integral(nBin, maxBin)/hS.Integral(), hB.Integral(nBin, maxBin)/hB.Integral()) for nBin in range(1, maxBin + 1) ]
  effsS = [hS.Integral(0, nBin)/hS.Integral(0, maxBin+1) for nBin in range(0, maxBin + 1) ]
  rejB = [1. - hB.Integral(0, nBin)/hB.Integral(0, maxBin+1) for nBin in range(0, maxBin + 1) ]

  rocCurve = ROOT.TGraph(maxBin, numpy.asarray(effsS), numpy.asarray(rejB))
  return rocCurve

def insertText(pad, insetText, lowX=0.18, lowY=0.77, colour=1, textSize=0.05):
    pad.cd()

    if not hasattr(pad, 'insertTexts'):
        pad.insertTexts = []
    insertText = ROOT.TPaveText(lowX, lowY, lowX+0.8, lowY+0.16, "NDC")
    insertText.SetBorderSize(0)
    insertText.SetFillStyle(0)
    insertText.SetTextAlign(12)
    insertText.SetTextSize(textSize)
    insertText.SetTextFont(42)
    insertText.SetTextColor (colour)
    insertText.AddText(insetText)
    insertText.Draw('same')
    pad.insertTexts.append(insertText)
    pad.Update()
    return insertText


def makePlotsROC(trees, setup):
    c = ROOT.TCanvas()
    variables = rocVars[setup]
    for var in variables:
        varDict = variables[var]
        hists = []
        rocGraphs = []
        for iTree, tree in enumerate(trees):
            histS = ROOT.TH1F(var+str(iTree)+'S', '', varDict['nbinsx'], varDict['xmin'], varDict['xmax'])
            histS.Sumw2()
            tree.Project(histS.GetName(), varDict['var'], selections[setup]['signal'])

            histB= ROOT.TH1F(var+str(iTree)+'B', '', varDict['nbinsx'], varDict['xmin'], varDict['xmax'])
            histB.Sumw2()
            tree.Project(histB.GetName(), varDict['var'], selections[setup]['background'])
            
            hists += [histS, histB]

            rocGraph = rocCurve(histS, histB)
            rocGraphs.append(rocGraph) # needs to be in scope

            rocGraph.SetLineColor(colours[iTree])
            rocGraph.SetLineStyle(iTree+1)
            rocGraph.SetMarkerColor(colours[iTree])
            rocGraph.GetYaxis().SetTitle('1-#varepsilon_{B} (bg. rejection)')
            rocGraph.GetXaxis().SetTitle('#varepsilon_{S} (signal efficiency)')

            # rocGraph.GetXaxis().SetRangeUser(0., 1.00001)
            if setup in ['muon', 'electron']:
                rocGraph.GetXaxis().SetRangeUser(0.5, 1.00001)
                rocGraph.GetYaxis().SetRangeUser(0.821, 1.)
            else:
                rocGraph.GetXaxis().SetRangeUser(0.3, 1.00001)
                rocGraph.GetYaxis().SetRangeUser(0.3, 1.)

            if iTree == 0:
                rocGraph.Draw('ACL')
            else:
                rocGraph.Draw('CL')
            sample = samples[iTree]

            insertText(c, sampleDict[sample]['name'], lowX=0.2, lowY=0.57 - 0.07*iTree, colour=colours[iTree])
        insertText(c, varDict['title'], lowX=0.2, lowY=0.88)

        # for iHist, hist in enumerate(hists):
        #     hist.SetLineColor(iHist+1)
        #     hist.SetLineStyle(iHist+1)
        #     hist.SetMarkerColor(iHist+1)
        #     if hist.Integral():
        #         hist.Scale(1./hist.Integral())
        #     hist.Draw('HIST E' + ('SAME' if iHist else ''))
        #     hist.GetYaxis().SetRangeUser(0., max(h.GetMaximum()/h.Integral() for h in hists if h.Integral()) * 1.1)
        #     insertText(c, fileLabels[iHist], 0.77 - 0.07*iHist, 0.67, iHist+1)
        # insertText(c, '{sel} selection'.format(sel=sel), 0.88, 0.15)
        ensureDir(OUTPUT_DIR+'/{setup}/'.format(setup=setup))
        c.Print(OUTPUT_DIR+'/{setup}/{title}.pdf'.format(setup=setup, title=var))
        c.Print(OUTPUT_DIR+'/{setup}/{title}.png'.format(setup=setup, title=var))


def makePlotsVars(trees, setup):
    c = ROOT.TCanvas()
    variables = plotVariables[setup]
    for var in variables:
        varDict = variables[var]
        hists = []
        for iTree, tree in enumerate(trees):
            hist = ROOT.TH1F(var+str(iTree), '', varDict['nbinsx'], varDict['xmin'], varDict['xmax'])
            hist.Sumw2()
            tree.Project(hist.GetName(), varDict['var'], selections[setup]['signal'])
           
            hists += [hist]
            if hist.Integral():
                hist.Scale(1./hist.Integral())

        maxY = max(hist.GetMaximum() for hist in hists)

        for iTree, hist in enumerate(hists):
            hist.SetLineColor(colours[iTree])
            hist.SetLineStyle(iTree+1)
            hist.SetMarkerColor(colours[iTree])

            hist.GetXaxis().SetTitle(varDict['title'])
            hist.GetYaxis().SetTitle('a.u.')

            hist.GetYaxis().SetRangeUser(0, maxY * 1.2)

            if iTree == 0:
                hist.Draw('HIST')
            else:
                hist.Draw('HIST SAME')
            sample = samples[iTree]

            if 'res' in var:
                mean = hist.GetMean()
                rms = hist.GetRMS()
                insertText(c, sampleDict[sample]['name'], lowX=0.63, lowY=0.81 - 0.08*iTree, colour=colours[iTree], textSize=0.03)
                insertText(c, 'm: {mean:.2f} RMS: {rms:.2f}'.format(mean=mean, rms=rms), lowX=0.63, lowY=0.77 - 0.08*iTree, colour=colours[iTree], textSize=0.03)
            else:
                insertText(c, sampleDict[sample]['name'], lowX=0.53, lowY=0.81 - 0.06*iTree, colour=colours[iTree])

        # for iHist, hist in enumerate(hists):
        #     hist.SetLineColor(iHist+1)
        #     hist.SetLineStyle(iHist+1)
        #     hist.SetMarkerColor(iHist+1)
        #     if hist.Integral():
        #         hist.Scale(1./hist.Integral())
        #     hist.Draw('HIST E' + ('SAME' if iHist else ''))
        #     hist.GetYaxis().SetRangeUser(0., max(h.GetMaximum()/h.Integral() for h in hists if h.Integral()) * 1.1)
        #     insertText(c, fileLabels[iHist], 0.77 - 0.07*iHist, 0.67, iHist+1)
        # insertText(c, '{sel} selection'.format(sel=sel), 0.88, 0.15)
        ensureDir(OUTPUT_DIR+'/{setup}/'.format(setup=setup))
        c.Print(OUTPUT_DIR+'/{setup}/{title}.pdf'.format(setup=setup, title=var))
        c.Print(OUTPUT_DIR+'/{setup}/{title}.png'.format(setup=setup, title=var))


def makePlotsEff(trees, setup):
    c = ROOT.TCanvas()
    variables = effVariables[setup]
    for var in variables:
        varDict = variables[var]
        for effSelection in effSelections[setup]:
            effSelectionDict = effSelections[setup][effSelection]
            hists = []
            for iTree, tree in enumerate(trees):
                histS = ROOT.TH1F(var+str(iTree)+'S', '', varDict['nbinsx'], varDict['xmin'], varDict['xmax'])
                histS.Sumw2()
                tree.Project(histS.GetName(), varDict['var'], effSelectionDict['signal'])

                histB = ROOT.TH1F(var+str(iTree)+'B', '', varDict['nbinsx'], varDict['xmin'], varDict['xmax'])
                histB.Sumw2()
                tree.Project(histB.GetName(), varDict['var'], effSelectionDict['background'])
                
                histS.Divide(histS, histB, 1., 1., 'B')
                hists.append(histS)

            # maxY = max(histS.GetMaximum() for hist in hists)

            for iTree, hist in enumerate(hists):
                hist.SetLineColor(colours[iTree])
                hist.SetLineStyle(iTree+1)
                hist.SetMarkerColor(colours[iTree])

                hist.GetXaxis().SetTitle(varDict['title'])
                hist.GetYaxis().SetTitle('Efficiency')

                hist.GetYaxis().SetRangeUser(0, 1.00001)

                if iTree == 0:
                    hist.Draw('HIST')
                else:
                    hist.Draw('HIST SAME')
                sample = samples[iTree]

                insertText(c, sampleDict[sample]['name'], lowX=0.53, lowY=0.81 - 0.06*iTree, colour=colours[iTree])
            insertText(c, varDict['title'], lowX=0.2, lowY=0.88)
            # for iHist, hist in enumerate(hists):
            #     hist.SetLineColor(iHist+1)
            #     hist.SetLineStyle(iHist+1)
            #     hist.SetMarkerColor(iHist+1)
            #     if hist.Integral():
            #         hist.Scale(1./hist.Integral())
            #     hist.Draw('HIST E' + ('SAME' if iHist else ''))
            #     hist.GetYaxis().SetRangeUser(0., max(h.GetMaximum()/h.Integral() for h in hists if h.Integral()) * 1.1)
            #     insertText(c, fileLabels[iHist], 0.77 - 0.07*iHist, 0.67, iHist+1)
            # insertText(c, '{sel} selection'.format(sel=sel), 0.88, 0.15)
            ensureDir(OUTPUT_DIR+'/{setup}/'.format(setup=setup))
            c.Print(OUTPUT_DIR+'/{setup}/{title}.pdf'.format(setup=setup, title=var+effSelection))
            c.Print(OUTPUT_DIR+'/{setup}/{title}.png'.format(setup=setup, title=var+effSelection))

if __name__ == '__main__':
    if mode == 'top':
        for setup in setupsROC:
            setupDict = setupsROC[setup]
            tfileNames = []
            for sample in samples:
                tfileNames.append('{dir}/{sample}/{treeLoc}'.format(dir=directory, sample=sample, treeLoc=setupDict['tree']))
            tfiles = [ROOT.TFile(tfileName) for tfileName in tfileNames]    
            trees = [tfile.Get(setupDict['treename']) for tfile in tfiles]
            makePlotsROC(trees, setup)

    for setup in setupsVars:
        setupDict = setupsVars[setup]
        tfileNames = []
        for sample in samples:
            tfileNames.append('{dir}/{sample}/{treeLoc}'.format(dir=directory, sample=sample, treeLoc=setupDict['tree']))
        tfiles = [ROOT.TFile(tfileName) for tfileName in tfileNames]    
        trees = [tfile.Get(setupDict['treename']) for tfile in tfiles]
        makePlotsVars(trees, setup)

    for setup in setupsEff:
        setupDict = setupsEff[setup]
        tfileNames = []
        for sample in samples:
            tfileNames.append('{dir}/{sample}/{treeLoc}'.format(dir=directory, sample=sample, treeLoc=setupDict['tree']))
        tfiles = [ROOT.TFile(tfileName) for tfileName in tfileNames]    
        trees = [tfile.Get(setupDict['treename']) for tfile in tfiles]
        makePlotsEff(trees, setup)        
