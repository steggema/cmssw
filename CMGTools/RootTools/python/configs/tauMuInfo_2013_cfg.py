import CMGTools.RootTools.fwlite.Config as cfg
from CMGTools.RootTools.fwlite.Config import printComps

from CMGTools.RootTools.RootTools import * 



vertexAna = cfg.Analyzer(
    'VertexAnalyzer',
    goodVertices = 'goodPVFilter',
    vertexWeight = None,
    fixedWeight = 1,
    verbose = False,
    )


pileUpAna = cfg.Analyzer(
    'PileUpAnalyzer',
    true = True
    )


muonAna = cfg.Analyzer(
    'MuonAnalyzer',
    genSrc = 'genParticles',
    absGenIds = [13, 5],
    muonSrc = 'muons',
    minPt = 20.,
    maxEta = 2.4,
    matchDeltaR = 0.3,
    )

electronAna = cfg.Analyzer(
    'ElectronAnalyzer',
    genSrc = 'genParticles',
    absGenIds = [11, 5],
    electronSrc = 'gedGsfElectrons',
    minPt = 20.,
    maxEta = 2.4,
    matchDeltaR = 0.3,
    )

tauAna = cfg.Analyzer(
    'TauAnalyzer',
    genSrc = 'genParticles',
    absGenIds = [15, 5, 4, 3, 2, 1], # all jets from top may fake a tau
    tauSrc = 'hpsPFTauProducer',
    minPt = 20.,
    maxEta = 2.4,
    matchDeltaR = 0.5,
    )

genMetAna = cfg.Analyzer(
    'GenMetAnalyzer',
    genSrc = 'genParticles',
    metSrc = 'pfMet'
    )

# defined for vbfAna and eventSorter
vbfKwargs = dict( Mjj = 500,
                  deltaEta = 3.5    
                  )

jetAna = cfg.Analyzer(
    'PFJetAnalyzer',
    # jetCol = 'ak5PFJets',
    jetCol = 'ak5PFJetsCHS',
    jetPt = 30.,
    jetEta = 3.0,
    btagSFseed = 123456,
    relaxJetId = False, 
    jerCorr = False,
    #jesCorr = 1.,
    )

# jetAnaCHS = cfg.Analyzer(
#     'JetAnalyzer',
#     jetCol = 'cmgPFJetSelCHS',
#     jetPt = 20.,
#     jetEta = 4.7,
#     btagSFseed = 123456,
#     relaxJetId = False, 
#     jerCorr = False,
#     #jesCorr = 1.,
#     )


treeProducer = cfg.Analyzer(
    'H2TauTauTreeProducerPFStudies'
    )

metTreeProducer = cfg.Analyzer(
    'METTreeProducer'
    )
jetTreeProducer = cfg.Analyzer(
    'JetTreeProducer'
    )
electronTreeProducer = cfg.Analyzer(
    'ElectronTreeProducer'
    )
muonTreeProducer = cfg.Analyzer(
    'MuonTreeProducer'
    )
partonTreeProducer = cfg.Analyzer(
    'PartonTreeProducer'
    )
tauTreeProducer = cfg.Analyzer(
    'TauTreeProducer'
    )
#########################################################################################

from CMGTools.RootTools.pf_samples import allsamples

#########################################################################################

TT_FromSim = cfg.MCComponent(
    name = 'TT_FromSim',
    files = ['/afs/cern.ch/work/s/steggema/PF71X/reco.root'],
    xSection = 1. ,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 
    )

TT_FromGen = cfg.MCComponent(
    name = 'TT_FromGen',
    files = ['/afs/cern.ch/work/s/steggema/PF71X/reco_from_gen.root'],
    xSection = 1. ,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 
    )

TT_62X = cfg.MCComponent(
    name = 'TT_62X',
    files = ['/afs/cern.ch/work/s/steggema/TT62X_pythia8_AODIM.root'],
    xSection = 1. ,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 
    )

TT_62X_20PU = cfg.MCComponent(
    name = 'TT_62X_20PU',
    files = ['/afs/cern.ch/work/s/steggema/TT62X_pythia8_AODIM_20PU.root'],
    xSection = 1. ,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 
    )


TTRelVal_FromGen = cfg.MCComponent(
    name = 'TTRelVal_FromGen',
    files = ['/afs/cern.ch/work/s/steggema/PF71X/reco_from_gen_relval.root'],
    xSection = 1. ,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 
    )

sequence = cfg.Sequence( [
    # eventSelector,
    # triggerAna,
    vertexAna, 
    # dyJetsFakeAna,
    # NJetsAna,
    # higgsWeighter, 
    jetAna,
    muonAna,
    electronAna,
    tauAna,
    genMetAna,
    # treeProducer,
    metTreeProducer,
    jetTreeProducer,
    electronTreeProducer,
    muonTreeProducer,
    partonTreeProducer,
    tauTreeProducer
   ] )

# allsamples = [TT_62X_20PU, TTRelVal_FromGen, TT_FromGen, TT_62X, TT_FromSim] + allsamples
allsamples = allsamples

selectedComponents = allsamples


test = 0
if test==1:
    # comp = embed_Run2012C_22Jan
    # comp = DYJets
    # comp = HiggsGGH125
    # comp = HiggsSUSYGluGlu1000
    comp = allsamples[0]
    comp = allsamples[5]
    print comp
    # comp = data_Run2012A
    selectedComponents = [comp]
    comp.splitFactor = 1
    # comp.files = comp.files[:10]
    # comp.files = ['tauMu_fullsel_tree_CMG.root']
elif test==2:
    selectedComponents = allsamples[:12]
    for comp in selectedComponents:
        comp.splitFactor = 1
        comp.files = comp.files[:5]


config = cfg.Config( components = selectedComponents,
                     sequence = sequence )

printComps(config.components, True)
