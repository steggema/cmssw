import itertools
from CMGTools.RootTools.fwlite.Config import printComps
from CMGTools.RootTools.utils.connect import connect
import CMGTools.RootTools.fwlite.Config as cfg

aliases = {
    '/TT_Tune4C_13TeV.*REALNOTIME.*':'TTNoTime',
    '/TT_Tune4C_13TeV.*_NOTIME.*':'TTChi2',
    '/TT_Tune4C_13TeV.*_NOCHI2.*':'TTNoChi2',
    '/TT_Tune4C_13TeV.*TIMEFROMSEED.*':'TTTimeFromSeed',
    '/TT_Tune4C_13TeV.*ONLYTIMECUT.*':'TTTimeCutOnly',
    '/TT_Tune4C_13TeV.*3SIGMACUT.*':'TT3SigmaCut',
    '/TT_Tune4C_13TeV.*3SIGMANEIGHBOUR.*':'TT3SigmaNeighbour',
    '/RelValTTbar.*NOTIME':'RelValTTNoTime',
    '/RelValTTbar.*TIMEFROMSEED.*':'RelValTTTimeFromSeed',
    '/DYJetsToLL_M-50_13TeV.*TIMEFROMSEEDMORE':'DYTimeFromSeed',
    '/DYJetsToLL_M-50_13TeV.*3SIGMACUT':'DY3SigmaCut',
    '/DYJetsToLL_M-50_13TeV.*3SIGMANEIGHBOUR':'DY3SigmaNeighbour',
    '/DYJetsToLL_M-50_13TeV.*ONLYTIMECUT':'DYTimeCutOnly',
    '/DYJetsToLL_M-50_13TeV.*NOTIME':'DYNoTime',
    }

DY3SigmaCut = cfg.MCComponent(
    name = 'DY3SigmaCut',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

DY3SigmaNeighbour = cfg.MCComponent(
    name = 'DY3SigmaNeighbour',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

DYTimeCutOnly = cfg.MCComponent(
    name = 'DYTimeCutOnly',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

DYNoTime = cfg.MCComponent(
    name = 'DYNoTime',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

DYTimeFromSeed = cfg.MCComponent(
    name = 'DYTimeFromSeed',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

RelValTTTimeFromSeed = cfg.MCComponent(
    name = 'RelValTTTimeFromSeed',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )


RelValTTNoTime = cfg.MCComponent(
    name = 'RelValTTNoTime',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TTNoTime = cfg.MCComponent(
    name = 'TTNoTime',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TTChi2 = cfg.MCComponent(
    name = 'TTChi2',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TTNoChi2 = cfg.MCComponent(
    name = 'TTNoChi2',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TTTimeFromSeed = cfg.MCComponent(
    name = 'TTTimeFromSeed',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TTTimeCutOnly = cfg.MCComponent(
    name = 'TTTimeCutOnly',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TT3SigmaCut = cfg.MCComponent(
    name = 'TT3SigmaCut',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

TT3SigmaNeighbour = cfg.MCComponent(
    name = 'TT3SigmaNeighbour',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )



allsamples = [DY3SigmaCut, DY3SigmaNeighbour, DYTimeCutOnly, DYNoTime, DYTimeFromSeed, RelValTTTimeFromSeed, TTNoTime, TTChi2, TTNoChi2, TTTimeFromSeed, TTTimeCutOnly, TT3SigmaCut, TT3SigmaNeighbour]

connect( allsamples, '%PF71X%', 'reco.*root', aliases, cache=True, verbose=True)
connect( [RelValTTNoTime], '%PFTIMERECO%', 'reco.*root', aliases, cache=True, verbose=True)
allsamples = allsamples + [RelValTTNoTime]
for c in allsamples:
    c.splitFactor = 20# splitFactor(c, 5e4)
