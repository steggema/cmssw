import itertools
from CMGTools.RootTools.fwlite.Config import printComps
from CMGTools.RootTools.utils.connect import connect
import CMGTools.RootTools.fwlite.Config as cfg

aliases = {
    '/TT_Tune4C_13TeV.*':'TTNoTime',
    }

TTNoTime = cfg.MCComponent(
    name = 'TTNoTime',
    files = [],
    xSection = 400.0,
    nGenEvents = 1,
    triggers = [],
    effCorrFactor = 1 )

allsamples = [TTNoTime]

connect( allsamples, '%PF71X_NOTIME_500', 'reco.*root', aliases, cache=True, verbose=False)
for c in allsamples:
    c.splitFactor = 100# splitFactor(c, 5e4)
