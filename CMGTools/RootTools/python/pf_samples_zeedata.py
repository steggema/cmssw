from CMGTools.RootTools.utils.connect import connect
import CMGTools.RootTools.fwlite.Config as cfg

aliases = {
    '/DoubleElectron.*PFTIME_DATA_NOTIME':'Zee_NoTime',
    '/DoubleElectron.*PFTIME_DATA_3DDEFAULT':'Zee_3DDefault',
    '/DoubleElectron.*PFTIME_DATA_TIMECUTS':'Zee_TimeCuts',
    '/DoubleElectron.*PFTIME_DATA_3DLooser':'Zee_3DLooser',
    }

Zee_NoTime = cfg.DataComponent(
    name = 'Zee_NoTime',
    files = [],
    intLumi = 1000., # dummy
    triggers = [])

Zee_3DDefault = cfg.DataComponent(
    name = 'Zee_3DDefault',
    files = [],
    intLumi = 1000., # dummy
    triggers = [])
Zee_3DLooser = cfg.DataComponent(
    name = 'Zee_3DLooser',
    files = [],
    intLumi = 1000., # dummy
    triggers = [])
Zee_TimeCuts = cfg.DataComponent(
    name = 'Zee_TimeCuts',
    files = [],
    intLumi = 1000., # dummy
    triggers = [])

allsamples = [Zee_NoTime, Zee_3DDefault, Zee_TimeCuts, Zee_3DLooser]

connect( allsamples, '%PFTIME%', 'reco.*root', aliases, cache=True, verbose=True)
for c in allsamples:
    c.splitFactor = 50# splitFactor(c, 5e4)
