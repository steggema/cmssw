import numpy as np

from DataFormats.FWLite import Events, Handle

events = Events('DeepMETTest.root')

products = {
    'gen':{'label':'genMetTrue', 'handle':Handle('vector<reco::GenMET>')},
    'pf':{'label':'slimmedMETs', 'handle':Handle('vector<pat::MET>')},
    'puppi':{'label':'slimmedMETsPuppi', 'handle':Handle('vector<pat::MET>')},
    'deep':{'label':'deepMETProducer', 'handle':Handle('vector<pat::MET>')},
}

rec_labels = ['pf', 'puppi', 'deep']

px = {label:[] for label in rec_labels}
py = {label:[] for label in rec_labels}

for ev in events:
    for prod in products.values():
        ev.getByLabel(prod['label'], prod['handle'])

    gen = products['gen']['handle'].product()[0]
    mets = {label:products[label]['handle'].product()[0] for label in rec_labels}

    for label in rec_labels:
        met = mets[label]
        px[label].append(met.px() - gen.px())
        py[label].append(met.py() - gen.py())

for label in rec_labels:
    print 'Variance px {:.2f} py {:.2f}'.format(np.var(px[label]), np.var(py[label]))
for label in rec_labels:
    print 'Mean px {:.2f} py {:.2f}'.format(np.mean(px[label]), np.mean(py[label]))
