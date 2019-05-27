import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

pfCandidatesTable = cms.EDProducer("SimpleCandidateFlatTableProducer",
    src = cms.InputTag("packedPFCandidates"),
    cut = cms.string(""), # filtered already above
    name = cms.string("PF"),
    doc  = cms.string("PF candidates"),
    singleton = cms.bool(False), # the number of entries is variable
    extension = cms.bool(False), # this is the main table for pf candidates
    variables = cms.PSet(P4Vars,
        dz = Var("dz",float,doc="dz (with sign) wrt first PV, in cm",precision=10),
        dxy = Var("dxy",float,doc="dxy (with sign) wrt first PV, in cm",precision=10),
        charge = Var("charge", int, doc="charge"),
        fromPV = Var("fromPV", int, doc="isolated track comes from PV"),
        puppiWeight = Var("puppiWeight", float, doc="PUPPI weight"),
        puppiWeightNoLep = Var("puppiWeightNoLep", float, doc="PUPPI weight (no lep)"),
        pdgId = Var("pdgId",int,doc="PDG id of PF cand"),
    ),
)
