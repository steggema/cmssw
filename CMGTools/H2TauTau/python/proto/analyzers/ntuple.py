#!/bin/env python
from CMGTools.RootTools.utils.DeltaR import deltaR, deltaPhi

def var( tree, varName, type=float ):
    tree.var(varName, type)

def fill( tree, varName, value ):
    tree.fill( varName, value )

# simple particle

def bookParticle( tree, pName ):
    var(tree, '{pName}_pt'.format(pName=pName))
    var(tree, '{pName}_eta'.format(pName=pName))
    var(tree, '{pName}_phi'.format(pName=pName))
    var(tree, '{pName}_charge'.format(pName=pName))

def fillParticle( tree, pName, particle ):
    fill(tree, '{pName}_pt'.format(pName=pName), particle.pt() )
    fill(tree, '{pName}_eta'.format(pName=pName), particle.eta() )
    fill(tree, '{pName}_phi'.format(pName=pName), particle.phi() )
    fill(tree, '{pName}_charge'.format(pName=pName), particle.charge() )

def bookGenParticle(tree, pName):
    bookParticle(tree, pName)
    var(tree, '{pName}_mass'.format(pName=pName))
    var(tree, '{pName}_pdgId'.format(pName=pName))
    
def fillGenParticle( tree, pName, particle ):
    fillParticle( tree, pName, particle )
    fill(tree, '{pName}_mass'.format(pName=pName), particle.mass() )
    fill(tree, '{pName}_pdgId'.format(pName=pName), particle.pdgId() )

# di-tau

def bookDiLepton(tree):
    var( tree, 'visMass')
    var( tree, 'svfitMass')
    var( tree, 'pZetaMET')
    var( tree, 'pZetaVis')
    var( tree, 'pZetaDisc')
    var( tree, 'mt')
    var( tree, 'mtleg1')
    var( tree, 'met')
    var( tree, 'metphi')
    var( tree, 'pthiggs')
    var( tree, 'deltaPhiL1L2')
    var( tree, 'deltaEtaL1L2')
    var( tree, 'deltaRL1L2')
    var( tree, 'deltaPhiL1MET')
    var( tree, 'deltaPhiL2MET')


def fillDiLepton(tree, diLepton):
    fill(tree, 'visMass', diLepton.mass())
    fill(tree, 'svfitMass', diLepton.massSVFit())
    fill(tree, 'pZetaMET', diLepton.pZetaMET())
    fill(tree, 'pZetaVis', diLepton.pZetaVis())
    fill(tree, 'pZetaDisc', diLepton.pZetaDisc())
    fill(tree, 'mt', diLepton.mTLeg2())
    fill(tree, 'mtleg1', diLepton.mTLeg1())
    fill(tree, 'met', diLepton.met().pt())
    fill(tree, 'metphi', diLepton.met().phi())

    pthiggs = (diLepton.leg1().p4()+diLepton.leg2().p4()+diLepton.met().p4()).pt()
    fill(tree, 'pthiggs', pthiggs)
    
    l1eta = diLepton.leg1().eta()
    l2eta = diLepton.leg2().eta()
    l1phi = diLepton.leg1().phi()
    l2phi = diLepton.leg2().phi()
    metphi = diLepton.met().phi()
    fill(tree, 'deltaPhiL1L2', deltaPhi(l1phi, l2phi))
    fill(tree, 'deltaEtaL1L2', abs(l1eta-l2eta))
    fill(tree, 'deltaRL1L2', deltaR(l1eta, l1phi, l2eta, l2phi))
    fill(tree, 'deltaPhiL1MET', deltaPhi(l1phi, metphi))
    fill(tree, 'deltaPhiL2MET', deltaPhi(l2phi, metphi))
    
# lepton

def bookLepton( tree, pName ):
    bookParticle(tree, pName )
    var(tree, '{pName}_relIso05'.format(pName=pName))
    var(tree, '{pName}_dxy'.format(pName=pName))
    var(tree, '{pName}_dz'.format(pName=pName))

    # var(tree, '{pName}_weight'.format(pName=pName))
    # var(tree, '{pName}_triggerWeight'.format(pName=pName))
    # var(tree, '{pName}_triggerEffData'.format(pName=pName))
    # var(tree, '{pName}_triggerEffMC'.format(pName=pName))
    # var(tree, '{pName}_recEffWeight'.format(pName=pName))

def fillLepton( tree, pName, lepton ):
    fillParticle(tree, pName, lepton )
    fill(tree, '{pName}_relIso05'.format(pName=pName), lepton.relIsoAllChargedDB05() )
    fill(tree, '{pName}_dxy'.format(pName=pName), lepton.dxy() )
    fill(tree, '{pName}_dz'.format(pName=pName), lepton.dz() )

    # fill(tree, '{pName}_weight'.format(pName=pName), lepton.weight )
    # fill(tree, '{pName}_triggerWeight'.format(pName=pName), lepton.triggerWeight )
    # fill(tree, '{pName}_triggerEffData'.format(pName=pName), lepton.triggerEffData )
    # fill(tree, '{pName}_triggerEffMC'.format(pName=pName), lepton.triggerEffMC )
    # fill(tree, '{pName}_recEffWeight'.format(pName=pName), lepton.recEffWeight )


# muon


def bookMuon( tree, pName ):
    bookLepton(tree, pName )
    var(tree, '{pName}_mvaIso'.format(pName=pName))
    var(tree, '{pName}_looseId'.format(pName=pName))
    var(tree, '{pName}_tightId'.format(pName=pName))
    var(tree, '{pName}_chargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_chargedAllIso04'.format(pName=pName))
    var(tree, '{pName}_puChargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_neutralHadronIso04'.format(pName=pName))
    var(tree, '{pName}_photonIso04'.format(pName=pName))

    var(tree, '{pName}_NchargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_NchargedAllIso04'.format(pName=pName))
    var(tree, '{pName}_NpuChargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_NneutralHadronIso04'.format(pName=pName))
    var(tree, '{pName}_NphotonIso04'.format(pName=pName))

    var(tree, '{pName}_MaxchargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_MaxchargedAllIso04'.format(pName=pName))
    var(tree, '{pName}_MaxpuChargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_MaxneutralHadronIso04'.format(pName=pName))
    var(tree, '{pName}_MaxphotonIso04'.format(pName=pName))

def fillMuon( tree, pName, muon ):
    fillLepton(tree, pName, muon)
    fill(tree, '{pName}_mvaIso'.format(pName=pName), muon.mvaIso() )
    fill(tree, '{pName}_looseId'.format(pName=pName), muon.looseId() )
    fill(tree, '{pName}_tightId'.format(pName=pName), muon.tightId() )

    fill(tree, '{pName}_chargedHadronIso04'.format(pName=pName), muon.chargedHadronIso(0.4))
    fill(tree, '{pName}_chargedAllIso04'.format(pName=pName), muon.chargedAllIso(0.4))
    fill(tree, '{pName}_puChargedHadronIso04'.format(pName=pName), muon.puChargedHadronIso(0.4))
    fill(tree, '{pName}_neutralHadronIso04'.format(pName=pName), muon.neutralHadronIso(0.4))
    fill(tree, '{pName}_photonIso04'.format(pName=pName), muon.photonIso(0.4))

    fill(tree, '{pName}_NchargedHadronIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(4)).countWithin(0.4, muon.chargedHadronVetos()))
    fill(tree, '{pName}_NchargedAllIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(13)).countWithin(0.4, muon.chargedAllVetos()))
    fill(tree, '{pName}_NpuChargedHadronIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(12)).countWithin(0.4, muon.puChargedHadronVetos()))
    fill(tree, '{pName}_NneutralHadronIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(5)).countWithin(0.4, muon.neutralHadronVetos()))
    fill(tree, '{pName}_NphotonIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(6)).countWithin(0.4, muon.photonVetos() ))

    fill(tree, '{pName}_MaxchargedHadronIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(4)).maxWithin(0.4, muon.chargedHadronVetos()))
    fill(tree, '{pName}_MaxchargedAllIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(13)).maxWithin(0.4, muon.chargedAllVetos()))
    fill(tree, '{pName}_MaxpuChargedHadronIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(12)).maxWithin(0.4, muon.puChargedHadronVetos()))
    fill(tree, '{pName}_MaxneutralHadronIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(5)).maxWithin(0.4, muon.neutralHadronVetos()))
    fill(tree, '{pName}_MaxphotonIso04'.format(pName=pName), (muon.sourcePtr().isoDeposit(6)).maxWithin(0.4, muon.photonVetos() ))

    # import pdb; pdb.set_trace()


# electron


def bookEle( tree, pName ):
    bookLepton(tree, pName )
    var(tree, '{pName}_mvaIso'.format(pName=pName))
    var(tree, '{pName}_mvaTrigV0'.format(pName=pName))
    var(tree, '{pName}_mvaNonTrigV0'.format(pName=pName))
    var(tree, '{pName}_looseId'.format(pName=pName))
    var(tree, '{pName}_tightId'.format(pName=pName))
    var(tree, '{pName}_numberOfMissingHits'.format(pName=pName))
    var(tree, '{pName}_passConversionVeto'.format(pName=pName))
    var(tree, '{pName}_chargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_chargedAllIso04'.format(pName=pName))
    var(tree, '{pName}_puChargedHadronIso04'.format(pName=pName))
    var(tree, '{pName}_neutralHadronIso04'.format(pName=pName))
    var(tree, '{pName}_photonIso04'.format(pName=pName))

def fillEle( tree, pName, ele ):
    fillLepton(tree, pName, ele)
    fill(tree, '{pName}_mvaIso'.format(pName=pName), ele.mvaIso() )
    fill(tree, '{pName}_mvaTrigV0'.format(pName=pName), ele.sourcePtr().electronID("mvaTrigV0") )
    fill(tree, '{pName}_mvaNonTrigV0'.format(pName=pName), ele.sourcePtr().electronID("mvaNonTrigV0") )
    fill(tree, '{pName}_looseId'.format(pName=pName), ele.looseIdForEleTau() )
    fill(tree, '{pName}_tightId'.format(pName=pName), ele.tightIdForEleTau() )
    fill(tree, '{pName}_numberOfMissingHits'.format(pName=pName), ele.numberOfHits() )
    fill(tree, '{pName}_passConversionVeto'.format(pName=pName), ele.passConversionVeto() )
    fill(tree, '{pName}_chargedHadronIso04'.format(pName=pName), ele.chargedHadronIso(0.4))
    fill(tree, '{pName}_chargedAllIso04'.format(pName=pName), ele.chargedAllIso()) # updated veto cones in HTauTauEle
    fill(tree, '{pName}_puChargedHadronIso04'.format(pName=pName), ele.puChargedHadronIso(0.4))
    fill(tree, '{pName}_neutralHadronIso04'.format(pName=pName), ele.neutralHadronIso(0.4))
    fill(tree, '{pName}_photonIso04'.format(pName=pName), ele.photonIso())# updated veto cones in HTauTauEle

# tau 

def bookTau( tree, pName ):
    bookLepton(tree, pName )
    var(tree, '{pName}_decayModeFinding'.format(pName=pName))
    var(tree, '{pName}_veryLooseIso'.format(pName=pName))
    var(tree, '{pName}_looseIso'.format(pName=pName))
    var(tree, '{pName}_mediumIso'.format(pName=pName))
    var(tree, '{pName}_tightIso'.format(pName=pName))

    var(tree, '{pName}_againstMuonTight'.format(pName=pName))
    var(tree, '{pName}_againstMuonTight2'.format(pName=pName))
    var(tree, '{pName}_againstElectronLoose'.format(pName=pName))    

    var(tree, '{pName}_byLooseIsoMVA'.format(pName=pName))    
    var(tree, '{pName}_againstElectronMVA'.format(pName=pName))    
    var(tree, '{pName}_againstElectronTightMVA2'.format(pName=pName))
    var(tree, '{pName}_againstElectronTightMVA3'.format(pName=pName))
    var(tree, '{pName}_againstElectronMedium'.format(pName=pName))
    var(tree, '{pName}_againstElectronMVA3Medium'.format(pName=pName))
    var(tree, '{pName}_againstElectronMVA3raw'.format(pName=pName))
    
    var(tree, '{pName}_againstMuonLoose'.format(pName=pName))
    var(tree, '{pName}_againstMuonLoose2'.format(pName=pName))

    var(tree, '{pName}_rawMvaIso'.format(pName=pName))
    var(tree, '{pName}_looseMvaIso'.format(pName=pName))
    var(tree, '{pName}_mediumMvaIso'.format(pName=pName))
    var(tree, '{pName}_tightMvaIso'.format(pName=pName))

    var(tree, '{pName}_threeHitIso'.format(pName=pName))
   
    var(tree, '{pName}_EOverp'.format(pName=pName))
    var(tree, '{pName}_decayMode'.format(pName=pName))
    var(tree, '{pName}_mass'.format(pName=pName))
    var(tree, '{pName}_zImpact'.format(pName=pName))

    var(tree, '{pName}_genVisPt'.format(pName=pName))
    var(tree, '{pName}_genVisEta'.format(pName=pName))
    var(tree, '{pName}_genVisPhi'.format(pName=pName))
    var(tree, '{pName}_genVisMass'.format(pName=pName))

def fillTau( tree, pName, tau ):
    fillLepton(tree, pName, tau)
    fill(tree, '{pName}_decayModeFinding'.format(pName=pName), tau.tauID("decayModeFinding"))
    fill(tree, '{pName}_veryLooseIso'.format(pName=pName),
         tau.tauID("byVLooseCombinedIsolationDeltaBetaCorr"))
    fill(tree, '{pName}_looseIso'.format(pName=pName),
         tau.tauID("byLooseCombinedIsolationDeltaBetaCorr"))
    fill(tree, '{pName}_mediumIso'.format(pName=pName),
         tau.tauID("byMediumCombinedIsolationDeltaBetaCorr"))
    fill(tree, '{pName}_tightIso'.format(pName=pName),
         tau.tauID("byTightCombinedIsolationDeltaBetaCorr"))

    fill(tree, '{pName}_againstMuonTight'.format(pName=pName),
         tau.tauID("againstMuonTight"))
    fill(tree, '{pName}_againstMuonTight2'.format(pName=pName),
         tau.tauID("againstMuonTight2"))
    fill(tree, '{pName}_againstElectronLoose'.format(pName=pName),
         tau.tauID("againstElectronLoose"))

    fill(tree, '{pName}_byLooseIsoMVA'.format(pName=pName),
         tau.tauID("byLooseIsoMVA"))
    fill(tree, '{pName}_againstElectronMVA'.format(pName=pName),
         tau.tauID("againstElectronMVA"))
    fill(tree, '{pName}_againstElectronTightMVA2'.format(pName=pName),
         tau.tauID("againstElectronTightMVA2"))
    fill(tree, '{pName}_againstElectronTightMVA3'.format(pName=pName),
         tau.tauID("againstElectronTightMVA3"))
    fill(tree, '{pName}_againstElectronMedium'.format(pName=pName),
         tau.tauID("againstElectronMedium"))
    fill(tree, '{pName}_againstElectronMVA3Medium'.format(pName=pName),
         tau.electronMVA3Medium())
    fill(tree, '{pName}_againstElectronMVA3raw'.format(pName=pName),
         tau.tauID('againstElectronMVA3raw') )
    

    fill(tree, '{pName}_againstMuonLoose'.format(pName=pName),
         tau.tauID("againstMuonLoose"))
    fill(tree, '{pName}_againstMuonLoose2'.format(pName=pName),
         tau.tauID("againstMuonLoose2"))

    fill(tree, '{pName}_rawMvaIso'.format(pName=pName),
         tau.tauID("byRawIsoMVA"))
    fill(tree, '{pName}_looseMvaIso'.format(pName=pName),
         tau.tauID("byLooseIsoMVA"))
    fill(tree, '{pName}_mediumMvaIso'.format(pName=pName),
         tau.tauID("byMediumIsoMVA"))
    fill(tree, '{pName}_tightMvaIso'.format(pName=pName),
         tau.tauID("byTightIsoMVA"))

    fill(tree, '{pName}_threeHitIso'.format(pName=pName),
         tau.tauID("byCombinedIsolationDeltaBetaCorrRaw3Hits"))

    fill(tree, '{pName}_rawMvaIso'.format(pName=pName),
         tau.tauID("byRawIsoMVA"))
    fill(tree, '{pName}_EOverp'.format(pName=pName),
         tau.calcEOverP())
    fill(tree, '{pName}_decayMode'.format(pName=pName),
         tau.decayMode())
    fill(tree, '{pName}_mass'.format(pName=pName),
         tau.mass())
    fill(tree, '{pName}_zImpact'.format(pName=pName),
         tau.zImpact())

    if hasattr(tau, 'genVisP4'):
        fill(tree, '{pName}_genVisPt'.format(pName=pName), tau.genVisP4.pt())
        fill(tree, '{pName}_genVisEta'.format(pName=pName), tau.genVisP4.eta())
        fill(tree, '{pName}_genVisPhi'.format(pName=pName), tau.genVisP4.phi())
        fill(tree, '{pName}_genVisMass'.format(pName=pName), tau.genVisP4.mass())

# jet

def bookJet( tree, pName ):
    bookParticle(tree, pName )
    var(tree, '{pName}_puMvaFull53X'.format(pName=pName))
    var(tree, '{pName}_puMvaFull'.format(pName=pName))
    var(tree, '{pName}_puMvaSimple'.format(pName=pName))
    var(tree, '{pName}_puMvaCutBased'.format(pName=pName))
    var(tree, '{pName}_looseJetId'.format(pName=pName))
    # var(tree, '{pName}_btagMVA'.format(pName=pName))
    var(tree, '{pName}_area'.format(pName=pName))
    var(tree, '{pName}_genJetPt'.format(pName=pName))
    var(tree, '{pName}_puJetIdPassed'.format(pName=pName))
    var(tree, '{pName}_pfJetIdPassed'.format(pName=pName))
    var(tree, '{pName}_partonFlavour'.format(pName=pName))
    var(tree, '{pName}_CHfraction'.format(pName=pName))
    var(tree, '{pName}_PHfraction'.format(pName=pName))
    var(tree, '{pName}_NHfraction'.format(pName=pName))

def fillJet( tree, pName, jet ):
    fillParticle(tree, pName, jet )
    fill(tree, '{pName}_puMvaFull53X'.format(pName=pName), jet.puMva('full53x') )
    fill(tree, '{pName}_puMvaFull'.format(pName=pName), jet.puMva('full') )
    fill(tree, '{pName}_puMvaSimple'.format(pName=pName), jet.puMva('simple'))
    fill(tree, '{pName}_puMvaCutBased'.format(pName=pName), jet.puMva('cut-based'))
    fill(tree, '{pName}_looseJetId'.format(pName=pName), jet.looseJetId())
    # fill(tree, '{pName}_btagMVA'.format(pName=pName), jet.btagMVA)
    fill(tree, '{pName}_area'.format(pName=pName), jet.jetArea())
    if hasattr(jet, 'genJet') and jet.genJet:
        fill(tree, '{pName}_genJetPt'.format(pName=pName), jet.genJet.pt())
    fill(tree, '{pName}_puJetIdPassed'.format(pName=pName), jet.puJetIdPassed)
    fill(tree, '{pName}_pfJetIdPassed'.format(pName=pName), jet.pfJetIdPassed)
    fill(tree, '{pName}_partonFlavour'.format(pName=pName), jet.partonFlavour())
    fill(tree, '{pName}_CHfraction'.format(pName=pName), jet.component(1).fraction())
    fill(tree, '{pName}_PHfraction'.format(pName=pName), jet.component(4).fraction())
    fill(tree, '{pName}_NHfraction'.format(pName=pName), jet.component(5).fraction())

# vbf

def bookVBF( tree, pName ):
    var(tree, '{pName}_mjj'.format(pName=pName))
    var(tree, '{pName}_deta'.format(pName=pName))
    var(tree, '{pName}_nCentral'.format(pName=pName))
    var(tree, '{pName}_mva'.format(pName=pName))
    var(tree, '{pName}_jdphi'.format(pName=pName))
    var(tree, '{pName}_dijetpt'.format(pName=pName))
    var(tree, '{pName}_dijetphi'.format(pName=pName))
    var(tree, '{pName}_hdijetphi'.format(pName=pName))
    var(tree, '{pName}_visjeteta'.format(pName=pName))
    var(tree, '{pName}_ptvis'.format(pName=pName))
    
def fillVBF( tree, pName, vbf ):
    fill(tree, '{pName}_mjj'.format(pName=pName), vbf.mjj )
    fill(tree, '{pName}_deta'.format(pName=pName), vbf.deta )
    fill(tree, '{pName}_nCentral'.format(pName=pName), len(vbf.centralJets) )
    fill(tree, '{pName}_mva'.format(pName=pName), vbf.mva )
    fill(tree, '{pName}_jdphi'.format(pName=pName), vbf.dphi )
    fill(tree, '{pName}_dijetpt'.format(pName=pName), vbf.dijetpt )
    fill(tree, '{pName}_dijetphi'.format(pName=pName), vbf.dijetphi )
    fill(tree, '{pName}_hdijetphi'.format(pName=pName), vbf.dphidijethiggs )
    fill(tree, '{pName}_visjeteta'.format(pName=pName), vbf.visjeteta )
    fill(tree, '{pName}_ptvis'.format(pName=pName), vbf.ptvis )
    
