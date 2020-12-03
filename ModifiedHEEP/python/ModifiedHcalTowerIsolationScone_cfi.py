import FWCore.ParameterSet.Config as cms
from ZprimeTo4l.ModifiedHEEP.ModifiedElectronTrackIsolations_cfi import trkIsol03CfgV2

ModifiedHcalTowerIsolationScone = cms.EDProducer("ModifiedTowerIsolationProducer",
    absolut = cms.bool(True),
    intRadius = cms.double(0.15), # to be orthogonal with the H/E ID cut
    extRadius = cms.double(0.3),
    towerProducer = cms.InputTag("towerMaker"),
    etMin = cms.double(0.0),
    Depth = cms.int32(-1),
    emObjectProducer = cms.InputTag("gedGsfElectrons")
)


ModifiedHcalDepth1TowerIsolationScone = cms.EDProducer("ModifiedTowerIsolationProducer",
    absolut = cms.bool(True),
    intRadius = cms.double(0.15), # to be orthogonal with the H/E ID cut
    extRadius = cms.double(0.3),
    towerProducer = cms.InputTag("towerMaker"),
    etMin = cms.double(0.0),
    Depth = cms.int32(1),
    emObjectProducer = cms.InputTag("gedGsfElectrons"),
    gsfTrks = cms.InputTag("electronGsfTracks"),
    trkIsoConfig= trkIsol03CfgV2
)

ModifiedHcalDepth2TowerIsolationScone = cms.EDProducer("ModifiedTowerIsolationProducer",
    absolut = cms.bool(True),
    intRadius = cms.double(0.15), # to be orthogonal with the H/E ID cut
    extRadius = cms.double(0.3),
    towerProducer = cms.InputTag("towerMaker"),
    etMin = cms.double(0.0),
    Depth = cms.int32(2),
    emObjectProducer = cms.InputTag("gedGsfElectrons")
)
