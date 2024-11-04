import FWCore.ParameterSet.Config as cms

mergedEleBkgFilter = cms.EDFilter("MergedEleBkgFilter",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("modifiedEcalRecHitIsolationScone2nd","EcalRecHitIso"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  lheEvent = cms.InputTag("externalLHEProducer"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  ptThres = cms.double(50.),
  select0J = cms.bool(False),
  selectHT = cms.bool(False),
  maxHT = cms.double(70.)
)
