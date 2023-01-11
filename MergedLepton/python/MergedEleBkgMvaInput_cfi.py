import FWCore.ParameterSet.Config as cms

mergedEleBkgMvaInput = cms.EDAnalyzer("MergedEleBkgMvaInput",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("modifiedEcalRecHitIsolationScone2nd","EcalRecHitIso"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  generator = cms.InputTag("generator"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  EBrecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
  EErecHits = cms.InputTag("reducedEgamma","reducedEERecHits"),
  ptThres=cms.double(50.),
  drThres=cms.double(0.3)
)
