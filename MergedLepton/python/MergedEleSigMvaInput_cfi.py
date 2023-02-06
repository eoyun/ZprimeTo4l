import FWCore.ParameterSet.Config as cms

mergedEleSigMvaInput = cms.EDAnalyzer("MergedEleSigMvaInput",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("modifiedEcalRecHitIsolationScone2nd","EcalRecHitIso"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  generator = cms.InputTag("generator"),
  EBrecHits = cms.InputTag("reducedEgamma","reducedEBRecHits"),
  EErecHits = cms.InputTag("reducedEgamma","reducedEERecHits"),
  ptThres=cms.double(20.),
  drThres=cms.double(0.3)
)
