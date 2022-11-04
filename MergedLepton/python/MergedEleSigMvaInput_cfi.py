import FWCore.ParameterSet.Config as cms

mergedEleSigMvaInput = cms.EDAnalyzer("MergedEleSigMvaInput",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
  nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  rho=cms.InputTag("fixedGridRhoFastjetAll"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  generator = cms.InputTag("generator"),
  ptThres=cms.double(20.),
  drThres=cms.double(0.3)
)
