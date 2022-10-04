import FWCore.ParameterSet.Config as cms

mergedFakeAnalyzer = cms.EDAnalyzer("MergedFakeAnalyzer",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  trkIsoMap=cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
  ecalIsoMap=cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
  nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  rho=cms.InputTag("fixedGridRhoFastjetAll"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  generator = cms.InputTag("generator"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  KFParameters = cms.PSet(
        maxDistance = cms.double(0.01),
        maxNbrOfIterations = cms.int32(10)
  ),
  ptThres=cms.double(20.),
  drThres=cms.double(0.3)
)
