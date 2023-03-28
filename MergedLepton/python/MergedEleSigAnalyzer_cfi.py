import FWCore.ParameterSet.Config as cms

mergedEleSigAnalyzer = cms.EDAnalyzer("MergedEleSigAnalyzer",
  srcEle=cms.InputTag("slimmedElectrons"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  conversions = cms.InputTag("reducedEgamma:reducedConversions"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  generator = cms.InputTag("generator"),
  ptThres = cms.double(20.),
  ptThres2nd = cms.double(10.),
  drThres = cms.double(0.3)
)
