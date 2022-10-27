import FWCore.ParameterSet.Config as cms

mergedMuonAnalyzer = cms.EDAnalyzer("MergedMuonAnalyzer",
  isMC = cms.untracked.bool(True),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcMET = cms.InputTag("slimmedMETsPuppi"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  ptThres = cms.double(200.),
  drThres = cms.double(0.3),
  ratioBarrel = cms.double(0.8),
  ratioEndcap = cms.double(0.5)
)
