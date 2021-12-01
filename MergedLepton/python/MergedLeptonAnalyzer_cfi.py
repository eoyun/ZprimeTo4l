import FWCore.ParameterSet.Config as cms

mergedLeptonAnalyzer = cms.EDAnalyzer("MergedLeptonAnalyzer",
  srcMuon=cms.InputTag("slimmedMuons"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  srcPv=cms.InputTag("offlineSlimmedPrimaryVertices"),
  srcPFMET=cms.InputTag("slimmedMETs"),
  srcPuppiMET=cms.InputTag("slimmedMETsPuppi"),
  ptThres=cms.double(20.),
  drThres=cms.double(0.3)
)
