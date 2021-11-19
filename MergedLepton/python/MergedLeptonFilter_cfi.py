import FWCore.ParameterSet.Config as cms

mergedLeptonFilter = cms.EDFilter("MergedLeptonFilter",
  srcMuon=cms.InputTag("muons"),
  srcGenPtc=cms.InputTag("genParticles"),
  srcPv=cms.InputTag("offlinePrimaryVertices"),
  ptThres=cms.double(20.),
  drThres=cms.double(0.3)
)
