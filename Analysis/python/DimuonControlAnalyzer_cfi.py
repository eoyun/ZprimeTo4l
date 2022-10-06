import FWCore.ParameterSet.Config as cms

dimuonControlAnalyzer = cms.EDAnalyzer("DimuonControlAnalyzer",
  isMC = cms.untracked.bool(True),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  ptThres = cms.double(25.)
)
