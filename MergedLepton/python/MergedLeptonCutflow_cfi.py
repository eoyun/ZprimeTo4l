import FWCore.ParameterSet.Config as cms

mergedLeptonCutflow = cms.EDAnalyzer("MergedLeptonCutflow",
  isMC = cms.untracked.bool(True),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  cutflow_modifiedHEEP = cms.InputTag("mergedLeptonIDProducer","cutflowModifiedHEEP"),
  cutflow_HEEP = cms.InputTag("mergedLeptonIDProducer","cutflowHEEP"),
  cutflow_mergedElectron = cms.InputTag("mergedLeptonIDProducer","cutflowMergedElectron"),
  mva_mergedElectron = cms.InputTag("mergedLeptonIDProducer","mvaMergedElectron"),
  openingAngle_mergedElectron = cms.InputTag("mergedLeptonIDProducer","openingAngleMergedElectron"),
  nresolvedElectron = cms.InputTag("mergedLeptonIDProducer","nresolvedElectron"),
  nmergedElectron = cms.InputTag("mergedLeptonIDProducer","nmergedElectron")
)