import FWCore.ParameterSet.Config as cms

mergedEleCRanalyzer = cms.EDAnalyzer("MergedEleCRanalyzer",
  isMC = cms.untracked.bool(True),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  cutflow_modifiedHEEP = cms.InputTag("mergedLeptonIDProducer","cutflowModifiedHEEP"),
  cutflow_HEEP = cms.InputTag("mergedLeptonIDProducer","cutflowHEEP"),
  status_mergedElectron = cms.InputTag("mergedLeptonIDProducer","statusMergedElectron"),
  mva_mergedElectron = cms.InputTag("mergedLeptonIDProducer","mvaMergedElectron"),
  GSFtype_mergedElectron = cms.InputTag("mergedLeptonIDProducer","GSFtypeMergedElectron")
)
