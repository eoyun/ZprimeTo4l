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
  GSFtype_mergedElectron = cms.InputTag("mergedLeptonIDProducer","GSFtypeMergedElectron"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  trigList = cms.vstring(
  # https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.)
)
