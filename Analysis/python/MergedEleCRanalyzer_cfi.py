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
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/FF_20UL16.root"),
  trigList = cms.vstring(
  # https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.),
  ssBoundary = cms.double(200.),
  osBoundary = cms.double(500.)
)

mergedEleCRanalyzer20UL16APV = mergedEleCRanalyzer.clone(
  cutflow_modifiedHEEP = cms.InputTag("mergedLeptonIDProducer20UL16APV","cutflowModifiedHEEP"),
  cutflow_HEEP = cms.InputTag("mergedLeptonIDProducer20UL16APV","cutflowHEEP"),
  status_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL16APV","statusMergedElectron"),
  mva_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL16APV","mvaMergedElectron"),
  GSFtype_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL16APV","GSFtypeMergedElectron"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16APV.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL17 = mergedEleCRanalyzer.clone(
  cutflow_modifiedHEEP = cms.InputTag("mergedLeptonIDProducer20UL17","cutflowModifiedHEEP"),
  cutflow_HEEP = cms.InputTag("mergedLeptonIDProducer20UL17","cutflowHEEP"),
  status_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL17","statusMergedElectron"),
  mva_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL17","mvaMergedElectron"),
  GSFtype_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL17","GSFtypeMergedElectron"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL17.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL18 = mergedEleCRanalyzer.clone(
  cutflow_modifiedHEEP = cms.InputTag("mergedLeptonIDProducer20UL18","cutflowModifiedHEEP"),
  cutflow_HEEP = cms.InputTag("mergedLeptonIDProducer20UL18","cutflowHEEP"),
  status_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL18","statusMergedElectron"),
  mva_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL18","mvaMergedElectron"),
  GSFtype_mergedElectron = cms.InputTag("mergedLeptonIDProducer20UL18","GSFtypeMergedElectron"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL18.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle25_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(28.)
)
