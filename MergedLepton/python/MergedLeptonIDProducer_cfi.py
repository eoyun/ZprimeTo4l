import FWCore.ParameterSet.Config as cms

mergedLeptonIDProducer = cms.EDProducer("MergedLeptonIDProducer", # 20UL16
  srcEle=cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  etThresEB = cms.double(200.),
  etThresEE = cms.double(150.),
  minEt = cms.double(50.),
  xgbPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL16.xml"),
  meanstdPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL16.csv"),
  xgbPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL16.xml"),
  meanstdPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL16.csv"),
  xgbPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL16.xml"),
  meanstdPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL16.csv"),
  xgbPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL16.xml"),
  meanstdPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL16.csv"),
  xgbPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL16.xml"),
  meanstdPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL16.csv"),
  xgbPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL16.xml"),
  meanstdPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL16.csv")
)

mergedLeptonIDProducer20UL16APV = mergedLeptonIDProducer.clone(
  xgbPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL16APV.xml"),
  meanstdPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL16APV.csv"),
  xgbPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL16APV.xml"),
  meanstdPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL16APV.csv"),
  xgbPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL16APV.xml"),
  meanstdPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL16APV.csv"),
  xgbPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL16APV.xml"),
  meanstdPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL16APV.csv"),
  xgbPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL16APV.xml"),
  meanstdPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL16APV.csv"),
  xgbPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL16APV.xml"),
  meanstdPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL16APV.csv")
)

mergedLeptonIDProducer20UL17 = mergedLeptonIDProducer.clone(
  xgbPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL17.xml"),
  meanstdPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL17.csv"),
  xgbPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL17.xml"),
  meanstdPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL17.csv"),
  xgbPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL17.xml"),
  meanstdPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL17.csv"),
  xgbPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL17.xml"),
  meanstdPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL17.csv"),
  xgbPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL17.xml"),
  meanstdPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL17.csv"),
  xgbPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL17.xml"),
  meanstdPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL17.csv")
)

mergedLeptonIDProducer20UL18 = mergedLeptonIDProducer.clone(
  xgbPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL18.xml"),
  meanstdPathDR1Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR1Et2_EB_20UL18.csv"),
  xgbPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL18.xml"),
  meanstdPathDR2Et1EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EB_20UL18.csv"),
  xgbPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL18.xml"),
  meanstdPathDR2Et2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EB_20UL18.csv"),
  xgbPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL18.xml"),
  meanstdPathDR2Et1EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et1_EE_20UL18.csv"),
  xgbPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL18.xml"),
  meanstdPathDR2Et2EE=cms.FileInPath("ZprimeTo4l/MergedLepton/data/DR2Et2_EE_20UL18.csv"),
  xgbPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL18.xml"),
  meanstdPathBkgEt2EB=cms.FileInPath("ZprimeTo4l/MergedLepton/data/NoneEt2_EB_20UL18.csv")
)
