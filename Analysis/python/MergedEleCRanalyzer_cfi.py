import FWCore.ParameterSet.Config as cms

mergedEleCRanalyzer = cms.EDAnalyzer("MergedEleCRanalyzer",
  isMC = cms.untracked.bool(True),
  isDY = cms.untracked.bool(False),
  selectHT = cms.untracked.bool(False),
  selectZpT = cms.untracked.bool(False),
  maxHT = cms.double(70.),
  maxVpT = cms.double(50.),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  lheEvent = cms.InputTag("externalLHEProducer"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  recoSFlowPtPath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptBelow20_EGM2D_UL2016postVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16.root"),
  trigList = cms.vstring(
  # https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.),
  ssBoundary = cms.double(200.)
)

mergedEleCRanalyzer20UL16APV = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  recoSFlowPtPath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptBelow20_EGM2D_UL2016preVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16APV.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL17 = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  recoSFlowPtPath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptBelow20_EGM2D_UL2017.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL17.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL18 = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  recoSFlowPtPath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptBelow20_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL18.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle25_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(28.)
)
