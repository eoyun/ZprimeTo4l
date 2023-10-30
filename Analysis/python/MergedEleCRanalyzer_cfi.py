import FWCore.ParameterSet.Config as cms

mergedEleCRanalyzer = cms.EDAnalyzer("MergedEleCRanalyzer",
  isMC = cms.bool(True),
  selectHT = cms.bool(False),
  selectZpT = cms.bool(False),
  maxHT = cms.double(70.),
  maxVpT = cms.double(50.),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  lheEvent = cms.InputTag("externalLHEProducer"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  addPackedCandMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddPackedCand"),
  union5x5dEtaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dEtaIn"),
  union5x5Energy = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5Energy"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  trigList = cms.vstring(
  # https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter"
  ),
  etThres1 = cms.double(35.),
  ssBoundary = cms.double(200.)
)

mergedEleCRanalyzer20UL16APV = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16APV.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16APV.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL17 = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL17.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*"
  ),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL18 = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL18.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle25_CaloIdL_MW_v*"
  ),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  etThres1 = cms.double(28.)
)
