import FWCore.ParameterSet.Config as cms

mergedMuCRanalyzer = cms.EDAnalyzer("MergedMuCRanalyzer",
  isMC = cms.bool(True),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  genptc = cms.InputTag("prunedGenParticles"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  srcMET = cms.InputTag("slimmedMETs"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
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
  trigList = cms.vstring(
  # https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonHLT
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  ),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL16.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL16.root"),
  ptThres = cms.double(50.),
  ptMuThres = cms.double(20.),
  drThres = cms.double(0.3),
  drThresCR = cms.double(0.6),
  ratioThresLo = cms.double(0.5),
  ratioThresHi = cms.double(1.5),
  year = cms.string("2016nonAPV")
)

mergedMuCRanalyzer20UL16APV = mergedMuCRanalyzer.clone(
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
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  ),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL16APV.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016aUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_preVFP_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_preVFP_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_preVFP_schemaV2.json"),
  year = cms.string("2016APV")
)

mergedMuCRanalyzer20UL17 = mergedMuCRanalyzer.clone(
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_OldMu100_v*",
    "HLT_TkMu100_v*"
  ),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL17.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2017UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2017_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2017_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2017_schemaV2.json"),
  year = cms.string("2017")
)

mergedMuCRanalyzer20UL18 = mergedMuCRanalyzer.clone(
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_OldMu100_v*",
    "HLT_TkMu100_v*"
  ),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL18.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2018UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2018_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2018_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2018_schemaV2.json"),
  year = cms.string("2018")
)
