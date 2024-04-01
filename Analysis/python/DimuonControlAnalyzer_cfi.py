import FWCore.ParameterSet.Config as cms

dimuonControlAnalyzer = cms.EDAnalyzer("DimuonControlAnalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  genptc = cms.InputTag("prunedGenParticles"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPhoton = cms.InputTag("slimmedPhotons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  ptThres = cms.double(20.),
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root")
)

dimuonControlAnalyzer20UL16APV = dimuonControlAnalyzer.clone(
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016aUL.txt")
)

dimuonControlAnalyzer20UL17 = dimuonControlAnalyzer.clone(
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2017UL.txt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root")
)

dimuonControlAnalyzer20UL18 = dimuonControlAnalyzer.clone(
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2018UL.txt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root")
)
