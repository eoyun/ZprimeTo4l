import FWCore.ParameterSet.Config as cms

resolvedFakeMuCRanalyzer20UL17 = cms.EDAnalyzer("ResolvedFakeMuCRanalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  srcMET = cms.InputTag("slimmedMETs"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL17.json"),
  PUname = cms.string("Collisions17_UltraLegacy_goldenJSON"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  trigList = cms.vstring(
    "HLT_Ele32_WPTight_Gsf_L1DoubleEG_v*",
  ),
  emulateEle32WPTightGsf = cms.bool(True),
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
)

resolvedFakeMuCRanalyzer20UL18 = resolvedFakeMuCRanalyzer20UL17.clone(
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL18.json"),
  PUname = cms.string("Collisions18_UltraLegacy_goldenJSON"),
  trigList = cms.vstring(
    "HLT_Ele32_WPTight_Gsf_v*",
  ),
  emulateEle32WPTightGsf = cms.bool(False),
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
)

resolvedFakeMuCRanalyzer20UL16 = resolvedFakeMuCRanalyzer20UL17.clone(
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  trigList = cms.vstring(
    "HLT_Ele27_WPTight_Gsf_v*"
  ),
  emulateEle32WPTightGsf = cms.bool(False),
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
)

resolvedFakeMuCRanalyzer20UL16APV = resolvedFakeMuCRanalyzer20UL17.clone(
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16APV.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  trigList = cms.vstring(
    "HLT_Ele27_WPTight_Gsf_v*"
  ),
  emulateEle32WPTightGsf = cms.bool(False),
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
)
