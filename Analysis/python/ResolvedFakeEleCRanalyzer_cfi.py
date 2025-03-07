import FWCore.ParameterSet.Config as cms

resolvedFakeEleCRanalyzer20UL17 = cms.EDAnalyzer("ResolvedFakeEleCRanalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  srcMET = cms.InputTag("slimmedMETs"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  modifiedTrkIso = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleTrkPtIso"),
  modifiedEcalIso = cms.InputTag("modifiedEcalRecHitIsolationScone2nd","EcalRecHitIso"),
  dPerpIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dPerpIn"),
  rho = cms.InputTag("fixedGridRhoFastjetAll"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL17.json"),
  PUname = cms.string("Collisions17_UltraLegacy_goldenJSON"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  trigList = cms.vstring(
    "HLT_IsoMu27_v*",
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
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
)

resolvedFakeEleCRanalyzer20UL18 = resolvedFakeEleCRanalyzer20UL17.clone(
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL18.json"),
  PUname = cms.string("Collisions18_UltraLegacy_goldenJSON"),
  trigList = cms.vstring(
    "HLT_IsoMu24_v*",
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
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
)

resolvedFakeEleCRanalyzer20UL16 = resolvedFakeEleCRanalyzer20UL17.clone(
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  trigList = cms.vstring(
    "HLT_IsoMu24_v*",
    "HLT_IsoTkMu24_v*"
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
)

resolvedFakeEleCRanalyzer20UL16APV = resolvedFakeEleCRanalyzer20UL17.clone(
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16APV.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  trigList = cms.vstring(
    "HLT_IsoMu24_v*",
    "HLT_IsoTkMu24_v*"
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
)
