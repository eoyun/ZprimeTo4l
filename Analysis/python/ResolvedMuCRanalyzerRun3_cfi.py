import FWCore.ParameterSet.Config as cms

# ResolvedMuCRanalyzerRun3 — Run-3 (2022 preEE) 4mu resolved CR analyzer.
# 물리 로직은 Run-2 ResolvedMuCRanalyzer 와 동일. prefiring 제거, egamma VID-only.
# payload(SF/PU/RoccoR/FF)는 아직 Run-3 버전 부재 → Run-2 UL placeholder + TODO.
resolvedMuCRanalyzerRun3 = cms.EDAnalyzer("ResolvedMuCRanalyzerRun3",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcMET = cms.InputTag("slimmedMETs"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  genptc = cms.InputTag("prunedGenParticles"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  # Run-3 MET filter (HBHENoise/HBHENoiseIso 제거, Flag_ 접두 통일, ecalBadCalib 추가)
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_eeBadScFilter",
    "Flag_ecalBadCalibFilter"
  ),
  # Run-3 single-mu / high-pT-mu triggers (Run-2 TkMu50 부재 -> Cascade/HighPtTk)
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_CascadeMu100_v*",
    "HLT_HighPtTkMu100_v*"
  ),
  # ---- payloads: Run-3 버전 부재 → Run-2 UL placeholder (TODO 교체) ----
  # TODO(Run-3 2022preEE): RoccoR Run-3 공식 파일 부재. 유지, muon momentum 마이그레이션은 후속.
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  # TODO(Run-3 2022preEE): Muon POG Run-3 highPt SF json 으로 교체.
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_schemaV2.json"),
  muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL16.root"),
  # TODO(Run-3 2022preEE): LUM POG Run-3 PU json 으로 교체.
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16.json"),
  PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  # TODO(Run-3 2022preEE): FF 재측정 payload 로 교체.
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_run2.root"),
  ptThresTrig = cms.double(52.),
  ptThres = cms.double(20.),
  ffSystCL95 = cms.double(0.0),
  # 고pT muon scale/smear (GE-method). TODO(Run-3): Run-3 재검증. 값은 Run-2 유지.
  muScaleBias = cms.vdouble(0.002, 0.049, -0.045, -0.061, 0.016, -0.108,
                            -0.005, -0.064, 0.023, -0.004, 0.041, -0.023,
                            0.026, -0.142, 0.039, 0.000, 0.003, -0.091 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.0102, 6.77e-05, -3.72e-08, 8.53e-12, 0.0129, 6.48e-05, -3.04e-08, 6.63e-12)
)
