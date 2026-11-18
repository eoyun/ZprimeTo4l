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
  # MC: BX=0 true/observed PU is saved unweighted before the analyzer selections.
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
  # ---- payloads ----
  # 2016/Run-2 payload 전부 비활성화 (analyzer 에서 SF/PU/FF/momentum 보정 주석 처리됨).
  # cms.FileInPath 는 C++ 사용 여부와 무관하게 config 로드 시 파일 존재를 검증하므로,
  # 없는 파일(muIso20UL16.root, RMFF_run2.root)까지 포함해 아래를 통째로 주석 처리한다.
  # Run-3 payload 준비되면: 아래 되살리고 analyzer 의 대응 주석도 해제.
  #
  # rochesterPath 만 유지 — checkIso(순수 geometry) 용 MuonCorrectionHelper 1-인자 생성에 필요.
  # (momentum 보정 자체는 analyzer 에서 미적용. RoccoR 파일은 load 되지만 물리에 반영 안 됨.)
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  # TODO(Run-3 2022preEE): Muon POG Run-3 highPt SF json 으로 교체 후 되살리기.
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2022_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2022_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2022_schemaV2.json"),
  # muonBoostIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/muIso20UL16.root"),   # 파일 없음
  # TODO(Run-3 2022preEE): LUM POG Run-3 PU json 으로 교체 후 되살리기.
  # PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16.json"),
  # PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  # TODO(Run-3 2022preEE): FF 재측정 payload 로 교체 후 되살리기.
  # FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_run2.root"),   # 파일 없음
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
