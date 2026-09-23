import FWCore.ParameterSet.Config as cms

# ResolvedEleCRanalyzerRun3 — Run-3 (2022 preEE) 4e resolved CR analyzer.
# 물리 로직은 Run-2 ResolvedEleCRanalyzer 와 동일. prefiring 제거, egamma VID-only.
# electron SF/PU/FF payload 는 아직 Run-3 버전 부재 + 파일도 없음 → 전부 비활성(주석).
# electron energy 보정도 OFF → analyzer 에서 raw polarP4 사용.
resolvedEleCRanalyzerRun3 = cms.EDAnalyzer("ResolvedEleCRanalyzerRun3",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
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
  # HLT_DoubleEle33_CaloIdL_MW filter labels (seeded leg / unseeded leg)
  trigFilters = cms.vstring(
    "hltEG33CaloIdLMWPMS2Filter",
    "hltEle33CaloIdLGsfTrkIdVLDPhiFilter"
  ),
  trigUnseededFilters = cms.vstring(
    "hltDiEle33CaloIdLMWPMS2UnseededFilter",
    "hltDiEle33CaloIdLGsfTrkIdVLDPhiUnseededFilter"
  ),
  etThres1 = cms.double(35.),
  etThres2 = cms.double(20.),
  ffSystCL95 = cms.double(0.0),

  # ---- payloads: 전부 비활성 (analyzer 에서 SF/PU/FF/energy 보정 주석 처리됨) ----
  # cms.FileInPath 는 config 로드 시 파일 존재를 검증하므로, 없는 파일들까지 포함해 통째로 주석.
  # Run-3 payload 준비되면: 아래 되살리고 analyzer 의 대응 주석도 해제.
  # TODO(Run-3 2022preEE): EGM Run-3 reco/modHeep/trig SF, FF 재측정, LUM PU json 으로 교체.
  # recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),   # 파일 없음
  # FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_run2.root"),                                        # 파일 없음
  # PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/puWeights20UL16.json"),                                  # 파일 없음
  # PUname = cms.string("Collisions16_UltraLegacy_goldenJSON"),
  # modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016postVFP.root"),    # 파일 없음
  # modHeepSFmuEB1 = cms.double(0.986),
  # modHeepSFmuEB2 = cms.double(1.002),
  # modHeepSFmuEE = cms.double(1.022),
  # modHeepSFcl95EB1 = cms.double(0.02226),
  # modHeepSFcl95EB2 = cms.double(0.03557),
  # modHeepSFcl95EE = cms.double(0.05408),
  # modHeepSFupperEB1 = cms.double(1.098),
  # modHeepSFupperEB2 = cms.double(1.139),
  # modHeepSFupperEE = cms.double(1.344),
  # modHeepPolEB1str = cms.string("2.805*0.000001*x+0.986"),
  # modHeepPolEB2str = cms.string("3.102*0.000001*x+1.001"),
  # modHeepPolEEstr = cms.string("0.0001103*x+1.010"),
  # trigSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_passingDoubleEle33_EGM2D_UL2016postVFP.root"),           # 파일 없음
  # trigUnseededSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_passingDoubleEle33unseeded_EGM2D_UL2016postVFP.root"),  # 파일 없음
  # trigSF = cms.vdouble(0.994,0.999,0.997),
  # trigSFcl95 = cms.vdouble(0.001724,0.001826,0.001514),
  # trigSFupper = cms.vdouble(1.002,1.002,1.001),
  # trigPol = cms.vstring("5.591*0.000001*x+0.994","-5.228*0.000001*x+1.000","2.738*0.000001*x+0.995"),
  # trigUnseededSF = cms.vdouble(1.001,1.001,0.999),
  # trigUnseededSFcl95 = cms.vdouble(0.001553,0.001571,0.001475),
  # trigUnseededSFupper = cms.vdouble(1.002,1.001,1.001),
  # trigUnseededPol = cms.vstring("1.482*0.000001*x+1.001","-5.783*0.000001*x+1.002","-1.141*0.0000001*x+0.999")
)
