import FWCore.ParameterSet.Config as cms

resolvedMuCRanalyzer = cms.EDAnalyzer("ResolvedMuCRanalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  genptc = cms.InputTag("prunedGenParticles"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16.root"),
  ptThresTrig = cms.double(52.),
  ptThres = cms.double(20.),
  ffSystCL95 = cms.double(0.5),
  muScaleBias = cms.vdouble(0.002, 0.049, -0.045, -0.061, 0.016, -0.108,
                            -0.005, -0.064, 0.023, -0.004, 0.041, -0.023,
                            0.026, -0.142, 0.039, 0.000, 0.003, -0.091 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.0102, 6.77e-05, -3.72e-08, 8.53e-12, 0.0129, 6.48e-05, -3.04e-08, 6.63e-12)
)

resolvedMuCRanalyzer20UL16APV = resolvedMuCRanalyzer.clone(
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016aUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2016_preVFP_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2016_preVFP_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2016_preVFP_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16.root"),
  muScaleBias = cms.vdouble(-0.075, -0.053, -0.020, -0.009, -0.010, 0.078,
                            -0.049, -0.100, 0.072, 0.057, -0.036, -0.073,
                            -0.007, 0.003, -0.015, -0.003, 0.060, -0.196 ),
  muSmearFactors = cms.vdouble(0.,0.),
  muSmearParams = cms.vdouble(0.011, 6.87e-05, -3.88e-08, 9.03e-12, 0.013, 6.93e-05, -3.46e-08, 7.72e-12)
)

resolvedMuCRanalyzer20UL17 = resolvedMuCRanalyzer.clone(
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
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2017_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2017_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2017_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL17.root"),
  muScaleBias = cms.vdouble(-0.049, -0.032, 0.027, -0.026, 0.012, -0.151,
                            -0.048, 0.001, 0.039, 0.027, -0.007, 0.058,
                            -0.025, -0.019, -0.049, 0.018, 0.000, 0.004 ),
  muSmearFactors = cms.vdouble(0.,0.3202),
  muSmearParams = cms.vdouble(0.0104, 6.11e-05, -3.31e-08, 6.73e-12, 0.0121, 5.92e-05, -2.61e-08, 5.11e-12)
)

resolvedMuCRanalyzer20UL18 = resolvedMuCRanalyzer.clone(
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
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_HLT_2018_schemaV2.json"),
  muonIdIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_IDISO_2018_schemaV2.json"),
  muonRecoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/ScaleFactors_Muon_highPt_RECO_2018_schemaV2.json"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL18.root"),
  muScaleBias = cms.vdouble(-0.089, -0.005, -0.014, 0.020, 0.042, -0.012,
                            0.057, -0.028, -0.009, 0.007, -0.068, -0.235,
                            0.149, 0.021, -0.028, 0.000, -0.059, 0.025),
  muSmearFactors = cms.vdouble(0.,0.46),
  muSmearParams = cms.vdouble(0.0108, 5.93e-05, -3.08e-08, 6.04e-12, 0.0136, 5.47e-05, -2.3e-08, 4.66e-12)
)
