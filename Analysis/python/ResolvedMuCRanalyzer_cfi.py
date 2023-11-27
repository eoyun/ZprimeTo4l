import FWCore.ParameterSet.Config as cms

resolvedMuCRanalyzer = cms.EDAnalyzer("ResolvedMuCRanalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
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
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16.root"),
  ptThresTrig = cms.double(52.),
  ptThres = cms.double(20.),
  ffSystCL95 = cms.double(0.5)
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
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16.root")
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
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL17.root"),
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
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL18.root"),
)
