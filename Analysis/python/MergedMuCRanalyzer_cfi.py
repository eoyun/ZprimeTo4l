import FWCore.ParameterSet.Config as cms

mergedMuCRanalyzer = cms.EDAnalyzer("MergedMuCRanalyzer",
  isMC = cms.bool(True),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  genptc = cms.InputTag("prunedGenParticles"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  srcMET = cms.InputTag("slimmedMETsPuppi"),
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt"),
  ptThres = cms.double(50.),
  ptMuThres = cms.double(20.),
  drThres = cms.double(0.3),
  drThresCR = cms.double(0.6),
  ratioThresLo = cms.double(0.5),
  ratioThresHi = cms.double(1.5)
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016aUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_HIPM_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt")
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2017UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2017_UL_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt")
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
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2018UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  muonIdSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ID.root"),
  idHistName = cms.string("NUM_HighPtID_DEN_TrackerMuons_abseta_pt"),
  idHistNameTrkHighPt = cms.string("NUM_TrkHighPtID_DEN_TrackerMuons_abseta_pt"),
  muonIsoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_ISO.root"),
  isoHistName = cms.string("NUM_LooseRelTkIso_DEN_HighPtIDandIPCut_abseta_pt"),
  isoHistNameTrkHighPt = cms.string("NUM_LooseRelTkIso_DEN_TrkHighPtIDandIPCut_abseta_pt")
)
