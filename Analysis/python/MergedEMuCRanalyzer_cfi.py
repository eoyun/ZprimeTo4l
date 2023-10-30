import FWCore.ParameterSet.Config as cms

mergedEMuCRanalyzer = cms.EDAnalyzer("MergedEMuCRanalyzer",
  isMC = cms.bool(True),
  srcEle = cms.InputTag("slimmedElectrons"),
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
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  addPackedCandMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddPackedCand"),
  union5x5dEtaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dEtaIn"),
  union5x5Energy = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5Energy"),
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
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16.root"),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL16.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  ptThres = cms.double(50.),
  ptMuThres = cms.double(20.),
  drThres = cms.double(0.3),
  drThresCR = cms.double(0.6),
  ratioThresLo = cms.double(0.5),
  ratioThresHi = cms.double(1.5)
)

mergedEMuCRanalyzer20UL18 = mergedEMuCRanalyzer.clone(
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
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL18.root"),
  MMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MMFF_20UL18.root"),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2018UL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2018_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_OldMu100_or_TkMu100_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root")
)
