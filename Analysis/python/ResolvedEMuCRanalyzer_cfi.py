import FWCore.ParameterSet.Config as cms

resolvedEMuCRanalyzer = cms.EDAnalyzer("ResolvedEMuCRanalyzer",
  isMC = cms.untracked.bool(True),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  genptc = cms.InputTag("prunedGenParticles"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  trigList = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  ),
  rochesterPath = cms.FileInPath("ZprimeTo4l/Analysis/data/RoccoR2016bUL.txt"),
  triggerSF = cms.FileInPath("ZprimeTo4l/Analysis/data/Efficiencies_muon_generalTracks_Z_Run2016_UL_SingleMuonTriggers.root"),
  trigHistName = cms.string("NUM_Mu50_or_TkMu50_DEN_CutBasedIdGlobalHighPt_and_TkIsoLoose_eta_pt"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  RMFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/RMFF_20UL16.root"),
  REFFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL16.root"),
  ptThresTrig = cms.double(52.),
  ptThres = cms.double(20.)
)
