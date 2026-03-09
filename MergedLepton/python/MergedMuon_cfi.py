import FWCore.ParameterSet.Config as cms

_muTkIsoDefaultCuts = cms.PSet(
  minPt = cms.double(1.0),
  maxDR = cms.double(0.3),
  minDR = cms.double(0.0),
  minDEta = cms.double(0.005),
  dEta2nd = cms.double(0.005),
  dPhi2nd = cms.double(0.05),
  maxDZ = cms.double(0.1),
  maxDPtPt = cms.double(0.1),
  addTrkMinPt = cms.double(10.0),
  addTrkDR2 = cms.double(0.4),
  addTrkREguard = cms.double(0.001),
  addTrkHoE = cms.double(0.1),
  minHits = cms.int32(8),
  minPixelHits = cms.int32(1),
  allowedQualities = cms.vstring(),
  algosToReject = cms.vstring("jetCoreRegionalStep")
)

mergedMuon = cms.EDAnalyzer("MergedMuon",
  isMC = cms.bool(True),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  trackCands = cms.VInputTag(
    cms.InputTag("packedPFCandidates"),
    cms.InputTag("lostTracks")
  ),
  trackCandsVetos = cms.vstring("NONE", "NONE"),
  muonTkIsoCalc = cms.PSet(
    cuts = _muTkIsoDefaultCuts.clone()
  ),
  packedPFcand = cms.InputTag("packedPFCandidates"),
  genptc = cms.InputTag("prunedGenParticles"),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  trigList = cms.vstring(
    "HLT_Mu9_IP6*",
    "HLT_Mu12_IP6*"
  ),
  PUrwgt = cms.FileInPath("ZprimeTo4l/MergedLepton/data/BPH_Mu9_or_Mu12_IP6_PUrwgt.root"),
  IPthresTag = cms.double(6.),
  dzThres = cms.double(0.1),
  d0Thres = cms.double(0.06),
  probThres = cms.double(10e-2),
  cosAlpha2dThres = cms.double(0.95),
  ptThresTag = cms.double(9.),
  ptThresK = cms.double(3.5),
  imageSize = cms.int32(7),
  ESimageSize = cms.int32(3)
)
