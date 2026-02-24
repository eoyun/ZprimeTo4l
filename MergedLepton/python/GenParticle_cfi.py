import FWCore.ParameterSet.Config as cms

genParticle = cms.EDAnalyzer("GenParticle",
  isMC = cms.bool(True),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  genptc = cms.InputTag("prunedGenParticles"),
  generator = cms.InputTag("generator"),
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
)
