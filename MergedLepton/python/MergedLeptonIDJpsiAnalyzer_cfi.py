import FWCore.ParameterSet.Config as cms

mergedLeptonIDJpsiAnalyzer = cms.EDAnalyzer("MergedLeptonIDJpsiAnalyzer",
  isMC = cms.bool(True),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcMuon = cms.InputTag("slimmedMuons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  beamSpot = cms.InputTag("offlineBeamSpot"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  addPackedCandMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddPackedCand"),
  trkIsoMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleTrkPtIso"),
  dPerpIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dPerpIn"),
  dEtaInSeed2nd = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dEtaInSeed2nd"),
  dPhiInSC2nd = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","dPhiInSC2nd"),
  alphaTrack = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","alphaTrack"),
  alphaCalo = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","alphaCalo"),
  normalizedDParaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","normalizedDParaIn"),
  union5x5covIeIe = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5covIeIe"),
  union5x5covIeIp = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5covIeIp"),
  union5x5dEtaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dEtaIn"),
  union5x5dPhiIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dPhiIn"),
  union5x5Energy = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5Energy"),
  packedPFcand = cms.InputTag("packedPFCandidates"),
  generator = cms.InputTag("generator"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  trigList = cms.vstring(
    "HLT_Mu9_IP6_part*"
  ),
  PUrwgt = cms.FileInPath("ZprimeTo4l/MergedLepton/data/BPH_Mu9_IP6_PUrwgt.root"),
  IPthresTag = cms.double(6.),
  dzThres = cms.double(0.1),
  d0Thres = cms.double(0.06),
  probThres = cms.double(10e-2),
  cosAlpha2dThres = cms.double(0.95),
  ptThresTag = cms.double(9.),
  ptThresK = cms.double(3.5)
)