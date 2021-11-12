import FWCore.ParameterSet.Config as cms

Modified2ndGsfValueMaps = cms.EDProducer("Modified2ndGsfValueMapProducer",
  eleSrc = cms.InputTag("slimmedElectrons"),
  addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
  pvSrc = cms.InputTag("offlineSlimmedPrimaryVertices")
)
