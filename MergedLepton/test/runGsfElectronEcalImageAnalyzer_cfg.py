import FWCore.ParameterSet.Config as cms

process = cms.Process("EcalImageAnalysis")

process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(10))
process.source = cms.Source("PoolSource",
    fileNames=cms.untracked.vstring('file:your_input_file.root')
)

process.ecalImageAnalyzer = cms.EDAnalyzer('GsfElectronEcalImageAnalyzer',
    electrons=cms.InputTag("gedGsfElectrons"),
    ecalRecHits=cms.InputTag("ecalRecHit", "EcalRecHitsEE"),
    imageSize=cms.int32(5)
)

process.p = cms.Path(process.ecalImageAnalyzer)

