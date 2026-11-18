import FWCore.ParameterSet.Config as cms

mergedLeptonNtuplizer = cms.EDAnalyzer(
    "MergedLeptonNtuplizer",
    isMC=cms.bool(True),
    srcMuon=cms.InputTag("slimmedMuons"),
    srcElectron=cms.InputTag("slimmedElectrons"),
    srcPV=cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamSpot=cms.InputTag("offlineBeamSpot"),
    rho=cms.InputTag("fixedGridRhoFastjetAll"),
    srcMET=cms.InputTag("slimmedMETs"),
    srcPuppiMET=cms.InputTag("slimmedMETsPuppi"),
    triggerResults=cms.InputTag("TriggerResults", "", "HLT"),
    eventFilters=cms.InputTag("TriggerResults", "", "PAT"),
    triggerObjects=cms.InputTag("slimmedPatTrigger"),
    generator=cms.InputTag("generator"),
    pileup=cms.InputTag("slimmedAddPileupInfo"),
    genParticles=cms.InputTag("prunedGenParticles"),
    # Run-2 Analysis acceptance. Resolved uses TuneP pt; merged uses reco pt.
    muPtMin=cms.double(20.),
    muEtaMax=cms.double(2.4),
    eleEtaMax=cms.double(2.5),
    eleGapLo=cms.double(1.4442),
    eleGapHi=cms.double(1.566),
    # No channel/multiplicity skim by default; also keeps electron-only events.
    minLeptons=cms.uint32(0),
)
