import FWCore.ParameterSet.Config as cms

# ResolvedNtuplizer 기본 설정. 실제 dataset/era/trigger path/GT 는 cfg 에서 override.
resolvedNtuplizer = cms.EDAnalyzer('ResolvedNtuplizer',
    isMC            = cms.bool(True),
    srcMuon         = cms.InputTag('slimmedMuons'),
    srcEle          = cms.InputTag('slimmedElectrons'),
    srcPv           = cms.InputTag('offlineSlimmedPrimaryVertices'),
    beamSpot        = cms.InputTag('offlineBeamSpot'),
    # modified-HEEP producer(2nd instance)가 만드는 추가 GSF track ValueMap
    addGsfTrk       = cms.InputTag('modifiedHEEPIDVarValueMaps2nd', 'eleAddGsfTrk'),
    generator       = cms.InputTag('generator'),
    pileupSummary   = cms.InputTag('slimmedAddPileupInfo'),
    triggerResults  = cms.InputTag('TriggerResults', '', 'HLT'),
    triggerObjects  = cms.InputTag('slimmedPatTrigger'),
    METfilters      = cms.InputTag('TriggerResults', '', 'PAT'),
    trigList        = cms.vstring(),   # 사용자 조사 (예: HLT_Mu50_v*, double-e path)
    METfilterList   = cms.vstring(),   # 사용자 조사

    # 저장/skim 기준 (느슨하게 — P/F/3-lepton 보존)
    muPtMinStore    = cms.double(20.),
    muEtaMaxStore   = cms.double(2.4),
    eleEtaMaxStore  = cms.double(2.5),
    eleGapLo        = cms.double(1.4442),
    eleGapHi        = cms.double(1.566),
    minLeptonsSkim  = cms.int32(2),
)
