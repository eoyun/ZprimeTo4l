import FWCore.ParameterSet.Config as cms

mergedMuRun3 = cms.EDAnalyzer(
    'mergedMuRun3',
    isMC=cms.bool(True),
    srcMuon=cms.InputTag('slimmedMuons'),
    srcPv=cms.InputTag('offlineSlimmedPrimaryVertices'),
    generator=cms.InputTag('generator'),
    triggerResults=cms.InputTag('TriggerResults', '', 'HLT'),
    METfilters=cms.InputTag('TriggerResults', '', 'PAT'),
    requireTrigger=cms.bool(True),
    triggerPaths=cms.vstring('HLT_Mu50_v*', 'HLT_CascadeMu100_v*', 'HLT_HighPtTkMu100_v*'),
    requireMETFilters=cms.bool(True),  # event quality flags; no MET object or kinematic cut
    METfilterNames=cms.vstring(
        'Flag_goodVertices', 'Flag_globalSuperTightHalo2016Filter',
        'Flag_EcalDeadCellTriggerPrimitiveFilter', 'Flag_BadPFMuonFilter',
        'Flag_BadPFMuonDzFilter', 'Flag_hfNoisyHitsFilter',
        'Flag_eeBadScFilter', 'Flag_ecalBadCalibFilter',
    ),
    minMuonPt=cms.double(20.),
    maxMuonEta=cms.double(2.4),
    requireExactlyThreeMuons=cms.bool(True),
    resolvedRelIsoMax=cms.double(0.1),
    # BDT deliberately unimplemented: bdtScore=-999, bdtValid=False; no ID cut.
)
