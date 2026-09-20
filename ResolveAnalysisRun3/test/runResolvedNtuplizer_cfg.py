# =============================================================================
# ResolvedNtuplizer 예시 cfg (Run 3)
#
# 이 파일은 "틀"입니다. 아래 [USER] 표시된 곳(GlobalTag, source 파일, era,
# trigger path, METfilter 목록)만 님이 채우면 됩니다. 나머지 producer 체인
# (modified-HEEP ValueMap + egamma post-reco)은 ResolvedNtuplizer 가 읽는
# userFloat("ecalTrkEnergyPostCorr") / electronID("modifiedHeepElectronID") /
# addGsfTrk ValueMap 을 만들기 위해 반드시 필요하므로 그대로 둡니다.
#
# 실행:  cmsRun runResolvedNtuplizer_cfg.py
# 출력:  Events.root  (TTree "Events" + norm/meta TH1)
# =============================================================================
import FWCore.ParameterSet.Config as cms

process = cms.Process('raRun3Ntuple')

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)   # [USER] 테스트 땐 예: 2000 로 줄여도 됨
)

# ---- [USER] 입력 MiniAOD 파일 ----
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        # 'root://cms-xrd-global.cern.ch//store/mc/Run3.../MINIAODSIM/....root'
    '/store/mc/Run3Summer22MiniAODv4/WtoLNu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/03d57826-8e8d-4b11-8683-c75d0310b16f.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# ---- 출력 ntuple ----
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('Events.root')
)

# ---- [USER] Global Tag / Geometry / MagneticField (Run-3 era 에 맞게) ----
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('130X_mcRun3_2022_realistic_v5')   # [USER] 예: '130X_mcRun3_2022_realistic_v5'
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')

# ---- modified-HEEP ValueMap producer (electron ID 입력·추가 GSF track) ----
process.load('ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi')
process.load('ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi')

# ---- VID + egamma post-reco: ecalTrkEnergyPostCorr(보정) + modifiedHeepElectronID 생성 ----
runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
]
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,   # → userFloat("ecalTrkEnergyPostCorr")
                       runVID=True,                 # → electronID/userInt("modifiedHeepElectronID")
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')           # [USER] Run-3 era 에 맞게 (예: 2022-Prompt / 2023-Prompt)

# egamma post-reco 뒤에 slimmedElectrons 로 다시 도는 2nd instance (addGsfTrk 여기서 나옴)
process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
    elesMiniAOD = cms.InputTag('slimmedElectrons')
)

# ---- ntuplizer ----
process.load('ZprimeTo4l.ResolveAnalysisRun3.ResolvedNtuplizer_cfi')
process.resolvedNtuplizer.isMC = cms.bool(True)   # [USER] data 면 False
process.resolvedNtuplizer.trigList = cms.vstring(
    # [USER] Run-3 trigger path (예시) — DAS/HLT menu 에서 확인해 채우기
    # 'HLT_Mu50_v*', 'HLT_CascadeMu100_v*', 'HLT_HighPtTkMu100_v*',
    # 'HLT_DoubleEle33_CaloIdL_MW_v*',
)
process.resolvedNtuplizer.METfilterList = cms.vstring(
    # [USER] 요구할 MET filter 이름들
    # 'Flag_goodVertices', 'Flag_globalSuperTightHalo2016Filter', ...
)
# addGsfTrk 은 cfi 기본값 modifiedHEEPIDVarValueMaps2nd:eleAddGsfTrk 사용 (위 producer 와 일치).

# ---- 실행 순서 (producer 먼저, ntuplizer 마지막) ----
process.p = cms.Path(
    process.ModifiedHEEPIDVarValueMaps
    + process.ModifiedEcalRecHitIsolationScone
    + process.egammaPostRecoSeq
    + process.modifiedHEEPIDVarValueMaps2nd
    + process.resolvedNtuplizer
)
