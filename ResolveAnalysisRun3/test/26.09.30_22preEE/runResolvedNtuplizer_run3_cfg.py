# =============================================================================
# ResolvedNtuplizer cfg (Run 3) — DATA 용 (2022 preEE)
#
# test/runResolvedNtuplizer_cfg.py 를 기반으로 한 data 전용 버전.
# CRAB 로 돌릴 때 source.fileNames 는 crab 가 덮어쓰므로 아래 예시는 로컬 테스트용.
#
# MC 버전(runResolvedNtuplizer_run3MC_cfg.py)과의 차이는 아래 [DATA] 표시 3곳뿐:
#   1) GlobalTag  → data GT (님이 직접 채움)
#   2) isMC       → False
#   3) METfilters InputTag process → data 는 보통 'RECO'/'PAT' 중 파일에 맞게 확인
# egammaPostRecoSeq(runVID=True)는 MC 와 동일하게 유지(modifiedHeepElectronID embed).
# energy 보정은 runEnergyCorrections=False 라 보정 파일 불필요.
#
# 실행:  cmsRun runResolvedNtuplizer_run3_cfg.py
# 출력:  Events.root
# =============================================================================
import FWCore.ParameterSet.Config as cms

from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('raRun3Ntuple', Run3)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)   # 테스트 땐 예: 2000
)

# ---- 입력 MiniAOD (crab 가 덮어씀; 로컬 테스트용 예시) ----
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        # [DATA] 2022 preEE data MiniAOD 예시로 교체
        '/store/data/Run2022C/Muon/MINIAOD/22Sep2023-v1/2520000/21e3c4bc-4477-46d8-ae3d-fc11cea5730a.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# ---- 출력 ntuple ----
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('Events.root')
)

# ---- [DATA] Global Tag / Geometry / MagneticField ----
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('130X_dataRun3_v2')   # [DATA] 예: '130X_dataRun3_v2' / ReReco GT 등 직접 채움
process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')

# ---- modified-HEEP ValueMap producer (electron ID 입력·추가 GSF track) ----
process.load('ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi')
process.load('ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi')
# Independent Run3 producers supply the Run3 VID inputs.
process.load('ZprimeTo4l.ModifiedHEEPRun3.ModifiedHEEPIdVarValueMapProducer_cfi')
process.load('ZprimeTo4l.ModifiedHEEPRun3.ModifiedEcalRecHitIsolationScone_cfi')

# ---- VID + egamma post-reco: modifiedHeepElectronID(VID) 생성 ----
# energy scale/smearing 보정은 후단 correctionlib → runEnergyCorrections=False (보정 파일 불필요)
# runVID=True 만 필요: ntuplizer 가 electronID("modifiedHeepElectronID") 를 읽음.
runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
    'ZprimeTo4l.ModifiedHEEPRun3.Identification.modifiedHeepElectronID_cff',
]
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')

# egamma post-reco 뒤에 slimmedElectrons 로 다시 도는 2nd instance (addGsfTrk 여기서 나옴)
process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
    elesMiniAOD = cms.InputTag('slimmedElectrons')
)

# seed 결정체 gain(12/6/1) ValueMap<int> — correctionlib scale 입력
process.seedGainEle = cms.EDProducer('ElectronSeedGainProducer',
    src = cms.InputTag('slimmedElectrons')
)

# ---- ntuplizer ----
process.load('ZprimeTo4l.ResolveAnalysisRun3.ResolvedNtuplizer_cfi')
process.resolvedNtuplizer.isMC = cms.bool(False)   # [DATA] data
process.resolvedNtuplizer.trigList = cms.vstring(
    'HLT_Mu50_v*', 'HLT_CascadeMu100_v*', 'HLT_HighPtTkMu100_v*',
    'HLT_DoubleEle33_CaloIdL_MW_v*',
)
process.resolvedNtuplizer.METfilterList = cms.vstring(
    'Flag_goodVertices',
    'Flag_globalSuperTightHalo2016Filter',
    'Flag_EcalDeadCellTriggerPrimitiveFilter',
    'Flag_BadPFMuonFilter',
    'Flag_BadPFMuonDzFilter',
    'Flag_hfNoisyHitsFilter',
    'Flag_eeBadScFilter',
    'Flag_ecalBadCalibFilter',
)

# ---- 실행 순서 (producer 먼저, ntuplizer 마지막) ----
process.p = cms.Path(
    process.ModifiedHEEPIDVarValueMaps
    + process.ModifiedEcalRecHitIsolationScone
    + process.ModifiedHEEPIDVarValueMapsRun3
    + process.ModifiedEcalRecHitIsolationSconeRun3
    + process.egammaPostRecoSeq
    + process.modifiedHEEPIDVarValueMaps2nd
    + process.seedGainEle
    + process.resolvedNtuplizer
)
