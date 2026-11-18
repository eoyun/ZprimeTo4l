# =============================================================================
# ResolvedNtuplizer cfg (Run 3) — MC 용 (2022 preEE)
#
# test/runResolvedNtuplizer_cfg.py 를 기반으로 한 MC 전용 버전.
# CRAB 로 돌릴 때 source.fileNames 는 crab 가 덮어쓰므로 아래 예시는 로컬 테스트용.
#
# egammaPostRecoSeq 는 유지: runEnergyCorrections=False 라 energy 보정 파일은
# 안 쓰지만(→ 없어도 됨), runVID=True 로 modifiedHeepElectronID 를 electron 에
# embed 해야 ResolvedNtuplizer 가 electronID("modifiedHeepElectronID") 를 읽는다.
#
# 실행:  cmsRun runResolvedNtuplizer_run3MC_cfg.py
# 출력:  Events.root  (TTree "Events" + norm/meta TH1)
# =============================================================================
import FWCore.ParameterSet.Config as cms

# Run-3 era modifier 필수 (run3_common → CaloGeometryBuilder 가 CASTOR 제외; Run-3 GT엔 PCastorRcd 없음)
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('raRun3Ntuple', Run3)

process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)   # 테스트 땐 예: 2000
)

# ---- 입력 MiniAODSIM (crab 가 덮어씀; 로컬 테스트용 예시) ----
process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        #'/store/mc/Run3Summer22MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/017afd05-e111-41dc-9802-86e708952417.root'
        '/store/user/yeo/HToAATo4L_H2000A1_TuneCP5_13p6TeV-pythia8/Run3Summer22EE_MiniAODv4/240917_120305/0000/MiniAODv4_1.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))

# ---- 출력 ntuple ----
process.TFileService = cms.Service('TFileService',
    fileName = cms.string('Events.root')
)

# ---- Global Tag / Geometry / MagneticField (2022 preEE MC) ----
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = cms.string('130X_mcRun3_2022_realistic_v5')   # 2022 preEE MC
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
# energy scale/smearing 보정은 후단 correctionlib 로 처리 → runEnergyCorrections=False
# (2022 energy 보정 파일 안 씀. ntuplizer 도 raw kinematics 저장). runVID=True 만 필요.
runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
    'ZprimeTo4l.ModifiedHEEPRun3.Identification.modifiedHeepElectronID_cff',
]
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,  # scale/smearing 은 후단 correctionlib (보정 파일 불필요)
                       runVID=True,                 # → electronID/userInt("modifiedHeepElectronID")
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')           # 2022 egamma VID era

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
process.resolvedNtuplizer.isMC = cms.bool(True)   # MC
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
