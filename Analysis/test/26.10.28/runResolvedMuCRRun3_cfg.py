# =============================================================================
# ResolvedMuCRanalyzerRun3 cfg (Run 3, 2022 preEE) — MC/Data 통합
#
# isMC 스위치 하나로 GT/era/analyzer.isMC 분기. egamma 는 VID-only
# (runEnergyCorrections=False) — mu 채널은 electron 을 veto ID 로만 사용.
# L1 timing weight (Run-2 전용) 미포함.
#
# 실행:  cmsRun runResolvedMuCRRun3_cfg.py
# 출력:  hists.root
# =============================================================================
import FWCore.ParameterSet.Config as cms

isMC = True   # Data 로 돌릴 땐 False 로 (아래 GT/METfilters 태그 분기)

# Run-3 era modifier 필수 (Run-3 GT 엔 PCastorRcd 없음 → run3_common 이 CASTOR 제외)
from Configuration.Eras.Era_Run3_cff import Run3
process = cms.Process('ResolvedMuCRanalysisRun3', Run3)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))

process.source = cms.Source('PoolSource',
    fileNames = cms.untracked.vstring(
        # 로컬 테스트용 예시 (Run3Summer22 DY MiniAODSIM). crab/실행 시 override.
        '/store/mc/Run3Summer22MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/017afd05-e111-41dc-9802-86e708952417.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(wantSummary = cms.untracked.bool(True))
process.TFileService = cms.Service('TFileService', fileName = cms.string('hists.root'))

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')
process.load('RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi')

# ---- GlobalTag (2022 preEE) ----
if isMC:
    process.GlobalTag.globaltag = cms.string('130X_mcRun3_2022_realistic_v5')
else:
    # TODO(Run-3 2022preEE Data): prompt/재처리 GT 확정치로 교체.
    process.GlobalTag.globaltag = cms.string('130X_dataRun3_PromptAnalysis_v1')

# ---- VID + egamma post-reco: modifiedHeepElectronID 생성 (electron veto 용) ----
# energy scale/smearing 보정 off → 후단 correctionlib. mu 채널은 energy 미사용.
runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
]
from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')

process.load('ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi')
process.load('ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi')
process.load('ZprimeTo4l.Analysis.ResolvedMuCRanalyzerRun3_cfi')
process.resolvedMuCRanalyzerRun3.isMC = cms.bool(isMC)
if not isMC:
    process.resolvedMuCRanalyzerRun3.METfilters = cms.InputTag('TriggerResults', '', 'RECO')

process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')
process.evtCounter.isMC = cms.bool(isMC)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = cms.vstring('HLT_Mu50_v*', 'HLT_CascadeMu100_v*', 'HLT_HighPtTkMu100_v*')
process.hltFilter.TriggerResultsTag = cms.InputTag('TriggerResults', '', 'HLT')

# MC must reach the analyzer before HLT rejection to fill the raw PU profile.
# The analyzer still applies its own trigger selection after filling PU.
hltSelection = cms.ignore(process.hltFilter) if isMC else process.hltFilter

process.p = cms.Path(
    process.evtCounter+
    hltSelection+
    process.ModifiedHEEPIDVarValueMaps+
    process.ModifiedEcalRecHitIsolationScone+
    process.egammaPostRecoSeq+
    process.resolvedMuCRanalyzerRun3
)
