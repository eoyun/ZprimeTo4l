import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process("MUONRECO",Run3)

# ---- Load Standard Services ----
process.load("Configuration.StandardSequences.Services_cff")
process.load("FWCore.MessageService.MessageLogger_cfi")

# ---- Set the Global Tag ----
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag.globaltag = cms.string("124X_mcRun3_2022_realistic_v12")

# ---- Load Geometry and Magnetic Field ----
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

# ---- Load Tracking and Muon Reconstruction ----
process.load("Configuration.StandardSequences.Reconstruction_cff")

# ---- Define Input Data ----
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:MiniAODv4_1.root'
        #'root://cms-xrd-global.cern.ch//store/data/Run2018D/SingleMuon/AOD/22Jan2019-v2/10000/FF6E5D32-E8C5-3D45-9A06-5F6BB48F5BC9.root'
    )
)

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(100))

# ---- Full Muon Reconstruction with Track Fitting ----
process.muonRecoWithTrackFit = cms.Sequence(
    process.offlineBeamSpot +         # Beam Spot 생성
    process.siPixelRecHits +          # 픽셀 히트 생성
    process.siStripMatchedRecHits +   # 스트립 히트 생성
    process.initialStepSeeds +        # 트랙 시드 생성
    process.ckfTrackCandidates +      # CKF 트랙 후보 생성
    process.initialStepTracks +       # Kalman Filter를 이용한 트랙 피팅
    process.generalTracks +           # 최적화된 피팅된 트랙
    process.muonlocalreco +           # 뮤온 검출기 로컬 재구성
    process.muonrecoComplete          # Standalone + Global Muon Reconstruction
)

# ---- Define Output Module ----
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("MuonRecoWithTrackFit.root"),
    outputCommands = cms.untracked.vstring("keep *")
)

# ---- Run the Sequence ----
process.p = cms.Path(process.muonRecoWithTrackFit)
process.e = cms.EndPath(process.out)
