import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3


process = cms.Process("MiniAODTrackFitterDemo",Run3)

# --- Standard CMSSW Services & Global Tag ---
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 100 # Print log every 100 events

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

from Configuration.AlCa.GlobalTag import GlobalTag
# You must replace this with the correct GlobalTag for your data/MC.
# Example: process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')
# Or for MC: process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')
# Use the GlobalTag appropriate for your MiniAOD file (e.g., 106X_mcRun2_asymptotic_v17)
#process.GlobalTag.globaltag = '106X_mcRun2_asymptotic_v17' # Example for UL2018 MC, modify for your file!
process.GlobalTag.globaltag = cms.string("124X_mcRun3_2022_realistic_v12")

# --- Input files ---
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # Enter the path to your actual miniAOD file here.
        # Example: '/store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/00000/YOUR_EXAMPLE_FILE.root'
        'file:MiniAODv4_1.root'
    )
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10) # Number of events to process
)

# --- Analyzer ---
process.trackHitAnalyzer = cms.EDAnalyzer("MiniAODTrackFitter",
    packedCandidates = cms.InputTag("packedPFCandidates"),
    dtSegments = cms.InputTag("dt4DSegments"),
    cscSegments = cms.InputTag("cscSegments"),
    rpcRecHits = cms.InputTag("rpcRecHits"),
    ttrhBuilder = cms.string("WithTrackAngle") # Name of the TransientTrackingRecHitBuilder
)

# --- Path ---
process.p = cms.Path(process.trackHitAnalyzer)

# --- Optional: TFileService for output ---
# process.TFileService = cms.Service("TFileService",
# fileName = cms.string("track_fitter_output.root")
# )
