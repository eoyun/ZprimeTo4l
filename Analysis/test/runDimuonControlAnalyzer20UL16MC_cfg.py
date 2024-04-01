import FWCore.ParameterSet.Config as cms

process = cms.Process('DimuonControlAnalyzer')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:MiniAOD.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hists.root')
)

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("106X_mcRun2_asymptotic_v17")

process.load("ZprimeTo4l.Analysis.DimuonControlAnalyzer_cfi")

process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')
process.evtCounter.isMC = cms.bool(True)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.p = cms.Path(
    process.evtCounter+
    process.hltFilter
    process.dimuonControlAnalyzer
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
