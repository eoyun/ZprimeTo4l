import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_cff import Run3

options = VarParsing('analysis')
options.register('isMC', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'MC or data')
options.register('requireTrigger', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Apply HLT OR')
options.register('requireMETFilters', True, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Apply event quality flags')
options.parseArguments()

process = cms.Process('MERGEDMURUN3', Run3)
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = cms.Source('PoolSource', fileNames=cms.untracked.vstring(options.inputFiles))
process.TFileService = cms.Service('TFileService', fileName=cms.string(options.outputFile))
process.load('ZprimeTo4l.Analysis.mergedMuRun3_cfi')
process.mergedMuRun3.isMC = options.isMC
process.mergedMuRun3.requireTrigger = options.requireTrigger
process.mergedMuRun3.requireMETFilters = options.requireMETFilters
if not options.isMC:
    process.mergedMuRun3.METfilters = cms.InputTag('TriggerResults', '', 'RECO')
process.p = cms.Path(process.mergedMuRun3)
# Uses embedded MiniAOD quantities only; no conditions-dependent correction modules.
