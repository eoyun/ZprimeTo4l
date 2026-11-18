import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from Configuration.Eras.Era_Run3_cff import Run3

options = VarParsing("analysis")
options.outputFile = "Leptons.root"
options.register("isMC", True, VarParsing.multiplicity.singleton,
                 VarParsing.varType.bool, "Read MC generator, pileup and gen leptons")
options.register("filterProcess", "", VarParsing.multiplicity.singleton,
                 VarParsing.varType.string, "Event-quality TriggerResults process; default PAT for MC, RECO for data")
options.register("minLeptons", 0, VarParsing.multiplicity.singleton,
                 VarParsing.varType.int, "Minimum stored muons + acceptance electrons; 0 disables skim")
options.parseArguments()
if options.minLeptons < 0:
    raise ValueError("minLeptons must be nonnegative")
if not options.inputFiles and options.maxEvents != 0:
    raise ValueError("Supply inputFiles=<MiniAOD> or use maxEvents=0 for the sample-free configuration smoke check")

process = cms.Process("MERGEDLEPTONS", Run3)
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))
process.source = (
    cms.Source("PoolSource", fileNames=cms.untracked.vstring(options.inputFiles))
    if options.inputFiles else cms.Source("EmptySource")
)
process.TFileService = cms.Service("TFileService", fileName=cms.string(options.outputFile))
process.load("ZprimeTo4l.MergedAnalysisRun3.mergedLeptonNtuplizer_cfi")
process.mergedLeptonNtuplizer.isMC = options.isMC
process.mergedLeptonNtuplizer.minLeptons = options.minLeptons
process.mergedLeptonNtuplizer.eventFilters = cms.InputTag(
    "TriggerResults", "", options.filterProcess or ("PAT" if options.isMC else "RECO")
)
process.p = cms.Path(process.mergedLeptonNtuplizer)
# Reads embedded MiniAOD quantities only. No GT, HEEP producers, BDT model,
# momentum/SF payload or external input file is needed for maxEvents=0.
