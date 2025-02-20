import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('mergedLeptonIDAnalyzer',Run3)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('FWCore.Modules.printContent_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    # 'file:MiniAOD.root'
    #'/store/data/Run2022C/ParkingSingleMuon0/MINIAOD/PromptReco-v1/000/356/489/00000/4144d69c-dc68-4443-b307-765fdaaef674.root'
    #'/store/data/Run2022C/ParkingSingleMuon0/MINIAOD/PromptReco-v1/000/356/488/00000/05236315-30ed-4c16-9288-2a5b7c8786e0.root'
    'file:MiniAODv4_1.root'
    #"/store/mc/Run3Summer22MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/017afd05-e111-41dc-9802-86e708952417.root"
    ),
    secondaryFileNames = cms.untracked.vstring()
)

#process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True),SkipEvent = cms.untracked.vstring('ProductNotFound') )
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hists_Image_signal.root')
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.printContent.verbose = True

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("124X_mcRun3_2022_realistic_v12")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("Geometry.CommonTopologies.bareGlobalTrackingGeometry_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi")
process.load("ZprimeTo4l.MergedLepton.MergedLeptonIDProducer_cfi")
process.load("ZprimeTo4l.MergedLepton.MergedLeptonIDImageEle_cfi")


runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
    'ZprimeTo4l.MergedLepton.Identification.mergedElectronID_20UL18_cff'
]

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')

process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
    elesMiniAOD=cms.InputTag("slimmedElectrons")
)


process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')
process.evtCounter.isMC = cms.bool(False)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(False)
process.hltFilter.HLTPaths = cms.vstring("HLT_Mu12_IP6*") # HLT_Mu9_IP6_part* # HLT_IsoMu24_v*
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.mergedLeptonIDImageEle.isMC = cms.bool(False)

process.p = cms.Path(
    process.evtCounter+
    #process.hltFilter+
    process.ModifiedHEEPIDVarValueMaps+
    process.ModifiedEcalRecHitIsolationScone+
    process.mergedLeptonIDProducer20UL18+
    process.egammaPostRecoSeq+
    process.modifiedHEEPIDVarValueMaps2nd+
    process.mergedLeptonIDImageEle
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
