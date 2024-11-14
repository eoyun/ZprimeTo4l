import FWCore.ParameterSet.Config as cms
from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('mergedEleBkgFilter',Run3)

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
        #fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Summer22MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/017afd05-e111-41dc-9802-86e708952417.root'), # root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/100000/0100885A-621D-A440-A951-96DE57DDB091.root
        fileNames = cms.untracked.vstring('root://cms-xrd-global.cern.ch//store/mc/Run3Summer22MiniAODv4/DYto2L-4Jets_MLL-50_TuneCP5_13p6TeV_madgraphMLM-pythia8/MINIAODSIM/130X_mcRun3_2022_realistic_v5-v2/2520000/02186c1f-c0c2-46d3-8583-2e9bb486a963.root'), # root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/100000/0100885A-621D-A440-A951-96DE57DDB091.root
    #fileNames = cms.untracked.vstring('file:MiniAODv4_1.root'), # root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL16MiniAODv2/DYJetsToLL_M-50_HT-100to200_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/MINIAODSIM/106X_mcRun2_asymptotic_v17-v2/100000/0100885A-621D-A440-A951-96DE57DDB091.root
    secondaryFileNames = cms.untracked.vstring()
)

process.out = cms.OutputModule("PoolOutputModule",
    compressionAlgorithm = cms.untracked.string('LZMA'),
    compressionLevel = cms.untracked.int32(4),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('MINIAODSIM'),
        filterName = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    eventAutoFlushCompressedSize = cms.untracked.int32(-900),
    fastCloning = cms.untracked.bool(False),
    fileName = cms.untracked.string('MiniAOD_skimmed.root'),
    SelectEvents = cms.untracked.PSet(SelectEvents = cms.vstring("p")),
    outputCommands = process.MINIAODSIMEventContent.outputCommands
)

process.outpath = cms.EndPath(process.out)
process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("124X_mcRun3_2022_realistic_v12")

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi")
process.load("ZprimeTo4l.MergedLepton.MergedEleBkgFilter_cfi")

runVIDmodules = [
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
]

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=False,
                       runVID=True,
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2022-Prompt')

process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
    #elesMiniAOD=cms.InputTag("slimmedElectrons")
)

process.modifiedEcalRecHitIsolationScone2nd = process.ModifiedEcalRecHitIsolationScone.clone(
    #emObjectProducer=cms.InputTag("slimmedElectrons"),
    #addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk")
)

process.p = cms.Path(
    process.ModifiedHEEPIDVarValueMaps +
    process.ModifiedEcalRecHitIsolationScone +
    process.egammaPostRecoSeq +
    process.modifiedHEEPIDVarValueMaps2nd +
    process.modifiedEcalRecHitIsolationScone2nd +
    process.mergedEleBkgFilter
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
