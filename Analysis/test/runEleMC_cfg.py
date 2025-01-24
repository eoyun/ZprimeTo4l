import FWCore.ParameterSet.Config as cms
#from Configuration.Eras.Era_Run3_cff import Run3

process = cms.Process('EleAnalyzer')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        # 'file:MiniAOD.root'
        # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/230000/0D0716D0-2537-4443-AFFE-BB3BF90C9E9E.root'
        # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/HToAATo4L_H2000A1_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/53014CCB-3FCA-4141-9B71-2DA7DC19F58E.root'
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/HToAATo4L_H250A1_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/60000/B2B6BECF-EF59-6D43-9C7B-406DA2535C4B.root'
        # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2500000/2BCCF1A8-BBA8-624A-82D0-61928D264815.root'
        #'root://cms-xrd-global.cern.ch//store/user/yeo/HToAATo4L_H2000A10_TuneCP5_13p6TeV-pythia8/Run3Summer22_MiniAODv4/240917_112751/0000/MiniAODv4_1.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hists_el.root')
)

process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("106X_upgrade2018_realistic_v16_L1v1")
#process.GlobalTag.globaltag = cms.string("124X_mcRun3_2022_realistic_v12")

runVIDmodules = [
    'RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff',
    'ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_cff',
    'ZprimeTo4l.MergedLepton.Identification.mergedElectronID_20UL18_cff'
]

from RecoEgamma.EgammaTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,
                       runEnergyCorrections=True,
                       runVID=True,
                       eleIDModules=runVIDmodules,
                       phoIDModules=[],
                       era='2018-UL')

from PhysicsTools.PatUtils.l1PrefiringWeightProducer_cfi import l1PrefiringWeightProducer
process.prefiringweight = l1PrefiringWeightProducer.clone(
    TheJets = cms.InputTag("slimmedJets",processName=cms.InputTag.skipCurrentProcess()), #this should be the slimmedJets collection with up to date JECs
    TheMuons = cms.InputTag('slimmedMuons',processName=cms.InputTag.skipCurrentProcess()),
    ThePhotons = cms.InputTag('slimmedPhotons',processName=cms.InputTag.skipCurrentProcess()),
    DataEraECAL = cms.string("None"), #Use 2016BtoH for 2016
    DataEraMuon = cms.string("20172018"), #Use 2016 for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUnctyECAL = cms.double(0.2),
    PrefiringRateSystematicUnctyMuon = cms.double(0.2)
)

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi")
process.load("ZprimeTo4l.MergedLepton.MergedLeptonIDProducer_cfi")
process.load("ZprimeTo4l.Analysis.MergedEleCRanalyzer_cfi")
process.load("ZprimeTo4l.Analysis.ResolvedEleCRanalyzer_cfi")

process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')
process.evtCounter.isMC = cms.bool(True)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = cms.vstring("HLT_DoubleEle25_CaloIdL_MW_v*")
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
    elesMiniAOD=cms.InputTag("slimmedElectrons")
)

process.p = cms.Path(
    process.evtCounter+
    process.hltFilter+
    process.prefiringweight+
    process.ModifiedHEEPIDVarValueMaps+
    process.ModifiedEcalRecHitIsolationScone+
    process.mergedLeptonIDProducer20UL18+
    process.egammaPostRecoSeq+
    process.modifiedHEEPIDVarValueMaps2nd+
    process.mergedEleCRanalyzer20UL18+
    process.resolvedEleCRanalyzer20UL18
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
