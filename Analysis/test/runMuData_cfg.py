import FWCore.ParameterSet.Config as cms

process = cms.Process('MuAnalyzer')

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
        #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/HToAATo4L_H2000A1_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/53014CCB-3FCA-4141-9B71-2DA7DC19F58E.root'
        'root://cms-xrd-global.cern.ch//store/data/Run2018A/SingleMuon/MINIAOD/UL2018_MiniAODv2_GT36-v1/2820000/000EE25A-A8E8-1444-8A0B-0DBEBE5634FB.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hists.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
process.GlobalTag.globaltag = cms.string("106X_dataRun2_v36")
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")

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

process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi")
process.load("ZprimeTo4l.MergedLepton.MergedLeptonIDProducer_cfi")
process.load("ZprimeTo4l.Analysis.MergedEMuCRanalyzer_cfi")
process.load("ZprimeTo4l.Analysis.MergedMuCRanalyzer_cfi")
process.load("ZprimeTo4l.Analysis.ResolvedEMuCRanalyzer_cfi")
process.load("ZprimeTo4l.Analysis.ResolvedMuCRanalyzer_cfi")

process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
    elesMiniAOD=cms.InputTag("slimmedElectrons")
)

process.evtCounter = cms.EDAnalyzer('SimpleEventCounter')
process.evtCounter.isMC = cms.bool(False)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_OldMu100_v*",
    "HLT_TkMu100_v*"
)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

from ZprimeTo4l.Analysis.MergedEMuCRanalyzer_cfi import mergedEMuCRanalyzer20UL18
process.mergedEMuCRanalyzerData = mergedEMuCRanalyzer20UL18.clone(
    isMC = cms.bool(False)
)

from ZprimeTo4l.Analysis.ResolvedEMuCRanalyzer_cfi import resolvedEMuCRanalyzer20UL18
process.resolvedEMuCRanalyzerData = resolvedEMuCRanalyzer20UL18.clone(
    isMC = cms.bool(False)
)

from ZprimeTo4l.Analysis.MergedMuCRanalyzer_cfi import mergedMuCRanalyzer20UL18
process.mergedMuCRanalyzerData = mergedMuCRanalyzer20UL18.clone(
    isMC = cms.bool(False)
)

from ZprimeTo4l.Analysis.ResolvedMuCRanalyzer_cfi import resolvedMuCRanalyzer20UL18
process.resolvedMuCRanalyzerData = resolvedMuCRanalyzer20UL18.clone(
    isMC = cms.bool(False)
)

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           electronColl=cms.InputTag("slimmedElectrons","","PAT"),
                           photonColl=cms.InputTag("slimmedPhotons","","PAT"),
                           isData=True)

process.p = cms.Path(
    process.evtCounter+
    process.hltFilter+
    process.prefiringweight+
    process.ModifiedHEEPIDVarValueMaps+
    process.ModifiedEcalRecHitIsolationScone+
    process.mergedLeptonIDProducer20UL18+
    process.egammaPostRecoSeq+
    process.modifiedHEEPIDVarValueMaps2nd+
    process.mergedMuCRanalyzerData+
    process.mergedEMuCRanalyzerData+
    process.resolvedMuCRanalyzerData+
    process.resolvedEMuCRanalyzerData,
    process.fullPatMetTask
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
