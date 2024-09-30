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
        # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/HToAATo4L_H2000A1_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/80000/23D8A666-0100-B342-8685-59FD0084097E.root'
        # 'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/HToAATo4L_H750A1_TuneCP5_13TeV-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/60000/3FEF02BA-9FF6-0149-9526-484E16FE6376.root'
        'root://cms-xrd-global.cern.ch//store/mc/RunIISummer20UL18MiniAODv2/WZTo3LNu_TuneCP5_13TeV-amcatnloFXFX-pythia8/MINIAODSIM/106X_upgrade2018_realistic_v16_L1v1-v2/2500000/228C9449-8510-454F-83F0-3E6807CE1C53.root'
    ),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('hists_mu.root')
)

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("106X_upgrade2018_realistic_v16_L1v1")
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

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
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
process.evtCounter.isMC = cms.bool(True)

from HLTrigger.HLTfilters.hltHighLevel_cfi import hltHighLevel
process.hltFilter = hltHighLevel.clone()
process.hltFilter.throw = cms.bool(True)
process.hltFilter.HLTPaths = cms.vstring(
    "HLT_Mu50_v*",
    "HLT_OldMu100_v*",
    "HLT_TkMu100_v*"
)
process.hltFilter.TriggerResultsTag = cms.InputTag("TriggerResults","","HLT")

from PhysicsTools.PatUtils.tools.runMETCorrectionsAndUncertainties import runMetCorAndUncFromMiniAOD
runMetCorAndUncFromMiniAOD(process,
                           electronColl=cms.InputTag("slimmedElectrons","","PAT"),
                           photonColl=cms.InputTag("slimmedPhotons","","PAT"),
                           isData=False)

process.p = cms.Path(
    process.evtCounter+
    process.hltFilter+
    process.fullPatMetSequence+
    process.prefiringweight+
    process.ModifiedHEEPIDVarValueMaps+
    process.ModifiedEcalRecHitIsolationScone+
    process.mergedLeptonIDProducer20UL18+
    process.egammaPostRecoSeq+
    process.modifiedHEEPIDVarValueMaps2nd+
    process.mergedMuCRanalyzer20UL18+
    process.mergedEMuCRanalyzer20UL18+
    process.resolvedMuCRanalyzer20UL18+
    process.resolvedEMuCRanalyzer20UL18
)

# Automatic addition of the customisation function from Configuration.DataProcessing.Utils
from Configuration.DataProcessing.Utils import addMonitoring

#call to customisation function addMonitoring imported from Configuration.DataProcessing.Utils
process = addMonitoring(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion
