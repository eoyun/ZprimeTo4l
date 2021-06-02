import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProduction")
isMC = True

## MessageLogger
process.load("FWCore.MessageLogger.MessageLogger_cfi")

## Options and Output Report
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(True)
)
process.MessageLogger.cerr.FwkReport.reportEvery = 1000

'''
from os import listdir
DirPath = '/pnfs/knu.ac.kr/data/cms/store/user/sako/ZZTo4L/ZZTo4L_13TeV_powheg_pythia8/ZZTo4L_13TeV_powheg_pythia8_reMiniAOD/180418_045641/0000/'
InputFilesList = listdir(DirPath)
InputFilesList = ['dcap://cluster142.knu.ac.kr/'+DirPath+x for x in InputFilesList if '.root' in x]
'''

## Source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
       ''#'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv3/ZZTo4L_13TeV_powheg_pythia8/MINIAODSIM/PUMoriond17_94X_mcRun2_asymptotic_v3-v1/100000/42BA2638-E9C6-E811-9BEB-001A649D47FD.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("94X_mcRun2_asymptotic_v3")#80X_mcRun2_asymptotic_2016_TrancheIV_v6

process.TFileService = cms.Service("TFileService",
    fileName = cms.string('Ntuple.root')
)

process.goodOfflinePrimaryVertices = cms.EDFilter("VertexSelector",
    src = cms.InputTag("offlineSlimmedPrimaryVertices"),
    cut = cms.string("!isFake && ndof > 4 && abs(z) < 24 && position.Rho < 2"), # tracksSize() > 3 for the older cut
    filter = cms.bool(True),   # otherwise it won't filter the events, just produce an empty vertex collection.
)

process.noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False),
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

from PhysicsTools.PatUtils.l1ECALPrefiringWeightProducer_cfi import l1ECALPrefiringWeightProducer
process.prefiringweight = l1ECALPrefiringWeightProducer.clone(
    TheJets = cms.InputTag("slimmedJets"), #this should be the slimmedJets collection with up to date JECs !
    DataEra = cms.string("2016BtoH"), #Use 2016BtoH for 2016
    UseJetEMPt = cms.bool(False),
    PrefiringRateSystematicUncty = cms.double(0.2),
    SkipWarnings = False
)

from EgammaUser.EgammaPostRecoTools.EgammaPostRecoTools import setupEgammaPostRecoSeq
setupEgammaPostRecoSeq(process,era='2016-Legacy')

newTask = cms.Task()
process.egammaPostRecoSeq.associate(newTask)

process.load("RecoLocalCalo.EcalRecAlgos.EcalSeverityLevelESProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi")
process.load("ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi")
newTask.add(process.ModifiedHEEPIDVarValueMaps)
newTask.add(process.ModifiedEcalRecHitIsolationScone)

process.tree = cms.EDAnalyzer("Minituplizer",
    isMC = cms.untracked.bool(True),
    TriggerResults = cms.InputTag("TriggerResults","","HLT"),
    TriggerSummary = cms.InputTag("hltTriggerSummaryAOD"),
    Vertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    PU = cms.InputTag("slimmedAddPileupInfo"),
    BeamSpot = cms.InputTag("offlineBeamSpot"),
    PFCandidates = cms.InputTag("packedPFCandidates"),
    Muons = cms.InputTag("slimmedMuons"),
    Electrons = cms.InputTag("slimmedElectrons"),
    Conversions = cms.InputTag("reducedEgamma:reducedConversions"),
    GenParticles = cms.InputTag("prunedGenParticles"),
    Generator = cms.InputTag("generator"),
    GenJets = cms.InputTag("slimmedGenJets"),
    Photons = cms.InputTag("slimmedPhotons"),
    rho = cms.InputTag("fixedGridRhoFastjetAll"),
    miniRho = cms.InputTag("fixedGridRhoFastjetCentralNeutral"),
    Jets = cms.InputTag("slimmedJets"),
    MET = cms.InputTag("slimmedMETs"),
    GsfTracks = cms.InputTag("reducedEgamma:reducedGsfTracks"),
    nrSatCrysMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNrSaturateIn5x5"),
    trkIsoMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
    lostTracks = cms.InputTag("lostTracks"),
    eleTracks = cms.InputTag("lostTracks:eleTracks"),
    addGsfTrkSelMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrkSel"),
    EcalRecHitIsoMap = cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
    noSelectedGsfTrk = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleNoSelectedGsfTrk"),
    addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"),
    KFParameters = cms.PSet(
        maxDistance = cms.double(0.01),
        maxNbrOfIterations = cms.int32(10)
    )
)

process.p = cms.Path(
    process.goodOfflinePrimaryVertices*
    process.prefiringweight*
    process.egammaPostRecoSeq*
    process.tree
)
