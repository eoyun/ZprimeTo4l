import FWCore.ParameterSet.Config as cms

process = cms.Process("NtupleProduction")
isMC = False

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
       ''#'root://cms-xrd-global.cern.ch//store/data/Run2016G/DoubleEG/MINIAOD/17Jul2018-v1/00000/52D5ADE9-C095-E811-82BF-34E6D7BDDEDB.root'
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("94X_dataRun2_v10")

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
    isMC = cms.untracked.bool(False),
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
    process.egammaPostRecoSeq*
    process.tree
)
