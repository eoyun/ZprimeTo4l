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
       ''
    )
)

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

## Geometry and Detector Conditions (needed for a few patTuple production steps)
process.load("Configuration.Geometry.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("TrackingTools/TransientTrack/TransientTrackBuilder_cfi")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = cms.string("80X_mcRun2_asymptotic_2016_TrancheIV_v6")

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

from PhysicsTools.SelectorUtils.tools.vid_id_tools import *

dataFormat = DataFormat.MiniAOD
switchOnVIDElectronIdProducer(process, dataFormat)
ElectronIDs = ['RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff']
for idmod in ElectronIDs:
    setupAllVIDIdsInModule(process,idmod,setupVIDElectronSelection)

newTask = cms.Task()
process.egmGsfElectronIDSequence.associate(newTask)

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
    HEEPVID = cms.InputTag("egmGsfElectronIDs:heepElectronID-HEEPV70"),
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
    process.egmGsfElectronIDSequence*
    process.tree
)
