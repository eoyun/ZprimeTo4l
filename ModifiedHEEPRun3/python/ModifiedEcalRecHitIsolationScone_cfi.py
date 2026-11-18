import FWCore.ParameterSet.Config as cms
import RecoEgamma.EgammaElectronProducers.gedGsfElectrons_cfi as gedGsfEle
from ZprimeTo4l.ModifiedHEEPRun3.ModifiedElectronTrackIsolations_cfi import trkIsol03CfgV2

stdGsfEle = gedGsfEle.gedGsfElectronsTmp

ModifiedEcalRecHitIsolationSconeRun3 = cms.EDProducer("ModifiedEcalRecHitIsolationProducerRun3",
    ecalBarrelRecHitCollection = cms.InputTag("reducedEgamma","reducedEBRecHits"),
    ecalEndcapRecHitCollection = cms.InputTag("reducedEgamma","reducedEERecHits"),

    useNumCrystals = cms.bool(True),
    intRadiusBarrel = cms.double(3.0),
    intRadiusEndcap = cms.double(3.0),
    jurassicWidth = cms.double(1.5), #dEta strip width
    extRadius = cms.double(0.3),
    etMinBarrel = cms.double(0.0),
    eMinBarrel = cms.double(0.095),
    etMinEndcap = cms.double(0.110),
    eMinEndcap = cms.double(0.0),

    useIsolEt = cms.bool(True),
    tryBoth   = cms.bool(True),
    subtract  = cms.bool(False),
    vetoClustered  = cms.bool(False),

    emObjectProducer = cms.InputTag("slimmedElectrons",processName=cms.InputTag.skipCurrentProcess()),

    recHitFlagsExclBarrel = stdGsfEle.recHitFlagsToBeExcludedBarrel,
    recHitFlagsExclEndcaps = stdGsfEle.recHitFlagsToBeExcludedEndcaps,
    recHitSeverityExclBarrel = stdGsfEle.recHitSeverityToBeExcludedBarrel,
    recHitSeverityExclEndcaps = stdGsfEle.recHitSeverityToBeExcludedEndcaps,

    addGsfTrkMap = cms.InputTag("ModifiedHEEPIDVarValueMapsRun3","eleAddGsfTrk"),
    addPackedCandMap = cms.InputTag("ModifiedHEEPIDVarValueMapsRun3","eleAddPackedCand"),
    dEtaInSeed2nd = cms.InputTag("ModifiedHEEPIDVarValueMapsRun3","dEtaInSeed2nd"),
    dPhiInSC2nd = cms.InputTag("ModifiedHEEPIDVarValueMapsRun3","dPhiInSC2nd")
)
