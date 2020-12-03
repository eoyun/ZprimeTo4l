import FWCore.ParameterSet.Config as cms
import RecoEgamma.EgammaElectronProducers.gedGsfElectrons_cfi as gedGsfEle
from ZprimeTo4l.ModifiedHEEP.ModifiedElectronTrackIsolations_cfi import trkIsol03CfgV2

stdGsfEle = gedGsfEle.gedGsfElectronsTmp

ModifiedEcalRecHitIsolationScone = cms.EDProducer("ModifiedEcalRecHitIsolationProducer",

    # ecalBarrelRecHitProducer = cms.InputTag("ecalRecHit"),
    # ecalBarrelRecHitCollection = cms.InputTag("reducedEcalRecHitsEB"),
    # ecalEndcapRecHitProducer = cms.InputTag("ecalRecHit"),
    # ecalEndcapRecHitCollection = cms.InputTag("reducedEcalRecHitsEE"),
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

    # emObjectProducer = cms.InputTag("gedGsfElectrons"),
    emObjectProducer = cms.InputTag("slimmedElectrons"),

    recHitFlagsExclBarrel = stdGsfEle.recHitFlagsToBeExcludedBarrel,
    recHitFlagsExclEndcaps = stdGsfEle.recHitFlagsToBeExcludedEndcaps,
    recHitSeverityExclBarrel = stdGsfEle.recHitSeverityToBeExcludedBarrel,
    recHitSeverityExclEndcaps = stdGsfEle.recHitSeverityToBeExcludedEndcaps,

    jurassicWidth2nd = cms.double(2.0),
    intRadius2nd = cms.double(4.0),

    trkIsoConfig = trkIsol03CfgV2,
    # gsfTrks = cms.InputTag("electronGsfTracks")
    gsfTrks = cms.InputTag("reducedEgamma:reducedGsfTracks")
)
