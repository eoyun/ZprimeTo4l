import FWCore.ParameterSet.Config as cms

from ZprimeTo4l.ModifiedHEEP.ModifiedElectronTrackIsolations_cfi import trkIsol03CfgV2

ModifiedHEEPIDVarValueMaps = cms.EDProducer("ModifiedHEEPIDValueMapProducer",
                                    beamSpot=cms.InputTag("offlineBeamSpot"),
                                    ebRecHitsAOD=cms.InputTag("reducedEcalRecHitsEB"),
                                    eeRecHitsAOD=cms.InputTag("reducedEcalRecHitsEE"),
                                    candsAOD=cms.VInputTag("packedCandsForTkIso",
                                                           "lostTracksForTkIso",
                                                           "lostTracksForTkIso:eleTracks"),
                                    #because GsfTracks of electrons are in "packedPFCandidates"
                                    #end KF tracks of electrons are in lostTracks:eleTracks, need to
                                    #tell producer to veto electrons in the first collection
                                    candVetosAOD=cms.vstring("ELES","NONE","NONELES"),
                                    gsfTrksAOD=cms.InputTag("electronGsfTracks"),
                                    elesAOD=cms.InputTag("gedGsfElectrons"),
                                    ebRecHitsMiniAOD=cms.InputTag("reducedEgamma","reducedEBRecHits"),
                                    eeRecHitsMiniAOD=cms.InputTag("reducedEgamma","reducedEERecHits"),
                                    candsMiniAOD=cms.VInputTag("packedPFCandidates",
                                                               "lostTracks",
                                                               "lostTracks:eleTracks"),
                                    candVetosMiniAOD=cms.vstring("ELES","NONE","NONELES"),
                                    gsfTrksMiniAOD=cms.InputTag("reducedEgamma:reducedGsfTracks"),
                                    elesMiniAOD=cms.InputTag("slimmedElectrons"),
                                    dataFormat=cms.int32(0),#0 = auto detection, 1 = AOD, 2 = miniAOD
                                    trkIsoConfig= trkIsol03CfgV2,
                                    trkIso04Config= trkIsol04CfgV2,
                                    makeTrkIso04 = cms.bool(True)
)
