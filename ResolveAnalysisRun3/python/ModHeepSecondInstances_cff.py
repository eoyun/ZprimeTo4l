import FWCore.ParameterSet.Config as cms

# Second instances of the modified-HEEP ValueMap producers, keyed on the post-reco
# slimmedElectrons (current process) that ResolvedNtuplizer reads.
#
# Why: the first instances (ModifiedHEEPIDVarValueMaps*, ModifiedEcalRecHitIsolationScone*)
# run on slimmedElectrons with skipCurrentProcess() so that VID can evaluate the IDs on the
# input-MiniAOD electrons. egamma post-reco then writes a NEW slimmedElectrons collection with
# the VID results embedded. Those maps cannot be looked up with the new electron refs, so the
# ntuplizer stored -999 for every modified ECAL iso value (schema <= v6).
#
# Requires in the cfg (before calling): the four base cfis loaded, i.e.
#   ZprimeTo4l.ModifiedHEEP.ModifiedHEEPIdVarValueMapProducer_cfi
#   ZprimeTo4l.ModifiedHEEP.ModifiedEcalRecHitIsolationScone_cfi
#   ZprimeTo4l.ModifiedHEEPRun3.ModifiedHEEPIdVarValueMapProducer_cfi
#   ZprimeTo4l.ModifiedHEEPRun3.ModifiedEcalRecHitIsolationScone_cfi
# Insert the returned sequence in the path AFTER process.egammaPostRecoSeq.

def setupModHeepSecondInstances(process, eleSrc='slimmedElectrons'):
    ele = cms.InputTag(eleSrc)

    # old ModifiedHEEP: track-side variables (eleTrkPtIso, dPerpIn, addGsfTrk, ...)
    if not hasattr(process, 'modifiedHEEPIDVarValueMaps2nd'):
        process.modifiedHEEPIDVarValueMaps2nd = process.ModifiedHEEPIDVarValueMaps.clone(
            elesMiniAOD = ele)
    # old ModifiedHEEP: ECAL iso, using the 2nd-instance track association
    process.ModifiedEcalRecHitIsolationScone2nd = process.ModifiedEcalRecHitIsolationScone.clone(
        emObjectProducer = ele,
        addGsfTrkMap     = cms.InputTag('modifiedHEEPIDVarValueMaps2nd', 'eleAddGsfTrk'),
        addPackedCandMap = cms.InputTag('modifiedHEEPIDVarValueMaps2nd', 'eleAddPackedCand'),
        dEtaInSeed2nd    = cms.InputTag('modifiedHEEPIDVarValueMaps2nd', 'dEtaInSeed2nd'),
        dPhiInSC2nd      = cms.InputTag('modifiedHEEPIDVarValueMaps2nd', 'dPhiInSC2nd'))

    # Run3 ModifiedHEEP: track-side variables, then ECAL iso using that association
    process.modifiedHEEPIDVarValueMapsRun32nd = process.ModifiedHEEPIDVarValueMapsRun3.clone(
        elesMiniAOD = ele)
    process.ModifiedEcalRecHitIsolationSconeRun32nd = process.ModifiedEcalRecHitIsolationSconeRun3.clone(
        emObjectProducer = ele,
        addGsfTrkMap     = cms.InputTag('modifiedHEEPIDVarValueMapsRun32nd', 'eleAddGsfTrk'),
        addPackedCandMap = cms.InputTag('modifiedHEEPIDVarValueMapsRun32nd', 'eleAddPackedCand'),
        dEtaInSeed2nd    = cms.InputTag('modifiedHEEPIDVarValueMapsRun32nd', 'dEtaInSeed2nd'),
        dPhiInSC2nd      = cms.InputTag('modifiedHEEPIDVarValueMapsRun32nd', 'dPhiInSC2nd'))

    return cms.Sequence(
        process.modifiedHEEPIDVarValueMaps2nd
        + process.ModifiedEcalRecHitIsolationScone2nd
        + process.modifiedHEEPIDVarValueMapsRun32nd
        + process.ModifiedEcalRecHitIsolationSconeRun32nd)
