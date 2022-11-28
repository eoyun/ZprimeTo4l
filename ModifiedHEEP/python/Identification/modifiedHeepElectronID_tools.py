import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.Identification.heepElectronID_tools import *

# copied from psetGsfEleTrkPtFall16IsoCut
def psetGsfModifiedEleTrkPtIsoCut(wpEB, wpEE):
    return cms.PSet(
        cutName = cms.string('GsfEleValueMapIsoRhoCut'),
        # Three constants for the GsfEleTrkPtIsoCut
        #     cut = constTerm if value < slopeStart
        #     cut = slopeTerm * (value - slopeStart) + constTerm if value >= slopeStart
        slopeTermEB = cms.double( wpEB.trkIsoSlopeTerm ),
        slopeTermEE = cms.double( wpEE.trkIsoSlopeTerm ),
        slopeStartEB = cms.double( wpEB.trkIsoSlopeStart ),
        slopeStartEE = cms.double( wpEE.trkIsoSlopeStart ),
        constTermEB = cms.double( wpEB.trkIsoConstTerm ),
        constTermEE = cms.double( wpEE.trkIsoConstTerm ),
        #no rho so we zero it out, if the input tag is empty, its ignored anyways
        rhoEtStartEB = cms.double( 999999. ),
        rhoEtStartEE = cms.double( 999999. ),
        rhoEAEB = cms.double( 0. ),
        rhoEAEE = cms.double( 0. ),
        rho = cms.InputTag(""),
        value = cms.InputTag("ModifiedHEEPIDVarValueMaps","eleTrkPtIso"),
        needsAdditionalProducts = cms.bool(True),
        isIgnored = cms.bool(False)
        )

def psetGsfModifiedEleEmHadD1IsoRhoCut(wpEB, wpEE, energyType="EcalTrk"):
    return cms.PSet(
        cutName = cms.string('GsfEleModifiedEmHadD1IsoRhoCut'),
        slopeTermEB = cms.double( wpEB.ehIsoSlopeTerm ),
        slopeTermEE = cms.double( wpEE.ehIsoSlopeTerm ),
        slopeStartEB = cms.double( wpEB.ehIsoSlopeStart ),
        slopeStartEE = cms.double( wpEE.ehIsoSlopeStart ),
        constTermEB = cms.double( wpEB.ehIsoConstTerm ),
        constTermEE = cms.double( wpEE.ehIsoConstTerm ),
        energyType = cms.string( energyType ),
        rhoConstant = cms.double( wpEB.effAreaForEHIso ), # expected to be the same for EB and EE
        rho = cms.InputTag("fixedGridRhoFastjetAll"),
        value = cms.InputTag("ModifiedEcalRecHitIsolationScone","EcalRecHitIso"),
        needsAdditionalProducts = cms.bool(True),
        isIgnored = cms.bool(False)
        )

def configureModifiedHEEPElectronID(idName, wpEB, wpEE):
    parameterSet = cms.PSet(
        idName = cms.string(idName),
        cutFlow = cms.VPSet(
            psetMinPtCut(cutValue=12.),                   #0
            psetGsfEleSCEtaMultiRangeCut(),               #1
            psetGsfEleDPhiInCut(wpEB,wpEE),               #2
            psetGsfEleHadronicOverEMLinearCut(wpEB,wpEE), #3
            psetGsfModifiedEleTrkPtIsoCut(wpEB,wpEE),     #4
            psetGsfModifiedEleEmHadD1IsoRhoCut(wpEB,wpEE),#5
            psetGsfEleDxyCut(wpEB,wpEE),                  #6
            psetGsfEleMissingHitsCut(wpEB,wpEE),          #7,
            psetGsfEleEcalDrivenCut(wpEB,wpEE)            #8
            )
        )

    return parameterSet

def configureHEEPElectronID_V70PtX(idName, wpEB, wpEE):
    parameterSet = cms.PSet(
        idName = cms.string(idName),
        cutFlow = cms.VPSet(
            psetMinPtCut(cutValue=12.),                   #0
            psetGsfEleSCEtaMultiRangeCut(),               #1
            psetGsfEleDEtaInSeedCut(wpEB,wpEE),           #2
            psetGsfEleDPhiInCut(wpEB,wpEE),               #3
            psetGsfEleFull5x5SigmaIEtaIEtaWithSatCut(wpEB,wpEE), #4
            psetGsfEleFull5x5E2x5OverE5x5WithSatCut(wpEB,wpEE),  #5
            psetGsfEleHadronicOverEMLinearCut(wpEB,wpEE), #6
            psetGsfEleTrkPtIsoCut(wpEB,wpEE,useHEEPIso=True),#7
            psetGsfEleEmHadD1IsoRhoCut(wpEB,wpEE),        #8
            psetGsfEleDxyCut(wpEB,wpEE),                  #9
            psetGsfEleMissingHitsCut(wpEB,wpEE),          #10,
            psetGsfEleEcalDrivenCut(wpEB,wpEE)            #11
            )
        )
    return parameterSet
