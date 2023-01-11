import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.Identification.heepElectronID_tools import *
from ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_tools import *

def configureMergedElectronID(mvaWP, wpEB, wpEE, moduleEra=""):
    pSet = cms.PSet(
        idName = cms.string( mvaWP.idName ),
        cutFlow = cms.VPSet(
            # merged electron MVA ID
            cms.PSet( cutName = cms.string("GsfEleMVACut"), #9
                      mvaCuts = cms.vstring( mvaWP.getCutStrings() ),
                      mvaValueMapName = cms.InputTag( "mergedLeptonIDProducer"+moduleEra+":mvaMergedElectronValues" ),
                      mvaCategoriesMapName = cms.InputTag( "mergedLeptonIDProducer"+moduleEra+":mvaMergedElectronCategories" ),
                      needsAdditionalProducts = cms.bool(True),
                      isIgnored = cms.bool(False) )
            )
        )

    return pSet
