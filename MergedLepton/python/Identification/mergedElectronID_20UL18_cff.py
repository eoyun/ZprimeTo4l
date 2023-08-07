import FWCore.ParameterSet.Config as cms
from RecoEgamma.ElectronIdentification.Identification.mvaElectronID_tools import EleMVA_WP
from RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff import WP_HEEP70_EB,WP_HEEP70_EE
from ZprimeTo4l.MergedLepton.Identification.mergedElectronID_tools import configureMergedElectronID
import copy

WP_modHEEP_EB = copy.deepcopy(WP_HEEP70_EB)
WP_modHEEP_EE = copy.deepcopy(WP_HEEP70_EE)

mvaMergedElectron_20UL18_container = EleMVA_WP(
    idName = "mvaMergedElectron",
    mvaTag = "", # never used since we hacked ValueMap names
    cutCategory0 = "0.595", # HasTrk
    cutCategory1 = "0.625", # NoTrkEt2
)

# let's make a chimera of cutflow & MVA ID
mvaMergedElectron_20UL18 = configureMergedElectronID( mvaMergedElectron_20UL18_container, WP_modHEEP_EB, WP_modHEEP_EE, "20UL18" )
mvaMergedElectron_20UL18.isPOGApproved = cms.untracked.bool(False)
