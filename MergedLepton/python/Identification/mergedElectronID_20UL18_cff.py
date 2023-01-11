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
    cutCategory0 = "0.666", # DR1Et2EB
    cutCategory1 = "0.490", # DR2Et1EB
    cutCategory2 = "0.202", # DR2Et2EB
    cutCategory3 = "0.553", # DR2Et1EE
    cutCategory4 = "0.467", # DR2Et2EE
    cutCategory5 = "0.802", # bkgEt2EB
)

# let's make a chimera of cutflow & MVA ID
mvaMergedElectron_20UL18 = configureMergedElectronID( mvaMergedElectron_20UL18_container, WP_modHEEP_EB, WP_modHEEP_EE, "20UL18" )
mvaMergedElectron_20UL18.isPOGApproved = cms.untracked.bool(False)
