import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.Identification.heepElectronID_tools import HEEP_WorkingPoint_V1
from ZprimeTo4l.ModifiedHEEP.Identification.modifiedHeepElectronID_tools import configureHEEPElectronID_V70PtX
from RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff import WP_HEEP70_EB,WP_HEEP70_EE
import copy

idName = "heepElectronID-HEEPV70PtX"
WP_ptxHEEP_EB = copy.deepcopy(WP_HEEP70_EB)
WP_ptxHEEP_EB.idName = str(idName)
WP_ptxHEEP_EE = copy.deepcopy(WP_HEEP70_EE)
WP_ptxHEEP_EE.idName = str(idName)

heepElectronID_HEEPV70PtX = configureHEEPElectronID_V70PtX(idName, WP_ptxHEEP_EB, WP_ptxHEEP_EE)
heepElectronID_HEEPV70PtX.isPOGApproved = cms.untracked.bool(False)
