import FWCore.ParameterSet.Config as cms

from RecoEgamma.ElectronIdentification.Identification.heepElectronID_tools import HEEP_WorkingPoint_V1
from ZprimeTo4l.ModifiedHEEPRun3.Identification.modifiedHeepElectronID_tools import configureModifiedHEEPElectronID
from RecoEgamma.ElectronIdentification.Identification.heepElectronID_HEEPV70_cff import WP_HEEP70_EB,WP_HEEP70_EE
import copy

idName = "ModifiedHEEPRun3ID"
WP_modHEEP_EB = copy.deepcopy(WP_HEEP70_EB)
WP_modHEEP_EB.idName = str(idName)
WP_modHEEP_EE = copy.deepcopy(WP_HEEP70_EE)
WP_modHEEP_EE.idName = str(idName)
WP_modHEEP_EE.dEtaInSeedCut = 0.01  # Common EE threshold for dEtaInSeed OR dPerpIn.
WP_modHEEP_EE.full5x5SigmaIEtaIEtaCut = 0.1

ModifiedHEEPRun3ID = configureModifiedHEEPElectronID(idName, WP_modHEEP_EB, WP_modHEEP_EE)
ModifiedHEEPRun3ID.isPOGApproved = cms.untracked.bool(False)
