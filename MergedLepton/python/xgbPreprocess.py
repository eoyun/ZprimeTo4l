import ROOT
import numpy as np
import pandas as pd
import math

def deltaPhi(phi1, phi2):
    dphi = phi2 - phi1
    if dphi > math.pi:
        dphi -= 2.*math.pi
    elif dphi <= -math.pi:
        dphi += 2.*math.pi

    return dphi

def deltaR2(eta1,phi1,eta2,phi2):
    return abs(eta1-eta2)**2 + abs(deltaPhi(phi1,phi2))**2

class DataframeInitializer(object):
    """docstring for DataframeInitializer."""
    def __init__(self, det, ang):
        super(DataframeInitializer, self).__init__()
        self.det_ = det
        self.ang_ = ang
        self.col_ = [
        'charge',
        'etaSCWidth',
        'phiSCWidth',
        'full5x5_sigmaIetaIeta',
        'full5x5_sigmaIphiIphi',
        'full5x5_E1x5/E5x5',
        'full5x5_E2x5/E5x5',
        'full5x5_r9',
        'dEtaIn',
        'dPhiIn',
        'dPhiSeed',
        'dEtaEle',
        'dPhiEle',
        'dEtaSeed',
        'EseedOverP',
        'EOverP',
        'dxy',
        'dz',
        'fbrem',
        'relModTrkIso',
        'relModEcalHcalD1Iso'
        ]

    def col_names(self):
        return self.col_

    def fill_arr(self,aTree,aGsfTree):
        anArr = np.zeros(shape=(aTree.GetEntries(),len(self.col_)))
        wgts = np.zeros(shape=(aTree.GetEntries(),))

        for iEl, el in enumerate(aTree):
            if self.det_=='EB':
                if abs(el.eta) > 1.5:
                    continue
            elif self.det_=='EE':
                if abs(el.eta) < 1.5:
                    continue
            else:
                raise NameError('Please check EB/EE argument, the current argument is %s' % self.det_)

            aGsfTree.GetEntry(iEl)

            if self.ang_=='dr1':
                if self.det_=='EB':
                    if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < 0.0174**2:
                        continue
                else:
                    if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < (0.00864*abs(math.sinh(el.eta)))**2:
                        continue
            elif self.ang_=='dr2':
                if self.det_=='EB':
                    if not ( deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > 0.0174**2 and
                             deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < (2*0.0174)**2 ):
                        continue
                else:
                    if not ( deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > (0.00864*abs(math.sinh(el.eta)))**2 and
                             deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < (2*0.00864*abs(math.sinh(el.eta)))**2 ):
                        continue
            elif self.ang_=='dr3':
                if self.det_=='EB':
                    if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > (2*0.0174)**2:
                        continue
                else:
                    if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > (2*0.00864*abs(math.sinh(el.eta)))**2:
                        continue
            elif self.ang_=='all':
                pass
            else:
                raise NameError('Please check dR(dr1/dr2/dr3) argument, the current argument is %s' % self.ang_)

            if float(el.passConversionVeto) < 0.5: # veto conversions
                continue

            idx = iEl
            wgts[idx] = el.prefiringweight # xgb does not support negative weights
            anArr[idx,0] = float(el.charge)
            anArr[idx,1] = el.etaSCWidth
            anArr[idx,2] = el.phiSCWidth
            anArr[idx,3] = el.full5x5_sigmaIetaIeta
            anArr[idx,4] = el.full5x5_sigmaIphiIphi
            anArr[idx,5] = el.full5x5_E1x5/el.full5x5_E5x5
            anArr[idx,6] = el.full5x5_E2x5/el.full5x5_E5x5
            anArr[idx,7] = el.full5x5_r9
            anArr[idx,8] = el.dEtaIn
            anArr[idx,9] = el.dPhiIn
            anArr[idx,10] = el.dPhiSeed
            anArr[idx,11] = el.dEtaEle
            anArr[idx,12] = el.dPhiEle
            anArr[idx,13] = el.dEtaSeed
            anArr[idx,14] = el.EseedOverP
            anArr[idx,15] = el.EOverP
            anArr[idx,16] = el.dxy
            anArr[idx,17] = el.dz
            anArr[idx,18] = el.fbrem
            anArr[idx,19] = (el.modTrkIso)/el.pt
            anArr[idx,20] = (el.modEcalIso+el.dr03HcalDepth1TowerSumEt)/el.pt

        # remove rows with only zeros
        mask = ~np.all(anArr==0.0,axis=1)
        anArr = anArr[mask]
        wgts = wgts[mask]

        # sanity check
        if anArr.shape[0] != wgts.shape[0]:
            raise IndexError("Index of the array and the weight does not match!")

        aDf = pd.DataFrame(anArr,columns=self.col_)

        return aDf, wgts

class SampleProcessor(object):
    """docstring for SampleProcessor."""
    def __init__(self, filename, xsec):
        super(SampleProcessor, self).__init__()
        self.filename_ = filename
        self.xsec_ = xsec

    def read(self,dfProducer):
        afile = ROOT.TFile.Open(self.filename_)
        atree = afile.Get("mergedFakeAnalyzer/fake_elTree")
        aGsfTree = afile.Get("mergedFakeAnalyzer/fakeGsf_addGsfTree")
        aWgtHist = afile.Get("mergedFakeAnalyzer/totWeightedSum")
        totWgtSum = aWgtHist.GetBinContent(1)

        df, wgts = dfProducer.fill_arr(atree,aGsfTree)
        wgts = wgts*(self.xsec_*1000./totWgtSum)

        return df, wgts
