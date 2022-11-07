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
    def __init__(self, det, ang, etrange):
        super(DataframeInitializer, self).__init__()
        self.det_ = det
        self.ang_ = ang
        self.et_ = etrange
        self.col_ = [
        'etaSCWidth',
        'phiSCWidth', #
        'E1x5/E5x5',
        'R9',
        '|dEtaIn|',
        '|dPhiIn|',
        'Eseed/p',
        'fbrem'
        ]

    def col_names(self):
        return self.col_

    def fill_arr(self,aTree,aGsfTree=None):
        anArr = np.zeros(shape=(aTree.GetEntries(),len(self.col_)))
        wgts = np.zeros(shape=(aTree.GetEntries(),))
        pts = np.zeros(shape=(aTree.GetEntries(),))

        etThres = 50.

        if self.det_=='EB':
            etThres = 200.
        elif self.det_=='EE':
            etThres = 150.
        else:
            raise NameError('Please check EB/EE argument, the current argument is %s' % self.det_)

        for iEl, el in enumerate(aTree):
            if self.det_=='EB':
                if abs(el.eta) > 1.5:
                    continue
            elif self.det_=='EE':
                if abs(el.eta) < 1.5:
                    continue
            else:
                raise NameError('Please check EB/EE argument, the current argument is %s' % self.det_)

            if aGsfTree is not None:
                aGsfTree.GetEntry(iEl)

                if self.ang_=='DR1':
                    if self.det_=='EB':
                        if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.Gsfeta,el.Gsfphi) < 0.0174**2:
                            continue
                    else:
                        if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.Gsfeta,el.Gsfphi) < (0.00864*abs(math.sinh(el.eta)))**2:
                            continue
                elif self.ang_=='DR2':
                    if self.det_=='EB':
                        if not ( deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.Gsfeta,el.Gsfphi) > 0.0174**2 ):
                            continue
                    else:
                        if not ( deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.Gsfeta,el.Gsfphi) > (0.00864*abs(math.sinh(el.eta)))**2 ):
                            continue
                elif self.ang_=='None':
                    pass
                else:
                    raise NameError('Please check dR(DR1/DR2) argument, the current argument is %s' % self.ang_)

            if self.et_=='Et1':
                if not ( el.et > 50. and el.et < etThres ):
                    continue
            elif self.et_=='Et2':
                if not ( el.et >= etThres ):
                    continue
            else:
                raise NameError('check Et argument (Et1/Et2)')

            if float(el.passConversionVeto) < 0.5: # veto conversions
                continue

            idx = iEl
            wgts[idx] = el.weight/abs(el.weight)
            pts[idx] = el.et
            anArr[idx,0] = el.etaSCWidth
            anArr[idx,1] = el.phiSCWidth
            anArr[idx,2] = el.full5x5_E1x5/el.full5x5_E5x5
            anArr[idx,3] = el.full5x5_r9
            anArr[idx,4] = abs(el.dEtaIn)
            anArr[idx,5] = abs(el.dPhiIn)
            anArr[idx,6] = el.EseedOverP
            anArr[idx,7] = el.fbrem

        # remove rows with only zeros
        mask = ~np.all(anArr==0.0,axis=1)
        anArr = anArr[mask]
        wgts = wgts[mask]
        pts = pts[mask]

        # sanity check
        if anArr.shape[0] != wgts.shape[0]:
            raise IndexError("Index of the array and the weight does not match!")

        aDf = pd.DataFrame(anArr,columns=self.col_)

        return aDf, wgts, pts

class SampleProcessor(object):
    """docstring for SampleProcessor."""
    def __init__(self, filename, xsec):
        super(SampleProcessor, self).__init__()
        self.filename_ = filename
        self.xsec_ = xsec
        self.df_ = None   # dataframe that contains ID variables
        self.wgts_ = None # np array of wgts for each electron
        self.pts = None   # np array of Et of each electron

    def read(self,dfProducer):
        afile = ROOT.TFile.Open(self.filename_)
        atree = None
        aGsfTree = None

        if dfProducer.ang_=="None":
            atree = afile.Get("mergedEleBkgMvaInput/bkg_elTree")
        else:
            atree = afile.Get("mergedEleBkgMvaInput/fake_elTree")
            aGsfTree = afile.Get("mergedEleBkgMvaInput/fakeGsf_addGsfTree")

        aWgtHist = afile.Get("mergedEleBkgMvaInput/totWeightedSum")
        totWgtSum = aWgtHist.GetBinContent(1)

        df, wgts, pts = dfProducer.fill_arr(atree,aGsfTree)
        wgts = wgts*(self.xsec_*1000./totWgtSum)

        self.df_ = df
        self.wgts_ = wgts
        self.pts_ = pts

        return
