#ifndef MergedMvaEstimator_H
#define MergedMvaEstimator_H 1

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "CommonTools/MVAUtils/interface/GBRForestTools.h"

#include <iostream>
#include <memory>

class MergedMvaEstimator {
public:
  enum mergedElectronVar {
    charge,
    etaSCWidth,
    phiSCWidth,
    full5x5_sigmaIetaIeta,
    full5x5_sigmaIphiIphi,
    full5x5_E1x5oE5x5,
    full5x5_E2x5oE5x5,
    full5x5_r9,
    dEtaIn,
    dPhiIn,
    dPhiSeed,
    dEtaEle,
    dPhiEle,
    dEtaSeed,
    EseedOverP,
    EOverP,
    dxy,
    dz,
    fbrem,
    relModTrkIso,
    relModEcalHcalD1Iso
  };

  MergedMvaEstimator(const edm::FileInPath& weightsfile, const edm::FileInPath& meanstdfile);
  ~MergedMvaEstimator() {}

  void setHas2ndGsf(const bool flag) { has2ndGsf_ = flag; }

  double computeMva(const edm::Ptr<pat::Electron>& el,
                    const reco::Vertex& pv,
                    const float& trkIso,
                    const float& ecalIso);

private:
  bool has2ndGsf_;
  std::unique_ptr<const GBRForest> gbrForest_;
  std::vector<double> scale_mean_;
  std::vector<double> scale_std_;
};

#endif
