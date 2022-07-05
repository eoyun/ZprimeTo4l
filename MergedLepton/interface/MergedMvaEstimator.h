#ifndef MergedMvaEstimator_H
#define MergedMvaEstimator_H 1

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include <iostream>
#include <memory>

class GBRForest;

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
    fbrem
  };

  MergedMvaEstimator(const edm::FileInPath& weightsfile, const edm::FileInPath& meanstdfile);
  ~MergedMvaEstimator() {}

  double computeMva(const edm::Ptr<pat::Electron>& el, const edm::Ptr<reco::Vertex>& pv);

private:
  std::unique_ptr<const GBRForest> gbrForest_;
  std::vector<double> scale_mean_;
  std::vector<double> scale_std_;
};

#endif
