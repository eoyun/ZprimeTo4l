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
    etaSCWidth,
    phiSCWidth,
    full5x5_E1x5oE5x5,
    full5x5_r9,
    dEtaIn,
    dPhiIn,
    EseedOverP,
    fbrem
  };

  MergedMvaEstimator(const edm::FileInPath& weightsfile, const edm::FileInPath& meanstdfile);
  ~MergedMvaEstimator() {}

  double computeMva(const pat::ElectronRef& el);

private:
  std::unique_ptr<const GBRForest> gbrForest_;
  std::vector<double> scale_mean_;
  std::vector<double> scale_std_;
};

#endif
