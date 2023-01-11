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
    dr03HcalDepth1TowerSumEt,
    full5x5_E1x5oE5x5,
    full5x5_r9,
    full5x5_r1,
    absDEtaIn,
    absDPhiIn,
    EOverP,
    EoverP_1st,
    EoverP_2nd
  };

  MergedMvaEstimator(const edm::FileInPath& weightsfile, const edm::FileInPath& meanstdfile);
  ~MergedMvaEstimator() {}

  double computeMva(const pat::ElectronRef& el,
                    const edm::Handle<EcalRecHitCollection>& EBrecHitHandle,
                    const edm::Handle<EcalRecHitCollection>& EErecHitHandle,
                    const reco::GsfTrackRef& addGsfTrk,
                    bool ignoreEoP2nd=false);

private:
  std::unique_ptr<const GBRForest> gbrForest_;
  std::vector<double> scale_mean_;
  std::vector<double> scale_std_;
};

#endif
