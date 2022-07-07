#include "ZprimeTo4l/MergedLepton/interface/MergedMvaEstimator.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedMvaEstimator.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

MergedMvaEstimator::MergedMvaEstimator(const edm::FileInPath& weightsfile, const edm::FileInPath& meanstdfile) {
  gbrForest_ = createGBRForest(weightsfile);
  std::ifstream fstr(meanstdfile.fullPath());

  if (fstr) {
    std::string line;
    std::getline(fstr, line);
    std::istringstream iss_0(line);
    for (std::string part; std::getline(iss_0,part,' '); scale_mean_.emplace_back(std::stod(part))) {}

    std::getline(fstr, line);
    std::istringstream iss_1(line);
    for (std::string part; std::getline(iss_1,part,' '); scale_std_.emplace_back(std::stod(part))) {}

  } else {
    throw cms::Exception("ConfigError") << "Error: cannot open a file " << meanstdfile.fullPath() << std::endl;
  }

}

double MergedMvaEstimator::computeMva(const edm::Ptr<pat::Electron>& el, const reco::Vertex& pv) {
  const unsigned int nvar = scale_mean_.size();
  float var[nvar];

  var[mergedElectronVar::charge] = static_cast<float>(el->charge());
  var[mergedElectronVar::etaSCWidth] = el->superCluster()->etaWidth();
  var[mergedElectronVar::phiSCWidth] = el->superCluster()->phiWidth();
  var[mergedElectronVar::full5x5_sigmaIetaIeta] = el->full5x5_sigmaIetaIeta();
  var[mergedElectronVar::full5x5_sigmaIphiIphi] = el->full5x5_sigmaIphiIphi();
  var[mergedElectronVar::full5x5_E1x5oE5x5] = el->full5x5_e1x5()/el->full5x5_e5x5();
  var[mergedElectronVar::full5x5_E2x5oE5x5] = el->full5x5_e2x5Max()/el->full5x5_e5x5();
  var[mergedElectronVar::full5x5_r9] = el->full5x5_r9();
  var[mergedElectronVar::dEtaIn] = el->deltaEtaSuperClusterTrackAtVtx();
  var[mergedElectronVar::dPhiIn] = el->deltaPhiSuperClusterTrackAtVtx();
  var[mergedElectronVar::dPhiSeed] = el->deltaPhiSeedClusterTrackAtCalo();
  var[mergedElectronVar::dEtaEle] = el->deltaEtaEleClusterTrackAtCalo();
  var[mergedElectronVar::dPhiEle] = el->deltaPhiEleClusterTrackAtCalo();
  var[mergedElectronVar::dEtaSeed] = el->deltaEtaSeedClusterTrackAtVtx();
  var[mergedElectronVar::EseedOverP] = el->eSeedClusterOverP();
  var[mergedElectronVar::EOverP] = el->eSuperClusterOverP();
  var[mergedElectronVar::dxy] = el->gsfTrack()->dxy(pv.position());
  var[mergedElectronVar::dz] = el->gsfTrack()->dz(pv.position());
  var[mergedElectronVar::fbrem] = el->fbrem();

  for (unsigned int idx = 0; idx < nvar; idx++)
    var[idx] = ( var[idx]-scale_mean_.at(idx) ) / scale_std_.at(idx);

  double rawScore = gbrForest_->GetResponse(var);

  return 1./( 1. + std::exp(-rawScore) ); // return sigmoid
}
