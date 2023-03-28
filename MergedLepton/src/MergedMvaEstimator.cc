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

double MergedMvaEstimator::computeMva(const pat::ElectronRef& el) {
  const unsigned int nvar = scale_mean_.size();
  float var[nvar];

  var[mergedElectronVar::dr03HcalDepth1TowerSumEt] = el->dr03HcalDepth1TowerSumEt();
  var[mergedElectronVar::full5x5_r9] = el->full5x5_r9();
  var[mergedElectronVar::full5x5_sigmaIetaIeta] = el->full5x5_sigmaIetaIeta();
  var[mergedElectronVar::absDEtaInSeed] = std::abs(el->deltaEtaSeedClusterTrackAtVtx());
  var[mergedElectronVar::absDPhiIn] = std::abs(el->deltaPhiSuperClusterTrackAtVtx());
  var[mergedElectronVar::fbrem] = el->fbrem();
  var[mergedElectronVar::EOverP] = el->eSuperClusterOverP();

  for (unsigned int idx = 0; idx < nvar; idx++)
    var[idx] = ( var[idx]-scale_mean_.at(idx) ) / scale_std_.at(idx);

  double rawScore = gbrForest_->GetResponse(var);

  return 1./( 1. + std::exp(-rawScore) ); // return sigmoid
}
