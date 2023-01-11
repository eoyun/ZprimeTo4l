#include "ZprimeTo4l/MergedLepton/interface/MergedMvaEstimator.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedMvaEstimator.h"

#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

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

double MergedMvaEstimator::computeMva(const pat::ElectronRef& el,
                                      const edm::Handle<EcalRecHitCollection>& EBrecHitHandle,
                                      const edm::Handle<EcalRecHitCollection>& EErecHitHandle,
                                      const reco::GsfTrackRef& addGsfTrk,
                                      bool ignoreEoP2nd) {
  const unsigned int nvar = scale_mean_.size();
  float var[nvar];

  var[mergedElectronVar::dr03HcalDepth1TowerSumEt] = el->dr03HcalDepth1TowerSumEt();
  var[mergedElectronVar::full5x5_E1x5oE5x5] = el->full5x5_e1x5()/el->full5x5_e5x5();
  var[mergedElectronVar::full5x5_r9] = el->full5x5_r9();
  var[mergedElectronVar::absDEtaIn] = std::abs(el->deltaEtaSuperClusterTrackAtVtx());
  var[mergedElectronVar::absDPhiIn] = std::abs(el->deltaPhiSuperClusterTrackAtVtx());
  var[mergedElectronVar::EOverP] = el->eSuperClusterOverP();

  std::vector<std::pair<DetId,float>> hitsAndEnergy;

  for (auto cluster = el->basicClustersBegin(); cluster != el->basicClustersEnd(); ++cluster) {
    for (auto& xtal : (*cluster)->hitsAndFractions()) {
      if ( xtal.first.subdetId() != EcalBarrel && xtal.first.subdetId() != EcalEndcap )
        continue;

      const auto* recHits = xtal.first.subdetId() == EcalBarrel ? &(*EBrecHitHandle) : &(*EErecHitHandle);
      const auto& theHit = recHits->find(xtal.first);

      if ( theHit != recHits->end() )
        hitsAndEnergy.push_back(std::make_pair( xtal.first, xtal.second*theHit->energy() ));
    }
  }

  auto sortByEnergy = [] (const std::pair<DetId,float>& a, const std::pair<DetId,float>& b) {
    return a.second > b.second;
  };

  std::sort(hitsAndEnergy.begin(),hitsAndEnergy.end(),sortByEnergy);

  var[mergedElectronVar::full5x5_r1] = hitsAndEnergy.front().second/el->full5x5_e5x5();
  var[mergedElectronVar::EoverP_1st] = hitsAndEnergy.front().second/el->gsfTrack()->p();

  if (!ignoreEoP2nd)
    var[mergedElectronVar::EoverP_2nd] = hitsAndEnergy.size() > 1 ? hitsAndEnergy.at(1).second/addGsfTrk->p() : 0.;

  for (unsigned int idx = 0; idx < nvar; idx++)
    var[idx] = ( var[idx]-scale_mean_.at(idx) ) / scale_std_.at(idx);

  double rawScore = gbrForest_->GetResponse(var);

  return 1./( 1. + std::exp(-rawScore) ); // return sigmoid
}
