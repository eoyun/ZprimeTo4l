#include "ZprimeTo4l/Analysis/interface/ElectronSystematicsHelper.h"

ElectronSystematicsHelper::ElectronSystematicsHelper(const edm::FileInPath& modHeepSFpath) {
  modHeepSFfile_ = std::make_unique<TFile>(modHeepSFpath.fullPath().c_str(),"READ");
  modHeepSF_ = static_cast<TH2D*>(modHeepSFfile_->Get("EGamma_SF2D"));
}

double ElectronSystematicsHelper::GetModifiedHeepSF(const pat::ElectronRef& aEle) {
  double sf = 1.;

  // if ( aEle->pt() < 35. )
  //   sf = modHeepSF_->GetBinContent( modHeepSF_->FindFixBin(aEle->superCluster()->eta(), aEle->pt()) );

  if ( std::abs(aEle->superCluster()->eta()) < 0.8 )
    sf = modHeepSFmuEB1_;
  else if ( std::abs(aEle->superCluster()->eta()) < 1.4442 )
    sf = modHeepSFmuEB2_;
  else
    sf = modHeepSFmuEE_;

  return sf;
}

std::pair<double,double> ElectronSystematicsHelper::GetModifiedHeepSFcl95UpDn(const pat::ElectronRef& aEle) {
  double sf = GetModifiedHeepSF(aEle);
  double cl95up = 2.;
  double cl95dn = 0.5;

  // if ( aEle->pt() < 35. ) {
  //   cl95up = 2.*modHeepSF_->GetBinError( modHeepSF_->FindFixBin(aEle->superCluster()->eta(), aEle->pt()) );
  //   cl95dn = cl95up;

  TF1* pol = nullptr;
  double upper = 2.;
  double cl95 = 0.5;

  if ( std::abs(aEle->superCluster()->eta()) < 0.8 ) {
    pol = modHeepPolEB1_.get();
    upper = modHeepSFupperEB1_;
    cl95 = modHeepSFcl95EB1_;
  } else if ( std::abs(aEle->superCluster()->eta()) < 1.4442 ) {
    pol = modHeepPolEB2_.get();
    upper = modHeepSFupperEB2_;
    cl95 = modHeepSFcl95EB2_;
  } else {
    pol = modHeepPolEE_.get();
    upper = modHeepSFupperEE_;
    cl95 = modHeepSFcl95EE_;
  }

  if ( aEle->pt() < 35. ) {
    double sfBelow35 = modHeepSF_->GetBinContent( modHeepSF_->FindFixBin(aEle->superCluster()->eta(), aEle->pt()) );
    cl95 = std::hypot( cl95, sfBelow35 - sf );
  }

  double diff = pol->Eval( aEle->pt() ) - sf;
  bool isUp = ( diff > 0. );

  if (isUp) {
    cl95dn = cl95;
    cl95up = std::min( std::hypot(cl95, diff), upper - sf );
  } else {
    cl95up = cl95;
    cl95dn = std::min( std::hypot(cl95, diff), sf );
  }

  return std::make_pair( (sf+cl95up)/sf, (sf-cl95dn)/sf );
}

double ElectronSystematicsHelper::GetMergedEleSF(const pat::ElectronRef& aEle) const {
  if ( aEle->userInt("mvaMergedElectronCategories")==1 )
    return mergedEleSFmuNoTrk_;

  return mergedEleSFmuHasTrk_;
}

std::pair<double,double> ElectronSystematicsHelper::GetMergedEleSFcl95UpDn(const pat::ElectronRef& aEle, const double u5x5Et) {
  double sf = GetMergedEleSF(aEle);
  double cl95up = 2.;
  double cl95dn = 0.5;

  const bool hasTrk = aEle->userInt("mvaMergedElectronCategories")!=1;

  const TF1* pol = hasTrk ? mergedElePolHasTrk_.get() : mergedElePolNoTrk_.get();
  const double upper = hasTrk ? mergedEleSFupperHasTrk_ : mergedEleSFupperNoTrk_;
  const double cl95 = hasTrk ? mergedEleSFcl95HasTrk_ : mergedEleSFcl95NoTrk_;

  double diff = pol->Eval( u5x5Et ) - sf;
  bool isUp = ( diff > 0. );

  if (isUp) {
    cl95dn = cl95;
    cl95up = std::min( std::hypot(cl95, diff), upper - sf );
  } else {
    cl95up = cl95;
    cl95dn = std::min( std::hypot(cl95, diff), sf );
  }

  return std::make_pair( (sf+cl95up)/sf, (sf-cl95dn)/sf );
}

void ElectronSystematicsHelper::SetModifiedHeepSF(const double eb1, const double eb2, const double ee) {
  modHeepSFmuEB1_ = eb1;
  modHeepSFmuEB2_ = eb2;
  modHeepSFmuEE_ = ee;
}

void ElectronSystematicsHelper::SetModifiedHeepSFcl95(const double eb1, const double eb2, const double ee) {
  modHeepSFcl95EB1_ = eb1;
  modHeepSFcl95EB2_ = eb2;
  modHeepSFcl95EE_ = ee;
}

void ElectronSystematicsHelper::SetModifiedHeepSFupper(const double eb1, const double eb2, const double ee) {
  modHeepSFupperEB1_ = eb1;
  modHeepSFupperEB2_ = eb2;
  modHeepSFupperEE_ = ee;
}

void ElectronSystematicsHelper::SetModifiedHeepPol1(const std::string& eb1, const std::string& eb2, const std::string& ee) {
  modHeepPolEB1_ = std::make_unique<TF1>("modHeepPolEB1",eb1.c_str(),35.,2000.);
  modHeepPolEB2_ = std::make_unique<TF1>("modHeepPolEB2",eb2.c_str(),35.,2000.);
  modHeepPolEE_ = std::make_unique<TF1>("modHeepPolEE",ee.c_str(),35.,2000.);
}

void ElectronSystematicsHelper::SetMergedEleSF(const double hasTrk, const double noTrk) {
  mergedEleSFmuHasTrk_ = hasTrk;
  mergedEleSFmuNoTrk_ = noTrk;
}

void ElectronSystematicsHelper::SetMergedEleSFcl95(const double hasTrk, const double noTrk) {
  mergedEleSFcl95HasTrk_ = hasTrk;
  mergedEleSFcl95NoTrk_ = noTrk;
}

void ElectronSystematicsHelper::SetMergedEleSFupper(const double hasTrk, const double noTrk) {
  mergedEleSFupperHasTrk_ = hasTrk;
  mergedEleSFupperNoTrk_ = noTrk;
}

void ElectronSystematicsHelper::SetMergedElePol1(const std::string& hasTrk, const std::string& noTrk) {
  mergedElePolHasTrk_ = std::make_unique<TF1>("mergedElePolHasTrk",hasTrk.c_str(),50.,2000.);
  mergedElePolNoTrk_ = std::make_unique<TF1>("mergedElePolNoTrk",noTrk.c_str(),50.,2000.);
}
