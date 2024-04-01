#include "ZprimeTo4l/Analysis/interface/ElectronSystematicsHelper.h"

ElectronSystematicsHelper::ElectronSystematicsHelper(const edm::FileInPath& modHeepSFpath,
                                                     const edm::FileInPath& recoSFpath) {
  modHeepSFfile_ = std::make_unique<TFile>(modHeepSFpath.fullPath().c_str(),"READ");
  modHeepSF_ = static_cast<TH2D*>(modHeepSFfile_->Get("EGamma_SF2D"));

  recoSFfile_ = std::make_unique<TFile>(recoSFpath.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));
}

ElectronSystematicsHelper::ElectronSystematicsHelper(const edm::FileInPath& modHeepSFpath,
                                                     const edm::FileInPath& recoSFpath,
                                                     const edm::FileInPath& trigSFpath,
                                                     const edm::FileInPath& trigUnseededSFpath) {
  modHeepSFfile_ = std::make_unique<TFile>(modHeepSFpath.fullPath().c_str(),"READ");
  modHeepSF_ = static_cast<TH2D*>(modHeepSFfile_->Get("EGamma_SF2D"));

  recoSFfile_ = std::make_unique<TFile>(recoSFpath.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));

  trigSFfile_ = std::make_unique<TFile>(trigSFpath.fullPath().c_str(),"READ");
  trigUnseededSFfile_ = std::make_unique<TFile>(trigUnseededSFpath.fullPath().c_str(),"READ");

  trigSF_ = static_cast<TH2D*>(trigSFfile_->Get("EGamma_SF2D"));
  trigUnseededSF_ = static_cast<TH2D*>(trigUnseededSFfile_->Get("EGamma_SF2D"));
}

double ElectronSystematicsHelper::GetRecoSF(const pat::ElectronRef& aEle) {
  double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();

  return recoSF_->GetBinContent( recoSF_->FindBin(aEle->superCluster()->eta(), apt ) );
}

double ElectronSystematicsHelper::GetRecoSFerr(const pat::ElectronRef& aEle) {
  double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();

  return recoSF_->GetBinError( recoSF_->FindBin(aEle->superCluster()->eta(), apt ) );
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

double ElectronSystematicsHelper::GetTrigSF(const pat::ElectronRef& aEle) {
  double sf = 1.;

  if ( std::abs(aEle->superCluster()->eta()) < 0.8 )
    sf = trigSFmuEB1_;
  else if ( std::abs(aEle->superCluster()->eta()) < 1.4442 )
    sf = trigSFmuEB2_;
  else
    sf = trigSFmuEE_;

  return sf;
}

double ElectronSystematicsHelper::GetTrigUnseededSF(const pat::ElectronRef& aEle) {
  double sf = 1.;

  if ( std::abs(aEle->superCluster()->eta()) < 0.8 )
    sf = trigUnseededSFmuEB1_;
  else if ( std::abs(aEle->superCluster()->eta()) < 1.4442 )
    sf = trigUnseededSFmuEB2_;
  else
    sf = trigUnseededSFmuEE_;

  return sf;
}

std::pair<double,double> ElectronSystematicsHelper::GetTrigSFcl95UpDn(const pat::ElectronRef& aEle) {
  double sf = GetTrigSF(aEle);
  double cl95up = 2.;
  double cl95dn = 0.5;

  TF1* pol = nullptr;
  double upper = 2.;
  double cl95 = 0.5;

  if ( std::abs(aEle->superCluster()->eta()) < 0.8 ) {
    pol = trigPolEB1_.get();
    upper = trigSFupperEB1_;
    cl95 = trigSFcl95EB1_;
  } else if ( std::abs(aEle->superCluster()->eta()) < 1.4442 ) {
    pol = trigPolEB2_.get();
    upper = trigSFupperEB2_;
    cl95 = trigSFcl95EB2_;
  } else {
    pol = trigPolEE_.get();
    upper = trigSFupperEE_;
    cl95 = trigSFcl95EE_;
  }

  double perBinSF = trigSF_->GetBinContent( trigSF_->FindFixBin(aEle->superCluster()->eta(), aEle->pt()) );
  cl95 = std::hypot( cl95, perBinSF - sf );

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

std::pair<double,double> ElectronSystematicsHelper::GetTrigUnseededSFcl95UpDn(const pat::ElectronRef& aEle) {
  double sf = GetTrigUnseededSF(aEle);
  double cl95up = 2.;
  double cl95dn = 0.5;

  TF1* pol = nullptr;
  double upper = 2.;
  double cl95 = 0.5;

  if ( std::abs(aEle->superCluster()->eta()) < 0.8 ) {
    pol = trigUnseededPolEB1_.get();
    upper = trigUnseededSFupperEB1_;
    cl95 = trigUnseededSFcl95EB1_;
  } else if ( std::abs(aEle->superCluster()->eta()) < 1.4442 ) {
    pol = trigUnseededPolEB2_.get();
    upper = trigUnseededSFupperEB2_;
    cl95 = trigUnseededSFcl95EB2_;
  } else {
    pol = trigUnseededPolEE_.get();
    upper = trigUnseededSFupperEE_;
    cl95 = trigUnseededSFcl95EE_;
  }

  double perBinSF = trigUnseededSF_->GetBinContent( trigUnseededSF_->FindFixBin(aEle->superCluster()->eta(), aEle->pt()) );
  cl95 = std::hypot( cl95, perBinSF - sf );

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

std::pair<double,double> ElectronSystematicsHelper::GetTrigSFMergedUpDn(const pat::ElectronRef& aEle) {
  double sf = GetTrigSF(aEle);
  const auto pair = GetTrigSFcl95UpDn(aEle);

  double errUp = std::hypot( (pair.first-1.)*sf, sf - trigMergedSFmu_);
  errUp = std::hypot(errUp, trigMergedSFerr_);

  double errDn = std::hypot( (pair.second-1.)*sf, sf - trigMergedSFmu_);
  errDn = std::hypot(errDn, trigMergedSFerr_);

  return std::make_pair( (sf+errUp)/sf, (sf-errDn)/sf );
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

void ElectronSystematicsHelper::SetTrigSF(const double eb1, const double eb2, const double ee) {
  trigSFmuEB1_ = eb1;
  trigSFmuEB2_ = eb2;
  trigSFmuEE_ = ee;
}

void ElectronSystematicsHelper::SetTrigSFcl95(const double eb1, const double eb2, const double ee) {
  trigSFcl95EB1_ = eb1;
  trigSFcl95EB2_ = eb2;
  trigSFcl95EE_ = ee;
}

void ElectronSystematicsHelper::SetTrigSFupper(const double eb1, const double eb2, const double ee) {
  trigSFupperEB1_ = eb1;
  trigSFupperEB2_ = eb2;
  trigSFupperEE_ = ee;
}

void ElectronSystematicsHelper::SetTrigPol1(const std::string& eb1, const std::string& eb2, const std::string& ee) {
  trigPolEB1_ = std::make_unique<TF1>("trigPolEB1",eb1.c_str(),28.,2000.);
  trigPolEB2_ = std::make_unique<TF1>("trigPolEB2",eb2.c_str(),28.,2000.);
  trigPolEE_ = std::make_unique<TF1>("trigPolEE",ee.c_str(),28.,2000.);
}

void ElectronSystematicsHelper::SetTrigUnseededSF(const double eb1, const double eb2, const double ee) {
  trigUnseededSFmuEB1_ = eb1;
  trigUnseededSFmuEB2_ = eb2;
  trigUnseededSFmuEE_ = ee;
}

void ElectronSystematicsHelper::SetTrigUnseededSFcl95(const double eb1, const double eb2, const double ee) {
  trigUnseededSFcl95EB1_ = eb1;
  trigUnseededSFcl95EB2_ = eb2;
  trigUnseededSFcl95EE_ = ee;
}

void ElectronSystematicsHelper::SetTrigUnseededSFupper(const double eb1, const double eb2, const double ee) {
  trigUnseededSFupperEB1_ = eb1;
  trigUnseededSFupperEB2_ = eb2;
  trigUnseededSFupperEE_ = ee;
}

void ElectronSystematicsHelper::SetTrigUnseededPol1(const std::string& eb1, const std::string& eb2, const std::string& ee) {
  trigUnseededPolEB1_ = std::make_unique<TF1>("trigUnseededPolEB1",eb1.c_str(),28.,2000.);
  trigUnseededPolEB2_ = std::make_unique<TF1>("trigUnseededPolEB2",eb2.c_str(),28.,2000.);
  trigUnseededPolEE_ = std::make_unique<TF1>("trigUnseededPolEE",ee.c_str(),28.,2000.);
}

void ElectronSystematicsHelper::SetTrigMergedSFandErr(const double sf, const double err) {
  trigMergedSFmu_ = sf;
  trigMergedSFerr_ = err;
}

void ElectronSystematicsHelper::SetAbcdScaleSmear(const double scaleAbove,
                                                  const double scaleBelow,
                                                  const double smear) {
  abcdScaleAbove_ = scaleAbove;
  abcdScaleBelow_ = scaleBelow;
  abcdSmear_ = smear;
}

double ElectronSystematicsHelper::mergedEleScale(const pat::ElectronRef& aEle) const {
  if ( aEle->userInt("mvaMergedElectronCategories")==1 )
    return aEle->userFloat("ecalEnergyPostCorr")/aEle->energy();

  return 0.992;
}

double ElectronSystematicsHelper::mergedEleSmear(const pat::ElectronRef& aEle, const double u5x5En) {
  if ( aEle->userInt("mvaMergedElectronCategories")==1 )
    return aEle->userFloat("ecalEnergyPostCorr")/aEle->energy();

  double corr = u5x5En + rng_.Gaus(0.,u5x5En*0.025);

  return corr/u5x5En;
}

double ElectronSystematicsHelper::GetSingleAbcdScaleSmear(const math::PtEtaPhiMLorentzVector& lvec) {
  double scale = abcdScaleBelow_;

  double smear = lvec.Et()*scale + rng_.Gaus(0.,lvec.Et()*abcdSmear_);

  return smear/lvec.Et();
}

std::pair<double,double> ElectronSystematicsHelper::GetAbcdScaleSmear(const math::PtEtaPhiMLorentzVector& lvec1,
                                                                      const math::PtEtaPhiMLorentzVector& lvec2) {
  const double ptll = (lvec1+lvec2).Pt();
  double scale = 1.;

  if (ptll > 50.)
    scale = abcdScaleAbove_;
  else
    scale = abcdScaleBelow_;

  double smear1 = lvec1.Et()*scale + rng_.Gaus(0.,lvec1.Et()*abcdSmear_);
  double smear2 = lvec2.Et()*scale + rng_.Gaus(0.,lvec2.Et()*abcdSmear_);

  return std::make_pair(smear1/lvec1.Et(),smear2/lvec2.Et());
}
