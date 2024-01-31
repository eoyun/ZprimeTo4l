#ifndef ElectronSystematicsHelper_h
#define ElectronSystematicsHelper_h 1

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TFile.h"
#include "TH2D.h"
#include "TF1.h"
#include "TRandom3.h"

class ElectronSystematicsHelper {
public:
  ElectronSystematicsHelper(const edm::FileInPath& modHeepSFpath);
  ~ElectronSystematicsHelper() = default;

  void SetModifiedHeepSF(const double eb1, const double eb2, const double ee);
  void SetModifiedHeepSFcl95(const double eb1, const double eb2, const double ee);
  void SetModifiedHeepSFupper(const double eb1, const double eb2, const double ee);
  void SetModifiedHeepPol1(const std::string& eb1, const std::string& eb2, const std::string& ee);

  void SetMergedEleSF(const double hasTrk, const double noTrk);
  void SetMergedEleSFcl95(const double hasTrk, const double noTrk);
  void SetMergedEleSFupper(const double hasTrk, const double noTrk);
  void SetMergedElePol1(const std::string& hasTrk, const std::string& noTrk);

  double GetModifiedHeepSF(const pat::ElectronRef& aEle);
  std::pair<double,double> GetModifiedHeepSFcl95UpDn(const pat::ElectronRef& aEle);

  double GetMergedEleSF(const pat::ElectronRef& aEle) const;
  std::pair<double,double> GetMergedEleSFcl95UpDn(const pat::ElectronRef& aEle, const double u5x5Et);

  double mergedEleScale(const pat::ElectronRef& aEle) const;
  double mergedEleSmear(const pat::ElectronRef& aEle, const double u5x5En);

private:
  double modHeepSFmuEB1_ = 0.;
  double modHeepSFmuEB2_ = 0.;
  double modHeepSFmuEE_ = 0.;
  double modHeepSFcl95EB1_ = 0.;
  double modHeepSFcl95EB2_ = 0.;
  double modHeepSFcl95EE_ = 0.;
  double modHeepSFupperEB1_ = 0.;
  double modHeepSFupperEB2_ = 0.;
  double modHeepSFupperEE_ = 0.;

  double mergedEleSFmuHasTrk_ = 0.;
  double mergedEleSFmuNoTrk_ = 0.;
  double mergedEleSFcl95HasTrk_ = 0.;
  double mergedEleSFcl95NoTrk_ = 0.;
  double mergedEleSFupperHasTrk_ = 0.;
  double mergedEleSFupperNoTrk_ = 0.;

  std::unique_ptr<TFile> modHeepSFfile_;
  TH2D* modHeepSF_;

  std::unique_ptr<TF1> modHeepPolEB1_;
  std::unique_ptr<TF1> modHeepPolEB2_;
  std::unique_ptr<TF1> modHeepPolEE_;

  std::unique_ptr<TF1> mergedElePolHasTrk_;
  std::unique_ptr<TF1> mergedElePolNoTrk_;

  TRandom3 rng_ = TRandom3(0);
};

#endif
