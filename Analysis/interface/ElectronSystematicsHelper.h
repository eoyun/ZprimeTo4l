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
  ElectronSystematicsHelper(const edm::FileInPath& modHeepSFpath, const edm::FileInPath& recoSFpath);
  ElectronSystematicsHelper(const edm::FileInPath& modHeepSFpath,
                            const edm::FileInPath& recoSFpath,
                            const edm::FileInPath& trigSFpath,
                            const edm::FileInPath& trigUnseededSFpath);
  ~ElectronSystematicsHelper() = default;

  void SetModifiedHeepSF(const double eb1, const double eb2, const double ee);
  void SetModifiedHeepSFcl95(const double eb1, const double eb2, const double ee);
  void SetModifiedHeepSFupper(const double eb1, const double eb2, const double ee);
  void SetModifiedHeepPol1(const std::string& eb1, const std::string& eb2, const std::string& ee);

  void SetMergedEleSF(const double hasTrk, const double noTrk);
  void SetMergedEleSFcl95(const double hasTrk, const double noTrk);
  void SetMergedEleSFupper(const double hasTrk, const double noTrk);
  void SetMergedElePol1(const std::string& hasTrk, const std::string& noTrk);

  void SetTrigSF(const double eb1, const double eb2, const double ee);
  void SetTrigSFcl95(const double eb1, const double eb2, const double ee);
  void SetTrigSFupper(const double eb1, const double eb2, const double ee);
  void SetTrigPol1(const std::string& eb1, const std::string& eb2, const std::string& ee);

  void SetTrigUnseededSF(const double eb1, const double eb2, const double ee);
  void SetTrigUnseededSFcl95(const double eb1, const double eb2, const double ee);
  void SetTrigUnseededSFupper(const double eb1, const double eb2, const double ee);
  void SetTrigUnseededPol1(const std::string& eb1, const std::string& eb2, const std::string& ee);

  void SetTrigMergedSFandErr(const double sf, const double err);

  void SetAbcdScaleSmear(const double scaleAbove, const double scaleBelow, const double smear);
  double GetSingleAbcdScaleSmear(const math::PtEtaPhiMLorentzVector& lvec);
  std::pair<double,double> GetAbcdScaleSmear(const math::PtEtaPhiMLorentzVector& lvec1,
                                             const math::PtEtaPhiMLorentzVector& lvec2);

  double GetModifiedHeepSF(const pat::ElectronRef& aEle);
  std::pair<double,double> GetModifiedHeepSFcl95UpDn(const pat::ElectronRef& aEle);

  double GetMergedEleSF(const pat::ElectronRef& aEle) const;
  std::pair<double,double> GetMergedEleSFcl95UpDn(const pat::ElectronRef& aEle, const double u5x5Et);

  double GetTrigSF(const pat::ElectronRef& aEle);
  std::pair<double,double> GetTrigSFcl95UpDn(const pat::ElectronRef& aEle);

  double GetTrigUnseededSF(const pat::ElectronRef& aEle);
  std::pair<double,double> GetTrigUnseededSFcl95UpDn(const pat::ElectronRef& aEle);

  std::pair<double,double> GetTrigSFMergedUpDn(const pat::ElectronRef& aEle);

  double GetRecoSF(const pat::ElectronRef& aEle);
  double GetRecoSFerr(const pat::ElectronRef& aEle);

  double mergedEleScale(const pat::ElectronRef& aEle, const bool isMC) const;
  double mergedEleSmear(const pat::ElectronRef& aEle, const double u5x5En, const bool isMC);

  bool isNonHeepEle(const pat::ElectronRef& aEle, const double trkIso, const double dPerpIn) const;

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

  double trigSFmuEB1_ = 0.;
  double trigSFmuEB2_ = 0.;
  double trigSFmuEE_ = 0.;
  double trigSFcl95EB1_ = 0.;
  double trigSFcl95EB2_ = 0.;
  double trigSFcl95EE_ = 0.;
  double trigSFupperEB1_ = 0.;
  double trigSFupperEB2_ = 0.;
  double trigSFupperEE_ = 0.;

  double trigUnseededSFmuEB1_ = 0.;
  double trigUnseededSFmuEB2_ = 0.;
  double trigUnseededSFmuEE_ = 0.;
  double trigUnseededSFcl95EB1_ = 0.;
  double trigUnseededSFcl95EB2_ = 0.;
  double trigUnseededSFcl95EE_ = 0.;
  double trigUnseededSFupperEB1_ = 0.;
  double trigUnseededSFupperEB2_ = 0.;
  double trigUnseededSFupperEE_ = 0.;

  double trigMergedSFmu_ = 0.;
  double trigMergedSFerr_ = 0.;

  double abcdScaleAbove_ = 1.;
  double abcdScaleBelow_ = 1.;
  double abcdSmear_ = 0.;

  std::unique_ptr<TFile> modHeepSFfile_;
  TH2D* modHeepSF_;

  std::unique_ptr<TF1> modHeepPolEB1_;
  std::unique_ptr<TF1> modHeepPolEB2_;
  std::unique_ptr<TF1> modHeepPolEE_;

  std::unique_ptr<TF1> mergedElePolHasTrk_;
  std::unique_ptr<TF1> mergedElePolNoTrk_;

  std::unique_ptr<TFile> recoSFfile_;
  TH2D* recoSF_;

  std::unique_ptr<TFile> trigSFfile_;
  TH2D* trigSF_;

  std::unique_ptr<TF1> trigPolEB1_;
  std::unique_ptr<TF1> trigPolEB2_;
  std::unique_ptr<TF1> trigPolEE_;

  std::unique_ptr<TFile> trigUnseededSFfile_;
  TH2D* trigUnseededSF_;

  std::unique_ptr<TF1> trigUnseededPolEB1_;
  std::unique_ptr<TF1> trigUnseededPolEB2_;
  std::unique_ptr<TF1> trigUnseededPolEE_;

  TRandom3 rng_ = TRandom3(0);
};

#endif
