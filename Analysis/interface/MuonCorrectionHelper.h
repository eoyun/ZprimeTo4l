#ifndef MuonCorrectionHelper_h
#define MuonCorrectionHelper_h 1

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TRandom3.h"
#include "TH2D.h"
#include "TFile.h"

#include "ZprimeTo4l/Analysis/interface/RoccoR.h"
#include "correction.h"

class MuonCorrectionHelper {
public:
  MuonCorrectionHelper(const edm::FileInPath& rochesterPath);
  MuonCorrectionHelper(const edm::FileInPath& rochesterPath,
                       const edm::FileInPath& muonTrigSFpath,
                       const edm::FileInPath& muonIdIsoSFpath,
                       const edm::FileInPath& muonBoostIsoSFpath,
                       const edm::FileInPath& muonRecoSFpath);
  ~MuonCorrectionHelper()=default;

  double rochesterData(const pat::MuonRef& mu) const;
  double rochesterMC(const pat::MuonRef& mu, const edm::Handle<edm::View<reco::GenParticle>>& ahandle);
  double nominalData(const pat::MuonRef&mu) const;
  double nominalMC(const pat::MuonRef& mu, const edm::Handle<edm::View<reco::GenParticle>>& ahandle);
  double altScale(const pat::MuonRef&mu, const std::vector<double>& kb, bool isMC=false) const;
  double smear(const pat::MuonRef& mu, const std::vector<double>& params, const std::vector<double>& factors, bool isMC=false);

  double recoSF(const pat::MuonRef& mu);
  double recoSFsyst(const pat::MuonRef& mu);
  double trigSFtracker(const pat::MuonRef& mu);
  double trigSFtrackerSyst(const pat::MuonRef& mu);
  double trigSFglobal(const pat::MuonRef& mu);
  double trigSFglobalSyst(const pat::MuonRef& mu);
  double highptIdSF(const pat::MuonRef& mu);
  double highptIdSFsyst(const pat::MuonRef& mu);
  double trkHighptIdSF(const pat::MuonRef& mu);
  double trkHighptIdSFsyst(const pat::MuonRef& mu);
  double looseIsoSF(const pat::MuonRef& mu);
  double looseIsoSFsyst(const pat::MuonRef& mu);
  double looseIsoSFtracker(const pat::MuonRef& mu);
  double looseIsoSFtrackerSyst(const pat::MuonRef& mu);

  double boostIsoSF(const pat::MuonRef& mu, const reco::Vertex& pv);
  std::pair<double,double> boostIsoSFupdn(const pat::MuonRef& mu, const reco::Vertex& pv);

  bool checkIso(const pat::MuonRef& mu, const reco::TrackRef& trk, const reco::BeamSpot& bs);

private:
  RoccoR rochester_;
  TRandom3 rng_ = TRandom3(0); // workaround when no GEN matching available

  std::unique_ptr<correction::CorrectionSet> recoSF_;
  std::unique_ptr<correction::CorrectionSet> trigSF_;
  std::unique_ptr<correction::CorrectionSet> idisoSF_;

  std::unique_ptr<TFile> boostIsoFile_;
  TH2D* boostIsoSF_;
  TH2D* boostIsoSFup_;
  TH2D* boostIsoSFdn_;
  TH2D* boostIsoTrackerSF_;
  TH2D* boostIsoTrackerSFup_;
  TH2D* boostIsoTrackerSFdn_;
};

#endif
