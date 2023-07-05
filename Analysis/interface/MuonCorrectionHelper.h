#ifndef MuonCorrectionHelper_h
#define MuonCorrectionHelper_h 1

#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

#include "TRandom3.h"
#include "TFile.h"
#include "TH2D.h"

#include "ZprimeTo4l/Analysis/interface/RoccoR.h"

class MuonCorrectionHelper {
public:
  MuonCorrectionHelper(const edm::FileInPath& rochesterPath, const edm::FileInPath& trigSFpath, const std::string& trigHistName);
  ~MuonCorrectionHelper();

  double rochesterData(const pat::MuonRef& mu) const;
  double rochesterMC(const pat::MuonRef& mu, const edm::Handle<edm::View<reco::GenParticle>>& ahandle);

  double triggerSF(const pat::MuonRef& mu);

private:
  RoccoR rochester_;
  TRandom3 rng_; // workaround when no GEN matching available

  std::unique_ptr<TFile> triggerSFfile_;
  TH2D* triggerSF_;
};

#endif
