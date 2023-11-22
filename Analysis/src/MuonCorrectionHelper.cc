#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"
#include "DataFormats/Math/interface/deltaR.h"

MuonCorrectionHelper::MuonCorrectionHelper(const edm::FileInPath& rochesterPath, const edm::FileInPath& trigSFpath, const std::string& trigHistName) {
  rochester_.init(rochesterPath.fullPath());
  triggerSFfile_ = std::make_unique<TFile>(trigSFpath.fullPath().c_str(),"READ");
  triggerSF_ = static_cast<TH2D*>(triggerSFfile_->Get(trigHistName.c_str()));
}

MuonCorrectionHelper::MuonCorrectionHelper(const edm::FileInPath& rochesterPath,
                                           const edm::FileInPath& trigSFpath,
                                           const edm::FileInPath& muonIdSFpath,
                                           const edm::FileInPath& muonIsoSFpath,
                                           const std::string& trigHistName) {
  rochester_.init(rochesterPath.fullPath());
  triggerSFfile_ = std::make_unique<TFile>(trigSFpath.fullPath().c_str(),"READ");
  triggerSF_ = static_cast<TH2D*>(triggerSFfile_->Get(trigHistName.c_str()));

  idSFfile_ = std::make_unique<TFile>(muonIdSFpath.fullPath().c_str(),"READ");
  isoSFfile_ = std::make_unique<TFile>(muonIsoSFpath.fullPath().c_str(),"READ");
}

MuonCorrectionHelper::~MuonCorrectionHelper() {
  triggerSFfile_->Close();
  idSFfile_->Close();
  isoSFfile_->Close();
}

double MuonCorrectionHelper::triggerSF(const pat::MuonRef& mu) {
  double apt = mu->pt() > 200. ? 199.9 : mu->pt();
  return triggerSF_->GetBinContent( triggerSF_->FindBin(mu->eta(),apt) );
}

double MuonCorrectionHelper::idSF(const pat::MuonRef& mu, const std::string& name) {
  double apt = mu->pt() > 120. ? 119.9 : mu->pt();
  TH2D* ahist = static_cast<TH2D*>(idSFfile_->Get(name.c_str()));

  return ahist->GetBinContent( ahist->FindBin(std::abs(mu->eta()),apt) );
}

double MuonCorrectionHelper::isoSF(const pat::MuonRef& mu, const std::string& name) {
  double apt = mu->pt() > 120. ? 119.9 : mu->pt();
  TH2D* ahist = static_cast<TH2D*>(isoSFfile_->Get(name.c_str()));

  return ahist->GetBinContent( ahist->FindBin(std::abs(mu->eta()),apt) );
}

double MuonCorrectionHelper::rochesterData(const pat::MuonRef& mu) const {
  return rochester_.kScaleDT(mu->charge(), mu->pt(), mu->eta(), mu->phi());
}

double MuonCorrectionHelper::rochesterMC(const pat::MuonRef& mu, const edm::Handle<edm::View<reco::GenParticle>>& ahandle) {
  std::vector<reco::GenParticleRef> matched;

  auto sortByClosestPt = [&mu](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
    return std::abs( mu->pt() - a->pt() ) < std::abs( mu->pt() - b->pt() );
  };

  // try prompt muon first
  for (unsigned int idx = 0; idx < ahandle->size(); ++idx) {
    const auto& genptc = ahandle->refAt(idx);

    if ( std::abs(genptc->pdgId())!=13 || genptc->status()!=1 )
      continue;

    if ( genptc->isPromptFinalState() ) {
      if ( reco::deltaR2(mu->eta(),mu->phi(),genptc->eta(),genptc->phi()) < 0.09 )
        matched.push_back( genptc.castTo<reco::GenParticleRef>() );
    }
  }

  // then try any final state muon
  if ( matched.empty() ) {
    for (unsigned int idx = 0; idx < ahandle->size(); ++idx) {
      const auto& genptc = ahandle->refAt(idx);

      if ( std::abs(genptc->pdgId())!=13 || genptc->status()!=1 )
        continue;

      if ( reco::deltaR2(mu->eta(),mu->phi(),genptc->eta(),genptc->phi()) < 0.09 )
        matched.push_back( genptc.castTo<reco::GenParticleRef>() );
    }
  }

  if ( !matched.empty() ) {
    // GEN matched
    std::sort(matched.begin(),matched.end(),sortByClosestPt);

    return rochester_.kSpreadMC(mu->charge(), mu->pt(), mu->eta(), mu->phi(), matched.front()->pt());
  }

  if (!mu->isTrackerMuon()) // ???
    return 1.;

  // no GEN matched muon found
  return rochester_.kSmearMC(mu->charge(),
                             mu->pt(),
                             mu->eta(),
                             mu->phi(),
                             mu->innerTrack()->hitPattern().trackerLayersWithMeasurement(),
                             rng_.Uniform());
}
