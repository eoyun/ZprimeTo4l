#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"
#include "DataFormats/Math/interface/deltaR.h"

#include <cmath>

MuonCorrectionHelper::MuonCorrectionHelper(const edm::FileInPath& rochesterPath) {
  rochester_.init(rochesterPath.fullPath());
}

MuonCorrectionHelper::MuonCorrectionHelper(const edm::FileInPath& rochesterPath,
                                           const edm::FileInPath& muonTrigSFpath,
                                           const edm::FileInPath& muonIdIsoSFpath,
                                           const edm::FileInPath& muonBoostIsoSFpath,
                                           const edm::FileInPath& muonRecoSFpath) {
  rochester_.init(rochesterPath.fullPath());

  recoSF_ = std::move(correction::CorrectionSet::from_file(muonRecoSFpath.fullPath()));
  trigSF_ = std::move(correction::CorrectionSet::from_file(muonTrigSFpath.fullPath()));
  idisoSF_ = std::move(correction::CorrectionSet::from_file(muonIdIsoSFpath.fullPath()));

  boostIsoFile_ = std::make_unique<TFile>(muonBoostIsoSFpath.fullPath().c_str(),"READ");
  boostIsoSF_ = (TH2D*)boostIsoFile_->Get("passHighPt_2D");
  boostIsoSFup_ = (TH2D*)boostIsoFile_->Get("passHighPt_2D_hi");
  boostIsoSFdn_ = (TH2D*)boostIsoFile_->Get("passHighPt_2D_lo");
  boostIsoTrackerSF_ = (TH2D*)boostIsoFile_->Get("passTrkHighPt_2D");
  boostIsoTrackerSFup_ = (TH2D*)boostIsoFile_->Get("passTrkHighPt_2D_hi");
  boostIsoTrackerSFdn_ = (TH2D*)boostIsoFile_->Get("passTrkHighPt_2D_lo");
}

double MuonCorrectionHelper::boostIsoSF(const pat::MuonRef& mu, const reco::Vertex& pv) {
  TH2D* sf = nullptr;

  if (muon::isHighPtMuon(*mu,pv))
    sf = boostIsoSF_;
  else if (muon::isTrackerHighPtMuon(*mu,pv))
    sf = boostIsoTrackerSF_;
  else
    return 1.;

  int ibin = sf->FindFixBin(std::abs(mu->tunePMuonBestTrack()->eta()),std::min(mu->tunePMuonBestTrack()->pt(),999.));

  return sf->GetBinContent(ibin);
}

std::pair<double,double> MuonCorrectionHelper::boostIsoSFupdn(const pat::MuonRef& mu, const reco::Vertex& pv) {
  TH2D* up = nullptr;
  TH2D* dn = nullptr;
  const double sf = boostIsoSF(mu,pv);

  if (muon::isHighPtMuon(*mu,pv)) {
    up = boostIsoSFup_;
    dn = boostIsoSFdn_;
  } else if (muon::isTrackerHighPtMuon(*mu,pv)) {
    up = boostIsoTrackerSFup_;
    dn = boostIsoTrackerSFdn_;
  } else
    return std::make_pair(0.,0.);

  int ibin = up->FindFixBin(std::abs(mu->tunePMuonBestTrack()->eta()),std::min(mu->tunePMuonBestTrack()->pt(),999.));

  return std::make_pair((sf+up->GetBinContent(ibin))/sf,(sf-dn->GetBinContent(ibin))/sf);
}

bool MuonCorrectionHelper::checkIso(const pat::MuonRef& mu,
                                    const reco::TrackRef& trk,
                                    const reco::BeamSpot& bs) {
  const double dr2 = reco::deltaR2(mu->eta(),mu->phi(),trk->eta(),trk->phi());

  if (dr2 > 0.09 || dr2 < 0.0001)
    return false;

  if ( std::abs(mu->vz() - trk->vz()) > 0.2 )
    return false;

  if ( trk->dxy( bs.position() ) > 0.1 )
    return false;

  return true;
}

void MuonCorrectionHelper::modifiedIso(std::vector<pat::MuonRef>& isolatedHighPtMuons,
                                       std::vector<pat::MuonRef>& isolatedHighPtTrackerMuons,
                                       std::vector<pat::MuonRef>& boostedMuons,
                                       const std::vector<pat::MuonRef>& highPtMuons,
                                       const std::vector<pat::MuonRef>& highPtTrackerMuons,
                                       const reco::BeamSpot& bs) {
  auto sortByTuneP = [] (const pat::MuonRef& a, const pat::MuonRef& b) {
    return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
  };

  for (unsigned idx = 0; idx < highPtMuons.size(); idx++) {
    const auto& firstMuon = highPtMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 ) {
      isolatedHighPtMuons.push_back(firstMuon);

      if (!isoCands.empty())
        boostedMuons.push_back(firstMuon);
    }
  }

  for (unsigned idx = 0; idx < highPtTrackerMuons.size(); idx++) {
    const auto& firstMuon = highPtTrackerMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 ) {
      isolatedHighPtTrackerMuons.push_back(firstMuon);

      if (!isoCands.empty())
        boostedMuons.push_back(firstMuon);
    }
  }
}

void MuonCorrectionHelper::nonHighPtMuonIso(std::vector<pat::MuonRef>& nonHighPtMuons,
                                            std::vector<pat::MuonRef>& nonHighPtMuonsVLiso,
                                            std::map<pat::MuonRef,float>& nonHighPtIsos,
                                            const edm::Handle<edm::View<pat::Muon>>& muonHandle,
                                            const std::vector<pat::MuonRef>& allHighPtMuons,
                                            const std::vector<pat::MuonRef>& highPtMuons,
                                            const std::vector<pat::MuonRef>& highPtTrackerMuons,
                                            const reco::Vertex& pv,
                                            const reco::BeamSpot& bs,
                                            const double ptThres) {
  auto sortByTuneP = [] (const pat::MuonRef& a, const pat::MuonRef& b) {
    return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
  };

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( !aMuon->isTrackerMuon() || aMuon->track().isNull() )
      continue;

    if ( aMuon->tunePMuonBestTrack()->pt() < ptThres || std::abs(aMuon->eta()) > 2.4 )
      continue;

    const auto castMu = aMuon.castTo<pat::MuonRef>();

    if ( std::find( allHighPtMuons.begin(), allHighPtMuons.end(), castMu ) != allHighPtMuons.end() )
      continue;

    bool hits = aMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
                && aMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

    bool ip = std::abs(aMuon->innerTrack()->dxy(pv.position())) < 0.2
              && std::abs(aMuon->innerTrack()->dz(pv.position())) < 0.5;

    if ( aMuon->numberOfMatchedStations() >= 1 && hits && ip )
      nonHighPtMuons.push_back(castMu);
  }

  for (unsigned idx = 0; idx < nonHighPtMuons.size(); idx++) {
    const auto& firstMuon = nonHighPtMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < nonHighPtMuons.size(); jdx++) {
      const auto& secondMuon = nonHighPtMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if ( checkIso(firstMuon,secondMuon->innerTrack(),bs) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    nonHighPtIsos.insert(std::pair<pat::MuonRef,float>(firstMuon,firstIso));

    if (firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.5)
      nonHighPtMuonsVLiso.push_back(firstMuon);
  }
}

double MuonCorrectionHelper::recoSF(const pat::MuonRef& mu) {
  double momentum = std::min(std::max(mu->tunePMuonBestTrack()->p(),50.001),3500.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return recoSF_->at("NUM_GlobalMuons_DEN_TrackerMuonProbes")->evaluate({std::abs(eta),momentum,"nominal"});
}

double MuonCorrectionHelper::recoSFsyst(const pat::MuonRef& mu) {
  double momentum = std::min(std::max(mu->tunePMuonBestTrack()->p(),50.001),3500.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = recoSF_->at("NUM_GlobalMuons_DEN_TrackerMuonProbes")->evaluate({std::abs(eta),momentum,"stat"});
  double syst = recoSF_->at("NUM_GlobalMuons_DEN_TrackerMuonProbes")->evaluate({std::abs(eta),momentum,"syst"});

  return std::hypot(stat,syst);
}

double MuonCorrectionHelper::trigSFtracker(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return trigSF_->at("NUM_HLT_DEN_TrkHighPtLooseRelIsoProbes")->evaluate({std::abs(eta),pt,"nominal"});
}

double MuonCorrectionHelper::trigSFtrackerSyst(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = trigSF_->at("NUM_HLT_DEN_TrkHighPtLooseRelIsoProbes")->evaluate({std::abs(eta),pt,"stat"});
  double syst = trigSF_->at("NUM_HLT_DEN_TrkHighPtLooseRelIsoProbes")->evaluate({std::abs(eta),pt,"syst"});

  return std::hypot(stat,syst);
}

double MuonCorrectionHelper::trigSFglobal(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return trigSF_->at("NUM_HLT_DEN_HighPtLooseRelIsoProbes")->evaluate({std::abs(eta),pt,"nominal"});
}

double MuonCorrectionHelper::trigSFglobalSyst(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = trigSF_->at("NUM_HLT_DEN_HighPtLooseRelIsoProbes")->evaluate({std::abs(eta),pt,"stat"});
  double syst = trigSF_->at("NUM_HLT_DEN_HighPtLooseRelIsoProbes")->evaluate({std::abs(eta),pt,"syst"});

  return std::hypot(stat,syst);
}

double MuonCorrectionHelper::highptIdSF(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return idisoSF_->at("NUM_HighPtID_DEN_GlobalMuonProbes")->evaluate({std::abs(eta),pt,"nominal"});
}

double MuonCorrectionHelper::highptIdSFsyst(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = idisoSF_->at("NUM_HighPtID_DEN_GlobalMuonProbes")->evaluate({std::abs(eta),pt,"stat"});
  double syst = idisoSF_->at("NUM_HighPtID_DEN_GlobalMuonProbes")->evaluate({std::abs(eta),pt,"syst"});

  return std::hypot(stat,syst);
}

double MuonCorrectionHelper::trkHighptIdSF(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return idisoSF_->at("NUM_TrkHighPtID_DEN_GlobalMuonProbes")->evaluate({std::abs(eta),pt,"nominal"});
}

double MuonCorrectionHelper::trkHighptIdSFsyst(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = idisoSF_->at("NUM_TrkHighPtID_DEN_GlobalMuonProbes")->evaluate({std::abs(eta),pt,"stat"});
  double syst = idisoSF_->at("NUM_TrkHighPtID_DEN_GlobalMuonProbes")->evaluate({std::abs(eta),pt,"syst"});

  return std::hypot(stat,syst);
}

double MuonCorrectionHelper::looseIsoSF(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return idisoSF_->at("NUM_probe_LooseRelTkIso_DEN_HighPtProbes")->evaluate({std::abs(eta),pt,"nominal"});
}

double MuonCorrectionHelper::looseIsoSFsyst(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = idisoSF_->at("NUM_probe_LooseRelTkIso_DEN_HighPtProbes")->evaluate({std::abs(eta),pt,"stat"});
  double syst = idisoSF_->at("NUM_probe_LooseRelTkIso_DEN_HighPtProbes")->evaluate({std::abs(eta),pt,"syst"});

  return std::hypot(stat,syst);
}

double MuonCorrectionHelper::looseIsoSFtracker(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  return idisoSF_->at("NUM_probe_LooseRelTkIso_DEN_TrkHighPtProbes")->evaluate({std::abs(eta),pt,"nominal"});
}

double MuonCorrectionHelper::looseIsoSFtrackerSyst(const pat::MuonRef& mu) {
  double pt = std::min(std::max(mu->tunePMuonBestTrack()->pt(),50.001),1000.-1e-3);
  double eta = mu->tunePMuonBestTrack()->eta();

  double stat = idisoSF_->at("NUM_probe_LooseRelTkIso_DEN_TrkHighPtProbes")->evaluate({std::abs(eta),pt,"stat"});
  double syst = idisoSF_->at("NUM_probe_LooseRelTkIso_DEN_TrkHighPtProbes")->evaluate({std::abs(eta),pt,"syst"});

  return std::hypot(stat,syst);
}


double MuonCorrectionHelper::rochesterData(const pat::MuonRef& mu) const {
  return rochester_.kScaleDT(mu->tunePMuonBestTrack()->charge(),
                             mu->tunePMuonBestTrack()->pt(),
                             mu->tunePMuonBestTrack()->eta(),
                             mu->tunePMuonBestTrack()->phi());
}

double MuonCorrectionHelper::rochesterMC(const pat::MuonRef& mu, const edm::Handle<edm::View<reco::GenParticle>>& ahandle) {
  std::vector<reco::GenParticleRef> matched;

  auto sortByClosestPt = [&mu](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
    return std::abs( mu->tunePMuonBestTrack()->pt() - a->pt() ) < std::abs( mu->tunePMuonBestTrack()->pt() - b->pt() );
  };

  // try prompt muon first
  for (unsigned int idx = 0; idx < ahandle->size(); ++idx) {
    const auto& genptc = ahandle->refAt(idx);

    if ( std::abs(genptc->pdgId())!=13 || genptc->status()!=1 )
      continue;

    if ( genptc->isPromptFinalState() ) {
      const double dr2 = reco::deltaR2(mu->tunePMuonBestTrack()->eta(),
                                       mu->tunePMuonBestTrack()->phi(),
                                       genptc->eta(),
                                       genptc->phi());
      if ( dr2 < 0.01 )
        matched.push_back( genptc.castTo<reco::GenParticleRef>() );
    }
  }

  // then try any final state muon
  if ( matched.empty() ) {
    for (unsigned int idx = 0; idx < ahandle->size(); ++idx) {
      const auto& genptc = ahandle->refAt(idx);

      if ( std::abs(genptc->pdgId())!=13 || genptc->status()!=1 )
        continue;

      const double dr2 = reco::deltaR2(mu->tunePMuonBestTrack()->eta(),
                                       mu->tunePMuonBestTrack()->phi(),
                                       genptc->eta(),
                                       genptc->phi());

      if ( dr2 < 0.01 )
        matched.push_back( genptc.castTo<reco::GenParticleRef>() );
    }
  }

  if ( !matched.empty() ) {
    // GEN matched
    std::sort(matched.begin(),matched.end(),sortByClosestPt);

    return rochester_.kSpreadMC(mu->tunePMuonBestTrack()->charge(),
                                mu->tunePMuonBestTrack()->pt(),
                                mu->tunePMuonBestTrack()->eta(),
                                mu->tunePMuonBestTrack()->phi(),
                                matched.front()->pt());
  }

  if (!mu->isTrackerMuon()) // ???
    return 1.;

  // no GEN matched muon found
  return rochester_.kSmearMC(mu->tunePMuonBestTrack()->charge(),
                             mu->tunePMuonBestTrack()->pt(),
                             mu->tunePMuonBestTrack()->eta(),
                             mu->tunePMuonBestTrack()->phi(),
                             mu->innerTrack()->hitPattern().trackerLayersWithMeasurement(),
                             rng_.Uniform());
}

double MuonCorrectionHelper::nominalData(const pat::MuonRef&mu) const {
  if (mu->tunePMuonBestTrack()->pt() <= 200.)
    return rochesterData(mu);

  return 1.;
}

double MuonCorrectionHelper::nominalMC(const pat::MuonRef& mu, const edm::Handle<edm::View<reco::GenParticle>>& ahandle) {
  if (mu->tunePMuonBestTrack()->pt() <= 200.)
    return rochesterMC(mu,ahandle);

  return 1.;
}

double MuonCorrectionHelper::altScale(const pat::MuonRef&mu, const std::vector<double>& kb, bool isMC) const {
  // if (mu->tunePMuonBestTrack()->pt() <= 200.)
  //   return 1.;

  if (!isMC) // alter MC instead of data
    return 1.;

  const auto tp = mu->tunePMuonBestTrack();
  int col = -1, row = -1;

  if (tp->eta() < -2.1)
    col = 0;
  else if (tp->eta() < -1.2)
    col = 1;
  else if (tp->eta() < 0.)
    col = 2;
  else if (tp->eta() < 1.2)
    col = 3;
  else if (tp->eta() < 2.1)
    col = 4;
  else
    col = 5;

  if (tp->phi() < -M_PI/3.)
    row = 0;
  else if (tp->phi() < M_PI/3.)
    row = 1;
  else
    row = 2;

  unsigned idx = 6*row + col;
  const double bias = kb.at(idx);

  // bias unit is charge/TeV...
  return std::max(1./(1.+static_cast<double>(tp->charge())*tp->pt()*0.001*bias),0.);
}

double MuonCorrectionHelper::smear(const pat::MuonRef& mu,
                                   const std::vector<double>& params,
                                   const std::vector<double>& factors,
                                   bool isMC) {
  if (!isMC)
    return 1.;

  const bool isMB = std::abs(mu->tunePMuonBestTrack()->eta()) < 1.2;
  const double mom = mu->tunePMuonBestTrack()->p();
  const double p0 = isMB ? params.at(0) : params.at(4);
  const double p1 = isMB ? params.at(1) : params.at(5);
  const double p2 = isMB ? params.at(2) : params.at(6);
  const double p3 = isMB ? params.at(3) : params.at(7);
  const double sigma = p0 + p1*mom + p2*mom*mom + p3*mom*mom*mom;
  const double factor = isMB ? factors.at(0) : factors.at(1);

  return (1. + rng_.Gaus(0.,sigma)*factor);
}
