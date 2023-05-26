#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedEleTkIsolFromCands.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

ModifiedEleTkIsolFromCands::TrkCuts::TrkCuts(const edm::ParameterSet& para) {
  auto sq = [](double val) { return val*val; };
  minPt = para.getParameter<double>("minPt");
  minDR2 = sq(para.getParameter<double>("minDR"));
  maxDR2 = sq(para.getParameter<double>("maxDR"));
  minDEta = para.getParameter<double>("minDEta");
  dEta2nd = para.getParameter<double>("dEta2nd");
  dPhi2nd = para.getParameter<double>("dPhi2nd");
  maxDZ = para.getParameter<double>("maxDZ");
  minHits = para.getParameter<int>("minHits");
  minPixelHits = para.getParameter<int>("minPixelHits");
  maxDPtPt = para.getParameter<double>("maxDPtPt");
  addTrkMinPt = para.getParameter<double>("addTrkMinPt");
  addTrkDR2 = sq(para.getParameter<double>("addTrkDR2"));
  addTrkREguard = para.getParameter<double>("addTrkREguard");

  auto qualNames = para.getParameter<std::vector<std::string>>("allowedQualities");
  auto algoNames = para.getParameter<std::vector<std::string>>("algosToReject");

  for (auto& qualName : qualNames)
    allowedQualities.push_back(reco::TrackBase::qualityByName(qualName));

  for(auto& algoName : algoNames)
    algosToReject.push_back(reco::TrackBase::algoByName(algoName));

  std::sort(algosToReject.begin(),algosToReject.end());
}

ModifiedEleTkIsolFromCands::ModifiedEleTkIsolFromCands(const edm::ParameterSet& para):
  barrelCuts_(para.getParameter<edm::ParameterSet>("barrelCuts")),
  endcapCuts_(para.getParameter<edm::ParameterSet>("endcapCuts"))
{}

double ModifiedEleTkIsolFromCands::calIsol(const reco::TrackBase& eleTrk,
                                           const edm::Handle<edm::View<pat::PackedCandidate>>& cands,
                                           const reco::TrackBase& addTrk,
                                           const PIDVeto pidVeto) const {
  return calIsol(eleTrk.eta(),eleTrk.phi(),eleTrk.vz(),cands,addTrk,pidVeto);
}

double ModifiedEleTkIsolFromCands::calIsol(const double eleEta,
                                           const double elePhi,
                                           const double eleVZ,
                                           const edm::Handle<edm::View<pat::PackedCandidate>>& cands,
                                           const reco::TrackBase& addTrk,
                                           const PIDVeto pidVeto) const {
  double ptSum=0.;

  const TrkCuts& cuts = std::abs(eleEta)<1.5 ? barrelCuts_ : endcapCuts_;

  for (unsigned idx = 0; idx < cands->size(); ++idx) {
    const auto& cand = cands->refAt(idx);

    if ( cand->hasTrackDetails() && cand->charge()!=0 && passPIDVeto(cand->pdgId(),pidVeto) ) { // 94X
      const reco::Track& trk = cand->pseudoTrack();

      if ( passTrkSel(trk,trk.pt(),cuts,eleEta,elePhi,eleVZ) ) {
        if( std::abs(addTrk.eta()-trk.eta()) < cuts.dEta2nd && std::abs( reco::deltaPhi(addTrk.phi(),trk.phi()) ) < cuts.dPhi2nd )
          continue;

      	ptSum += trk.pt();
      }
    }
  }

  return ptSum;
}

double ModifiedEleTkIsolFromCands::calIsol(const reco::TrackBase& eleTrk,
                                           const reco::TrackCollection& tracks,
                                           const reco::TrackBase& addTrk) const {
  return calIsol(eleTrk.eta(),eleTrk.phi(),eleTrk.vz(),tracks,addTrk);
}

double ModifiedEleTkIsolFromCands::calIsol(const double eleEta,
                                           const double elePhi,
                                           const double eleVZ,
                                           const reco::TrackCollection& tracks,
                                           const reco::TrackBase& addTrk) const {
  double ptSum=0.;

  const TrkCuts& cuts = std::abs(eleEta)<1.5 ? barrelCuts_ : endcapCuts_;

  for (auto& trk : tracks) {
    if ( passTrkSel(trk,trk.pt(),cuts,eleEta,elePhi,eleVZ) ) {
      if ( std::abs(addTrk.eta()-trk.eta()) < cuts.dEta2nd && std::abs( reco::deltaPhi(addTrk.phi(),trk.phi()) ) < cuts.dPhi2nd )
        continue;

      ptSum += trk.pt();
    }
  }

  return ptSum;
}

bool ModifiedEleTkIsolFromCands::passPIDVeto(const int pdgId, const ModifiedEleTkIsolFromCands::PIDVeto veto) {
  int pidAbs = std::abs(pdgId);

  switch (veto) {
  case PIDVeto::NONE:
    return true;
  case PIDVeto::ELES:
    if (pidAbs==11)
      return false;
    else
      return true;
  case PIDVeto::NONELES:
    if (pidAbs==11)
      return true;
    else
      return false;
  }

  throw cms::Exception("CodeError") <<
    "invalid PIDVeto " << static_cast<int>(veto) << ", " <<
    "this is likely due to some static casting of invalid ints somewhere";
}

ModifiedEleTkIsolFromCands::PIDVeto ModifiedEleTkIsolFromCands::pidVetoFromStr(const std::string& vetoStr) {
  if (vetoStr=="NONE")
    return PIDVeto::NONE;
  else if (vetoStr=="ELES")
    return PIDVeto::ELES;
  else if (vetoStr=="NONELES")
    return PIDVeto::NONELES;
  else {
    throw cms::Exception("CodeError") << "unrecognised string " << vetoStr << ", either a typo or this function needs to be updated";
  }
}

bool ModifiedEleTkIsolFromCands::passTrkSel(const reco::TrackBase& trk,
                                            const double trkPt,
                                            const TrkCuts& cuts,
                                            const double eleEta,
                                            const double elePhi,
                                            const double eleVZ) {
  const float dR2 = reco::deltaR2(eleEta,elePhi,trk.eta(),trk.phi());
  const float dEta = trk.eta()-eleEta;
  const float dZ = eleVZ - trk.vz();

  return dR2 >= cuts.minDR2 && dR2 <= cuts.maxDR2 &&
    std::abs(dEta) >= cuts.minDEta &&
    std::abs(dZ) < cuts.maxDZ &&
    trk.hitPattern().numberOfValidHits() >= cuts.minHits &&
    trk.hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits &&
    (trk.ptError()/trkPt < cuts.maxDPtPt || cuts.maxDPtPt<0) &&
    passQual(trk,cuts.allowedQualities) &&
    passAlgo(trk,cuts.algosToReject) &&
    trkPt > cuts.minPt;
}

bool ModifiedEleTkIsolFromCands::additionalTrkSel(const reco::TrackBase& addTrk,
                                                  const reco::TrackBase& eleTrk,
                                                  const TrkCuts& cuts) {
  const float dR2 = reco::deltaR2(eleTrk.eta(),eleTrk.phi(),addTrk.eta(),addTrk.phi());
  const float dZ = eleTrk.vz() - addTrk.vz();

  return dR2 <= cuts.addTrkDR2 &&
    std::abs(dZ) < cuts.maxDZ &&
    addTrk.hitPattern().numberOfValidHits() >= cuts.minHits &&
    addTrk.hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits &&
    passQual(addTrk,cuts.allowedQualities) &&
    passAlgo(addTrk,cuts.algosToReject) &&
    addTrk.pt() > cuts.addTrkMinPt;
}

bool ModifiedEleTkIsolFromCands::additionalTrkSel(const edm::RefToBase<pat::PackedCandidate>& cand,
                                                  const reco::TrackBase& eleTrk,
                                                  const TrkCuts& cuts) {
  const reco::Track* addTrk = cand->bestTrack();
  const float dR2 = reco::deltaR2(eleTrk.eta(),eleTrk.phi(),addTrk->eta(),addTrk->phi());
  const float dZ = eleTrk.vz() - addTrk->vz();

  return dR2 <= cuts.addTrkDR2 &&
    std::abs(dZ) < cuts.maxDZ &&
    addTrk->hitPattern().numberOfValidHits() >= cuts.minHits &&
    addTrk->hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits &&
    cand->trackHighPurity() &&
    passAlgo(*addTrk,cuts.algosToReject) &&
    addTrk->pt() > cuts.addTrkMinPt;
}

const reco::GsfTrackRef ModifiedEleTkIsolFromCands::additionalGsfTrkSelector(const reco::GsfElectron& ele, const edm::Handle<edm::View<reco::GsfTrack>>& gsfTrks) {
  std::vector<reco::GsfTrackRef> additionalGsfTrks;
  const reco::GsfTrackRef eleTrk = ele.gsfTrack();

  for (unsigned iGsf = 0; iGsf < gsfTrks->size(); iGsf++) {
    const auto& gsfTrk = gsfTrks->refAt(iGsf);
    const reco::GsfTrackRef trkRef = gsfTrk.castTo<reco::GsfTrackRef>();

    // skip the electron's Gsf track (obviously)
    if ( trkRef==eleTrk )
      continue;

    const ModifiedEleTkIsolFromCands::TrkCuts& cuts = std::abs(trkRef->eta()) < 1.5 ? barrelCuts_ : endcapCuts_;

    if ( additionalTrkSel(*trkRef,*eleTrk,cuts) )
      additionalGsfTrks.push_back(trkRef);
  }

  std::sort(additionalGsfTrks.begin(), additionalGsfTrks.end(), [](const reco::GsfTrackRef& a, const reco::GsfTrackRef& b) {return a->pt() > b->pt();} );

  if (additionalGsfTrks.empty())
    return eleTrk;

  return additionalGsfTrks.front();
}

const pat::PackedCandidateRef ModifiedEleTkIsolFromCands::additionalPackedCandSelector(const pat::Electron& ele,
                                                                                       const std::vector<edm::Handle<edm::View<pat::PackedCandidate>>>& cands,
                                                                                       const std::vector<ModifiedEleTkIsolFromCands::PIDVeto>& candVetos) {
  std::vector<pat::PackedCandidateRef> additionalCands;

  for (unsigned iHandle = 0; iHandle < cands.size(); iHandle++) {
    const auto& ahandle = cands.at(iHandle);
    const auto& pidVeto = candVetos.at(iHandle);

    for (unsigned icand = 0; icand < ahandle->size(); icand++) {
      const auto& acand = ahandle->refAt(icand);

      if ( !acand->hasTrackDetails() )
        continue;

      if ( !passPIDVeto(acand->pdgId(),pidVeto) )
        continue;

      if (pidVeto!=PIDVeto::NONE) {
        // skip PF electron candidates & electron CTF tracks since we use GSF tracks instead of them
        // essentially reject all lostTracks:eleTracks products
        if ( std::abs(acand->pdgId())==11 )
          continue;
      }

      const reco::Track* atrack = acand->bestTrack();
      const ModifiedEleTkIsolFromCands::TrkCuts& cuts = std::abs(atrack->eta()) < 1.5 ? barrelCuts_ : endcapCuts_;

      // avoid double-counting of the closest CTF track
      // note the rounding error due to the compression of packedPFCandidates
      if ( ele.closestCtfTrackRef().isNonnull() &&
           std::abs( atrack->eta() - ele.closestCtfTrackRef()->eta() ) < cuts.addTrkREguard &&
           std::abs( atrack->phi() - ele.closestCtfTrackRef()->phi() ) < cuts.addTrkREguard )
        continue;

      if ( additionalTrkSel(acand,*(ele.gsfTrack()),cuts) )
        additionalCands.push_back(acand.castTo<pat::PackedCandidateRef>());
    }
  }

  if (additionalCands.empty())
    return pat::PackedCandidateRef();

  auto sortByTrackPt = [] (const pat::PackedCandidateRef& acand, const pat::PackedCandidateRef& bcand) {
    return acand->bestTrack()->pt() > bcand->bestTrack()->pt();
  };

  std::sort(additionalCands.begin(), additionalCands.end(), sortByTrackPt);

  return additionalCands.front();
}

bool ModifiedEleTkIsolFromCands::passQual(const reco::TrackBase& trk,
                                          const std::vector<reco::TrackBase::TrackQuality>& quals) {
  if (quals.empty()) return true;

  for (auto qual : quals) {
    if (trk.quality(qual))
      return true;
  }

  return false;
}

bool ModifiedEleTkIsolFromCands::passAlgo(const reco::TrackBase& trk,
                                          const std::vector<reco::TrackBase::TrackAlgorithm>& algosToRej) {
  return algosToRej.empty() || !std::binary_search(algosToRej.begin(),algosToRej.end(),trk.algo());
}
