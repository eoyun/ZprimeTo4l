#include "ZprimeTo4l/MergedLepton/interface/MergedMuonTkIsolFromCands.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "TMath.h"

MergedMuonTkIsolFromCands::TrkCuts::TrkCuts(const edm::ParameterSet& para) {
  auto sq = [](double val) { return val * val; };
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
  addTrkHoE = para.getParameter<double>("addTrkHoE");

  auto qualNames = para.getParameter<std::vector<std::string>>("allowedQualities");
  auto algoNames = para.getParameter<std::vector<std::string>>("algosToReject");

  for (auto& qualName : qualNames)
    allowedQualities.push_back(reco::TrackBase::qualityByName(qualName));

  for (auto& algoName : algoNames)
    algosToReject.push_back(reco::TrackBase::algoByName(algoName));

  std::sort(algosToReject.begin(), algosToReject.end());
}

MergedMuonTkIsolFromCands::MergedMuonTkIsolFromCands(const edm::ParameterSet& para)
    : barrelCuts_(para.getParameter<edm::ParameterSet>("barrelCuts")),
      endcapCuts_(para.getParameter<edm::ParameterSet>("endcapCuts")) {}

MergedMuonTkIsolFromCands::MergedMuonTkIsolFromCands(const edm::ParameterSet& para, edm::ConsumesCollector iC)
    : MergedMuonTkIsolFromCands(para) {
  (void)iC;
}

double MergedMuonTkIsolFromCands::calIsol(const reco::TrackBase& muTrk,
                                          const edm::Handle<edm::View<pat::PackedCandidate>>& cands,
                                          const reco::TrackBase& addTrk,
                                          const PIDVeto pidVeto) const {
  return calIsol(muTrk.eta(), muTrk.phi(), muTrk.vz(), cands, addTrk, pidVeto);
}

double MergedMuonTkIsolFromCands::calIsol(const double muEta,
                                          const double muPhi,
                                          const double muVZ,
                                          const edm::Handle<edm::View<pat::PackedCandidate>>& cands,
                                          const reco::TrackBase& addTrk,
                                          const PIDVeto pidVeto) const {
  double ptSum = 0.;

  const TrkCuts& cuts = std::abs(muEta) < 1.5 ? barrelCuts_ : endcapCuts_;

  for (unsigned idx = 0; idx < cands->size(); ++idx) {
    const auto& cand = cands->refAt(idx);

    if (cand->hasTrackDetails() && cand->charge() != 0 && passPIDVeto(cand->pdgId(), pidVeto)) {
      const reco::Track& trk = cand->pseudoTrack();

      if (passTrkSel(trk, trk.pt(), cuts, muEta, muPhi, muVZ)) {
        if (std::abs(addTrk.eta() - trk.eta()) < cuts.dEta2nd &&
            std::abs(reco::deltaPhi(addTrk.phi(), trk.phi())) < cuts.dPhi2nd)
          continue;

        ptSum += trk.pt();
      }
    }
  }

  return ptSum;
}

double MergedMuonTkIsolFromCands::calIsol(const reco::TrackBase& muTrk,
                                          const reco::TrackCollection& tracks,
                                          const reco::TrackBase& addTrk) const {
  return calIsol(muTrk.eta(), muTrk.phi(), muTrk.vz(), tracks, addTrk);
}

double MergedMuonTkIsolFromCands::calIsol(const double muEta,
                                          const double muPhi,
                                          const double muVZ,
                                          const reco::TrackCollection& tracks,
                                          const reco::TrackBase& addTrk) const {
  double ptSum = 0.;

  const TrkCuts& cuts = std::abs(muEta) < 1.5 ? barrelCuts_ : endcapCuts_;

  for (auto& trk : tracks) {
    if (passTrkSel(trk, trk.pt(), cuts, muEta, muPhi, muVZ)) {
      if (std::abs(addTrk.eta() - trk.eta()) < cuts.dEta2nd &&
          std::abs(reco::deltaPhi(addTrk.phi(), trk.phi())) < cuts.dPhi2nd)
        continue;

      ptSum += trk.pt();
    }
  }

  return ptSum;
}

bool MergedMuonTkIsolFromCands::passPIDVeto(const int pdgId, const MergedMuonTkIsolFromCands::PIDVeto veto) {
  int pidAbs = std::abs(pdgId);

  switch (veto) {
    case PIDVeto::NONE:
      return true;
    case PIDVeto::MUNS:
      return pidAbs != 13;
    case PIDVeto::NONMUNS:
      return pidAbs == 13;
  }

  throw cms::Exception("CodeError") << "invalid PIDVeto " << static_cast<int>(veto)
                                    << ", this is likely due to static casting of invalid ints somewhere";
}

MergedMuonTkIsolFromCands::PIDVeto MergedMuonTkIsolFromCands::pidVetoFromStr(const std::string& vetoStr) {
  if (vetoStr == "NONE")
    return PIDVeto::NONE;
  if (vetoStr == "MUNS")
    return PIDVeto::MUNS;
  if (vetoStr == "NONMUNS")
    return PIDVeto::NONMUNS;

  throw cms::Exception("CodeError") << "unrecognised string " << vetoStr;
}

bool MergedMuonTkIsolFromCands::passTrkSel(const reco::TrackBase& trk,
                                           const double trkPt,
                                           const TrkCuts& cuts,
                                           const double muEta,
                                           const double muPhi,
                                           const double muVZ) {
  const float dR2 = reco::deltaR2(muEta, muPhi, trk.eta(), trk.phi());
  const float dEta = trk.eta() - muEta;
  const float dZ = muVZ - trk.vz();

  return dR2 >= cuts.minDR2 && dR2 <= cuts.maxDR2 && std::abs(dEta) >= cuts.minDEta && std::abs(dZ) < cuts.maxDZ &&
         trk.hitPattern().numberOfValidHits() >= cuts.minHits &&
         trk.hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits &&
         (trk.ptError() / trkPt < cuts.maxDPtPt || cuts.maxDPtPt < 0) && passQual(trk, cuts.allowedQualities) &&
         passAlgo(trk, cuts.algosToReject) && trkPt > cuts.minPt;
}

bool MergedMuonTkIsolFromCands::additionalTrkSel(const reco::TrackBase& addTrk,
                                                 const reco::TrackBase& muTrk,
                                                 const TrkCuts& cuts) {
  const float dR2 = reco::deltaR2(muTrk.eta(), muTrk.phi(), addTrk.eta(), addTrk.phi());
  const float dZ = muTrk.vz() - addTrk.vz();

  return dR2 <= cuts.addTrkDR2 && std::abs(dZ) < cuts.maxDZ &&
         addTrk.hitPattern().numberOfValidHits() >= cuts.minHits &&
         addTrk.hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits && passQual(addTrk, cuts.allowedQualities) &&
         passAlgo(addTrk, cuts.algosToReject) && addTrk.pt() > cuts.addTrkMinPt;
}

bool MergedMuonTkIsolFromCands::additionalTrkSel(const edm::RefToBase<pat::PackedCandidate>& cand,
                                                 const reco::TrackBase& muTrk,
                                                 const TrkCuts& cuts) {
  const reco::TrackRef bestTrkRef = cand->bestTrackRef();
  const reco::Track addTrk = (bestTrkRef.isNonnull() && bestTrkRef.isAvailable()) ? *bestTrkRef : cand->pseudoTrack();
  const float dR2 = reco::deltaR2(muTrk.eta(), muTrk.phi(), addTrk.eta(), addTrk.phi());
  const float dZ = muTrk.vz() - addTrk.vz();

  return dR2 <= cuts.addTrkDR2 && std::abs(dZ) < cuts.maxDZ &&
         addTrk.hitPattern().numberOfValidHits() >= cuts.minHits &&
         addTrk.hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits && cand->trackHighPurity() &&
         cand->hcalFraction() < cuts.addTrkHoE / (1. + cuts.addTrkHoE) && passAlgo(addTrk, cuts.algosToReject) &&
         addTrk.pt() > cuts.addTrkMinPt;
}

const pat::PackedCandidateRef MergedMuonTkIsolFromCands::additionalPackedCandSelector(
    const pat::Muon& mu,
    const std::vector<edm::Handle<edm::View<pat::PackedCandidate>>>& cands,
    const std::vector<MergedMuonTkIsolFromCands::PIDVeto>& candVetos,
    const TransientTrackBuilder& TTbuilder) {
  std::vector<std::pair<pat::PackedCandidateRef, double>> additionalCands;

  const reco::TrackRef muTrkRef = mu.muonBestTrack();
  if (muTrkRef.isNull() || !muTrkRef.isAvailable())
    return pat::PackedCandidateRef();

  auto fitter = KalmanVertexFitter();
  auto firstMu = TTbuilder.build(muTrkRef);

  for (unsigned iHandle = 0; iHandle < cands.size(); iHandle++) {
    const auto& ahandle = cands.at(iHandle);
    const auto& pidVeto = candVetos.at(iHandle);

    for (unsigned icand = 0; icand < ahandle->size(); icand++) {
      const auto& acand = ahandle->refAt(icand);

      if (!acand->hasTrackDetails())
        continue;
      if (!passPIDVeto(acand->pdgId(), pidVeto))
        continue;

      const reco::TrackRef bestTrkRef = acand->bestTrackRef();
      const reco::Track atrack = (bestTrkRef.isNonnull() && bestTrkRef.isAvailable()) ? *bestTrkRef : acand->pseudoTrack();
      const TrkCuts& cuts = std::abs(atrack.eta()) < 1.5 ? barrelCuts_ : endcapCuts_;

      if (reco::deltaR2(atrack.eta(), atrack.phi(), muTrkRef->eta(), muTrkRef->phi()) <
          cuts.addTrkREguard * cuts.addTrkREguard)
        continue;

      if (additionalTrkSel(acand, *muTrkRef, cuts)) {
        if (std::isnan(atrack.dzError()) || std::isinf(atrack.dzError()) || std::isnan(atrack.dxyError()) ||
            std::isinf(atrack.dxyError()) || std::isnan(atrack.d0Error()) || std::isinf(atrack.d0Error()))
          continue;

        std::vector<reco::TransientTrack> trackPair;
        trackPair.push_back(firstMu);
        trackPair.push_back(TTbuilder.build(atrack));

        auto aVtx = fitter.vertex(trackPair);
        if (aVtx.isValid()) {
          double prob = TMath::Prob(aVtx.totalChiSquared(), static_cast<int>(std::rint(aVtx.degreesOfFreedom())));
          additionalCands.push_back(std::make_pair(acand.castTo<pat::PackedCandidateRef>(), prob));
        }
      }
    }
  }

  if (additionalCands.empty())
    return pat::PackedCandidateRef();

  auto sortByProb = [](const std::pair<pat::PackedCandidateRef, double>& a,
                       const std::pair<pat::PackedCandidateRef, double>& b) { return a.second > b.second; };
  std::sort(additionalCands.begin(), additionalCands.end(), sortByProb);

  return additionalCands.front().first;
}

bool MergedMuonTkIsolFromCands::passQual(const reco::TrackBase& trk,
                                         const std::vector<reco::TrackBase::TrackQuality>& quals) {
  if (quals.empty())
    return true;

  for (auto qual : quals) {
    if (trk.quality(qual))
      return true;
  }

  return false;
}

bool MergedMuonTkIsolFromCands::passAlgo(const reco::TrackBase& trk,
                                         const std::vector<reco::TrackBase::TrackAlgorithm>& algosToRej) {
  return algosToRej.empty() || !std::binary_search(algosToRej.begin(), algosToRej.end(), trk.algo());
}
