#ifndef ZprimeTo4l_MergedLepton_MergedMuonTkIsoFromCands_H
#define ZprimeTo4l_MergedLepton_MergedMuonTkIsoFromCands_H 1

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

class MergedMuonTkIsolFromCands {
public:
  enum class PIDVeto {
    NONE = 0,
    MUNS,
    NONMUNS,
  };

  struct TrkCuts {
    float minPt;
    float minDR2;
    float maxDR2;
    float minDEta;
    float dEta2nd;
    float dPhi2nd;
    float maxDZ;
    float minHits;
    float minPixelHits;
    float maxDPtPt;
    float addTrkMinPt;
    float addTrkDR2;
    float addTrkREguard;
    float addTrkHoE;
    std::vector<reco::TrackBase::TrackQuality> allowedQualities;
    std::vector<reco::TrackBase::TrackAlgorithm> algosToReject;
    explicit TrkCuts(const edm::ParameterSet& para);
  };

  TrkCuts cuts_;

  explicit MergedMuonTkIsolFromCands(const edm::ParameterSet& para, edm::ConsumesCollector iC);
  MergedMuonTkIsolFromCands(const MergedMuonTkIsolFromCands&) = default;
  ~MergedMuonTkIsolFromCands() = default;
  MergedMuonTkIsolFromCands& operator=(const MergedMuonTkIsolFromCands&) = default;

  double calIsol(const reco::TrackBase& trk,
                 const edm::Handle<edm::View<pat::PackedCandidate>>& cands,
                 const reco::TrackBase& addTrk,
                 const PIDVeto = PIDVeto::NONE) const;
  double calIsol(const double muEta,
                 const double muPhi,
                 const double muVZ,
                 const edm::Handle<edm::View<pat::PackedCandidate>>& cands,
                 const reco::TrackBase& addTrk,
                 const PIDVeto = PIDVeto::NONE) const;

  double calIsol(const reco::TrackBase& trk, const reco::TrackCollection& tracks, const reco::TrackBase& addTrk) const;
  double calIsol(const double muEta,
                 const double muPhi,
                 const double muVZ,
                 const reco::TrackCollection& tracks,
                 const reco::TrackBase& addTrk) const;

  template <typename... Args>
  double calIsolPt(Args&&... args) const {
    return calIsol(std::forward<Args>(args)...);
  }

  static PIDVeto pidVetoFromStr(const std::string& vetoStr);
  static bool passPIDVeto(const int pdgId, const MergedMuonTkIsolFromCands::PIDVeto pidVeto);

  bool additionalTrkSel(const reco::TrackBase& addTrk, const reco::TrackBase& muTrk, const TrkCuts& cuts);
  bool additionalTrkSel(const edm::RefToBase<pat::PackedCandidate>& cand, const reco::TrackBase& muTrk, const TrkCuts& cuts);

  const pat::PackedCandidateRef additionalPackedCandSelector(
      const pat::Muon& mu,
      const std::vector<edm::Handle<edm::View<pat::PackedCandidate>>>& cands,
      const std::vector<MergedMuonTkIsolFromCands::PIDVeto>& candVetos,
      const edm::EventSetup& iSetup);

private:
  edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;

  static bool passTrkSel(const reco::TrackBase& trk,
                         const double trkPt,
                         const TrkCuts& cuts,
                         const double muEta,
                         const double muPhi,
                         const double muVZ);
  static bool passQual(const reco::TrackBase& trk, const std::vector<reco::TrackBase::TrackQuality>& quals);
  static bool passAlgo(const reco::TrackBase& trk, const std::vector<reco::TrackBase::TrackAlgorithm>& algosToRej);
};

#endif
