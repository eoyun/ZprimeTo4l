#ifndef MergedLeptonIDs_H
#define MergedLeptonIDs_H 1

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

namespace MergedLeptonIDs {
  enum cutflowElectron {
    baseline = -1,
    minEnergy,
    ecalDriven,
    dPhiIn,
    saturatedXtal,
    missingInnerHits,
    trackIso,
    HoE,
    caloIso,
    dxy, // passed modified HEEP
    dEtaIn,
    showerShape, // passed HEEP
    failedHEEP = 20,
    no2ndGsf,
    has2ndGsf,
    hasMatchedElectron, // has 2nd GSF track but it's from another electron

  };

  inline bool isSameGsfTrack(const reco::GsfTrackRef& aTrk, const reco::GsfTrackRef& bTrk) {
    return ( aTrk->pt()==bTrk->pt() &&
             aTrk->eta()==bTrk->eta() &&
             aTrk->phi()==bTrk->phi() );
  }

  bool isModifiedHEEP(const reco::GsfElectron& el,
                      const reco::Vertex& primaryVertex,
                      const float& trkIso,
                      const float& ecalIso,
                      const int& nrSatCrys,
                      const double& rho,
                      cutflowElectron& cutflow);

  bool hasPassedHEEP(const reco::GsfElectron& el,
                     const reco::Vertex& primaryVertex,
                     const int& nrSatCrys,
                     const double& rho,
                     cutflowElectron& cutflow);

  bool isHighPtTrackerMuon(const reco::Muon& muon,
                           const reco::Vertex& vtx);
} /* MergedLeptonIDs */

#endif
