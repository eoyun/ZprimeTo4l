#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

bool MergedLeptonIDs::isHighPtTrackerMuon(const reco::Muon& muon, const reco::Vertex& vtx) {
  if (!muon.isTrackerMuon())
    return false;

  bool muMatchedSt = muon.numberOfMatchedStations() > 1;
  if (!muMatchedSt) {
    if (muon.isTrackerMuon() && muon.numberOfMatchedStations() == 1) {
      if (muon.expectedNnumberOfMatchedStations() < 2 || !(muon.stationMask() == 1 || muon.stationMask() == 16) ||
          muon.numberOfMatchedRPCLayers() > 2)
        muMatchedSt = true;
    }
  }

  bool muID = muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
              muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError() / muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip =
      std::abs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && std::abs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;
}

bool MergedLeptonIDs::isModifiedHEEP(const reco::GsfElectron& el,
                                     const reco::Vertex& primaryVertex,
                                     const float& trkIso,
                                     const float& ecalIso,
                                     const int& nrSatCrys,
                                     const double& rho,
                                     int& cutflow) {
  auto sq = [](const double& val) { return val*val; };
  double R = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()) + sq(el.superCluster()->z()));
  double Rt = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()));
  double etSC = (el.superCluster()->energy())*(Rt/R);
  double etaSC = el.superCluster()->eta();
  double caloIsoVal = ecalIso + el.dr03HcalDepth1TowerSumEt();
  if ( etSC < 20. )
    return false;
  cutflow++;
  if ( !el.ecalDriven() )
    return false;
  cutflow++;
  if ( !(std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.06) )
    return false;
  cutflow++;
  if ( nrSatCrys > 0 )
    return false;
  cutflow++;
  if ( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1 )
    return false;
  cutflow++;
  if ( trkIso > 5. )
    return false;
  cutflow++;
  if ( std::abs(etaSC) < 1.4442 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (1./el.energy() + 0.05) );
    bool caloIso = ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.02;

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut );
  } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (5./el.energy() + 0.05) );
    bool caloIso = ( etSC < 50. ) ? ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho )
                                  : ( caloIsoVal < 2 + 0.03*(etSC-50.) + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.05;

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut );
  }

  return false;
}

bool MergedLeptonIDs::hasPassedHEEP(const reco::GsfElectron& el,
                                    const reco::Vertex& primaryVertex,
                                    const int& nrSatCrys,
                                    const double& rho,
                                    int& cutflow) {
  auto sq = [](const double& val) { return val*val; };
  double R = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()) + sq(el.superCluster()->z()));
  double Rt = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()));
  double etSC = (el.superCluster()->energy())*(Rt/R);
  double etaSC = el.superCluster()->eta();
  double caloIsoVal = el.dr03EcalRecHitSumEt() + el.dr03HcalDepth1TowerSumEt();
  if ( etSC < 20. )
    return false;
  cutflow++;
  if ( !el.ecalDriven() )
    return false;
  cutflow++;
  if ( !(std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.06) )
    return false;
  cutflow++;
  if ( nrSatCrys > 0 )
    return false;
  cutflow++;
  if ( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1 )
    return false;
  cutflow++;
  if ( el.dr03TkSumPtHEEP() > 5. )
    return false;
  cutflow++;
  if ( std::abs(etaSC) < 1.4442 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (1./el.energy() + 0.05) );
    bool caloIso = ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.02;
    bool dEtaInCut = std::abs( el.deltaEtaSeedClusterTrackAtVtx() ) < 0.004;
    bool ssCut = ( el.full5x5_e2x5Max()/el.full5x5_e5x5() > 0.94 ||
                   el.full5x5_e1x5()/el.full5x5_e5x5() > 0.83 );

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;
    if (dEtaInCut)
      cutflow++;
    if (ssCut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut && dEtaInCut && ssCut );
  } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (5./el.energy() + 0.05) );
    bool caloIso = ( etSC < 50. ) ? ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho )
                                  : ( caloIsoVal < 2 + 0.03*(etSC-50.) + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.05;
    bool dEtaInCut = std::abs( el.deltaEtaSeedClusterTrackAtVtx() ) < 0.006;
    bool ssCut = ( el.full5x5_sigmaIetaIeta() < 0.03 );

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;
    if (dEtaInCut)
      cutflow++;
    if (ssCut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut && dEtaInCut && ssCut );
  }

  return false;
}