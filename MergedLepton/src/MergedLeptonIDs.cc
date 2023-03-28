#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"

MergedLeptonIDs::GSFtype MergedLeptonIDs::checkElectronGSFtype( const pat::ElectronRef& aEle,
                                                                const reco::GsfTrackRef& orgGsfTrk,
                                                                const reco::GsfTrackRef& addGsfTrk,
                                                                const double etThresEB,
                                                                const double etThresEE,
                                                                const double minEt ) {
  double dr2 = reco::deltaR2(orgGsfTrk->eta(),orgGsfTrk->phi(),addGsfTrk->eta(),addGsfTrk->phi());
  auto square = [](double input) { return input*input; };

  if ( aEle->et() < minEt )
    return GSFtype::nulltype;

  if ( aEle->isEB() ) {
    if ( dr2 < square(0.0174) && aEle->et() > etThresEB )
      return GSFtype::DR1Et2EB;
    else if ( dr2 < square(0.0174) )
      return GSFtype::extendedCR;
    else if ( dr2 > square(0.0174) ) {
      if ( aEle->et() > etThresEB )
        return GSFtype::DR2Et2EB;
      else
        return GSFtype::DR2Et1EB;
    }
  } else {
    if ( dr2 < square( 0.00864*std::abs(std::sinh(aEle->eta())) ) && aEle->et() > etThresEE )
      return GSFtype::DR1Et2EE;
    else if ( dr2 > square( 0.00864*std::abs(std::sinh(aEle->eta())) ) ) {
      if ( aEle->et() > etThresEE )
        return GSFtype::DR2Et2EE;
      else
        return GSFtype::DR2Et1EE;
    }
  }

  return GSFtype::nulltype;
}

void MergedLeptonIDs::fillCutflow(TH1* ahist, const int acutflow, const double aWeight) {
  float acutfloat = static_cast<float>(acutflow) + 0.5;
  int bin = ahist->FindBin(acutfloat);

  for (int idx = 1; idx <= bin; idx++)
    ahist->Fill( ahist->GetBinCenter(idx), aWeight );
}

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
                                     const double& rho,
                                     cutflowElectron& cutflow) {
  cutflow = cutflowElectron::baseline;
  double etaSC = el.superCluster()->eta();
  double caloIsoVal = ecalIso + el.dr03HcalDepth1TowerSumEt();
  if ( el.et() < 20. )
    return false;
  cutflow = cutflowElectron::minEnergy;
  if ( !el.ecalDriven() )
    return false;
  cutflow = cutflowElectron::ecalDriven;
  if ( !(std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.06) )
    return false;
  cutflow = cutflowElectron::dPhiIn;
  if ( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1 )
    return false;
  cutflow = cutflowElectron::missingInnerHits;
  if ( trkIso > 5. )
    return false;
  cutflow = cutflowElectron::trackIso;
  if ( std::abs(etaSC) < 1.4442 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (1./el.energy() + 0.05) );
    bool caloIso = ( caloIsoVal < 2 + 0.03*el.et() + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.02;

    if (HoEcut)
      cutflow = cutflowElectron::HoE;
    if (caloIso)
      cutflow = cutflowElectron::caloIso;
    if (dxycut)
      cutflow = cutflowElectron::dxy;

    return ( caloIso && HoEcut && dxycut );
  } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (5./el.energy() + 0.05) );
    bool caloIso = ( el.et() < 50. ) ? ( caloIsoVal < 2 + 0.03*el.et() + 0.28*rho )
                                  : ( caloIsoVal < 2 + 0.03*(el.et()-50.) + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.05;

    if (HoEcut)
      cutflow = cutflowElectron::HoE;
    if (caloIso)
      cutflow = cutflowElectron::caloIso;
    if (dxycut)
      cutflow = cutflowElectron::dxy;

    return ( caloIso && HoEcut && dxycut );
  }

  return false;
}

bool MergedLeptonIDs::hasPassedHEEP(const reco::GsfElectron& el,
                                    const reco::Vertex& primaryVertex,
                                    const double& rho,
                                    cutflowElectron& cutflow) {
  cutflow = cutflowElectron::baseline;
  double etaSC = el.superCluster()->eta();
  double caloIsoVal = el.dr03EcalRecHitSumEt() + el.dr03HcalDepth1TowerSumEt();
  if ( el.et() < 20. )
    return false;
  cutflow = cutflowElectron::minEnergy;
  if ( !el.ecalDriven() )
    return false;
  cutflow = cutflowElectron::ecalDriven;
  if ( !(std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.06) )
    return false;
  cutflow = cutflowElectron::dPhiIn;
  if ( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1 )
    return false;
  cutflow = cutflowElectron::missingInnerHits;
  if ( el.dr03TkSumPtHEEP() > 5. )
    return false;
  cutflow = cutflowElectron::trackIso;
  if ( std::abs(etaSC) < 1.4442 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (1./el.energy() + 0.05) );
    bool caloIso = ( caloIsoVal < 2 + 0.03*el.et() + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.02;
    bool dEtaInCut = std::abs( el.deltaEtaSeedClusterTrackAtVtx() ) < 0.004;
    bool ssCut = ( el.full5x5_e2x5Max()/el.full5x5_e5x5() > 0.94 ||
                   el.full5x5_e1x5()/el.full5x5_e5x5() > 0.83 );

    if (HoEcut)
      cutflow = cutflowElectron::HoE;
    if (caloIso)
      cutflow = cutflowElectron::caloIso;
    if (dxycut)
      cutflow = cutflowElectron::dxy;
    if (dEtaInCut)
      cutflow = cutflowElectron::dEtaIn;
    if (ssCut)
      cutflow = cutflowElectron::showerShape;

    return ( caloIso && HoEcut && dxycut && dEtaInCut && ssCut );
  } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (5./el.energy() + 0.05) );
    bool caloIso = ( el.et() < 50. ) ? ( caloIsoVal < 2 + 0.03*el.et() + 0.28*rho )
                                  : ( caloIsoVal < 2 + 0.03*(el.et()-50.) + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.05;
    bool dEtaInCut = std::abs( el.deltaEtaSeedClusterTrackAtVtx() ) < 0.006;
    bool ssCut = ( el.full5x5_sigmaIetaIeta() < 0.03 );

    if (HoEcut)
      cutflow = cutflowElectron::HoE;
    if (caloIso)
      cutflow = cutflowElectron::caloIso;
    if (dxycut)
      cutflow = cutflowElectron::dxy;
    if (dEtaInCut)
      cutflow = cutflowElectron::dEtaIn;
    if (ssCut)
      cutflow = cutflowElectron::showerShape;

    return ( caloIso && HoEcut && dxycut && dEtaInCut && ssCut );
  }

  return false;
}
