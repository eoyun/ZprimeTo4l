#include "ZprimeTo4l/Analysis/interface/PairingHelper.h"

bool PairingHelper::pairBoostedElectrons(const pat::ElectronRef& a,
                                         const pat::ElectronRef& b,
                                         const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkMap) {
  return ( a->gsfTrack()==(*addGsfTrkMap)[b] && b->gsfTrack()==(*addGsfTrkMap)[a] );
}

bool PairingHelper::scanBoostedElectrons(const std::vector<pat::ElectronRef>& eles,
                                         const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkMap,
                                         std::pair<pat::ElectronRef,pat::ElectronRef>& pair) {
  bool paired = false;

  for (unsigned idx = 0; idx < eles.size()-1; idx++) {
    const auto& aEle1 = eles.at(idx);

    if ( aEle1->gsfTrack()==(*addGsfTrkMap)[aEle1] ) // there is nearby GSF track (2)
      continue;

    for (unsigned jdx = idx+1; jdx < eles.size(); jdx++) {
      const auto& aEle2 = eles.at(jdx);

      if ( aEle2->gsfTrack()==(*addGsfTrkMap)[aEle2] ) // there is nearby GSF track (2)
        continue;

      if ( pairBoostedElectrons(aEle1, aEle2, addGsfTrkMap) ) {
        pair = std::make_pair(aEle1,aEle2);
        paired = true;
        break;
      }
    }

    if (paired)
      break;
  }

  return paired;
}

bool PairingHelper::pair4E(const std::vector<pat::ElectronRef>& eles,
                           const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkMap,
                           std::pair<pat::ElectronRef,pat::ElectronRef>& pair1,
                           std::pair<pat::ElectronRef,pat::ElectronRef>& pair2) {
  if ( eles.size()!=4 )
    return false; // supposed to work only with 4E

  bool boosted = scanBoostedElectrons(eles,addGsfTrkMap,pair1);

  // 1st & 2nd boosted
  if (boosted) {
    std::vector<pat::ElectronRef> leftover;
    for (unsigned idx = 0; idx < eles.size(); idx++) {
      if ( eles.at(idx)!=pair1.first && eles.at(idx)!=pair1.second )
        leftover.push_back(eles.at(idx));
    }

    pat::ElectronRef aEle3 = leftover.front();
    pat::ElectronRef aEle4 = leftover.at(1);

    if ( aEle3->gsfTrack()==(*addGsfTrkMap)[aEle3] && aEle4->gsfTrack()==(*addGsfTrkMap)[aEle4] ) {
      // 3rd & 4th not boosted
      pair2 = std::make_pair(aEle3,aEle4);
      return true;
    } else {
      // 3rd & 4th boosted
      if ( pairBoostedElectrons(aEle3,aEle4,addGsfTrkMap) ) {
        pair2 = std::make_pair(aEle3,aEle4);
        return true;
      }

      return false;
    }
  } // 1st & 2nd boosted

  // now consider the case that nothing boosted
  if ( eles.at(0)->gsfTrack()!=(*addGsfTrkMap)[eles.at(0)] ||
       eles.at(1)->gsfTrack()!=(*addGsfTrkMap)[eles.at(1)] ||
       eles.at(2)->gsfTrack()!=(*addGsfTrkMap)[eles.at(2)] ||
       eles.at(3)->gsfTrack()!=(*addGsfTrkMap)[eles.at(3)] )
    return false;

  return pairByDR(eles,pair1,pair2);
}

bool PairingHelper::pairByInvM(const std::vector<pat::ElectronRef>& eles,
                               std::pair<pat::ElectronRef,pat::ElectronRef>& pair1,
                               std::pair<pat::ElectronRef,pat::ElectronRef>& pair2) {
  // what a ugly function isn't it
  std::vector< std::tuple< std::pair<pat::ElectronRef,pat::ElectronRef>,
                           std::pair<pat::ElectronRef,pat::ElectronRef> > > casecontainer;

  casecontainer.push_back( std::make_tuple( std::make_pair(eles.at(0),eles.at(1)), std::make_pair(eles.at(2),eles.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(eles.at(0),eles.at(2)), std::make_pair(eles.at(1),eles.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(eles.at(0),eles.at(3)), std::make_pair(eles.at(1),eles.at(2)) ) );

  // i love lambda
  auto sortByMassDiff = [] (const std::tuple<std::pair<pat::ElectronRef,pat::ElectronRef>,std::pair<pat::ElectronRef,pat::ElectronRef>>& a,
                            const std::tuple<std::pair<pat::ElectronRef,pat::ElectronRef>,std::pair<pat::ElectronRef,pat::ElectronRef>>& b) {
    const auto lvecE1a = std::get<0>(a).first->polarP4()*std::get<0>(a).first->userFloat("ecalTrkEnergyPostCorr")/std::get<0>(a).first->energy();
    const auto lvecE2a = std::get<0>(a).second->polarP4()*std::get<0>(a).second->userFloat("ecalTrkEnergyPostCorr")/std::get<0>(a).second->energy();
    const auto lvecE3a = std::get<1>(a).first->polarP4()*std::get<1>(a).first->userFloat("ecalTrkEnergyPostCorr")/std::get<1>(a).first->energy();
    const auto lvecE4a = std::get<1>(a).second->polarP4()*std::get<1>(a).second->userFloat("ecalTrkEnergyPostCorr")/std::get<1>(a).second->energy();
    const auto lvecA1a = lvecE1a + lvecE2a;
    const auto lvecA2a = lvecE3a + lvecE4a;

    const auto lvecE1b = std::get<0>(b).first->polarP4()*std::get<0>(b).first->userFloat("ecalTrkEnergyPostCorr")/std::get<0>(b).first->energy();
    const auto lvecE2b = std::get<0>(b).second->polarP4()*std::get<0>(b).second->userFloat("ecalTrkEnergyPostCorr")/std::get<0>(b).second->energy();
    const auto lvecE3b = std::get<1>(b).first->polarP4()*std::get<1>(b).first->userFloat("ecalTrkEnergyPostCorr")/std::get<1>(b).first->energy();
    const auto lvecE4b = std::get<1>(b).second->polarP4()*std::get<1>(b).second->userFloat("ecalTrkEnergyPostCorr")/std::get<1>(b).second->energy();
    const auto lvecA1b = lvecE1b + lvecE2b;
    const auto lvecA2b = lvecE3b + lvecE4b;

    return std::abs( lvecA1a.M()-lvecA2a.M() ) < std::abs( lvecA1b.M()-lvecA2b.M() );
  };

  // choose pairing with the smallest |M(A1)-M(A2)|
  std::sort(casecontainer.begin(),casecontainer.end(),sortByMassDiff);

  pair1 = std::get<0>(casecontainer.front());
  pair2 = std::get<1>(casecontainer.front());

  return true;
}

bool PairingHelper::pairByDR(const std::vector<pat::ElectronRef>& eles,
                             std::pair<pat::ElectronRef,pat::ElectronRef>& pair1,
                             std::pair<pat::ElectronRef,pat::ElectronRef>& pair2) {
  std::vector< std::tuple< std::pair<pat::ElectronRef,pat::ElectronRef>,
                           std::pair<pat::ElectronRef,pat::ElectronRef> > > casecontainer;

  casecontainer.push_back( std::make_tuple( std::make_pair(eles.at(0),eles.at(1)), std::make_pair(eles.at(2),eles.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(eles.at(0),eles.at(2)), std::make_pair(eles.at(1),eles.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(eles.at(0),eles.at(3)), std::make_pair(eles.at(1),eles.at(2)) ) );

  auto sortByDRDiff = [] (const std::tuple<std::pair<pat::ElectronRef,pat::ElectronRef>,std::pair<pat::ElectronRef,pat::ElectronRef>>& a,
                          const std::tuple<std::pair<pat::ElectronRef,pat::ElectronRef>,std::pair<pat::ElectronRef,pat::ElectronRef>>& b) {
    const double dR1a = reco::deltaR2(std::get<0>(a).first->eta(),std::get<0>(a).first->phi(),std::get<0>(a).second->eta(),std::get<0>(a).second->phi());
    const double dR2a = reco::deltaR2(std::get<1>(a).first->eta(),std::get<1>(a).first->phi(),std::get<1>(a).second->eta(),std::get<1>(a).second->phi());
    const double dR1b = reco::deltaR2(std::get<0>(b).first->eta(),std::get<0>(b).first->phi(),std::get<0>(b).second->eta(),std::get<0>(b).second->phi());
    const double dR2b = reco::deltaR2(std::get<1>(b).first->eta(),std::get<1>(b).first->phi(),std::get<1>(b).second->eta(),std::get<1>(b).second->phi());

    return std::min( dR1a, dR2a ) < std::min( dR1b, dR2b );
  };

  std::sort(casecontainer.begin(),casecontainer.end(),sortByDRDiff);

  pair1 = std::get<0>(casecontainer.front());
  pair2 = std::get<1>(casecontainer.front());

  return true;
}

bool PairingHelper::pair4M(const std::vector<pat::MuonRef>& muons,
                           std::pair<pat::MuonRef,pat::MuonRef>& pair1,
                           std::pair<pat::MuonRef,pat::MuonRef>& pair2) {
  if ( muons.size()!=4 )
    return false; // supposed to work only with 4M

  return pairByDR(muons,pair1,pair2);
}

bool PairingHelper::pairByDR(const std::vector<pat::MuonRef>& muons,
                             std::pair<pat::MuonRef,pat::MuonRef>& pair1,
                             std::pair<pat::MuonRef,pat::MuonRef>& pair2) {
  std::vector< std::tuple< std::pair<pat::MuonRef,pat::MuonRef>,
                           std::pair<pat::MuonRef,pat::MuonRef> > > casecontainer;

  casecontainer.push_back( std::make_tuple( std::make_pair(muons.at(0),muons.at(1)), std::make_pair(muons.at(2),muons.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(muons.at(0),muons.at(2)), std::make_pair(muons.at(1),muons.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(muons.at(0),muons.at(3)), std::make_pair(muons.at(1),muons.at(2)) ) );

  auto sortByDRDiff = [] (const std::tuple<std::pair<pat::MuonRef,pat::MuonRef>,std::pair<pat::MuonRef,pat::MuonRef>>& a,
                          const std::tuple<std::pair<pat::MuonRef,pat::MuonRef>,std::pair<pat::MuonRef,pat::MuonRef>>& b) {
    const double dR1a = reco::deltaR2(std::get<0>(a).first->eta(),std::get<0>(a).first->phi(),std::get<0>(a).second->eta(),std::get<0>(a).second->phi());
    const double dR2a = reco::deltaR2(std::get<1>(a).first->eta(),std::get<1>(a).first->phi(),std::get<1>(a).second->eta(),std::get<1>(a).second->phi());
    const double dR1b = reco::deltaR2(std::get<0>(b).first->eta(),std::get<0>(b).first->phi(),std::get<0>(b).second->eta(),std::get<0>(b).second->phi());
    const double dR2b = reco::deltaR2(std::get<1>(b).first->eta(),std::get<1>(b).first->phi(),std::get<1>(b).second->eta(),std::get<1>(b).second->phi());

    return std::min( dR1a, dR2a ) < std::min( dR1b, dR2b );
  };

  std::sort(casecontainer.begin(),casecontainer.end(),sortByDRDiff);

  pair1 = std::get<0>(casecontainer.front());
  pair2 = std::get<1>(casecontainer.front());

  return true;
}

bool PairingHelper::pairByDR(const std::vector<edm::Ptr<reco::RecoCandidate>>& leptons,
                             std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>& pair1,
                             std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>& pair2) {
  std::vector< std::tuple< std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>,
                           std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>> > > casecontainer;

  casecontainer.push_back( std::make_tuple( std::make_pair(leptons.at(0),leptons.at(1)), std::make_pair(leptons.at(2),leptons.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(leptons.at(0),leptons.at(2)), std::make_pair(leptons.at(1),leptons.at(3)) ) );
  casecontainer.push_back( std::make_tuple( std::make_pair(leptons.at(0),leptons.at(3)), std::make_pair(leptons.at(1),leptons.at(2)) ) );

  auto sortByDRDiff = [] (const std::tuple<std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>,
                                           std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>>& a,
                          const std::tuple<std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>,
                                           std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>>>& b) {
    const double dR1a = reco::deltaR2(std::get<0>(a).first->eta(),std::get<0>(a).first->phi(),std::get<0>(a).second->eta(),std::get<0>(a).second->phi());
    const double dR2a = reco::deltaR2(std::get<1>(a).first->eta(),std::get<1>(a).first->phi(),std::get<1>(a).second->eta(),std::get<1>(a).second->phi());
    const double dR1b = reco::deltaR2(std::get<0>(b).first->eta(),std::get<0>(b).first->phi(),std::get<0>(b).second->eta(),std::get<0>(b).second->phi());
    const double dR2b = reco::deltaR2(std::get<1>(b).first->eta(),std::get<1>(b).first->phi(),std::get<1>(b).second->eta(),std::get<1>(b).second->phi());

    return std::min( dR1a, dR2a ) < std::min( dR1b, dR2b );
  };

  std::sort(casecontainer.begin(),casecontainer.end(),sortByDRDiff);

  pair1 = std::get<0>(casecontainer.front());
  pair2 = std::get<1>(casecontainer.front());

  if ( reco::deltaR2(pair1.first->eta(),pair1.first->phi(),pair1.second->eta(),pair1.second->phi())
       > reco::deltaR2(pair2.first->eta(),pair2.first->phi(),pair2.second->eta(),pair2.second->phi()) )
    std::swap(pair1,pair2);

  return true;
}
