#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

MergedLeptonHelper::MergedLeptonHelper() :
pFS_(nullptr) {
  elstr_ = TString("weight/F:pt:eta:phi:en:et:charge/I:")
  + "enSC/F:etSC:etaSC:phiSC:etaSCWidth:phiSCWidth:"
  + "full5x5_sigmaIetaIeta:full5x5_sigmaIphiIphi:"
  + "full5x5_E1x5:full5x5_E2x5:full5x5_E5x5:full5x5_hOverE:full5x5_r9:"
  + "dEtaIn:dPhiIn:dPhiSeed:dEtaEle:dPhiEle:dEtaSeed:"
  + "EseedOverP:EOverP:"
  + "ecalEn:ecalErr:trkErr:combErr:PFcombErr:"
  + "dr03TkSumPtHEEP:dr03EcalRecHitSumEt:dr03HcalDepth1TowerSumEt:"
  + "nrSatCrys/I:"
  + "modTrkIso/F:modEcalIso:"
  + "lostHits/I:nValidHits:nValidPixelHits:GsfHits:"
  + "chi2/F:d0:d0Err:dxyErr:vz:dzErr:dxy:dz:"
  + "Gsfpt/F:Gsfeta:Gsfphi:"
  + "expectedMissingInnerHits/I:"
  + "convVtxFitProb/F:convVtxChi2:convDist:convDcot:convRadius:"
  + "passConversionVeto/I:nbrem:"
  + "fbrem/F:fbremSC";

  addgsfstr_ = TString("weight/F:Gsfpt:Gsfeta:Gsfphi:")
  + "lostHits/I:nValidHits:nValidPixelHits:"
  + "chi2/F:d0:d0Err:dxyErr:vz:dzErr:dxy:dz:ptErr";

  mustr_ = TString("weight/F:numberOfValidTrackerHits/I:numberOfValidPixelHits:numberOfValidStripHits:")
  + "trackerLayersWithMeasurement:pixelLayersWithMeasurement:stripLayersWithMeasurement:"
  + "trackerLayersWithoutMeasurement:pixelLayersWithoutMeasurement:stripLayersWithoutMeasurement:"
  + "isHighPt:isTrackerHighPt:isPFloose:isPFmedium:isPFtight:"
  + "isGlobal:isTracker:isStdAlone:bestType:tunePtype:nShower:nMatchedStations:nExpectedStations:"
  + "trackerVoM/F:pixelVoM:stripVoM:pfPt:tunepPt:genPt:genEta:genPhi:PFoGen:TPoGen:trackerPt:TRKoGen:"
  + "pairPt:pairEta:pairPhi:"
  + "globalChi2:globalNormChi2:trackerChi2:trackerNormChi2:tunepChi2:tunepNormChi2";

  metstr_ = TString("weight/F:PFphi:PFdPhi:PFpt:PFoGen:PFoMu:MuEta:PFSumEt:")
  + "TPphi:TPdPhi:TPpt:TPoGen:TPoMu:TPSumEt:dSumEt:sumEtRatio:sumEtRatioTP:nLepton/I:"
  + "bitmask/i";
}

void MergedLeptonHelper::initElectronTree(const std::string& name,
                                          const std::string& prefix,
                                          const std::string& postfix) {
  elvalues_[prefix+"_"+postfix] = ElectronStruct();
  tree_[prefix+"_"+postfix+"Tree"] = (*pFS_)->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(elvalues_[prefix+"_"+postfix]),elstr_);
}

void MergedLeptonHelper::initAddGsfTree(const std::string& name,
                                        const std::string& prefix,
                                        const std::string& postfix) {
  gsfvalues_[prefix+"_"+postfix] = AddGsfStruct();
  tree_[prefix+"_"+postfix+"Tree"] = (*pFS_)->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(gsfvalues_[prefix+"_"+postfix]),addgsfstr_);
}

void MergedLeptonHelper::initMuonTree(const std::string& name,
                                      const std::string& prefix,
                                      const std::string& postfix) {
  histo1d_[prefix+"_muType"] = (*pFS_)->make<TH1D>(TString(prefix)+"_muType","muon type",3,0.,3.);
  histo1d_[prefix+"_algo"] = (*pFS_)->make<TH1D>(TString(prefix)+"_algo","muon algo",46,0.,46.);
  histo1d_[prefix+"_orgAlgo"] = (*pFS_)->make<TH1D>(TString(prefix)+"_orgAlgo","muon original algo",46,0.,46.);
  histo1d_[prefix+"_algoMask"] = (*pFS_)->make<TH1D>(TString(prefix)+"_algoMask","muon algo mask",46,0.,46.);
  histo1d_[prefix+"_quality"] = (*pFS_)->make<TH1D>(TString(prefix)+"_quality","muon quality",8,0.,8.);
  histo1d_[prefix+"_passIDs"] = (*pFS_)->make<TH1D>(TString(prefix)+"_passIDs","muon passed IDs",5,0.,5.);

  muonvalues_[prefix+"_"+postfix] = MuonStruct();
  tree_[prefix+"_"+postfix+"Tree"] = (*pFS_)->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(muonvalues_[prefix+"_"+postfix]),mustr_);
}

void MergedLeptonHelper::initMETTree(const std::string& name,
                                     const std::string& prefix,
                                     const std::string& postfix) {
  metvalues_[prefix+"_"+postfix] = METstruct();
  tree_[prefix+"_"+postfix+"Tree"] = (*pFS_)->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(metvalues_[prefix+"_"+postfix]),metstr_);
}

void MergedLeptonHelper::fillElectrons(const edm::Ptr<reco::GsfElectron>& el,
                                       const float& trkIso,
                                       const float& ecalIso,
                                       const int& nrSatCrys,
                                       const edm::Handle<reco::ConversionCollection>& conversions,
                                       const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                                       const std::string& prefix) {
  auto square = [](const double& val) { return val*val; };
  double rad = std::sqrt(square(el->superCluster()->x()) + square(el->superCluster()->y()) + square(el->superCluster()->z()));
  double radTrans = std::sqrt(square(el->superCluster()->x()) + square(el->superCluster()->y()));

  elvalues_[prefix+"_el"].weight = mcweight_;

  elvalues_[prefix+"_el"].pt = el->pt();
  elvalues_[prefix+"_el"].eta = el->eta();
  elvalues_[prefix+"_el"].phi = el->phi();
  elvalues_[prefix+"_el"].en = el->energy();
  elvalues_[prefix+"_el"].et = el->et();
  elvalues_[prefix+"_el"].charge = el->charge();

  elvalues_[prefix+"_el"].enSC = el->superCluster()->energy();
  elvalues_[prefix+"_el"].etSC = (el->superCluster()->energy())*(radTrans/rad);
  elvalues_[prefix+"_el"].etaSC = el->superCluster()->eta();
  elvalues_[prefix+"_el"].phiSC = el->superCluster()->phi();
  elvalues_[prefix+"_el"].etaSCWidth = el->superCluster()->etaWidth();
  elvalues_[prefix+"_el"].phiSCWidth = el->superCluster()->phiWidth();

  elvalues_[prefix+"_el"].full5x5_sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
  elvalues_[prefix+"_el"].full5x5_sigmaIphiIphi = el->full5x5_sigmaIphiIphi();
  elvalues_[prefix+"_el"].full5x5_E1x5 = el->full5x5_e1x5();
  elvalues_[prefix+"_el"].full5x5_E2x5 = el->full5x5_e2x5Max();
  elvalues_[prefix+"_el"].full5x5_E5x5 = el->full5x5_e5x5();
  elvalues_[prefix+"_el"].full5x5_hOverE = el->full5x5_hcalOverEcal();
  elvalues_[prefix+"_el"].full5x5_r9 = el->full5x5_r9();

  elvalues_[prefix+"_el"].dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
  elvalues_[prefix+"_el"].dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
  elvalues_[prefix+"_el"].dPhiSeed = el->deltaPhiSeedClusterTrackAtCalo();
  elvalues_[prefix+"_el"].dEtaEle = el->deltaEtaEleClusterTrackAtCalo();
  elvalues_[prefix+"_el"].dPhiEle = el->deltaPhiEleClusterTrackAtCalo();
  elvalues_[prefix+"_el"].dEtaSeed = el->deltaEtaSeedClusterTrackAtVtx();
  elvalues_[prefix+"_el"].EseedOverP = el->eSeedClusterOverP();
  elvalues_[prefix+"_el"].EOverP = el->eSuperClusterOverP();

  elvalues_[prefix+"_el"].ecalEn = el->ecalEnergy();
  elvalues_[prefix+"_el"].ecalErr = el->ecalEnergyError();
  elvalues_[prefix+"_el"].trkErr = el->trackMomentumError();
  elvalues_[prefix+"_el"].combErr = el->p4Error(reco::GsfElectron::P4_COMBINATION);
  elvalues_[prefix+"_el"].PFcombErr = el->p4Error(reco::GsfElectron::P4_PFLOW_COMBINATION);

  elvalues_[prefix+"_el"].dr03TkSumPtHEEP = el->dr03TkSumPtHEEP();
  elvalues_[prefix+"_el"].dr03EcalRecHitSumEt = el->dr03EcalRecHitSumEt();
  elvalues_[prefix+"_el"].dr03HcalDepth1TowerSumEt = el->dr03HcalDepth1TowerSumEt();
  elvalues_[prefix+"_el"].nrSatCrys = nrSatCrys;
  elvalues_[prefix+"_el"].modTrkIso = trkIso;
  elvalues_[prefix+"_el"].modEcalIso = ecalIso;

  reco::GsfTrackRef TrackRef = el->gsfTrack();
  elvalues_[prefix+"_el"].lostHits = TrackRef->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  elvalues_[prefix+"_el"].nValidHits = TrackRef->hitPattern().numberOfValidHits();
  elvalues_[prefix+"_el"].nValidPixelHits = TrackRef->hitPattern().numberOfValidPixelHits();
  elvalues_[prefix+"_el"].chi2 = TrackRef->normalizedChi2();
  elvalues_[prefix+"_el"].GsfHits = TrackRef->hitPattern().trackerLayersWithMeasurement();

  elvalues_[prefix+"_el"].d0 = TrackRef->d0();
  elvalues_[prefix+"_el"].d0Err = TrackRef->d0Error();
  elvalues_[prefix+"_el"].dxyErr = TrackRef->dxyError();
  elvalues_[prefix+"_el"].vz = TrackRef->vz();
  elvalues_[prefix+"_el"].dzErr = TrackRef->dzError();

  elvalues_[prefix+"_el"].dxy = TrackRef->dxy(pPV_->position());
  elvalues_[prefix+"_el"].dz = TrackRef->dz(pPV_->position());

  elvalues_[prefix+"_el"].Gsfpt = TrackRef->pt();
  elvalues_[prefix+"_el"].Gsfeta = TrackRef->eta();
  elvalues_[prefix+"_el"].Gsfphi = TrackRef->phi();

  elvalues_[prefix+"_el"].expectedMissingInnerHits = el->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS); // 94X

  elvalues_[prefix+"_el"].convVtxFitProb = -1.;
  elvalues_[prefix+"_el"].convVtxChi2 = -1.;
  elvalues_[prefix+"_el"].passConversionVeto = !ConversionTools::hasMatchedConversion(*el, *conversions, beamSpotHandle->position());
  const reco::Conversion* convRef = ConversionTools::matchedConversion(*el, *conversions, beamSpotHandle->position());

  if ( convRef!=nullptr ) {
    reco::Vertex convVtx = convRef->conversionVertex();
    if (convVtx.isValid()) {
      elvalues_[prefix+"_el"].convVtxFitProb = TMath::Prob(convVtx.chi2(), convVtx.ndof());
      elvalues_[prefix+"_el"].convVtxChi2 = convVtx.normalizedChi2();
    }
  }

  elvalues_[prefix+"_el"].convDist = el->convDist();
  elvalues_[prefix+"_el"].convDcot = el->convDcot();
  elvalues_[prefix+"_el"].convRadius = el->convRadius();

  elvalues_[prefix+"_el"].fbrem = el->fbrem();
  elvalues_[prefix+"_el"].nbrem = el->numberOfBrems();
  elvalues_[prefix+"_el"].fbremSC = el->superClusterFbrem();

  tree_[prefix+"_elTree"]->Fill();
}

void MergedLeptonHelper::fillGsfTracks(const reco::GsfTrackRef& addGsfTrk,
                                       const reco::GsfTrackRef& orgGsfTrk,
                                       const std::string& prefix) {
  gsfvalues_[prefix+"_addGsf"].weight = mcweight_;

  gsfvalues_[prefix+"_addGsf"].Gsfpt = addGsfTrk->pt();
  gsfvalues_[prefix+"_addGsf"].Gsfeta = addGsfTrk->eta();
  gsfvalues_[prefix+"_addGsf"].Gsfphi = addGsfTrk->phi();
  gsfvalues_[prefix+"_addGsf"].lostHits = addGsfTrk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  gsfvalues_[prefix+"_addGsf"].nValidHits = addGsfTrk->hitPattern().numberOfValidHits();
  gsfvalues_[prefix+"_addGsf"].nValidPixelHits = addGsfTrk->hitPattern().numberOfValidPixelHits();
  gsfvalues_[prefix+"_addGsf"].chi2 = addGsfTrk->normalizedChi2();
  gsfvalues_[prefix+"_addGsf"].d0 = addGsfTrk->d0();
  gsfvalues_[prefix+"_addGsf"].d0Err = addGsfTrk->d0Error();
  gsfvalues_[prefix+"_addGsf"].dxyErr = addGsfTrk->dxyError();
  gsfvalues_[prefix+"_addGsf"].vz = addGsfTrk->vz();
  gsfvalues_[prefix+"_addGsf"].dzErr = addGsfTrk->dzError();

  gsfvalues_[prefix+"_addGsf"].dxy = addGsfTrk->dxy(pPV_->position());
  gsfvalues_[prefix+"_addGsf"].dz = addGsfTrk->dz(pPV_->position());
  gsfvalues_[prefix+"_addGsf"].ptErr = addGsfTrk->ptError();

  tree_[prefix+"_addGsfTree"]->Fill();
}

void MergedLeptonHelper::fillMuons(const edm::Ptr<reco::Muon>& aMuon,
                                   const edm::Ptr<reco::GenParticle>& matched,
                                   const edm::Ptr<reco::GenParticle>& probe,
                                   std::string prefix) {
  if (aMuon->isGlobalMuon())
    histo1d_[prefix+"_muType"]->Fill(0.5);
  if (aMuon->isTrackerMuon())
    histo1d_[prefix+"_muType"]->Fill(1.5);
  if (aMuon->isStandAloneMuon())
    histo1d_[prefix+"_muType"]->Fill(2.5);

  if (muon::isHighPtMuon(*aMuon,*pPV_))
    histo1d_[prefix+"_passIDs"]->Fill(0.5);
  if (MergedLeptonIDs::isHighPtTrackerMuon(*aMuon,*pPV_))
    histo1d_[prefix+"_passIDs"]->Fill(1.5);
  if (muon::isLooseMuon(*aMuon))
    histo1d_[prefix+"_passIDs"]->Fill(2.5);
  if (muon::isMediumMuon(*aMuon))
    histo1d_[prefix+"_passIDs"]->Fill(3.5);
  if (muon::isTightMuon(*aMuon,*pPV_))
    histo1d_[prefix+"_passIDs"]->Fill(4.5);

  muonvalues_[prefix+"_muon"].weight = mcweight_;

  muonvalues_[prefix+"_muon"].isHighPt = muon::isHighPtMuon(*aMuon,*pPV_);
  muonvalues_[prefix+"_muon"].isTrackerHighPt = MergedLeptonIDs::isHighPtTrackerMuon(*aMuon,*pPV_);
  muonvalues_[prefix+"_muon"].isPFloose = muon::isLooseMuon(*aMuon);
  muonvalues_[prefix+"_muon"].isPFmedium = muon::isMediumMuon(*aMuon);
  muonvalues_[prefix+"_muon"].isPFtight = muon::isTightMuon(*aMuon,*pPV_);

  muonvalues_[prefix+"_muon"].isGlobal = aMuon->isGlobalMuon();
  muonvalues_[prefix+"_muon"].isTracker = aMuon->isTrackerMuon();
  muonvalues_[prefix+"_muon"].isStdAlone = aMuon->isStandAloneMuon();
  muonvalues_[prefix+"_muon"].bestType = static_cast<int>(aMuon->muonBestTrackType());
  muonvalues_[prefix+"_muon"].tunePtype = static_cast<int>(aMuon->tunePMuonBestTrackType());
  muonvalues_[prefix+"_muon"].nShower = aMuon->numberOfShowers();
  muonvalues_[prefix+"_muon"].nMatchedStations = aMuon->numberOfMatchedStations();
  muonvalues_[prefix+"_muon"].nExpectedStations = aMuon->expectedNnumberOfMatchedStations();

  muonvalues_[prefix+"_muon"].pfPt = aMuon->pt();
  muonvalues_[prefix+"_muon"].tunepPt = aMuon->tunePMuonBestTrack()->pt();
  muonvalues_[prefix+"_muon"].genPt = matched->pt();
  muonvalues_[prefix+"_muon"].genEta = matched->eta();
  muonvalues_[prefix+"_muon"].genPhi = matched->phi();
  muonvalues_[prefix+"_muon"].PFoGen = aMuon->pt()/matched->pt();
  muonvalues_[prefix+"_muon"].TPoGen = aMuon->tunePMuonBestTrack()->pt()/matched->pt();
  muonvalues_[prefix+"_muon"].trackerPt = aMuon->isTrackerMuon() ? aMuon->innerTrack()->pt() : -1.;
  muonvalues_[prefix+"_muon"].TRKoGen = aMuon->isTrackerMuon() ? aMuon->innerTrack()->pt()/matched->pt() : -1.;
  muonvalues_[prefix+"_muon"].pairPt = probe->pt();
  muonvalues_[prefix+"_muon"].pairEta = probe->eta();
  muonvalues_[prefix+"_muon"].pairPhi = probe->phi();
  muonvalues_[prefix+"_muon"].globalChi2 = aMuon->isGlobalMuon() ? aMuon->globalTrack()->chi2() : -1.;
  muonvalues_[prefix+"_muon"].globalNormChi2 = aMuon->isGlobalMuon() ? aMuon->globalTrack()->normalizedChi2() : -1.;
  muonvalues_[prefix+"_muon"].trackerChi2 = aMuon->isTrackerMuon() ? aMuon->innerTrack()->chi2() : -1.;
  muonvalues_[prefix+"_muon"].trackerNormChi2 = aMuon->isTrackerMuon() ? aMuon->innerTrack()->normalizedChi2() : -1.;
  muonvalues_[prefix+"_muon"].tunepChi2 = aMuon->isGlobalMuon() ? aMuon->tunePMuonBestTrack()->chi2() : -1.;
  muonvalues_[prefix+"_muon"].tunepNormChi2 = aMuon->isGlobalMuon() ? aMuon->tunePMuonBestTrack()->normalizedChi2() : -1.;

  if (!aMuon->isTrackerMuon()) {
    muonvalues_[prefix+"_muon"].numberOfValidTrackerHits = -1;
    muonvalues_[prefix+"_muon"].numberOfValidPixelHits = -1;
    muonvalues_[prefix+"_muon"].numberOfValidStripHits = -1;
    muonvalues_[prefix+"_muon"].trackerLayersWithMeasurement = -1;
    muonvalues_[prefix+"_muon"].pixelLayersWithMeasurement = -1;
    muonvalues_[prefix+"_muon"].stripLayersWithMeasurement = -1;
    muonvalues_[prefix+"_muon"].trackerLayersWithoutMeasurement = -1;
    muonvalues_[prefix+"_muon"].pixelLayersWithoutMeasurement = -1;
    muonvalues_[prefix+"_muon"].stripLayersWithoutMeasurement = -1;
    muonvalues_[prefix+"_muon"].trackerVoM = -1.;
    muonvalues_[prefix+"_muon"].pixelVoM = -1.;
    muonvalues_[prefix+"_muon"].stripVoM = -1.;
    tree_[prefix+"_muonTree"]->Fill();

    return;
  }

  auto innerTrack = aMuon->innerTrack();
  muonvalues_[prefix+"_muon"].numberOfValidTrackerHits = innerTrack->hitPattern().numberOfValidTrackerHits();
  muonvalues_[prefix+"_muon"].numberOfValidPixelHits = innerTrack->hitPattern().numberOfValidPixelHits();
  muonvalues_[prefix+"_muon"].numberOfValidStripHits = innerTrack->hitPattern().numberOfValidStripHits();
  muonvalues_[prefix+"_muon"].trackerLayersWithMeasurement = innerTrack->hitPattern().trackerLayersWithMeasurement();
  muonvalues_[prefix+"_muon"].pixelLayersWithMeasurement = innerTrack->hitPattern().pixelLayersWithMeasurement();
  muonvalues_[prefix+"_muon"].stripLayersWithMeasurement = innerTrack->hitPattern().stripLayersWithMeasurement();
  muonvalues_[prefix+"_muon"].trackerLayersWithoutMeasurement = innerTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  muonvalues_[prefix+"_muon"].pixelLayersWithoutMeasurement = innerTrack->hitPattern().pixelLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  muonvalues_[prefix+"_muon"].stripLayersWithoutMeasurement = innerTrack->hitPattern().stripLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  muonvalues_[prefix+"_muon"].trackerVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidTrackerHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
  muonvalues_[prefix+"_muon"].pixelVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidPixelHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());
  muonvalues_[prefix+"_muon"].stripVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidStripHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().stripLayersWithMeasurement());
  tree_[prefix+"_muonTree"]->Fill();

  histo1d_[prefix+"_algo"]->Fill(static_cast<float>(innerTrack->algo())+0.5);
  histo1d_[prefix+"_orgAlgo"]->Fill(static_cast<float>(innerTrack->originalAlgo())+0.5);

  for (int en = reco::TrackBase::undefAlgorithm; en != reco::TrackBase::algoSize; ++en) {
    if ( (innerTrack->algoMask())[en] )
      histo1d_[prefix+"_algoMask"]->Fill(static_cast<float>(en)+0.5);
  }

  for (int en = reco::TrackBase::loose; en != reco::TrackBase::qualitySize; ++en) {
    if ( innerTrack->quality( static_cast<reco::TrackBase::TrackQuality>(en) ) )
      histo1d_[prefix+"_quality"]->Fill(static_cast<float>(en)+0.5);
  }
}

void MergedLeptonHelper::fillMETs(const edm::Ptr<reco::Muon>& aMuon,
                                  const edm::Ptr<reco::GenParticle>& matched,
                                  const pat::MET& met,
                                  const std::vector<edm::Ptr<reco::Muon>>& finalMuons,
                                  const unsigned int& bitmask,
                                  const std::string& prefix) {
  metvalues_[prefix+"_MET"].weight = mcweight_;
  metvalues_[prefix+"_MET"].PFphi = met.phi();
  metvalues_[prefix+"_MET"].PFdPhi = reco::deltaPhi(met.phi(),matched->phi());
  metvalues_[prefix+"_MET"].PFpt = met.pt();
  metvalues_[prefix+"_MET"].PFoGen = met.pt()/matched->pt();
  metvalues_[prefix+"_MET"].PFoMu = met.pt()/aMuon->pt();
  metvalues_[prefix+"_MET"].MuEta = aMuon->eta();
  metvalues_[prefix+"_MET"].PFSumEt = met.sumEt();
  metvalues_[prefix+"_MET"].bitmask = bitmask;

  pat::MET::Vector2 tpvec;
  tpvec.px = met.px() + aMuon->px() - aMuon->tunePMuonBestTrack()->px();
  tpvec.py = met.py() + aMuon->py() - aMuon->tunePMuonBestTrack()->py();
  metvalues_[prefix+"_MET"].TPphi = tpvec.phi();
  metvalues_[prefix+"_MET"].TPdPhi = reco::deltaPhi(tpvec.phi(),matched->phi());
  metvalues_[prefix+"_MET"].TPpt = tpvec.pt();
  metvalues_[prefix+"_MET"].TPoGen = tpvec.pt()/matched->pt();
  metvalues_[prefix+"_MET"].TPoMu = tpvec.pt()/aMuon->tunePMuonBestTrack()->pt();

  float sumpt = 0.;
  float sumTPpt = 0.;

  for (auto mu : finalMuons) {
    sumpt += mu->pt();
    sumTPpt += mu->tunePMuonBestTrack()->pt();
  }

  metvalues_[prefix+"_MET"].TPSumEt = met.sumEt() + sumTPpt - sumpt;
  metvalues_[prefix+"_MET"].dSumEt = met.sumEt() - sumpt;
  metvalues_[prefix+"_MET"].sumEtRatio = met.sumEt()/sumpt;
  metvalues_[prefix+"_MET"].sumEtRatioTP = metvalues_[prefix+"_MET"].TPSumEt/sumTPpt;
  metvalues_[prefix+"_MET"].nLepton = static_cast<int>(finalMuons.size());

  tree_[prefix+"_METTree"]->Fill();
}

void MergedLeptonHelper::fillMETs(const edm::Ptr<reco::Candidate>& aLepton,
                                  const pat::MET& met,
                                  const std::vector<edm::Ptr<reco::Candidate>>& finalLeptons,
                                  const unsigned int& bitmask,
                                  const std::string& prefix) {
  metvalues_[prefix+"_MET"].weight = mcweight_;
  metvalues_[prefix+"_MET"].PFphi = met.phi();
  metvalues_[prefix+"_MET"].PFdPhi = reco::deltaPhi(met.phi(),aLepton->phi());
  metvalues_[prefix+"_MET"].PFpt = met.pt();
  metvalues_[prefix+"_MET"].PFoGen = 1.;
  metvalues_[prefix+"_MET"].PFoMu = met.pt()/aLepton->pt();
  metvalues_[prefix+"_MET"].MuEta = aLepton->eta();
  metvalues_[prefix+"_MET"].PFSumEt = met.sumEt();
  metvalues_[prefix+"_MET"].bitmask = bitmask;

  pat::MET::Vector2 tpvec;
  tpvec.px = met.px();
  tpvec.py = met.py();
  metvalues_[prefix+"_MET"].TPphi = tpvec.phi();
  metvalues_[prefix+"_MET"].TPdPhi = reco::deltaPhi(met.phi(),aLepton->phi());
  metvalues_[prefix+"_MET"].TPpt = tpvec.pt();
  metvalues_[prefix+"_MET"].TPoGen = 1.;
  metvalues_[prefix+"_MET"].TPoMu = tpvec.pt()/aLepton->pt();

  float sumpt = 0.;
  float sumTPpt = 0.;

  for (auto lep : finalLeptons) {
    sumpt += lep->pt();
    sumTPpt += lep->pt();
  }

  metvalues_[prefix+"_MET"].TPSumEt = met.sumEt() + sumTPpt - sumpt;
  metvalues_[prefix+"_MET"].dSumEt = met.sumEt() - sumpt;
  metvalues_[prefix+"_MET"].sumEtRatio = met.sumEt()/sumpt;
  metvalues_[prefix+"_MET"].sumEtRatioTP = metvalues_[prefix+"_MET"].TPSumEt/sumTPpt;
  metvalues_[prefix+"_MET"].nLepton = static_cast<int>(finalLeptons.size());

  tree_[prefix+"_METTree"]->Fill();
}
