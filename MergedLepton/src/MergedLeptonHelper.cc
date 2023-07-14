#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

bool MergedLeptonHelperFct::isNotMerged(const pat::ElectronRef& aEle,
                                        const edm::Handle<edm::View<pat::Electron>>& eleHandle,
                                        const reco::GsfTrackRef& addGsfTrk) {
  // additional GSF track must not be an electron
  bool notMerged = false;

  // additional GSF track is not selected
  if (aEle->gsfTrack()==addGsfTrk)
    return notMerged;

  // GSF track must be ambiguous
  for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
    const auto& secEle = eleHandle->refAt(jdx);
    const auto& secGsfTrk = secEle->gsfTrack();

    if (aEle==secEle.castTo<pat::ElectronRef>())
      continue;

    if ( addGsfTrk==secGsfTrk ) {
      notMerged = true;
      break;
    }
  }

  return notMerged;
}

MergedLeptonHelper::MergedLeptonHelper()
: pFS_(nullptr) {
  elstr_ = TString("weight/F:pt:eta:phi:en:et:charge/I:")
  + "enSC/F:etSC:etaSC:phiSC:etaSeed:phiSeed:etaSCWidth:phiSCWidth:"
  + "full5x5_sigmaIetaIeta:full5x5_sigmaIphiIphi:"
  + "full5x5_E1x5:full5x5_E2x5:full5x5_E5x5:full5x5_hOverE:full5x5_r9:"
  + "dEtaIn:dPhiIn:dPhiSeed:dEtaEle:dPhiEle:dEtaSeed:"
  + "EseedOverP:EOverP:"
  + "ecalEn:ecalErr:trkErr:combErr:PFcombErr:"
  + "dr03TkSumPtHEEP:dr03EcalRecHitSumEt:dr03HcalDepth1TowerSumEt:"
  + "modTrkIso/F:modEcalIso:"
  + "lostHits/I:nValidHits:nValidPixelHits:GsfHits:"
  + "chi2/F:d0:d0Err:dxyErr:vz:dzErr:dxy:dz:"
  + "Gsfpt/F:Gsfeta:Gsfphi:GsfPtErr:"
  + "expectedMissingInnerHits/I:"
  + "convDist/F:convDcot:convRadius:"
  + "passConversionVeto/I:nbrem:"
  + "fbrem/F:fbremSC:"
  + "full5x5_e2x5Left:full5x5_e2x5Right:full5x5_e2x5Top:full5x5_e2x5Bottom:"
  + "full5x5_eLeft:full5x5_eRight:full5x5_eTop:full5x5_eBottom:"
  + "full5x5_eMax:full5x5_e2nd:"
  + "clus2ndMoment_sMaj:clus2ndMoment_sMin:clus2ndMoment_alpha:"
  + "dEtaSeedMiniAOD:dPhiInMiniAOD:sigIeIeMiniAOD:"
  + "union5x5dEtaIn:union5x5dPhiIn:"
  + "union5x5Energy:union5x5covIeIe:union5x5covIeIp:union5x5covIpIp:"
  + "union5x5covMaj:union5x5covMin:"
  + "alphaCalo";

  addtrkstr_ = TString("weight/F:pt:eta:phi:")
  + "lostHits/I:nValidHits:nValidPixelHits:"
  + "chi2/F:d0:d0Err:dxyErr:vz:dzErr:dxy:dz:ptErr:"
  + "deltaEtaSeedClusterAtVtx:deltaPhiSuperClusterAtVtx:"
  + "dPerpIn:normalizedDParaIn:alphaTrack";
}

void MergedLeptonHelper::initElectronTree(const std::string& name,
                                          const std::string& prefix,
                                          const std::string& postfix) {
  elvalues_[prefix+"_"+postfix] = ElectronStruct();
  tree_[prefix+"_"+postfix+"Tree"] = (*pFS_)->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(elvalues_[prefix+"_"+postfix]),elstr_);
}

void MergedLeptonHelper::initAddTrkTree(const std::string& name,
                                        const std::string& prefix,
                                        const std::string& postfix) {
  trkvalues_[prefix+"_"+postfix] = AddTrkStruct();
  tree_[prefix+"_"+postfix+"Tree"] = (*pFS_)->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(trkvalues_[prefix+"_"+postfix]),addtrkstr_);
}

void MergedLeptonHelper::fillElectrons(const pat::ElectronRef& el,
                                       const float& trkIso,
                                       const float& ecalIso,
                                       const ModifiedShowerShape::variables& variables,
                                       const EcalRecHitCollection* ecalRecHits,
                                       const edm::EventSetup& iSetup,
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
  elvalues_[prefix+"_el"].etaSeed = el->superCluster()->seed()->eta();
  elvalues_[prefix+"_el"].phiSeed = el->superCluster()->seed()->phi();
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
  elvalues_[prefix+"_el"].modTrkIso = trkIso;
  elvalues_[prefix+"_el"].modEcalIso = ecalIso;

  const auto GSFtrack = el->gsfTrack();
  elvalues_[prefix+"_el"].lostHits = GSFtrack->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  elvalues_[prefix+"_el"].nValidHits = GSFtrack->hitPattern().numberOfValidHits();
  elvalues_[prefix+"_el"].nValidPixelHits = GSFtrack->hitPattern().numberOfValidPixelHits();
  elvalues_[prefix+"_el"].chi2 = GSFtrack->normalizedChi2();
  elvalues_[prefix+"_el"].GsfHits = GSFtrack->hitPattern().trackerLayersWithMeasurement();

  elvalues_[prefix+"_el"].d0 = GSFtrack->d0();
  elvalues_[prefix+"_el"].d0Err = GSFtrack->d0Error();
  elvalues_[prefix+"_el"].dxyErr = GSFtrack->dxyError();
  elvalues_[prefix+"_el"].vz = GSFtrack->vz();
  elvalues_[prefix+"_el"].dzErr = GSFtrack->dzError();

  elvalues_[prefix+"_el"].dxy = GSFtrack->dxy(pPV_->position());
  elvalues_[prefix+"_el"].dz = GSFtrack->dz(pPV_->position());

  elvalues_[prefix+"_el"].Gsfpt = GSFtrack->pt();
  elvalues_[prefix+"_el"].Gsfeta = GSFtrack->eta();
  elvalues_[prefix+"_el"].Gsfphi = GSFtrack->phi();
  elvalues_[prefix+"_el"].GsfPtErr = GSFtrack->ptError();

  elvalues_[prefix+"_el"].expectedMissingInnerHits = el->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS); // 94X

  elvalues_[prefix+"_el"].passConversionVeto = el->passConversionVeto();
  elvalues_[prefix+"_el"].convDist = el->convDist();
  elvalues_[prefix+"_el"].convDcot = el->convDcot();
  elvalues_[prefix+"_el"].convRadius = el->convRadius();

  elvalues_[prefix+"_el"].fbrem = el->fbrem();
  elvalues_[prefix+"_el"].nbrem = el->numberOfBrems();
  elvalues_[prefix+"_el"].fbremSC = el->superClusterFbrem();

  elvalues_[prefix+"_el"].full5x5_eMax = el->full5x5_showerShape().eMax;
  elvalues_[prefix+"_el"].full5x5_e2nd = el->full5x5_showerShape().e2nd;

  elvalues_[prefix+"_el"].full5x5_e2x5Left = el->full5x5_e2x5Left();
  elvalues_[prefix+"_el"].full5x5_e2x5Right = el->full5x5_e2x5Right();
  elvalues_[prefix+"_el"].full5x5_e2x5Top = el->full5x5_e2x5Top();
  elvalues_[prefix+"_el"].full5x5_e2x5Bottom = el->full5x5_e2x5Bottom();
  elvalues_[prefix+"_el"].full5x5_eLeft = el->full5x5_eLeft();
  elvalues_[prefix+"_el"].full5x5_eRight = el->full5x5_eRight();
  elvalues_[prefix+"_el"].full5x5_eTop = el->full5x5_eTop();
  elvalues_[prefix+"_el"].full5x5_eBottom = el->full5x5_eBottom();

  std::vector<std::pair<const EcalRecHit*,float>> recHitsAndFractions;

  for (auto& xtal : el->superCluster()->hitsAndFractions()) {
    EcalRecHitCollection::const_iterator aRecHit = ecalRecHits->find(xtal.first);

    if ( aRecHit!=ecalRecHits->end() )
      recHitsAndFractions.push_back( std::make_pair(&(*aRecHit),xtal.second) );
  }

  const reco::CaloClusterPtr seedClus = el->superCluster()->seed();
  auto cluster2ndMoments = noZS::EcalClusterTools::cluster2ndMoments(recHitsAndFractions,1.0);
  elvalues_[prefix+"_el"].clus2ndMoment_sMaj = cluster2ndMoments.sMaj;
  elvalues_[prefix+"_el"].clus2ndMoment_sMin = cluster2ndMoments.sMin;
  elvalues_[prefix+"_el"].clus2ndMoment_alpha = cluster2ndMoments.alpha;

  elvalues_[prefix+"_el"].dEtaSeedMiniAOD = std::numeric_limits<float>::max();
  elvalues_[prefix+"_el"].dPhiInMiniAOD = std::numeric_limits<float>::max();

  auto dEtaInSeedCalc = ModifiedDEtaInSeed(posCalcLog_);
  auto scAtVtx = EleRelPointPair(math::XYZPoint(),math::XYZPoint(),pBS_->position());
  auto seedAtCalo = EleRelPointPair(math::XYZPoint(),math::XYZPoint(),pBS_->position());

  if ( dEtaInSeedCalc.extrapolate(*el,*(el->gsfTrack()),pBS_->position(),iSetup,scAtVtx,seedAtCalo) ) {
    elvalues_[prefix+"_el"].dEtaSeedMiniAOD = scAtVtx.dEta() - el->superCluster()->eta() + el->superCluster()->seed()->eta();
    elvalues_[prefix+"_el"].dPhiInMiniAOD = scAtVtx.dPhi();
  }

  elvalues_[prefix+"_el"].union5x5dEtaIn = variables.dEtaInUnion5x5;
  elvalues_[prefix+"_el"].union5x5dPhiIn = variables.dPhiInUnion5x5;
  elvalues_[prefix+"_el"].union5x5Energy = variables.union5x5Energy;
  elvalues_[prefix+"_el"].union5x5covIeIe = variables.covEE;
  elvalues_[prefix+"_el"].union5x5covIeIp = variables.covEP;
  elvalues_[prefix+"_el"].union5x5covIpIp = variables.covPP;
  elvalues_[prefix+"_el"].alphaCalo = variables.alpha;

  double covEE = variables.covEE;
  double covEP = variables.covEP;
  double covPP = variables.covPP;
  double covMaj = ((covEE + covPP) + std::sqrt((covEE - covPP)*(covEE - covPP) + 4.*covEP*covEP)) / 2.;
  double covMin = ((covEE + covPP) - std::sqrt((covEE - covPP)*(covEE - covPP) + 4.*covEP*covEP)) / 2.;
  elvalues_[prefix+"_el"].union5x5covMaj = covMaj;
  elvalues_[prefix+"_el"].union5x5covMin = covMin;

  // Get Calo Topology
  edm::ESHandle<CaloTopology> caloTopoHandle;
  iSetup.get<CaloTopologyRecord>().get(caloTopoHandle);
  const CaloTopology* caloTopo = caloTopoHandle.product();

  auto matrix5x5of1stGSF = noZS::EcalClusterTools::matrixDetId(caloTopo, el->superCluster()->seed()->seed(), 2 );
  std::vector<DetId> unionMatrix5x5(matrix5x5of1stGSF);
  std::vector<std::pair<DetId,float>> hitFracUnion5x5;

  for (const auto& adetId : unionMatrix5x5)
    hitFracUnion5x5.push_back(std::make_pair(adetId,1.));

  auto showerShapeCalc = ModifiedShowerShape(posCalcLog_);

  double sumE = noZS::EcalClusterTools::e5x5(*(el->superCluster()->seed()),ecalRecHits,caloTopo);
  auto sigmas = showerShapeCalc.calcSigmas(*el,hitFracUnion5x5,ecalRecHits,sumE);
  double sigEE = std::get<0>(sigmas);
  elvalues_[prefix+"_el"].sigIeIeMiniAOD = 0.01745*std::sqrt(sigEE);

  tree_[prefix+"_elTree"]->Fill();
}

void MergedLeptonHelper::fillAddTracks(const pat::ElectronRef& el,
                                       const reco::TrackBase* addTrk,
                                       const ModifiedDEtaInSeed::variables& variables,
                                       const EcalRecHitCollection* ecalRecHits,
                                       const edm::EventSetup& iSetup,
                                       const std::string& prefix) {
  trkvalues_[prefix+"_addTrk"].weight = mcweight_;

  trkvalues_[prefix+"_addTrk"].pt = addTrk->pt();
  trkvalues_[prefix+"_addTrk"].eta = addTrk->eta();
  trkvalues_[prefix+"_addTrk"].phi = addTrk->phi();
  trkvalues_[prefix+"_addTrk"].lostHits = addTrk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  trkvalues_[prefix+"_addTrk"].nValidHits = addTrk->hitPattern().numberOfValidHits();
  trkvalues_[prefix+"_addTrk"].nValidPixelHits = addTrk->hitPattern().numberOfValidPixelHits();
  trkvalues_[prefix+"_addTrk"].chi2 = addTrk->normalizedChi2();
  trkvalues_[prefix+"_addTrk"].d0 = addTrk->d0();
  trkvalues_[prefix+"_addTrk"].d0Err = addTrk->d0Error();
  trkvalues_[prefix+"_addTrk"].dxyErr = addTrk->dxyError();
  trkvalues_[prefix+"_addTrk"].vz = addTrk->vz();
  trkvalues_[prefix+"_addTrk"].dzErr = addTrk->dzError();

  trkvalues_[prefix+"_addTrk"].dxy = addTrk->dxy(pPV_->position());
  trkvalues_[prefix+"_addTrk"].dz = addTrk->dz(pPV_->position());
  trkvalues_[prefix+"_addTrk"].ptErr = addTrk->ptError();

  trkvalues_[prefix+"_addTrk"].deltaEtaSeedClusterAtVtx = variables.dEtaInSeed2nd;
  trkvalues_[prefix+"_addTrk"].deltaPhiSuperClusterAtVtx = variables.dPhiInSC2nd;
  trkvalues_[prefix+"_addTrk"].dPerpIn = variables.dPerpIn;
  trkvalues_[prefix+"_addTrk"].alphaTrack = variables.alphaTrack;
  trkvalues_[prefix+"_addTrk"].normalizedDParaIn = variables.normalizedDParaIn;

  tree_[prefix+"_addTrkTree"]->Fill();
}
