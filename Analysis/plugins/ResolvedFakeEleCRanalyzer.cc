#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"

#include "correction.h"

class ResolvedFakeEleCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedFakeEleCRanalyzer(const edm::ParameterSet&);
  virtual ~ResolvedFakeEleCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> modifiedTrkIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<double> rhoToken_;

  std::unique_ptr<correction::CorrectionSet> purwgt_;
  const std::string puname_;

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  std::map<std::string,TH1*> histo1d_;

  TTree* tree_1m_denom_ = nullptr;
  TTree* tree_1m_numer_ = nullptr;
  TTree* tree_2m_denom_ = nullptr;
  TTree* tree_2m_numer_ = nullptr;
  float var_weight_ = 0.;
  float var_1m_muPt_ = -1.;
  float var_1m_muEta_ = std::numeric_limits<float>::max();
  float var_1m_muPhi_ = std::numeric_limits<float>::max();
  float var_1m_met_ = -1.;
  float var_1m_metPhi_ = std::numeric_limits<float>::max();
  float var_2m_muPt_ = -1.;
  float var_2m_muEta_ = std::numeric_limits<float>::max();
  float var_2m_muPhi_ = std::numeric_limits<float>::max();
  float var_e1Pt_ = -1.;
  float var_e1scEn_ = -1.;
  float var_e1scEta_ = std::numeric_limits<float>::max();
  float var_e1Eta_ = std::numeric_limits<float>::max();
  float var_e1Phi_ = std::numeric_limits<float>::max();
  float var_e1DPerpIn_ = std::numeric_limits<float>::max();
  float var_e1DEtaInSeed_ = std::numeric_limits<float>::max();
  float var_e1MvaScore_ = -1.;
  int var_e1MvaCategory_ = -1;
  int var_e1PassMva_ = -1;
  float var_e1TrkIso_ = -1.;
  float var_e1EcalIso_ = -1.;
  float var_e1HcalIso_ = -1.;
  float var_e1HoE_ = -1.;
  float var_e2Pt_ = -1.;
  float var_e2scEn_ = -1.;
  float var_e2scEta_ = std::numeric_limits<float>::max();
  float var_e2Eta_ = std::numeric_limits<float>::max();
  float var_e2Phi_ = std::numeric_limits<float>::max();
  float var_e2DPerpIn_ = std::numeric_limits<float>::max();
  float var_e2DEtaInSeed_ = std::numeric_limits<float>::max();
  float var_e2MvaScore_ = -1.;
  int var_e2MvaCategory_ = -1;
  int var_e2PassMva_ = -1;
  float var_e2TrkIso_ = -1.;
  float var_e2EcalIso_ = -1.;
  float var_e2HcalIso_ = -1.;
  float var_e2HoE_ = -1.;
  float var_e1e2Rho_ = -1.;
  float var_e1e2InvM_ = -1.;
  float var_e1e2DR_ = -1.;
  float var_e1e2scDR_ = -1.;
};

ResolvedFakeEleCRanalyzer::ResolvedFakeEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
modifiedTrkIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("modifiedTrkIso"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("modifiedEcalIso"))),
dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
purwgt_(std::move(correction::CorrectionSet::from_file((iConfig.getParameter<edm::FileInPath>("PUrwgt")).fullPath()))),
puname_(iConfig.getParameter<std::string>("PUname")),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")) {
  usesResource("TFileService");
}

void ResolvedFakeEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  tree_1m_denom_ = fs->make<TTree>("denomTree1m","denomTree1m");
  tree_1m_denom_->Branch("wgt",&var_weight_,"wgt/F");
  tree_1m_denom_->Branch("muPt",&var_1m_muPt_,"muPt/F");
  tree_1m_denom_->Branch("muEta",&var_1m_muEta_,"muEta/F");
  tree_1m_denom_->Branch("muPhi",&var_1m_muPhi_,"muPhi/F");
  tree_1m_denom_->Branch("met",&var_1m_met_,"met/F");
  tree_1m_denom_->Branch("metPhi",&var_1m_metPhi_,"metPhi/F");
  tree_1m_numer_ = fs->make<TTree>("numerTree1m","numerTree1m");
  tree_1m_numer_->Branch("wgt",&var_weight_,"wgt/F");
  tree_1m_numer_->Branch("muPt",&var_1m_muPt_,"muPt/F");
  tree_1m_numer_->Branch("muEta",&var_1m_muEta_,"muEta/F");
  tree_1m_numer_->Branch("muPhi",&var_1m_muPhi_,"muPhi/F");
  tree_1m_numer_->Branch("met",&var_1m_met_,"met/F");
  tree_1m_numer_->Branch("metPhi",&var_1m_metPhi_,"metPhi/F");

  auto setEleBranch = [this] (TTree* atree) {
    atree->Branch("e1Pt",&var_e1Pt_,"e1Pt/F");
    atree->Branch("e1scEn",&var_e1scEn_,"e1scEn/F");
    atree->Branch("e1scEta",&var_e1scEta_,"e1scEta/F");
    atree->Branch("e1Eta",&var_e1Eta_,"e1Eta/F");
    atree->Branch("e1Phi",&var_e1Phi_,"e1Phi/F");
    atree->Branch("e1DPerpIn",&var_e1DPerpIn_,"e1DPerpIn/F");
    atree->Branch("e1DEtaInSeed",&var_e1DEtaInSeed_,"e1DEtaInSeed/F");
    atree->Branch("e1MvaScore",&var_e1MvaScore_,"e1MvaScore/F");
    atree->Branch("e1MvaCategory",&var_e1MvaCategory_,"e1MvaCategory/I");
    atree->Branch("e1PassMva",&var_e1PassMva_,"e1PassMva/I");
    atree->Branch("e1TrkIso",&var_e1TrkIso_,"e1TrkIso/F");
    atree->Branch("e1EcalIso",&var_e1EcalIso_,"e1EcalIso/F");
    atree->Branch("e1HcalIso",&var_e1HcalIso_,"e1HcalIso/F");
    atree->Branch("e1HoE",&var_e1HoE_,"e1HoE/F");

    atree->Branch("e2Pt",&var_e2Pt_,"e2Pt/F");
    atree->Branch("e2scEn",&var_e2scEn_,"e2scEn/F");
    atree->Branch("e2scEta",&var_e2scEta_,"e2scEta/F");
    atree->Branch("e2Eta",&var_e2Eta_,"e2Eta/F");
    atree->Branch("e2Phi",&var_e2Phi_,"e2Phi/F");
    atree->Branch("e2DPerpIn",&var_e2DPerpIn_,"e2DPerpIn/F");
    atree->Branch("e2DEtaInSeed",&var_e2DEtaInSeed_,"e2DEtaInSeed/F");
    atree->Branch("e2MvaScore",&var_e2MvaScore_,"e2MvaScore/F");
    atree->Branch("e2MvaCategory",&var_e2MvaCategory_,"e2MvaCategory/I");
    atree->Branch("e2PassMva",&var_e2PassMva_,"e2PassMva/I");
    atree->Branch("e2TrkIso",&var_e2TrkIso_,"e2TrkIso/F");
    atree->Branch("e2EcalIso",&var_e2EcalIso_,"e2EcalIso/F");
    atree->Branch("e2HcalIso",&var_e2HcalIso_,"e2HcalIso/F");
    atree->Branch("e2HoE",&var_e2HoE_,"e2HoE/F");

    atree->Branch("e1e2Rho",&var_e1e2Rho_,"e1e2Rho/F");
    atree->Branch("e1e2InvM",&var_e1e2InvM_,"e1e2InvM/F");
    atree->Branch("e1e2DR",&var_e1e2DR_,"e1e2DR/F");
    atree->Branch("e1e2scDR",&var_e1e2scDR_,"e1e2scDR/F");
  };

  setEleBranch(tree_1m_denom_);
  setEleBranch(tree_1m_numer_);

  tree_2m_denom_ = fs->make<TTree>("denomTree2m","denomTree2m");
  tree_2m_denom_->Branch("wgt",&var_weight_,"wgt/F");
  tree_2m_denom_->Branch("mu1Pt",&var_1m_muPt_,"mu1Pt/F");
  tree_2m_denom_->Branch("mu1Eta",&var_1m_muEta_,"mu1Eta/F");
  tree_2m_denom_->Branch("mu1Phi",&var_1m_muPhi_,"mu1Phi/F");
  tree_2m_denom_->Branch("mu2Pt",&var_2m_muPt_,"mu2Pt/F");
  tree_2m_denom_->Branch("mu2Eta",&var_2m_muEta_,"mu2Eta/F");
  tree_2m_denom_->Branch("mu2Phi",&var_2m_muPhi_,"mu2Phi/F");
  tree_2m_numer_ = fs->make<TTree>("numerTree2m","numerTree2m");
  tree_2m_numer_->Branch("wgt",&var_weight_,"wgt/F");
  tree_2m_numer_->Branch("mu1Pt",&var_1m_muPt_,"mu1Pt/F");
  tree_2m_numer_->Branch("mu1Eta",&var_1m_muEta_,"mu1Eta/F");
  tree_2m_numer_->Branch("mu1Phi",&var_1m_muPhi_,"mu1Phi/F");
  tree_2m_numer_->Branch("mu2Pt",&var_2m_muPt_,"mu2Pt/F");
  tree_2m_numer_->Branch("mu2Eta",&var_2m_muEta_,"mu2Eta/F");
  tree_2m_numer_->Branch("mu2Phi",&var_2m_muPhi_,"mu2Phi/F");

  setEleBranch(tree_2m_denom_);
  setEleBranch(tree_2m_numer_);
}

void ResolvedFakeEleCRanalyzer::endJob() {}

void ResolvedFakeEleCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  double aWeight = 1.;
  double purwgtNo = 1.;
  double prefireNo = 1.;

  if (isMC_) {
    edm::Handle<double> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    prefireNo = *theprefweight;

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    aWeight = prefireNo*mcweight/std::abs(mcweight);

    edm::Handle<edm::View<PileupSummaryInfo>> pusummary;
    iEvent.getByToken(pileupToken_, pusummary);

    for (unsigned int idx = 0; idx < pusummary->size(); ++idx) {
      const auto& apu = pusummary->refAt(idx);

      int bx = apu->getBunchCrossing();

      if (bx==0) { // in-time PU only
        purwgtNo = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"nominal"});

        aWeight *= purwgtNo;

        break;
      }
    }
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  std::vector<std::pair<std::string, int>> indices;
  edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    std::string trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
        if (trigResultHandle.product()->accept(iTrig))
          isFired = true;
      }
    } // wanted triggers
  } // fired triggers

  if (!isFired)
    return;

  edm::Handle<edm::TriggerResults> METfilterHandle;
  iEvent.getByToken(METfilterToken_,METfilterHandle);
  edm::TriggerNames METfilters = iEvent.triggerNames(*METfilterHandle);

  unsigned int nPassedFilters = 0;

  for (unsigned int iTrig = 0; iTrig < METfilterHandle.product()->size(); iTrig++) {
    const std::string trigname = METfilters.triggerName(iTrig);

    if (METfilterHandle.product()->accept(iTrig)) {
      for (const auto& filterName : METfilterList_) {
        if (trigname.find(filterName) != std::string::npos)
          nPassedFilters++;
      }
    }
  }

  if (nPassedFilters!=METfilterList_.size())
    return;

  std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjs;

  for (unsigned iTrig = 0; iTrig < trigObjHandle->size(); iTrig++) {
    const auto& trigObj = trigObjHandle->refAt(iTrig);
    auto trigObjInst = trigObjHandle->at(iTrig); // workaround for copy
    trigObjInst.unpackPathNames(trigList);
    const auto& pathNames = trigObjInst.pathNames();

    for (const auto name : pathNames) {
      for (unsigned int jTrig = 0; jTrig < trigList_.size(); jTrig++) {
        if ( name.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos &&
             trigObjInst.hasPathName(name,true,true) ) {
          trigObjs.push_back(trigObj);
        }
      } // wanted triggers
    } // fired triggers
  } // trigger objs

  if ( !pvHandle.isValid() || pvHandle->empty() )
    return;

  const reco::Vertex* primaryVertex = pvHandle->ptrAt(0).get();

  std::vector<pat::MuonRef> isolatedHighPtMuons;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->tunePMuonBestTrack()->pt() < 20. || std::abs(aMuon->eta()) > 2.4 )
      continue;

    if ( muon::isHighPtMuon(*aMuon,*primaryVertex) ) {
      if ( aMuon->trackIso()/aMuon->tunePMuonBestTrack()->pt() < 0.1 )
        isolatedHighPtMuons.push_back(aMuon.castTo<pat::MuonRef>());
    }
  }

  if ( isolatedHighPtMuons.empty() )
    return;

  auto sortByTuneP = [](const pat::MuonRef& a, const pat::MuonRef& b) {
    return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
  };

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);

  bool trigMatched = false;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    const auto& leadMu = isolatedHighPtMuons.front();

    if ( leadMu->tunePMuonBestTrack()->pt() > 30. &&
         reco::deltaR2(trigObj->eta(),
                       trigObj->phi(),
                       leadMu->tunePMuonBestTrack()->eta(),
                       leadMu->tunePMuonBestTrack()->phi()) < 0.01 ) {
      trigMatched = true;

      // if (isMC_)
      //   aWeight *= mucorrHelper_.trigSFglobal(leadMu);

      break;
    }
  }

  if ( !trigMatched )
    return;

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<float>> modifiedTrkIsoHandle;
  iEvent.getByToken(modifiedTrkIsoToken_, modifiedTrkIsoHandle);

  edm::Handle<edm::ValueMap<float>> modifiedEcalIsoHandle;
  iEvent.getByToken(ecalIsoToken_, modifiedEcalIsoHandle);

  edm::Handle<edm::ValueMap<float>> dPerpInHandle;
  iEvent.getByToken(dPerpInToken_, dPerpInHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);

  std::vector<pat::ElectronRef> acceptEles;
  std::vector<pat::ElectronRef> nonHeepEles;

  auto sortByEt = [](const pat::ElectronRef& a, const pat::ElectronRef& b) {
    return a->et() > b->et();
  };

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( std::abs(aEle->superCluster()->eta()) > 2.5 )
      continue;

    // veto EBEE gap
    if ( std::abs(aEle->superCluster()->eta()) > 1.4442 && std::abs(aEle->superCluster()->eta()) < 1.566 )
      continue;

    bool selected = false;

    if ( aEle->electronID("modifiedHeepElectronID") ) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      acceptEles.push_back(castEle);
      selected = true;
    }

    if (!selected) {
      auto castEle = aEle.castTo<pat::ElectronRef>();

      int32_t bitmap = aEle->userInt("modifiedHeepElectronID");
      int32_t mask = 0x000001F4; // 0001 1111 0100
      int32_t pass = bitmap | mask;
      bool passMaskedId = pass==0x00000FFF; // HEEP ID has 12 cuts
      bool passDEtaIn = (std::abs((*dPerpInHandle)[aEle]) < 0.01 || std::abs(aEle->deltaEtaSeedClusterTrackAtVtx()) < 0.01);

      if (passMaskedId && passDEtaIn)
        nonHeepEles.push_back(castEle);
    }
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByEt);
  std::sort(nonHeepEles.begin(),nonHeepEles.end(),sortByEt);

  auto fillEleBranch = [&,this] (const pat::ElectronRef& e1, const pat::ElectronRef& e2) {
    var_e1Pt_ = e1->pt();
    var_e1scEn_ = e1->superCluster()->energy();
    var_e1scEta_ = e1->superCluster()->eta();
    var_e1Eta_ = e1->eta();
    var_e1Phi_ = e1->phi();
    var_e1DPerpIn_ = (*dPerpInHandle)[e1];
    var_e1DEtaInSeed_ = e1->deltaEtaSeedClusterTrackAtVtx();
    var_e1MvaScore_ = e1->userFloat("mvaMergedElectronValues");
    var_e1MvaCategory_ = e1->userInt("mvaMergedElectronCategories");
    var_e1PassMva_ = e1->electronID("mvaMergedElectron");
    var_e1TrkIso_ = (*modifiedTrkIsoHandle)[e1];
    var_e1EcalIso_ = (*modifiedEcalIsoHandle)[e1];
    var_e1HcalIso_ = e1->dr03HcalDepth1TowerSumEt();
    var_e1HoE_ = e1->hadronicOverEm();
    var_e2Pt_ = e2->pt();
    var_e2scEn_ = e2->superCluster()->energy();
    var_e2scEta_ = e2->superCluster()->eta();
    var_e2Eta_ = e2->eta();
    var_e2Phi_ = e2->phi();
    var_e2DPerpIn_ = (*dPerpInHandle)[e2];
    var_e2DEtaInSeed_ = e2->deltaEtaSeedClusterTrackAtVtx();
    var_e2MvaScore_ = e2->userFloat("mvaMergedElectronValues");
    var_e2MvaCategory_ = e2->userInt("mvaMergedElectronCategories");
    var_e2PassMva_ = e2->electronID("mvaMergedElectron");
    var_e2TrkIso_ = (*modifiedTrkIsoHandle)[e2];
    var_e2EcalIso_ = (*modifiedEcalIsoHandle)[e2];
    var_e2HcalIso_ = e2->dr03HcalDepth1TowerSumEt();
    var_e2HoE_ = e2->hadronicOverEm();
    var_e1e2Rho_ = *rhoHandle;
    var_e1e2InvM_ = (e1->polarP4()+e2->polarP4()).M();
    var_e1e2DR_ = reco::deltaR(e1->eta(),e1->phi(),e2->eta(),e2->phi());
    var_e1e2scDR_ = reco::deltaR(e1->superCluster()->eta(),e1->superCluster()->phi(),
                                 e2->superCluster()->eta(),e2->superCluster()->phi());
  };

  var_weight_ = aWeight;

  if (isolatedHighPtMuons.size()==1) { // W+jet CR
    edm::Handle<edm::View<pat::MET>> metHandle;
    iEvent.getByToken(metToken_, metHandle);

    const auto& aMET = metHandle->at(0);

    if ( aMET.pt() < 50. )
      return;

    var_1m_muPt_ = isolatedHighPtMuons.front()->tunePMuonBestTrack()->pt();
    var_1m_muEta_ = isolatedHighPtMuons.front()->tunePMuonBestTrack()->eta();
    var_1m_muPhi_ = isolatedHighPtMuons.front()->tunePMuonBestTrack()->phi();
    var_1m_met_ = aMET.pt();
    var_1m_metPhi_ = aMET.phi();

    if ( acceptEles.empty() && nonHeepEles.size()==2 ) {
      if (reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),
                        nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi()) < 0.09) {
        fillEleBranch(nonHeepEles.front(),nonHeepEles.at(1));
        tree_1m_denom_->Fill();
      }
    } else if ( acceptEles.size()==2 && nonHeepEles.empty() ) {
      if (reco::deltaR2(acceptEles.front()->eta(),acceptEles.front()->phi(),
                        acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) < 0.09) {
        fillEleBranch(acceptEles.front(),acceptEles.at(1));
        tree_1m_numer_->Fill();
      }
    }
  }

  if (isolatedHighPtMuons.size()==2) {

    var_1m_muPt_ = isolatedHighPtMuons.front()->tunePMuonBestTrack()->pt();
    var_1m_muEta_ = isolatedHighPtMuons.front()->tunePMuonBestTrack()->eta();
    var_1m_muPhi_ = isolatedHighPtMuons.front()->tunePMuonBestTrack()->phi();
    var_2m_muPt_ = isolatedHighPtMuons.at(1)->tunePMuonBestTrack()->pt();
    var_2m_muEta_ = isolatedHighPtMuons.at(1)->tunePMuonBestTrack()->eta();
    var_2m_muPhi_ = isolatedHighPtMuons.at(1)->tunePMuonBestTrack()->phi();

    if ( acceptEles.empty() && nonHeepEles.size()==2 ) {
      if (reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),
                        nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi()) < 0.09) {
        fillEleBranch(nonHeepEles.front(),nonHeepEles.at(1));
        tree_2m_denom_->Fill();
      }
    } else if ( acceptEles.size()==2 && nonHeepEles.empty() ) {
      if (reco::deltaR2(acceptEles.front()->eta(),acceptEles.front()->phi(),
                        acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) < 0.09) {
        fillEleBranch(acceptEles.front(),acceptEles.at(1));
        tree_2m_numer_->Fill();
      }
    }
  }
}

DEFINE_FWK_MODULE(ResolvedFakeEleCRanalyzer);
