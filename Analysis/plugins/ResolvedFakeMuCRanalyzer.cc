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

#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"

#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"

#include "correction.h"

class ResolvedFakeMuCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedFakeMuCRanalyzer(const edm::ParameterSet&);
  virtual ~ResolvedFakeMuCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool passEle32WPTight(const edm::Event& iEvent, std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>>& trigObjs);

  math::PtEtaPhiMLorentzVector lvecFromTuneP(const pat::MuonRef& aMu);

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  std::unique_ptr<correction::CorrectionSet> purwgt_;
  const std::string puname_;

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const bool emulateEle32WPTightGsf_;
  const double mumass_ = 0.1056583745;

  MuonCorrectionHelper mucorrHelper_;

  std::map<std::string,TH1*> histo1d_;

  TTree* tree_1e_denom_ = nullptr;
  TTree* tree_1e_numer_ = nullptr;
  TTree* tree_2e_denom_ = nullptr;
  TTree* tree_2e_numer_ = nullptr;
  float var_weight_ = 0.;
  float var_1e_elPt_ = -1.;
  float var_1e_elEta_ = std::numeric_limits<float>::max();
  float var_1e_elPhi_ = std::numeric_limits<float>::max();
  float var_1e_met_ = -1.;
  float var_1e_metPhi_ = std::numeric_limits<float>::max();
  float var_2e_elPt_ = -1.;
  float var_2e_elEta_ = std::numeric_limits<float>::max();
  float var_2e_elPhi_ = std::numeric_limits<float>::max();
  float var_m1Pt_ = -1.;
  float var_m1Eta_ = std::numeric_limits<float>::max();
  float var_m1Phi_ = std::numeric_limits<float>::max();
  float var_m1TrkIso_ = -1.;
  float var_m2Pt_ = -1.;
  float var_m2Eta_ = std::numeric_limits<float>::max();
  float var_m2Phi_ = std::numeric_limits<float>::max();
  float var_m2TrkIso_ = -1.;
  float var_m1m2InvM_ = -1.;
  float var_m1m2Pt_ = -1.;
  float var_m1m2DR_ = -1.;
};

ResolvedFakeMuCRanalyzer::ResolvedFakeMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
purwgt_(std::move(correction::CorrectionSet::from_file((iConfig.getParameter<edm::FileInPath>("PUrwgt")).fullPath()))),
puname_(iConfig.getParameter<std::string>("PUname")),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
emulateEle32WPTightGsf_(iConfig.getParameter<bool>("emulateEle32WPTightGsf")),
mucorrHelper_() {
  usesResource("TFileService");
}

void ResolvedFakeMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  tree_1e_denom_ = fs->make<TTree>("denomTree1e","denomTree1e");
  tree_1e_denom_->Branch("wgt",&var_weight_,"wgt/F");
  tree_1e_denom_->Branch("elPt",&var_1e_elPt_,"elPt/F");
  tree_1e_denom_->Branch("elEta",&var_1e_elEta_,"elEta/F");
  tree_1e_denom_->Branch("elPhi",&var_1e_elPhi_,"elPhi/F");
  tree_1e_denom_->Branch("met",&var_1e_met_,"met/F");
  tree_1e_denom_->Branch("metPhi",&var_1e_metPhi_,"metPhi/F");
  tree_1e_numer_ = fs->make<TTree>("numerTree1e","numerTree1e");
  tree_1e_numer_->Branch("wgt",&var_weight_,"wgt/F");
  tree_1e_numer_->Branch("elPt",&var_1e_elPt_,"elPt/F");
  tree_1e_numer_->Branch("elEta",&var_1e_elEta_,"elEta/F");
  tree_1e_numer_->Branch("elPhi",&var_1e_elPhi_,"elPhi/F");
  tree_1e_numer_->Branch("met",&var_1e_met_,"met/F");
  tree_1e_numer_->Branch("metPhi",&var_1e_metPhi_,"metPhi/F");

  auto setEleBranch = [this] (TTree* atree) {
    atree->Branch("m1Pt",&var_m1Pt_,"m1Pt/F");
    atree->Branch("m1Eta",&var_m1Eta_,"m1Eta/F");
    atree->Branch("m1Phi",&var_m1Phi_,"m1Phi/F");
    atree->Branch("m1TrkIso",&var_m1TrkIso_,"m1TrkIso/F");
    atree->Branch("m2Pt",&var_m2Pt_,"m2Pt/F");
    atree->Branch("m2Eta",&var_m2Eta_,"m2Eta/F");
    atree->Branch("m2Phi",&var_m2Phi_,"m2Phi/F");
    atree->Branch("m2TrkIso",&var_m2TrkIso_,"m2TrkIso/F");
    atree->Branch("m1m2InvM",&var_m1m2InvM_,"m1m2InvM/F");
    atree->Branch("m1m2Pt",&var_m1m2Pt_,"m1m2Pt/F");
    atree->Branch("m1m2DR",&var_m1m2DR_,"m1m2DR/F");
  };

  setEleBranch(tree_1e_denom_);
  setEleBranch(tree_1e_numer_);

  tree_2e_denom_ = fs->make<TTree>("denomTree2e","denomTree2e");
  tree_2e_denom_->Branch("wgt",&var_weight_,"wgt/F");
  tree_2e_denom_->Branch("el1Pt",&var_1e_elPt_,"el1Pt/F");
  tree_2e_denom_->Branch("el1Eta",&var_1e_elEta_,"el1Eta/F");
  tree_2e_denom_->Branch("el1Phi",&var_1e_elPhi_,"el1Phi/F");
  tree_2e_denom_->Branch("el2Pt",&var_2e_elPt_,"el2Pt/F");
  tree_2e_denom_->Branch("el2Eta",&var_2e_elEta_,"el2Eta/F");
  tree_2e_denom_->Branch("el2Phi",&var_2e_elPhi_,"el2Phi/F");
  tree_2e_numer_ = fs->make<TTree>("numerTree2e","numerTree2e");
  tree_2e_numer_->Branch("wgt",&var_weight_,"wgt/F");
  tree_2e_numer_->Branch("el1Pt",&var_1e_elPt_,"el1Pt/F");
  tree_2e_numer_->Branch("el1Eta",&var_1e_elEta_,"el1Eta/F");
  tree_2e_numer_->Branch("el1Phi",&var_1e_elPhi_,"el1Phi/F");
  tree_2e_numer_->Branch("el2Pt",&var_2e_elPt_,"el2Pt/F");
  tree_2e_numer_->Branch("el2Eta",&var_2e_elEta_,"el2Eta/F");
  tree_2e_numer_->Branch("el2Phi",&var_2e_elPhi_,"el2Phi/F");

  setEleBranch(tree_2e_denom_);
  setEleBranch(tree_2e_numer_);
}

void ResolvedFakeMuCRanalyzer::endJob() {}

bool ResolvedFakeMuCRanalyzer::passEle32WPTight(const edm::Event& iEvent,
                                                std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>>& trigObjs) {
  edm::Handle<edm::View<pat::Electron>> electronHandle;
  iEvent.getByToken(srcEle_, electronHandle);

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  // Check if there's an electron matched to a trigger object which passes the two filters
  for (const auto& aEle : *electronHandle) {
    for (unsigned iTrig = 0; iTrig < trigObjHandle->size(); iTrig++) {
      const auto& trigObj = trigObjHandle->refAt(iTrig);
      auto trigObjInst = trigObjHandle->at(iTrig); // workaround for copy
      trigObjInst.unpackFilterLabels(iEvent, *trigResultHandle);

      if (!trigObjInst.hasFilterLabel("hltEle32L1DoubleEGWPTightGsfTrackIsoFilter"))
        continue;

      if (!trigObjInst.hasFilterLabel("hltEGL1SingleEGOrFilter"))
        continue;

      if (reco::deltaR2(aEle.eta(),aEle.phi(),trigObj->eta(),trigObj->phi()) > 0.01)
        continue;

      trigObjs.push_back(trigObj);

      return true;
    }
  }

  return false;
}

math::PtEtaPhiMLorentzVector ResolvedFakeMuCRanalyzer::lvecFromTuneP(const pat::MuonRef& aMu) {
  return math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                      aMu->tunePMuonBestTrack()->eta(),
                                      aMu->tunePMuonBestTrack()->phi(),
                                      mumass_);
}

void ResolvedFakeMuCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

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

  if (emulateEle32WPTightGsf_) {
    trigObjs.clear();
    isFired = passEle32WPTight(iEvent,trigObjs);
  }

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

  if ( !pvHandle.isValid() || pvHandle->empty() )
    return;

  const reco::Vertex* primaryVertex = pvHandle->ptrAt(0).get();

  std::vector<pat::ElectronRef> mediumEles;

  for (unsigned idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( !aEle->electronID("cutBasedElectronID-Fall17-94X-V2-medium") )
      continue;

    if (aEle->pt() < 20.)
      continue;

    if ( std::abs(aEle->superCluster()->eta()) > 1.4442 && std::abs(aEle->superCluster()->eta()) < 1.566 )
      continue;

    mediumEles.push_back( aEle.castTo<pat::ElectronRef>() );
  }

  if (mediumEles.empty())
    return;

  auto sortByEt = [](const pat::ElectronRef& a, const pat::ElectronRef& b) {
    return a->et() > b->et();
  };

  std::sort(mediumEles.begin(),mediumEles.end(),sortByEt);

  bool trigMatched = false;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    const auto& leadEle = mediumEles.front();

    if ( leadEle->et() > 35. &&
         reco::deltaR2(trigObj->eta(),
                       trigObj->phi(),
                       leadEle->eta(),
                       leadEle->phi()) < 0.01 ) {
      trigMatched = true;

      break;
    }
  }

  if ( !trigMatched )
    return;

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons; // but not highPtMuon
  std::vector<pat::MuonRef> nonHighPtMuons; // pass acceptance but not ID
  std::vector<pat::MuonRef> nonHighPtMuonsVLiso; // pass acceptance but not ID
  std::map<pat::MuonRef,float> nonHighPtIsos;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->tunePMuonBestTrack()->pt() < 20. || std::abs(aMuon->eta()) > 2.4 )
      continue;

    if ( muon::isHighPtMuon(*aMuon,*primaryVertex) )
      highPtMuons.push_back(aMuon.castTo<pat::MuonRef>());
    else if ( muon::isTrackerHighPtMuon(*aMuon,*primaryVertex) )
      highPtTrackerMuons.push_back(aMuon.castTo<pat::MuonRef>());
  }

  std::vector<pat::MuonRef> isolatedHighPtMuons;
  std::vector<pat::MuonRef> isolatedHighPtTrackerMuons;
  std::vector<pat::MuonRef> boostedMuons;

  auto sortByTuneP = [](const pat::MuonRef& a, const pat::MuonRef& b) {
    return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
  };

  mucorrHelper_.modifiedIso(isolatedHighPtMuons,isolatedHighPtTrackerMuons,boostedMuons,
                            highPtMuons,highPtTrackerMuons,*beamSpotHandle);

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  // concatenate(-ish) highPtMuons & highPtTrackerMuons - order is important!
  std::vector<pat::MuonRef> allHighPtMuons(isolatedHighPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), isolatedHighPtTrackerMuons.begin(), isolatedHighPtTrackerMuons.end() );
  std::sort(allHighPtMuons.begin(),allHighPtMuons.end(),sortByTuneP);

  mucorrHelper_.nonHighPtMuonIso(nonHighPtMuons,
                                 nonHighPtMuonsVLiso,
                                 nonHighPtIsos,
                                 muonHandle,
                                 allHighPtMuons,
                                 highPtMuons,
                                 highPtTrackerMuons,
                                 *primaryVertex,
                                 *beamSpotHandle,
                                 20.);

  auto checkTrackerMuPair = [&isolatedHighPtMuons] (const pat::MuonRef& m1, const pat::MuonRef& m2, bool checkId=true) -> bool {
    if (!checkId)
      return !m1->isGlobalMuon() && !m2->isGlobalMuon();

    bool check1st = std::find(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),m1) != isolatedHighPtMuons.end();
    bool check2nd = std::find(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),m2) != isolatedHighPtMuons.end();

    return !check1st && !check2nd;
  };

  auto fillMuBranch = [&,this] (const pat::MuonRef& m1, const pat::MuonRef& m2) {
    var_m1Pt_ = m1->tunePMuonBestTrack()->pt();
    var_m1Eta_ = m1->tunePMuonBestTrack()->eta();
    var_m1Phi_ = m1->tunePMuonBestTrack()->phi();
    var_m1TrkIso_ = nonHighPtIsos.find(m1)!=nonHighPtIsos.end() ? nonHighPtIsos.at(m1) : -1.;
    var_m2Pt_ = m2->tunePMuonBestTrack()->pt();
    var_m2Eta_ = m2->tunePMuonBestTrack()->eta();
    var_m2Phi_ = m2->tunePMuonBestTrack()->phi();
    var_m2TrkIso_ = nonHighPtIsos.find(m2)!=nonHighPtIsos.end() ? nonHighPtIsos.at(m2) : -1.;

    const auto lvec1 = lvecFromTuneP(m1);
    const auto lvec2 = lvecFromTuneP(m2);
    const auto lvecll = lvec1+lvec2;
    var_m1m2InvM_ = lvecll.M();
    var_m1m2Pt_ = lvecll.Pt();
    var_m1m2DR_ = reco::deltaR(lvec1.eta(),lvec1.phi(),lvec2.eta(),lvec2.phi());
  };

  var_weight_ = aWeight;

  if (mediumEles.size()==1) { // W+jet CR
    edm::Handle<edm::View<pat::MET>> metHandle;
    iEvent.getByToken(metToken_, metHandle);

    const auto& aMET = metHandle->at(0);

    if ( aMET.pt() < 50. )
      return;

    var_1e_elPt_ = mediumEles.front()->pt();
    var_1e_elEta_ = mediumEles.front()->eta();
    var_1e_elPhi_ = mediumEles.front()->phi();
    var_1e_met_ = aMET.pt();
    var_1e_metPhi_ = aMET.phi();

    if ( allHighPtMuons.empty() && nonHighPtMuonsVLiso.size()==2 ) {
      if (reco::deltaR2(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(),
                        nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(),
                        nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->eta(),
                        nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->phi()) < 0.09) {
        if (!checkTrackerMuPair(nonHighPtMuonsVLiso.front(),nonHighPtMuonsVLiso.at(1),false)) {
          fillMuBranch(nonHighPtMuonsVLiso.front(),nonHighPtMuonsVLiso.at(1));
          tree_1e_denom_->Fill();
        }
      }
    } else if ( allHighPtMuons.size()==2 && nonHighPtMuonsVLiso.empty() ) {
      if (reco::deltaR2(allHighPtMuons.front()->tunePMuonBestTrack()->eta(),
                        allHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                        allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),
                        allHighPtMuons.at(1)->tunePMuonBestTrack()->phi()) < 0.09) {
        if (!checkTrackerMuPair(allHighPtMuons.front(),allHighPtMuons.at(1),true)) {
          fillMuBranch(allHighPtMuons.front(),allHighPtMuons.at(1));
          tree_1e_numer_->Fill();
        }
      }
    }
  }

  if (mediumEles.size()==2) {

    var_1e_elPt_ = mediumEles.front()->pt();
    var_1e_elEta_ = mediumEles.front()->eta();
    var_1e_elPhi_ = mediumEles.front()->phi();
    var_2e_elPt_ = mediumEles.at(1)->pt();
    var_2e_elEta_ = mediumEles.at(1)->eta();
    var_2e_elPhi_ = mediumEles.at(1)->phi();

    if ( allHighPtMuons.empty() && nonHighPtMuonsVLiso.size()==2 ) {
      if (reco::deltaR2(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(),
                        nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(),
                        nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->eta(),
                        nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->phi()) < 0.09) {
        if (!checkTrackerMuPair(nonHighPtMuonsVLiso.front(),nonHighPtMuonsVLiso.at(1),false)) {
          fillMuBranch(nonHighPtMuonsVLiso.front(),nonHighPtMuonsVLiso.at(1));
          tree_2e_denom_->Fill();
        }
      }
    } else if ( allHighPtMuons.size()==2 && nonHighPtMuonsVLiso.empty() ) {
      if (reco::deltaR2(allHighPtMuons.front()->tunePMuonBestTrack()->eta(),
                        allHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                        allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),
                        allHighPtMuons.at(1)->tunePMuonBestTrack()->phi()) < 0.09) {
        if (!checkTrackerMuPair(allHighPtMuons.front(),allHighPtMuons.at(1),true)) {
          fillMuBranch(allHighPtMuons.front(),allHighPtMuons.at(1));
          tree_2e_numer_->Fill();
        }
      }
    }
  }
}

DEFINE_FWK_MODULE(ResolvedFakeMuCRanalyzer);
