// 무엇: MiniAOD → Events.root ntuple 생산 (resolved 4l 후단용).
// 어떻게: acceptance 수준(ID 무관)의 muon/electron 을 NtupleSchema 이름의 가변길이 branch 로
//   저장. correction 전/후 저장(electron 은 userFloat 로 raw/corr 둘 다 진짜값; muon 은
//   raw=TuneP, corr 는 Run-3 payload 확보 전까지 raw 와 동일 + muonCorrApplied=false 플래그).
//   느슨한 skim(>=minLeptons)만 적용, 정규화(skim 이전 sumw)는 별도 TH1 에 기록.
// 의존: pat::Muon/Electron, muon::isHighPtMuon, addGsfTrk ValueMap, GenEventInfo/Pileup,
//   TriggerResults/Objects, NtupleSchema.h(branch 이름 공유 계약).
//
// 재현 근거(기존 analyzer): TuneP=tunePMuonBestTrack, ID=muon::isHighPtMuon/isTrackerHighPtMuon,
//   modified-HEEP=electronID/userInt("modifiedHeepElectronID"),
//   electron mass corr = polarP4*ecalTrkEnergyPostCorr/energy.

#include <memory>
#include <vector>
#include <string>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"   // muon::isHighPtMuon 등
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "ZprimeTo4l/ResolveAnalysisRun3/interface/NtupleSchema.h"

#include "TTree.h"
#include "TH1D.h"

namespace sch = raRun3::schema;

class ResolvedNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedNtuplizer(const edm::ParameterSet&);
  ~ResolvedNtuplizer() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}

  void clearBuffers();

  // ---- 입력 ----
  const bool isMC_;
  const edm::EDGetTokenT<edm::View<pat::Muon>>     muonToken_;
  const edm::EDGetTokenT<edm::View<pat::Electron>> eleToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>>  pvToken_;
  const edm::EDGetTokenT<reco::BeamSpot>           bsToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfToken_;
  const edm::EDGetTokenT<GenEventInfoProduct>          genToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> puToken_;
  const edm::EDGetTokenT<edm::TriggerResults>         trigResToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> trigObjToken_;
  const edm::EDGetTokenT<edm::TriggerResults>         metFilterToken_;
  const std::vector<std::string> trigList_;
  const std::vector<std::string> metFilterList_;

  // ---- skim / 저장 기준 ----
  const double muPtMinStore_;   // 저장할 muon TuneP pt 하한
  const double muEtaMaxStore_;  // 저장할 muon |eta| 상한
  const double eleEtaMaxStore_; // 저장할 electron |etaSC| 상한
  const double eleGapLo_;
  const double eleGapHi_;
  const int    minLeptonsSkim_; // (muon+electron) 저장 개수 하한

  const double mumass_ = 0.1056583745;

  // ---- 출력 ----
  TTree* tree_ = nullptr;
  TH1D*  norm_ = nullptr;   // 정규화: bin1=Nevents, bin2=sumGenWeights, bin3=sumGenWeights2 (skim 이전)
  TH1D*  meta_ = nullptr;   // metadata 플래그: bin1=schemaVersion, bin2=muonCorrApplied

  // event 스칼라
  unsigned int b_run_ = 0, b_lumi_ = 0;
  unsigned long long b_event_ = 0;
  float b_genWeight_ = 1.f, b_puTrue_ = -1.f;
  int   b_hltFired_ = 0, b_passMET_ = 0;

  // muon 배열
  std::vector<int>   mu_index_, mu_charge_, mu_isHighPt_, mu_isTrackerHighPt_, mu_isTrackerMuon_,
                     mu_trkLayers_, mu_pixelHits_, mu_matchedStations_;
  std::vector<float> mu_corrTunePpt_, mu_rawTunePpt_, mu_tunePeta_, mu_tunePphi_,
                     mu_trackIso_, mu_innerPt_, mu_recoEta_, mu_recoPhi_, mu_recoVz_,
                     mu_innerEta_, mu_innerPhi_, mu_innerVz_, mu_innerDxyBS_,
                     mu_relPtErr_, mu_dxy_, mu_dz_;

  // electron 배열
  std::vector<int>   ele_index_, ele_charge_, ele_passModHeep_, ele_modHeepBitmap_, ele_addGsfIdx_;
  std::vector<float> ele_selEt_, ele_etaSC_, ele_rawPt_, ele_rawEta_, ele_rawPhi_, ele_rawEnergy_,
                     ele_corrPt_, ele_corrEta_, ele_corrPhi_, ele_corrM_;

  // trigger object 배열
  std::vector<float> trig_pt_, trig_eta_, trig_phi_;
  std::vector<int>   trig_filterBits_;
};

ResolvedNtuplizer::ResolvedNtuplizer(const edm::ParameterSet& iConfig) :
  isMC_(iConfig.getParameter<bool>("isMC")),
  muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
  eleToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
  pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  addGsfToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrk"))),
  genToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
  puToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
  trigResToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  trigObjToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  metFilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
  trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
  metFilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
  muPtMinStore_(iConfig.getParameter<double>("muPtMinStore")),
  muEtaMaxStore_(iConfig.getParameter<double>("muEtaMaxStore")),
  eleEtaMaxStore_(iConfig.getParameter<double>("eleEtaMaxStore")),
  eleGapLo_(iConfig.getParameter<double>("eleGapLo")),
  eleGapHi_(iConfig.getParameter<double>("eleGapHi")),
  minLeptonsSkim_(iConfig.getParameter<int>("minLeptonsSkim")) {
  usesResource("TFileService");
}

void ResolvedNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Events", "resolved 4l ntuple");

  // event 스칼라
  tree_->Branch(sch::ev::run, &b_run_);
  tree_->Branch(sch::ev::lumi, &b_lumi_);
  tree_->Branch(sch::ev::event, &b_event_);
  tree_->Branch(sch::ev::genWeight, &b_genWeight_);
  tree_->Branch(sch::ev::puTrue, &b_puTrue_);
  tree_->Branch(sch::ev::hltFired, &b_hltFired_);
  tree_->Branch(sch::ev::passMETfilters, &b_passMET_);

  // muon
  tree_->Branch(sch::mu::index, &mu_index_);
  tree_->Branch(sch::mu::charge, &mu_charge_);
  tree_->Branch(sch::mu::corrTunePpt, &mu_corrTunePpt_);
  tree_->Branch(sch::mu::rawTunePpt, &mu_rawTunePpt_);
  tree_->Branch(sch::mu::tunePeta, &mu_tunePeta_);
  tree_->Branch(sch::mu::tunePphi, &mu_tunePphi_);
  tree_->Branch(sch::mu::isHighPt, &mu_isHighPt_);
  tree_->Branch(sch::mu::isTrackerHighPt, &mu_isTrackerHighPt_);
  tree_->Branch(sch::mu::isTrackerMuon, &mu_isTrackerMuon_);
  tree_->Branch(sch::mu::trackIso, &mu_trackIso_);
  tree_->Branch(sch::mu::innerPt, &mu_innerPt_);
  tree_->Branch(sch::mu::recoEta, &mu_recoEta_);
  tree_->Branch(sch::mu::recoPhi, &mu_recoPhi_);
  tree_->Branch(sch::mu::recoVz, &mu_recoVz_);
  tree_->Branch(sch::mu::innerEta, &mu_innerEta_);
  tree_->Branch(sch::mu::innerPhi, &mu_innerPhi_);
  tree_->Branch(sch::mu::innerVz, &mu_innerVz_);
  tree_->Branch(sch::mu::innerDxyBS, &mu_innerDxyBS_);
  tree_->Branch(sch::mu::relPtErr, &mu_relPtErr_);
  tree_->Branch(sch::mu::trkLayers, &mu_trkLayers_);
  tree_->Branch(sch::mu::pixelHits, &mu_pixelHits_);
  tree_->Branch(sch::mu::matchedStations, &mu_matchedStations_);
  tree_->Branch(sch::mu::dxy, &mu_dxy_);
  tree_->Branch(sch::mu::dz, &mu_dz_);

  // electron
  tree_->Branch(sch::ele::index, &ele_index_);
  tree_->Branch(sch::ele::charge, &ele_charge_);
  tree_->Branch(sch::ele::selEt, &ele_selEt_);
  tree_->Branch(sch::ele::etaSC, &ele_etaSC_);
  tree_->Branch(sch::ele::rawPt, &ele_rawPt_);
  tree_->Branch(sch::ele::rawEta, &ele_rawEta_);
  tree_->Branch(sch::ele::rawPhi, &ele_rawPhi_);
  tree_->Branch(sch::ele::rawEnergy, &ele_rawEnergy_);
  tree_->Branch(sch::ele::corrPt, &ele_corrPt_);
  tree_->Branch(sch::ele::corrEta, &ele_corrEta_);
  tree_->Branch(sch::ele::corrPhi, &ele_corrPhi_);
  tree_->Branch(sch::ele::corrM, &ele_corrM_);
  tree_->Branch(sch::ele::passModHeep, &ele_passModHeep_);
  tree_->Branch(sch::ele::modHeepBitmap, &ele_modHeepBitmap_);
  tree_->Branch(sch::ele::addGsfIdx, &ele_addGsfIdx_);

  // trigger objects
  tree_->Branch(sch::trig::pt, &trig_pt_);
  tree_->Branch(sch::trig::eta, &trig_eta_);
  tree_->Branch(sch::trig::phi, &trig_phi_);
  tree_->Branch(sch::trig::filterBits, &trig_filterBits_);

  // 정규화(skim 이전): bin1=Nevents, bin2=sumGenWeights, bin3=sumGenWeights2
  norm_ = fs->make<TH1D>("norm", "normalization (pre-skim)", 3, 0., 3.);
  // metadata 플래그: bin1=schemaVersion, bin2=muonCorrApplied(0=미적용)
  meta_ = fs->make<TH1D>("meta", "metadata", 2, 0., 2.);
  meta_->SetBinContent(1, static_cast<double>(sch::kVersion));
  meta_->SetBinContent(2, 0.);  // muonCorrApplied=false (Run-3 payload 확보 전)
}

void ResolvedNtuplizer::clearBuffers() {
  mu_index_.clear(); mu_charge_.clear(); mu_isHighPt_.clear(); mu_isTrackerHighPt_.clear();
  mu_isTrackerMuon_.clear(); mu_trkLayers_.clear(); mu_pixelHits_.clear(); mu_matchedStations_.clear();
  mu_corrTunePpt_.clear(); mu_rawTunePpt_.clear(); mu_tunePeta_.clear(); mu_tunePphi_.clear();
  mu_trackIso_.clear(); mu_innerPt_.clear(); mu_recoEta_.clear(); mu_recoPhi_.clear(); mu_recoVz_.clear();
  mu_innerEta_.clear(); mu_innerPhi_.clear(); mu_innerVz_.clear(); mu_innerDxyBS_.clear();
  mu_relPtErr_.clear(); mu_dxy_.clear(); mu_dz_.clear();
  ele_index_.clear(); ele_charge_.clear(); ele_passModHeep_.clear(); ele_modHeepBitmap_.clear();
  ele_addGsfIdx_.clear();
  ele_selEt_.clear(); ele_etaSC_.clear(); ele_rawPt_.clear(); ele_rawEta_.clear(); ele_rawPhi_.clear();
  ele_rawEnergy_.clear(); ele_corrPt_.clear(); ele_corrEta_.clear(); ele_corrPhi_.clear(); ele_corrM_.clear();
  trig_pt_.clear(); trig_eta_.clear(); trig_phi_.clear(); trig_filterBits_.clear();
}

void ResolvedNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  clearBuffers();

  // ---- weight (정규화는 skim 이전에 채운다) ----
  b_genWeight_ = 1.f;
  b_puTrue_ = -1.f;
  if (isMC_) {
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(genToken_, genInfo);
    b_genWeight_ = static_cast<float>(genInfo->weight());

    edm::Handle<edm::View<PileupSummaryInfo>> pu;
    iEvent.getByToken(puToken_, pu);
    for (unsigned i = 0; i < pu->size(); ++i) {
      if (pu->at(i).getBunchCrossing() == 0) {  // in-time PU
        b_puTrue_ = static_cast<float>(pu->at(i).getTrueNumInteractions());
        break;
      }
    }
  }
  const double w = b_genWeight_;
  norm_->AddBinContent(1, 1.0);      // Nevents
  norm_->AddBinContent(2, w);        // sumGenWeights
  norm_->AddBinContent(3, w * w);    // sumGenWeights2

  b_run_ = iEvent.id().run();
  b_lumi_ = iEvent.id().luminosityBlock();
  b_event_ = iEvent.id().event();

  // ---- HLT decision (wanted path OR) ----
  edm::Handle<edm::TriggerResults> trigRes;
  iEvent.getByToken(trigResToken_, trigRes);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigRes);
  b_hltFired_ = 0;
  for (unsigned i = 0; i < trigRes->size(); ++i) {
    if (!trigRes->accept(i)) continue;
    const std::string name = trigNames.triggerName(i);
    for (const auto& want : trigList_) {
      const std::string stem = want.substr(0, want.find("*"));
      if (name.find(stem) != std::string::npos) { b_hltFired_ = 1; break; }
    }
    if (b_hltFired_) break;
  }

  // ---- MET filters (모든 요구 filter 통과?) ----
  edm::Handle<edm::TriggerResults> metRes;
  iEvent.getByToken(metFilterToken_, metRes);
  const edm::TriggerNames& metNames = iEvent.triggerNames(*metRes);
  unsigned nPassed = 0;
  for (unsigned i = 0; i < metRes->size(); ++i) {
    if (!metRes->accept(i)) continue;
    const std::string name = metNames.triggerName(i);
    for (const auto& f : metFilterList_)
      if (name.find(f) != std::string::npos) ++nPassed;
  }
  b_passMET_ = (nPassed == metFilterList_.size()) ? 1 : 0;

  // ---- trigger objects (wanted path 에 매칭된 것 저장) ----
  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjs;
  iEvent.getByToken(trigObjToken_, trigObjs);
  for (unsigned i = 0; i < trigObjs->size(); ++i) {
    auto obj = trigObjs->at(i);      // copy (unpack 위해)
    obj.unpackPathNames(trigNames);
    bool matched = false;
    for (const auto& pathName : obj.pathNames()) {
      for (const auto& want : trigList_) {
        const std::string stem = want.substr(0, want.find("*"));
        if (pathName.find(stem) != std::string::npos && obj.hasPathName(pathName, true, true)) {
          matched = true; break;
        }
      }
      if (matched) break;
    }
    if (!matched) continue;
    trig_pt_.push_back(obj.pt());
    trig_eta_.push_back(obj.eta());
    trig_phi_.push_back(obj.phi());
    trig_filterBits_.push_back(1);
  }

  // ---- PV / BeamSpot ----
  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  if (!pvHandle.isValid() || pvHandle->empty()) return;  // PV 없으면 저장 안 함
  const reco::Vertex& pv = pvHandle->at(0);

  edm::Handle<reco::BeamSpot> bsHandle;
  iEvent.getByToken(bsToken_, bsHandle);

  // ---- muon (acceptance 수준, ID 무관 저장) ----
  edm::Handle<edm::View<pat::Muon>> muHandle;
  iEvent.getByToken(muonToken_, muHandle);
  for (unsigned i = 0; i < muHandle->size(); ++i) {
    const auto& m = muHandle->at(i);
    if (m.tunePMuonBestTrack().isNull() || m.innerTrack().isNull()) continue;
    const double tunePpt = m.tunePMuonBestTrack()->pt();
    if (tunePpt < muPtMinStore_ || std::abs(m.eta()) > muEtaMaxStore_) continue;

    mu_index_.push_back(static_cast<int>(i));
    mu_charge_.push_back(m.charge());
    mu_rawTunePpt_.push_back(tunePpt);
    mu_corrTunePpt_.push_back(tunePpt);  // Run-3 payload 전까지 raw=corr (meta muonCorrApplied=false)
    mu_tunePeta_.push_back(m.tunePMuonBestTrack()->eta());
    mu_tunePphi_.push_back(m.tunePMuonBestTrack()->phi());
    mu_isHighPt_.push_back(muon::isHighPtMuon(m, pv) ? 1 : 0);
    mu_isTrackerHighPt_.push_back(muon::isTrackerHighPtMuon(m, pv) ? 1 : 0);
    mu_isTrackerMuon_.push_back(m.isTrackerMuon() ? 1 : 0);
    mu_trackIso_.push_back(m.trackIso());
    mu_innerPt_.push_back(m.innerTrack()->pt());
    mu_recoEta_.push_back(m.eta());
    mu_recoPhi_.push_back(m.phi());
    mu_recoVz_.push_back(m.vz());
    mu_innerEta_.push_back(m.innerTrack()->eta());
    mu_innerPhi_.push_back(m.innerTrack()->phi());
    mu_innerVz_.push_back(m.innerTrack()->vz());
    mu_innerDxyBS_.push_back(bsHandle.isValid() ? m.innerTrack()->dxy(bsHandle->position()) : 0.f);
    mu_relPtErr_.push_back(m.tunePMuonBestTrack()->ptError() / m.tunePMuonBestTrack()->pt());
    mu_trkLayers_.push_back(m.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    mu_pixelHits_.push_back(m.innerTrack()->hitPattern().numberOfValidPixelHits());
    mu_matchedStations_.push_back(m.numberOfMatchedStations());
    mu_dxy_.push_back(m.innerTrack()->dxy(pv.position()));
    mu_dz_.push_back(m.innerTrack()->dz(pv.position()));
  }

  // ---- electron (acceptance 수준, ID 무관 저장) ----
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(eleToken_, eleHandle);
  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfMap;
  iEvent.getByToken(addGsfToken_, addGsfMap);
  const bool haveAddGsf = addGsfMap.isValid();

  // 먼저 저장할 electron 의 (원본 index, gsfTrack, addGsfTrk) 를 모아 addGsfIdx 계산에 쓴다.
  std::vector<unsigned> storeIdx;
  for (unsigned i = 0; i < eleHandle->size(); ++i) {
    const auto& e = eleHandle->at(i);
    const double aeta = std::abs(e.superCluster()->eta());
    if (aeta > eleEtaMaxStore_) continue;                       // |etaSC|>2.5
    if (aeta > eleGapLo_ && aeta < eleGapHi_) continue;         // EB-EE gap
    storeIdx.push_back(i);
  }

  for (unsigned s = 0; s < storeIdx.size(); ++s) {
    const unsigned i = storeIdx[s];
    const auto& e = eleHandle->at(i);
    const auto eref = eleHandle->refAt(i);

    // addGsfIdx: e 의 추가 GSF track 이 다른 저장 electron 의 primary gsfTrack 과 같으면 그 index.
    int addGsfIdx = -1;
    if (haveAddGsf) {
      const reco::GsfTrackRef addTrk = (*addGsfMap)[eref];
      if (addTrk != e.gsfTrack()) {  // 추가 GSF track 이 있음
        for (unsigned t = 0; t < storeIdx.size(); ++t) {
          if (t == s) continue;
          const auto& other = eleHandle->at(storeIdx[t]);
          if (other.gsfTrack() == addTrk) { addGsfIdx = static_cast<int>(storeIdx[t]); break; }
        }
      }
    }

    // correction 전/후: raw=polarP4/energy(), corr=ecalTrkEnergyPostCorr 적용 (기존 재현).
    const double post = e.userFloat("ecalTrkEnergyPostCorr");
    const auto rawP4 = e.polarP4();
    const auto corrP4 = rawP4 * post / e.energy();

    ele_index_.push_back(static_cast<int>(i));
    ele_charge_.push_back(e.charge());
    ele_selEt_.push_back(e.et());
    ele_etaSC_.push_back(e.superCluster()->eta());
    ele_rawPt_.push_back(rawP4.pt());
    ele_rawEta_.push_back(rawP4.eta());
    ele_rawPhi_.push_back(rawP4.phi());
    ele_rawEnergy_.push_back(e.energy());
    ele_corrPt_.push_back(corrP4.pt());
    ele_corrEta_.push_back(corrP4.eta());
    ele_corrPhi_.push_back(corrP4.phi());
    ele_corrM_.push_back(corrP4.M());
    ele_passModHeep_.push_back(e.electronID("modifiedHeepElectronID") ? 1 : 0);
    ele_modHeepBitmap_.push_back(e.userInt("modifiedHeepElectronID"));
    ele_addGsfIdx_.push_back(addGsfIdx);
  }

  // ---- 느슨한 skim: (muon+electron) 저장 개수 >= minLeptons ----
  const int nLep = static_cast<int>(mu_index_.size() + ele_index_.size());
  if (nLep < minLeptonsSkim_) return;  // 정규화는 이미 위에서 기록됨(skim 이전)

  tree_->Fill();
}

DEFINE_FWK_MODULE(ResolvedNtuplizer);
