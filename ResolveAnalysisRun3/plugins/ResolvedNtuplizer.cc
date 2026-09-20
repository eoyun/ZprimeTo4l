// 무엇: MiniAOD → Events.root (flat ntuple). 외부 코드/condor 가 struct 없이 바로 읽게
//   모든 branch 를 std::vector<float/int> flat 배열(muon_pt, electron_pt ...)로 저장.
// 어떻게: acceptance 수준(ID 무관) muon/electron 저장. correction 전/후. modified-HEEP ID
//   study 를 위해 electron 의 표준 HEEP 입력변수 + modified ValueMap 변수 전부 저장.
//   느슨한 skim(>=minLeptons), skim 이전 정규화 TH1.
// 의존: pat::Muon/Electron, muon selectors, ModifiedHEEP ValueMaps, NtupleSchema.h.

#include <memory>
#include <vector>
#include <string>
#include <map>

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
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "ZprimeTo4l/ResolveAnalysisRun3/interface/NtupleSchema.h"

#include "TTree.h"
#include "TH1D.h"

namespace sch = raRun3::schema;

namespace {
// modified-HEEP ValueMap: (저장 branch 이름, producer instance label) 쌍.
struct ModVM { const char* branch; const char* label; };
const std::vector<ModVM> kModVMs = {
  {sch::ele::modTrkIso,       "eleTrkPtIso"},
  {sch::ele::union5x5covIeIe, "union5x5covIeIe"},
  {sch::ele::union5x5covIeIp, "union5x5covIeIp"},
  {sch::ele::union5x5covIpIp, "union5x5covIpIp"},
  {sch::ele::union5x5dEtaIn,  "union5x5dEtaIn"},
  {sch::ele::union5x5dPhiIn,  "union5x5dPhiIn"},
  {sch::ele::union5x5Energy,  "union5x5Energy"},
  {sch::ele::dEtaInSeed2nd,   "dEtaInSeed2nd"},
  {sch::ele::dPhiInSC2nd,     "dPhiInSC2nd"},
  {sch::ele::dPerpIn,         "dPerpIn"},
  {sch::ele::alphaTrack,      "alphaTrack"},
  {sch::ele::alphaCalo,       "alphaCalo"},
  {sch::ele::normDParaIn,     "normalizedDParaIn"},
};
}  // namespace

class ResolvedNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedNtuplizer(const edm::ParameterSet&);
  ~ResolvedNtuplizer() override = default;

private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void endJob() override {}
  void clearBuffers();

  // 입력
  const bool isMC_;
  const edm::EDGetTokenT<edm::View<pat::Muon>>     muonToken_;
  const edm::EDGetTokenT<edm::View<pat::Electron>> eleToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>>  pvToken_;
  const edm::EDGetTokenT<reco::BeamSpot>           bsToken_;
  const edm::EDGetTokenT<double>                   rhoToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfToken_;
  const edm::EDGetTokenT<GenEventInfoProduct>          genToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> puToken_;
  const edm::EDGetTokenT<edm::TriggerResults>         trigResToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> trigObjToken_;
  const edm::EDGetTokenT<edm::TriggerResults>         metFilterToken_;
  const std::vector<std::string> trigList_;
  const std::vector<std::string> metFilterList_;
  // modified-HEEP ValueMap 토큰 (map: branch 이름 → token)
  std::map<std::string, edm::EDGetTokenT<edm::ValueMap<float>>> modTokens_;
  edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;

  // skim / 저장 기준
  const double muPtMinStore_, muEtaMaxStore_, eleEtaMaxStore_, eleGapLo_, eleGapHi_;
  const int    minLeptonsSkim_;

  // 출력
  TTree* tree_ = nullptr;
  TH1D*  norm_ = nullptr;
  TH1D*  meta_ = nullptr;

  // event 스칼라
  unsigned int b_run_ = 0, b_lumi_ = 0, b_nPV_ = 0;
  unsigned long long b_event_ = 0;
  float b_genWeight_ = 1.f, b_puTrue_ = -1.f, b_rho_ = -1.f;
  int   b_hltFired_ = 0, b_passMET_ = 0;

  // muon
  std::vector<int>   mu_index_, mu_charge_, mu_isHighPt_, mu_isTrackerHighPt_, mu_isTrackerMuon_,
                     mu_trkLayers_, mu_pixelHits_, mu_matchedStations_;
  std::vector<float> mu_pt_, mu_ptRaw_, mu_eta_, mu_phi_, mu_recoEta_, mu_recoPhi_, mu_recoVz_,
                     mu_innerPt_, mu_innerEta_, mu_innerPhi_, mu_innerVz_, mu_innerDxyBS_,
                     mu_trackIso_, mu_relPtErr_, mu_dxy_, mu_dz_;

  // electron 기본
  std::vector<int>   ele_index_, ele_charge_, ele_passModHeep_, ele_modHeepBitmap_, ele_addGsfIdx_,
                     ele_missingHits_, ele_ecalDriven_;
  std::vector<float> ele_et_, ele_ptRaw_, ele_etaRaw_, ele_phiRaw_, ele_energyRaw_,
                     ele_pt_, ele_eta_, ele_phi_, ele_mass_, ele_energyCorr_, ele_etaSC_,
                     ele_hOverE_, ele_sigmaIeta_, ele_dEtaInSeed_, ele_dPhiIn_, ele_e2x5_, ele_e5x5_,
                     ele_dr03TkSumPt_, ele_dr03EcalIso_, ele_dr03HcalIso_, ele_dxy_, ele_dz_;
  // electron modified ValueMap (map: branch 이름 → 값 배열) + ecal iso
  std::map<std::string, std::vector<float>> ele_modBuf_;
  std::vector<float> ele_modEcalHcalIso_;

  // trigger
  std::vector<float> trig_pt_, trig_eta_, trig_phi_;
  std::vector<int>   trig_filterBits_;
};

ResolvedNtuplizer::ResolvedNtuplizer(const edm::ParameterSet& iConfig) :
  isMC_(iConfig.getParameter<bool>("isMC")),
  muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
  eleToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
  pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
  bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  addGsfToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrk"))),
  genToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
  puToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
  trigResToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
  trigObjToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
  metFilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
  trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
  metFilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
  ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("modEcalIso"))),
  muPtMinStore_(iConfig.getParameter<double>("muPtMinStore")),
  muEtaMaxStore_(iConfig.getParameter<double>("muEtaMaxStore")),
  eleEtaMaxStore_(iConfig.getParameter<double>("eleEtaMaxStore")),
  eleGapLo_(iConfig.getParameter<double>("eleGapLo")),
  eleGapHi_(iConfig.getParameter<double>("eleGapHi")),
  minLeptonsSkim_(iConfig.getParameter<int>("minLeptonsSkim")) {
  usesResource("TFileService");
  // modified-HEEP ValueMap 토큰들: producer instance(module) + label 로 구성.
  const std::string modModule = iConfig.getParameter<std::string>("modHeepModule");
  for (const auto& mv : kModVMs) {
    modTokens_[mv.branch] = consumes<edm::ValueMap<float>>(edm::InputTag(modModule, mv.label));
    ele_modBuf_[mv.branch] = {};  // 버퍼 자리 확보(주소 고정)
  }
}

void ResolvedNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Events", "resolved 4l flat ntuple");

  tree_->Branch(sch::ev::run, &b_run_);
  tree_->Branch(sch::ev::lumi, &b_lumi_);
  tree_->Branch(sch::ev::event, &b_event_);
  tree_->Branch(sch::ev::genWeight, &b_genWeight_);
  tree_->Branch(sch::ev::puTrue, &b_puTrue_);
  tree_->Branch(sch::ev::nPV, &b_nPV_);
  tree_->Branch(sch::ev::rho, &b_rho_);
  tree_->Branch(sch::ev::hltFired, &b_hltFired_);
  tree_->Branch(sch::ev::passMETfilters, &b_passMET_);

  tree_->Branch(sch::mu::index, &mu_index_);
  tree_->Branch(sch::mu::charge, &mu_charge_);
  tree_->Branch(sch::mu::pt, &mu_pt_);
  tree_->Branch(sch::mu::ptRaw, &mu_ptRaw_);
  tree_->Branch(sch::mu::eta, &mu_eta_);
  tree_->Branch(sch::mu::phi, &mu_phi_);
  tree_->Branch(sch::mu::recoEta, &mu_recoEta_);
  tree_->Branch(sch::mu::recoPhi, &mu_recoPhi_);
  tree_->Branch(sch::mu::recoVz, &mu_recoVz_);
  tree_->Branch(sch::mu::innerPt, &mu_innerPt_);
  tree_->Branch(sch::mu::innerEta, &mu_innerEta_);
  tree_->Branch(sch::mu::innerPhi, &mu_innerPhi_);
  tree_->Branch(sch::mu::innerVz, &mu_innerVz_);
  tree_->Branch(sch::mu::innerDxyBS, &mu_innerDxyBS_);
  tree_->Branch(sch::mu::isHighPt, &mu_isHighPt_);
  tree_->Branch(sch::mu::isTrackerHighPt, &mu_isTrackerHighPt_);
  tree_->Branch(sch::mu::isTrackerMuon, &mu_isTrackerMuon_);
  tree_->Branch(sch::mu::trackIso, &mu_trackIso_);
  tree_->Branch(sch::mu::relPtErr, &mu_relPtErr_);
  tree_->Branch(sch::mu::trkLayers, &mu_trkLayers_);
  tree_->Branch(sch::mu::pixelHits, &mu_pixelHits_);
  tree_->Branch(sch::mu::matchedStations, &mu_matchedStations_);
  tree_->Branch(sch::mu::dxy, &mu_dxy_);
  tree_->Branch(sch::mu::dz, &mu_dz_);

  tree_->Branch(sch::ele::index, &ele_index_);
  tree_->Branch(sch::ele::charge, &ele_charge_);
  tree_->Branch(sch::ele::et, &ele_et_);
  tree_->Branch(sch::ele::ptRaw, &ele_ptRaw_);
  tree_->Branch(sch::ele::etaRaw, &ele_etaRaw_);
  tree_->Branch(sch::ele::phiRaw, &ele_phiRaw_);
  tree_->Branch(sch::ele::energyRaw, &ele_energyRaw_);
  tree_->Branch(sch::ele::pt, &ele_pt_);
  tree_->Branch(sch::ele::eta, &ele_eta_);
  tree_->Branch(sch::ele::phi, &ele_phi_);
  tree_->Branch(sch::ele::mass, &ele_mass_);
  tree_->Branch(sch::ele::energyCorr, &ele_energyCorr_);
  tree_->Branch(sch::ele::etaSC, &ele_etaSC_);
  tree_->Branch(sch::ele::passModHeep, &ele_passModHeep_);
  tree_->Branch(sch::ele::modHeepBitmap, &ele_modHeepBitmap_);
  tree_->Branch(sch::ele::addGsfIdx, &ele_addGsfIdx_);
  tree_->Branch(sch::ele::hOverE, &ele_hOverE_);
  tree_->Branch(sch::ele::sigmaIeta, &ele_sigmaIeta_);
  tree_->Branch(sch::ele::dEtaInSeed, &ele_dEtaInSeed_);
  tree_->Branch(sch::ele::dPhiIn, &ele_dPhiIn_);
  tree_->Branch(sch::ele::e2x5, &ele_e2x5_);
  tree_->Branch(sch::ele::e5x5, &ele_e5x5_);
  tree_->Branch(sch::ele::dr03TkSumPt, &ele_dr03TkSumPt_);
  tree_->Branch(sch::ele::dr03EcalIso, &ele_dr03EcalIso_);
  tree_->Branch(sch::ele::dr03HcalIso, &ele_dr03HcalIso_);
  tree_->Branch(sch::ele::missingHits, &ele_missingHits_);
  tree_->Branch(sch::ele::dxy, &ele_dxy_);
  tree_->Branch(sch::ele::dz, &ele_dz_);
  tree_->Branch(sch::ele::ecalDriven, &ele_ecalDriven_);
  tree_->Branch(sch::ele::modEcalHcalIso, &ele_modEcalHcalIso_);
  for (const auto& mv : kModVMs)
    tree_->Branch(mv.branch, &ele_modBuf_[mv.branch]);  // map 값 주소 고정

  tree_->Branch(sch::trig::pt, &trig_pt_);
  tree_->Branch(sch::trig::eta, &trig_eta_);
  tree_->Branch(sch::trig::phi, &trig_phi_);
  tree_->Branch(sch::trig::filterBits, &trig_filterBits_);

  norm_ = fs->make<TH1D>("norm", "normalization (pre-skim)", 3, 0., 3.);  // Nevents/sumw/sumw2
  meta_ = fs->make<TH1D>("meta", "metadata", 2, 0., 2.);
  meta_->SetBinContent(1, static_cast<double>(sch::kVersion));
  meta_->SetBinContent(2, 0.);  // muonCorrApplied=false (Run-3 muon payload 확보 전)
}

void ResolvedNtuplizer::clearBuffers() {
  mu_index_.clear(); mu_charge_.clear(); mu_isHighPt_.clear(); mu_isTrackerHighPt_.clear();
  mu_isTrackerMuon_.clear(); mu_trkLayers_.clear(); mu_pixelHits_.clear(); mu_matchedStations_.clear();
  mu_pt_.clear(); mu_ptRaw_.clear(); mu_eta_.clear(); mu_phi_.clear(); mu_recoEta_.clear();
  mu_recoPhi_.clear(); mu_recoVz_.clear(); mu_innerPt_.clear(); mu_innerEta_.clear();
  mu_innerPhi_.clear(); mu_innerVz_.clear(); mu_innerDxyBS_.clear(); mu_trackIso_.clear();
  mu_relPtErr_.clear(); mu_dxy_.clear(); mu_dz_.clear();

  ele_index_.clear(); ele_charge_.clear(); ele_passModHeep_.clear(); ele_modHeepBitmap_.clear();
  ele_addGsfIdx_.clear(); ele_missingHits_.clear(); ele_ecalDriven_.clear();
  ele_et_.clear(); ele_ptRaw_.clear(); ele_etaRaw_.clear(); ele_phiRaw_.clear(); ele_energyRaw_.clear();
  ele_pt_.clear(); ele_eta_.clear(); ele_phi_.clear(); ele_mass_.clear(); ele_energyCorr_.clear();
  ele_etaSC_.clear(); ele_hOverE_.clear(); ele_sigmaIeta_.clear(); ele_dEtaInSeed_.clear();
  ele_dPhiIn_.clear(); ele_e2x5_.clear(); ele_e5x5_.clear(); ele_dr03TkSumPt_.clear();
  ele_dr03EcalIso_.clear(); ele_dr03HcalIso_.clear(); ele_dxy_.clear(); ele_dz_.clear();
  ele_modEcalHcalIso_.clear();
  for (auto& kv : ele_modBuf_) kv.second.clear();

  trig_pt_.clear(); trig_eta_.clear(); trig_phi_.clear(); trig_filterBits_.clear();
}

void ResolvedNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup&) {
  clearBuffers();

  // weight (정규화는 skim 이전)
  b_genWeight_ = 1.f; b_puTrue_ = -1.f;
  if (isMC_) {
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(genToken_, genInfo);
    b_genWeight_ = static_cast<float>(genInfo->weight());
    edm::Handle<edm::View<PileupSummaryInfo>> pu;
    iEvent.getByToken(puToken_, pu);
    for (unsigned i = 0; i < pu->size(); ++i)
      if (pu->at(i).getBunchCrossing() == 0) { b_puTrue_ = pu->at(i).getTrueNumInteractions(); break; }
  }
  const double w = b_genWeight_;
  norm_->AddBinContent(1, 1.0);
  norm_->AddBinContent(2, w);
  norm_->AddBinContent(3, w * w);

  b_run_ = iEvent.id().run();
  b_lumi_ = iEvent.id().luminosityBlock();
  b_event_ = iEvent.id().event();

  edm::Handle<double> rhoH; iEvent.getByToken(rhoToken_, rhoH);
  b_rho_ = rhoH.isValid() ? static_cast<float>(*rhoH) : -1.f;

  // HLT decision
  edm::Handle<edm::TriggerResults> trigRes; iEvent.getByToken(trigResToken_, trigRes);
  const edm::TriggerNames& trigNames = iEvent.triggerNames(*trigRes);
  b_hltFired_ = 0;
  for (unsigned i = 0; i < trigRes->size() && !b_hltFired_; ++i) {
    if (!trigRes->accept(i)) continue;
    const std::string name = trigNames.triggerName(i);
    for (const auto& want : trigList_)
      if (name.find(want.substr(0, want.find("*"))) != std::string::npos) { b_hltFired_ = 1; break; }
  }

  // MET filters
  edm::Handle<edm::TriggerResults> metRes; iEvent.getByToken(metFilterToken_, metRes);
  const edm::TriggerNames& metNames = iEvent.triggerNames(*metRes);
  unsigned nPassed = 0;
  for (unsigned i = 0; i < metRes->size(); ++i) {
    if (!metRes->accept(i)) continue;
    const std::string name = metNames.triggerName(i);
    for (const auto& f : metFilterList_) if (name.find(f) != std::string::npos) ++nPassed;
  }
  b_passMET_ = (nPassed == metFilterList_.size()) ? 1 : 0;

  // trigger objects
  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjs; iEvent.getByToken(trigObjToken_, trigObjs);
  for (unsigned i = 0; i < trigObjs->size(); ++i) {
    auto obj = trigObjs->at(i);
    obj.unpackPathNames(trigNames);
    bool matched = false;
    for (const auto& pathName : obj.pathNames()) {
      for (const auto& want : trigList_)
        if (pathName.find(want.substr(0, want.find("*"))) != std::string::npos &&
            obj.hasPathName(pathName, true, true)) { matched = true; break; }
      if (matched) break;
    }
    if (!matched) continue;
    trig_pt_.push_back(obj.pt()); trig_eta_.push_back(obj.eta());
    trig_phi_.push_back(obj.phi()); trig_filterBits_.push_back(1);
  }

  // PV / BeamSpot
  edm::Handle<edm::View<reco::Vertex>> pvHandle; iEvent.getByToken(pvToken_, pvHandle);
  if (!pvHandle.isValid() || pvHandle->empty()) return;
  const reco::Vertex& pv = pvHandle->at(0);
  b_nPV_ = pvHandle->size();
  edm::Handle<reco::BeamSpot> bsHandle; iEvent.getByToken(bsToken_, bsHandle);

  // muon (acceptance 수준, ID 무관)
  edm::Handle<edm::View<pat::Muon>> muHandle; iEvent.getByToken(muonToken_, muHandle);
  for (unsigned i = 0; i < muHandle->size(); ++i) {
    const auto& m = muHandle->at(i);
    if (m.tunePMuonBestTrack().isNull() || m.innerTrack().isNull()) continue;
    const double tunePpt = m.tunePMuonBestTrack()->pt();
    if (tunePpt < muPtMinStore_ || std::abs(m.eta()) > muEtaMaxStore_) continue;

    mu_index_.push_back((int)i); mu_charge_.push_back(m.charge());
    mu_ptRaw_.push_back(tunePpt);
    mu_pt_.push_back(tunePpt);  // 현재 raw=corr (meta muonCorrApplied=0)
    mu_eta_.push_back(m.tunePMuonBestTrack()->eta());
    mu_phi_.push_back(m.tunePMuonBestTrack()->phi());
    mu_recoEta_.push_back(m.eta()); mu_recoPhi_.push_back(m.phi()); mu_recoVz_.push_back(m.vz());
    mu_innerPt_.push_back(m.innerTrack()->pt());
    mu_innerEta_.push_back(m.innerTrack()->eta()); mu_innerPhi_.push_back(m.innerTrack()->phi());
    mu_innerVz_.push_back(m.innerTrack()->vz());
    mu_innerDxyBS_.push_back(bsHandle.isValid() ? m.innerTrack()->dxy(bsHandle->position()) : 0.f);
    mu_isHighPt_.push_back(muon::isHighPtMuon(m, pv) ? 1 : 0);
    mu_isTrackerHighPt_.push_back(muon::isTrackerHighPtMuon(m, pv) ? 1 : 0);
    mu_isTrackerMuon_.push_back(m.isTrackerMuon() ? 1 : 0);
    mu_trackIso_.push_back(m.trackIso());
    mu_relPtErr_.push_back(m.tunePMuonBestTrack()->ptError() / m.tunePMuonBestTrack()->pt());
    mu_trkLayers_.push_back(m.innerTrack()->hitPattern().trackerLayersWithMeasurement());
    mu_pixelHits_.push_back(m.innerTrack()->hitPattern().numberOfValidPixelHits());
    mu_matchedStations_.push_back(m.numberOfMatchedStations());
    mu_dxy_.push_back(m.innerTrack()->dxy(pv.position()));
    mu_dz_.push_back(m.innerTrack()->dz(pv.position()));
  }

  // electron (acceptance 수준, ID 무관) + modified-HEEP study 변수
  edm::Handle<edm::View<pat::Electron>> eleHandle; iEvent.getByToken(eleToken_, eleHandle);
  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfMap; iEvent.getByToken(addGsfToken_, addGsfMap);
  const bool haveAddGsf = addGsfMap.isValid();

  // modified ValueMap handle 들
  std::map<std::string, edm::Handle<edm::ValueMap<float>>> modH;
  for (const auto& mv : kModVMs) iEvent.getByToken(modTokens_[mv.branch], modH[mv.branch]);
  edm::Handle<edm::ValueMap<float>> ecalIsoH; iEvent.getByToken(ecalIsoToken_, ecalIsoH);

  // ValueMap 안전 읽기: 없거나 ref 미포함이면 sentinel(-999).
  auto readVM = [](const edm::Handle<edm::ValueMap<float>>& h,
                   const edm::RefToBase<pat::Electron>& ref) -> float {
    if (!h.isValid() || !h->contains(ref.id())) return sch::kMissing;
    return (*h)[ref];
  };

  // 저장할 electron index(acceptance) 먼저 수집 (addGsfIdx 계산용)
  std::vector<unsigned> storeIdx;
  for (unsigned i = 0; i < eleHandle->size(); ++i) {
    const double a = std::abs(eleHandle->at(i).superCluster()->eta());
    if (a > eleEtaMaxStore_) continue;
    if (a > eleGapLo_ && a < eleGapHi_) continue;
    storeIdx.push_back(i);
  }

  for (unsigned s = 0; s < storeIdx.size(); ++s) {
    const unsigned i = storeIdx[s];
    const auto& e = eleHandle->at(i);
    const auto eref = eleHandle->refAt(i);

    // addGsfIdx (reciprocal-GSF 상대)
    int addGsfIdx = -1;
    if (haveAddGsf) {
      const reco::GsfTrackRef addTrk = (*addGsfMap)[eref];
      if (addTrk != e.gsfTrack()) {
        for (unsigned t = 0; t < storeIdx.size(); ++t) {
          if (t == s) continue;
          if (eleHandle->at(storeIdx[t]).gsfTrack() == addTrk) { addGsfIdx = (int)storeIdx[t]; break; }
        }
      }
    }

    // correction 전/후 (electron 은 둘 다 진짜값)
    const double post = e.userFloat("ecalTrkEnergyPostCorr");
    const auto rawP4 = e.polarP4();
    const auto corrP4 = rawP4 * post / e.energy();

    ele_index_.push_back((int)i); ele_charge_.push_back(e.charge());
    ele_et_.push_back(e.et());
    ele_ptRaw_.push_back(rawP4.pt()); ele_etaRaw_.push_back(rawP4.eta()); ele_phiRaw_.push_back(rawP4.phi());
    ele_energyRaw_.push_back(e.energy());
    ele_pt_.push_back(corrP4.pt()); ele_eta_.push_back(corrP4.eta()); ele_phi_.push_back(corrP4.phi());
    ele_mass_.push_back(corrP4.M()); ele_energyCorr_.push_back(post);
    ele_etaSC_.push_back(e.superCluster()->eta());
    ele_passModHeep_.push_back(e.electronID("modifiedHeepElectronID") ? 1 : 0);
    ele_modHeepBitmap_.push_back(e.userInt("modifiedHeepElectronID"));
    ele_addGsfIdx_.push_back(addGsfIdx);

    // --- 표준 HEEP 입력 (electron 메서드) ---
    ele_hOverE_.push_back(e.hadronicOverEm());
    ele_sigmaIeta_.push_back(e.full5x5_sigmaIetaIeta());
    ele_dEtaInSeed_.push_back(e.deltaEtaSeedClusterTrackAtVtx());
    ele_dPhiIn_.push_back(e.deltaPhiSuperClusterTrackAtVtx());
    ele_e2x5_.push_back(e.full5x5_e2x5Max());
    ele_e5x5_.push_back(e.full5x5_e5x5());
    ele_dr03TkSumPt_.push_back(e.dr03TkSumPt());
    ele_dr03EcalIso_.push_back(e.dr03EcalRecHitSumEt());
    ele_dr03HcalIso_.push_back(e.dr03HcalTowerSumEt());
    ele_missingHits_.push_back(e.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
    ele_dxy_.push_back(e.gsfTrack()->dxy(pv.position()));
    ele_dz_.push_back(e.gsfTrack()->dz(pv.position()));
    ele_ecalDriven_.push_back(e.ecalDrivenSeed() ? 1 : 0);

    // --- modified 변수 (ValueMap; 없으면 sentinel) ---
    ele_modEcalHcalIso_.push_back(readVM(ecalIsoH, eref));
    for (const auto& mv : kModVMs)
      ele_modBuf_[mv.branch].push_back(readVM(modH[mv.branch], eref));
  }

  // skim
  if ((int)(mu_index_.size() + ele_index_.size()) < minLeptonsSkim_) return;
  tree_->Fill();
}

DEFINE_FWK_MODULE(ResolvedNtuplizer);
