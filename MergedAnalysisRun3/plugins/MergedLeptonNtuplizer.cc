// MiniAOD -> flat lepton ntuple. No BDT, candidate arbitration, mass or region cuts.
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include <vector>

#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "ZprimeTo4l/MergedAnalysisRun3/interface/NtupleSchema.h"
#include "ZprimeTo4l/MergedAnalysisRun3/interface/MuonIsolation.h"
#include "TH1D.h"
#include "TNamed.h"
#include "TTree.h"

namespace {
namespace sch = mergedRun3::schema;
template <typename Ref> bool available(const Ref& ref) { return ref.isNonnull() && ref.isAvailable(); }

// Fixed maps own stable branch addresses. Starting every row with sentinels keeps
// optional track/ID columns aligned and prevents stale values between objects.
struct Columns {
  std::map<std::string, std::vector<float>> floats;
  std::map<std::string, std::vector<int>> ints;
  template <std::size_t NF, std::size_t NI>
  void book(TTree* tree, const std::string& prefix, const char* const (&fs)[NF], const char* const (&is)[NI]) {
    for (auto name : fs) tree->Branch((prefix + name).c_str(), &floats[name]);
    for (auto name : is) tree->Branch((prefix + name).c_str(), &ints[name]);
  }
  void clear() {
    for (auto& item : floats) item.second.clear();
    for (auto& item : ints) item.second.clear();
  }
  void row() {
    for (auto& item : floats) item.second.push_back(sch::missing);
    for (auto& item : ints) item.second.push_back(-1);
  }
  void f(const std::string& key, double value) { floats.at(key).back() = value; }
  void i(const std::string& key, int value) { ints.at(key).back() = value; }
};

struct MetColumns {
  bool valid = false;
  float pt = sch::missing, phi = sch::missing, sumEt = sch::missing;
  void book(TTree* tree, const std::string& prefix) {
    tree->Branch((prefix + "valid").c_str(), &valid);
    tree->Branch((prefix + "pt").c_str(), &pt);
    tree->Branch((prefix + "phi").c_str(), &phi);
    tree->Branch((prefix + "sumEt").c_str(), &sumEt);
  }
  void read(const edm::Handle<edm::View<pat::MET>>& handle) {
    valid = handle.isValid() && !handle->empty();
    pt = phi = sumEt = sch::missing;
    if (valid) {
      pt = handle->at(0).pt(); phi = handle->at(0).phi(); sumEt = handle->at(0).sumEt();
    }
  }
};

// Preserve the complete menu and status, so all electron/muon channels can choose
// their triggers downstream. -1=not run, 0=reject, 1=accept, 2=error.
void triggerStatus(const edm::Event& event, const edm::TriggerResults& results,
                   std::vector<std::string>& names, std::vector<int>& status) {
  names.clear(); status.clear();
  const auto& menu = event.triggerNames(results);
  for (unsigned i = 0; i < results.size(); ++i) {
    names.push_back(menu.triggerName(i));
    status.push_back(results.error(i) ? 2 : (!results.wasrun(i) ? -1 : int(results.accept(i))));
  }
}
}  // namespace

class MergedLeptonNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonNtuplizer(const edm::ParameterSet& cfg);
private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  void fillMuons(const edm::View<pat::Muon>&, const reco::Vertex&, const reco::BeamSpot*);
  void fillElectrons(const edm::View<pat::Electron>&, const reco::Vertex&);

  const edm::ParameterSet config_;
  const bool isMC_;
  const double muPtMin_, muEtaMax_, eleEtaMax_, eleGapLo_, eleGapHi_;
  const unsigned minLeptons_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  const edm::EDGetTokenT<edm::View<pat::Electron>> electrons_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> vertices_;
  const edm::EDGetTokenT<reco::BeamSpot> beamSpot_;
  const edm::EDGetTokenT<double> rhoToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_, puppiMetToken_;
  const edm::EDGetTokenT<edm::TriggerResults> hltToken_, filterToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerObjects_;
  const edm::EDGetTokenT<GenEventInfoProduct> generator_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileup_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genParticles_;
  TTree* tree_ = nullptr;
  TH1D* norm_ = nullptr;
  TH1D* cutflow_ = nullptr;
  Columns mu_, ele_;
  MetColumns met_, puppiMet_;
  unsigned run_ = 0, lumi_ = 0, nPV_ = 0, nMuon_ = 0, nElectron_ = 0;
  unsigned long long event_ = 0;
  double genWeight_ = 1.;
  float rho_ = sch::missing, puTrue_ = sch::missing;
  bool beamSpotValid_ = false, puValid_ = false;
  std::vector<std::string> hltNames_, filterNames_, triggerPaths_, triggerLabels_;
  std::vector<int> hltStatus_, filterStatus_, triggerId_, triggerPathOffsets_, triggerLabelOffsets_;
  std::vector<int> triggerTypes_, triggerTypeOffsets_;
  std::vector<float> triggerPt_, triggerEta_, triggerPhi_;
  std::vector<int> genIndex_, genPdgId_, genMotherPdgId_;
  std::vector<float> genPt_, genEta_, genPhi_, genMass_;
};

MergedLeptonNtuplizer::MergedLeptonNtuplizer(const edm::ParameterSet& cfg)
    : config_(cfg), isMC_(cfg.getParameter<bool>("isMC")),
      muPtMin_(cfg.getParameter<double>("muPtMin")), muEtaMax_(cfg.getParameter<double>("muEtaMax")),
      eleEtaMax_(cfg.getParameter<double>("eleEtaMax")),
      eleGapLo_(cfg.getParameter<double>("eleGapLo")), eleGapHi_(cfg.getParameter<double>("eleGapHi")),
      minLeptons_(cfg.getParameter<unsigned>("minLeptons")),
      muons_(consumes<edm::View<pat::Muon>>(cfg.getParameter<edm::InputTag>("srcMuon"))),
      electrons_(consumes<edm::View<pat::Electron>>(cfg.getParameter<edm::InputTag>("srcElectron"))),
      vertices_(consumes<edm::View<reco::Vertex>>(cfg.getParameter<edm::InputTag>("srcPV"))),
      beamSpot_(consumes<reco::BeamSpot>(cfg.getParameter<edm::InputTag>("beamSpot"))),
      rhoToken_(consumes<double>(cfg.getParameter<edm::InputTag>("rho"))),
      metToken_(consumes<edm::View<pat::MET>>(cfg.getParameter<edm::InputTag>("srcMET"))),
      puppiMetToken_(consumes<edm::View<pat::MET>>(cfg.getParameter<edm::InputTag>("srcPuppiMET"))),
      hltToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"))),
      filterToken_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("eventFilters"))),
      triggerObjects_(consumes<edm::View<pat::TriggerObjectStandAlone>>(cfg.getParameter<edm::InputTag>("triggerObjects"))),
      generator_(consumes<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("generator"))),
      pileup_(consumes<edm::View<PileupSummaryInfo>>(cfg.getParameter<edm::InputTag>("pileup"))),
      genParticles_(consumes<edm::View<reco::GenParticle>>(cfg.getParameter<edm::InputTag>("genParticles"))) {
  usesResource("TFileService");
  if (!(muPtMin_ > 0. && muEtaMax_ > 0. && eleGapLo_ >= 0. && eleGapHi_ > eleGapLo_ && eleEtaMax_ > eleGapHi_))
    throw cms::Exception("Configuration") << "Invalid lepton acceptance parameters";
}

void MergedLeptonNtuplizer::beginJob() {
  edm::Service<TFileService> fs;
  tree_ = fs->make<TTree>("Events", "HighPt muons and acceptance electrons; merged ID deferred");
  tree_->Branch("run", &run_); tree_->Branch("lumi", &lumi_); tree_->Branch("event", &event_);
  tree_->Branch("nPV", &nPV_); tree_->Branch("rho", &rho_);
  tree_->Branch("genWeight", &genWeight_); tree_->Branch("puTrue", &puTrue_); tree_->Branch("puValid", &puValid_);
  tree_->Branch("beamSpotValid", &beamSpotValid_);
  tree_->Branch("nMuon", &nMuon_); tree_->Branch("nElectron", &nElectron_);
  mu_.book(tree_, "muon_", sch::muonFloats, sch::muonInts);
  ele_.book(tree_, "electron_", sch::electronFloats, sch::electronInts);
  met_.book(tree_, "met_"); puppiMet_.book(tree_, "puppiMet_");
  tree_->Branch("hlt_names", &hltNames_); tree_->Branch("hlt_status", &hltStatus_);
  tree_->Branch("filter_names", &filterNames_); tree_->Branch("filter_status", &filterStatus_);
  tree_->Branch("trigObj_pt", &triggerPt_); tree_->Branch("trigObj_eta", &triggerEta_);
  tree_->Branch("trigObj_phi", &triggerPhi_); tree_->Branch("trigObj_pdgId", &triggerId_);
  tree_->Branch("trigObj_types", &triggerTypes_); tree_->Branch("trigObj_typeOffsets", &triggerTypeOffsets_);
  tree_->Branch("trigObj_pathOffsets", &triggerPathOffsets_); tree_->Branch("trigObj_paths", &triggerPaths_);
  tree_->Branch("trigObj_labelOffsets", &triggerLabelOffsets_); tree_->Branch("trigObj_labels", &triggerLabels_);
  tree_->Branch("gen_index", &genIndex_); tree_->Branch("gen_pdgId", &genPdgId_);
  tree_->Branch("gen_motherPdgId", &genMotherPdgId_);
  tree_->Branch("gen_pt", &genPt_); tree_->Branch("gen_eta", &genEta_);
  tree_->Branch("gen_phi", &genPhi_); tree_->Branch("gen_mass", &genMass_);
  norm_ = fs->make<TH1D>("norm", "Pre-skim normalization", 3, 0., 3.);
  const char* normLabels[] = {"Nevents", "sumGenWeights", "sumGenWeights2"};
  for (int i = 0; i < 3; ++i) norm_->GetXaxis()->SetBinLabel(i + 1, normLabels[i]);
  cutflow_ = fs->make<TH1D>("cutflow", "Unweighted event counts", 3, 0., 3.);
  cutflow_->GetXaxis()->SetBinLabel(1, "all");
  cutflow_->GetXaxis()->SetBinLabel(2, "has PV");
  cutflow_->GetXaxis()->SetBinLabel(3, "stored");
  auto* meta = fs->make<TH1D>("meta", "Production metadata", 5, 0., 5.);
  const char* metaLabels[] = {"schemaVersion", "isMC", "muonCorrectionsApplied", "mergedIdApplied", "electronIdApplied"};
  const double metaValues[] = {double(sch::version), double(isMC_), 0., 0., 0.};
  for (int i = 0; i < 5; ++i) { meta->GetXaxis()->SetBinLabel(i + 1, metaLabels[i]); meta->SetBinContent(i + 1, metaValues[i]); }
  fs->make<TNamed>("configuration", config_.dump().c_str());
}

void MergedLeptonNtuplizer::fillMuons(const edm::View<pat::Muon>& muons,
                                   const reco::Vertex& pv, const reco::BeamSpot* bs) {
  mu_.clear();
  std::vector<mergedRun3::IsolationMuon> inputs(muons.size());
  std::vector<bool> highPt(muons.size()), trackerHighPt(muons.size()), mergedBaseline(muons.size());
  for (unsigned i = 0; i < muons.size(); ++i) {
    const auto& m = muons[i];
    const auto inner = m.innerTrack(), tune = m.tunePMuonBestTrack();
    auto& in = inputs[i];
    in.eta = m.eta(); in.phi = m.phi(); in.vz = m.vz();
    in.isoValid = m.isIsolationValid() && std::isfinite(m.isolationR03().sumPt);
    if (in.isoValid) in.rawIso = m.isolationR03().sumPt;
    in.tracksValid = available(inner) && available(tune) && tune->pt() > 0.;
    if (!in.tracksValid) continue;
    in.innerEta = inner->eta(); in.innerPhi = inner->phi(); in.innerVz = inner->vz();
    in.innerPt = inner->pt(); in.tunePt = tune->pt();
    in.innerDxyBS = bs ? inner->dxy(bs->position()) : 0.;
    highPt[i] = available(m.globalTrack()) && muon::isHighPtMuon(m, pv);
    trackerHighPt[i] = muon::isTrackerHighPtMuon(m, pv);
    const bool etaAccept = std::abs(m.eta()) <= muEtaMax_;
    // ResolvedMuCRanalyzer uses TuneP pt; MergedMuCRanalyzer uses reco pt.
    in.inResolvedPool = etaAccept && tune->pt() >= muPtMin_ && (highPt[i] || trackerHighPt[i]);
    mergedBaseline[i] = etaAccept && m.pt() >= muPtMin_ && highPt[i];
  }
  for (unsigned i = 0; i < muons.size(); ++i) {
    const auto& in = inputs[i];
    if (!in.inResolvedPool && !mergedBaseline[i]) continue;
    const auto& m = muons[i];
    mu_.row();
    mu_.i("index", i); mu_.i("charge", m.charge());
    mu_.f("pt", m.pt()); mu_.f("eta", m.eta()); mu_.f("phi", m.phi());
    mu_.f("mass", m.mass()); mu_.f("energy", m.energy()); mu_.f("vz", m.vz());
    mu_.i("isHighPt", highPt[i]); mu_.i("isTrackerHighPt", trackerHighPt[i]);
    mu_.i("passMergedBaseline", mergedBaseline[i]); mu_.i("passResolvedAcceptance", in.inResolvedPool);
    mu_.i("isGlobal", m.isGlobalMuon()); mu_.i("isTracker", m.isTrackerMuon()); mu_.i("isPF", m.isPFMuon());
    mu_.i("isLoose", muon::isLooseMuon(m)); mu_.i("isMedium", muon::isMediumMuon(m));
    mu_.i("isTight", available(m.globalTrack()) && available(m.muonBestTrack()) && muon::isTightMuon(m, pv));
    auto track = [&](const std::string& prefix, const reco::TrackRef& ref) {
      const bool valid = available(ref);
      mu_.i(prefix + "_valid", valid);
      if (valid) { mu_.f(prefix + "_pt", ref->pt()); mu_.f(prefix + "_eta", ref->eta()); mu_.f(prefix + "_phi", ref->phi()); }
    };
    track("tuneP", m.tunePMuonBestTrack()); track("inner", m.innerTrack());
    track("outer", m.outerTrack()); track("global", m.globalTrack()); track("best", m.muonBestTrack());
    const auto inner = m.innerTrack(), tune = m.tunePMuonBestTrack();
    mu_.f("tuneP_relPtError", tune->ptError() / tune->pt());
    mu_.f("inner_vz", inner->vz());
    if (bs) mu_.f("inner_dxyBS", inner->dxy(bs->position()));
    mu_.f("dxy", inner->dxy(pv.position())); mu_.f("dz", inner->dz(pv.position()));
    mu_.i("nValidPixelHits", inner->hitPattern().numberOfValidPixelHits());
    mu_.i("trkLayers", inner->hitPattern().trackerLayersWithMeasurement());
    mu_.i("pixelLayers", inner->hitPattern().pixelLayersWithMeasurement());
    mu_.f("normChi2_inner", inner->normalizedChi2()); mu_.f("validFraction", inner->validFraction());
    if (available(m.globalTrack())) {
      mu_.i("nValidMuonHits", m.globalTrack()->hitPattern().numberOfValidMuonHits());
      mu_.f("normChi2_global", m.globalTrack()->normalizedChi2());
    }
    mu_.i("numOfMatchedStations", m.numberOfMatchedStations()); mu_.i("numOfChambers", m.numberOfChambers());
    mu_.i("nMatches", m.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
    mu_.i("matchesValid", m.isMatchesValid()); mu_.i("qualityValid", m.isQualityValid());
    mu_.f("segCompatibility", muon::segmentCompatibility(m)); mu_.f("caloCompatibility", m.caloCompatibility());
    mu_.f("chi2LocalPosition", m.combinedQuality().chi2LocalPosition);
    mu_.f("trkKink", m.combinedQuality().trkKink); mu_.f("glbKink", m.combinedQuality().glbKink);
    mu_.i("energyValid", m.isEnergyValid());
    // Keep the same values as MergedMuon.cc; validity is provided separately.
    const auto& e = m.calEnergy();
    mu_.f("em", e.em); mu_.f("emS9", e.emS9); mu_.f("emS25", e.emS25); mu_.f("emMax", e.emMax);
    mu_.f("had", e.had); mu_.f("hadS9", e.hadS9); mu_.f("hadMax", e.hadMax); mu_.f("ho", e.ho); mu_.f("hoS9", e.hoS9);
    mu_.i("trackIsoValid", in.isoValid);
    if (in.isoValid) mu_.f("trackiso", in.rawIso);
    const auto studyIso = mergedRun3::modifiedIsolation(inputs, i, bs != nullptr, false);
    const auto run2Iso = mergedRun3::modifiedIsolation(inputs, i, bs != nullptr, true);
    mu_.i("modifiedIsoValid", run2Iso.valid);
    mu_.f("modifiedtrackiso", studyIso.sumPt); mu_.f("modifiedIsoRun2", run2Iso.sumPt);
    mu_.i("isoNeighborStudy", studyIso.neighborIndex); mu_.i("isoNeighborRun2", run2Iso.neighborIndex);
    if (run2Iso.valid) mu_.f("relIsoRun2", run2Iso.sumPt / tune->pt());
    mu_.i("passResolved", in.inResolvedPool && run2Iso.valid && run2Iso.sumPt / tune->pt() < 0.1);
    const bool pfValid = m.isPFMuon() && m.isPFIsolationValid() && m.pt() > 0.;
    mu_.i("pfIsoValid", pfValid);
    if (pfValid) {
      const auto& iso = m.pfIsolationR04();
      mu_.f("pfrelIso04", (iso.sumChargedHadronPt + std::max(0.f, iso.sumNeutralHadronEt + iso.sumPhotonEt - 0.5f * iso.sumPUPt)) / m.pt());
    }
  }
  nMuon_ = mu_.ints.at("index").size();
}

void MergedLeptonNtuplizer::fillElectrons(const edm::View<pat::Electron>& electrons, const reco::Vertex& pv) {
  ele_.clear();
  for (unsigned i = 0; i < electrons.size(); ++i) {
    const auto& e = electrons[i];
    if (!available(e.superCluster())) continue;
    const double etaSC = std::abs(e.superCluster()->eta());
    if (etaSC > eleEtaMax_ || (etaSC > eleGapLo_ && etaSC < eleGapHi_)) continue;
    // Run-2 pre-ID acceptance. No electron ID or ET hard cut until Run-3 HEEP is integrated.
    ele_.row(); ele_.i("index", i); ele_.i("charge", e.charge());
    ele_.f("pt", e.pt()); ele_.f("eta", e.eta()); ele_.f("phi", e.phi());
    ele_.f("mass", e.mass()); ele_.f("energy", e.energy()); ele_.f("et", e.et());
    ele_.f("etaSC", e.superCluster()->eta()); ele_.f("r9", e.full5x5_r9());
    ele_.f("hOverE", e.hadronicOverEm()); ele_.f("sigmaIetaIeta", e.full5x5_sigmaIetaIeta());
    if (available(e.superCluster()->seed()))
      ele_.f("dEtaInSeed", e.deltaEtaSuperClusterTrackAtVtx() - e.superCluster()->eta() + e.superCluster()->seed()->eta());
    ele_.f("dPhiIn", e.deltaPhiSuperClusterTrackAtVtx());
    ele_.f("e2x5", e.full5x5_e2x5Max()); ele_.f("e5x5", e.full5x5_e5x5());
    ele_.f("trackIso", e.dr03TkSumPt()); ele_.f("ecalIso", e.dr03EcalRecHitSumEt());
    ele_.f("hcalIso", e.dr03HcalTowerSumEt(1)); ele_.i("ecalDriven", e.ecalDrivenSeed());
    ele_.i("passRun2Pt", e.pt() > 20.);  // modified HEEP's MinPtCut, not ET.
    const bool trackValid = available(e.gsfTrack());
    ele_.i("gsfTrackValid", trackValid);
    if (trackValid) {
      ele_.i("missingHits", e.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS));
      ele_.f("dxy", e.gsfTrack()->dxy(pv.position())); ele_.f("dz", e.gsfTrack()->dz(pv.position()));
    }
    ele_.i("modifiedHeepRun3Valid", 0);  // passModifiedHeepRun3 remains -1 (not evaluated).
  }
  nElectron_ = ele_.ints.at("index").size();
}

void MergedLeptonNtuplizer::analyze(const edm::Event& event, const edm::EventSetup&) {
  genWeight_ = isMC_ ? event.get(generator_).weight() : 1.;
  norm_->AddBinContent(1, 1.); norm_->AddBinContent(2, genWeight_); norm_->AddBinContent(3, genWeight_ * genWeight_);
  cutflow_->Fill(0.5);
  const auto& vertices = event.get(vertices_);
  if (vertices.empty()) return;
  cutflow_->Fill(1.5);
  const auto& pv = vertices[0];  // Exactly the PV convention used in Run-2 Analysis.
  run_ = event.id().run(); lumi_ = event.luminosityBlock(); event_ = event.id().event(); nPV_ = vertices.size();
  rho_ = event.get(rhoToken_);
  const auto bs = event.getHandle(beamSpot_); beamSpotValid_ = bs.isValid();
  fillMuons(event.get(muons_), pv, beamSpotValid_ ? bs.product() : nullptr);
  fillElectrons(event.get(electrons_), pv);
  if (nMuon_ + nElectron_ < minLeptons_) return;
  met_.read(event.getHandle(metToken_)); puppiMet_.read(event.getHandle(puppiMetToken_));
  const auto& hlt = event.get(hltToken_);
  triggerStatus(event, hlt, hltNames_, hltStatus_);
  triggerStatus(event, event.get(filterToken_), filterNames_, filterStatus_);
  triggerPt_.clear(); triggerEta_.clear(); triggerPhi_.clear(); triggerId_.clear();
  triggerPaths_.clear(); triggerLabels_.clear(); triggerPathOffsets_ = {0}; triggerLabelOffsets_ = {0};
  triggerTypes_.clear(); triggerTypeOffsets_ = {0};
  for (auto obj : event.get(triggerObjects_)) {
    obj.unpackPathNames(event.triggerNames(hlt));
    obj.unpackFilterLabels(event, hlt);
    triggerPt_.push_back(obj.pt()); triggerEta_.push_back(obj.eta()); triggerPhi_.push_back(obj.phi()); triggerId_.push_back(obj.pdgId());
    for (int type : obj.triggerObjectTypes()) triggerTypes_.push_back(type);
    triggerTypeOffsets_.push_back(triggerTypes_.size());
    for (const auto& path : obj.pathNames())
      if (obj.hasPathName(path, true, true)) triggerPaths_.push_back(path);
    for (const auto& label : obj.filterLabels()) triggerLabels_.push_back(label);
    triggerPathOffsets_.push_back(triggerPaths_.size()); triggerLabelOffsets_.push_back(triggerLabels_.size());
  }
  puTrue_ = sch::missing; puValid_ = false;
  genIndex_.clear(); genPdgId_.clear(); genMotherPdgId_.clear();
  genPt_.clear(); genEta_.clear(); genPhi_.clear(); genMass_.clear();
  if (isMC_) {
    const auto pu = event.getHandle(pileup_);
    if (pu.isValid()) for (const auto& entry : *pu) if (entry.getBunchCrossing() == 0) {
      puTrue_ = entry.getTrueNumInteractions(); puValid_ = true; break;
    }
    const auto& gen = event.get(genParticles_);
    for (unsigned i = 0; i < gen.size(); ++i) {
      const auto& p = gen[i];
      if ((std::abs(p.pdgId()) != 11 && std::abs(p.pdgId()) != 13) || !p.fromHardProcessFinalState()) continue;
      genIndex_.push_back(i); genPdgId_.push_back(p.pdgId());
      genMotherPdgId_.push_back(p.numberOfMothers() ? p.mother()->pdgId() : 0);
      genPt_.push_back(p.pt()); genEta_.push_back(p.eta()); genPhi_.push_back(p.phi()); genMass_.push_back(p.mass());
    }
  }
  cutflow_->Fill(2.5); tree_->Fill();
}

DEFINE_FWK_MODULE(MergedLeptonNtuplizer);
