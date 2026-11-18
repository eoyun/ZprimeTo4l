#include <algorithm>
#include <cmath>
#include <limits>
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
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "TH1D.h"
#include "TTree.h"

// Three reconstructed objects: two resolved muons and one merged-muon candidate.
// No missing-muon/MET reconstruction. All valid assignments are retained until
// a trained merged ID and an event-level candidate arbitration are supplied.
class mergedMuRun3 : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit mergedMuRun3(const edm::ParameterSet&);
private:
  void beginJob() override;
  void analyze(const edm::Event&, const edm::EventSetup&) override;
  struct IdResult { float score; bool valid; };
  IdResult mergedId(const pat::Muon&) const { return {-999.f, false}; } // TODO: XGBoost inference
  bool passesPaths(const edm::Event&, const edm::TriggerResults&,
                   const std::vector<std::string>&, bool prefix) const;
  bool isResolved(const pat::Muon&, const reco::Vertex&) const;

  const bool isMC_, requireTrigger_, requireFilters_, requireExactlyThree_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> muons_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> vertices_;
  const edm::EDGetTokenT<GenEventInfoProduct> generator_;
  const edm::EDGetTokenT<edm::TriggerResults> triggers_, filters_;
  const std::vector<std::string> triggerPaths_, filterNames_;
  const double minPt_, maxEta_, resolvedRelIsoMax_;
  TH1D* cutflow_ = nullptr;
  TTree* events_ = nullptr;
  TTree* candidates_ = nullptr;
  unsigned run_ = 0, lumi_ = 0;
  unsigned long long event_ = 0;
  double weight_ = 1.;
  unsigned nMuons_ = 0, nCandidates_ = 0;
  unsigned mergedIndex_ = 0, resolved1Index_ = 0, resolved2Index_ = 0;
  bool bdtValid_ = false, trackisoValid_ = false;
  float bdtScore_ = -999.f;
  float mergedPt_, mergedEta_, mergedPhi_, mergedMass_, trackiso_, trackRelIso_;
  float em_, emS9_, emS25_, emMax_, had_, hadS9_, hadMax_, ho_, hoS9_;
  float resolved1Pt_, resolved1Eta_, resolved1Phi_;
  float resolved2Pt_, resolved2Eta_, resolved2Phi_;
  float resolvedMass_, resolvedDR_, visibleMass_, visiblePt_;
};

mergedMuRun3::mergedMuRun3(const edm::ParameterSet& cfg)
    : isMC_(cfg.getParameter<bool>("isMC")),
      requireTrigger_(cfg.getParameter<bool>("requireTrigger")),
      requireFilters_(cfg.getParameter<bool>("requireMETFilters")),
      requireExactlyThree_(cfg.getParameter<bool>("requireExactlyThreeMuons")),
      muons_(consumes<edm::View<pat::Muon>>(cfg.getParameter<edm::InputTag>("srcMuon"))),
      vertices_(consumes<edm::View<reco::Vertex>>(cfg.getParameter<edm::InputTag>("srcPv"))),
      generator_(consumes<GenEventInfoProduct>(cfg.getParameter<edm::InputTag>("generator"))),
      triggers_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("triggerResults"))),
      filters_(consumes<edm::TriggerResults>(cfg.getParameter<edm::InputTag>("METfilters"))),
      triggerPaths_(cfg.getParameter<std::vector<std::string>>("triggerPaths")),
      filterNames_(cfg.getParameter<std::vector<std::string>>("METfilterNames")),
      minPt_(cfg.getParameter<double>("minMuonPt")),
      maxEta_(cfg.getParameter<double>("maxMuonEta")),
      resolvedRelIsoMax_(cfg.getParameter<double>("resolvedRelIsoMax")) {
  usesResource("TFileService");
  if (requireTrigger_ && triggerPaths_.empty())
    throw cms::Exception("Configuration") << "requireTrigger needs triggerPaths";
  if (requireFilters_ && filterNames_.empty())
    throw cms::Exception("Configuration") << "requireMETFilters needs METfilterNames";
}

void mergedMuRun3::beginJob() {
  edm::Service<TFileService> fs;
  cutflow_ = fs->make<TH1D>("cutflow", "Event cutflow", 6, 0., 6.);
  const char* labels[] = {"all", "trigger", "event filters", "PV", "muon multiplicity", "has candidate"};
  for (int i = 0; i < 6; ++i) cutflow_->GetXaxis()->SetBinLabel(i + 1, labels[i]);
  cutflow_->Sumw2();
  events_ = fs->make<TTree>("events", "One row per event passing PV and multiplicity; no candidate double counting");
  candidates_ = fs->make<TTree>("candidates", "All two-resolved plus one-merged assignments; BDT pending");
  for (auto* tree : {events_, candidates_}) {
    tree->Branch("run", &run_, "run/i");
    tree->Branch("lumi", &lumi_, "lumi/i");
    tree->Branch("event", &event_, "event/l");
    tree->Branch("weight", &weight_, "weight/D");
    tree->Branch("nMuons", &nMuons_, "nMuons/i");
    tree->Branch("nCandidates", &nCandidates_, "nCandidates/i");
  }
  candidates_->Branch("mergedIndex", &mergedIndex_, "mergedIndex/i");
  candidates_->Branch("resolved1Index", &resolved1Index_, "resolved1Index/i");
  candidates_->Branch("resolved2Index", &resolved2Index_, "resolved2Index/i");
  candidates_->Branch("bdtValid", &bdtValid_, "bdtValid/O");
  candidates_->Branch("bdtScore", &bdtScore_, "bdtScore/F");
  candidates_->Branch("trackisoValid", &trackisoValid_, "trackisoValid/O");
  candidates_->Branch("mergedPt", &mergedPt_, "mergedPt/F");
  candidates_->Branch("mergedEta", &mergedEta_, "mergedEta/F");
  candidates_->Branch("mergedPhi", &mergedPhi_, "mergedPhi/F");
  candidates_->Branch("mergedMass", &mergedMass_, "mergedMass/F");
  candidates_->Branch("trackiso", &trackiso_, "trackiso/F");
  candidates_->Branch("trackRelIso", &trackRelIso_, "trackRelIso/F");
  candidates_->Branch("em", &em_, "em/F");
  candidates_->Branch("emS9", &emS9_, "emS9/F");
  candidates_->Branch("emS25", &emS25_, "emS25/F");
  candidates_->Branch("emMax", &emMax_, "emMax/F");
  candidates_->Branch("had", &had_, "had/F");
  candidates_->Branch("hadS9", &hadS9_, "hadS9/F");
  candidates_->Branch("hadMax", &hadMax_, "hadMax/F");
  candidates_->Branch("ho", &ho_, "ho/F");
  candidates_->Branch("hoS9", &hoS9_, "hoS9/F");
  candidates_->Branch("resolved1Pt", &resolved1Pt_, "resolved1Pt/F");
  candidates_->Branch("resolved1Eta", &resolved1Eta_, "resolved1Eta/F");
  candidates_->Branch("resolved1Phi", &resolved1Phi_, "resolved1Phi/F");
  candidates_->Branch("resolved2Pt", &resolved2Pt_, "resolved2Pt/F");
  candidates_->Branch("resolved2Eta", &resolved2Eta_, "resolved2Eta/F");
  candidates_->Branch("resolved2Phi", &resolved2Phi_, "resolved2Phi/F");
  candidates_->Branch("resolvedMass", &resolvedMass_, "resolvedMass/F");
  candidates_->Branch("resolvedDR", &resolvedDR_, "resolvedDR/F");
  candidates_->Branch("visibleMass", &visibleMass_, "visibleMass/F");
  candidates_->Branch("visiblePt", &visiblePt_, "visiblePt/F");
}

bool mergedMuRun3::passesPaths(const edm::Event& event, const edm::TriggerResults& results,
                               const std::vector<std::string>& paths, bool prefix) const {
  const auto& names = event.triggerNames(results);
  for (const auto& pattern : paths) {
    const auto key = prefix ? pattern.substr(0, pattern.find('*')) : pattern;
    bool accepted = false;
    for (unsigned i = 0; i < results.size(); ++i) {
      const auto& name = names.triggerName(i);
      if ((prefix ? name.compare(0, key.size(), key) == 0 : name == key) && results.accept(i))
        accepted = true;
    }
    if (prefix && accepted) return true; // HLT OR
    if (!prefix && !accepted) return false; // filter AND; absent flags fail
  }
  return !prefix;
}

bool mergedMuRun3::isResolved(const pat::Muon& mu, const reco::Vertex& pv) const {
  // Initial working selection; merged candidate never passes through this cut.
  return (muon::isHighPtMuon(mu, pv) || muon::isTrackerHighPtMuon(mu, pv)) &&
         mu.isIsolationValid() && mu.isolationR03().sumPt / mu.pt() < resolvedRelIsoMax_;
}

void mergedMuRun3::analyze(const edm::Event& event, const edm::EventSetup&) {
  run_ = event.id().run();
  lumi_ = event.luminosityBlock();
  event_ = event.id().event();
  weight_ = isMC_ ? event.get(generator_).weight() : 1.; // full generator weight, no Run-2 SF
  nCandidates_ = 0;
  cutflow_->Fill(0.5, weight_);
  if (requireTrigger_ && !passesPaths(event, event.get(triggers_), triggerPaths_, true)) return;
  cutflow_->Fill(1.5, weight_);
  if (requireFilters_ && !passesPaths(event, event.get(filters_), filterNames_, false)) return;
  cutflow_->Fill(2.5, weight_);
  const auto& vertices = event.get(vertices_);
  const auto pv = std::find_if(vertices.begin(), vertices.end(), [](const reco::Vertex& v) {
    return !v.isFake() && v.ndof() > 4 && std::abs(v.z()) < 24. && v.position().rho() < 2.;
  });
  if (pv == vertices.end()) return;
  cutflow_->Fill(3.5, weight_);
  const auto& muons = event.get(muons_);
  std::vector<unsigned> baseline;
  for (unsigned i = 0; i < muons.size(); ++i) {
    const auto& mu = muons[i];
    // No high-pT ID, isolation, additional-track, or BDT requirement on merged candidates.
    if (mu.pt() >= minPt_ && std::abs(mu.eta()) <= maxEta_ && muon::isLooseMuon(mu))
      baseline.push_back(i);
  }
  nMuons_ = baseline.size();
  if (nMuons_ < 3 || (requireExactlyThree_ && nMuons_ != 3)) return;
  cutflow_->Fill(4.5, weight_);
  struct Assignment { unsigned merged, first, second; };
  std::vector<Assignment> assignments;
  for (auto m : baseline) {
    for (unsigned a = 0; a < baseline.size(); ++a) {
      const auto i = baseline[a];
      if (i == m || !isResolved(muons[i], *pv)) continue;
      for (unsigned b = a + 1; b < baseline.size(); ++b) {
        const auto j = baseline[b];
        if (j == m || !isResolved(muons[j], *pv)) continue;
        assignments.push_back(muons[i].pt() >= muons[j].pt() ? Assignment{m,i,j} : Assignment{m,j,i});
      }
    }
  }
  nCandidates_ = assignments.size();
  events_->Fill();
  if (assignments.empty()) return;
  cutflow_->Fill(5.5, weight_);
  for (const auto& assignment : assignments) {
    mergedIndex_ = assignment.merged;
    resolved1Index_ = assignment.first;
    resolved2Index_ = assignment.second;
    const auto& merged = muons[mergedIndex_];
    const auto& first = muons[resolved1Index_];
    const auto& second = muons[resolved2Index_];
    const auto id = mergedId(merged);
    bdtScore_ = id.score;
    bdtValid_ = id.valid;
    mergedPt_ = merged.pt(); mergedEta_ = merged.eta(); mergedPhi_ = merged.phi();
    mergedMass_ = merged.mass();
    trackisoValid_ = merged.isIsolationValid();
    trackiso_ = trackisoValid_ ? merged.isolationR03().sumPt : -999.f;
    trackRelIso_ = trackisoValid_ ? trackiso_ / merged.pt() : -999.f;
    const auto& energy = merged.calEnergy();
    em_ = energy.em; emS9_ = energy.emS9; emS25_ = energy.emS25; emMax_ = energy.emMax;
    had_ = energy.had; hadS9_ = energy.hadS9; hadMax_ = energy.hadMax;
    ho_ = energy.ho; hoS9_ = energy.hoS9;
    resolved1Pt_ = first.pt(); resolved1Eta_ = first.eta(); resolved1Phi_ = first.phi();
    resolved2Pt_ = second.pt(); resolved2Eta_ = second.eta(); resolved2Phi_ = second.phi();
    const auto resolved = first.p4() + second.p4();
    const auto visible = resolved + merged.p4();
    resolvedMass_ = resolved.mass();
    resolvedDR_ = reco::deltaR(first, second);
    // Visible mass only: a single muon p4 is NOT a calibrated dimuon/merged p4.
    visibleMass_ = visible.mass(); visiblePt_ = visible.pt();
    candidates_->Fill();
  }
}

DEFINE_FWK_MODULE(mergedMuRun3);
