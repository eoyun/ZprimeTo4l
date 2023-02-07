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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"
#include "ZprimeTo4l/Analysis/interface/PairingHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"

class ResolvedEleCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedEleCRanalyzer(const edm::ParameterSet&);
  virtual ~ResolvedEleCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;

  const edm::FileInPath recoSFpath_;
  const edm::FileInPath recoSFlowPtPath_;
  const edm::FileInPath FFpath_;
  const std::vector<std::string> trigList_;
  const double etThres1_;
  const double etThres2_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  std::unique_ptr<TFile> recoSFfile_;
  std::unique_ptr<TFile> recoSFlowPtFile_;
  TH2D* recoSF_;
  TH2D* recoSFlowPt_;

  std::unique_ptr<TFile> FFfile_;
  TF1* ff_EB_below_;
  TF1* ff_EE_below_;
  TF1* ff_EB_above_;
  TF1* ff_EE_above_;
  TF1* ff_PPFF_EB_below_;
  TF1* ff_PPFF_EE_below_;
  TF1* ff_PPFF_EB_above_;
  TF1* ff_PPFF_EE_above_;
};

ResolvedEleCRanalyzer::ResolvedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
recoSFpath_(iConfig.getParameter<edm::FileInPath>("recoSFpath")),
recoSFlowPtPath_(iConfig.getParameter<edm::FileInPath>("recoSFlowPtPath")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
etThres1_(iConfig.getParameter<double>("etThres1")),
etThres2_(iConfig.getParameter<double>("etThres2")) {
  usesResource("TFileService");
}

void ResolvedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  recoSFfile_ = std::make_unique<TFile>(recoSFpath_.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));

  recoSFlowPtFile_ = std::make_unique<TFile>(recoSFlowPtPath_.fullPath().c_str(),"READ");
  recoSFlowPt_ = static_cast<TH2D*>(recoSFlowPtFile_->Get("EGamma_SF2D"));

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  ff_EB_below_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_dr_below_EB"));
  ff_EE_below_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_dr_below_EE"));
  ff_EB_above_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_dr_above_EB"));
  ff_EE_above_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_dr_above_EE"));
  ff_PPFF_EB_below_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_PPFF_dr_below_EB"));
  ff_PPFF_EE_below_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_PPFF_dr_below_EE"));
  ff_PPFF_EB_above_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_PPFF_dr_above_EB"));
  ff_PPFF_EE_above_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_PPFF_dr_above_EE"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  // 3P0F
  histo1d_["3P0F_P1_Et"] = fs->make<TH1D>("3P0F_P1_Et","3P0F Pass E1 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["3P0F_P1_eta"] = fs->make<TH1D>("3P0F_P1_eta","3P0F Pass E1 #eta",200,-2.5,2.5);
  histo1d_["3P0F_P1_phi"] = fs->make<TH1D>("3P0F_P1_phi","3P0F Pass E1 #phi",128,-3.2,3.2);

  histo1d_["3P0F_P2_Et"] = fs->make<TH1D>("3P0F_P2_Et","3P0F Pass E2 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["3P0F_P2_eta"] = fs->make<TH1D>("3P0F_P2_eta","3P0F Pass E2 #eta",200,-2.5,2.5);
  histo1d_["3P0F_P2_phi"] = fs->make<TH1D>("3P0F_P2_phi","3P0F Pass E2 #phi",128,-3.2,3.2);

  histo1d_["3P0F_P3_Et"] = fs->make<TH1D>("3P0F_P3_Et","3P0F Pass E3 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["3P0F_P3_eta"] = fs->make<TH1D>("3P0F_P3_eta","3P0F Pass E3 #eta",200,-2.5,2.5);
  histo1d_["3P0F_P3_phi"] = fs->make<TH1D>("3P0F_P3_phi","3P0F Pass E3 #phi",128,-3.2,3.2);

  histo1d_["3P0F_ll_invM"] = fs->make<TH1D>("3P0F_ll_invM","3P0F P1P2 M(ll);M [GeV];",120,60.,120.);
  histo1d_["3P0F_ll_pt"] = fs->make<TH1D>("3P0F_ll_pt","3P0F P1P2 p_{T}(ll);p_{T} [GeV];",200,0.,200.);
  histo1d_["3P0F_ll_dr"] = fs->make<TH1D>("3P0F_ll_dr","3P0F P1P2 dR(ll)",128,0.,6.4);

  histo1d_["3P0F_Et_EB"] = fs->make<TH1D>("3P0F_Et_EB","3P0F Pass E_{T} (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["3P0F_Et_EE"] = fs->make<TH1D>("3P0F_Et_EE","3P0F Pass E_{T} (EE);E_{T} [GeV];",200,0.,1000.);

  histo1d_["3P0F_dr_EB"] = fs->make<TH1D>("3P0F_dr_EB","3P0F Pass min dR(ll) (EB)",128,0.,6.4);
  histo1d_["3P0F_dr_EE"] = fs->make<TH1D>("3P0F_dr_EE","3P0F Pass min dR(ll) (EE)",128,0.,6.4);

  // 2P1F
  histo1d_["2P1F_P1_Et"] = fs->make<TH1D>("2P1F_P1_Et","2P1F Pass E1 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P1F_P1_eta"] = fs->make<TH1D>("2P1F_P1_eta","2P1F Pass E1 #eta",200,-2.5,2.5);
  histo1d_["2P1F_P1_phi"] = fs->make<TH1D>("2P1F_P1_phi","2P1F Pass E1 #phi",128,-3.2,3.2);

  histo1d_["2P1F_P2_Et"] = fs->make<TH1D>("2P1F_P2_Et","2P1F Pass E2 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P1F_P2_eta"] = fs->make<TH1D>("2P1F_P2_eta","2P1F Pass E2 #eta",200,-2.5,2.5);
  histo1d_["2P1F_P2_phi"] = fs->make<TH1D>("2P1F_P2_phi","2P1F Pass E2 #phi",128,-3.2,3.2);

  histo1d_["2P1F_F1_Et"] = fs->make<TH1D>("2P1F_F1_Et","2P1F Fail E1 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P1F_F1_eta"] = fs->make<TH1D>("2P1F_F1_eta","2P1F Fail E1 #eta",200,-2.5,2.5);
  histo1d_["2P1F_F1_phi"] = fs->make<TH1D>("2P1F_F1_phi","2P1F Fail E1 #phi",128,-3.2,3.2);

  histo1d_["2P1F_ll_invM"] = fs->make<TH1D>("2P1F_ll_invM","2P1F P1P2 M(ll);M [GeV];",120,60.,120.);
  histo1d_["2P1F_ll_pt"] = fs->make<TH1D>("2P1F_ll_pt","2P1F P1P2 p_{T}(ll);p_{T} [GeV];",200,0.,200.);
  histo1d_["2P1F_ll_dr"] = fs->make<TH1D>("2P1F_ll_dr","2P1F P1P2 dR(ll)",128,0.,6.4);

  histo1d_["2P1F_Et_EB"] = fs->make<TH1D>("2P1F_Et_EB","2P1F Fail E_{T} (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P1F_Et_EE"] = fs->make<TH1D>("2P1F_Et_EE","2P1F Fail E_{T} (EE);E_{T} [GeV];",200,0.,1000.);

  histo1d_["2P1F_dr_EB"] = fs->make<TH1D>("2P1F_dr_EB","2P1F Fail min dR(ll) (EB)",128,0.,6.4);
  histo1d_["2P1F_dr_EE"] = fs->make<TH1D>("2P1F_dr_EE","2P1F Fail min dR(ll) (EE)",128,0.,6.4);

  // 2P2F
  histo1d_["2P2F_P1_Et"] = fs->make<TH1D>("2P2F_P1_Et","2P2F Pass E1 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_P1_eta"] = fs->make<TH1D>("2P2F_P1_eta","2P2F Pass E1 #eta",200,-2.5,2.5);
  histo1d_["2P2F_P1_phi"] = fs->make<TH1D>("2P2F_P1_phi","2P2F Pass E1 #phi",128,-3.2,3.2);

  histo1d_["2P2F_P2_Et"] = fs->make<TH1D>("2P2F_P2_Et","2P2F Pass E2 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_P2_eta"] = fs->make<TH1D>("2P2F_P2_eta","2P2F Pass E2 #eta",200,-2.5,2.5);
  histo1d_["2P2F_P2_phi"] = fs->make<TH1D>("2P2F_P2_phi","2P2F Pass E2 #phi",128,-3.2,3.2);

  histo1d_["2P2F_F1_Et"] = fs->make<TH1D>("2P2F_F1_Et","2P2F Fail E1 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_F1_eta"] = fs->make<TH1D>("2P2F_F1_eta","2P2F Fail E1 #eta",200,-2.5,2.5);
  histo1d_["2P2F_F1_phi"] = fs->make<TH1D>("2P2F_F1_phi","2P2F Fail E1 #phi",128,-3.2,3.2);

  histo1d_["2P2F_F2_Et"] = fs->make<TH1D>("2P2F_F2_Et","2P2F Fail E2 E_{T};E_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_F2_eta"] = fs->make<TH1D>("2P2F_F2_eta","2P2F Fail E2 #eta",200,-2.5,2.5);
  histo1d_["2P2F_F2_phi"] = fs->make<TH1D>("2P2F_F2_phi","2P2F Fail E2 #phi",128,-3.2,3.2);

  histo1d_["2P2F_llll_invM"] = fs->make<TH1D>("2P2F_llll_invM","2P2F M(4l);M [GeV];",100,0.,1000.);
  histo1d_["2P2F_llll_pt"] = fs->make<TH1D>("2P2F_llll_pt","2P2F p_{T}(4l);p_{T} [GeV];",100,0.,1000.);
  histo1d_["2P2F_ll1ll2_dr"] = fs->make<TH1D>("2P2F_ll1ll2_dr","2P2F dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_ll1_invM"] = fs->make<TH1D>("2P2F_ll1_invM","2P2F M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_ll1_pt"] = fs->make<TH1D>("2P2F_ll1_pt","2P2F p_{T}(ll1);p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_ll1_dr"] = fs->make<TH1D>("2P2F_ll1_dr","2P2F dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_ll2_invM"] = fs->make<TH1D>("2P2F_ll2_invM","2P2F M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_ll2_pt"] = fs->make<TH1D>("2P2F_ll2_pt","2P2F p_{T}(ll2);p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_ll2_dr"] = fs->make<TH1D>("2P2F_ll2_dr","2P2F dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_Et_P1"] = fs->make<TH1D>("2P2F_CR_Et_P1","2P2F CR Pass E1 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_CR_Et_P2"] = fs->make<TH1D>("2P2F_CR_Et_P2","2P2F CR Pass E2 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_CR_Et_F1"] = fs->make<TH1D>("2P2F_CR_Et_F1","2P2F CR Fail E1 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_CR_Et_F2"] = fs->make<TH1D>("2P2F_CR_Et_F2","2P2F CR Fail E2 E_{T};E_{T} [GeV];",100,0.,500.);

  histo1d_["2P2F_CR_PPFF_Et_EB"] = fs->make<TH1D>("2P2F_CR_PPFF_Et_EB","2P2F CR PPFF Fail E_{T} (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_PPFF_Et_EE"] = fs->make<TH1D>("2P2F_CR_PPFF_Et_EE","2P2F CR PPFF Fail E_{T} (EE);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_PPFF_dr_EB"] = fs->make<TH1D>("2P2F_CR_PPFF_dr_EB","2P2F CR PPFF Fail min dR(ll) (EB)",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_dr_EE"] = fs->make<TH1D>("2P2F_CR_PPFF_dr_EE","2P2F CR PPFF Fail min dR(ll) (EE)",128,0.,6.4);

  histo1d_["2P2F_CR_PFPF_Et_EB_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_Et_EB_xFF","2P2F CR PFPF Fail E_{T} x Fake Factor (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_PFPF_Et_EE_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_Et_EE_xFF","2P2F CR PFPF Fail E_{T} x Fake Factor (EE);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_PFPF_dr_EB_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_dr_EB_xFF","2P2F CR PFPF Fail min dR(ll) x Fake Factor (EB)",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_dr_EE_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_dr_EE_xFF","2P2F CR PFPF Fail min dR(ll) x Fake Factor (EE)",128,0.,6.4);

  histo1d_["2P2F_CR_PFPF_llll_invM"] = fs->make<TH1D>("2P2F_CR_PFPF_llll_invM","2P2F CR PFPF M(4l);M [GeV];",150,50.,200.);
  histo1d_["2P2F_CR_PFPF_llll_pt"] = fs->make<TH1D>("2P2F_CR_PFPF_llll_pt","2P2F CR PFPF p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll1ll2_dr"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1ll2_dr","2P2F CR PFPF dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_ll1_invM"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_invM","2P2F CR PFPF M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll1_pt"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_pt","2P2F CR PFPF p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll1_dr"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_dr","2P2F CR PFPF dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_ll2_invM"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_invM","2P2F CR PFPF M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll2_pt"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_pt","2P2F CR PFPF p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll2_dr"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_dr","2P2F CR PFPF dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_PPFF_llll_invM"] = fs->make<TH1D>("2P2F_CR_PPFF_llll_invM","2P2F CR PPFF M(4l);M [GeV];",150,50.,200.);
  histo1d_["2P2F_CR_PPFF_llll_pt"] = fs->make<TH1D>("2P2F_CR_PPFF_llll_pt","2P2F CR PPFF p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll1ll2_dr"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1ll2_dr","2P2F CR PPFF dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_ll1_invM"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_invM","2P2F CR PPFF M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll1_pt"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_pt","2P2F CR PPFF p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll1_dr"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_dr","2P2F CR PPFF dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_ll2_invM"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_invM","2P2F CR PPFF M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll2_pt"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_pt","2P2F CR PPFF p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll2_dr"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_dr","2P2F CR PPFF dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_PFPF_llll_invM_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_llll_invM_xFF","2P2F CR PFPF M(4l) x Fake factor;M [GeV];",150,50.,200.);
  histo1d_["2P2F_CR_PFPF_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1ll2_dr_xFF","2P2F CR PFPF dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_invM_xFF","2P2F CR PFPF M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_dr_xFF","2P2F CR PFPF dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_invM_xFF","2P2F CR PFPF M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_dr_xFF","2P2F CR PFPF dR(ll2) x Fake factor",128,0.,6.4);

  histo1d_["2P2F_CR_PFPF_llll_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_PFPF_llll_invM_xFF2","2P2F CR PFPF M(4l) x Fake factor^2;M [GeV];",150,50.,200.);
  histo1d_["2P2F_CR_PFPF_ll1ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1ll2_dr_xFF2","2P2F CR PFPF dR(ll1ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_ll1_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_invM_xFF2","2P2F CR PFPF M(ll1) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll1_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_PFPF_ll1_dr_xFF2","2P2F CR PFPF dR(ll1) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_PFPF_ll2_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_invM_xFF2","2P2F CR PFPF M(ll2) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PFPF_ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_PFPF_ll2_dr_xFF2","2P2F CR PFPF dR(ll2) x Fake factor^2",128,0.,6.4);

  histo1d_["2P2F_CR_PPFF_llll_invM_xFF"] = fs->make<TH1D>("2P2F_CR_PPFF_llll_invM_xFF","2P2F CR PPFF M(4l) x Fake factor;M [GeV];",150,50.,200.);
  histo1d_["2P2F_CR_PPFF_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1ll2_dr_xFF","2P2F CR PPFF dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_invM_xFF","2P2F CR PPFF M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_dr_xFF","2P2F CR PPFF dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_invM_xFF","2P2F CR PPFF M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_dr_xFF","2P2F CR PPFF dR(ll2) x Fake factor",128,0.,6.4);

  histo1d_["2P2F_CR_PPFF_llll_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_PPFF_llll_invM_xFF2","2P2F CR PPFF M(4l) x Fake factor^2;M [GeV];",150,50.,200.);
  histo1d_["2P2F_CR_PPFF_ll1ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1ll2_dr_xFF2","2P2F CR PPFF dR(ll1ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_ll1_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_invM_xFF2","2P2F CR PPFF M(ll1) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll1_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_PPFF_ll1_dr_xFF2","2P2F CR PPFF dR(ll1) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_PPFF_ll2_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_invM_xFF2","2P2F CR PPFF M(ll2) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_PPFF_ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_PPFF_ll2_dr_xFF2","2P2F CR PPFF dR(ll2) x Fake factor^2",128,0.,6.4);

  // 3P1F
  histo1d_["3P1F_P1_Et"] = fs->make<TH1D>("3P1F_P1_Et","3P1F Pass E1 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_P1_eta"] = fs->make<TH1D>("3P1F_P1_eta","3P1F Pass E1 #eta",100,-2.5,2.5);
  histo1d_["3P1F_P1_phi"] = fs->make<TH1D>("3P1F_P1_phi","3P1F Pass E1 #phi",64,-3.2,3.2);

  histo1d_["3P1F_P2_Et"] = fs->make<TH1D>("3P1F_P2_Et","3P1F Pass E2 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_P2_eta"] = fs->make<TH1D>("3P1F_P2_eta","3P1F Pass E2 #eta",100,-2.5,2.5);
  histo1d_["3P1F_P2_phi"] = fs->make<TH1D>("3P1F_P2_phi","3P1F Pass E2 #phi",64,-3.2,3.2);

  histo1d_["3P1F_P3_Et"] = fs->make<TH1D>("3P1F_P3_Et","3P1F Pass E3 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_P3_eta"] = fs->make<TH1D>("3P1F_P3_eta","3P1F Pass E3 #eta",100,-2.5,2.5);
  histo1d_["3P1F_P3_phi"] = fs->make<TH1D>("3P1F_P3_phi","3P1F Pass E3 #phi",64,-3.2,3.2);

  histo1d_["3P1F_F1_Et"] = fs->make<TH1D>("3P1F_F1_Et","3P1F Fail E1 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_F1_eta"] = fs->make<TH1D>("3P1F_F1_eta","3P1F Fail E1 #eta",100,-2.5,2.5);
  histo1d_["3P1F_F1_phi"] = fs->make<TH1D>("3P1F_F1_phi","3P1F Fail E1 #phi",64,-3.2,3.2);

  histo1d_["3P1F_llll_invM"] = fs->make<TH1D>("3P1F_llll_invM","3P1F M(4l);M [GeV];",100,0.,500.);
  histo1d_["3P1F_llll_pt"] = fs->make<TH1D>("3P1F_llll_pt","3P1F p_{T}(4l);p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_ll1ll2_dr"] = fs->make<TH1D>("3P1F_ll1ll2_dr","3P1F dR(ll1ll2)",64,0.,6.4);
  histo1d_["3P1F_ll1_invM"] = fs->make<TH1D>("3P1F_ll1_invM","3P1F M(ll1);M [GeV];",200,0.,200.);
  histo1d_["3P1F_ll1_pt"] = fs->make<TH1D>("3P1F_ll1_pt","3P1F p_{T}(ll1);p_{T} [GeV];",200,0.,500.);
  histo1d_["3P1F_ll1_dr"] = fs->make<TH1D>("3P1F_ll1_dr","3P1F dR(ll1)",64,0.,6.4);
  histo1d_["3P1F_ll2_invM"] = fs->make<TH1D>("3P1F_ll2_invM","3P1F M(ll2);M [GeV];",200,0.,200.);
  histo1d_["3P1F_ll2_pt"] = fs->make<TH1D>("3P1F_ll2_pt","3P1F p_{T}(ll2);p_{T} [GeV];",200,0.,500.);
  histo1d_["3P1F_ll2_dr"] = fs->make<TH1D>("3P1F_ll2_dr","3P1F dR(ll2)",64,0.,6.4);

  histo1d_["3P1F_CR_Et_P1"] = fs->make<TH1D>("3P1F_CR_Et_P1","3P1F CR Pass E1 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_CR_Et_P2"] = fs->make<TH1D>("3P1F_CR_Et_P2","3P1F CR Pass E2 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_CR_Et_P3"] = fs->make<TH1D>("3P1F_CR_Et_P3","3P1F CR Pass E3 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_CR_Et_F1"] = fs->make<TH1D>("3P1F_CR_Et_F1","3P1F CR Fail E1 E_{T};E_{T} [GeV];",100,0.,500.);

  histo1d_["3P1F_CR_Et_EB"] = fs->make<TH1D>("3P1F_CR_Et_EB","3P1F CR Pass E_{T} (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["3P1F_CR_Et_EE"] = fs->make<TH1D>("3P1F_CR_Et_EE","3P1F CR Pass E_{T} (EE);E_{T} [GeV];",200,0.,1000.);
  histo1d_["3P1F_CR_dr_EB"] = fs->make<TH1D>("3P1F_CR_dr_EB","3P1F CR Pass min dR(ll) (EB)",128,0.,6.4);
  histo1d_["3P1F_CR_dr_EE"] = fs->make<TH1D>("3P1F_CR_dr_EE","3P1F CR Pass min dR(ll) (EE)",128,0.,6.4);

  histo1d_["3P1F_CR_llll_invM"] = fs->make<TH1D>("3P1F_CR_llll_invM","3P1F CR M(4l);M [GeV];",75,50.,200.);
  histo1d_["3P1F_CR_llll_pt"] = fs->make<TH1D>("3P1F_CR_llll_pt","3P1F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1ll2_dr"] = fs->make<TH1D>("3P1F_CR_ll1ll2_dr","3P1F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["3P1F_CR_ll1_invM"] = fs->make<TH1D>("3P1F_CR_ll1_invM","3P1F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_pt"] = fs->make<TH1D>("3P1F_CR_ll1_pt","3P1F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_dr"] = fs->make<TH1D>("3P1F_CR_ll1_dr","3P1F CR dR(ll1)",64,0.,6.4);
  histo1d_["3P1F_CR_ll2_invM"] = fs->make<TH1D>("3P1F_CR_ll2_invM","3P1F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_pt"] = fs->make<TH1D>("3P1F_CR_ll2_pt","3P1F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_dr"] = fs->make<TH1D>("3P1F_CR_ll2_dr","3P1F CR dR(ll2)",64,0.,6.4);

  histo1d_["3P1F_CR_llll_invM_xFF"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF","3P1F CR M(4l) x Fake factor;M [GeV];",75,50.,200.);
  histo1d_["3P1F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll1ll2_dr_xFF","3P1F CR dR(ll1ll2) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_ll1_invM_xFF"] = fs->make<TH1D>("3P1F_CR_ll1_invM_xFF","3P1F CR M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll1_dr_xFF","3P1F CR dR(ll1) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_ll2_invM_xFF"] = fs->make<TH1D>("3P1F_CR_ll2_invM_xFF","3P1F CR M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll2_dr_xFF","3P1F CR dR(ll2) x Fake factor",64,0.,6.4);

  // 4P0F
  histo1d_["4P0F_CR_P1_eta"] = fs->make<TH1D>("4P0F_CR_P1_eta","4P0F CR Pass E1 #eta",100,-2.5,2.5);
  histo1d_["4P0F_CR_P1_phi"] = fs->make<TH1D>("4P0F_CR_P1_phi","4P0F CR Pass E1 #phi",64,-3.2,3.2);

  histo1d_["4P0F_CR_P2_eta"] = fs->make<TH1D>("4P0F_CR_P2_eta","4P0F CR Pass E2 #eta",100,-2.5,2.5);
  histo1d_["4P0F_CR_P2_phi"] = fs->make<TH1D>("4P0F_CR_P2_phi","4P0F CR Pass E2 #phi",64,-3.2,3.2);

  histo1d_["4P0F_CR_P3_eta"] = fs->make<TH1D>("4P0F_CR_P3_eta","4P0F CR Pass E3 #eta",100,-2.5,2.5);
  histo1d_["4P0F_CR_P3_phi"] = fs->make<TH1D>("4P0F_CR_P3_phi","4P0F CR Pass E3 #phi",64,-3.2,3.2);

  histo1d_["4P0F_CR_P4_eta"] = fs->make<TH1D>("4P0F_CR_P4_eta","4P0F CR Pass E4 #eta",100,-2.5,2.5);
  histo1d_["4P0F_CR_P4_phi"] = fs->make<TH1D>("4P0F_CR_P4_phi","4P0F CR Pass E4 #phi",64,-3.2,3.2);

  histo1d_["4P0F_CR_Et_P1"] = fs->make<TH1D>("4P0F_CR_Et_P1","4P0F CR Pass E1 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_Et_P2"] = fs->make<TH1D>("4P0F_CR_Et_P2","4P0F CR Pass E2 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_Et_P3"] = fs->make<TH1D>("4P0F_CR_Et_P3","4P0F CR Pass E3 E_{T};E_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_Et_P4"] = fs->make<TH1D>("4P0F_CR_Et_P4","4P0F CR Pass E4 E_{T};E_{T} [GeV];",100,0.,500.);

  histo1d_["4P0F_CR_llll_invM"] = fs->make<TH1D>("4P0F_CR_llll_invM","4P0F CR M(4l);M [GeV];",75,50.,200.);
  histo1d_["4P0F_CR_llll_pt"] = fs->make<TH1D>("4P0F_CR_llll_pt","4P0F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll1ll2_dr","4P0F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["4P0F_CR_ll1_invM"] = fs->make<TH1D>("4P0F_CR_ll1_invM","4P0F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_pt"] = fs->make<TH1D>("4P0F_CR_ll1_pt","4P0F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_dr"] = fs->make<TH1D>("4P0F_CR_ll1_dr","4P0F CR dR(ll1)",64,0.,6.4);
  histo1d_["4P0F_CR_ll2_invM"] = fs->make<TH1D>("4P0F_CR_ll2_invM","4P0F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_pt"] = fs->make<TH1D>("4P0F_CR_ll2_pt","4P0F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll2_dr","4P0F CR dR(ll2)",64,0.,6.4);
}

void ResolvedEleCRanalyzer::endJob() {
  recoSFfile_->Close();
  recoSFlowPtFile_->Close();
}

void ResolvedEleCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  double aWeight = 1.;

  if (isMC_) {
    edm::Handle<double> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    double prefiringweight = *theprefweight;

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    aWeight = prefiringweight*mcweight/std::abs(mcweight);
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

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->ptrAt(0);
  else
    return;

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

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
    if ( std::abs(aEle->superCluster()->eta()) > 1.44 && std::abs(aEle->superCluster()->eta()) < 1.57 )
      continue;

    bool selected = false;

    if ( aEle->electronID("modifiedHeepElectronID") ) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      acceptEles.push_back(castEle);
      selected = true;
    }

    if (!selected && aEle->et() > etThres2_) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      nonHeepEles.push_back(castEle);
    }
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByEt);
  std::sort(nonHeepEles.begin(),nonHeepEles.end(),sortByEt);

  bool accepted = false;

  if ( acceptEles.size() > 1 ) {
    if ( acceptEles.at(1)->et() > etThres1_ )
      accepted = true;
  }

  if (!accepted)
    return;

  bool trigMatched1 = false;
  bool trigMatched2 = false;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()) < 0.09 )
      trigMatched1 = true;

    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) < 0.09 )
      trigMatched2 = true;
  }

  if ( !trigMatched1 || !trigMatched2 )
    return;

  std::vector<pat::ElectronRef> mergedEls;
  for (unsigned int idx = 0; idx < acceptEles.size(); ++idx) {
    auto aEle = acceptEles.at(idx);

    // check whether there are ME
    if ( aEle->electronID("mvaMergedElectron") )
      mergedEls.push_back(aEle);
  }

  // veto ME
  if (mergedEls.size()>0)
    return;

  // recoSF
  if (isMC_) {
    for (const auto& aEle : acceptEles) {
      TH2D* aSF = aEle->pt() >= 20. ? recoSF_ : recoSFlowPt_;
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= aSF->GetBinContent( aSF->FindBin(aEle->superCluster()->eta(), apt ) );
    }

    for (const auto& aEle : nonHeepEles) {
      TH2D* aSF = aEle->pt() >= 20. ? recoSF_ : recoSFlowPt_;
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= aSF->GetBinContent( aSF->FindBin(aEle->superCluster()->eta(), apt ) );
    }
  }

  switch (acceptEles.size()) {
    case 4:
      // 4P0F
      if ( nonHeepEles.size()==0 ) {
        std::vector<pat::ElectronRef> allEles(acceptEles);
        std::pair<pat::ElectronRef,pat::ElectronRef> pair1, pair2;

        bool paired = PairingHelper::pair4E(allEles,addGsfTrkMap,pair1,pair2);

        if (paired) {
          const auto lvecE1 = pair1.first->polarP4()*pair1.first->userFloat("ecalTrkEnergyPostCorr")/pair1.first->energy();
          const auto lvecE2 = pair1.second->polarP4()*pair1.second->userFloat("ecalTrkEnergyPostCorr")/pair1.second->energy();
          const auto lvecE3 = pair2.first->polarP4()*pair2.first->userFloat("ecalTrkEnergyPostCorr")/pair2.first->energy();
          const auto lvecE4 = pair2.second->polarP4()*pair2.second->userFloat("ecalTrkEnergyPostCorr")/pair2.second->energy();

          const auto lvecA1 = lvecE1 + lvecE2;
          const auto lvecA2 = lvecE3 + lvecE4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
          const double dr2ll2 = reco::deltaR2(lvecE3.eta(),lvecE3.phi(),lvecE4.eta(),lvecE4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();

          if ( m4l > 50. && m4l < 200. ) {
            histo1d_["4P0F_CR_P1_eta"]->Fill(acceptEles.front()->eta(), aWeight);
            histo1d_["4P0F_CR_P1_phi"]->Fill(acceptEles.front()->phi(), aWeight);
            histo1d_["4P0F_CR_P2_eta"]->Fill(acceptEles.at(1)->eta(), aWeight);
            histo1d_["4P0F_CR_P2_phi"]->Fill(acceptEles.at(1)->phi(), aWeight);
            histo1d_["4P0F_CR_P3_eta"]->Fill(acceptEles.at(2)->eta(), aWeight);
            histo1d_["4P0F_CR_P3_phi"]->Fill(acceptEles.at(2)->phi(), aWeight);
            histo1d_["4P0F_CR_P4_eta"]->Fill(acceptEles.at(3)->eta(), aWeight);
            histo1d_["4P0F_CR_P4_phi"]->Fill(acceptEles.at(3)->phi(), aWeight);

            histo1d_["4P0F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["4P0F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
            histo1d_["4P0F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
            histo1d_["4P0F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
            histo1d_["4P0F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
            histo1d_["4P0F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
            histo1d_["4P0F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
            histo1d_["4P0F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
            histo1d_["4P0F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

            histo1d_["4P0F_CR_Et_P1"]->Fill(acceptEles.front()->et(), aWeight);
            histo1d_["4P0F_CR_Et_P2"]->Fill(acceptEles.at(1)->et(), aWeight);
            histo1d_["4P0F_CR_Et_P3"]->Fill(acceptEles.at(2)->et(), aWeight);
            histo1d_["4P0F_CR_Et_P4"]->Fill(acceptEles.at(3)->et(), aWeight);
          } // CR
        } // paired
      } // 4P0F

      break;
    case 3:
      // 3P0F - FF numerator
      if ( nonHeepEles.size()==0 ) {
        if ( acceptEles.at(2)->et() < etThres1_ ) {
          auto& first = acceptEles.front();
          auto& second = acceptEles.at(1);
          auto& probe = acceptEles.at(2);

          const auto lvecE1 = first->polarP4()*first->userFloat("ecalTrkEnergyPostCorr")/first->energy();
          const auto lvecE2 = second->polarP4()*second->userFloat("ecalTrkEnergyPostCorr")/second->energy();
          const auto lvecll = lvecE1 + lvecE2;
          const double mll = lvecll.M();
          const double dr2ll = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());

          if ( mll > 84.19 && mll < 98.19 && first->charge()*second->charge() < 0 ) {
            histo1d_["3P0F_P1_Et"]->Fill(first->et(), aWeight);
            histo1d_["3P0F_P1_eta"]->Fill(first->eta(), aWeight);
            histo1d_["3P0F_P1_phi"]->Fill(first->phi(), aWeight);
            histo1d_["3P0F_P2_Et"]->Fill(second->et(), aWeight);
            histo1d_["3P0F_P2_eta"]->Fill(second->eta(), aWeight);
            histo1d_["3P0F_P2_phi"]->Fill(second->phi(), aWeight);

            histo1d_["3P0F_P3_Et"]->Fill(probe->et(), aWeight);
            histo1d_["3P0F_P3_eta"]->Fill(probe->eta(), aWeight);
            histo1d_["3P0F_P3_phi"]->Fill(probe->phi(), aWeight);

            histo1d_["3P0F_ll_invM"]->Fill(lvecll.M(), aWeight);
            histo1d_["3P0F_ll_pt"]->Fill(lvecll.pt(), aWeight);
            histo1d_["3P0F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

            const double mindr2 = std::min( reco::deltaR2(probe->eta(),probe->phi(),first->eta(),first->phi()),
                                            reco::deltaR2(probe->eta(),probe->phi(),second->eta(),second->phi()) );

            if ( std::abs( probe->superCluster()->eta() ) < 1.5 ) {
              histo1d_["3P0F_Et_EB"]->Fill(probe->et(), aWeight);
              histo1d_["3P0F_dr_EB"]->Fill( std::sqrt(mindr2), aWeight );
            } else {
              histo1d_["3P0F_Et_EE"]->Fill(probe->et(), aWeight);
              histo1d_["3P0F_dr_EE"]->Fill( std::sqrt(mindr2), aWeight );
            }
          } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
        } else {
          for (unsigned idx = 0; idx < acceptEles.size(); idx++) {
            const auto& first = acceptEles.at(idx);
            const auto& second = acceptEles.at( (idx+1)%acceptEles.size() );
            const auto& probe = acceptEles.at( (idx+2)%acceptEles.size() );

            const auto lvecE1 = first->polarP4()*first->userFloat("ecalTrkEnergyPostCorr")/first->energy();
            const auto lvecE2 = second->polarP4()*second->userFloat("ecalTrkEnergyPostCorr")/second->energy();
            const auto lvecll = lvecE1 + lvecE2;
            const double mll = lvecll.M();
            const double dr2ll = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());

            if ( mll > 84.19 && mll < 98.19 && first->charge()*second->charge() < 0 ) {
              histo1d_["3P0F_P1_Et"]->Fill(first->et(), aWeight);
              histo1d_["3P0F_P1_eta"]->Fill(first->eta(), aWeight);
              histo1d_["3P0F_P1_phi"]->Fill(first->phi(), aWeight);
              histo1d_["3P0F_P2_Et"]->Fill(second->et(), aWeight);
              histo1d_["3P0F_P2_eta"]->Fill(second->eta(), aWeight);
              histo1d_["3P0F_P2_phi"]->Fill(second->phi(), aWeight);

              histo1d_["3P0F_P3_Et"]->Fill(probe->et(), aWeight);
              histo1d_["3P0F_P3_eta"]->Fill(probe->eta(), aWeight);
              histo1d_["3P0F_P3_phi"]->Fill(probe->phi(), aWeight);

              histo1d_["3P0F_ll_invM"]->Fill(lvecll.M(), aWeight);
              histo1d_["3P0F_ll_pt"]->Fill(lvecll.pt(), aWeight);
              histo1d_["3P0F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

              const double mindr2 = std::min( reco::deltaR2(probe->eta(),probe->phi(),first->eta(),first->phi()),
                                              reco::deltaR2(probe->eta(),probe->phi(),second->eta(),second->phi()) );

              if ( std::abs( probe->superCluster()->eta() ) < 1.5 ) {
                histo1d_["3P0F_Et_EB"]->Fill(probe->et(), aWeight);
                histo1d_["3P0F_dr_EB"]->Fill( std::sqrt(mindr2), aWeight );
              } else {
                histo1d_["3P0F_Et_EE"]->Fill(probe->et(), aWeight);
                histo1d_["3P0F_dr_EE"]->Fill( std::sqrt(mindr2), aWeight );
              }
            } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
          } // Z candidate combination
        }
      } // 3P0F

      // 3P1F
      if ( nonHeepEles.size()==1 ) {
        histo1d_["3P1F_P1_Et"]->Fill(acceptEles.front()->et(), aWeight);
        histo1d_["3P1F_P1_eta"]->Fill(acceptEles.front()->eta(), aWeight);
        histo1d_["3P1F_P1_phi"]->Fill(acceptEles.front()->phi(), aWeight);
        histo1d_["3P1F_P2_Et"]->Fill(acceptEles.at(1)->et(), aWeight);
        histo1d_["3P1F_P2_eta"]->Fill(acceptEles.at(1)->eta(), aWeight);
        histo1d_["3P1F_P2_phi"]->Fill(acceptEles.at(1)->phi(), aWeight);
        histo1d_["3P1F_P3_Et"]->Fill(acceptEles.at(2)->et(), aWeight);
        histo1d_["3P1F_P3_eta"]->Fill(acceptEles.at(2)->eta(), aWeight);
        histo1d_["3P1F_P3_phi"]->Fill(acceptEles.at(2)->phi(), aWeight);

        histo1d_["3P1F_F1_Et"]->Fill(nonHeepEles.front()->et(), aWeight);
        histo1d_["3P1F_F1_eta"]->Fill(nonHeepEles.front()->eta(), aWeight);
        histo1d_["3P1F_F1_phi"]->Fill(nonHeepEles.front()->phi(), aWeight);

        std::vector<pat::ElectronRef> allEles(acceptEles);
        allEles.insert( allEles.end(), nonHeepEles.begin(), nonHeepEles.end() );
        std::pair<pat::ElectronRef,pat::ElectronRef> pair1, pair2;

        bool paired = PairingHelper::pair4E(allEles,addGsfTrkMap,pair1,pair2);

        if (paired) {
          const auto lvecE1 = pair1.first->polarP4()*pair1.first->userFloat("ecalTrkEnergyPostCorr")/pair1.first->energy();
          const auto lvecE2 = pair1.second->polarP4()*pair1.second->userFloat("ecalTrkEnergyPostCorr")/pair1.second->energy();
          const auto lvecE3 = pair2.first->polarP4()*pair2.first->userFloat("ecalTrkEnergyPostCorr")/pair2.first->energy();
          const auto lvecE4 = pair2.second->polarP4()*pair2.second->userFloat("ecalTrkEnergyPostCorr")/pair2.second->energy();

          const auto lvecA1 = lvecE1 + lvecE2;
          const auto lvecA2 = lvecE3 + lvecE4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
          const double dr2ll2 = reco::deltaR2(lvecE3.eta(),lvecE3.phi(),lvecE4.eta(),lvecE4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();

          histo1d_["3P1F_llll_invM"]->Fill(m4l, aWeight);
          histo1d_["3P1F_llll_pt"]->Fill(lvec4l.pt(), aWeight);
          histo1d_["3P1F_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
          histo1d_["3P1F_ll1_invM"]->Fill(lvecA1.M(), aWeight);
          histo1d_["3P1F_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
          histo1d_["3P1F_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
          histo1d_["3P1F_ll2_invM"]->Fill(lvecA2.M(), aWeight);
          histo1d_["3P1F_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
          histo1d_["3P1F_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

          if ( m4l > 50. && m4l < 200. ) {
            const double drFake = std::sqrt( std::min({ reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(2)->eta(),acceptEles.at(2)->phi()),
                                                        reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                                        reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) }) );
            const TF1* func = nonHeepEles.front()->isEB() ? ( drFake < 0.3 ? ff_EB_below_ : ff_EB_above_ ) : ( drFake < 0.3 ? ff_EE_below_ : ff_EE_above_ );
            const double ff = func->Eval(drFake);

            const auto& fakePair = ( nonHeepEles.front()==pair1.first || nonHeepEles.front()==pair1.second ) ? pair1 : pair2;
            const auto& passPartner = ( nonHeepEles.front()==fakePair.first ) ? fakePair.second : fakePair.first;

            histo1d_["3P1F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["3P1F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
            histo1d_["3P1F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
            histo1d_["3P1F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
            histo1d_["3P1F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
            histo1d_["3P1F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
            histo1d_["3P1F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
            histo1d_["3P1F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
            histo1d_["3P1F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

            histo1d_["3P1F_CR_Et_P1"]->Fill(acceptEles.front()->et(), aWeight);
            histo1d_["3P1F_CR_Et_P2"]->Fill(acceptEles.at(1)->et(), aWeight);
            histo1d_["3P1F_CR_Et_P3"]->Fill(acceptEles.at(2)->et(), aWeight);
            histo1d_["3P1F_CR_Et_F1"]->Fill(nonHeepEles.front()->et(), aWeight);

            const double drPass = std::sqrt( reco::deltaR2(passPartner->eta(),passPartner->phi(),nonHeepEles.front()->eta(),nonHeepEles.front()->phi()) );

            if ( std::abs( passPartner->superCluster()->eta() ) < 1.5 ) {
              histo1d_["3P1F_CR_Et_EB"]->Fill(passPartner->et(), aWeight);
              histo1d_["3P1F_CR_dr_EB"]->Fill(drPass, aWeight);
            } else {
              histo1d_["3P1F_CR_Et_EE"]->Fill(nonHeepEles.front()->et(), aWeight);
              histo1d_["3P1F_CR_dr_EE"]->Fill(drPass, aWeight);
            }

            histo1d_["3P1F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff);
            histo1d_["3P1F_CR_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff);

            histo1d_["3P1F_CR_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff);
            histo1d_["3P1F_CR_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff);
            histo1d_["3P1F_CR_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff);
            histo1d_["3P1F_CR_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff);
          } // CR
        } // paired
      } // 3P1F

      break;
    case 2:
      // 2P1F - FF denominator
      if ( nonHeepEles.size()==1 && acceptEles.front()->charge()*acceptEles.at(1)->charge() < 0 ) {
        const auto lvecE1 = acceptEles.front()->polarP4()*acceptEles.front()->userFloat("ecalTrkEnergyPostCorr")/acceptEles.front()->energy();
        const auto lvecE2 = acceptEles.at(1)->polarP4()*acceptEles.at(1)->userFloat("ecalTrkEnergyPostCorr")/acceptEles.at(1)->energy();
        const auto lvecll = lvecE1 + lvecE2;
        const double mll = lvecll.M();
        const double dr2ll = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());

        if ( mll > 84.19 && mll < 98.19 && acceptEles.front()->charge()*acceptEles.at(1)->charge() < 0 ) {
          histo1d_["2P1F_P1_Et"]->Fill(acceptEles.front()->et(), aWeight);
          histo1d_["2P1F_P1_eta"]->Fill(acceptEles.front()->eta(), aWeight);
          histo1d_["2P1F_P1_phi"]->Fill(acceptEles.front()->phi(), aWeight);
          histo1d_["2P1F_P2_Et"]->Fill(acceptEles.at(1)->et(), aWeight);
          histo1d_["2P1F_P2_eta"]->Fill(acceptEles.at(1)->eta(), aWeight);
          histo1d_["2P1F_P2_phi"]->Fill(acceptEles.at(1)->phi(), aWeight);

          histo1d_["2P1F_F1_Et"]->Fill(nonHeepEles.front()->et(), aWeight);
          histo1d_["2P1F_F1_eta"]->Fill(nonHeepEles.front()->eta(), aWeight);
          histo1d_["2P1F_F1_phi"]->Fill(nonHeepEles.front()->phi(), aWeight);

          histo1d_["2P1F_ll_invM"]->Fill(mll, aWeight);
          histo1d_["2P1F_ll_pt"]->Fill(lvecll.pt(), aWeight);
          histo1d_["2P1F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

          const double mindr2 = std::min( reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                          reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) );

          if ( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 ) {
            histo1d_["2P1F_Et_EB"]->Fill(nonHeepEles.front()->et(), aWeight);
            histo1d_["2P1F_dr_EB"]->Fill(std::sqrt(mindr2), aWeight);
          } else {
            histo1d_["2P1F_Et_EE"]->Fill(nonHeepEles.front()->et(), aWeight);
            histo1d_["2P1F_dr_EE"]->Fill(std::sqrt(mindr2), aWeight);
          }
        } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
      } // 2P1F

      // 2P2F
      if ( nonHeepEles.size()==2 ) {
        histo1d_["2P2F_P1_Et"]->Fill(acceptEles.front()->et(), aWeight);
        histo1d_["2P2F_P1_eta"]->Fill(acceptEles.front()->eta(), aWeight);
        histo1d_["2P2F_P1_phi"]->Fill(acceptEles.front()->phi(), aWeight);
        histo1d_["2P2F_P2_Et"]->Fill(acceptEles.at(1)->et(), aWeight);
        histo1d_["2P2F_P2_eta"]->Fill(acceptEles.at(1)->eta(), aWeight);
        histo1d_["2P2F_P2_phi"]->Fill(acceptEles.at(1)->phi(), aWeight);

        histo1d_["2P2F_F1_Et"]->Fill(nonHeepEles.front()->et(), aWeight);
        histo1d_["2P2F_F1_eta"]->Fill(nonHeepEles.front()->eta(), aWeight);
        histo1d_["2P2F_F1_phi"]->Fill(nonHeepEles.front()->phi(), aWeight);
        histo1d_["2P2F_F2_Et"]->Fill(nonHeepEles.at(1)->et(), aWeight);
        histo1d_["2P2F_F2_eta"]->Fill(nonHeepEles.at(1)->eta(), aWeight);
        histo1d_["2P2F_F2_phi"]->Fill(nonHeepEles.at(1)->phi(), aWeight);

        std::vector<pat::ElectronRef> allEles(acceptEles);
        allEles.insert( allEles.end(), nonHeepEles.begin(), nonHeepEles.end() );
        std::pair<pat::ElectronRef,pat::ElectronRef> pair1, pair2;

        bool paired = PairingHelper::pair4E(allEles,addGsfTrkMap,pair1,pair2);
        bool pairedFakes = ( ( pair1.first==nonHeepEles.front() && pair1.second==nonHeepEles.at(1) ) ||
                             ( pair1.second==nonHeepEles.front() && pair1.first==nonHeepEles.at(1) ) ||
                             ( pair2.first==nonHeepEles.front() && pair2.second==nonHeepEles.at(1) ) ||
                             ( pair2.second==nonHeepEles.front() && pair2.first==nonHeepEles.at(1) ) );

        if (paired) {
          const auto lvecE1 = pair1.first->polarP4()*pair1.first->userFloat("ecalTrkEnergyPostCorr")/pair1.first->energy();
          const auto lvecE2 = pair1.second->polarP4()*pair1.second->userFloat("ecalTrkEnergyPostCorr")/pair1.second->energy();
          const auto lvecE3 = pair2.first->polarP4()*pair2.first->userFloat("ecalTrkEnergyPostCorr")/pair2.first->energy();
          const auto lvecE4 = pair2.second->polarP4()*pair2.second->userFloat("ecalTrkEnergyPostCorr")/pair2.second->energy();

          const auto lvecA1 = lvecE1 + lvecE2;
          const auto lvecA2 = lvecE3 + lvecE4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
          const double dr2ll2 = reco::deltaR2(lvecE3.eta(),lvecE3.phi(),lvecE4.eta(),lvecE4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();

          histo1d_["2P2F_llll_invM"]->Fill(m4l, aWeight);
          histo1d_["2P2F_llll_pt"]->Fill(lvec4l.pt(), aWeight);
          histo1d_["2P2F_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
          histo1d_["2P2F_ll1_invM"]->Fill(lvecA1.M(), aWeight);
          histo1d_["2P2F_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
          histo1d_["2P2F_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
          histo1d_["2P2F_ll2_invM"]->Fill(lvecA2.M(), aWeight);
          histo1d_["2P2F_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
          histo1d_["2P2F_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

          if ( m4l > 50. && m4l < 200. ) {
            const double drFake1 = std::sqrt( std::min({ reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi()),
                                                         reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                                         reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) }) );
            const double drFake2 = std::sqrt( std::min({ reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi()),
                                                         reco::deltaR2(nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                                         reco::deltaR2(nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) }) );
            const TF1* func1 = nonHeepEles.front()->isEB() ? ( drFake1 < 0.3 ? ff_EB_below_ : ff_EB_above_ ) : ( drFake1 < 0.3 ? ff_EE_below_ : ff_EE_above_ );
            const TF1* func2 = nonHeepEles.at(1)->isEB() ? ( drFake2 < 0.3 ? ff_EB_below_ : ff_EB_above_ ) : ( drFake2 < 0.3 ? ff_EE_below_ : ff_EE_above_ );
            const double ff1 = func1->Eval(drFake1);
            const double ff2 = func2->Eval(drFake2);
            const double drFake = std::sqrt( reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi() ) );
            const TF1* func3 = nonHeepEles.front()->isEB() ? ( drFake < 0.3 ? ff_PPFF_EB_below_ : ff_PPFF_EB_above_ ) : ( drFake < 0.3 ? ff_PPFF_EE_below_ : ff_PPFF_EE_above_ );
            const double gg = func3->Eval(drFake);

            if (!pairedFakes) {
              histo1d_["2P2F_CR_PFPF_llll_invM"]->Fill(m4l, aWeight);
              histo1d_["2P2F_CR_PFPF_llll_pt"]->Fill(lvec4l.pt(), aWeight);
              histo1d_["2P2F_CR_PFPF_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
              histo1d_["2P2F_CR_PFPF_ll1_invM"]->Fill(lvecA1.M(), aWeight);
              histo1d_["2P2F_CR_PFPF_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
              histo1d_["2P2F_CR_PFPF_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
              histo1d_["2P2F_CR_PFPF_ll2_invM"]->Fill(lvecA2.M(), aWeight);
              histo1d_["2P2F_CR_PFPF_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
              histo1d_["2P2F_CR_PFPF_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

              histo1d_["2P2F_CR_PFPF_llll_invM_xFF"]->Fill(m4l, aWeight*ff1);
              histo1d_["2P2F_CR_PFPF_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff1);
              histo1d_["2P2F_CR_PFPF_llll_invM_xFF"]->Fill(m4l, aWeight*ff2);
              histo1d_["2P2F_CR_PFPF_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff2);

              histo1d_["2P2F_CR_PFPF_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff1);
              histo1d_["2P2F_CR_PFPF_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff1);
              histo1d_["2P2F_CR_PFPF_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff2);
              histo1d_["2P2F_CR_PFPF_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff2);

              histo1d_["2P2F_CR_PFPF_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff1);
              histo1d_["2P2F_CR_PFPF_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff1);
              histo1d_["2P2F_CR_PFPF_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff2);
              histo1d_["2P2F_CR_PFPF_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff2);

              histo1d_["2P2F_CR_PFPF_llll_invM_xFF2"]->Fill(m4l, aWeight*ff1*ff2);
              histo1d_["2P2F_CR_PFPF_ll1ll2_dr_xFF2"]->Fill(std::sqrt(dr2A12), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_PFPF_ll1_invM_xFF2"]->Fill(lvecA1.M(), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_PFPF_ll1_dr_xFF2"]->Fill(std::sqrt(dr2ll1), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_PFPF_ll2_invM_xFF2"]->Fill(lvecA2.M(), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_PFPF_ll2_dr_xFF2"]->Fill(std::sqrt(dr2ll2), aWeight*ff1*ff2);

              if ( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 ) {
                histo1d_["2P2F_CR_PFPF_Et_EB_xFF"]->Fill(nonHeepEles.front()->et(), aWeight*ff1);
                histo1d_["2P2F_CR_PFPF_dr_EB_xFF"]->Fill(drFake, aWeight*ff1);
              } else {
                histo1d_["2P2F_CR_PFPF_Et_EE_xFF"]->Fill(nonHeepEles.front()->et(), aWeight*ff1);
                histo1d_["2P2F_CR_PFPF_dr_EE_xFF"]->Fill(drFake, aWeight*ff1);
              }

              if ( std::abs( nonHeepEles.at(1)->superCluster()->eta() ) < 1.5 ) {
                histo1d_["2P2F_CR_PFPF_Et_EB_xFF"]->Fill(nonHeepEles.at(1)->et(), aWeight*ff2);
                histo1d_["2P2F_CR_PFPF_dr_EB_xFF"]->Fill(drFake, aWeight*ff2);
              } else {
                histo1d_["2P2F_CR_PFPF_Et_EE_xFF"]->Fill(nonHeepEles.at(1)->et(), aWeight*ff2);
                histo1d_["2P2F_CR_PFPF_dr_EE_xFF"]->Fill(drFake, aWeight*ff2);
              }
            }

            if (pairedFakes) {
              histo1d_["2P2F_CR_PPFF_llll_invM"]->Fill(m4l, aWeight);
              histo1d_["2P2F_CR_PPFF_llll_pt"]->Fill(lvec4l.pt(), aWeight);
              histo1d_["2P2F_CR_PPFF_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
              histo1d_["2P2F_CR_PPFF_ll1_invM"]->Fill(lvecA1.M(), aWeight);
              histo1d_["2P2F_CR_PPFF_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
              histo1d_["2P2F_CR_PPFF_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
              histo1d_["2P2F_CR_PPFF_ll2_invM"]->Fill(lvecA2.M(), aWeight);
              histo1d_["2P2F_CR_PPFF_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
              histo1d_["2P2F_CR_PPFF_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

              histo1d_["2P2F_CR_PPFF_llll_invM_xFF"]->Fill(m4l, aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_llll_invM_xFF"]->Fill(m4l, aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*gg);

              histo1d_["2P2F_CR_PPFF_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*gg);

              histo1d_["2P2F_CR_PPFF_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*gg);
              histo1d_["2P2F_CR_PPFF_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*gg);

              histo1d_["2P2F_CR_PPFF_llll_invM_xFF2"]->Fill(m4l, aWeight*ff1*gg);
              histo1d_["2P2F_CR_PPFF_ll1ll2_dr_xFF2"]->Fill(std::sqrt(dr2A12), aWeight*ff1*gg);
              histo1d_["2P2F_CR_PPFF_ll1_invM_xFF2"]->Fill(lvecA1.M(), aWeight*ff1*gg);
              histo1d_["2P2F_CR_PPFF_ll1_dr_xFF2"]->Fill(std::sqrt(dr2ll1), aWeight*ff1*gg);
              histo1d_["2P2F_CR_PPFF_ll2_invM_xFF2"]->Fill(lvecA2.M(), aWeight*ff1*gg);
              histo1d_["2P2F_CR_PPFF_ll2_dr_xFF2"]->Fill(std::sqrt(dr2ll2), aWeight*ff1*gg);

              if ( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 ) {
                histo1d_["2P2F_CR_PPFF_Et_EB"]->Fill(nonHeepEles.front()->et(), aWeight);
                histo1d_["2P2F_CR_PPFF_dr_EB"]->Fill(drFake, aWeight);
              } else {
                histo1d_["2P2F_CR_PPFF_Et_EE"]->Fill(nonHeepEles.front()->et(), aWeight);
                histo1d_["2P2F_CR_PPFF_dr_EE"]->Fill(drFake, aWeight);
              }

              if ( std::abs( nonHeepEles.at(1)->superCluster()->eta() ) < 1.5 ) {
                histo1d_["2P2F_CR_PPFF_Et_EB"]->Fill(nonHeepEles.at(1)->et(), aWeight);
                histo1d_["2P2F_CR_PPFF_dr_EB"]->Fill(drFake, aWeight);
              } else {
                histo1d_["2P2F_CR_PPFF_Et_EE"]->Fill(nonHeepEles.at(1)->et(), aWeight);
                histo1d_["2P2F_CR_PPFF_dr_EE"]->Fill(drFake, aWeight);
              }
            } // pairedFakes

            histo1d_["2P2F_CR_Et_P1"]->Fill(acceptEles.front()->et(), aWeight);
            histo1d_["2P2F_CR_Et_P2"]->Fill(acceptEles.at(1)->et(), aWeight);
            histo1d_["2P2F_CR_Et_F1"]->Fill(nonHeepEles.front()->et(), aWeight);
            histo1d_["2P2F_CR_Et_F2"]->Fill(nonHeepEles.at(1)->et(), aWeight);
          } // CR
        } // paired
      } // 2P2F

      break;
      // case acceptEles.size()==2
  } // switch acceptEles.size()
}

DEFINE_FWK_MODULE(ResolvedEleCRanalyzer);
