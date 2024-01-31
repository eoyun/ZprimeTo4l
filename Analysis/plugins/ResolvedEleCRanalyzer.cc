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
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "ZprimeTo4l/Analysis/interface/PairingHelper.h"
#include "ZprimeTo4l/Analysis/interface/ElectronSystematicsHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"

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
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;

  const edm::FileInPath recoSFpath_;
  const edm::FileInPath FFpath_;
  const edm::FileInPath purwgtPath_;
  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;
  const double etThres1_;
  const double etThres2_;
  const double ffSystCL95_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  std::unique_ptr<TFile> recoSFfile_;
  TH2D* recoSF_;

  std::unique_ptr<TFile> FFfile_;
  TF1* ff_dr_;
  TFitResultPtr ffFitResult_;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  ElectronSystematicsHelper systHelperEle_;
};

ResolvedEleCRanalyzer::ResolvedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
recoSFpath_(iConfig.getParameter<edm::FileInPath>("recoSFpath")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
etThres1_(iConfig.getParameter<double>("etThres1")),
etThres2_(iConfig.getParameter<double>("etThres2")),
ffSystCL95_(iConfig.getParameter<double>("ffSystCL95")),
systHelperEle_(ElectronSystematicsHelper(iConfig.getParameter<edm::FileInPath>("modHeepSFpath"))) {
  usesResource("TFileService");

  systHelperEle_.SetModifiedHeepSF(iConfig.getParameter<double>("modHeepSFmuEB1"),
                                   iConfig.getParameter<double>("modHeepSFmuEB2"),
                                   iConfig.getParameter<double>("modHeepSFmuEE"));
  systHelperEle_.SetModifiedHeepSFcl95(iConfig.getParameter<double>("modHeepSFcl95EB1"),
                                       iConfig.getParameter<double>("modHeepSFcl95EB2"),
                                       iConfig.getParameter<double>("modHeepSFcl95EE"));
  systHelperEle_.SetModifiedHeepSFupper(iConfig.getParameter<double>("modHeepSFupperEB1"),
                                        iConfig.getParameter<double>("modHeepSFupperEB2"),
                                        iConfig.getParameter<double>("modHeepSFupperEE"));
  systHelperEle_.SetModifiedHeepPol1(iConfig.getParameter<std::string>("modHeepPolEB1str"),
                                     iConfig.getParameter<std::string>("modHeepPolEB2str"),
                                     iConfig.getParameter<std::string>("modHeepPolEEstr"));
}

void ResolvedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  recoSFfile_ = std::make_unique<TFile>(recoSFpath_.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  ff_dr_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFF_dr_all"));
  ffFitResult_ = (static_cast<TH1D*>(FFfile_->Get("all_dr_3P0F_rebin")))->Fit(ff_dr_,"RS");

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow_4E"] = fs->make<TH1D>("cutflow_4E","cutflow (4E)",10,0.,10.);

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
  histo1d_["3P0F_dr"] = fs->make<TH1D>("3P0F_dr","3P0F Pass min dR(ll)",128,0.,6.4);

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
  histo1d_["2P1F_dr"] = fs->make<TH1D>("2P1F_dr","2P1F Fail min dR(ll)",128,0.,6.4);

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

  histo1d_["2P2F_CR_Et_EB"] = fs->make<TH1D>("2P2F_CR_Et_EB","2P2F CR Fail E_{T} (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_Et_EE"] = fs->make<TH1D>("2P2F_CR_Et_EE","2P2F CR Fail E_{T} (EE);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_dr_EB"] = fs->make<TH1D>("2P2F_CR_dr_EB","2P2F CR Fail min dR(ll) (EB)",128,0.,6.4);
  histo1d_["2P2F_CR_dr_EE"] = fs->make<TH1D>("2P2F_CR_dr_EE","2P2F CR Fail min dR(ll) (EE)",128,0.,6.4);

  histo1d_["2P2F_CR_Et_EB_xFF"] = fs->make<TH1D>("2P2F_CR_Et_EB_xFF","2P2F CR Fail E_{T} x Fake Factor (EB);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_Et_EE_xFF"] = fs->make<TH1D>("2P2F_CR_Et_EE_xFF","2P2F CR Fail E_{T} x Fake Factor (EE);E_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_dr_EB_xFF"] = fs->make<TH1D>("2P2F_CR_dr_EB_xFF","2P2F CR Fail min dR(ll) x Fake Factor (EB)",128,0.,6.4);
  histo1d_["2P2F_CR_dr_EE_xFF"] = fs->make<TH1D>("2P2F_CR_dr_EE_xFF","2P2F CR Fail min dR(ll) x Fake Factor (EE)",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM"] = fs->make<TH1D>("2P2F_CR_llll_invM","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_idUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_idUp","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_idDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_idDn","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_pt"] = fs->make<TH1D>("2P2F_CR_llll_pt","2P2F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1ll2_dr"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr","2P2F CR dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM"] = fs->make<TH1D>("2P2F_CR_ll1_invM","2P2F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_pt"] = fs->make<TH1D>("2P2F_CR_ll1_pt","2P2F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr"] = fs->make<TH1D>("2P2F_CR_ll1_dr","2P2F CR dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM"] = fs->make<TH1D>("2P2F_CR_ll2_invM","2P2F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_pt"] = fs->make<TH1D>("2P2F_CR_ll2_pt","2P2F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr"] = fs->make<TH1D>("2P2F_CR_ll2_dr","2P2F CR dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM_xFF"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_idUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_idUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_idDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_idDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF","2P2F CR dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF","2P2F CR M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF","2P2F CR dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF","2P2F CR M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF","2P2F CR dR(ll2) x Fake factor",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffUp","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffDn","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_idUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_idUp","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_idDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_idDn","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF2","2P2F CR dR(ll1ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF2","2P2F CR M(ll1) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF2","2P2F CR dR(ll1) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF2","2P2F CR M(ll2) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF2","2P2F CR dR(ll2) x Fake factor^2",128,0.,6.4);

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

  histo1d_["3P1F_CR_llll_invM"] = fs->make<TH1D>("3P1F_CR_llll_invM","3P1F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_idUp"] = fs->make<TH1D>("3P1F_CR_llll_invM_idUp","3P1F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_idDn"] = fs->make<TH1D>("3P1F_CR_llll_invM_idDn","3P1F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_pt"] = fs->make<TH1D>("3P1F_CR_llll_pt","3P1F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1ll2_dr"] = fs->make<TH1D>("3P1F_CR_ll1ll2_dr","3P1F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["3P1F_CR_ll1_invM"] = fs->make<TH1D>("3P1F_CR_ll1_invM","3P1F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_pt"] = fs->make<TH1D>("3P1F_CR_ll1_pt","3P1F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_dr"] = fs->make<TH1D>("3P1F_CR_ll1_dr","3P1F CR dR(ll1)",64,0.,6.4);
  histo1d_["3P1F_CR_ll2_invM"] = fs->make<TH1D>("3P1F_CR_ll2_invM","3P1F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_pt"] = fs->make<TH1D>("3P1F_CR_ll2_pt","3P1F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_dr"] = fs->make<TH1D>("3P1F_CR_ll2_dr","3P1F CR dR(ll2)",64,0.,6.4);

  histo1d_["3P1F_CR_llll_invM_xFF"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffUp"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffUp","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffDn"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffDn","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_idUp"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_idUp","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_idDn"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_idDn","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
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

  histo1d_["4P0F_CR_llll_invM"] = fs->make<TH1D>("4P0F_CR_llll_invM","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_scaleUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_scaleUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_scaleDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_scaleDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_sigmaUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_sigmaUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_sigmaDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_sigmaDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_idUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_idUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_idDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_idDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
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
  FFfile_->Close();
  purwgtFile_->Close();
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

    edm::Handle<edm::View<PileupSummaryInfo>> pusummary;
    iEvent.getByToken(pileupToken_, pusummary);

    for (unsigned int idx = 0; idx < pusummary->size(); ++idx) {
      const auto& apu = pusummary->refAt(idx);

      int bx = apu->getBunchCrossing();

      if (bx==0) { // in-time PU only
        aWeight *= purwgt_->GetBinContent( purwgt_->FindBin(apu->getTrueNumInteractions()) );

        break;
      }
    }
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);
  histo1d_["cutflow_4E"]->Fill(0.5,aWeight);

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

  histo1d_["cutflow_4E"]->Fill(1.5,aWeight);

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

  histo1d_["cutflow_4E"]->Fill(2.5,aWeight);

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
    if ( std::abs(aEle->superCluster()->eta()) > 1.4442 && std::abs(aEle->superCluster()->eta()) < 1.566 )
      continue;

    bool selected = false;

    if ( aEle->electronID("modifiedHeepElectronID") ) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      acceptEles.push_back(castEle);
      selected = true;
    }

    if (!selected && aEle->et() > etThres2_ && aEle->ecalDriven()) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      nonHeepEles.push_back(castEle);
    }
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByEt);
  std::sort(nonHeepEles.begin(),nonHeepEles.end(),sortByEt);

  histo1d_["cutflow_4E"]->Fill(3.5,aWeight);

  bool accepted = false;

  if ( acceptEles.size() > 1 ) {
    if ( acceptEles.at(1)->et() > etThres1_ )
      accepted = true;
  }

  if (!accepted)
    return;

  histo1d_["cutflow_4E"]->Fill(4.5,aWeight);

  bool trigMatched1 = false;
  bool trigMatched2 = false;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()) < 0.01 )
      trigMatched1 = true;

    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) < 0.01 )
      trigMatched2 = true;
  }

  if ( !trigMatched1 || !trigMatched2 )
    return;

  histo1d_["cutflow_4E"]->Fill(5.5,aWeight);

  if ( acceptEles.size()==4 )
    histo1d_["cutflow_4E"]->Fill(6.5,aWeight);

  // SF
  if (isMC_) {
    for (const auto& aEle : acceptEles) {
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= recoSF_->GetBinContent( recoSF_->FindBin(aEle->superCluster()->eta(), apt ) );
      aWeight *= systHelperEle_.GetModifiedHeepSF(aEle);
    }

    for (const auto& aEle : nonHeepEles) {
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= recoSF_->GetBinContent( recoSF_->FindBin(aEle->superCluster()->eta(), apt ) );
    }
  }

  auto modHeepSFcl95 = [this,&acceptEles] (const double wgt) -> std::pair<double,double> {
    if (!isMC_)
      return std::make_pair(wgt,wgt);

    double wgtUp = wgt;
    double wgtDn = wgt;

    for (const auto& aEle : acceptEles) {
      auto apair = systHelperEle_.GetModifiedHeepSFcl95UpDn(aEle);
      wgtUp *= apair.first;
      wgtDn *= apair.second;
    } // acceptEles

    return std::make_pair(wgtUp,wgtDn);
  };

  auto ciFF = [this] (const double dr) {
    const double ff = ff_dr_->Eval(dr);
    const double xval[1] = {dr};
    double ci[1];
    ffFitResult_->GetConfidenceIntervals(1,1,0,xval,ci,0.95,false);

    return std::hypot(ci[0],ff*ffSystCL95_);
  };

  switch (acceptEles.size()) {
    case 4:
      // 4P0F
      histo1d_["cutflow_4E"]->Fill(7.5,aWeight);

      if (true) {
        std::vector<pat::ElectronRef> allEles(acceptEles);
        std::pair<pat::ElectronRef,pat::ElectronRef> pair1, pair2;

        bool paired = PairingHelper::pair4E(allEles,addGsfTrkMap,pair1,pair2);

        if (paired) {
          histo1d_["cutflow_4E"]->Fill(8.5,aWeight);

          auto lvecCorr = [] (const pat::ElectronRef& aEle, const std::string& opt) {
            return aEle->polarP4()*aEle->userFloat(opt)/aEle->energy();
          };

          const auto lvecE1 = lvecCorr(pair1.first,"ecalTrkEnergyPostCorr");
          const auto lvecE2 = lvecCorr(pair1.second,"ecalTrkEnergyPostCorr");
          const auto lvecE3 = lvecCorr(pair2.first,"ecalTrkEnergyPostCorr");
          const auto lvecE4 = lvecCorr(pair2.second,"ecalTrkEnergyPostCorr");

          const auto lvecE1_scaleUp = lvecCorr(pair1.first,"energyScaleUp");
          const auto lvecE2_scaleUp = lvecCorr(pair1.second,"energyScaleUp");
          const auto lvecE3_scaleUp = lvecCorr(pair2.first,"energyScaleUp");
          const auto lvecE4_scaleUp = lvecCorr(pair2.second,"energyScaleUp");
          const auto lvecE1_scaleDn = lvecCorr(pair1.first,"energyScaleDown");
          const auto lvecE2_scaleDn = lvecCorr(pair1.second,"energyScaleDown");
          const auto lvecE3_scaleDn = lvecCorr(pair2.first,"energyScaleDown");
          const auto lvecE4_scaleDn = lvecCorr(pair2.second,"energyScaleDown");
          const auto lvecE1_sigmaUp = lvecCorr(pair1.first,"energySigmaUp");
          const auto lvecE2_sigmaUp = lvecCorr(pair1.second,"energySigmaUp");
          const auto lvecE3_sigmaUp = lvecCorr(pair2.first,"energySigmaUp");
          const auto lvecE4_sigmaUp = lvecCorr(pair2.second,"energySigmaUp");
          const auto lvecE1_sigmaDn = lvecCorr(pair1.first,"energySigmaDown");
          const auto lvecE2_sigmaDn = lvecCorr(pair1.second,"energySigmaDown");
          const auto lvecE3_sigmaDn = lvecCorr(pair2.first,"energySigmaDown");
          const auto lvecE4_sigmaDn = lvecCorr(pair2.second,"energySigmaDown");

          const auto lvecA1 = lvecE1 + lvecE2;
          const auto lvecA2 = lvecE3 + lvecE4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
          const double dr2ll2 = reco::deltaR2(lvecE3.eta(),lvecE3.phi(),lvecE4.eta(),lvecE4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();
          const double m4l_scaleUp = (lvecE1_scaleUp+lvecE2_scaleUp+lvecE3_scaleUp+lvecE4_scaleUp).M();
          const double m4l_scaleDn = (lvecE1_scaleDn+lvecE2_scaleDn+lvecE3_scaleDn+lvecE4_scaleDn).M();
          const double m4l_sigmaUp = (lvecE1_sigmaUp+lvecE2_sigmaUp+lvecE3_sigmaUp+lvecE4_sigmaUp).M();
          const double m4l_sigmaDn = (lvecE1_sigmaDn+lvecE2_sigmaDn+lvecE3_sigmaDn+lvecE4_sigmaDn).M();

          if ( m4l > 50. && m4l < 200. )
            histo1d_["cutflow_4E"]->Fill(9.5,aWeight);

          if ( m4l > 50. && (m4l < 200. || isMC_) && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["4P0F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["4P0F_CR_llll_invM_scaleUp"]->Fill(m4l_scaleUp, aWeight);
            histo1d_["4P0F_CR_llll_invM_scaleDn"]->Fill(m4l_scaleDn, aWeight);
            histo1d_["4P0F_CR_llll_invM_sigmaUp"]->Fill(m4l_sigmaUp, aWeight);
            histo1d_["4P0F_CR_llll_invM_sigmaDn"]->Fill(m4l_sigmaDn, aWeight);
            histo1d_["4P0F_CR_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
            histo1d_["4P0F_CR_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
          }

          if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["cutflow_4E"]->Fill(10.5,aWeight);

            histo1d_["4P0F_CR_P1_eta"]->Fill(acceptEles.front()->eta(), aWeight);
            histo1d_["4P0F_CR_P1_phi"]->Fill(acceptEles.front()->phi(), aWeight);
            histo1d_["4P0F_CR_P2_eta"]->Fill(acceptEles.at(1)->eta(), aWeight);
            histo1d_["4P0F_CR_P2_phi"]->Fill(acceptEles.at(1)->phi(), aWeight);
            histo1d_["4P0F_CR_P3_eta"]->Fill(acceptEles.at(2)->eta(), aWeight);
            histo1d_["4P0F_CR_P3_phi"]->Fill(acceptEles.at(2)->phi(), aWeight);
            histo1d_["4P0F_CR_P4_eta"]->Fill(acceptEles.at(3)->eta(), aWeight);
            histo1d_["4P0F_CR_P4_phi"]->Fill(acceptEles.at(3)->phi(), aWeight);

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
      }

      break;
    case 3:
      // 3P0F - FF numerator
      if ( nonHeepEles.size()==0 ) {
        if ( acceptEles.at(2)->et() < etThres1_ ) {
          auto& first = acceptEles.front();
          auto& second = acceptEles.at(1);
          auto& probe = acceptEles.at(2);

          const auto lvecE1 = first->polarP4();
          const auto lvecE2 = second->polarP4();
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

            histo1d_["3P0F_dr"]->Fill( std::sqrt(mindr2), aWeight );

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

            const auto lvecE1 = first->polarP4();
            const auto lvecE2 = second->polarP4();
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

              histo1d_["3P0F_dr"]->Fill( std::sqrt(mindr2), aWeight );

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

          const double drFake = std::sqrt( std::min({ reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(2)->eta(),acceptEles.at(2)->phi()),
                                                      reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                                      reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) }) );
          const double ff = std::max(ff_dr_->Eval(drFake),0.);
          const double ci = ciFF(drFake);

          if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["3P1F_llll_pt"]->Fill(lvec4l.pt(), aWeight);
            histo1d_["3P1F_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
            histo1d_["3P1F_ll1_invM"]->Fill(lvecA1.M(), aWeight);
            histo1d_["3P1F_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
            histo1d_["3P1F_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
            histo1d_["3P1F_ll2_invM"]->Fill(lvecA2.M(), aWeight);
            histo1d_["3P1F_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
            histo1d_["3P1F_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

            histo1d_["3P1F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["3P1F_CR_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
            histo1d_["3P1F_CR_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
            histo1d_["3P1F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff);
            histo1d_["3P1F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff+ci));
            histo1d_["3P1F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff-std::min(ff,ci)));
            histo1d_["3P1F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff);
            histo1d_["3P1F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff);
          }

          if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            const auto& fakePair = ( nonHeepEles.front()==pair1.first || nonHeepEles.front()==pair1.second ) ? pair1 : pair2;
            const auto& passPartner = ( nonHeepEles.front()==fakePair.first ) ? fakePair.second : fakePair.first;

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
      if ( nonHeepEles.size()==1 ) {
        const auto lvecE1 = acceptEles.front()->polarP4();
        const auto lvecE2 = acceptEles.at(1)->polarP4();
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

          histo1d_["2P1F_dr"]->Fill(std::sqrt(mindr2), aWeight);

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

        if (paired) {
          bool pairedFakes = ( ( pair1.first==nonHeepEles.front() && pair1.second==nonHeepEles.at(1) ) ||
                               ( pair1.second==nonHeepEles.front() && pair1.first==nonHeepEles.at(1) ) ||
                               ( pair2.first==nonHeepEles.front() && pair2.second==nonHeepEles.at(1) ) ||
                               ( pair2.second==nonHeepEles.front() && pair2.first==nonHeepEles.at(1) ) );

          if (!pairedFakes) {
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

            const double drFake1 = std::sqrt( std::min({ reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi()),
                                                         reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                                         reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) }) );
            const double drFake2 = std::sqrt( std::min({ reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi()),
                                                         reco::deltaR2(nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()),
                                                         reco::deltaR2(nonHeepEles.at(1)->eta(),nonHeepEles.at(1)->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) }) );
            const double ff1 = std::max(ff_dr_->Eval(drFake1),0.);
            const double ff2 = std::max(ff_dr_->Eval(drFake2),0.);
            const double ci1 = ciFF(drFake1);
            const double ci2 = ciFF(drFake2);

            if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
              histo1d_["2P2F_llll_pt"]->Fill(lvec4l.pt(), aWeight);
              histo1d_["2P2F_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
              histo1d_["2P2F_ll1_invM"]->Fill(lvecA1.M(), aWeight);
              histo1d_["2P2F_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
              histo1d_["2P2F_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
              histo1d_["2P2F_ll2_invM"]->Fill(lvecA2.M(), aWeight);
              histo1d_["2P2F_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
              histo1d_["2P2F_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

              histo1d_["2P2F_CR_llll_invM"]->Fill(m4l, aWeight);
              histo1d_["2P2F_CR_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
              histo1d_["2P2F_CR_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
              histo1d_["2P2F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff1+ci1));
              histo1d_["2P2F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff2+ci2));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff1-std::min(ff1,ci1)));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff2-std::min(ff2,ci2)));
              histo1d_["2P2F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2"]->Fill(m4l, aWeight*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2_ffUp"]->Fill(m4l, aWeight*(ff1+ci1)*(ff2+ci2));
              histo1d_["2P2F_CR_llll_invM_xFF2_ffDn"]->Fill(m4l, aWeight*(ff1-std::min(ff1,ci1))*(ff2-std::min(ff2,ci2)));
              histo1d_["2P2F_CR_llll_invM_xFF2_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff1*ff2);
            }

            if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
              histo1d_["2P2F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
              histo1d_["2P2F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
              histo1d_["2P2F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
              histo1d_["2P2F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
              histo1d_["2P2F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
              histo1d_["2P2F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
              histo1d_["2P2F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
              histo1d_["2P2F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

              histo1d_["2P2F_CR_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff1);
              histo1d_["2P2F_CR_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff2);

              histo1d_["2P2F_CR_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff1);
              histo1d_["2P2F_CR_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff1);
              histo1d_["2P2F_CR_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff2);
              histo1d_["2P2F_CR_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff2);

              histo1d_["2P2F_CR_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff1);
              histo1d_["2P2F_CR_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff1);
              histo1d_["2P2F_CR_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff2);
              histo1d_["2P2F_CR_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff2);

              histo1d_["2P2F_CR_ll1ll2_dr_xFF2"]->Fill(std::sqrt(dr2A12), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_ll1_invM_xFF2"]->Fill(lvecA1.M(), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_ll1_dr_xFF2"]->Fill(std::sqrt(dr2ll1), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_ll2_invM_xFF2"]->Fill(lvecA2.M(), aWeight*ff1*ff2);
              histo1d_["2P2F_CR_ll2_dr_xFF2"]->Fill(std::sqrt(dr2ll2), aWeight*ff1*ff2);

              if ( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 ) {
                histo1d_["2P2F_CR_Et_EB"]->Fill(nonHeepEles.front()->et(), aWeight);
                histo1d_["2P2F_CR_dr_EB"]->Fill(drFake1, aWeight);
                histo1d_["2P2F_CR_Et_EB_xFF"]->Fill(nonHeepEles.front()->et(), aWeight*ff1);
                histo1d_["2P2F_CR_dr_EB_xFF"]->Fill(drFake1, aWeight*ff1);
              } else {
                histo1d_["2P2F_CR_Et_EE"]->Fill(nonHeepEles.front()->et(), aWeight);
                histo1d_["2P2F_CR_dr_EE"]->Fill(drFake1, aWeight);
                histo1d_["2P2F_CR_Et_EE_xFF"]->Fill(nonHeepEles.front()->et(), aWeight*ff1);
                histo1d_["2P2F_CR_dr_EE_xFF"]->Fill(drFake1, aWeight*ff1);
              }

              if ( std::abs( nonHeepEles.at(1)->superCluster()->eta() ) < 1.5 ) {
                histo1d_["2P2F_CR_Et_EB"]->Fill(nonHeepEles.at(1)->et(), aWeight);
                histo1d_["2P2F_CR_dr_EB"]->Fill(drFake2, aWeight);
                histo1d_["2P2F_CR_Et_EB_xFF"]->Fill(nonHeepEles.at(1)->et(), aWeight*ff2);
                histo1d_["2P2F_CR_dr_EB_xFF"]->Fill(drFake2, aWeight*ff2);
              } else {
                histo1d_["2P2F_CR_Et_EE"]->Fill(nonHeepEles.at(1)->et(), aWeight);
                histo1d_["2P2F_CR_dr_EE"]->Fill(drFake2, aWeight);
                histo1d_["2P2F_CR_Et_EE_xFF"]->Fill(nonHeepEles.at(1)->et(), aWeight*ff2);
                histo1d_["2P2F_CR_dr_EE_xFF"]->Fill(drFake2, aWeight*ff2);
              }

              histo1d_["2P2F_CR_Et_P1"]->Fill(acceptEles.front()->et(), aWeight);
              histo1d_["2P2F_CR_Et_P2"]->Fill(acceptEles.at(1)->et(), aWeight);
              histo1d_["2P2F_CR_Et_F1"]->Fill(nonHeepEles.front()->et(), aWeight);
              histo1d_["2P2F_CR_Et_F2"]->Fill(nonHeepEles.at(1)->et(), aWeight);
            } // CR
          } // !pairedFakes
        } // paired
      } // 2P2F

      break;
      // case acceptEles.size()==2
  } // switch acceptEles.size()
}

DEFINE_FWK_MODULE(ResolvedEleCRanalyzer);
