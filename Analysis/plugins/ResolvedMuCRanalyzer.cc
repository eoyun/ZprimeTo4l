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
#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"

class ResolvedMuCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedMuCRanalyzer(const edm::ParameterSet&);
  virtual ~ResolvedMuCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  math::PtEtaPhiMLorentzVector lvecFromTuneP(const pat::MuonRef& aMu);

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const std::vector<std::string> trigList_;

  const std::string trigHistName_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath purwgtPath_;

  const edm::FileInPath FFpath_;

  const double ptThres_;
  const double ptThresTrig_;

  MuonCorrectionHelper mucorrHelper_;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::unique_ptr<TFile> FFfile_;
  TF1* drFF_;
  TFitResultPtr ffFit_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  const double mumass_ = 0.1056583745;
};

ResolvedMuCRanalyzer::ResolvedMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
trigHistName_(iConfig.getParameter<std::string>("trigHistName")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptThresTrig_(iConfig.getParameter<double>("ptThresTrig")),
mucorrHelper_(rochesterPath_,triggerSFpath_,trigHistName_) {
  usesResource("TFileService");
}

math::PtEtaPhiMLorentzVector ResolvedMuCRanalyzer::lvecFromTuneP(const pat::MuonRef& aMu) {
  return math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                      aMu->tunePMuonBestTrack()->eta(),
                                      aMu->tunePMuonBestTrack()->phi(),
                                      mumass_);
}

void ResolvedMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  drFF_ = static_cast<TF1*>(FFfile_->FindObjectAny("RMFF_dr_MB"));
  ffFit_ = (static_cast<TH1D*>(FFfile_->Get("dr_3P0F_MB_rebin")))->Fit(drFF_,"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",10,0.,10.);

  // 3P0F
  histo1d_["3P0F_P1_pt"] = fs->make<TH1D>("3P0F_P1_pt","3P0F Pass M1 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["3P0F_P1_eta"] = fs->make<TH1D>("3P0F_P1_eta","3P0F Pass M1 #eta",200,-2.5,2.5);
  histo1d_["3P0F_P1_phi"] = fs->make<TH1D>("3P0F_P1_phi","3P0F Pass M1 #phi",128,-3.2,3.2);

  histo1d_["3P0F_P2_pt"] = fs->make<TH1D>("3P0F_P2_pt","3P0F Pass M2 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["3P0F_P2_eta"] = fs->make<TH1D>("3P0F_P2_eta","3P0F Pass M2 #eta",200,-2.5,2.5);
  histo1d_["3P0F_P2_phi"] = fs->make<TH1D>("3P0F_P2_phi","3P0F Pass M2 #phi",128,-3.2,3.2);

  histo1d_["3P0F_P3_pt"] = fs->make<TH1D>("3P0F_P3_pt","3P0F Pass M3 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["3P0F_P3_eta"] = fs->make<TH1D>("3P0F_P3_eta","3P0F Pass M3 #eta",200,-2.5,2.5);
  histo1d_["3P0F_P3_phi"] = fs->make<TH1D>("3P0F_P3_phi","3P0F Pass M3 #phi",128,-3.2,3.2);

  histo1d_["3P0F_ll_invM"] = fs->make<TH1D>("3P0F_ll_invM","3P0F P1P2 M(ll);M [GeV];",120,60.,120.);
  histo1d_["3P0F_ll_pt"] = fs->make<TH1D>("3P0F_ll_pt","3P0F P1P2 p_{T}(ll);p_{T} [GeV];",200,0.,200.);
  histo1d_["3P0F_ll_dr"] = fs->make<TH1D>("3P0F_ll_dr","3P0F P1P2 dR(ll)",128,0.,6.4);

  histo1d_["3P0F_pt_MB"] = fs->make<TH1D>("3P0F_pt_MB","3P0F Pass p_{T} (MB);p_{T} [GeV];",200,0.,1000.);
  histo1d_["3P0F_pt_ME"] = fs->make<TH1D>("3P0F_pt_ME","3P0F Pass p_{T} (ME);p_{T} [GeV];",200,0.,1000.);

  histo1d_["3P0F_dr_MB"] = fs->make<TH1D>("3P0F_dr_MB","3P0F Pass min dR(ll) (MB)",128,0.,6.4);
  histo1d_["3P0F_dr_ME"] = fs->make<TH1D>("3P0F_dr_ME","3P0F Pass min dR(ll) (ME)",128,0.,6.4);
  histo1d_["3P0F_dr"] = fs->make<TH1D>("3P0F_dr","3P0F Pass min dR(ll)",128,0.,6.4);

  // 2P1F
  histo1d_["2P1F_P1_pt"] = fs->make<TH1D>("2P1F_P1_pt","2P1F Pass M1 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P1F_P1_eta"] = fs->make<TH1D>("2P1F_P1_eta","2P1F Pass M1 #eta",200,-2.5,2.5);
  histo1d_["2P1F_P1_phi"] = fs->make<TH1D>("2P1F_P1_phi","2P1F Pass M1 #phi",128,-3.2,3.2);

  histo1d_["2P1F_P2_pt"] = fs->make<TH1D>("2P1F_P2_pt","2P1F Pass M2 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P1F_P2_eta"] = fs->make<TH1D>("2P1F_P2_eta","2P1F Pass M2 #eta",200,-2.5,2.5);
  histo1d_["2P1F_P2_phi"] = fs->make<TH1D>("2P1F_P2_phi","2P1F Pass M2 #phi",128,-3.2,3.2);

  histo1d_["2P1F_F1_pt"] = fs->make<TH1D>("2P1F_F1_pt","2P1F Fail M1 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P1F_F1_eta"] = fs->make<TH1D>("2P1F_F1_eta","2P1F Fail M1 #eta",200,-2.5,2.5);
  histo1d_["2P1F_F1_phi"] = fs->make<TH1D>("2P1F_F1_phi","2P1F Fail M1 #phi",128,-3.2,3.2);

  histo1d_["2P1F_ll_invM"] = fs->make<TH1D>("2P1F_ll_invM","2P1F P1P2 M(ll);M [GeV];",120,60.,120.);
  histo1d_["2P1F_ll_pt"] = fs->make<TH1D>("2P1F_ll_pt","2P1F P1P2 p_{T}(ll);p_{T} [GeV];",200,0.,200.);
  histo1d_["2P1F_ll_dr"] = fs->make<TH1D>("2P1F_ll_dr","2P1F P1P2 dR(ll)",128,0.,6.4);

  histo1d_["2P1F_pt_MB"] = fs->make<TH1D>("2P1F_pt_MB","2P1F Fail p_{T} (MB);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P1F_pt_ME"] = fs->make<TH1D>("2P1F_pt_ME","2P1F Fail p_{T} (ME);p_{T} [GeV];",200,0.,1000.);

  histo1d_["2P1F_dr_MB"] = fs->make<TH1D>("2P1F_dr_MB","2P1F Fail min dR(ll) (MB)",128,0.,6.4);
  histo1d_["2P1F_dr_ME"] = fs->make<TH1D>("2P1F_dr_ME","2P1F Fail min dR(ll) (ME)",128,0.,6.4);
  histo1d_["2P1F_dr"] = fs->make<TH1D>("2P1F_dr","2P1F Fail min dR(ll)",128,0.,6.4);

  // 2P2F
  histo1d_["2P2F_P1_pt"] = fs->make<TH1D>("2P2F_P1_pt","2P2F Pass M1 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_P1_eta"] = fs->make<TH1D>("2P2F_P1_eta","2P2F Pass M1 #eta",200,-2.5,2.5);
  histo1d_["2P2F_P1_phi"] = fs->make<TH1D>("2P2F_P1_phi","2P2F Pass M1 #phi",128,-3.2,3.2);

  histo1d_["2P2F_P2_pt"] = fs->make<TH1D>("2P2F_P2_pt","2P2F Pass M2 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_P2_eta"] = fs->make<TH1D>("2P2F_P2_eta","2P2F Pass M2 #eta",200,-2.5,2.5);
  histo1d_["2P2F_P2_phi"] = fs->make<TH1D>("2P2F_P2_phi","2P2F Pass M2 #phi",128,-3.2,3.2);

  histo1d_["2P2F_F1_pt"] = fs->make<TH1D>("2P2F_F1_pt","2P2F Fail M1 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_F1_eta"] = fs->make<TH1D>("2P2F_F1_eta","2P2F Fail M1 #eta",200,-2.5,2.5);
  histo1d_["2P2F_F1_phi"] = fs->make<TH1D>("2P2F_F1_phi","2P2F Fail M1 #phi",128,-3.2,3.2);

  histo1d_["2P2F_F2_pt"] = fs->make<TH1D>("2P2F_F2_pt","2P2F Fail M2 p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_F2_eta"] = fs->make<TH1D>("2P2F_F2_eta","2P2F Fail M2 #eta",200,-2.5,2.5);
  histo1d_["2P2F_F2_phi"] = fs->make<TH1D>("2P2F_F2_phi","2P2F Fail M2 #phi",128,-3.2,3.2);

  histo1d_["2P2F_llll_pt"] = fs->make<TH1D>("2P2F_llll_pt","2P2F p_{T}(4l);p_{T} [GeV];",100,0.,1000.);
  histo1d_["2P2F_ll1ll2_dr"] = fs->make<TH1D>("2P2F_ll1ll2_dr","2P2F dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_ll1_invM"] = fs->make<TH1D>("2P2F_ll1_invM","2P2F M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_ll1_pt"] = fs->make<TH1D>("2P2F_ll1_pt","2P2F p_{T}(ll1);p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_ll1_dr"] = fs->make<TH1D>("2P2F_ll1_dr","2P2F dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_ll2_invM"] = fs->make<TH1D>("2P2F_ll2_invM","2P2F M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_ll2_pt"] = fs->make<TH1D>("2P2F_ll2_pt","2P2F p_{T}(ll2);p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_ll2_dr"] = fs->make<TH1D>("2P2F_ll2_dr","2P2F dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_pt_P1"] = fs->make<TH1D>("2P2F_CR_pt_P1","2P2F CR Pass M1 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_CR_pt_P2"] = fs->make<TH1D>("2P2F_CR_pt_P2","2P2F CR Pass M2 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_CR_pt_F1"] = fs->make<TH1D>("2P2F_CR_pt_F1","2P2F CR Fail M1 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_CR_pt_F2"] = fs->make<TH1D>("2P2F_CR_pt_F2","2P2F CR Fail M2 p_{T};p_{T} [GeV];",100,0.,500.);

  histo1d_["2P2F_CR_pt_MB"] = fs->make<TH1D>("2P2F_CR_pt_MB","2P2F CR Fail p_{T} (MB);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_pt_ME"] = fs->make<TH1D>("2P2F_CR_pt_ME","2P2F CR Fail p_{T} (ME);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_dr_MB"] = fs->make<TH1D>("2P2F_CR_dr_MB","2P2F CR Fail min dR(ll) (MB)",128,0.,6.4);
  histo1d_["2P2F_CR_dr_ME"] = fs->make<TH1D>("2P2F_CR_dr_ME","2P2F CR Fail min dR(ll) (ME)",128,0.,6.4);

  histo1d_["2P2F_CR_pt_MB_xFF"] = fs->make<TH1D>("2P2F_CR_pt_MB_xFF","2P2F CR Fail p_{T} x Fake Factor (MB);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_pt_ME_xFF"] = fs->make<TH1D>("2P2F_CR_pt_ME_xFF","2P2F CR Fail p_{T} x Fake Factor (ME);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_dr_MB_xFF"] = fs->make<TH1D>("2P2F_CR_dr_MB_xFF","2P2F CR Fail min dR(ll) x Fake Factor (MB)",128,0.,6.4);
  histo1d_["2P2F_CR_dr_ME_xFF"] = fs->make<TH1D>("2P2F_CR_dr_ME_xFF","2P2F CR Fail min dR(ll) x Fake Factor (ME)",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM"] = fs->make<TH1D>("2P2F_CR_llll_invM","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_pt"] = fs->make<TH1D>("2P2F_CR_llll_pt","2P2F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1ll2_dr"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr","2P2F CR dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM"] = fs->make<TH1D>("2P2F_CR_ll1_invM","2P2F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_pt"] = fs->make<TH1D>("2P2F_CR_ll1_pt","2P2F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr"] = fs->make<TH1D>("2P2F_CR_ll1_dr","2P2F CR dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM"] = fs->make<TH1D>("2P2F_CR_ll2_invM","2P2F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_pt"] = fs->make<TH1D>("2P2F_CR_ll2_pt","2P2F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr"] = fs->make<TH1D>("2P2F_CR_ll2_dr","2P2F CR dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM_xFF"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF","2P2F CR dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF","2P2F CR M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF","2P2F CR dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF","2P2F CR M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF","2P2F CR dR(ll2) x Fake factor",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF2","2P2F CR dR(ll1ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF2","2P2F CR M(ll1) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF2","2P2F CR dR(ll1) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF2","2P2F CR M(ll2) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF2","2P2F CR dR(ll2) x Fake factor^2",128,0.,6.4);

  // 3P1F
  histo1d_["3P1F_P1_pt"] = fs->make<TH1D>("3P1F_P1_pt","3P1F Pass M1 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_P1_eta"] = fs->make<TH1D>("3P1F_P1_eta","3P1F Pass M1 #eta",100,-2.5,2.5);
  histo1d_["3P1F_P1_phi"] = fs->make<TH1D>("3P1F_P1_phi","3P1F Pass M1 #phi",64,-3.2,3.2);

  histo1d_["3P1F_P2_pt"] = fs->make<TH1D>("3P1F_P2_pt","3P1F Pass M2 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_P2_eta"] = fs->make<TH1D>("3P1F_P2_eta","3P1F Pass M2 #eta",100,-2.5,2.5);
  histo1d_["3P1F_P2_phi"] = fs->make<TH1D>("3P1F_P2_phi","3P1F Pass M2 #phi",64,-3.2,3.2);

  histo1d_["3P1F_P3_pt"] = fs->make<TH1D>("3P1F_P3_pt","3P1F Pass M3 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_P3_eta"] = fs->make<TH1D>("3P1F_P3_eta","3P1F Pass M3 #eta",100,-2.5,2.5);
  histo1d_["3P1F_P3_phi"] = fs->make<TH1D>("3P1F_P3_phi","3P1F Pass M3 #phi",64,-3.2,3.2);

  histo1d_["3P1F_F1_pt"] = fs->make<TH1D>("3P1F_F1_pt","3P1F Fail M1 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_F1_eta"] = fs->make<TH1D>("3P1F_F1_eta","3P1F Fail M1 #eta",100,-2.5,2.5);
  histo1d_["3P1F_F1_phi"] = fs->make<TH1D>("3P1F_F1_phi","3P1F Fail M1 #phi",64,-3.2,3.2);

  histo1d_["3P1F_llll_pt"] = fs->make<TH1D>("3P1F_llll_pt","3P1F p_{T}(4l);p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_ll1ll2_dr"] = fs->make<TH1D>("3P1F_ll1ll2_dr","3P1F dR(ll1ll2)",64,0.,6.4);
  histo1d_["3P1F_ll1_invM"] = fs->make<TH1D>("3P1F_ll1_invM","3P1F M(ll1);M [GeV];",200,0.,200.);
  histo1d_["3P1F_ll1_pt"] = fs->make<TH1D>("3P1F_ll1_pt","3P1F p_{T}(ll1);p_{T} [GeV];",200,0.,500.);
  histo1d_["3P1F_ll1_dr"] = fs->make<TH1D>("3P1F_ll1_dr","3P1F dR(ll1)",64,0.,6.4);
  histo1d_["3P1F_ll2_invM"] = fs->make<TH1D>("3P1F_ll2_invM","3P1F M(ll2);M [GeV];",200,0.,200.);
  histo1d_["3P1F_ll2_pt"] = fs->make<TH1D>("3P1F_ll2_pt","3P1F p_{T}(ll2);p_{T} [GeV];",200,0.,500.);
  histo1d_["3P1F_ll2_dr"] = fs->make<TH1D>("3P1F_ll2_dr","3P1F dR(ll2)",64,0.,6.4);

  histo1d_["3P1F_CR_pt_P1"] = fs->make<TH1D>("3P1F_CR_pt_P1","3P1F CR Pass M1 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_CR_pt_P2"] = fs->make<TH1D>("3P1F_CR_pt_P2","3P1F CR Pass M2 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_CR_pt_P3"] = fs->make<TH1D>("3P1F_CR_pt_P3","3P1F CR Pass M3 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_CR_pt_F1"] = fs->make<TH1D>("3P1F_CR_pt_F1","3P1F CR Fail M1 p_{T};p_{T} [GeV];",100,0.,500.);

  histo1d_["3P1F_CR_pt_MB"] = fs->make<TH1D>("3P1F_CR_pt_MB","3P1F CR Pass p_{T} (EB);p_{T} [GeV];",200,0.,1000.);
  histo1d_["3P1F_CR_pt_ME"] = fs->make<TH1D>("3P1F_CR_pt_ME","3P1F CR Pass p_{T} (EE);p_{T} [GeV];",200,0.,1000.);
  histo1d_["3P1F_CR_dr_MB"] = fs->make<TH1D>("3P1F_CR_dr_MB","3P1F CR Pass min dR(ll) (MB)",128,0.,6.4);
  histo1d_["3P1F_CR_dr_ME"] = fs->make<TH1D>("3P1F_CR_dr_ME","3P1F CR Pass min dR(ll) (ME)",128,0.,6.4);
  histo1d_["3P1F_CR_dr"] = fs->make<TH1D>("3P1F_CR_dr","3P1F CR Pass min dR(ll)",128,0.,6.4);

  histo1d_["3P1F_CR_llll_invM"] = fs->make<TH1D>("3P1F_CR_llll_invM","3P1F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_pt"] = fs->make<TH1D>("3P1F_CR_llll_pt","3P1F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1ll2_dr"] = fs->make<TH1D>("3P1F_CR_ll1ll2_dr","3P1F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["3P1F_CR_ll1_invM"] = fs->make<TH1D>("3P1F_CR_ll1_invM","3P1F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_pt"] = fs->make<TH1D>("3P1F_CR_ll1_pt","3P1F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_dr"] = fs->make<TH1D>("3P1F_CR_ll1_dr","3P1F CR dR(ll1)",64,0.,6.4);
  histo1d_["3P1F_CR_ll2_invM"] = fs->make<TH1D>("3P1F_CR_ll2_invM","3P1F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_pt"] = fs->make<TH1D>("3P1F_CR_ll2_pt","3P1F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_dr"] = fs->make<TH1D>("3P1F_CR_ll2_dr","3P1F CR dR(ll2)",64,0.,6.4);

  histo1d_["3P1F_CR_llll_invM_xFF"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
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

  histo1d_["4P0F_CR_pt_P1"] = fs->make<TH1D>("4P0F_CR_pt_P1","4P0F CR Pass M1 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_pt_P2"] = fs->make<TH1D>("4P0F_CR_pt_P2","4P0F CR Pass M2 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_pt_P3"] = fs->make<TH1D>("4P0F_CR_pt_P3","4P0F CR Pass M3 p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_pt_P4"] = fs->make<TH1D>("4P0F_CR_pt_P4","4P0F CR Pass M4 p_{T};p_{T} [GeV];",100,0.,500.);

  histo1d_["4P0F_CR_llll_invM"] = fs->make<TH1D>("4P0F_CR_llll_invM","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_pt"] = fs->make<TH1D>("4P0F_CR_llll_pt","4P0F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll1ll2_dr","4P0F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["4P0F_CR_ll1_invM"] = fs->make<TH1D>("4P0F_CR_ll1_invM","4P0F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_pt"] = fs->make<TH1D>("4P0F_CR_ll1_pt","4P0F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_dr"] = fs->make<TH1D>("4P0F_CR_ll1_dr","4P0F CR dR(ll1)",64,0.,6.4);
  histo1d_["4P0F_CR_ll2_invM"] = fs->make<TH1D>("4P0F_CR_ll2_invM","4P0F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_pt"] = fs->make<TH1D>("4P0F_CR_ll2_pt","4P0F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll2_dr","4P0F CR dR(ll2)",64,0.,6.4);
}

void ResolvedMuCRanalyzer::endJob() {
  purwgtFile_->Close();
  FFfile_->Close();
}

void ResolvedMuCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

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
  histo1d_["cutflow"]->Fill( 0.5, aWeight );

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

  if ( !pvHandle.isValid() || pvHandle->empty() )
    return;

  const reco::Vertex* primaryVertex = pvHandle->ptrAt(0).get();

  histo1d_["cutflow"]->Fill( 1.5, aWeight );

  std::vector<pat::ElectronRef> acceptEles;

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( std::abs(aEle->superCluster()->eta()) > 2.5 )
      continue;

    // veto EBEE gap
    if ( std::abs(aEle->superCluster()->eta()) > 1.4442 && std::abs(aEle->superCluster()->eta()) < 1.566 )
      continue;

    if ( aEle->electronID("modifiedHeepElectronID") )
      acceptEles.push_back(aEle.castTo<pat::ElectronRef>());
  }

  if ( !acceptEles.empty() )
    return;

  histo1d_["cutflow"]->Fill( 2.5, aWeight );

  if (muonHandle->empty())
    return;

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons; // but not highPtMuon
  std::vector<pat::MuonRef> nonHighPtMuons; // pass acceptance but not ID

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->tunePMuonBestTrack()->pt() < ptThres_ || std::abs(aMuon->eta()) > 2.4 )
      continue;

    if ( muon::isHighPtMuon(*aMuon,*primaryVertex) )
      highPtMuons.push_back(aMuon.castTo<pat::MuonRef>());
    else if ( muon::isTrackerHighPtMuon(*aMuon,*primaryVertex) )
      highPtTrackerMuons.push_back(aMuon.castTo<pat::MuonRef>());
  }

  std::vector<pat::MuonRef> isolatedHighPtMuons;
  std::vector<pat::MuonRef> isolatedHighPtTrackerMuons;

  auto sortByTuneP = [](const pat::MuonRef& a, const pat::MuonRef& b) {
    return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
  };

  for (unsigned idx = 0; idx < highPtMuons.size(); idx++) {
    const auto& firstMuon = highPtMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( reco::deltaR2(firstMuon->eta(),firstMuon->phi(),secondMuon->eta(),secondMuon->phi()) < 0.09 )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if ( reco::deltaR2(firstMuon->eta(),firstMuon->phi(),secondMuon->eta(),secondMuon->phi()) < 0.09 )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      if ( reco::deltaR2(firstMuon->innerTrack()->eta(),firstMuon->innerTrack()->phi(),
                         isoCands.front()->innerTrack()->eta(),isoCands.front()->innerTrack()->phi()) > 0.0001 )
        firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.05 )
      isolatedHighPtMuons.push_back(firstMuon);
  }

  for (unsigned idx = 0; idx < highPtTrackerMuons.size(); idx++) {
    const auto& firstMuon = highPtTrackerMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( reco::deltaR2(firstMuon->eta(),firstMuon->phi(),secondMuon->eta(),secondMuon->phi()) < 0.09 )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if ( reco::deltaR2(firstMuon->eta(),firstMuon->phi(),secondMuon->eta(),secondMuon->phi()) < 0.09 )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      if ( reco::deltaR2(firstMuon->innerTrack()->eta(),firstMuon->innerTrack()->phi(),
                         isoCands.front()->innerTrack()->eta(),isoCands.front()->innerTrack()->phi()) > 0.0001 )
        firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.05 )
      isolatedHighPtTrackerMuons.push_back(firstMuon);
  }

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  if ( isolatedHighPtMuons.size() < 2 )
    return;

  if ( isolatedHighPtMuons.front()->tunePMuonBestTrack()->pt() < ptThresTrig_ )
    return;

  histo1d_["cutflow"]->Fill( 3.5, aWeight );

  bool trigMatched = false;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),isolatedHighPtMuons.front()->eta(),isolatedHighPtMuons.front()->phi()) < 0.01 ) {
      trigMatched = true;

      if (isMC_)
        aWeight *= mucorrHelper_.triggerSF(isolatedHighPtMuons.front());

      break;
    }
  }

  if ( !trigMatched )
    return;

  histo1d_["cutflow"]->Fill( 4.5, aWeight );

  // concatenate(-ish) highPtMuons & highPtTrackerMuons - order is important!
  std::vector<pat::MuonRef> allHighPtMuons(isolatedHighPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), isolatedHighPtTrackerMuons.begin(), isolatedHighPtTrackerMuons.end() );
  std::sort(allHighPtMuons.begin(),allHighPtMuons.end(),sortByTuneP);

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( !aMuon->isTrackerMuon() || aMuon->track().isNull() )
      continue;

    if ( aMuon->tunePMuonBestTrack()->pt() < ptThres_ || std::abs(aMuon->eta()) > 2.4 )
      continue;

    const auto castMu = aMuon.castTo<pat::MuonRef>();

    if ( std::find( allHighPtMuons.begin(), allHighPtMuons.end(), castMu ) != allHighPtMuons.end() )
      continue;

    if ( castMu->numberOfMatchedStations() >= 1 )
      nonHighPtMuons.push_back(castMu);
  }

  auto checkTrackerHighPtPair = [&isolatedHighPtTrackerMuons] (const std::pair<pat::MuonRef,pat::MuonRef>& apair) -> bool {
    bool check1st = std::find(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),apair.first) != isolatedHighPtTrackerMuons.end();
    bool check2nd = std::find(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),apair.second) != isolatedHighPtTrackerMuons.end();

    return check1st && check2nd;
  };

  switch (allHighPtMuons.size()) {
    case 4:
      // 4P0F
      if ( nonHighPtMuons.size()==0 ) {
        std::vector<pat::MuonRef> allMuons(allHighPtMuons);
        std::pair<pat::MuonRef,pat::MuonRef> pair1, pair2;

        bool paired = PairingHelper::pair4M(allHighPtMuons,pair1,pair2);

        if (paired && !checkTrackerHighPtPair(pair1) && !checkTrackerHighPtPair(pair2)) {
          auto lvecM1 = lvecFromTuneP(pair1.first);
          auto lvecM2 = lvecFromTuneP(pair1.second);
          auto lvecM3 = lvecFromTuneP(pair2.first);
          auto lvecM4 = lvecFromTuneP(pair2.second);

          if (isMC_) {
            edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
            iEvent.getByToken(genptcToken_, genptcHandle);

            lvecM1 *= mucorrHelper_.rochesterMC(pair1.first,genptcHandle);
            lvecM2 *= mucorrHelper_.rochesterMC(pair1.second,genptcHandle);
            lvecM3 *= mucorrHelper_.rochesterMC(pair2.first,genptcHandle);
            lvecM4 *= mucorrHelper_.rochesterMC(pair2.second,genptcHandle);
          } else {
            lvecM1 *= mucorrHelper_.rochesterData(pair1.first);
            lvecM2 *= mucorrHelper_.rochesterData(pair1.second);
            lvecM3 *= mucorrHelper_.rochesterData(pair2.first);
            lvecM4 *= mucorrHelper_.rochesterData(pair2.second);
          }

          const auto lvecA1 = lvecM1 + lvecM2;
          const auto lvecA2 = lvecM3 + lvecM4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
          const double dr2ll2 = reco::deltaR2(lvecM3.eta(),lvecM3.phi(),lvecM4.eta(),lvecM4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();

          if ( m4l > 50. && (m4l < 200. || isMC_) && lvecA1.M() > 1. && lvecA2.M() > 1. )
            histo1d_["4P0F_CR_llll_invM"]->Fill(m4l, aWeight);

          if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["4P0F_CR_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["4P0F_CR_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
            histo1d_["4P0F_CR_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["4P0F_CR_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);
            histo1d_["4P0F_CR_P3_eta"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["4P0F_CR_P3_phi"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->phi(), aWeight);
            histo1d_["4P0F_CR_P4_eta"]->Fill(allHighPtMuons.at(3)->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["4P0F_CR_P4_phi"]->Fill(allHighPtMuons.at(3)->tunePMuonBestTrack()->phi(), aWeight);

            histo1d_["4P0F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
            histo1d_["4P0F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
            histo1d_["4P0F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
            histo1d_["4P0F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
            histo1d_["4P0F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
            histo1d_["4P0F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
            histo1d_["4P0F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
            histo1d_["4P0F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

            histo1d_["4P0F_CR_pt_P1"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["4P0F_CR_pt_P2"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["4P0F_CR_pt_P3"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["4P0F_CR_pt_P4"]->Fill(allHighPtMuons.at(3)->tunePMuonBestTrack()->pt(), aWeight);
          } // CR
        } // paired
      } // 4P0F

      break;
    case 3:
      // 3P0F - FF numerator
      if ( nonHighPtMuons.size()==0 ) {
        if ( allHighPtMuons.at(2)->tunePMuonBestTrack()->pt() < ptThresTrig_ ) {
          auto& first = allHighPtMuons.front();
          auto& second = allHighPtMuons.at(1);
          auto& probe = allHighPtMuons.at(2);

          auto lvecM1 = lvecFromTuneP(first);
          auto lvecM2 = lvecFromTuneP(second);

          if (isMC_) {
            edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
            iEvent.getByToken(genptcToken_, genptcHandle);

            lvecM1 *= mucorrHelper_.rochesterMC(first,genptcHandle);
            lvecM2 *= mucorrHelper_.rochesterMC(second,genptcHandle);
          } else {
            lvecM1 *= mucorrHelper_.rochesterData(first);
            lvecM2 *= mucorrHelper_.rochesterData(second);
          }

          const auto lvecll = lvecM1 + lvecM2;
          const double mll = lvecll.M();
          const double dr2ll = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());

          if ( mll > 84.19 && mll < 98.19 && first->charge()*second->charge() < 0 ) {
            histo1d_["3P0F_P1_pt"]->Fill(first->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["3P0F_P1_eta"]->Fill(first->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["3P0F_P1_phi"]->Fill(first->tunePMuonBestTrack()->phi(), aWeight);
            histo1d_["3P0F_P2_pt"]->Fill(second->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["3P0F_P2_eta"]->Fill(second->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["3P0F_P2_phi"]->Fill(second->tunePMuonBestTrack()->phi(), aWeight);

            histo1d_["3P0F_P3_pt"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["3P0F_P3_eta"]->Fill(probe->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["3P0F_P3_phi"]->Fill(probe->tunePMuonBestTrack()->phi(), aWeight);

            histo1d_["3P0F_ll_invM"]->Fill(lvecll.M(), aWeight);
            histo1d_["3P0F_ll_pt"]->Fill(lvecll.pt(), aWeight);
            histo1d_["3P0F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

            const double mindr2 = std::min( reco::deltaR2(probe->tunePMuonBestTrack()->eta(),probe->tunePMuonBestTrack()->phi(),
                                                          first->tunePMuonBestTrack()->eta(),first->tunePMuonBestTrack()->phi()),
                                            reco::deltaR2(probe->tunePMuonBestTrack()->eta(),probe->tunePMuonBestTrack()->phi(),
                                                          second->tunePMuonBestTrack()->eta(),second->tunePMuonBestTrack()->phi()) );

            histo1d_["3P0F_dr"]->Fill( std::sqrt(mindr2), aWeight );

            if ( std::abs( probe->tunePMuonBestTrack()->eta() ) < 1.4 ) {
              histo1d_["3P0F_pt_MB"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_dr_MB"]->Fill( std::sqrt(mindr2), aWeight );
            } else {
              histo1d_["3P0F_pt_ME"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_dr_ME"]->Fill( std::sqrt(mindr2), aWeight );
            }
          } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
        } else {
          for (unsigned idx = 0; idx < allHighPtMuons.size(); idx++) {
            const auto& first = allHighPtMuons.at(idx);
            const auto& second = allHighPtMuons.at( (idx+1)%allHighPtMuons.size() );
            const auto& probe = allHighPtMuons.at( (idx+2)%allHighPtMuons.size() );

            auto lvecM1 = lvecFromTuneP(first);
            auto lvecM2 = lvecFromTuneP(second);

            if (isMC_) {
              edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
              iEvent.getByToken(genptcToken_, genptcHandle);

              lvecM1 *= mucorrHelper_.rochesterMC(first,genptcHandle);
              lvecM2 *= mucorrHelper_.rochesterMC(second,genptcHandle);
            } else {
              lvecM1 *= mucorrHelper_.rochesterData(first);
              lvecM2 *= mucorrHelper_.rochesterData(second);
            }

            const auto lvecll = lvecM1 + lvecM2;
            const double mll = lvecll.M();
            const double dr2ll = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());

            if ( mll > 84.19 && mll < 98.19 && first->charge()*second->charge() < 0 ) {
              histo1d_["3P0F_P1_pt"]->Fill(first->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_P1_eta"]->Fill(first->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["3P0F_P1_phi"]->Fill(first->tunePMuonBestTrack()->phi(), aWeight);
              histo1d_["3P0F_P2_pt"]->Fill(second->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_P2_eta"]->Fill(second->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["3P0F_P2_phi"]->Fill(second->tunePMuonBestTrack()->phi(), aWeight);

              histo1d_["3P0F_P3_pt"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_P3_eta"]->Fill(probe->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["3P0F_P3_phi"]->Fill(probe->tunePMuonBestTrack()->phi(), aWeight);

              histo1d_["3P0F_ll_invM"]->Fill(lvecll.M(), aWeight);
              histo1d_["3P0F_ll_pt"]->Fill(lvecll.pt(), aWeight);
              histo1d_["3P0F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

              const double mindr2 = std::min( reco::deltaR2(probe->tunePMuonBestTrack()->eta(),probe->tunePMuonBestTrack()->phi(),
                                                            first->tunePMuonBestTrack()->eta(),first->tunePMuonBestTrack()->phi()),
                                              reco::deltaR2(probe->tunePMuonBestTrack()->eta(),probe->tunePMuonBestTrack()->phi(),
                                                            second->tunePMuonBestTrack()->eta(),second->tunePMuonBestTrack()->phi()) );

              histo1d_["3P0F_dr"]->Fill( std::sqrt(mindr2), aWeight );

              if ( std::abs( probe->tunePMuonBestTrack()->eta() ) < 1.4 ) {
                histo1d_["3P0F_pt_MB"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
                histo1d_["3P0F_dr_MB"]->Fill( std::sqrt(mindr2), aWeight );
              } else {
                histo1d_["3P0F_pt_ME"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
                histo1d_["3P0F_dr_ME"]->Fill( std::sqrt(mindr2), aWeight );
              }
            } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
          } // Z candidate combination
        }
      } // 3P0F

      // 3P1F
      if ( nonHighPtMuons.size()==1 ) {
        histo1d_["3P1F_P1_pt"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
        histo1d_["3P1F_P2_pt"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);
        histo1d_["3P1F_P3_pt"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_P3_eta"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_P3_phi"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->phi(), aWeight);

        histo1d_["3P1F_F1_pt"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_F1_eta"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_F1_phi"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);

        std::vector<pat::MuonRef> allMuons(allHighPtMuons);
        allMuons.insert( allMuons.end(), nonHighPtMuons.begin(), nonHighPtMuons.end() );
        std::pair<pat::MuonRef,pat::MuonRef> pair1, pair2;

        bool paired = PairingHelper::pair4M(allMuons,pair1,pair2);

        if (paired && !checkTrackerHighPtPair(pair1) && !checkTrackerHighPtPair(pair2)) {
          auto lvecM1 = lvecFromTuneP(pair1.first);
          auto lvecM2 = lvecFromTuneP(pair1.second);
          auto lvecM3 = lvecFromTuneP(pair2.first);
          auto lvecM4 = lvecFromTuneP(pair2.second);

          if (isMC_) {
            edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
            iEvent.getByToken(genptcToken_, genptcHandle);

            lvecM1 *= mucorrHelper_.rochesterMC(pair1.first,genptcHandle);
            lvecM2 *= mucorrHelper_.rochesterMC(pair1.second,genptcHandle);
            lvecM3 *= mucorrHelper_.rochesterMC(pair2.first,genptcHandle);
            lvecM4 *= mucorrHelper_.rochesterMC(pair2.second,genptcHandle);
          } else {
            lvecM1 *= mucorrHelper_.rochesterData(pair1.first);
            lvecM2 *= mucorrHelper_.rochesterData(pair1.second);
            lvecM3 *= mucorrHelper_.rochesterData(pair2.first);
            lvecM4 *= mucorrHelper_.rochesterData(pair2.second);
          }

          const auto lvecA1 = lvecM1 + lvecM2;
          const auto lvecA2 = lvecM3 + lvecM4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
          const double dr2ll2 = reco::deltaR2(lvecM3.eta(),lvecM3.phi(),lvecM4.eta(),lvecM4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();

          const double drFake = std::sqrt( std::min({ reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                    allHighPtMuons.at(2)->tunePMuonBestTrack()->eta(),allHighPtMuons.at(2)->tunePMuonBestTrack()->phi()),
                                                      reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                    allHighPtMuons.front()->tunePMuonBestTrack()->eta(),allHighPtMuons.front()->tunePMuonBestTrack()->phi()),
                                                      reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                    allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),allHighPtMuons.at(1)->tunePMuonBestTrack()->phi()) }) );
          const double ff = drFF_->Eval(drFake);
          const double xval[1] = {drFake};
          double ci[1];
          ffFit_->GetConfidenceIntervals(1,1,0,xval,ci,0.6827,false);

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
            histo1d_["3P1F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff);
          }

          if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            const auto& fakePair = ( nonHighPtMuons.front()==pair1.first || nonHighPtMuons.front()==pair1.second ) ? pair1 : pair2;
            const auto& passPartner = ( nonHighPtMuons.front()==fakePair.first ) ? fakePair.second : fakePair.first;

            histo1d_["3P1F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
            histo1d_["3P1F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
            histo1d_["3P1F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
            histo1d_["3P1F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
            histo1d_["3P1F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
            histo1d_["3P1F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
            histo1d_["3P1F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
            histo1d_["3P1F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

            histo1d_["3P1F_CR_pt_P1"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["3P1F_CR_pt_P2"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["3P1F_CR_pt_P3"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["3P1F_CR_pt_F1"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);

            const double drPass = std::sqrt( reco::deltaR2(passPartner->tunePMuonBestTrack()->eta(),passPartner->tunePMuonBestTrack()->phi(),
                                                           nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi()) );

            histo1d_["3P1F_CR_dr"]->Fill(drPass, aWeight);

            if ( std::abs( passPartner->tunePMuonBestTrack()->eta() ) < 1.4 ) {
              histo1d_["3P1F_CR_pt_MB"]->Fill(passPartner->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P1F_CR_dr_MB"]->Fill(drPass, aWeight);
            } else {
              histo1d_["3P1F_CR_pt_ME"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P1F_CR_dr_ME"]->Fill(drPass, aWeight);
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
      if ( nonHighPtMuons.size()==1 && allHighPtMuons.front()->charge()*allHighPtMuons.at(1)->charge() < 0 ) {
        auto lvecM1 = lvecFromTuneP(allHighPtMuons.front());
        auto lvecM2 = lvecFromTuneP(allHighPtMuons.at(1));

        if (isMC_) {
          edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
          iEvent.getByToken(genptcToken_, genptcHandle);

          lvecM1 *= mucorrHelper_.rochesterMC(allHighPtMuons.front(),genptcHandle);
          lvecM2 *= mucorrHelper_.rochesterMC(allHighPtMuons.at(1),genptcHandle);
        } else {
          lvecM1 *= mucorrHelper_.rochesterData(allHighPtMuons.front());
          lvecM2 *= mucorrHelper_.rochesterData(allHighPtMuons.at(1));
        }

        const auto lvecll = lvecM1 + lvecM2;
        const double mll = lvecll.M();
        const double dr2ll = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());

        if ( mll > 84.19 && mll < 98.19 && allHighPtMuons.front()->charge()*allHighPtMuons.at(1)->charge() < 0 ) {
          histo1d_["2P1F_P1_pt"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["2P1F_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["2P1F_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
          histo1d_["2P1F_P2_pt"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["2P1F_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["2P1F_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);

          histo1d_["2P1F_F1_pt"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["2P1F_F1_eta"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["2P1F_F1_phi"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);

          histo1d_["2P1F_ll_invM"]->Fill(mll, aWeight);
          histo1d_["2P1F_ll_pt"]->Fill(lvecll.pt(), aWeight);
          histo1d_["2P1F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

          const double mindr2 = std::min( reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                        allHighPtMuons.front()->tunePMuonBestTrack()->eta(),allHighPtMuons.front()->tunePMuonBestTrack()->phi()),
                                          reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                        allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),allHighPtMuons.at(1)->tunePMuonBestTrack()->phi()) );

          histo1d_["2P1F_dr"]->Fill(std::sqrt(mindr2), aWeight);

          if ( std::abs( nonHighPtMuons.front()->tunePMuonBestTrack()->eta() ) < 1.4 ) {
            histo1d_["2P1F_pt_MB"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P1F_dr_MB"]->Fill(std::sqrt(mindr2), aWeight);
          } else {
            histo1d_["2P1F_pt_ME"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P1F_dr_ME"]->Fill(std::sqrt(mindr2), aWeight);
          }
        } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
      } // 2P1F

      // 2P2F
      if ( nonHighPtMuons.size()==2 ) {
        std::vector<pat::MuonRef> allMuons(allHighPtMuons);
        allMuons.insert( allMuons.end(), nonHighPtMuons.begin(), nonHighPtMuons.end() );
        std::pair<pat::MuonRef,pat::MuonRef> pair1, pair2;

        bool paired = PairingHelper::pair4M(allMuons,pair1,pair2);

        if (paired && !checkTrackerHighPtPair(pair1) && !checkTrackerHighPtPair(pair2)) {
          auto lvecM1 = lvecFromTuneP(pair1.first);
          auto lvecM2 = lvecFromTuneP(pair1.second);
          auto lvecM3 = lvecFromTuneP(pair2.first);
          auto lvecM4 = lvecFromTuneP(pair2.second);

          if (isMC_) {
            edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
            iEvent.getByToken(genptcToken_, genptcHandle);

            lvecM1 *= mucorrHelper_.rochesterMC(pair1.first,genptcHandle);
            lvecM2 *= mucorrHelper_.rochesterMC(pair1.second,genptcHandle);
            lvecM3 *= mucorrHelper_.rochesterMC(pair2.first,genptcHandle);
            lvecM4 *= mucorrHelper_.rochesterMC(pair2.second,genptcHandle);
          } else {
            lvecM1 *= mucorrHelper_.rochesterData(pair1.first);
            lvecM2 *= mucorrHelper_.rochesterData(pair1.second);
            lvecM3 *= mucorrHelper_.rochesterData(pair2.first);
            lvecM4 *= mucorrHelper_.rochesterData(pair2.second);
          }

          const auto lvecA1 = lvecM1 + lvecM2;
          const auto lvecA2 = lvecM3 + lvecM4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
          const double dr2ll2 = reco::deltaR2(lvecM3.eta(),lvecM3.phi(),lvecM4.eta(),lvecM4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();

          const double drFake1 = std::sqrt( std::min({ reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                     nonHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),nonHighPtMuons.at(1)->tunePMuonBestTrack()->phi()),
                                                       reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                     allHighPtMuons.front()->tunePMuonBestTrack()->eta(),allHighPtMuons.front()->tunePMuonBestTrack()->phi()),
                                                       reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                     allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),allHighPtMuons.at(1)->tunePMuonBestTrack()->phi()) }) );
          const double drFake2 = std::sqrt( std::min({ reco::deltaR2(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(),nonHighPtMuons.front()->tunePMuonBestTrack()->phi(),
                                                                     nonHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),nonHighPtMuons.at(1)->tunePMuonBestTrack()->phi()),
                                                       reco::deltaR2(nonHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),nonHighPtMuons.at(1)->tunePMuonBestTrack()->phi(),
                                                                     allHighPtMuons.front()->tunePMuonBestTrack()->eta(),allHighPtMuons.front()->tunePMuonBestTrack()->phi()),
                                                       reco::deltaR2(nonHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),nonHighPtMuons.at(1)->tunePMuonBestTrack()->phi(),
                                                                     allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),allHighPtMuons.at(1)->tunePMuonBestTrack()->phi()) }) );

          const double ff1 = drFF_->Eval(drFake1);
          const double xval1[1] = {drFake1};
          double ci1[1];
          ffFit_->GetConfidenceIntervals(1,1,0,xval1,ci1,0.6827,false);

          const double ff2 = drFF_->Eval(drFake2);
          const double xval2[1] = {drFake2};
          double ci2[1];
          ffFit_->GetConfidenceIntervals(1,1,0,xval2,ci2,0.6827,false);

          if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["2P2F_P1_pt"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["2P2F_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
            histo1d_["2P2F_P2_pt"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["2P2F_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);

            histo1d_["2P2F_F1_pt"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_F1_eta"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["2P2F_F1_phi"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
            histo1d_["2P2F_F2_pt"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_F2_eta"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
            histo1d_["2P2F_F2_phi"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);

            histo1d_["2P2F_llll_pt"]->Fill(lvec4l.pt(), aWeight);
            histo1d_["2P2F_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
            histo1d_["2P2F_ll1_invM"]->Fill(lvecA1.M(), aWeight);
            histo1d_["2P2F_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
            histo1d_["2P2F_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
            histo1d_["2P2F_ll2_invM"]->Fill(lvecA2.M(), aWeight);
            histo1d_["2P2F_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
            histo1d_["2P2F_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

            histo1d_["2P2F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["2P2F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff1);
            histo1d_["2P2F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff2);
            histo1d_["2P2F_CR_llll_invM_xFF2"]->Fill(m4l, aWeight*ff1*ff2);
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

            if ( std::abs( nonHighPtMuons.front()->tunePMuonBestTrack()->eta() ) < 1.4 ) {
              histo1d_["2P2F_CR_pt_MB"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_dr_MB"]->Fill(drFake1, aWeight);
              histo1d_["2P2F_CR_pt_MB_xFF"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight*ff1);
              histo1d_["2P2F_CR_dr_MB_xFF"]->Fill(drFake1, aWeight*ff1);
            } else {
              histo1d_["2P2F_CR_pt_ME"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_dr_ME"]->Fill(drFake1, aWeight);
              histo1d_["2P2F_CR_pt_ME_xFF"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight*ff1);
              histo1d_["2P2F_CR_dr_ME_xFF"]->Fill(drFake1, aWeight*ff1);
            }

            if ( std::abs( nonHighPtMuons.at(1)->tunePMuonBestTrack()->eta() ) < 1.4 ) {
              histo1d_["2P2F_CR_pt_MB"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_dr_MB"]->Fill(drFake2, aWeight);
              histo1d_["2P2F_CR_pt_MB_xFF"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight*ff2);
              histo1d_["2P2F_CR_dr_MB_xFF"]->Fill(drFake2, aWeight*ff2);
            } else {
              histo1d_["2P2F_CR_pt_ME"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_dr_ME"]->Fill(drFake2, aWeight);
              histo1d_["2P2F_CR_pt_ME_xFF"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight*ff2);
              histo1d_["2P2F_CR_dr_ME_xFF"]->Fill(drFake2, aWeight*ff2);
            }

            histo1d_["2P2F_CR_pt_P1"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_CR_pt_P2"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_CR_pt_F1"]->Fill(nonHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P2F_CR_pt_F2"]->Fill(nonHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
          } // CR
        } // paired
      } // 2P2F

      break;
  } // switch allHighPtMuons.size()
}

DEFINE_FWK_MODULE(ResolvedMuCRanalyzer);
