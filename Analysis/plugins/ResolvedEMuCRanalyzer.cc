#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Common/interface/RefToPtr.h"
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
#include "ZprimeTo4l/Analysis/interface/ElectronSystematicsHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"

class ResolvedEMuCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedEMuCRanalyzer(const edm::ParameterSet&);
  virtual ~ResolvedEMuCRanalyzer() {}

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

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const std::string trigHistName_;
  const std::string idHistName_;
  const std::string idHistNameTrkHighPt_;
  const std::string isoHistName_;
  const std::string isoHistNameTrkHighPt_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath muonIdSFpath_;
  const edm::FileInPath muonIsoSFpath_;
  const edm::FileInPath purwgtPath_;
  const edm::FileInPath recoEleSFpath_;

  const edm::FileInPath RMFFpath_;
  const edm::FileInPath REFFpath_;

  const double ptThres_;
  const double ptThresTrig_;

  const double ffSystCL95_;

  MuonCorrectionHelper mucorrHelper_;
  ElectronSystematicsHelper systHelperEle_;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::unique_ptr<TFile> recoEleSFfile_;
  TH2D* recoEleSF_;

  std::unique_ptr<TFile> RMFFfile_;
  TF1* drRMFF_;
  TFitResultPtr ffMuonFit_;

  std::unique_ptr<TFile> REFFfile_;
  TF1* drREFF_;
  TFitResultPtr ffEleFit_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  const double mumass_ = 0.1056583745;
};

ResolvedEMuCRanalyzer::ResolvedEMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
trigHistName_(iConfig.getParameter<std::string>("trigHistName")),
idHistName_(iConfig.getParameter<std::string>("idHistName")),
idHistNameTrkHighPt_(iConfig.getParameter<std::string>("idHistNameTrkHighPt")),
isoHistName_(iConfig.getParameter<std::string>("isoHistName")),
isoHistNameTrkHighPt_(iConfig.getParameter<std::string>("isoHistNameTrkHighPt")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
muonIdSFpath_(iConfig.getParameter<edm::FileInPath>("muonIdSFpath")),
muonIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonIsoSFpath")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
recoEleSFpath_(iConfig.getParameter<edm::FileInPath>("recoEleSFpath")),
RMFFpath_(iConfig.getParameter<edm::FileInPath>("RMFFpath")),
REFFpath_(iConfig.getParameter<edm::FileInPath>("REFFpath")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptThresTrig_(iConfig.getParameter<double>("ptThresTrig")),
ffSystCL95_(iConfig.getParameter<double>("ffSystCL95")),
mucorrHelper_(rochesterPath_,triggerSFpath_,muonIdSFpath_,muonIsoSFpath_,trigHistName_),
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

math::PtEtaPhiMLorentzVector ResolvedEMuCRanalyzer::lvecFromTuneP(const pat::MuonRef& aMu) {
  return math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                      aMu->tunePMuonBestTrack()->eta(),
                                      aMu->tunePMuonBestTrack()->phi(),
                                      mumass_);
}

void ResolvedEMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  recoEleSFfile_ = std::make_unique<TFile>(recoEleSFpath_.fullPath().c_str(),"READ");
  recoEleSF_ = static_cast<TH2D*>(recoEleSFfile_->Get("EGamma_SF2D"));

  RMFFfile_ = std::make_unique<TFile>(RMFFpath_.fullPath().c_str(),"READ");
  drRMFF_ = static_cast<TF1*>(RMFFfile_->FindObjectAny("RMFF_dr_all"));
  ffMuonFit_ = (static_cast<TH1D*>(RMFFfile_->Get("all_dr_3P0F_rebin")))->Fit(drRMFF_,"RS");

  REFFfile_ = std::make_unique<TFile>(REFFpath_.fullPath().c_str(),"READ");
  drREFF_ = static_cast<TF1*>(REFFfile_->FindObjectAny("REFF_dr_all"));
  ffEleFit_ = (static_cast<TH1D*>(REFFfile_->Get("all_dr_3P0F_rebin")))->Fit(drREFF_,"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",10,0.,10.);

  // 2P2F
  histo1d_["2P2F_EP_pt"] = fs->make<TH1D>("2P2F_EP_pt","2P2F Pass E p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_EP_eta"] = fs->make<TH1D>("2P2F_EP_eta","2P2F Pass E #eta",200,-2.5,2.5);
  histo1d_["2P2F_EP_phi"] = fs->make<TH1D>("2P2F_EP_phi","2P2F Pass E #phi",128,-3.2,3.2);

  histo1d_["2P2F_MP_pt"] = fs->make<TH1D>("2P2F_MP_pt","2P2F Pass M p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_MP_eta"] = fs->make<TH1D>("2P2F_MP_eta","2P2F Pass M #eta",200,-2.5,2.5);
  histo1d_["2P2F_MP_phi"] = fs->make<TH1D>("2P2F_MP_phi","2P2F Pass M #phi",128,-3.2,3.2);

  histo1d_["2P2F_EF_pt"] = fs->make<TH1D>("2P2F_EF_pt","2P2F Fail E p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_EF_eta"] = fs->make<TH1D>("2P2F_EF_eta","2P2F Fail E #eta",200,-2.5,2.5);
  histo1d_["2P2F_EF_phi"] = fs->make<TH1D>("2P2F_EF_phi","2P2F Fail E #phi",128,-3.2,3.2);

  histo1d_["2P2F_MF_pt"] = fs->make<TH1D>("2P2F_MF_pt","2P2F Fail M p_{T};p_{T} [GeV];",200,0.,500.);
  histo1d_["2P2F_MF_eta"] = fs->make<TH1D>("2P2F_MF_eta","2P2F Fail M #eta",200,-2.5,2.5);
  histo1d_["2P2F_MF_phi"] = fs->make<TH1D>("2P2F_MF_phi","2P2F Fail M #phi",128,-3.2,3.2);

  histo1d_["2P2F_llll_pt"] = fs->make<TH1D>("2P2F_llll_pt","2P2F p_{T}(4l);p_{T} [GeV];",100,0.,1000.);
  histo1d_["2P2F_ll1ll2_dr"] = fs->make<TH1D>("2P2F_ll1ll2_dr","2P2F dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_ll1_invM"] = fs->make<TH1D>("2P2F_ll1_invM","2P2F M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_ll1_pt"] = fs->make<TH1D>("2P2F_ll1_pt","2P2F p_{T}(ll1);p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_ll1_dr"] = fs->make<TH1D>("2P2F_ll1_dr","2P2F dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_ll2_invM"] = fs->make<TH1D>("2P2F_ll2_invM","2P2F M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_ll2_pt"] = fs->make<TH1D>("2P2F_ll2_pt","2P2F p_{T}(ll2);p_{T} [GeV];",100,0.,500.);
  histo1d_["2P2F_ll2_dr"] = fs->make<TH1D>("2P2F_ll2_dr","2P2F dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CR_pt_E"] = fs->make<TH1D>("2P2F_CR_pt_E","2P2F CR Fail p_{T} (E);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_pt_M"] = fs->make<TH1D>("2P2F_CR_pt_M","2P2F CR Fail p_{T} (M);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_dr_E"] = fs->make<TH1D>("2P2F_CR_dr_E","2P2F CR Fail min dR(ll) (E)",128,0.,6.4);
  histo1d_["2P2F_CR_dr_M"] = fs->make<TH1D>("2P2F_CR_dr_M","2P2F CR Fail min dR(ll) (M)",128,0.,6.4);

  histo1d_["2P2F_CR_pt_E_xFF"] = fs->make<TH1D>("2P2F_CR_pt_E_xFF","2P2F CR Fail p_{T} x Fake Factor (E);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_pt_M_xFF"] = fs->make<TH1D>("2P2F_CR_pt_M_xFF","2P2F CR Fail p_{T} x Fake Factor (M);p_{T} [GeV];",200,0.,1000.);
  histo1d_["2P2F_CR_dr_E_xFF"] = fs->make<TH1D>("2P2F_CR_dr_E_xFF","2P2F CR Fail min dR(ll) x Fake Factor (E)",128,0.,6.4);
  histo1d_["2P2F_CR_dr_M_xFF"] = fs->make<TH1D>("2P2F_CR_dr_M_xFF","2P2F CR Fail min dR(ll) x Fake Factor (M)",128,0.,6.4);

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
  histo1d_["2P2F_CR_llll_invM_xFF_ffUpE"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffUpE","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffDnE"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffDnE","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffUpM"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffUpM","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffDnM"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffDnM","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_idUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_idUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_idDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_idDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF","2P2F CR dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF","2P2F CR M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF","2P2F CR dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF","2P2F CR M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF","2P2F CR dR(ll2) x Fake factor",128,0.,6.4);

  histo1d_["2P2F_CR_llll_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffUpE"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffUpE","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffDnE"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffDnE","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffUpM"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffUpM","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffDnM"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffDnM","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_idUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_idUp","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_idDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_idDn","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF2","2P2F CR dR(ll1ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF2","2P2F CR M(ll1) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF2","2P2F CR dR(ll1) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF2","2P2F CR M(ll2) x Fake factor^2;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF2","2P2F CR dR(ll2) x Fake factor^2",128,0.,6.4);

  // 3P1F
  histo1d_["3P1F_EP_pt"] = fs->make<TH1D>("3P1F_EP_pt","3P1F Pass E p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_EP_eta"] = fs->make<TH1D>("3P1F_EP_eta","3P1F Pass E #eta",100,-2.5,2.5);
  histo1d_["3P1F_EP_phi"] = fs->make<TH1D>("3P1F_EP_phi","3P1F Pass E #phi",64,-3.2,3.2);

  histo1d_["3P1F_MP_pt"] = fs->make<TH1D>("3P1F_MP_pt","3P1F Pass M p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_MP_eta"] = fs->make<TH1D>("3P1F_MP_eta","3P1F Pass M #eta",100,-2.5,2.5);
  histo1d_["3P1F_MP_phi"] = fs->make<TH1D>("3P1F_MP_phi","3P1F Pass M #phi",64,-3.2,3.2);

  histo1d_["3P1F_EF_pt"] = fs->make<TH1D>("3P1F_EF_pt","3P1F Fail E p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_EF_eta"] = fs->make<TH1D>("3P1F_EF_eta","3P1F Fail E #eta",100,-2.5,2.5);
  histo1d_["3P1F_EF_phi"] = fs->make<TH1D>("3P1F_EF_phi","3P1F Fail E #phi",64,-3.2,3.2);

  histo1d_["3P1F_MF_pt"] = fs->make<TH1D>("3P1F_MF_pt","3P1F Fail M p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_MF_eta"] = fs->make<TH1D>("3P1F_MF_eta","3P1F Fail M #eta",100,-2.5,2.5);
  histo1d_["3P1F_MF_phi"] = fs->make<TH1D>("3P1F_MF_phi","3P1F Fail M #phi",64,-3.2,3.2);

  histo1d_["3P1F_llll_pt"] = fs->make<TH1D>("3P1F_llll_pt","3P1F p_{T}(4l);p_{T} [GeV];",100,0.,500.);
  histo1d_["3P1F_ll1ll2_dr"] = fs->make<TH1D>("3P1F_ll1ll2_dr","3P1F dR(ll1ll2)",64,0.,6.4);
  histo1d_["3P1F_ll1_invM"] = fs->make<TH1D>("3P1F_ll1_invM","3P1F M(ll1);M [GeV];",200,0.,200.);
  histo1d_["3P1F_ll1_pt"] = fs->make<TH1D>("3P1F_ll1_pt","3P1F p_{T}(ll1);p_{T} [GeV];",200,0.,500.);
  histo1d_["3P1F_ll1_dr"] = fs->make<TH1D>("3P1F_ll1_dr","3P1F dR(ll1)",64,0.,6.4);
  histo1d_["3P1F_ll2_invM"] = fs->make<TH1D>("3P1F_ll2_invM","3P1F M(ll2);M [GeV];",200,0.,200.);
  histo1d_["3P1F_ll2_pt"] = fs->make<TH1D>("3P1F_ll2_pt","3P1F p_{T}(ll2);p_{T} [GeV];",200,0.,500.);
  histo1d_["3P1F_ll2_dr"] = fs->make<TH1D>("3P1F_ll2_dr","3P1F dR(ll2)",64,0.,6.4);

  histo1d_["3P1F_CR_pt_E"] = fs->make<TH1D>("3P1F_CR_pt_E","3P1F CR Pass p_{T} (E);p_{T} [GeV];",200,0.,1000.);
  histo1d_["3P1F_CR_pt_M"] = fs->make<TH1D>("3P1F_CR_pt_M","3P1F CR Pass p_{T} (M);p_{T} [GeV];",200,0.,1000.);
  histo1d_["3P1F_CR_dr_E"] = fs->make<TH1D>("3P1F_CR_dr_E","3P1F CR Pass min dR(ll) (E)",128,0.,6.4);
  histo1d_["3P1F_CR_dr_M"] = fs->make<TH1D>("3P1F_CR_dr_M","3P1F CR Pass min dR(ll) (M)",128,0.,6.4);

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
  histo1d_["3P1F_CR_llll_invM_xFF_ffUpE"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffUpE","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffDnE"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffDnE","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffUpM"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffUpM","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffDnM"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffDnM","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_idUp"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_idUp","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_idDn"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_idDn","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll1ll2_dr_xFF","3P1F CR dR(ll1ll2) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_ll1_invM_xFF"] = fs->make<TH1D>("3P1F_CR_ll1_invM_xFF","3P1F CR M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll1_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll1_dr_xFF","3P1F CR dR(ll1) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_ll2_invM_xFF"] = fs->make<TH1D>("3P1F_CR_ll2_invM_xFF","3P1F CR M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["3P1F_CR_ll2_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll2_dr_xFF","3P1F CR dR(ll2) x Fake factor",64,0.,6.4);

  // 4P0F
  histo1d_["4P0F_CR_E_eta"] = fs->make<TH1D>("4P0F_CR_E_eta","4P0F CR Pass E #eta",100,-2.5,2.5);
  histo1d_["4P0F_CR_E_phi"] = fs->make<TH1D>("4P0F_CR_E_phi","4P0F CR Pass E #phi",64,-3.2,3.2);

  histo1d_["4P0F_CR_M_eta"] = fs->make<TH1D>("4P0F_CR_M_eta","4P0F CR Pass M #eta",100,-2.5,2.5);
  histo1d_["4P0F_CR_M_phi"] = fs->make<TH1D>("4P0F_CR_M_phi","4P0F CR Pass M #phi",64,-3.2,3.2);

  histo1d_["4P0F_CR_pt_E"] = fs->make<TH1D>("4P0F_CR_pt_E","4P0F CR Pass E p_{T};p_{T} [GeV];",100,0.,500.);
  histo1d_["4P0F_CR_pt_M"] = fs->make<TH1D>("4P0F_CR_pt_M","4P0F CR Pass M p_{T};p_{T} [GeV];",100,0.,500.);

  histo1d_["4P0F_CR_llll_invM"] = fs->make<TH1D>("4P0F_CR_llll_invM","4P0F CR M(4l);M [GeV];",500,0.,2500.);
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

void ResolvedEMuCRanalyzer::endJob() {
  purwgtFile_->Close();
  RMFFfile_->Close();
  REFFfile_->Close();
}

void ResolvedEMuCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  histo1d_["cutflow"]->Fill( 1.5, aWeight );

  std::vector<pat::ElectronRef> acceptEles;
  std::vector<pat::ElectronRef> nonHeepEles;

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

    if (!selected && aEle->et() > ptThres_ && aEle->ecalDriven()) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      nonHeepEles.push_back(castEle);
    }
  }

  if ( acceptEles.size() + nonHeepEles.size() != 2 )
    return;

  histo1d_["cutflow"]->Fill( 2.5, aWeight );

  if (muonHandle->empty())
    return;

  histo1d_["cutflow"]->Fill( 3.5, aWeight );

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

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 )
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

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 )
      isolatedHighPtTrackerMuons.push_back(firstMuon);
  }

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  if ( isolatedHighPtMuons.empty() )
    return;

  if ( isolatedHighPtMuons.front()->tunePMuonBestTrack()->pt() < ptThresTrig_ )
    return;

  histo1d_["cutflow"]->Fill( 4.5, aWeight );

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

  histo1d_["cutflow"]->Fill( 5.5, aWeight );

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

  if ( allHighPtMuons.size() + nonHighPtMuons.size() != 2 )
    return;

  histo1d_["cutflow"]->Fill( 6.5, aWeight );

  // SF
  if (isMC_) {
    for (const auto& aEle : acceptEles) {
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= recoEleSF_->GetBinContent( recoEleSF_->FindBin(aEle->superCluster()->eta(), apt ) );
      aWeight *= systHelperEle_.GetModifiedHeepSF(aEle);
    }

    for (const auto& aEle : nonHeepEles) {
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= recoEleSF_->GetBinContent( recoEleSF_->FindBin(aEle->superCluster()->eta(), apt ) );
    }

    for (const auto& aMu : isolatedHighPtMuons) {
      aWeight *= mucorrHelper_.idSF(aMu,idHistName_);
      aWeight *= mucorrHelper_.isoSF(aMu,isoHistName_);
    }

    for (const auto& aMu : isolatedHighPtTrackerMuons) {
      aWeight *= mucorrHelper_.idSF(aMu,idHistNameTrkHighPt_);
      aWeight *= mucorrHelper_.isoSF(aMu,isoHistNameTrkHighPt_);
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

  auto ciFF = [this] (const double dr, TF1* ffFit, TFitResultPtr& fitResult) {
    const double ff = ffFit->Eval(dr);
    const double xval[1] = {dr};
    double ci[1];
    fitResult->GetConfidenceIntervals(1,1,0,xval,ci,0.95,false);

    return std::hypot(ci[0],ff*ffSystCL95_);
  };

  switch ( allHighPtMuons.size() + acceptEles.size() ) {
    case 4:
      // 4P0F
      if ( nonHighPtMuons.empty() && nonHeepEles.empty() ) {
        auto lvecM1 = lvecFromTuneP( allHighPtMuons.at(0) );
        auto lvecM2 = lvecFromTuneP( allHighPtMuons.at(1) );
        const auto lvecE1 = acceptEles.at(0)->polarP4()*acceptEles.at(0)->userFloat("ecalTrkEnergyPostCorr")/acceptEles.at(0)->energy();
        const auto lvecE2 = acceptEles.at(1)->polarP4()*acceptEles.at(1)->userFloat("ecalTrkEnergyPostCorr")/acceptEles.at(1)->energy();

        if (isMC_) {
          edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
          iEvent.getByToken(genptcToken_, genptcHandle);

          lvecM1 *= mucorrHelper_.rochesterMC(allHighPtMuons.at(0),genptcHandle);
          lvecM2 *= mucorrHelper_.rochesterMC(allHighPtMuons.at(1),genptcHandle);
        } else {
          lvecM1 *= mucorrHelper_.rochesterData(allHighPtMuons.at(0));
          lvecM2 *= mucorrHelper_.rochesterData(allHighPtMuons.at(1));
        }

        const auto lvecA1 = lvecM1 + lvecM2;
        const auto lvecA2 = lvecE1 + lvecE2;
        const auto lvec4l = lvecA1 + lvecA2;
        const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
        const double dr2ll2 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
        const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
        const double m4l = lvec4l.M();

        if ( m4l > 50. && (m4l < 200. || isMC_) && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
          histo1d_["4P0F_CR_llll_invM"]->Fill(m4l, aWeight);
          histo1d_["4P0F_CR_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
          histo1d_["4P0F_CR_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
        }

        if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
          histo1d_["4P0F_CR_E_eta"]->Fill( acceptEles.front()->eta(), aWeight );
          histo1d_["4P0F_CR_E_phi"]->Fill( acceptEles.front()->phi(), aWeight );
          histo1d_["4P0F_CR_E_eta"]->Fill( acceptEles.at(1)->eta(), aWeight );
          histo1d_["4P0F_CR_E_phi"]->Fill( acceptEles.at(1)->phi(), aWeight );

          histo1d_["4P0F_CR_M_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["4P0F_CR_M_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
          histo1d_["4P0F_CR_M_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["4P0F_CR_M_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);

          histo1d_["4P0F_CR_pt_E"]->Fill( acceptEles.front()->et(), aWeight );
          histo1d_["4P0F_CR_pt_E"]->Fill( acceptEles.at(1)->et(), aWeight );
          histo1d_["4P0F_CR_pt_M"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["4P0F_CR_pt_M"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);

          histo1d_["4P0F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
          histo1d_["4P0F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
          histo1d_["4P0F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
          histo1d_["4P0F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
          histo1d_["4P0F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
          histo1d_["4P0F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
          histo1d_["4P0F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
          histo1d_["4P0F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);
        } // CR
      } // 4P0F

      break;
    case 3:
      // 3P1F
      if ( nonHighPtMuons.size()==1 || nonHeepEles.size()==1 ) {
        auto lvecM1 = lvecFromTuneP(allHighPtMuons.front());
        const auto secondMuon = nonHighPtMuons.empty() ? allHighPtMuons.at(1) : nonHighPtMuons.front();
        auto lvecM2 = lvecFromTuneP(secondMuon);
        const auto lvecE1 = acceptEles.at(0)->polarP4()*acceptEles.at(0)->userFloat("ecalTrkEnergyPostCorr")/acceptEles.at(0)->energy();
        const auto secondEle = nonHeepEles.empty() ? acceptEles.at(1) : nonHeepEles.front();
        const auto lvecE2 = secondEle->polarP4()*secondEle->userFloat("ecalTrkEnergyPostCorr")/secondEle->energy();
        const auto lvecFake = nonHighPtMuons.empty() ? lvecE2 : lvecM2;
        const auto lvecPair = nonHighPtMuons.empty() ? lvecE1 : lvecM1;

        histo1d_["3P1F_EP_pt"]->Fill(lvecE1.pt(), aWeight);
        histo1d_["3P1F_EP_eta"]->Fill(lvecE1.eta(), aWeight);
        histo1d_["3P1F_EP_phi"]->Fill(lvecE1.phi(), aWeight);

        histo1d_["3P1F_MP_pt"]->Fill(lvecM1.pt(), aWeight);
        histo1d_["3P1F_MP_eta"]->Fill(lvecM1.eta(), aWeight);
        histo1d_["3P1F_MP_phi"]->Fill(lvecM1.phi(), aWeight);

        if ( nonHighPtMuons.empty() ) {
          histo1d_["3P1F_MP_pt"]->Fill(lvecM2.pt(), aWeight);
          histo1d_["3P1F_MP_eta"]->Fill(lvecM2.eta(), aWeight);
          histo1d_["3P1F_MP_phi"]->Fill(lvecM2.phi(), aWeight);

          histo1d_["3P1F_EF_pt"]->Fill(lvecE2.pt(), aWeight);
          histo1d_["3P1F_EF_eta"]->Fill(lvecE2.eta(), aWeight);
          histo1d_["3P1F_EF_phi"]->Fill(lvecE2.phi(), aWeight);
        } else {
          histo1d_["3P1F_EP_pt"]->Fill(lvecE2.pt(), aWeight);
          histo1d_["3P1F_EP_eta"]->Fill(lvecE2.eta(), aWeight);
          histo1d_["3P1F_EP_phi"]->Fill(lvecE2.phi(), aWeight);

          histo1d_["3P1F_MF_pt"]->Fill(lvecM2.pt(), aWeight);
          histo1d_["3P1F_MF_eta"]->Fill(lvecM2.eta(), aWeight);
          histo1d_["3P1F_MF_phi"]->Fill(lvecM2.phi(), aWeight);
        }

        if (isMC_) {
          edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
          iEvent.getByToken(genptcToken_, genptcHandle);

          lvecM1 *= mucorrHelper_.rochesterMC( allHighPtMuons.front() ,genptcHandle);
          lvecM2 *= mucorrHelper_.rochesterMC( secondMuon ,genptcHandle);
        } else {
          lvecM1 *= mucorrHelper_.rochesterData(allHighPtMuons.front());
          lvecM2 *= mucorrHelper_.rochesterData(secondMuon);
        }

        const auto lvecA1 = lvecM1 + lvecM2;
        const auto lvecA2 = lvecE1 + lvecE2;
        const auto lvec4l = lvecA1 + lvecA2;
        const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
        const double dr2ll2 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
        const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
        const double m4l = lvec4l.M();

        const double drFake = reco::deltaR(lvecFake.eta(),lvecFake.phi(),lvecPair.eta(),lvecPair.phi());
        const double ff = nonHighPtMuons.empty() ? drREFF_->Eval(drFake) : drRMFF_->Eval(drFake);
        const double ciE = nonHighPtMuons.empty() ? ciFF(drFake, drREFF_, ffEleFit_) : 0.;
        const double ciM = nonHighPtMuons.empty() ? 0. : ciFF(drFake, drRMFF_, ffMuonFit_);

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
          histo1d_["3P1F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff);
          histo1d_["3P1F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff);
          histo1d_["3P1F_CR_llll_invM_xFF_ffUpE"]->Fill(m4l, aWeight*(ff+ciE));
          histo1d_["3P1F_CR_llll_invM_xFF_ffDnE"]->Fill(m4l, aWeight*(ff-std::min(ciE,ff)));
          histo1d_["3P1F_CR_llll_invM_xFF_ffUpM"]->Fill(m4l, aWeight*(ff+ciM));
          histo1d_["3P1F_CR_llll_invM_xFF_ffDnM"]->Fill(m4l, aWeight*(ff-std::min(ciM,ff)));
        }

        if ( m4l > 50. && m4l < 200. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
          histo1d_["3P1F_CR_llll_pt"]->Fill(lvec4l.pt(), aWeight);
          histo1d_["3P1F_CR_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);
          histo1d_["3P1F_CR_ll1_invM"]->Fill(lvecA1.M(), aWeight);
          histo1d_["3P1F_CR_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
          histo1d_["3P1F_CR_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
          histo1d_["3P1F_CR_ll2_invM"]->Fill(lvecA2.M(), aWeight);
          histo1d_["3P1F_CR_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
          histo1d_["3P1F_CR_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

          if ( nonHighPtMuons.empty() ) {
            histo1d_["3P1F_CR_pt_E"]->Fill( lvecE1.pt(), aWeight );
            histo1d_["3P1F_CR_dr_E"]->Fill( drFake, aWeight );
          } else {
            histo1d_["3P1F_CR_pt_M"]->Fill( lvecM1.pt(), aWeight );
            histo1d_["3P1F_CR_dr_M"]->Fill( drFake, aWeight );
          }

          histo1d_["3P1F_CR_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff);

          histo1d_["3P1F_CR_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff);
          histo1d_["3P1F_CR_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff);
          histo1d_["3P1F_CR_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff);
          histo1d_["3P1F_CR_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff);
        } // CR
      } // 3P1F

      break;
    case 2:
      // 2P2F
      if ( nonHighPtMuons.size() + nonHeepEles.size()==2 ) {
        auto lvecM1 = lvecFromTuneP(allHighPtMuons.front());
        const auto secondMuon = nonHighPtMuons.empty() ? allHighPtMuons.at(1) : nonHighPtMuons.front();
        auto lvecM2 = lvecFromTuneP(secondMuon);
        const auto firstEle = acceptEles.empty() ? nonHeepEles.at(0) : acceptEles.front();
        const auto lvecE1 = firstEle->polarP4()*firstEle->userFloat("ecalTrkEnergyPostCorr")/firstEle->energy();
        const auto secondEle = acceptEles.empty() ? nonHeepEles.at(1) : nonHeepEles.at(0);
        const auto lvecE2 = secondEle->polarP4()*secondEle->userFloat("ecalTrkEnergyPostCorr")/secondEle->energy();

        std::vector<edm::Ptr<reco::RecoCandidate>> allLeptons;
        const edm::Ptr<reco::RecoCandidate> firstElePtr = edm::refToPtr(firstEle);
        const edm::Ptr<reco::RecoCandidate> secondElePtr = edm::refToPtr(secondEle);
        allLeptons.push_back( edm::refToPtr(allHighPtMuons.front()) );
        allLeptons.push_back( edm::refToPtr(secondMuon) );
        allLeptons.push_back( firstElePtr );
        allLeptons.push_back( secondElePtr );
        std::pair<edm::Ptr<reco::RecoCandidate>,edm::Ptr<reco::RecoCandidate>> pair1, pair2;
        bool paired = PairingHelper::pairByDR(allLeptons,pair1,pair2);

        if (paired) {
          bool pairedEles = ( ( pair1.first==firstElePtr && pair1.second==secondElePtr ) ||
                              ( pair1.second==firstElePtr && pair1.first==secondElePtr ) );

          if ( !(pairedEles && acceptEles.empty()) ) {
            if (isMC_) {
              edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
              iEvent.getByToken(genptcToken_, genptcHandle);

              lvecM1 *= mucorrHelper_.rochesterMC(allHighPtMuons.front(),genptcHandle);
              lvecM2 *= mucorrHelper_.rochesterMC(secondMuon,genptcHandle);
            } else {
              lvecM1 *= mucorrHelper_.rochesterData(allHighPtMuons.front());
              lvecM2 *= mucorrHelper_.rochesterData(secondMuon);
            }

            const auto lvecA1 = lvecM1 + lvecM2;
            const auto lvecA2 = lvecE1 + lvecE2;
            const auto lvec4l = lvecA1 + lvecA2;
            const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
            const double dr2ll2 = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());
            const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
            const double m4l = lvec4l.M();

            const double drFakeE = reco::deltaR( lvecE1.eta(), lvecE1.phi(), lvecE2.eta(), lvecE2.phi() );
            const double drFake2 = acceptEles.empty() ? drFakeE : reco::deltaR( lvecM1.eta(), lvecM1.phi(), lvecM2.eta(), lvecM2.phi() );

            const double ff1 = drREFF_->Eval(drFakeE);
            const double ff2 = acceptEles.empty() ? drREFF_->Eval(drFake2) : drRMFF_->Eval(drFake2);
            const double ci1 = ciFF(drFakeE, drREFF_, ffEleFit_);
            const double ci2E = acceptEles.empty() ? ciFF(drFake2, drREFF_, ffEleFit_) : 0.;
            const double ci2M = acceptEles.empty() ? 0. : ciFF(drFake2, drRMFF_, ffMuonFit_);

            if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
              histo1d_["2P2F_MP_pt"]->Fill( lvecM1.pt(), aWeight );
              histo1d_["2P2F_MP_eta"]->Fill( lvecM1.eta(), aWeight );
              histo1d_["2P2F_MP_phi"]->Fill( lvecM1.phi(), aWeight );

              histo1d_["2P2F_EF_pt"]->Fill( lvecE2.pt(), aWeight );
              histo1d_["2P2F_EF_eta"]->Fill( lvecE2.eta(), aWeight );
              histo1d_["2P2F_EF_phi"]->Fill( lvecE2.phi(), aWeight );

              if (acceptEles.empty()) {
                histo1d_["2P2F_EF_pt"]->Fill( lvecE1.pt(), aWeight );
                histo1d_["2P2F_EF_eta"]->Fill( lvecE1.eta(), aWeight );
                histo1d_["2P2F_EF_phi"]->Fill( lvecE1.phi(), aWeight );

                histo1d_["2P2F_MP_pt"]->Fill( lvecM2.pt(), aWeight );
                histo1d_["2P2F_MP_eta"]->Fill( lvecM2.eta(), aWeight );
                histo1d_["2P2F_MP_phi"]->Fill( lvecM2.phi(), aWeight );
              } else {
                histo1d_["2P2F_EP_pt"]->Fill( lvecE1.pt(), aWeight );
                histo1d_["2P2F_EP_eta"]->Fill( lvecE1.eta(), aWeight );
                histo1d_["2P2F_EP_phi"]->Fill( lvecE1.phi(), aWeight );

                histo1d_["2P2F_MF_pt"]->Fill( lvecM2.pt(), aWeight );
                histo1d_["2P2F_MF_eta"]->Fill( lvecM2.eta(), aWeight );
                histo1d_["2P2F_MF_phi"]->Fill( lvecM2.phi(), aWeight );
              }

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
              histo1d_["2P2F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF_ffUpE"]->Fill(m4l, aWeight*(ff1+ci1));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDnE"]->Fill(m4l, aWeight*(ff1-std::min(ff1,ci1)));
              histo1d_["2P2F_CR_llll_invM_xFF_ffUpM"]->Fill(m4l, aWeight*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF_ffDnM"]->Fill(m4l, aWeight*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF_ffUpE"]->Fill(m4l, aWeight*(ff2+ci2E));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDnE"]->Fill(m4l, aWeight*(ff2-std::min(ff2,ci2E)));
              histo1d_["2P2F_CR_llll_invM_xFF_ffUpM"]->Fill(m4l, aWeight*(ff2+ci2M));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDnM"]->Fill(m4l, aWeight*(ff2-std::min(ff2,ci2M)));
              histo1d_["2P2F_CR_llll_invM_xFF2"]->Fill(m4l, aWeight*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2_ffUpE"]->Fill(m4l, aWeight*(ff1+ci1)*(ff2+ci2E));
              histo1d_["2P2F_CR_llll_invM_xFF2_ffDnE"]->Fill(m4l, aWeight*(ff1-std::min(ff1,ci1))*(ff2-std::min(ff2,ci2E)));
              histo1d_["2P2F_CR_llll_invM_xFF2_ffUpM"]->Fill(m4l, aWeight*ff1*(ff2+ci2M));
              histo1d_["2P2F_CR_llll_invM_xFF2_ffDnM"]->Fill(m4l, aWeight*ff1*(ff2-std::min(ff2,ci2M)));
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

              histo1d_["2P2F_CR_pt_E"]->Fill( lvecE2.pt(), aWeight );
              histo1d_["2P2F_CR_dr_E"]->Fill( drFakeE, aWeight );
              histo1d_["2P2F_CR_pt_E_xFF"]->Fill( lvecE2.pt(), aWeight*ff2 );
              histo1d_["2P2F_CR_dr_E_xFF"]->Fill( drFakeE, aWeight*ff2 );

              if (acceptEles.empty()) {
                histo1d_["2P2F_CR_pt_E"]->Fill( lvecE1.pt(), aWeight );
                histo1d_["2P2F_CR_dr_E"]->Fill( drFakeE, aWeight );
                histo1d_["2P2F_CR_pt_E_xFF"]->Fill( lvecE1.pt(), aWeight*ff1 );
                histo1d_["2P2F_CR_dr_E_xFF"]->Fill( drFakeE, aWeight*ff1 );
              } else {
                histo1d_["2P2F_CR_pt_M"]->Fill( lvecM2.pt(), aWeight );
                histo1d_["2P2F_CR_dr_M"]->Fill( drFake2, aWeight );
                histo1d_["2P2F_CR_pt_M_xFF"]->Fill( lvecM2.pt(), aWeight*ff2 );
                histo1d_["2P2F_CR_dr_M_xFF"]->Fill( drFake2, aWeight*ff2 );
              }
            } // CR
          }
        } // paired
      } // 2P2F

      break;
  } // switch allHighPtMuons.size()
}

DEFINE_FWK_MODULE(ResolvedEMuCRanalyzer);
