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

#include "ZprimeTo4l/Analysis/interface/PairingHelper.h"
#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TFitResult.h"

#include "correction.h"

class ResolvedMuCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit ResolvedMuCRanalyzer(const edm::ParameterSet&);
  virtual ~ResolvedMuCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  math::PtEtaPhiMLorentzVector lvecFromTuneP(const pat::MuonRef& aMu);
  std::pair<double,double> pairByInvM(const math::PtEtaPhiMLorentzVector& a,
                                      const math::PtEtaPhiMLorentzVector& b,
                                      const math::PtEtaPhiMLorentzVector& c,
                                      const math::PtEtaPhiMLorentzVector& d) const;

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<double> prefweightUp_token_;
  const edm::EDGetTokenT<double> prefweightDn_token_;

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath muonIdIsoSFpath_;
  const edm::FileInPath muonBoostIsoSFpath_;
  const edm::FileInPath muonRecoSFpath_;

  const edm::FileInPath FFpath_;
  const edm::FileInPath FFdr03path_;

  std::unique_ptr<correction::CorrectionSet> purwgt_;
  const std::string puname_;

  const std::vector<double> muScaleBias_;
  const std::vector<double> muSmearFactors_;
  const std::vector<double> muSmearParams_;

  const double ptThres_;
  const double ptThresTrig_;
  const double ffSystCL95_;

  MuonCorrectionHelper mucorrHelper_;

  std::unique_ptr<TFile> FFfile_;
  // TF1* drFF_;
  TF1* ptFF_;
  TFitResultPtr ffFit_;

  std::unique_ptr<TFile> FFdr03file_;
  TF1* ffdr03_;
  TFitResultPtr ffdr03FitResult_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  const double mumass_ = 0.1056583745;

  TTree* numerTree_ = nullptr;
  float numerPt_ = -1.;
  float numerDr_ = -1.;
  float numerWgt_ = 0.;
  float numerInvM_ = -1.;
  int numerIsMB_ = -1;

  TTree* denomTree_ = nullptr;
  float denomPt_ = -1.;
  float denomDr_ = -1.;
  float denomIso_ = -1.;
  float denomWgt_ = 0.;
  float denomInvM_ = -1.;
  int denomIsMB_ = -1;
  int denomPassHighPt_ = -1;

  TTree* tree_3P1F_ = nullptr;
  float pt_3P1F_ = -1.;
  float dr_3P1F_ = -1.;
  float iso_3P1F_ = -1.;
  float mll_3P1F_ = -1.;
  float m4l_3P1F_ = -1.;
  float wgt_3P1F_ = 0.;
  int passHighPt_3P1F_ = -1;

  TTree* tree_2P2F_ = nullptr;
  float pt1_2P2F_ = -1.;
  float pt2_2P2F_ = -1.;
  float dr1_2P2F_ = -1.;
  float dr2_2P2F_ = -1.;
  float iso1_2P2F_ = -1.;
  float iso2_2P2F_ = -1.;
  float mll1_2P2F_ = -1.;
  float mll2_2P2F_ = -1.;
  float m4l_2P2F_ = -1.;
  float wgt_2P2F_ = 0.;
  int pairedFakes_2P2F_ = -1;
  int passHighPt_2P2F_ = -1;

  TTree* tree_2P2F_dr03_ = nullptr;
  float var_weight_ = 0.;
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

ResolvedMuCRanalyzer::ResolvedMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
prefweightUp_token_(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
prefweightDn_token_(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
muonIdIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonIdIsoSFpath")),
muonBoostIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonBoostIsoSFpath")),
muonRecoSFpath_(iConfig.getParameter<edm::FileInPath>("muonRecoSFpath")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
FFdr03path_(iConfig.getParameter<edm::FileInPath>("FFdr03path")),
purwgt_(std::move(correction::CorrectionSet::from_file((iConfig.getParameter<edm::FileInPath>("PUrwgt")).fullPath()))),
puname_(iConfig.getParameter<std::string>("PUname")),
muScaleBias_(iConfig.getParameter<std::vector<double>>("muScaleBias")),
muSmearFactors_(iConfig.getParameter<std::vector<double>>("muSmearFactors")),
muSmearParams_(iConfig.getParameter<std::vector<double>>("muSmearParams")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptThresTrig_(iConfig.getParameter<double>("ptThresTrig")),
ffSystCL95_(iConfig.getParameter<double>("ffSystCL95")),
mucorrHelper_(rochesterPath_,triggerSFpath_,muonIdIsoSFpath_,muonBoostIsoSFpath_,muonRecoSFpath_) {
  usesResource("TFileService");
}

math::PtEtaPhiMLorentzVector ResolvedMuCRanalyzer::lvecFromTuneP(const pat::MuonRef& aMu) {
  return math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                      aMu->tunePMuonBestTrack()->eta(),
                                      aMu->tunePMuonBestTrack()->phi(),
                                      mumass_);
}

std::pair<double,double> ResolvedMuCRanalyzer::pairByInvM(const math::PtEtaPhiMLorentzVector& a,
                                                          const math::PtEtaPhiMLorentzVector& b,
                                                          const math::PtEtaPhiMLorentzVector& c,
                                                          const math::PtEtaPhiMLorentzVector& d) const {
  const double val11 = (a+b).M();
  const double val12 = (c+d).M();
  auto cand1 = std::make_pair( std::max(val11,val12), std::min(val11,val12) );

  const double val21 = (a+c).M();
  const double val22 = (b+d).M();
  auto cand2 = std::make_pair( std::max(val21,val22), std::min(val21,val22) );

  const double val31 = (a+d).M();
  const double val32 = (b+c).M();
  auto cand3 = std::make_pair( std::max(val31,val32), std::min(val31,val32) );

  std::vector<std::pair<double,double>> vec = {cand1,cand2,cand3};

  auto sortByInvM = [] (const std::pair<double,double>& i, const std::pair<double,double>& j) {
    return i.second < j.second;
  };

  std::sort(vec.begin(),vec.end(),sortByInvM);

  return vec.front();
}

void ResolvedMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  // drFF_ = static_cast<TF1*>(FFfile_->FindObjectAny("RMFF_dr_all"));
  ptFF_ = static_cast<TF1*>(FFfile_->FindObjectAny("RMFFfunc"));
  ffFit_ = (static_cast<TH1D*>(FFfile_->Get("RMFF_pt")))->Fit(ptFF_,"RS");

  FFdr03file_ = std::make_unique<TFile>(FFdr03path_.fullPath().c_str(),"READ");
  ffdr03_ = static_cast<TF1*>(FFdr03file_->FindObjectAny("RMFFdr03"));
  ffdr03FitResult_ = (static_cast<TH1D*>(FFdr03file_->Get("RMFFhist")))->Fit(ffdr03_,"RS");

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
  histo1d_["2P2F_CR_llll_invM_altMuScale_xFF"] = fs->make<TH1D>("2P2F_CR_llll_invM_altMuScale_xFF",";M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF_ffDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF_ffDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF","2P2F CR dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF","2P2F CR M(ll1) x Fake factor;M [GeV];",400,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF","2P2F CR dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF","2P2F CR M(ll2) x Fake factor;M [GeV];",400,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF","2P2F CR dR(ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CR_llHiM_invM_xFF"] = fs->make<TH1D>("2P2F_CR_llHiM_invM_xFF",";M [GeV];",400,0.,200.);
  histo1d_["2P2F_CR_llLoM_invM_xFF"] = fs->make<TH1D>("2P2F_CR_llLoM_invM_xFF",";M [GeV];",400,0.,200.);

  histo1d_["2P2F_CR_llll_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_altMuScale_xFF2"] = fs->make<TH1D>("2P2F_CR_llll_invM_altMuScale_xFF2",";M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffUp"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffUp","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_llll_invM_xFF2_ffDn"] = fs->make<TH1D>("2P2F_CR_llll_invM_xFF2_ffDn","2P2F CR M(4l) x Fake factor^2;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CR_ll1ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1ll2_dr_xFF2","2P2F CR dR(ll1ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll1_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_invM_xFF2","2P2F CR M(ll1) x Fake factor^2;M [GeV];",400,0.,200.);
  histo1d_["2P2F_CR_ll1_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll1_dr_xFF2","2P2F CR dR(ll1) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_ll2_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_invM_xFF2","2P2F CR M(ll2) x Fake factor^2;M [GeV];",400,0.,200.);
  histo1d_["2P2F_CR_ll2_dr_xFF2"] = fs->make<TH1D>("2P2F_CR_ll2_dr_xFF2","2P2F CR dR(ll2) x Fake factor^2",128,0.,6.4);
  histo1d_["2P2F_CR_llHiM_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_llHiM_invM_xFF2",";M [GeV];",400,0.,200.);
  histo1d_["2P2F_CR_llLoM_invM_xFF2"] = fs->make<TH1D>("2P2F_CR_llLoM_invM_xFF2",";M [GeV];",400,0.,200.);

  // 2P2F noniso
  histo1d_["2P2F_CRdr03_llll_invM"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_pt"] = fs->make<TH1D>("2P2F_CRdr03_llll_pt","2P2F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1ll2_dr"] = fs->make<TH1D>("2P2F_CRdr03_ll1ll2_dr","2P2F CR dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll1_invM"] = fs->make<TH1D>("2P2F_CRdr03_ll1_invM","2P2F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1_pt"] = fs->make<TH1D>("2P2F_CRdr03_ll1_pt","2P2F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1_dr"] = fs->make<TH1D>("2P2F_CRdr03_ll1_dr","2P2F CR dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll2_invM"] = fs->make<TH1D>("2P2F_CRdr03_ll2_invM","2P2F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_pt"] = fs->make<TH1D>("2P2F_CRdr03_ll2_pt","2P2F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_dr"] = fs->make<TH1D>("2P2F_CRdr03_ll2_dr","2P2F CR dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CRdr03_llll_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_altMuScale_xFF"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_altMuScale_xFF",";M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_xFF_ffUp"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF_ffUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_xFF_ffDn"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF_ffDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1ll2_dr_xFF","2P2F CR dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1_invM_xFF","2P2F CR M(ll1) x Fake factor;M [GeV];",400,0.,200.);
  histo1d_["2P2F_CRdr03_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1_dr_xFF","2P2F CR dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll2_invM_xFF","2P2F CR M(ll2) x Fake factor;M [GeV];",400,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll2_dr_xFF","2P2F CR dR(ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CRdr03_llHiM_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_llHiM_invM_xFF",";M [GeV];",400,0.,200.);
  histo1d_["2P2F_CRdr03_llLoM_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_llLoM_invM_xFF",";M [GeV];",400,0.,200.);

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
  histo1d_["3P1F_CR_llll_invM_altMuScale_xFF"] = fs->make<TH1D>("3P1F_CR_llll_invM_altMuScale_xFF",";M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffUp"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffUp","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_llll_invM_xFF_ffDn"] = fs->make<TH1D>("3P1F_CR_llll_invM_xFF_ffDn","3P1F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["3P1F_CR_ll1ll2_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll1ll2_dr_xFF","3P1F CR dR(ll1ll2) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_ll1_invM_xFF"] = fs->make<TH1D>("3P1F_CR_ll1_invM_xFF","3P1F CR M(ll1) x Fake factor;M [GeV];",400,0.,200.);
  histo1d_["3P1F_CR_ll1_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll1_dr_xFF","3P1F CR dR(ll1) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_ll2_invM_xFF"] = fs->make<TH1D>("3P1F_CR_ll2_invM_xFF","3P1F CR M(ll2) x Fake factor;M [GeV];",400,0.,200.);
  histo1d_["3P1F_CR_ll2_dr_xFF"] = fs->make<TH1D>("3P1F_CR_ll2_dr_xFF","3P1F CR dR(ll2) x Fake factor",64,0.,6.4);
  histo1d_["3P1F_CR_llHiM_invM_xFF"] = fs->make<TH1D>("3P1F_CR_llHiM_invM_xFF",";M [GeV];",400,0.,200.);
  histo1d_["3P1F_CR_llLoM_invM_xFF"] = fs->make<TH1D>("3P1F_CR_llLoM_invM_xFF",";M [GeV];",400,0.,200.);

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
  histo1d_["4P0F_CR_llll_invM_altMuScale"] = fs->make<TH1D>("4P0F_CR_llll_invM_altMuScale","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_altMuSmear"] = fs->make<TH1D>("4P0F_CR_llll_invM_altMuSmear","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_idUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_idUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_idDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_idDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_isoUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_isoUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_isoDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_isoDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_muBoostIsoUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_muBoostIsoUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_muBoostIsoDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_muBoostIsoDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_trigUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_trigUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_trigDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_trigDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_recoUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_recoUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_recoDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_recoDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_PUrwgtUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_PUrwgtUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_PUrwgtDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_PUrwgtDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_prefireUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_prefireUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_prefireDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_prefireDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_pt"] = fs->make<TH1D>("4P0F_CR_llll_pt","4P0F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll1ll2_dr","4P0F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["4P0F_CR_ll1_invM"] = fs->make<TH1D>("4P0F_CR_ll1_invM","4P0F CR M(ll1);M [GeV];",400,0.,200.);
  histo1d_["4P0F_CR_ll1_pt"] = fs->make<TH1D>("4P0F_CR_ll1_pt","4P0F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_dr"] = fs->make<TH1D>("4P0F_CR_ll1_dr","4P0F CR dR(ll1)",64,0.,6.4);
  histo1d_["4P0F_CR_ll2_invM"] = fs->make<TH1D>("4P0F_CR_ll2_invM","4P0F CR M(ll2);M [GeV];",400,0.,200.);
  histo1d_["4P0F_CR_ll2_pt"] = fs->make<TH1D>("4P0F_CR_ll2_pt","4P0F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll2_dr","4P0F CR dR(ll2)",64,0.,6.4);
  histo1d_["4P0F_CR_llHiM_invM"] = fs->make<TH1D>("4P0F_CR_llHiM_invM",";M [GeV];",400,0.,200.);
  histo1d_["4P0F_CR_llLoM_invM"] = fs->make<TH1D>("4P0F_CR_llLoM_invM",";M [GeV];",400,0.,200.);

  numerTree_ = fs->make<TTree>("numerTree","numerTree");
  numerTree_->Branch("pt",&numerPt_,"pt/F");
  numerTree_->Branch("dr",&numerDr_,"dr/F");
  numerTree_->Branch("wgt",&numerWgt_,"wgt/F");
  numerTree_->Branch("invM",&numerInvM_,"invM/F");
  numerTree_->Branch("isMB",&numerIsMB_,"isMB/I");

  denomTree_ = fs->make<TTree>("denomTree","denomTree");
  denomTree_->Branch("pt",&denomPt_,"pt/F");
  denomTree_->Branch("dr",&denomDr_,"dr/F");
  denomTree_->Branch("iso",&denomIso_,"iso/F");
  denomTree_->Branch("wgt",&denomWgt_,"wgt/F");
  denomTree_->Branch("invM",&denomInvM_,"invM/F");
  denomTree_->Branch("isMB",&denomIsMB_,"isMB/I");
  denomTree_->Branch("passHighPt",&denomPassHighPt_,"passHighPt/I");

  tree_3P1F_ = fs->make<TTree>("tree_3P1F","tree_3P1F");
  tree_3P1F_->Branch("pt",&pt_3P1F_,"pt/F");
  tree_3P1F_->Branch("dr",&dr_3P1F_,"dr/F");
  tree_3P1F_->Branch("iso",&iso_3P1F_,"iso/F");
  tree_3P1F_->Branch("mll",&mll_3P1F_,"mll/F");
  tree_3P1F_->Branch("m4l",&m4l_3P1F_,"m4l/F");
  tree_3P1F_->Branch("wgt",&wgt_3P1F_,"wgt/F");
  tree_3P1F_->Branch("passHighPt",&passHighPt_3P1F_,"passHighPt/I");

  tree_2P2F_ = fs->make<TTree>("tree_2P2F","tree_2P2F");
  tree_2P2F_->Branch("pt1",&pt1_2P2F_,"pt1/F");
  tree_2P2F_->Branch("pt2",&pt2_2P2F_,"pt2/F");
  tree_2P2F_->Branch("dr1",&dr1_2P2F_,"dr1/F");
  tree_2P2F_->Branch("dr2",&dr2_2P2F_,"dr2/F");
  tree_2P2F_->Branch("iso1",&iso1_2P2F_,"iso1/F");
  tree_2P2F_->Branch("iso2",&iso2_2P2F_,"iso2/F");
  tree_2P2F_->Branch("mll1",&mll1_2P2F_,"mll1/F");
  tree_2P2F_->Branch("mll2",&mll2_2P2F_,"mll2/F");
  tree_2P2F_->Branch("m4l",&m4l_2P2F_,"m4l/F");
  tree_2P2F_->Branch("wgt",&wgt_2P2F_,"wgt/F");
  tree_2P2F_->Branch("pairedFakes",&pairedFakes_2P2F_,"pairedFakes/I");
  tree_2P2F_->Branch("passHighPt",&passHighPt_2P2F_,"passHighPt/I");

  tree_2P2F_dr03_ = fs->make<TTree>("tree_2P2F_dr03","tree_2P2F_dr03");
  tree_2P2F_dr03_->Branch("wgt",&var_weight_,"wgt/F");
  tree_2P2F_dr03_->Branch("m1Pt",&var_m1Pt_,"m1Pt/F");
  tree_2P2F_dr03_->Branch("m1Eta",&var_m1Eta_,"m1Eta/F");
  tree_2P2F_dr03_->Branch("m1Phi",&var_m1Phi_,"m1Phi/F");
  tree_2P2F_dr03_->Branch("m1TrkIso",&var_m1TrkIso_,"m1TrkIso/F");
  tree_2P2F_dr03_->Branch("m2Pt",&var_m2Pt_,"m2Pt/F");
  tree_2P2F_dr03_->Branch("m2Eta",&var_m2Eta_,"m2Eta/F");
  tree_2P2F_dr03_->Branch("m2Phi",&var_m2Phi_,"m2Phi/F");
  tree_2P2F_dr03_->Branch("m2TrkIso",&var_m2TrkIso_,"m2TrkIso/F");
  tree_2P2F_dr03_->Branch("m1m2InvM",&var_m1m2InvM_,"m1m2InvM/F");
  tree_2P2F_dr03_->Branch("m1m2Pt",&var_m1m2Pt_,"m1m2Pt/F");
  tree_2P2F_dr03_->Branch("m1m2DR",&var_m1m2DR_,"m1m2DR/F");
}

void ResolvedMuCRanalyzer::endJob() {
  FFfile_->Close();
}

void ResolvedMuCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<pat::MET>> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  double aWeight = 1.;
  double purwgtUp = 1.;
  double purwgtDn = 1.;
  double purwgtNo = 1.;
  double prefireNo = 1.;
  double prefireUp = 1.;
  double prefireDn = 1.;

  if (isMC_) {
    edm::Handle<double> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    prefireNo = *theprefweight;

    edm::Handle<double> theprefweightUp;
    iEvent.getByToken(prefweightUp_token_, theprefweightUp);
    prefireUp = *theprefweightUp;

    edm::Handle<double> theprefweightDn;
    iEvent.getByToken(prefweightDn_token_, theprefweightDn);
    prefireDn = *theprefweightDn;

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
        purwgtUp = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"up"});
        purwgtDn = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"down"});
        purwgtNo = purwgt_->at(puname_)->evaluate({apu->getTrueNumInteractions(),"nominal"});

        aWeight *= purwgtNo;

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

  if (muonHandle->empty())
    return;

  histo1d_["cutflow"]->Fill( 2.5, aWeight );

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons; // but not highPtMuon
  std::vector<pat::MuonRef> nonHighPtMuons; // pass acceptance but not ID
  std::vector<pat::MuonRef> nonHighPtMuonsVLiso; // pass acceptance but not ID
  std::map<pat::MuonRef,float> nonHighPtIsos;

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

  if ( allHighPtMuons.empty() || allHighPtMuons.front()->tunePMuonBestTrack()->pt() < ptThresTrig_ )
    return;

  histo1d_["cutflow"]->Fill( 3.5, aWeight );

  bool trigMatched = false;
  std::vector<double> trigSyst;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    const auto& leadMu = allHighPtMuons.front();
    const bool isGlobal = std::find(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),leadMu)!=isolatedHighPtMuons.end();

    if ( reco::deltaR2(trigObj->eta(),
                       trigObj->phi(),
                       leadMu->tunePMuonBestTrack()->eta(),
                       leadMu->tunePMuonBestTrack()->phi()) < 0.01 ) {
      trigMatched = true;

      if (isMC_) {
        if (isGlobal) {
          aWeight *= mucorrHelper_.trigSFglobal(leadMu);
          trigSyst.push_back( mucorrHelper_.trigSFglobalSyst(leadMu) /
                              mucorrHelper_.trigSFglobal(leadMu) );
        } else {
          aWeight *= mucorrHelper_.trigSFtracker(leadMu);
          trigSyst.push_back( mucorrHelper_.trigSFtrackerSyst(leadMu) /
                              mucorrHelper_.trigSFtracker(leadMu) );
        }
      }

      break;
    }
  }

  if ( !trigMatched )
    return;

  histo1d_["cutflow"]->Fill( 4.5, aWeight );

  mucorrHelper_.nonHighPtMuonIso(nonHighPtMuons,
                                 nonHighPtMuonsVLiso,
                                 nonHighPtIsos,
                                 muonHandle,
                                 allHighPtMuons,
                                 highPtMuons,
                                 highPtTrackerMuons,
                                 *primaryVertex,
                                 *beamSpotHandle,
                                 ptThres_);

  bool hasEle = false;

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( std::abs(aEle->superCluster()->eta()) > 2.5 )
      continue;

    // veto EBEE gap
    if ( std::abs(aEle->superCluster()->eta()) > 1.4442 && std::abs(aEle->superCluster()->eta()) < 1.566 )
      continue;

    if ( aEle->electronID("modifiedHeepElectronID") )
      hasEle = true;
  }

  if (hasEle)
    return;

  histo1d_["cutflow"]->Fill( 5.5, aWeight );

  std::vector<double> idSyst;
  std::vector<double> isoSyst;
  std::vector<double> recoSyst;

  auto systSFratio = [] (const std::vector<double>& vec) -> std::pair<double,double> {
    double up = 1., dn = 1.;

    for (const auto& systOverSF : vec) {
      up *= (1.+systOverSF);
      dn *= std::max(1.-systOverSF,0.);
    }

    return std::make_pair(up,dn); // ratio to the nominal
  };

  auto boostIsoSFupdn = [this,&boostedMuons,&primaryVertex] (const double wgt) -> std::pair<double,double> {
    if (!isMC_)
      return std::make_pair(wgt,wgt);

    double wgtUp = wgt;
    double wgtDn = wgt;

    for (const auto& aMu : boostedMuons) {
      auto apair = mucorrHelper_.boostIsoSFupdn(aMu,*primaryVertex);
      wgtUp *= apair.first;
      wgtDn *= apair.second;
    }

    return std::make_pair(wgtUp,wgtDn);
  };

  if (isMC_) {
    for (const auto& aMu : isolatedHighPtMuons) {
      aWeight *= mucorrHelper_.highptIdSF(aMu);
      aWeight *= mucorrHelper_.looseIsoSF(aMu);
      aWeight *= mucorrHelper_.recoSF(aMu);
      idSyst.push_back( mucorrHelper_.highptIdSFsyst(aMu)/mucorrHelper_.highptIdSF(aMu) );
      recoSyst.push_back( mucorrHelper_.recoSFsyst(aMu)/mucorrHelper_.recoSF(aMu) );

      if (std::find(boostedMuons.begin(),boostedMuons.end(),aMu)==boostedMuons.end())
        isoSyst.push_back( mucorrHelper_.looseIsoSFsyst(aMu)/mucorrHelper_.looseIsoSF(aMu) );
    }

    for (const auto& aMu : isolatedHighPtTrackerMuons) {
      aWeight *= mucorrHelper_.trkHighptIdSF(aMu);
      aWeight *= mucorrHelper_.looseIsoSFtracker(aMu);
      aWeight *= mucorrHelper_.recoSF(aMu);
      idSyst.push_back( mucorrHelper_.trkHighptIdSFsyst(aMu)/mucorrHelper_.trkHighptIdSF(aMu) );

      if (std::find(boostedMuons.begin(),boostedMuons.end(),aMu)==boostedMuons.end())
        isoSyst.push_back( mucorrHelper_.looseIsoSFtrackerSyst(aMu)/mucorrHelper_.looseIsoSFtracker(aMu) );
    }

    for (const auto& aMu : nonHighPtMuonsVLiso) {
      aWeight *= mucorrHelper_.recoSF(aMu);
      recoSyst.push_back( mucorrHelper_.recoSFsyst(aMu)/mucorrHelper_.recoSF(aMu) );
    }
  }

  auto checkTrackerMuPair = [&isolatedHighPtMuons] (const std::pair<pat::MuonRef,pat::MuonRef>& apair, bool checkId=true) -> bool {
    if (!checkId)
      return !apair.first->isGlobalMuon() && !apair.second->isGlobalMuon();

    bool check1st = std::find(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),apair.first) != isolatedHighPtMuons.end();
    bool check2nd = std::find(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),apair.second) != isolatedHighPtMuons.end();

    return !check1st && !check2nd;
  };

  auto valFF = [this] (const double dr, const double pt) {
    // const double ff1 = drFF_->Eval(dr);
    // const double ff2 = ptFF_->Eval(pt);

    // return (ff1+ff2)/2.;
    return ptFF_->Eval(pt);
  };

  auto ciFF = [this] (const double dr, const double pt) {
    // const double ff1 = drFF_->Eval(dr);
    // const double ff2 = ptFF_->Eval(pt);
    //
    // return std::abs(ff1-ff2)/2.;

    const double xval[1] = {pt};
    double ci[1];
    ffFit_->GetConfidenceIntervals(1,1,0,xval,ci,0.68,false);

    return std::abs(ci[0]);
  };

  auto ciFFdr03 = [this] (const double ptll) {
    const double xval[1] = {ptll};
    double ci[1];
    ffdr03FitResult_->GetConfidenceIntervals(1,1,0,xval,ci,0.68,false);

    return ci[0];
  };

  switch (allHighPtMuons.size()) {
    case 4:
      if (true) {
        std::vector<pat::MuonRef> allMuons(allHighPtMuons);
        std::pair<pat::MuonRef,pat::MuonRef> pair1, pair2;

        bool paired = PairingHelper::pair4M(allHighPtMuons,pair1,pair2);

        if (paired && !checkTrackerMuPair(pair1) && !checkTrackerMuPair(pair2)) {
          auto lvecM1 = lvecFromTuneP(pair1.first);
          auto lvecM2 = lvecFromTuneP(pair1.second);
          auto lvecM3 = lvecFromTuneP(pair2.first);
          auto lvecM4 = lvecFromTuneP(pair2.second);

          lvecM1 *= mucorrHelper_.smear(pair1.first,muSmearParams_,muSmearFactors_,isMC_);
          lvecM2 *= mucorrHelper_.smear(pair1.second,muSmearParams_,muSmearFactors_,isMC_);
          lvecM3 *= mucorrHelper_.smear(pair2.first,muSmearParams_,muSmearFactors_,isMC_);
          lvecM4 *= mucorrHelper_.smear(pair2.second,muSmearParams_,muSmearFactors_,isMC_);
          const auto lvecM1alt = lvecM1*mucorrHelper_.altScale(pair1.first,muScaleBias_,isMC_);
          const auto lvecM2alt = lvecM2*mucorrHelper_.altScale(pair1.second,muScaleBias_,isMC_);
          const auto lvecM3alt = lvecM3*mucorrHelper_.altScale(pair2.first,muScaleBias_,isMC_);
          const auto lvecM4alt = lvecM4*mucorrHelper_.altScale(pair2.second,muScaleBias_,isMC_);

          // if (isMC_) {
          //   edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
          //   iEvent.getByToken(genptcToken_, genptcHandle);
          //
          //   lvecM1 *= mucorrHelper_.nominalMC(pair1.first,genptcHandle);
          //   lvecM2 *= mucorrHelper_.nominalMC(pair1.second,genptcHandle);
          //   lvecM3 *= mucorrHelper_.nominalMC(pair2.first,genptcHandle);
          //   lvecM4 *= mucorrHelper_.nominalMC(pair2.second,genptcHandle);
          // } else {
          //   lvecM1 *= mucorrHelper_.nominalData(pair1.first);
          //   lvecM2 *= mucorrHelper_.nominalData(pair1.second);
          //   lvecM3 *= mucorrHelper_.nominalData(pair2.first);
          //   lvecM4 *= mucorrHelper_.nominalData(pair2.second);
          // }

          const auto lvecM1smear = lvecFromTuneP(pair1.first)*mucorrHelper_.smear(pair1.first,muSmearParams_,{0.46,0.46},isMC_);
          const auto lvecM2smear = lvecFromTuneP(pair1.second)*mucorrHelper_.smear(pair1.second,muSmearParams_,{0.46,0.46},isMC_);
          const auto lvecM3smear = lvecFromTuneP(pair2.first)*mucorrHelper_.smear(pair2.first,muSmearParams_,{0.46,0.46},isMC_);
          const auto lvecM4smear = lvecFromTuneP(pair2.second)*mucorrHelper_.smear(pair2.second,muSmearParams_,{0.46,0.46},isMC_);

          const auto lvecA1 = lvecM1 + lvecM2;
          const auto lvecA2 = lvecM3 + lvecM4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
          const double dr2ll2 = reco::deltaR2(lvecM3.eta(),lvecM3.phi(),lvecM4.eta(),lvecM4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = std::min(lvec4l.M(),2499.9); // truncate at 2.5 TeV
          const double m4lAlt = std::min((lvecM1alt+lvecM2alt+lvecM3alt+lvecM4alt).M(),2499.9);
          const double m4lsmear = std::min((lvecM1smear+lvecM2smear+lvecM3smear+lvecM4smear).M(),2499.9);

          if ( m4l > 50. /*&& (m4l < 500. || isMC_) (unblinded)*/ && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["4P0F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["4P0F_CR_llll_invM_altMuScale"]->Fill(m4lAlt, aWeight);
            histo1d_["4P0F_CR_llll_invM_altMuSmear"]->Fill(m4lsmear, aWeight);
            histo1d_["4P0F_CR_llll_invM_idUp"]->Fill(m4l, aWeight*systSFratio(idSyst).first);
            histo1d_["4P0F_CR_llll_invM_idDn"]->Fill(m4l, aWeight*systSFratio(idSyst).second);
            histo1d_["4P0F_CR_llll_invM_isoUp"]->Fill(m4l, aWeight*systSFratio(isoSyst).first);
            histo1d_["4P0F_CR_llll_invM_isoDn"]->Fill(m4l, aWeight*systSFratio(isoSyst).second);
            histo1d_["4P0F_CR_llll_invM_muBoostIsoUp"]->Fill(m4l, boostIsoSFupdn(aWeight).first);
            histo1d_["4P0F_CR_llll_invM_muBoostIsoDn"]->Fill(m4l, boostIsoSFupdn(aWeight).second);
            histo1d_["4P0F_CR_llll_invM_trigUp"]->Fill(m4l, aWeight*systSFratio(trigSyst).first);
            histo1d_["4P0F_CR_llll_invM_trigDn"]->Fill(m4l, aWeight*systSFratio(trigSyst).second);
            histo1d_["4P0F_CR_llll_invM_recoUp"]->Fill(m4l, aWeight*systSFratio(recoSyst).first);
            histo1d_["4P0F_CR_llll_invM_recoDn"]->Fill(m4l, aWeight*systSFratio(recoSyst).second);
            histo1d_["4P0F_CR_llll_invM_PUrwgtUp"]->Fill(m4l, aWeight*purwgtUp/purwgtNo);
            histo1d_["4P0F_CR_llll_invM_PUrwgtDn"]->Fill(m4l, aWeight*purwgtDn/purwgtNo);
            histo1d_["4P0F_CR_llll_invM_prefireUp"]->Fill(m4l, aWeight*prefireUp/prefireNo);
            histo1d_["4P0F_CR_llll_invM_prefireDn"]->Fill(m4l, aWeight*prefireDn/prefireNo);
          }

          if ( m4l > 50. /*&& m4l < 500. (unblinded)*/ && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
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

            auto invMpair = pairByInvM(lvecM1,lvecM2,lvecM3,lvecM4);
            histo1d_["4P0F_CR_llHiM_invM"]->Fill(invMpair.first, aWeight);
            histo1d_["4P0F_CR_llLoM_invM"]->Fill(invMpair.second, aWeight);

            histo1d_["4P0F_CR_pt_P1"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["4P0F_CR_pt_P2"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["4P0F_CR_pt_P3"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["4P0F_CR_pt_P4"]->Fill(allHighPtMuons.at(3)->tunePMuonBestTrack()->pt(), aWeight);
          } // CR
        } // paired
      }

      break;
    case 3:
      // 3P0F - FF numerator
      if ( nonHighPtMuonsVLiso.size()==0 ) {
        for (unsigned idx = 0; idx < allHighPtMuons.size(); idx++) {
          if ( idx!=0 && allHighPtMuons.at(2)->tunePMuonBestTrack()->pt() < ptThresTrig_ )
            continue;

          const auto& first = allHighPtMuons.at(idx);
          const auto& second = allHighPtMuons.at( (idx+1)%allHighPtMuons.size() );
          const auto& probe = allHighPtMuons.at( (idx+2)%allHighPtMuons.size() );

          auto lvecM1 = lvecFromTuneP(first);
          auto lvecM2 = lvecFromTuneP(second);
          auto lvecProbe = lvecFromTuneP(probe);

          const auto lvecll = lvecM1 + lvecM2;
          const double mll = lvecll.M();
          const double dr2ll = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());

          const auto tpMET = metHandle->at(0).p4() - lvecM1 - lvecM2 - lvecProbe
                           + first->p4() + second->p4() + probe->p4();

          bool isClose = std::abs(reco::deltaPhi(tpMET.phi(),lvecM1.phi())) < 0.3;
          isClose = isClose || std::abs(reco::deltaPhi(tpMET.phi(),lvecM2.phi())) < 0.3;
          isClose = isClose || std::abs(reco::deltaPhi(tpMET.phi(),lvecProbe.phi())) < 0.3;
          bool isMergedMuSR = tpMET.pt() > 50. && isClose;

          if ( mll > 84.19 && mll < 98.19 && first->charge()*second->charge() < 0 && !isMergedMuSR ) {
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

            const double dr2_1 = reco::deltaR2(probe->tunePMuonBestTrack()->eta(),probe->tunePMuonBestTrack()->phi(),
                                               first->tunePMuonBestTrack()->eta(),first->tunePMuonBestTrack()->phi());
            const double dr2_2 = reco::deltaR2(probe->tunePMuonBestTrack()->eta(),probe->tunePMuonBestTrack()->phi(),
                                               second->tunePMuonBestTrack()->eta(),second->tunePMuonBestTrack()->phi());
            const auto lvecClosest = dr2_1 < dr2_2 ? lvecM1 : lvecM2;
            const double mindr2 = std::min( dr2_1, dr2_2 );

            histo1d_["3P0F_dr"]->Fill( std::sqrt(mindr2), aWeight );

            if ( std::abs( probe->tunePMuonBestTrack()->eta() ) < 1.2 ) {
              histo1d_["3P0F_pt_MB"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_dr_MB"]->Fill( std::sqrt(mindr2), aWeight );
            } else {
              histo1d_["3P0F_pt_ME"]->Fill(probe->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P0F_dr_ME"]->Fill( std::sqrt(mindr2), aWeight );
            }

            numerPt_ = probe->tunePMuonBestTrack()->pt();
            numerDr_ = std::sqrt(mindr2);
            numerWgt_ = aWeight;
            numerInvM_ = (lvecClosest + lvecProbe).M();
            numerIsMB_ = static_cast<int>( std::abs( probe->tunePMuonBestTrack()->eta() ) < 1.2 );
            numerTree_->Fill();
          } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
        } // Z candidate combination
      } // 3P0F

      // 3P1F
      if ( nonHighPtMuonsVLiso.size()==1 ) {
        histo1d_["3P1F_P1_pt"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
        histo1d_["3P1F_P2_pt"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);
        histo1d_["3P1F_P3_pt"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_P3_eta"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_P3_phi"]->Fill(allHighPtMuons.at(2)->tunePMuonBestTrack()->phi(), aWeight);

        histo1d_["3P1F_F1_pt"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
        histo1d_["3P1F_F1_eta"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(), aWeight);
        histo1d_["3P1F_F1_phi"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(), aWeight);

        std::vector<pat::MuonRef> allMuons(allHighPtMuons);
        allMuons.insert( allMuons.end(), nonHighPtMuonsVLiso.begin(), nonHighPtMuonsVLiso.end() );
        std::pair<pat::MuonRef,pat::MuonRef> pair1, pair2;

        bool paired = PairingHelper::pair4M(allMuons,pair1,pair2);

        if (paired && !checkTrackerMuPair(pair1,false) && !checkTrackerMuPair(pair2,false)) {
          auto lvecM1 = lvecFromTuneP(pair1.first);
          auto lvecM2 = lvecFromTuneP(pair1.second);
          auto lvecM3 = lvecFromTuneP(pair2.first);
          auto lvecM4 = lvecFromTuneP(pair2.second);
          auto lvecM1alt = lvecM1*mucorrHelper_.altScale(pair1.first,muScaleBias_,isMC_);
          auto lvecM2alt = lvecM2*mucorrHelper_.altScale(pair1.second,muScaleBias_,isMC_);
          auto lvecM3alt = lvecM3*mucorrHelper_.altScale(pair2.first,muScaleBias_,isMC_);
          auto lvecM4alt = lvecM4*mucorrHelper_.altScale(pair2.second,muScaleBias_,isMC_);

          // if (isMC_) {
          //   edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
          //   iEvent.getByToken(genptcToken_, genptcHandle);
          //
          //   lvecM1 *= mucorrHelper_.nominalMC(pair1.first,genptcHandle);
          //   lvecM2 *= mucorrHelper_.nominalMC(pair1.second,genptcHandle);
          //   lvecM3 *= mucorrHelper_.nominalMC(pair2.first,genptcHandle);
          //   lvecM4 *= mucorrHelper_.nominalMC(pair2.second,genptcHandle);
          // } else {
          //   lvecM1 *= mucorrHelper_.nominalData(pair1.first);
          //   lvecM2 *= mucorrHelper_.nominalData(pair1.second);
          //   lvecM3 *= mucorrHelper_.nominalData(pair2.first);
          //   lvecM4 *= mucorrHelper_.nominalData(pair2.second);
          // }

          const auto lvecA1 = lvecM1 + lvecM2;
          const auto lvecA2 = lvecM3 + lvecM4;
          const auto lvec4l = lvecA1 + lvecA2;
          const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
          const double dr2ll2 = reco::deltaR2(lvecM3.eta(),lvecM3.phi(),lvecM4.eta(),lvecM4.phi());
          const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
          const double m4l = lvec4l.M();
          const double m4lAlt = (lvecM1alt+lvecM2alt+lvecM3alt+lvecM4alt).M();

          const auto& fakePair = ( nonHighPtMuonsVLiso.front()==pair1.first || nonHighPtMuonsVLiso.front()==pair1.second ) ? pair1 : pair2;
          const auto& passPartner = ( nonHighPtMuonsVLiso.front()==fakePair.first ) ? fakePair.second : fakePair.first;
          const double drFake = reco::deltaR(passPartner->tunePMuonBestTrack()->eta(),passPartner->tunePMuonBestTrack()->phi(),
                                             nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(),nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi());

          const double ff = std::max( valFF(drFake,nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt()), 0.);
          const double ci = ciFF(drFake,nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt());

          if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. && drFake > 0.3 ) {
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
            histo1d_["3P1F_CR_llll_invM_altMuScale_xFF"]->Fill(m4lAlt, aWeight*ff);
            histo1d_["3P1F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff+ci));
            histo1d_["3P1F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff-std::min(ff,ci)));
          }

          if ( m4l > 50. /*&& m4l < 500. (unblinded)*/ && lvecA1.M() > 1. && lvecA2.M() > 1. && drFake > 0.3 ) {
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
            histo1d_["3P1F_CR_pt_F1"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);

            histo1d_["3P1F_CR_dr"]->Fill(drFake, aWeight);

            if ( std::abs( nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta() ) < 1.2 ) {
              histo1d_["3P1F_CR_pt_MB"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P1F_CR_dr_MB"]->Fill(drFake, aWeight);
            } else {
              histo1d_["3P1F_CR_pt_ME"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["3P1F_CR_dr_ME"]->Fill(drFake, aWeight);
            }

            histo1d_["3P1F_CR_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff);

            histo1d_["3P1F_CR_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff);
            histo1d_["3P1F_CR_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff);
            histo1d_["3P1F_CR_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff);
            histo1d_["3P1F_CR_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff);

            auto invMpair = pairByInvM(lvecM1,lvecM2,lvecM3,lvecM4);
            histo1d_["3P1F_CR_llHiM_invM_xFF"]->Fill(invMpair.first, aWeight*ff);
            histo1d_["3P1F_CR_llLoM_invM_xFF"]->Fill(invMpair.second, aWeight*ff);

            auto lvecFake = lvecFromTuneP(nonHighPtMuonsVLiso.front());
            auto lvecPartner = lvecFromTuneP(passPartner);
            const bool globHighPt = std::find( highPtMuons.begin(), highPtMuons.end(), nonHighPtMuonsVLiso.front() ) != highPtMuons.end();
            const bool trkHighPt = std::find( highPtTrackerMuons.begin(), highPtTrackerMuons.end(), nonHighPtMuonsVLiso.front() ) != highPtTrackerMuons.end();

            pt_3P1F_ = nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt();
            dr_3P1F_ = drFake;
            iso_3P1F_ = nonHighPtIsos.at(nonHighPtMuonsVLiso.front());
            mll_3P1F_ = (lvecFake + lvecPartner).M();
            m4l_3P1F_ = m4l;
            wgt_3P1F_ = aWeight;
            passHighPt_3P1F_ = static_cast<int>( globHighPt || trkHighPt );
            tree_3P1F_->Fill();
          } // CR
        } // paired
      } // 3P1F

      break;
    case 2:
      // 2P1F - FF denominator
      if ( nonHighPtMuonsVLiso.size()==1 ) {
        auto lvecM1 = lvecFromTuneP(allHighPtMuons.front());
        auto lvecM2 = lvecFromTuneP(allHighPtMuons.at(1));
        auto lvecProbe = lvecFromTuneP(nonHighPtMuonsVLiso.front());

        const auto lvecll = lvecM1 + lvecM2;
        const double mll = lvecll.M();
        const double dr2ll = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());

        const auto tpMET = metHandle->at(0).p4() - lvecM1 - lvecM2 - lvecProbe
                           + allHighPtMuons.front()->p4() + allHighPtMuons.at(1)->p4() + nonHighPtMuonsVLiso.front()->p4();

        bool isClose = std::abs(reco::deltaPhi(tpMET.phi(),lvecM1.phi())) < 0.3;
        isClose = isClose || std::abs(reco::deltaPhi(tpMET.phi(),lvecM2.phi())) < 0.3;
        isClose = isClose || std::abs(reco::deltaPhi(tpMET.phi(),lvecProbe.phi())) < 0.3;
        bool isMergedMuSR = tpMET.pt() > 50. && isClose;

        if ( mll > 84.19 && mll < 98.19 && allHighPtMuons.front()->charge()*allHighPtMuons.at(1)->charge() < 0 && !isMergedMuSR ) {
          histo1d_["2P1F_P1_pt"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["2P1F_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["2P1F_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
          histo1d_["2P1F_P2_pt"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["2P1F_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["2P1F_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);

          histo1d_["2P1F_F1_pt"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
          histo1d_["2P1F_F1_eta"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(), aWeight);
          histo1d_["2P1F_F1_phi"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(), aWeight);

          histo1d_["2P1F_ll_invM"]->Fill(mll, aWeight);
          histo1d_["2P1F_ll_pt"]->Fill(lvecll.pt(), aWeight);
          histo1d_["2P1F_ll_dr"]->Fill(std::sqrt(dr2ll), aWeight);

          const double dr2_1 = reco::deltaR2(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(),nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(),
                                             allHighPtMuons.front()->tunePMuonBestTrack()->eta(),allHighPtMuons.front()->tunePMuonBestTrack()->phi());
          const double dr2_2 = reco::deltaR2(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(),nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(),
                                             allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(),allHighPtMuons.at(1)->tunePMuonBestTrack()->phi());
          const auto lvecClosest = dr2_1 < dr2_2 ? lvecM1 : lvecM2;
          const double mindr2 = std::min( dr2_1, dr2_2 );

          histo1d_["2P1F_dr"]->Fill(std::sqrt(mindr2), aWeight);

          if ( std::abs( nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta() ) < 1.2 ) {
            histo1d_["2P1F_pt_MB"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P1F_dr_MB"]->Fill(std::sqrt(mindr2), aWeight);
          } else {
            histo1d_["2P1F_pt_ME"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
            histo1d_["2P1F_dr_ME"]->Fill(std::sqrt(mindr2), aWeight);
          }

          const bool globHighPt = std::find( highPtMuons.begin(), highPtMuons.end(), nonHighPtMuonsVLiso.front() ) != highPtMuons.end();
          const bool trkHighPt = std::find( highPtTrackerMuons.begin(), highPtTrackerMuons.end(), nonHighPtMuonsVLiso.front() ) != highPtTrackerMuons.end();

          denomPt_ = nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt();
          denomDr_ = std::sqrt(mindr2);
          denomIso_ = nonHighPtIsos.at(nonHighPtMuonsVLiso.front());
          denomWgt_ = aWeight;
          denomInvM_ = (lvecClosest + lvecProbe).M();
          denomIsMB_ = static_cast<int>( std::abs( nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta() ) < 1.2 );
          denomPassHighPt_ = static_cast<int>(globHighPt || trkHighPt);
          denomTree_->Fill();
        } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
      } // 2P1F

      // 2P2F
      if ( nonHighPtMuonsVLiso.size()==2 ) {
        std::vector<pat::MuonRef> allMuons(allHighPtMuons);
        allMuons.insert( allMuons.end(), nonHighPtMuonsVLiso.begin(), nonHighPtMuonsVLiso.end() );
        std::pair<pat::MuonRef,pat::MuonRef> pair1, pair2;

        bool paired = PairingHelper::pair4M(allMuons,pair1,pair2);

        if (paired && !checkTrackerMuPair(pair1,false) && !checkTrackerMuPair(pair2,false)) {
          bool pairedFakes = ( ( pair1.first==nonHighPtMuonsVLiso.front() && pair1.second==nonHighPtMuonsVLiso.at(1) ) ||
                               ( pair1.second==nonHighPtMuonsVLiso.front() && pair1.first==nonHighPtMuonsVLiso.at(1) ) ||
                               ( pair2.first==nonHighPtMuonsVLiso.front() && pair2.second==nonHighPtMuonsVLiso.at(1) ) ||
                               ( pair2.second==nonHighPtMuonsVLiso.front() && pair2.first==nonHighPtMuonsVLiso.at(1) ) );
          const auto firstFakePair = ( pair1.first==nonHighPtMuonsVLiso.front() || pair1.second==nonHighPtMuonsVLiso.front() )
                                     ? pair1 : pair2;
          const auto secondFakePair = ( pair1.first==nonHighPtMuonsVLiso.at(1) || pair1.second==nonHighPtMuonsVLiso.at(1) )
                                     ? pair1 : pair2;

          if (true) { // !pairedFakes
            auto lvecM1 = lvecFromTuneP(pair1.first);
            auto lvecM2 = lvecFromTuneP(pair1.second);
            auto lvecM3 = lvecFromTuneP(pair2.first);
            auto lvecM4 = lvecFromTuneP(pair2.second);
            auto lvecM1alt = lvecM1*mucorrHelper_.altScale(pair1.first,muScaleBias_,isMC_);
            auto lvecM2alt = lvecM2*mucorrHelper_.altScale(pair1.second,muScaleBias_,isMC_);
            auto lvecM3alt = lvecM3*mucorrHelper_.altScale(pair2.first,muScaleBias_,isMC_);
            auto lvecM4alt = lvecM4*mucorrHelper_.altScale(pair2.second,muScaleBias_,isMC_);

            // if (isMC_) {
            //   edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
            //   iEvent.getByToken(genptcToken_, genptcHandle);
            //
            //   lvecM1 *= mucorrHelper_.nominalMC(pair1.first,genptcHandle);
            //   lvecM2 *= mucorrHelper_.nominalMC(pair1.second,genptcHandle);
            //   lvecM3 *= mucorrHelper_.nominalMC(pair2.first,genptcHandle);
            //   lvecM4 *= mucorrHelper_.nominalMC(pair2.second,genptcHandle);
            // } else {
            //   lvecM1 *= mucorrHelper_.nominalData(pair1.first);
            //   lvecM2 *= mucorrHelper_.nominalData(pair1.second);
            //   lvecM3 *= mucorrHelper_.nominalData(pair2.first);
            //   lvecM4 *= mucorrHelper_.nominalData(pair2.second);
            // }

            const auto lvecA1 = lvecM1 + lvecM2;
            const auto lvecA2 = lvecM3 + lvecM4;
            const auto lvec4l = lvecA1 + lvecA2;
            const double dr2ll1 = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());
            const double dr2ll2 = reco::deltaR2(lvecM3.eta(),lvecM3.phi(),lvecM4.eta(),lvecM4.phi());
            const double dr2A12 = reco::deltaR2(lvecA1.eta(),lvecA1.phi(),lvecA2.eta(),lvecA2.phi());
            const double m4l = lvec4l.M();
            const double m4lAlt = (lvecM1alt+lvecM2alt+lvecM3alt+lvecM4alt).M();

            const double drFake1 = reco::deltaR( firstFakePair.first->tunePMuonBestTrack()->eta(), firstFakePair.first->tunePMuonBestTrack()->phi(),
                                                 firstFakePair.second->tunePMuonBestTrack()->eta(), firstFakePair.second->tunePMuonBestTrack()->phi() );
            const double drFake2 = reco::deltaR( secondFakePair.first->tunePMuonBestTrack()->eta(), secondFakePair.first->tunePMuonBestTrack()->phi(),
                                                 secondFakePair.second->tunePMuonBestTrack()->eta(), secondFakePair.second->tunePMuonBestTrack()->phi() );

            const double ff1 = std::max( valFF(drFake1,nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt()), 0.);
            const double ci1 = ciFF(drFake1,nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt());

            const double ff2 = std::max( valFF(drFake2,nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt()), 0.);
            const double ci2 = ciFF(drFake2,nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt());

            if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. && drFake1 > 0.3 && drFake2 > 0.3 ) {
              histo1d_["2P2F_P1_pt"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_P1_eta"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["2P2F_P1_phi"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->phi(), aWeight);
              histo1d_["2P2F_P2_pt"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_P2_eta"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["2P2F_P2_phi"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->phi(), aWeight);

              histo1d_["2P2F_F1_pt"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_F1_eta"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["2P2F_F1_phi"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->phi(), aWeight);
              histo1d_["2P2F_F2_pt"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_F2_eta"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->eta(), aWeight);
              histo1d_["2P2F_F2_phi"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->phi(), aWeight);

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
              histo1d_["2P2F_CR_llll_invM_altMuScale_xFF"]->Fill(m4lAlt, aWeight*ff1);
              histo1d_["2P2F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff1+ci1));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff1-std::min(ff1,ci1)));
              histo1d_["2P2F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff2);
              histo1d_["2P2F_CR_llll_invM_altMuScale_xFF"]->Fill(m4lAlt, aWeight*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff2+ci2));
              histo1d_["2P2F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff2-std::min(ff2,ci2)));
              histo1d_["2P2F_CR_llll_invM_xFF2"]->Fill(m4l, aWeight*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_altMuScale_xFF2"]->Fill(m4lAlt, aWeight*ff1*ff2);
              histo1d_["2P2F_CR_llll_invM_xFF2_ffUp"]->Fill(m4l, aWeight*(ff1+ci1)*(ff2+ci2));
              histo1d_["2P2F_CR_llll_invM_xFF2_ffDn"]->Fill(m4l, aWeight*(ff1-std::min(ff1,ci1))*(ff2-std::min(ff2,ci2)));
            }

            if ( m4l > 50. /*&& m4l < 500. (unblinded)*/ && lvecA1.M() > 1. && lvecA2.M() > 1. && drFake1 > 0.3 && drFake2 > 0.3 ) {
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

              auto invMpair = pairByInvM(lvecM1,lvecM2,lvecM3,lvecM4);
              histo1d_["2P2F_CR_llHiM_invM_xFF"]->Fill(invMpair.first, aWeight*ff1);
              histo1d_["2P2F_CR_llHiM_invM_xFF"]->Fill(invMpair.first, aWeight*ff2);
              histo1d_["2P2F_CR_llLoM_invM_xFF"]->Fill(invMpair.second, aWeight*ff1);
              histo1d_["2P2F_CR_llLoM_invM_xFF"]->Fill(invMpair.second, aWeight*ff2);
              histo1d_["2P2F_CR_llHiM_invM_xFF2"]->Fill(invMpair.first, aWeight*ff1*ff2);
              histo1d_["2P2F_CR_llLoM_invM_xFF2"]->Fill(invMpair.second, aWeight*ff1*ff2);

              if ( std::abs( nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->eta() ) < 1.2 ) {
                histo1d_["2P2F_CR_pt_MB"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
                histo1d_["2P2F_CR_dr_MB"]->Fill(drFake1, aWeight);
                histo1d_["2P2F_CR_pt_MB_xFF"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight*ff1);
                histo1d_["2P2F_CR_dr_MB_xFF"]->Fill(drFake1, aWeight*ff1);
              } else {
                histo1d_["2P2F_CR_pt_ME"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
                histo1d_["2P2F_CR_dr_ME"]->Fill(drFake1, aWeight);
                histo1d_["2P2F_CR_pt_ME_xFF"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight*ff1);
                histo1d_["2P2F_CR_dr_ME_xFF"]->Fill(drFake1, aWeight*ff1);
              }

              if ( std::abs( nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->eta() ) < 1.2 ) {
                histo1d_["2P2F_CR_pt_MB"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt(), aWeight);
                histo1d_["2P2F_CR_dr_MB"]->Fill(drFake2, aWeight);
                histo1d_["2P2F_CR_pt_MB_xFF"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt(), aWeight*ff2);
                histo1d_["2P2F_CR_dr_MB_xFF"]->Fill(drFake2, aWeight*ff2);
              } else {
                histo1d_["2P2F_CR_pt_ME"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt(), aWeight);
                histo1d_["2P2F_CR_dr_ME"]->Fill(drFake2, aWeight);
                histo1d_["2P2F_CR_pt_ME_xFF"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt(), aWeight*ff2);
                histo1d_["2P2F_CR_dr_ME_xFF"]->Fill(drFake2, aWeight*ff2);
              }

              histo1d_["2P2F_CR_pt_P1"]->Fill(allHighPtMuons.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_pt_P2"]->Fill(allHighPtMuons.at(1)->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_pt_F1"]->Fill(nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt(), aWeight);
              histo1d_["2P2F_CR_pt_F2"]->Fill(nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt(), aWeight);

              const auto lvecFirstFakeM1 = lvecFromTuneP(firstFakePair.first);
              const auto lvecFirstFakeM2 = lvecFromTuneP(firstFakePair.second);
              const auto lvecSecondFakeM1 = lvecFromTuneP(secondFakePair.first);
              const auto lvecSecondFakeM2 = lvecFromTuneP(secondFakePair.second);

              const bool globHighPt1 = std::find( highPtMuons.begin(), highPtMuons.end(), nonHighPtMuonsVLiso.front() ) != highPtMuons.end();
              const bool trkHighPt1 = std::find( highPtTrackerMuons.begin(), highPtTrackerMuons.end(), nonHighPtMuonsVLiso.front() ) != highPtTrackerMuons.end();
              const bool globHighPt2 = std::find( highPtMuons.begin(), highPtMuons.end(), nonHighPtMuonsVLiso.at(1) ) != highPtMuons.end();
              const bool trkHighPt2 = std::find( highPtTrackerMuons.begin(), highPtTrackerMuons.end(), nonHighPtMuonsVLiso.at(1) ) != highPtTrackerMuons.end();

              pt1_2P2F_ = nonHighPtMuonsVLiso.front()->tunePMuonBestTrack()->pt();
              pt2_2P2F_ = nonHighPtMuonsVLiso.at(1)->tunePMuonBestTrack()->pt();
              dr1_2P2F_ = drFake1;
              dr2_2P2F_ = drFake2;
              iso1_2P2F_ = nonHighPtIsos.at(nonHighPtMuonsVLiso.front());
              iso2_2P2F_ = nonHighPtIsos.at(nonHighPtMuonsVLiso.at(1));
              mll1_2P2F_ = (lvecFirstFakeM1 + lvecFirstFakeM2).M();
              mll2_2P2F_ = (lvecSecondFakeM1 + lvecSecondFakeM2).M();
              m4l_2P2F_ = m4l;
              wgt_2P2F_ = aWeight;
              pairedFakes_2P2F_ = static_cast<int>(pairedFakes);
              passHighPt_2P2F_ = static_cast<int>( (globHighPt1 || trkHighPt1) && (globHighPt2 || trkHighPt2) );
              tree_2P2F_->Fill();
            } // CR

            if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. && pairedFakes && drFake1 < 0.3 ) {
              const auto lvecllFake = firstFakePair.first->polarP4() + firstFakePair.second->polarP4();
              const double ptll = lvecllFake.Pt();
              const double valFFdr03 = ffdr03_->Eval(ptll);
              const double cidr03 = ciFFdr03(ptll);

              histo1d_["2P2F_CRdr03_llll_invM"]->Fill(m4l, aWeight);
              histo1d_["2P2F_CRdr03_llll_invM_xFF"]->Fill(m4l, aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_llll_invM_altMuScale_xFF"]->Fill(m4lAlt, aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(valFFdr03+cidr03));
              histo1d_["2P2F_CRdr03_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*std::max(0.,valFFdr03-cidr03));
              histo1d_["2P2F_CRdr03_llll_pt"]->Fill(lvec4l.pt(), aWeight);
              histo1d_["2P2F_CRdr03_ll1ll2_dr"]->Fill(std::sqrt(dr2A12), aWeight);

              histo1d_["2P2F_CRdr03_ll1_invM"]->Fill(lvecA1.M(), aWeight);
              histo1d_["2P2F_CRdr03_ll1_pt"]->Fill(lvecA1.pt(), aWeight);
              histo1d_["2P2F_CRdr03_ll1_dr"]->Fill(std::sqrt(dr2ll1), aWeight);
              histo1d_["2P2F_CRdr03_ll2_invM"]->Fill(lvecA2.M(), aWeight);
              histo1d_["2P2F_CRdr03_ll2_pt"]->Fill(lvecA2.pt(), aWeight);
              histo1d_["2P2F_CRdr03_ll2_dr"]->Fill(std::sqrt(dr2ll2), aWeight);

              histo1d_["2P2F_CRdr03_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*valFFdr03);

              auto invMpair = pairByInvM(lvecM1,lvecM2,lvecM3,lvecM4);
              histo1d_["2P2F_CRdr03_llHiM_invM_xFF"]->Fill(invMpair.first, aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_llLoM_invM_xFF"]->Fill(invMpair.second, aWeight*valFFdr03);

              const auto m1 = firstFakePair.first;
              const auto m2 = firstFakePair.second;
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
              tree_2P2F_dr03_->Fill();
            } // nonisolated CR
          } // !pairedFakes
        } // paired
      } // 2P2F

      break;
  } // switch allHighPtMuons.size()
}

DEFINE_FWK_MODULE(ResolvedMuCRanalyzer);
