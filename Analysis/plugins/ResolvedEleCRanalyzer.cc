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
#include "ZprimeTo4l/Analysis/interface/ElectronSystematicsHelper.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"
#include "TFitResult.h"

#include "correction.h"

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
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<double> prefweightUp_token_;
  const edm::EDGetTokenT<double> prefweightDn_token_;

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> modifiedTrkIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<double> rhoToken_;

  const edm::FileInPath FFpath_;
  const edm::FileInPath FFdr03path_;

  std::unique_ptr<correction::CorrectionSet> purwgt_;
  const std::string puname_;

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigFilters_;
  const std::vector<std::string> trigUnseededFilters_;
  const double etThres1_;
  const double etThres2_;
  const double ffSystCL95_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  std::unique_ptr<TFile> FFfile_;
  // TH1* ff_dr_;
  TF1* ff_pt_;
  TFitResultPtr ffFitResult_;

  std::unique_ptr<TFile> FFdr03file_;
  TF1* ffdr03_;
  TFitResultPtr ffdr03FitResult_;

  ElectronSystematicsHelper systHelperEle_;

  TTree* numerTree_ = nullptr;
  float numerPt_ = -1.;
  float numerDr_ = -1.;
  float numerWgt_ = 0.;
  float numerInvM_ = -1.;
  int numerIsEB_ = -1;

  TTree* denomTree_ = nullptr;
  float denomPt_ = -1.;
  float denomDr_ = -1.;
  float denomWgt_ = 0.;
  float denomDPerpIn_ = std::numeric_limits<float>::max();
  float denomDEtaInSeed_ = std::numeric_limits<float>::max();
  float denomMvaScore_ = -1.;
  int denomMvaCategory_ = -1;
  int denomPassMva_ = -1;
  float denomIso_ = -1.;
  float denomDxy_ = std::numeric_limits<float>::max();
  float denomEcalIso_ = -1.;
  float denomHcalIso_ = -1.;
  float denomHoE_ = -1.;
  float denomRho_ = -1.;
  float denomE_ = -1.;
  int denomMIH_ = -1;
  float denomInvM_ = -1.;
  int denomIsEB_ = -1;

  TTree* tree_3P1F_ = nullptr;
  float pt_3P1F_ = -1.;
  float dr_3P1F_ = -1.;
  float dPerpIn_3P1F_ = std::numeric_limits<float>::max();
  float dEtaInSeed_3P1F_ = std::numeric_limits<float>::max();
  float mvaScore_3P1F_ = -1.;
  int mvaCategory_3P1F_ = -1;
  int passMva_3P1F_ = -1;
  float iso_3P1F_ = -1.;
  float dxy_3P1F_ = std::numeric_limits<float>::max();
  float ecalIso_3P1F_ = -1.;
  float hcalIso_3P1F_ = -1.;
  float HoE_3P1F_ = -1.;
  float rho_3P1F_ = -1.;
  float enSC_3P1F_ = -1.;
  int mih_3P1F_ = -1;
  float mll_3P1F_ = -1.;
  float m4l_3P1F_ = -1.;
  float wgt_3P1F_ = 0.;

  TTree* tree_2P2F_ = nullptr;
  float pt1_2P2F_ = -1.;
  float pt2_2P2F_ = -1.;
  float dr1_2P2F_ = -1.;
  float dr2_2P2F_ = -1.;
  float dPerpIn1_2P2F_ = std::numeric_limits<float>::max();
  float dPerpIn2_2P2F_ = std::numeric_limits<float>::max();
  float dEtaInSeed1_2P2F_ = std::numeric_limits<float>::max();
  float dEtaInSeed2_2P2F_ = std::numeric_limits<float>::max();
  float mvaScore1_2P2F_ = -1.;
  float mvaScore2_2P2F_ = -1.;
  int mvaCategory1_2P2F_ = -1;
  int mvaCategory2_2P2F_ = -1;
  int passMva1_2P2F_ = -1;
  int passMva2_2P2F_ = -1;
  float iso1_2P2F_ = -1.;
  float iso2_2P2F_ = -1.;
  float dxy1_2P2F_ = std::numeric_limits<float>::max();
  float dxy2_2P2F_ = std::numeric_limits<float>::max();
  float ecalIso1_2P2F_ = -1.;
  float ecalIso2_2P2F_ = -1.;
  float hcalIso1_2P2F_ = -1.;
  float hcalIso2_2P2F_ = -1.;
  float HoE1_2P2F_ = -1.;
  float HoE2_2P2F_ = -1.;
  float enSC1_2P2F_ = -1.;
  float enSC2_2P2F_ = -1.;
  float rho_2P2F_ = -1.;
  int mih1_2P2F_ = -1;
  int mih2_2P2F_ = -1;
  float mll1_2P2F_ = -1.;
  float mll2_2P2F_ = -1.;
  float m4l_2P2F_ = -1.;
  float wgt_2P2F_ = 0.;
  int pairedFakes_2P2F_ = -1;

  TTree* tree_2P2F_dr03_ = nullptr;
  float var_weight_ = 0.;
  float var_e1Pt_ = -1.;
  float var_e1scEn_ = -1.;
  float var_e1scEta_ = std::numeric_limits<float>::max();
  float var_e1Eta_ = std::numeric_limits<float>::max();
  float var_e1Phi_ = std::numeric_limits<float>::max();
  float var_e1DPerpIn_ = std::numeric_limits<float>::max();
  float var_e1DEtaInSeed_ = std::numeric_limits<float>::max();
  float var_e1MvaScore_ = -1.;
  int var_e1MvaCategory_ = -1;
  int var_e1PassMva_ = -1;
  float var_e1TrkIso_ = -1.;
  float var_e1EcalIso_ = -1.;
  float var_e1HcalIso_ = -1.;
  float var_e1HoE_ = -1.;
  float var_e2Pt_ = -1.;
  float var_e2scEn_ = -1.;
  float var_e2scEta_ = std::numeric_limits<float>::max();
  float var_e2Eta_ = std::numeric_limits<float>::max();
  float var_e2Phi_ = std::numeric_limits<float>::max();
  float var_e2DPerpIn_ = std::numeric_limits<float>::max();
  float var_e2DEtaInSeed_ = std::numeric_limits<float>::max();
  float var_e2MvaScore_ = -1.;
  int var_e2MvaCategory_ = -1;
  int var_e2PassMva_ = -1;
  float var_e2TrkIso_ = -1.;
  float var_e2EcalIso_ = -1.;
  float var_e2HcalIso_ = -1.;
  float var_e2HoE_ = -1.;
  float var_e1e2Rho_ = -1.;
  float var_e1e2InvM_ = -1.;
  float var_e1e2DR_ = -1.;
  float var_e1e2scDR_ = -1.;
};

ResolvedEleCRanalyzer::ResolvedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
prefweightUp_token_(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
prefweightDn_token_(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
modifiedTrkIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("modifiedTrkIso"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("modifiedEcalIso"))),
dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
FFdr03path_(iConfig.getParameter<edm::FileInPath>("FFdr03path")),
purwgt_(std::move(correction::CorrectionSet::from_file((iConfig.getParameter<edm::FileInPath>("PUrwgt")).fullPath()))),
puname_(iConfig.getParameter<std::string>("PUname")),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigFilters_(iConfig.getParameter<std::vector<std::string>>("trigFilters")),
trigUnseededFilters_(iConfig.getParameter<std::vector<std::string>>("trigUnseededFilters")),
etThres1_(iConfig.getParameter<double>("etThres1")),
etThres2_(iConfig.getParameter<double>("etThres2")),
ffSystCL95_(iConfig.getParameter<double>("ffSystCL95")),
systHelperEle_(ElectronSystematicsHelper(iConfig.getParameter<edm::FileInPath>("modHeepSFpath"),
                                         iConfig.getParameter<edm::FileInPath>("recoSFpath"),
                                         iConfig.getParameter<edm::FileInPath>("trigSFpath"),
                                         iConfig.getParameter<edm::FileInPath>("trigUnseededSFpath"))) {
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

  std::vector<double> trigSF = iConfig.getParameter<std::vector<double>>("trigSF");
  std::vector<double> trigSFcl95 = iConfig.getParameter<std::vector<double>>("trigSFcl95");
  std::vector<double> trigSFupper = iConfig.getParameter<std::vector<double>>("trigSFupper");
  std::vector<std::string> trigPol = iConfig.getParameter<std::vector<std::string>>("trigPol");
  systHelperEle_.SetTrigSF(trigSF.at(0), trigSF.at(1), trigSF.at(2));
  systHelperEle_.SetTrigSFcl95(trigSFcl95.at(0), trigSFcl95.at(1), trigSFcl95.at(2));
  systHelperEle_.SetTrigSFupper(trigSFupper.at(0), trigSFupper.at(1), trigSFupper.at(2));
  systHelperEle_.SetTrigPol1(trigPol.at(0), trigPol.at(1), trigPol.at(2));

  std::vector<double> trigUnseededSF = iConfig.getParameter<std::vector<double>>("trigUnseededSF");
  std::vector<double> trigUnseededSFcl95 = iConfig.getParameter<std::vector<double>>("trigUnseededSFcl95");
  std::vector<double> trigUnseededSFupper = iConfig.getParameter<std::vector<double>>("trigUnseededSFupper");
  std::vector<std::string> trigUnseededPol = iConfig.getParameter<std::vector<std::string>>("trigUnseededPol");
  systHelperEle_.SetTrigUnseededSF(trigUnseededSF.at(0), trigUnseededSF.at(1), trigUnseededSF.at(2));
  systHelperEle_.SetTrigUnseededSFcl95(trigUnseededSFcl95.at(0), trigUnseededSFcl95.at(1), trigUnseededSFcl95.at(2));
  systHelperEle_.SetTrigUnseededSFupper(trigUnseededSFupper.at(0), trigUnseededSFupper.at(1), trigUnseededSFupper.at(2));
  systHelperEle_.SetTrigUnseededPol1(trigUnseededPol.at(0), trigUnseededPol.at(1), trigUnseededPol.at(2));
}

void ResolvedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  // ff_dr_ = static_cast<TH1*>(FFfile_->Get("all_dr_3P0F_rebin"));
  ff_pt_ = static_cast<TF1*>(FFfile_->FindObjectAny("REFFfunc"));
  ffFitResult_ = (static_cast<TH1D*>(FFfile_->Get("REFF_pt")))->Fit(ff_pt_,"RS");

  FFdr03file_ = std::make_unique<TFile>(FFdr03path_.fullPath().c_str(),"READ");
  ffdr03_ = static_cast<TF1*>(FFdr03file_->FindObjectAny("REFFdr03"));
  ffdr03FitResult_ = (static_cast<TH1D*>(FFdr03file_->Get("REFFhist")))->Fit(ffdr03_,"RS");

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
  histo1d_["3P0F_dr_EEneg_preHEM"] = fs->make<TH1D>("3P0F_dr_EEneg_preHEM","3P0F Pass min dR(ll) (EE)",128,0.,6.4);
  histo1d_["3P0F_dr_EEneg_postHEM"] = fs->make<TH1D>("3P0F_dr_EEneg_postHEM","3P0F Pass min dR(ll) (EE)",128,0.,6.4);
  histo1d_["3P0F_dr_preHEM"] = fs->make<TH1D>("3P0F_dr_preHEM","3P0F Pass min dR(ll)",128,0.,6.4);
  histo1d_["3P0F_dr_postHEM"] = fs->make<TH1D>("3P0F_dr_postHEM","3P0F Pass min dR(ll)",128,0.,6.4);
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
  histo1d_["2P1F_dr_EEneg_preHEM"] = fs->make<TH1D>("2P1F_dr_EEneg_preHEM","2P1F Fail min dR(ll) (EE)",128,0.,6.4);
  histo1d_["2P1F_dr_EEneg_postHEM"] = fs->make<TH1D>("2P1F_dr_EEneg_postHEM","2P1F Fail min dR(ll) (EE)",128,0.,6.4);
  histo1d_["2P1F_dr_preHEM"] = fs->make<TH1D>("2P1F_dr_preHEM","2P1F Fail min dR(ll)",128,0.,6.4);
  histo1d_["2P1F_dr_postHEM"] = fs->make<TH1D>("2P1F_dr_postHEM","2P1F Fail min dR(ll)",128,0.,6.4);
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

  // 2P2F noniso
  histo1d_["2P2F_CRdr03_llll_invM"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_idUp"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_idUp","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_idDn"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_idDn","2P2F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_pt"] = fs->make<TH1D>("2P2F_CRdr03_llll_pt","2P2F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1ll2_dr"] = fs->make<TH1D>("2P2F_CRdr03_ll1ll2_dr","2P2F CR dR(ll1ll2)",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll1_invM"] = fs->make<TH1D>("2P2F_CRdr03_ll1_invM","2P2F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1_pt"] = fs->make<TH1D>("2P2F_CRdr03_ll1_pt","2P2F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1_dr"] = fs->make<TH1D>("2P2F_CRdr03_ll1_dr","2P2F CR dR(ll1)",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll2_invM"] = fs->make<TH1D>("2P2F_CRdr03_ll2_invM","2P2F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_pt"] = fs->make<TH1D>("2P2F_CRdr03_ll2_pt","2P2F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_dr"] = fs->make<TH1D>("2P2F_CRdr03_ll2_dr","2P2F CR dR(ll2)",128,0.,6.4);

  histo1d_["2P2F_CRdr03_llll_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_xFF_ffUp"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF_ffUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_xFF_ffDn"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF_ffDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_xFF_idUp"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF_idUp","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_llll_invM_xFF_idDn"] = fs->make<TH1D>("2P2F_CRdr03_llll_invM_xFF_idDn","2P2F CR M(4l) x Fake factor;M [GeV];",500,0.,2500.);
  histo1d_["2P2F_CRdr03_ll1ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1ll2_dr_xFF","2P2F CR dR(ll1ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll1_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1_invM_xFF","2P2F CR M(ll1) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll1_dr_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1_dr_xFF","2P2F CR dR(ll1) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll2_invM_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll2_invM_xFF","2P2F CR M(ll2) x Fake factor;M [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_dr_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll2_dr_xFF","2P2F CR dR(ll2) x Fake factor",128,0.,6.4);
  histo1d_["2P2F_CRdr03_ll1_pt_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll1_pt_xFF","2P2F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["2P2F_CRdr03_ll2_pt_xFF"] = fs->make<TH1D>("2P2F_CRdr03_ll2_pt_xFF","2P2F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);

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
  histo1d_["4P0F_CR_llll_invM_elRecoUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_elRecoUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_elRecoDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_elRecoDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_elTrigUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_elTrigUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_elTrigDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_elTrigDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_PUrwgtUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_PUrwgtUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_PUrwgtDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_PUrwgtDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_prefireUp"] = fs->make<TH1D>("4P0F_CR_llll_invM_prefireUp","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_invM_prefireDn"] = fs->make<TH1D>("4P0F_CR_llll_invM_prefireDn","4P0F CR M(4l);M [GeV];",500,0.,2500.);
  histo1d_["4P0F_CR_llll_pt"] = fs->make<TH1D>("4P0F_CR_llll_pt","4P0F CR p_{T}(4l);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll1ll2_dr","4P0F CR dR(ll1ll2)",64,0.,6.4);
  histo1d_["4P0F_CR_ll1_invM"] = fs->make<TH1D>("4P0F_CR_ll1_invM","4P0F CR M(ll1);M [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_pt"] = fs->make<TH1D>("4P0F_CR_ll1_pt","4P0F CR p_{T}(ll1);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll1_dr"] = fs->make<TH1D>("4P0F_CR_ll1_dr","4P0F CR dR(ll1)",64,0.,6.4);
  histo1d_["4P0F_CR_ll2_invM"] = fs->make<TH1D>("4P0F_CR_ll2_invM","4P0F CR M(ll2);M [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_pt"] = fs->make<TH1D>("4P0F_CR_ll2_pt","4P0F CR p_{T}(ll2);p_{T} [GeV];",100,0.,200.);
  histo1d_["4P0F_CR_ll2_dr"] = fs->make<TH1D>("4P0F_CR_ll2_dr","4P0F CR dR(ll2)",64,0.,6.4);

  numerTree_ = fs->make<TTree>("numerTree","numerTree");
  numerTree_->Branch("pt",&numerPt_,"pt/F");
  numerTree_->Branch("dr",&numerDr_,"dr/F");
  numerTree_->Branch("wgt",&numerWgt_,"wgt/F");
  numerTree_->Branch("invM",&numerInvM_,"invM/F");
  numerTree_->Branch("isEB",&numerIsEB_,"isEB/I");

  denomTree_ = fs->make<TTree>("denomTree","denomTree");
  denomTree_->Branch("pt",&denomPt_,"pt/F");
  denomTree_->Branch("dr",&denomDr_,"dr/F");
  denomTree_->Branch("wgt",&denomWgt_,"wgt/F");
  denomTree_->Branch("dPerpIn",&denomDPerpIn_,"dPerpIn/F");
  denomTree_->Branch("dEtaInSeed",&denomDEtaInSeed_,"dEtaInSeed/F");
  denomTree_->Branch("mvaScore",&denomMvaScore_,"mvaScore/F");
  denomTree_->Branch("mvaCategory",&denomMvaCategory_,"mvaCategory/I");
  denomTree_->Branch("passMva",&denomPassMva_,"passMva/I");
  denomTree_->Branch("iso",&denomIso_,"iso/F");
  denomTree_->Branch("dxy",&denomDxy_,"dxy/F");
  denomTree_->Branch("ecalIso",&denomEcalIso_,"ecalIso/F");
  denomTree_->Branch("hcalIso",&denomHcalIso_,"hcalIso/F");
  denomTree_->Branch("HoE",&denomHoE_,"HoE/F");
  denomTree_->Branch("rho",&denomRho_,"rho/F");
  denomTree_->Branch("enSC",&denomE_,"enSC/F");
  denomTree_->Branch("mih",&denomMIH_,"mih/I");
  denomTree_->Branch("invM",&denomInvM_,"invM/F");
  denomTree_->Branch("isEB",&denomIsEB_,"isEB/I");

  tree_3P1F_ = fs->make<TTree>("tree_3P1F","tree_3P1F");
  tree_3P1F_->Branch("pt",&pt_3P1F_,"pt/F");
  tree_3P1F_->Branch("dr",&dr_3P1F_,"dr/F");
  tree_3P1F_->Branch("dPerpIn",&dPerpIn_3P1F_,"dPerpIn/F");
  tree_3P1F_->Branch("dEtaInSeed",&dEtaInSeed_3P1F_,"dEtaInSeed/F");
  tree_3P1F_->Branch("mvaScore",&mvaScore_3P1F_,"mvaScore/F");
  tree_3P1F_->Branch("mvaCategory",&mvaCategory_3P1F_,"mvaCategory/I");
  tree_3P1F_->Branch("passMva",&passMva_3P1F_,"passMva/I");
  tree_3P1F_->Branch("iso",&iso_3P1F_,"iso/F");
  tree_3P1F_->Branch("dxy",&dxy_3P1F_,"dxy/F");
  tree_3P1F_->Branch("ecalIso",&ecalIso_3P1F_,"ecalIso/F");
  tree_3P1F_->Branch("hcalIso",&hcalIso_3P1F_,"hcalIso/F");
  tree_3P1F_->Branch("HoE",&HoE_3P1F_,"HoE/F");
  tree_3P1F_->Branch("rho",&rho_3P1F_,"rho/F");
  tree_3P1F_->Branch("enSC",&enSC_3P1F_,"enSC/F");
  tree_3P1F_->Branch("mih",&mih_3P1F_,"mih/I");
  tree_3P1F_->Branch("mll",&mll_3P1F_,"mll/F");
  tree_3P1F_->Branch("m4l",&m4l_3P1F_,"m4l/F");
  tree_3P1F_->Branch("wgt",&wgt_3P1F_,"wgt/F");

  tree_2P2F_ = fs->make<TTree>("tree_2P2F","tree_2P2F");
  tree_2P2F_->Branch("pt1",&pt1_2P2F_,"pt1/F");
  tree_2P2F_->Branch("pt2",&pt2_2P2F_,"pt2/F");
  tree_2P2F_->Branch("dr1",&dr1_2P2F_,"dr1/F");
  tree_2P2F_->Branch("dr2",&dr2_2P2F_,"dr2/F");
  tree_2P2F_->Branch("dPerpIn1",&dPerpIn1_2P2F_,"dPerpIn1/F");
  tree_2P2F_->Branch("dPerpIn2",&dPerpIn2_2P2F_,"dPerpIn2/F");
  tree_2P2F_->Branch("dEtaInSeed1",&dEtaInSeed1_2P2F_,"dEtaInSeed1/F");
  tree_2P2F_->Branch("dEtaInSeed2",&dEtaInSeed2_2P2F_,"dEtaInSeed2/F");
  tree_2P2F_->Branch("mvaScore1",&mvaScore1_2P2F_,"mvaScore1/F");
  tree_2P2F_->Branch("mvaScore2",&mvaScore2_2P2F_,"mvaScore2/F");
  tree_2P2F_->Branch("mvaCategory1",&mvaCategory1_2P2F_,"mvaCategory1/I");
  tree_2P2F_->Branch("mvaCategory2",&mvaCategory2_2P2F_,"mvaCategory2/I");
  tree_2P2F_->Branch("passMva1",&passMva1_2P2F_,"passMva1/I");
  tree_2P2F_->Branch("passMva2",&passMva2_2P2F_,"passMva2/I");
  tree_2P2F_->Branch("iso1",&iso1_2P2F_,"iso1/F");
  tree_2P2F_->Branch("iso2",&iso2_2P2F_,"iso2/F");
  tree_2P2F_->Branch("dxy1",&dxy1_2P2F_,"dxy1/F");
  tree_2P2F_->Branch("dxy2",&dxy2_2P2F_,"dxy2/F");
  tree_2P2F_->Branch("ecalIso1",&ecalIso1_2P2F_,"ecalIso1/F");
  tree_2P2F_->Branch("ecalIso2",&ecalIso2_2P2F_,"ecalIso2/F");
  tree_2P2F_->Branch("hcalIso1",&hcalIso1_2P2F_,"hcalIso1/F");
  tree_2P2F_->Branch("hcalIso2",&hcalIso2_2P2F_,"hcalIso2/F");
  tree_2P2F_->Branch("HoE1",&HoE1_2P2F_,"HoE1/F");
  tree_2P2F_->Branch("HoE2",&HoE2_2P2F_,"HoE2/F");
  tree_2P2F_->Branch("enSC1",&enSC1_2P2F_,"enSC1/F");
  tree_2P2F_->Branch("enSC2",&enSC2_2P2F_,"enSC2/F");
  tree_2P2F_->Branch("rho",&rho_2P2F_,"rho/F");
  tree_2P2F_->Branch("mih1",&mih1_2P2F_,"mih1/I");
  tree_2P2F_->Branch("mih2",&mih2_2P2F_,"mih2/I");
  tree_2P2F_->Branch("mll1",&mll1_2P2F_,"mll1/F");
  tree_2P2F_->Branch("mll2",&mll2_2P2F_,"mll2/F");
  tree_2P2F_->Branch("m4l",&m4l_2P2F_,"m4l/F");
  tree_2P2F_->Branch("wgt",&wgt_2P2F_,"wgt/F");
  tree_2P2F_->Branch("pairedFakes",&pairedFakes_2P2F_,"pairedFakes/I");

  tree_2P2F_dr03_ = fs->make<TTree>("tree_2P2F_dr03","tree_2P2F_dr03");
  tree_2P2F_dr03_->Branch("wgt",&var_weight_,"wgt/F");
  tree_2P2F_dr03_->Branch("e1Pt",&var_e1Pt_,"e1Pt/F");
  tree_2P2F_dr03_->Branch("e1scEn",&var_e1scEn_,"e1scEn/F");
  tree_2P2F_dr03_->Branch("e1scEta",&var_e1scEta_,"e1scEta/F");
  tree_2P2F_dr03_->Branch("e1Eta",&var_e1Eta_,"e1Eta/F");
  tree_2P2F_dr03_->Branch("e1Phi",&var_e1Phi_,"e1Phi/F");
  tree_2P2F_dr03_->Branch("e1DPerpIn",&var_e1DPerpIn_,"e1DPerpIn/F");
  tree_2P2F_dr03_->Branch("e1DEtaInSeed",&var_e1DEtaInSeed_,"e1DEtaInSeed/F");
  tree_2P2F_dr03_->Branch("e1MvaScore",&var_e1MvaScore_,"e1MvaScore/F");
  tree_2P2F_dr03_->Branch("e1MvaCategory",&var_e1MvaCategory_,"e1MvaCategory/I");
  tree_2P2F_dr03_->Branch("e1PassMva",&var_e1PassMva_,"e1PassMva/I");
  tree_2P2F_dr03_->Branch("e1TrkIso",&var_e1TrkIso_,"e1TrkIso/F");
  tree_2P2F_dr03_->Branch("e1EcalIso",&var_e1EcalIso_,"e1EcalIso/F");
  tree_2P2F_dr03_->Branch("e1HcalIso",&var_e1HcalIso_,"e1HcalIso/F");
  tree_2P2F_dr03_->Branch("e1HoE",&var_e1HoE_,"e1HoE/F");

  tree_2P2F_dr03_->Branch("e2Pt",&var_e2Pt_,"e2Pt/F");
  tree_2P2F_dr03_->Branch("e2scEn",&var_e2scEn_,"e2scEn/F");
  tree_2P2F_dr03_->Branch("e2scEta",&var_e2scEta_,"e2scEta/F");
  tree_2P2F_dr03_->Branch("e2Eta",&var_e2Eta_,"e2Eta/F");
  tree_2P2F_dr03_->Branch("e2Phi",&var_e2Phi_,"e2Phi/F");
  tree_2P2F_dr03_->Branch("e2DPerpIn",&var_e2DPerpIn_,"e2DPerpIn/F");
  tree_2P2F_dr03_->Branch("e2DEtaInSeed",&var_e2DEtaInSeed_,"e2DEtaInSeed/F");
  tree_2P2F_dr03_->Branch("e2MvaScore",&var_e2MvaScore_,"e2MvaScore/F");
  tree_2P2F_dr03_->Branch("e2MvaCategory",&var_e2MvaCategory_,"e2MvaCategory/I");
  tree_2P2F_dr03_->Branch("e2PassMva",&var_e2PassMva_,"e2PassMva/I");
  tree_2P2F_dr03_->Branch("e2TrkIso",&var_e2TrkIso_,"e2TrkIso/F");
  tree_2P2F_dr03_->Branch("e2EcalIso",&var_e2EcalIso_,"e2EcalIso/F");
  tree_2P2F_dr03_->Branch("e2HcalIso",&var_e2HcalIso_,"e2HcalIso/F");
  tree_2P2F_dr03_->Branch("e2HoE",&var_e2HoE_,"e2HoE/F");

  tree_2P2F_dr03_->Branch("e1e2Rho",&var_e1e2Rho_,"e1e2Rho/F");
  tree_2P2F_dr03_->Branch("e1e2InvM",&var_e1e2InvM_,"e1e2InvM/F");
  tree_2P2F_dr03_->Branch("e1e2DR",&var_e1e2DR_,"e1e2DR/F");
  tree_2P2F_dr03_->Branch("e1e2scDR",&var_e1e2scDR_,"e1e2scDR/F");
}

void ResolvedEleCRanalyzer::endJob() {
  FFfile_->Close();
}

void ResolvedEleCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

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
  histo1d_["cutflow_4E"]->Fill(0.5,aWeight);

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

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
  std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjsUnseeded;

  for (unsigned iTrig = 0; iTrig < trigObjHandle->size(); iTrig++) {
    const auto& trigObj = trigObjHandle->refAt(iTrig);
    auto trigObjInst = trigObjHandle->at(iTrig); // workaround for copy
    trigObjInst.unpackFilterLabels(iEvent, *trigResultHandle);

    for (const auto& aname : trigFilters_) {
      if (trigObjInst.hasFilterLabel(aname))
        trigObjs.push_back(trigObj);
    }

    for (const auto& aname : trigUnseededFilters_) {
      if (trigObjInst.hasFilterLabel(aname))
        trigObjsUnseeded.push_back(trigObj);
    }
  } // trigger objs

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->ptrAt(0);
  else
    return;

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<float>> modifiedTrkIsoHandle;
  iEvent.getByToken(modifiedTrkIsoToken_, modifiedTrkIsoHandle);

  edm::Handle<edm::ValueMap<float>> modifiedEcalIsoHandle;
  iEvent.getByToken(ecalIsoToken_, modifiedEcalIsoHandle);

  edm::Handle<edm::ValueMap<float>> dPerpInHandle;
  iEvent.getByToken(dPerpInToken_, dPerpInHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_, rhoHandle);

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

    if (!selected) {
      auto castEle = aEle.castTo<pat::ElectronRef>();

      bool isNonHeepEle = systHelperEle_.isNonHeepEle(castEle,(*modifiedTrkIsoHandle)[castEle],(*dPerpInHandle)[castEle]);

      if (isNonHeepEle)
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

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->tunePMuonBestTrack()->pt() < 20. || std::abs(aMuon->eta()) > 2.4 )
      continue;

    if ( muon::isHighPtMuon(*aMuon,*primaryVertex) )
      highPtMuons.push_back(aMuon.castTo<pat::MuonRef>());
    else if ( muon::isTrackerHighPtMuon(*aMuon,*primaryVertex) )
      highPtTrackerMuons.push_back(aMuon.castTo<pat::MuonRef>());
  }

  if ( !highPtMuons.empty() || !highPtTrackerMuons.empty() )
    return;

  histo1d_["cutflow_4E"]->Fill(4.5,aWeight);

  bool trigMatched1 = false;
  bool trigMatched2 = false;

  for (const auto& trigObj : trigObjs) {
    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),acceptEles.front()->eta(),acceptEles.front()->phi()) < 0.01 )
      trigMatched1 = true;
  }

  for (const auto& trigObj : trigObjsUnseeded) {
    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi()) < 0.01 )
      trigMatched2 = true;
  }

  if ( !trigMatched1 || !trigMatched2 )
    return;

  // SF
  if (isMC_) {
    aWeight *= systHelperEle_.GetTrigSF(acceptEles.front());
    aWeight *= systHelperEle_.GetTrigUnseededSF(acceptEles.at(1));

    for (const auto& aEle : acceptEles) {
      aWeight *= systHelperEle_.GetRecoSF(aEle);
      aWeight *= systHelperEle_.GetModifiedHeepSF(aEle);
    }

    for (const auto& aEle : nonHeepEles)
      aWeight *= systHelperEle_.GetRecoSF(aEle);
  }

  histo1d_["cutflow_4E"]->Fill(5.5,aWeight);

  if ( acceptEles.size()==4 )
    histo1d_["cutflow_4E"]->Fill(6.5,aWeight);

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

  auto recoSFcl95 = [this,&acceptEles,&nonHeepEles] (const double wgt, bool ignore=true) -> std::pair<double,double> {
    if (!isMC_)
      return std::make_pair(wgt,wgt);

    double wgtUp = wgt;
    double wgtDn = wgt;

    for (const auto& aEle : acceptEles) {
      double nominal = systHelperEle_.GetRecoSF(aEle);
      double err = systHelperEle_.GetRecoSFerr(aEle);
      wgtUp *= (nominal+err)/nominal;
      wgtDn *= (nominal-err)/nominal;
    } // acceptEles

    if (!ignore) {
      for (const auto& aEle : nonHeepEles) {
        double nominal = systHelperEle_.GetRecoSF(aEle);
        double err = systHelperEle_.GetRecoSFerr(aEle);
        wgtUp *= (nominal+err)/nominal;
        wgtDn *= (nominal-err)/nominal;
      }
    }

    return std::make_pair(wgtUp,wgtDn);
  };

  auto valFF = [this] (const double dr, const double pt) {
    // const double ff1 = ff_dr_->Eval(dr);
    // const double ff2 = ff_pt_->Eval(pt);

    // return (ff1+ff2)/2.;
    // return ff_dr_->GetBinContent( ff_dr_->FindFixBin(dr) );
    return ff_pt_->Eval(pt);
  };

  auto ciFF = [this] (const double dr, const double pt) {
    // const double ff1 = ff_dr_->Eval(dr);
    // const double ff2 = ff_pt_->Eval(pt);
    //
    // return std::abs(ff1-ff2)/2.;

    const double xval[1] = {pt};
    double ci[1];
    ffFitResult_->GetConfidenceIntervals(1,1,0,xval,ci,0.68,false);

    return ci[0];
  };

  auto ciFFdr03 = [this] (const double ptll) {
    const double xval[1] = {ptll};
    double ci[1];
    ffdr03FitResult_->GetConfidenceIntervals(1,1,0,xval,ci,0.68,false);

    return ci[0];
  };

  auto isMergedEle = [this,&eleHandle,&addGsfTrkMap] (const pat::ElectronRef& aEle) {
    bool isNotMerged = MergedLeptonHelperFct::isNotMerged(aEle,eleHandle,(*addGsfTrkMap)[aEle]);

    return aEle->electronID("mvaMergedElectron") && !isNotMerged;
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
          const double m4l = std::min(lvec4l.M(),2499.9);
          const double m4l_scaleUp = std::min((lvecE1_scaleUp+lvecE2_scaleUp+lvecE3_scaleUp+lvecE4_scaleUp).M(),2499.9);
          const double m4l_scaleDn = std::min((lvecE1_scaleDn+lvecE2_scaleDn+lvecE3_scaleDn+lvecE4_scaleDn).M(),2499.9);
          const double m4l_sigmaUp = std::min((lvecE1_sigmaUp+lvecE2_sigmaUp+lvecE3_sigmaUp+lvecE4_sigmaUp).M(),2499.9);
          const double m4l_sigmaDn = std::min((lvecE1_sigmaDn+lvecE2_sigmaDn+lvecE3_sigmaDn+lvecE4_sigmaDn).M(),2499.9);

          if ( m4l > 50. /*&& (m4l < 500. || isMC_) (unblinded)*/ && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
            histo1d_["cutflow_4E"]->Fill(9.5,aWeight);

            histo1d_["4P0F_CR_llll_invM"]->Fill(m4l, aWeight);
            histo1d_["4P0F_CR_llll_invM_scaleUp"]->Fill(m4l_scaleUp, aWeight);
            histo1d_["4P0F_CR_llll_invM_scaleDn"]->Fill(m4l_scaleDn, aWeight);
            histo1d_["4P0F_CR_llll_invM_sigmaUp"]->Fill(m4l_sigmaUp, aWeight);
            histo1d_["4P0F_CR_llll_invM_sigmaDn"]->Fill(m4l_sigmaDn, aWeight);
            histo1d_["4P0F_CR_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
            histo1d_["4P0F_CR_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
            histo1d_["4P0F_CR_llll_invM_elRecoUp"]->Fill(m4l, recoSFcl95(aWeight).first);
            histo1d_["4P0F_CR_llll_invM_elRecoDn"]->Fill(m4l, recoSFcl95(aWeight).second);
            histo1d_["4P0F_CR_llll_invM_elTrigUp"]->Fill(m4l, aWeight*systHelperEle_.GetTrigSFcl95UpDn(acceptEles.front()).first
                                                                     *systHelperEle_.GetTrigUnseededSFcl95UpDn(acceptEles.front()).first);
            histo1d_["4P0F_CR_llll_invM_elTrigDn"]->Fill(m4l, aWeight*systHelperEle_.GetTrigSFcl95UpDn(acceptEles.front()).second
                                                                     *systHelperEle_.GetTrigUnseededSFcl95UpDn(acceptEles.front()).second);
            histo1d_["4P0F_CR_llll_invM_PUrwgtUp"]->Fill(m4l, aWeight*purwgtUp/purwgtNo);
            histo1d_["4P0F_CR_llll_invM_PUrwgtDn"]->Fill(m4l, aWeight*purwgtDn/purwgtNo);
            histo1d_["4P0F_CR_llll_invM_prefireUp"]->Fill(m4l, aWeight*prefireUp/prefireNo);
            histo1d_["4P0F_CR_llll_invM_prefireDn"]->Fill(m4l, aWeight*prefireDn/prefireNo);
          }

          if ( m4l > 50. /*&& m4l < 500. (unblinded)*/ && lvecA1.M() > 1. && lvecA2.M() > 1. ) {
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
        for (unsigned idx = 0; idx < acceptEles.size(); idx++) {
          if ( idx!=0 && acceptEles.at(2)->et() < etThres1_ )
            continue;

          const auto& first = acceptEles.at(idx);
          const auto& second = acceptEles.at( (idx+1)%acceptEles.size() );
          const auto& probe = acceptEles.at( (idx+2)%acceptEles.size() );

          const auto lvecE1 = first->polarP4();
          const auto lvecE2 = second->polarP4();
          const auto lvecll = lvecE1 + lvecE2;
          const double mll = lvecll.M();
          const double dr2ll = reco::deltaR2(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi());

          bool hasMergedEle = isMergedEle(first) || isMergedEle(second) || isMergedEle(probe);

          if ( mll > 84.19 && mll < 98.19 && first->charge()*second->charge() < 0 && !hasMergedEle ) {
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

            const double dr2_1 = reco::deltaR2(probe->eta(),probe->phi(),first->eta(),first->phi());
            const double dr2_2 = reco::deltaR2(probe->eta(),probe->phi(),second->eta(),second->phi());
            const auto lvecClosest = dr2_1 < dr2_2 ? lvecE1 : lvecE2;
            const double mindr2 = std::min( dr2_1, dr2_2 );

            histo1d_["3P0F_dr"]->Fill( std::sqrt(mindr2), aWeight );

            if (iEvent.id().run() < 319077)
              histo1d_["3P0F_dr_preHEM"]->Fill( std::sqrt(mindr2), aWeight );
            else
              histo1d_["3P0F_dr_postHEM"]->Fill( std::sqrt(mindr2), aWeight );

            if ( std::abs( probe->superCluster()->eta() ) < 1.5 ) {
              histo1d_["3P0F_Et_EB"]->Fill(probe->et(), aWeight);
              histo1d_["3P0F_dr_EB"]->Fill( std::sqrt(mindr2), aWeight );
            } else {
              histo1d_["3P0F_Et_EE"]->Fill(probe->et(), aWeight);
              histo1d_["3P0F_dr_EE"]->Fill( std::sqrt(mindr2), aWeight );
            }

            if ( probe->superCluster()->eta() < -1.566 ) {
              if (iEvent.id().run() < 319077)
                histo1d_["3P0F_dr_EEneg_preHEM"]->Fill( std::sqrt(mindr2), aWeight );
              else
                histo1d_["3P0F_dr_EEneg_postHEM"]->Fill( std::sqrt(mindr2), aWeight );
            }

            numerPt_ = probe->et();
            numerDr_ = std::sqrt(mindr2);
            numerWgt_ = aWeight;
            numerInvM_ = (lvecClosest + probe->polarP4()).M();
            numerIsEB_ = static_cast<int>( std::abs( probe->superCluster()->eta() ) < 1.5 );
            numerTree_->Fill();
          } // Z tagging ( M_Z - 7 < M < M_Z + 7 )
        } // Z candidate combination
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

          const auto& fakePair = ( nonHeepEles.front()==pair1.first || nonHeepEles.front()==pair1.second ) ? pair1 : pair2;
          const auto& passPartner = ( nonHeepEles.front()==fakePair.first ) ? fakePair.second : fakePair.first;
          const double drFake = reco::deltaR(passPartner->eta(),passPartner->phi(),nonHeepEles.front()->eta(),nonHeepEles.front()->phi());

          const double ff = std::max( valFF(drFake,nonHeepEles.front()->et()), 0.);
          const double ci = ciFF(drFake,nonHeepEles.front()->et());

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
            histo1d_["3P1F_CR_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
            histo1d_["3P1F_CR_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
            histo1d_["3P1F_CR_llll_invM_xFF"]->Fill(m4l, aWeight*ff);
            histo1d_["3P1F_CR_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(ff+ci));
            histo1d_["3P1F_CR_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*(ff-std::min(ff,ci)));
            histo1d_["3P1F_CR_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*ff);
            histo1d_["3P1F_CR_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*ff);
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

            histo1d_["3P1F_CR_Et_P1"]->Fill(acceptEles.front()->et(), aWeight);
            histo1d_["3P1F_CR_Et_P2"]->Fill(acceptEles.at(1)->et(), aWeight);
            histo1d_["3P1F_CR_Et_P3"]->Fill(acceptEles.at(2)->et(), aWeight);
            histo1d_["3P1F_CR_Et_F1"]->Fill(nonHeepEles.front()->et(), aWeight);

            if ( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 ) {
              histo1d_["3P1F_CR_Et_EB"]->Fill(nonHeepEles.front()->et(), aWeight);
              histo1d_["3P1F_CR_dr_EB"]->Fill(drFake, aWeight);
            } else {
              histo1d_["3P1F_CR_Et_EE"]->Fill(nonHeepEles.front()->et(), aWeight);
              histo1d_["3P1F_CR_dr_EE"]->Fill(drFake, aWeight);
            }

            histo1d_["3P1F_CR_ll1ll2_dr_xFF"]->Fill(std::sqrt(dr2A12), aWeight*ff);

            histo1d_["3P1F_CR_ll1_invM_xFF"]->Fill(lvecA1.M(), aWeight*ff);
            histo1d_["3P1F_CR_ll1_dr_xFF"]->Fill(std::sqrt(dr2ll1), aWeight*ff);
            histo1d_["3P1F_CR_ll2_invM_xFF"]->Fill(lvecA2.M(), aWeight*ff);
            histo1d_["3P1F_CR_ll2_dr_xFF"]->Fill(std::sqrt(dr2ll2), aWeight*ff);

            pt_3P1F_ = nonHeepEles.front()->et();
            dr_3P1F_ = drFake;
            dPerpIn_3P1F_ = (*dPerpInHandle)[nonHeepEles.front()];
            dEtaInSeed_3P1F_ = nonHeepEles.front()->deltaEtaSeedClusterTrackAtVtx();
            mvaScore_3P1F_ = nonHeepEles.front()->userFloat("mvaMergedElectronValues");
            mvaCategory_3P1F_ = nonHeepEles.front()->userInt("mvaMergedElectronCategories");
            passMva_3P1F_ = nonHeepEles.front()->electronID("mvaMergedElectron");
            iso_3P1F_ = (*modifiedTrkIsoHandle)[nonHeepEles.front()];
            dxy_3P1F_ = nonHeepEles.front()->gsfTrack()->dxy(primaryVertex->position());
            ecalIso_3P1F_ = (*modifiedEcalIsoHandle)[nonHeepEles.front()];
            hcalIso_3P1F_ = nonHeepEles.front()->dr03HcalDepth1TowerSumEt();
            HoE_3P1F_ = nonHeepEles.front()->hadronicOverEm();
            rho_3P1F_ = *rhoHandle;
            enSC_3P1F_ = nonHeepEles.front()->superCluster()->energy();
            mih_3P1F_ = nonHeepEles.front()->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
            mll_3P1F_ = (nonHeepEles.front()->polarP4() + passPartner->polarP4()).M();
            m4l_3P1F_ = m4l;
            wgt_3P1F_ = aWeight;
            tree_3P1F_->Fill();
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

        bool hasMergedEle = isMergedEle(acceptEles.at(0)) || isMergedEle(acceptEles.at(1)) || isMergedEle(nonHeepEles.front());

        if ( mll > 84.19 && mll < 98.19 && acceptEles.front()->charge()*acceptEles.at(1)->charge() < 0 && !hasMergedEle ) {
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

          const double dr2_1 = reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.front()->eta(),acceptEles.front()->phi());
          const double dr2_2 = reco::deltaR2(nonHeepEles.front()->eta(),nonHeepEles.front()->phi(),acceptEles.at(1)->eta(),acceptEles.at(1)->phi());
          const auto lvecClosest = dr2_1 < dr2_2 ? lvecE1 : lvecE2;
          const double mindr2 = std::min( dr2_1, dr2_2 );

          histo1d_["2P1F_dr"]->Fill(std::sqrt(mindr2), aWeight);

          if (iEvent.id().run() < 319077)
            histo1d_["2P1F_dr_preHEM"]->Fill( std::sqrt(mindr2), aWeight );
          else
            histo1d_["2P1F_dr_postHEM"]->Fill( std::sqrt(mindr2), aWeight );

          if ( nonHeepEles.front()->superCluster()->eta() < -1.566 ) {
            if (iEvent.id().run() < 319077)
              histo1d_["2P1F_dr_EEneg_preHEM"]->Fill( std::sqrt(mindr2), aWeight );
            else
              histo1d_["2P1F_dr_EEneg_postHEM"]->Fill( std::sqrt(mindr2), aWeight );
          }

          if ( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 ) {
            histo1d_["2P1F_Et_EB"]->Fill(nonHeepEles.front()->et(), aWeight);
            histo1d_["2P1F_dr_EB"]->Fill(std::sqrt(mindr2), aWeight);
          } else {
            histo1d_["2P1F_Et_EE"]->Fill(nonHeepEles.front()->et(), aWeight);
            histo1d_["2P1F_dr_EE"]->Fill(std::sqrt(mindr2), aWeight);
          }

          denomPt_ = nonHeepEles.front()->et();
          denomDr_ = std::sqrt(mindr2);
          denomWgt_ = aWeight;
          denomDPerpIn_ = (*dPerpInHandle)[nonHeepEles.front()];
          denomDEtaInSeed_ = nonHeepEles.front()->deltaEtaSeedClusterTrackAtVtx();
          denomMvaScore_ = nonHeepEles.front()->userFloat("mvaMergedElectronValues");
          denomMvaCategory_ = nonHeepEles.front()->userInt("mvaMergedElectronCategories");
          denomPassMva_ = nonHeepEles.front()->electronID("mvaMergedElectron");
          denomIso_ = (*modifiedTrkIsoHandle)[nonHeepEles.front()];
          denomDxy_ = nonHeepEles.front()->gsfTrack()->dxy(primaryVertex->position());
          denomEcalIso_ = (*modifiedEcalIsoHandle)[nonHeepEles.front()];
          denomHcalIso_ = nonHeepEles.front()->dr03HcalDepth1TowerSumEt();
          denomHoE_ = nonHeepEles.front()->hadronicOverEm();
          denomRho_ = *rhoHandle;
          denomE_ = nonHeepEles.front()->superCluster()->energy();
          denomMIH_ = nonHeepEles.front()->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
          denomInvM_ = (lvecClosest + nonHeepEles.front()->polarP4()).M();
          denomIsEB_ = static_cast<int>( std::abs( nonHeepEles.front()->superCluster()->eta() ) < 1.5 );
          denomTree_->Fill();
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
          const auto firstFakePair = ( pair1.first==nonHeepEles.front() || pair1.second==nonHeepEles.front() )
                                     ? pair1 : pair2;
          const auto secondFakePair = ( pair1.first==nonHeepEles.at(1) || pair1.second==nonHeepEles.at(1) )
                                     ? pair1 : pair2;

          if (true) { // !pairedFakes
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

            const double drFake1 = reco::deltaR( firstFakePair.first->eta(), firstFakePair.first->phi(),
                                                 firstFakePair.second->eta(), firstFakePair.second->phi() );
            const double drFake2 = reco::deltaR( secondFakePair.first->eta(), secondFakePair.first->phi(),
                                                 secondFakePair.second->eta(), secondFakePair.second->phi() );
            const double ff1 = std::max( valFF(drFake1,nonHeepEles.front()->et()), 0.);
            const double ff2 = std::max( valFF(drFake2,nonHeepEles.at(1)->et()), 0.);
            const double ci1 = ciFF(drFake1,nonHeepEles.front()->et());
            const double ci2 = ciFF(drFake2,nonHeepEles.at(1)->et());

            if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. && drFake1 > 0.3 && drFake2 > 0.3 ) {
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

              pt1_2P2F_ = nonHeepEles.front()->et();
              pt2_2P2F_ = nonHeepEles.at(1)->et();
              dr1_2P2F_ = drFake1;
              dr2_2P2F_ = drFake2;
              dPerpIn1_2P2F_ = (*dPerpInHandle)[nonHeepEles.front()];
              dPerpIn2_2P2F_ = (*dPerpInHandle)[nonHeepEles.at(1)];
              dEtaInSeed1_2P2F_ = nonHeepEles.front()->deltaEtaSeedClusterTrackAtVtx();
              dEtaInSeed2_2P2F_ = nonHeepEles.at(1)->deltaEtaSeedClusterTrackAtVtx();
              mvaScore1_2P2F_ = nonHeepEles.front()->userFloat("mvaMergedElectronValues");
              mvaScore2_2P2F_ = nonHeepEles.at(1)->userFloat("mvaMergedElectronValues");
              mvaCategory1_2P2F_ = nonHeepEles.front()->userInt("mvaMergedElectronCategories");
              mvaCategory2_2P2F_ = nonHeepEles.at(1)->userInt("mvaMergedElectronCategories");
              passMva1_2P2F_ = nonHeepEles.front()->electronID("mvaMergedElectron");
              passMva2_2P2F_ = nonHeepEles.at(1)->electronID("mvaMergedElectron");
              iso1_2P2F_ = (*modifiedTrkIsoHandle)[nonHeepEles.front()];
              iso2_2P2F_ = (*modifiedTrkIsoHandle)[nonHeepEles.at(1)];
              dxy1_2P2F_ = nonHeepEles.front()->gsfTrack()->dxy(primaryVertex->position());
              dxy2_2P2F_ = nonHeepEles.at(1)->gsfTrack()->dxy(primaryVertex->position());
              ecalIso1_2P2F_ = (*modifiedEcalIsoHandle)[nonHeepEles.front()];
              ecalIso2_2P2F_ = (*modifiedEcalIsoHandle)[nonHeepEles.at(1)];
              hcalIso1_2P2F_ = nonHeepEles.front()->dr03HcalDepth1TowerSumEt();
              hcalIso2_2P2F_ = nonHeepEles.at(1)->dr03HcalDepth1TowerSumEt();
              HoE1_2P2F_ = nonHeepEles.front()->hadronicOverEm();
              HoE2_2P2F_ = nonHeepEles.at(1)->hadronicOverEm();
              enSC1_2P2F_ = nonHeepEles.front()->superCluster()->energy();
              enSC2_2P2F_ = nonHeepEles.at(1)->superCluster()->energy();
              rho_2P2F_ = *rhoHandle;
              mih1_2P2F_ = nonHeepEles.front()->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
              mih2_2P2F_ = nonHeepEles.at(1)->gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
              mll1_2P2F_ = (firstFakePair.first->polarP4() + firstFakePair.second->polarP4()).M();
              mll2_2P2F_ = (secondFakePair.first->polarP4() + secondFakePair.second->polarP4()).M();
              m4l_2P2F_ = m4l;
              wgt_2P2F_ = aWeight;
              pairedFakes_2P2F_ = static_cast<int>(pairedFakes);
              tree_2P2F_->Fill();
            } // CR

            if ( m4l > 50. && lvecA1.M() > 1. && lvecA2.M() > 1. && pairedFakes && drFake1 < 0.3 ) {
              const auto lvecllFake = firstFakePair.first->polarP4() + firstFakePair.second->polarP4();
              const double ptll = lvecllFake.Pt();
              const double valFFdr03 = ffdr03_->Eval(ptll);
              const double cidr03 = ciFFdr03(ptll);
              histo1d_["2P2F_CRdr03_llll_invM"]->Fill(m4l, aWeight);
              histo1d_["2P2F_CRdr03_llll_invM_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first);
              histo1d_["2P2F_CRdr03_llll_invM_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second);
              histo1d_["2P2F_CRdr03_llll_invM_xFF"]->Fill(m4l, aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_llll_invM_xFF_ffUp"]->Fill(m4l, aWeight*(valFFdr03+cidr03));
              histo1d_["2P2F_CRdr03_llll_invM_xFF_ffDn"]->Fill(m4l, aWeight*std::max(0.,valFFdr03-cidr03));
              histo1d_["2P2F_CRdr03_llll_invM_xFF_idUp"]->Fill(m4l, modHeepSFcl95(aWeight).first*valFFdr03);
              histo1d_["2P2F_CRdr03_llll_invM_xFF_idDn"]->Fill(m4l, modHeepSFcl95(aWeight).second*valFFdr03);

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
              histo1d_["2P2F_CRdr03_ll1_pt_xFF"]->Fill(lvecA1.pt(), aWeight*valFFdr03);
              histo1d_["2P2F_CRdr03_ll2_pt_xFF"]->Fill(lvecA2.pt(), aWeight*valFFdr03);

              const auto e1 = firstFakePair.first;
              const auto e2 = firstFakePair.second;
              var_weight_ = aWeight;
              var_e1Pt_ = e1->pt();
              var_e1scEn_ = e1->superCluster()->energy();
              var_e1scEta_ = e1->superCluster()->eta();
              var_e1Eta_ = e1->eta();
              var_e1Phi_ = e1->phi();
              var_e1DPerpIn_ = (*dPerpInHandle)[e1];
              var_e1DEtaInSeed_ = e1->deltaEtaSeedClusterTrackAtVtx();
              var_e1MvaScore_ = e1->userFloat("mvaMergedElectronValues");
              var_e1MvaCategory_ = e1->userInt("mvaMergedElectronCategories");
              var_e1PassMva_ = e1->electronID("mvaMergedElectron");
              var_e1TrkIso_ = (*modifiedTrkIsoHandle)[e1];
              var_e1EcalIso_ = (*modifiedEcalIsoHandle)[e1];
              var_e1HcalIso_ = e1->dr03HcalDepth1TowerSumEt();
              var_e1HoE_ = e1->hadronicOverEm();
              var_e2Pt_ = e2->pt();
              var_e2scEn_ = e2->superCluster()->energy();
              var_e2scEta_ = e2->superCluster()->eta();
              var_e2Eta_ = e2->eta();
              var_e2Phi_ = e2->phi();
              var_e2DPerpIn_ = (*dPerpInHandle)[e2];
              var_e2DEtaInSeed_ = e2->deltaEtaSeedClusterTrackAtVtx();
              var_e2MvaScore_ = e2->userFloat("mvaMergedElectronValues");
              var_e2MvaCategory_ = e2->userInt("mvaMergedElectronCategories");
              var_e2PassMva_ = e2->electronID("mvaMergedElectron");
              var_e2TrkIso_ = (*modifiedTrkIsoHandle)[e2];
              var_e2EcalIso_ = (*modifiedEcalIsoHandle)[e2];
              var_e2HcalIso_ = e2->dr03HcalDepth1TowerSumEt();
              var_e2HoE_ = e2->hadronicOverEm();
              var_e1e2Rho_ = *rhoHandle;
              var_e1e2InvM_ = (e1->polarP4()+e2->polarP4()).M();
              var_e1e2DR_ = reco::deltaR(e1->eta(),e1->phi(),e2->eta(),e2->phi());
              var_e1e2scDR_ = reco::deltaR(e1->superCluster()->eta(),e1->superCluster()->phi(),
                                           e2->superCluster()->eta(),e2->superCluster()->phi());
              tree_2P2F_dr03_->Fill();
            } // nonisolated CR
          } // !pairedFakes
        } // paired
      } // 2P2F

      break;
      // case acceptEles.size()==2
  } // switch acceptEles.size()
}

DEFINE_FWK_MODULE(ResolvedEleCRanalyzer);
