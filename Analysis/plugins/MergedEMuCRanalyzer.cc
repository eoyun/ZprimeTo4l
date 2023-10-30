#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "TF1.h"
#include "TFitResult.h"

class MergedEMuCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEMuCRanalyzer(const edm::ParameterSet&);
  ~MergedEMuCRanalyzer() override {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;

  const edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;

  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dEtaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5EnergyToken_;

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const edm::FileInPath recoSFpath_;
  const edm::FileInPath FFpath_;
  const edm::FileInPath MMFFpath_;

  const std::string trigHistName_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath purwgtPath_;

  const double ptThres_;
  const double ptMuThres_;
  const double drThres_;
  const double drThresCR_;
  const double ratioThresLo_;
  const double ratioThresHi_;
  const double mumass_ = 0.1056583745;

  MuonCorrectionHelper mucorrHelper_;

  std::unique_ptr<TFile> recoSFfile_;
  TH2D* recoSF_;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::unique_ptr<TFile> FFfile_;
  TF1* ssboth_;
  TF1* osboth_;
  TFitResultPtr ssFit_;
  TFitResultPtr osFit_;
  TH1D* osPar_;

  std::unique_ptr<TFile> MMFFfile_;
  TF1* ffFunc_;
  TFitResultPtr fitResult_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;
};

MergedEMuCRanalyzer::MergedEMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
union5x5dEtaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dEtaIn"))),
union5x5EnergyToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5Energy"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
recoSFpath_(iConfig.getParameter<edm::FileInPath>("recoSFpath")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
MMFFpath_(iConfig.getParameter<edm::FileInPath>("MMFFpath")),
trigHistName_(iConfig.getParameter<std::string>("trigHistName")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptMuThres_(iConfig.getParameter<double>("ptMuThres")),
drThres_(iConfig.getParameter<double>("drThres")),
drThresCR_(iConfig.getParameter<double>("drThresCR")),
ratioThresLo_(iConfig.getParameter<double>("ratioThresLo")),
ratioThresHi_(iConfig.getParameter<double>("ratioThresHi")),
mucorrHelper_(rochesterPath_,triggerSFpath_,trigHistName_) {
  usesResource("TFileService");
}

void MergedEMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  recoSFfile_ = std::make_unique<TFile>(recoSFpath_.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  ssboth_ = static_cast<TF1*>(FFfile_->FindObjectAny("ssboth"));
  osboth_ = static_cast<TF1*>(FFfile_->FindObjectAny("osboth"));
  ssFit_ = (static_cast<TH1D*>(FFfile_->Get("2E_Et_SSCR_EB_mixedME")))->Fit(ssboth_,"RS");
  osFit_ = (static_cast<TH1D*>(FFfile_->Get("2E_Et_OSCR_EB_mixedME")))->Fit(osboth_,"RS");
  osPar_ = static_cast<TH1D*>(FFfile_->Get("2E_eta_OSCR_EB_mixedME"));

  MMFFfile_ = std::make_unique<TFile>(MMFFpath_.fullPath().c_str(),"READ");
  ffFunc_ = static_cast<TF1*>(MMFFfile_->FindObjectAny("MMFF"));
  fitResult_ = (static_cast<TH1D*>(MMFFfile_->Get("1M_MMFF_numer_rebin")))->Fit(ffFunc_,"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["totWeightedSum_2E2M"] = fs->make<TH1D>("totWeightedSum_2E2M","totWeightedSum_2E2M",1,0.,1.);

  histo1d_["cutflow_2M"] = fs->make<TH1D>("cutflow_2M","cutflow",10,0.,10.);
  histo1d_["cutflow_1M"] = fs->make<TH1D>("cutflow_1M","cutflow",15,0.,15.);

  histo1d_["2M_CRME_Et"] = fs->make<TH1D>("2M_CRME_Et",";E_{T};",40,0.,200.);
  histo1d_["2M_CRME_eta"] = fs->make<TH1D>("2M_CRME_eta",";#eta;",50,-2.5,2.5);
  histo1d_["2M_CRME_phi"] = fs->make<TH1D>("2M_CRME_phi",";#phi;",64,-3.2,3.2);

  histo1d_["2M_CRME_noGsf_Et"] = fs->make<TH1D>("2M_CRME_noGsf_Et",";E_{T};",40,0.,200.);
  histo1d_["2M_CRME_noGsf_eta"] = fs->make<TH1D>("2M_CRME_noGsf_eta",";#eta;",50,-2.5,2.5);
  histo1d_["2M_CRME_noGsf_phi"] = fs->make<TH1D>("2M_CRME_noGsf_phi",";#phi;",64,-3.2,3.2);

  histo1d_["2M_CRME_all_eta"] = fs->make<TH1D>("2M_CRME_all_eta",";#eta;",50,-2.5,2.5);

  histo1d_["2M_CRME_M1_pt"] = fs->make<TH1D>("2M_CRME_M1_pt",";p_{T};",40,0.,200.);
  histo1d_["2M_CRME_M1_eta"] = fs->make<TH1D>("2M_CRME_M1_eta",";#eta;",50,-2.5,2.5);
  histo1d_["2M_CRME_M1_phi"] = fs->make<TH1D>("2M_CRME_M1_phi",";#phi;",64,-3.2,3.2);

  histo1d_["2M_CRME_M2_pt"] = fs->make<TH1D>("2M_CRME_M2_pt",";p_{T};",40,0.,200.);
  histo1d_["2M_CRME_M2_eta"] = fs->make<TH1D>("2M_CRME_M2_eta",";#eta;",50,-2.5,2.5);
  histo1d_["2M_CRME_M2_phi"] = fs->make<TH1D>("2M_CRME_M2_phi",";#phi;",64,-3.2,3.2);

  histo1d_["2M_CRME_lll_invM"] = fs->make<TH1D>("2M_CRME_lll_invM",";M(3l);",1000,0.,2500.);
  histo1d_["2M_CRME_lll_pt"] = fs->make<TH1D>("2M_CRME_lll_pt",";p_{T}(3l);",100,0.,500.);
  histo1d_["2M_CRME_lll_dr_l1ME"] = fs->make<TH1D>("2M_CRME_lll_dr_l1ME","dR(l1ME)",64,0.,6.4);
  histo1d_["2M_CRME_lll_dr_l2ME"] = fs->make<TH1D>("2M_CRME_lll_dr_l2ME","dR(l2ME)",64,0.,6.4);
  histo1d_["2M_CRME_lll_dr_l1l2"] = fs->make<TH1D>("2M_CRME_lll_dr_l1l2","dR(l1l2)",64,0.,6.4);

  histo1d_["2M_antiME_Et"] = fs->make<TH1D>("2M_antiME_Et",";E_{T};",200,0.,500.);
  histo1d_["2M_antiME_eta"] = fs->make<TH1D>("2M_antiME_eta",";#eta;",200,-2.5,2.5);
  histo1d_["2M_antiME_phi"] = fs->make<TH1D>("2M_antiME_phi",";#phi;",128,-3.2,3.2);

  histo1d_["2M_antiME_noGsf_Et"] = fs->make<TH1D>("2M_antiME_noGsf_Et",";E_{T};",200,0.,500.);
  histo1d_["2M_antiME_noGsf_eta"] = fs->make<TH1D>("2M_antiME_noGsf_eta",";#eta;",200,-2.5,2.5);
  histo1d_["2M_antiME_noGsf_phi"] = fs->make<TH1D>("2M_antiME_noGsf_phi",";#phi;",128,-3.2,3.2);

  histo1d_["2M_antiME_M1_pt"] = fs->make<TH1D>("2M_antiME_M1_Et",";p_{T};",200,0.,500.);
  histo1d_["2M_antiME_M1_eta"] = fs->make<TH1D>("2M_antiME_M1_eta",";#eta;",200,-2.5,2.5);
  histo1d_["2M_antiME_M1_phi"] = fs->make<TH1D>("2M_antiME_M1_phi",";#phi;",128,-3.2,3.2);

  histo1d_["2M_antiME_M2_pt"] = fs->make<TH1D>("2M_antiME_M2_Et",";p_{T};",200,0.,500.);
  histo1d_["2M_antiME_M2_eta"] = fs->make<TH1D>("2M_antiME_M2_eta",";#eta;",200,-2.5,2.5);
  histo1d_["2M_antiME_M2_phi"] = fs->make<TH1D>("2M_antiME_M2_phi",";#phi;",128,-3.2,3.2);

  histo1d_["2M_antiME_lll_invM_CR"] = fs->make<TH1D>("2M_antiME_lll_invM_CR",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_invM_CR_xSSFF"] = fs->make<TH1D>("2M_antiME_lll_invM_CR_xSSFF",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_invM_CR_xOSFF"] = fs->make<TH1D>("2M_antiME_lll_invM_CR_xOSFF",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_invM_CR_xSSFF_up"] = fs->make<TH1D>("2M_antiME_lll_invM_CR_xSSFF_up",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_invM_CR_xOSFF_up"] = fs->make<TH1D>("2M_antiME_lll_invM_CR_xOSFF_up",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_invM_CR_xSSFF_dn"] = fs->make<TH1D>("2M_antiME_lll_invM_CR_xSSFF_dn",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_invM_CR_xOSFF_dn"] = fs->make<TH1D>("2M_antiME_lll_invM_CR_xOSFF_dn",";M(3l);",1000,0.,2500.);
  histo1d_["2M_antiME_lll_pt"] = fs->make<TH1D>("2M_antiME_lll_pt",";p_{T}(3l);",200,0.,1000.);
  histo1d_["2M_antiME_lll_dr_l1ME"] = fs->make<TH1D>("2M_antiME_lll_dr_l1ME","dR(l1ME)",128,0.,6.4);
  histo1d_["2M_antiME_lll_dr_l2ME"] = fs->make<TH1D>("2M_antiME_lll_dr_l2ME","dR(l2ME)",128,0.,6.4);
  histo1d_["2M_antiME_lll_dr_l1l2"] = fs->make<TH1D>("2M_antiME_lll_dr_l1l2","dR(l1l2)",128,0.,6.4);
  histo1d_["2M_antiME_eta_xSSFF"] = fs->make<TH1D>("2M_antiME_eta_xSSFF","CR anti ME #eta x SS Fake factor",100,-2.5,2.5);
  histo1d_["2M_antiME_eta_xOSFF"] = fs->make<TH1D>("2M_antiME_eta_xOSFF","CR anti ME #eta x OS Fake factor",100,-2.5,2.5);
  histo1d_["2M_antiME_eta_xSSFF_up"] = fs->make<TH1D>("2M_antiME_eta_xSSFF_up",";#eta;",100,-2.5,2.5);
  histo1d_["2M_antiME_eta_xOSFF_up"] = fs->make<TH1D>("2M_antiME_eta_xOSFF_up",";#eta;",100,-2.5,2.5);
  histo1d_["2M_antiME_eta_xSSFF_dn"] = fs->make<TH1D>("2M_antiME_eta_xSSFF_dn",";#eta;",100,-2.5,2.5);
  histo1d_["2M_antiME_eta_xOSFF_dn"] = fs->make<TH1D>("2M_antiME_eta_xOSFF_dn",";#eta;",100,-2.5,2.5);

  histo1d_["2M_Et_CR_EB_antiME"] = fs->make<TH1D>("2M_Et_CR_EB_antiME",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_antiME_xSSFF"] = fs->make<TH1D>("2M_Et_CR_EB_antiME_xSSFF",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_antiME_xSSFF_up"] = fs->make<TH1D>("2M_Et_CR_EB_antiME_xSSFF_up",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_antiME_xSSFF_dn"] = fs->make<TH1D>("2M_Et_CR_EB_antiME_xSSFF_dn",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_antiME_xOSFF"] = fs->make<TH1D>("2M_Et_CR_EB_antiME_xOSFF",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_antiME_xOSFF_up"] = fs->make<TH1D>("2M_Et_CR_EB_antiME_xOSFF_up",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_antiME_xOSFF_dn"] = fs->make<TH1D>("2M_Et_CR_EB_antiME_xOSFF_dn",";E_{T};",200,0.,1000.);
  histo1d_["2M_Et_CR_EB_CRME"] = fs->make<TH1D>("2M_Et_CR_EB_CRME",";E_{T};",200,0.,1000.);

  histo1d_["1M_CRME_MM_pt"] = fs->make<TH1D>("1M_CRME_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_MM_eta"] = fs->make<TH1D>("1M_CRME_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["1M_CRME_MM_phi"] = fs->make<TH1D>("1M_CRME_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["1M_CRME_Et"] = fs->make<TH1D>("1M_CRME_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_eta"] = fs->make<TH1D>("1M_CRME_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_phi"] = fs->make<TH1D>("1M_CRME_phi",";#phi;",64,-3.2,3.2);

  histo1d_["1M_CRME_noGsf_Et"] = fs->make<TH1D>("1M_CRME_noGsf_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_noGsf_eta"] = fs->make<TH1D>("1M_CRME_noGsf_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_noGsf_phi"] = fs->make<TH1D>("1M_CRME_noGsf_phi",";#phi;",64,-3.2,3.2);

  histo1d_["1M_CRME_all_Et"] = fs->make<TH1D>("1M_CRME_all_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_all_eta"] = fs->make<TH1D>("1M_CRME_all_eta",";#eta;",50,-2.5,2.5);

  histo1d_["1M_CRME_MET_pt"] = fs->make<TH1D>("1M_CRME_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_MET_phi"] = fs->make<TH1D>("1M_CRME_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_MET_dphi"] = fs->make<TH1D>("1M_CRME_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_ratioPt"] = fs->make<TH1D>("1M_CRME_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);
  histo1d_["1M_CRME_mt"] = fs->make<TH1D>("1M_CRME_mt","m_{T};m_{T};",500,0.,2500.);

  histo1d_["1M_CRME_antiRpt_MM_pt"] = fs->make<TH1D>("1M_CRME_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("1M_CRME_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiRpt_MM_eta"] = fs->make<TH1D>("1M_CRME_antiRpt_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["1M_CRME_antiRpt_all_Et"] = fs->make<TH1D>("1M_CRME_antiRpt_all_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiRpt_all_Et_xFF"] = fs->make<TH1D>("1M_CRME_antiRpt_all_Et_xFF",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiRpt_all_eta"] = fs->make<TH1D>("1M_CRME_antiRpt_all_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_antiRpt_MET_pt"] = fs->make<TH1D>("1M_CRME_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("1M_CRME_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiRpt_MET_phi"] = fs->make<TH1D>("1M_CRME_antiRpt_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_antiRpt_MET_dphi"] = fs->make<TH1D>("1M_CRME_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_antiRpt_ratioPt"] = fs->make<TH1D>("1M_CRME_antiRpt_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);
  histo1d_["1M_CRME_antiRpt_mt"] = fs->make<TH1D>("1M_CRME_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_antiRpt_mt_xFF"] = fs->make<TH1D>("1M_CRME_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_antiRpt_mt_xFF_up"] = fs->make<TH1D>("1M_CRME_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("1M_CRME_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);

  histo1d_["1M_CRME_CRdphi_MM_pt"] = fs->make<TH1D>("1M_CRME_CRdphi_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_CRdphi_MM_eta"] = fs->make<TH1D>("1M_CRME_CRdphi_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["1M_CRME_CRdphi_MM_phi"] = fs->make<TH1D>("1M_CRME_CRdphi_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["1M_CRME_CRdphi_Et"] = fs->make<TH1D>("1M_CRME_CRdphi_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_CRdphi_eta"] = fs->make<TH1D>("1M_CRME_CRdphi_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_CRdphi_phi"] = fs->make<TH1D>("1M_CRME_CRdphi_phi",";#phi;",64,-3.2,3.2);

  histo1d_["1M_CRME_CRdphi_noGsf_Et"] = fs->make<TH1D>("1M_CRME_CRdphi_noGsf_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_CRdphi_noGsf_eta"] = fs->make<TH1D>("1M_CRME_CRdphi_noGsf_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_CRdphi_noGsf_phi"] = fs->make<TH1D>("1M_CRME_CRdphi_noGsf_phi",";#phi;",64,-3.2,3.2);

  histo1d_["1M_CRME_CRdphi_all_Et"] = fs->make<TH1D>("1M_CRME_CRdphi_all_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_CRdphi_all_eta"] = fs->make<TH1D>("1M_CRME_CRdphi_all_eta",";#eta;",50,-2.5,2.5);

  histo1d_["1M_CRME_CRdphi_MET_pt"] = fs->make<TH1D>("1M_CRME_CRdphi_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_CRdphi_MET_phi"] = fs->make<TH1D>("1M_CRME_CRdphi_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_CRdphi_MET_dphi"] = fs->make<TH1D>("1M_CRME_CRdphi_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_CRdphi_ratioPt"] = fs->make<TH1D>("1M_CRME_CRdphi_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);
  histo1d_["1M_CRME_CRdphi_mt"] = fs->make<TH1D>("1M_CRME_CRdphi_mt","m_{T};m_{T};",500,0.,2500.);

  histo1d_["1M_CRME_CRdphi_antiRpt_MM_pt"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_MM_eta"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["1M_CRME_CRdphi_antiRpt_all_Et"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_all_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_CRdphi_antiRpt_all_Et_xFF"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_all_Et_xFF",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_CRdphi_antiRpt_all_eta"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_all_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_CRdphi_antiRpt_MET_pt"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_MET_phi"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_CRdphi_antiRpt_MET_dphi"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_CRdphi_antiRpt_ratioPt"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);
  histo1d_["1M_CRME_CRdphi_antiRpt_mt"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_mt_xFF"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_mt_xFF_up"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_CRdphi_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("1M_CRME_CRdphi_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);

  histo1d_["1M_CRME_antiDphi_antiRpt_MM_pt"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_MM_eta"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["1M_CRME_antiDphi_antiRpt_all_Et"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_all_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiDphi_antiRpt_all_Et_xFF"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_all_Et_xFF",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiDphi_antiRpt_all_eta"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_all_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_antiDphi_antiRpt_MET_pt"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_MET_phi"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_antiDphi_antiRpt_MET_dphi"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_antiDphi_antiRpt_ratioPt"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);
  histo1d_["1M_CRME_antiDphi_antiRpt_mt"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_mt_xFF"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_mt_xFF_up"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["1M_CRME_antiDphi_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("1M_CRME_antiDphi_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);

  histo1d_["1M_CRME_antiDphi_MM_pt"] = fs->make<TH1D>("1M_CRME_antiDphi_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiDphi_MM_eta"] = fs->make<TH1D>("1M_CRME_antiDphi_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["1M_CRME_antiDphi_MM_phi"] = fs->make<TH1D>("1M_CRME_antiDphi_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["1M_CRME_antiDphi_Et"] = fs->make<TH1D>("1M_CRME_antiDphi_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiDphi_eta"] = fs->make<TH1D>("1M_CRME_antiDphi_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_antiDphi_phi"] = fs->make<TH1D>("1M_CRME_antiDphi_phi",";#phi;",64,-3.2,3.2);

  histo1d_["1M_CRME_antiDphi_noGsf_Et"] = fs->make<TH1D>("1M_CRME_antiDphi_noGsf_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiDphi_noGsf_eta"] = fs->make<TH1D>("1M_CRME_antiDphi_noGsf_eta",";#eta;",50,-2.5,2.5);
  histo1d_["1M_CRME_antiDphi_noGsf_phi"] = fs->make<TH1D>("1M_CRME_antiDphi_noGsf_phi",";#phi;",64,-3.2,3.2);

  histo1d_["1M_CRME_antiDphi_all_Et"] = fs->make<TH1D>("1M_CRME_antiDphi_all_Et",";E_{T};",40,0.,200.);
  histo1d_["1M_CRME_antiDphi_all_eta"] = fs->make<TH1D>("1M_CRME_antiDphi_all_eta",";#eta;",50,-2.5,2.5);

  histo1d_["1M_CRME_antiDphi_MET_pt"] = fs->make<TH1D>("1M_CRME_antiDphi_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["1M_CRME_antiDphi_MET_phi"] = fs->make<TH1D>("1M_CRME_antiDphi_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_antiDphi_MET_dphi"] = fs->make<TH1D>("1M_CRME_antiDphi_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_CRME_antiDphi_ratioPt"] = fs->make<TH1D>("1M_CRME_antiDphi_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);
  histo1d_["1M_CRME_antiDphi_mt"] = fs->make<TH1D>("1M_CRME_antiDphi_mt","m_{T};m_{T};",500,0.,2500.);

  histo1d_["1M_FF_CRME_antiDphi_MMMET_pt"] = fs->make<TH1D>("1M_FF_CRME_antiDphi_MMMET_pt",";p_{T}+MET;",200,0.,1000.);
  histo1d_["1M_FF_CRME_MMMET_pt"] = fs->make<TH1D>("1M_FF_CRME_MMMET_pt",";p_{T}+MET;",200,0.,1000.);
  histo1d_["1M_FF_CRME_antiDphi_MET_pt"] = fs->make<TH1D>("1M_FF_CRME_antiDphi_MET_pt",";MET;",200,0.,500.);
  histo1d_["1M_FF_CRME_MET_pt"] = fs->make<TH1D>("1M_FF_CRME_MET_pt",";MET;",200,0.,500.);

  histo1d_["1M_ABCD_MET_dphi"] = fs->make<TH1D>("1M_ABCD_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["1M_ABCD_MET_ratioPt"] = fs->make<TH1D>("1M_ABCD_MET_ratioPt",";#Sigma p_{T}^{l}/MET;",200,0.,5.);

  histo2d_["1M_mt_dphi"] = fs->make<TH2D>("1M_mt_dphi",";m_{T};#Delta#phi",100,0.,500.,128,-3.2,3.2);
  histo2d_["1M_mt_ratioPt"] = fs->make<TH2D>("1M_mt_ratioPt",";m_{T};#Sigma p_{T}^{l}/MET",100,0.,500.,100,0.,5.);
  histo2d_["1M_dphi_ratioPt"] = fs->make<TH2D>("1M_dphi_ratioPt",";#Delta#phi;#Sigma p_{T}^{l}/MET",128,-3.2,3.2,100,0.,5.);
}

void MergedEMuCRanalyzer::endJob() {
  recoSFfile_->Close();
  FFfile_->Close();
  MMFFfile_->Close();
  purwgtFile_->Close();
}

void MergedEMuCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::View<pat::MET>> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

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

    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(genptcToken_, genptcHandle);

    std::vector<reco::GenParticleRef> promptMuons;
    std::vector<reco::GenParticleRef> promptElectrons;

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( ( std::abs(genptc->pdgId())==13 ) && genptc->fromHardProcessFinalState() )
        promptMuons.push_back(genptc.castTo<reco::GenParticleRef>());

      if ( ( std::abs(genptc->pdgId())==11 ) && genptc->fromHardProcessFinalState() )
        promptElectrons.push_back(genptc.castTo<reco::GenParticleRef>());
    }

    if ( promptMuons.size()==2 && promptElectrons.size()==2 )
      histo1d_["totWeightedSum_2E2M"]->Fill(0.5,aWeight);
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  const edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    const std::string& trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
        if (trigResultHandle.product()->accept(iTrig))
          isFired = true;
      }
    } // wanted triggers
  } // fired triggers

  std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjs;

  for (unsigned iTrig=0; iTrig < trigObjHandle->size(); iTrig++) {
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

  histo1d_["cutflow_2M"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_1M"]->Fill( 0.5, aWeight );

  if (!isFired)
    return;

  histo1d_["cutflow_2M"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_1M"]->Fill( 1.5, aWeight );

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

  histo1d_["cutflow_2M"]->Fill( 2.5, aWeight );
  histo1d_["cutflow_1M"]->Fill( 2.5, aWeight );

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons; // but not highPtMuon

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->pt() < ptMuThres_ || std::abs(aMuon->eta()) > 2.4 )
      continue;

    const auto casted = aMuon.castTo<pat::MuonRef>();

    if ( muon::isHighPtMuon(*aMuon,primaryVertex) )
      highPtMuons.push_back(casted);
    else if ( muon::isTrackerHighPtMuon(*aMuon,primaryVertex) )
      highPtTrackerMuons.push_back(casted);
    else {}
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

    if ( firstIso/firstMuon->pt() < 0.05 )
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

    if ( firstIso/firstMuon->pt() < 0.05 )
      isolatedHighPtTrackerMuons.push_back(firstMuon);
  }

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  if ( isolatedHighPtMuons.empty() )
    return;

  if ( isolatedHighPtMuons.front()->pt() < 52. )
    return;

  histo1d_["cutflow_2M"]->Fill( 3.5, aWeight );
  histo1d_["cutflow_1M"]->Fill( 3.5, aWeight );

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

  histo1d_["cutflow_2M"]->Fill( 4.5, aWeight );
  histo1d_["cutflow_1M"]->Fill( 4.5, aWeight );

  // concatenate(-ish) highPtMuons & highPtTrackerMuons - order is important!
  std::vector<pat::MuonRef> allHighPtMuons(isolatedHighPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), isolatedHighPtTrackerMuons.begin(), isolatedHighPtTrackerMuons.end() );

  unsigned int nHighPtMuons = allHighPtMuons.size();

  // do electron selection
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  std::vector<pat::ElectronRef> acceptEles;

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

    if ( aEle->electronID("modifiedHeepElectronID") ) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      acceptEles.push_back(castEle);
    }
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByEt);

  if ( acceptEles.empty() )
    return;

  std::vector<pat::ElectronRef> CRMEs;
  std::vector<pat::ElectronRef> antiMEs;

  auto estimateU5x5Et = [&union5x5dEtaInHandle,&union5x5EnergyHandle] (const pat::ElectronRef& aEle) -> float {
    const float eta1stGSF = -(aEle->deltaEtaSeedClusterTrackAtVtx() - aEle->superCluster()->seed()->eta());
    const float u5x5Eta = (*union5x5dEtaInHandle)[aEle] + eta1stGSF;
    const float u5x5Et = (*union5x5EnergyHandle)[aEle]/std::cosh(u5x5Eta);

    return u5x5Et;
  };

  for (unsigned int idx = 0; idx < acceptEles.size(); ++idx) {
    const auto aEle = acceptEles.at(idx);

    const auto& orgGsfTrk = aEle->gsfTrack();
    const auto& addGsfTrk = (*addGsfTrkMap)[aEle];

    const float u5x5Et = estimateU5x5Et(aEle);

    if ( u5x5Et < 50. )
      continue;

    bool isNotMerged = MergedLeptonHelperFct::isNotMerged(aEle,eleHandle,addGsfTrk);

    if ( aEle->electronID("mvaMergedElectron") && !isNotMerged )
      CRMEs.push_back(aEle);

    // electrons with mva score = -1 are out of our interest
    if ( !static_cast<bool>(aEle->electronID("mvaMergedElectron")) && aEle->userFloat("mvaMergedElectronValues")!=-1. && !isNotMerged )
      antiMEs.push_back(aEle);
  } // electron loop

  std::sort(CRMEs.begin(),CRMEs.end(),sortByEt);
  std::sort(antiMEs.begin(),antiMEs.end(),sortByEt);

  // recoSF
  if (isMC_) {
    for (const auto& aEle : acceptEles) {
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= recoSF_->GetBinContent( recoSF_->FindBin(aEle->superCluster()->eta(), apt ) );
    }
  }

  if ( nHighPtMuons==2 && acceptEles.size()==1 ) {
    histo1d_["cutflow_2M"]->Fill( 5.5, aWeight );

    const pat::MuonRef& aMu = allHighPtMuons.at(0);
    const pat::MuonRef& bMu = allHighPtMuons.at(1);

    double momcorr1 = 1.;
    double momcorr2 = 1.;

    if (isMC_) {
      edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
      iEvent.getByToken(genptcToken_, genptcHandle);

      momcorr1 = mucorrHelper_.rochesterMC(aMu,genptcHandle);
      momcorr2 = mucorrHelper_.rochesterMC(bMu,genptcHandle);
    } else {
      momcorr1 = mucorrHelper_.rochesterData(aMu);
      momcorr2 = mucorrHelper_.rochesterData(bMu);
    }

    const auto lvecMa = math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                                     aMu->tunePMuonBestTrack()->eta(),
                                                     aMu->tunePMuonBestTrack()->phi(),
                                                     mumass_)*momcorr1;
    const auto lvecMb = math::PtEtaPhiMLorentzVector(bMu->tunePMuonBestTrack()->pt(),
                                                     bMu->tunePMuonBestTrack()->eta(),
                                                     bMu->tunePMuonBestTrack()->phi(),
                                                     mumass_)*momcorr2;

    if ( CRMEs.size()==1 && antiMEs.size()==0 ) {
      const auto lvecCRME = CRMEs.front()->polarP4()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy();
      const auto lvecCRllME = lvecCRME + lvecMa + lvecMb;
      const double dr2l1ME = reco::deltaR2(CRMEs.front()->eta(),CRMEs.front()->phi(),lvecMa.eta(),lvecMa.phi());
      const double dr2l2ME = reco::deltaR2(CRMEs.front()->eta(),CRMEs.front()->phi(),lvecMb.eta(),lvecMb.phi());
      const double dr2l1l2 = reco::deltaR2(lvecMa.eta(),lvecMa.phi(),lvecMb.eta(),lvecMb.phi());
      const double mlll = lvecCRllME.M();

      const float u5x5Et = estimateU5x5Et(CRMEs.front());

      if ( mlll > 50. && (mlll < 200. || isMC_) )
        histo1d_["2M_CRME_lll_invM"]->Fill(mlll,aWeight);

      if ( mlll > 50. && mlll < 200. ) {
        if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
          histo1d_["2M_CRME_Et"]->Fill(CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(),aWeight);
          histo1d_["2M_CRME_eta"]->Fill(CRMEs.front()->eta(),aWeight);
          histo1d_["2M_CRME_phi"]->Fill(CRMEs.front()->phi(),aWeight);
        } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
          histo1d_["2M_CRME_noGsf_Et"]->Fill(CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(),aWeight);
          histo1d_["2M_CRME_noGsf_eta"]->Fill(CRMEs.front()->eta(),aWeight);
          histo1d_["2M_CRME_noGsf_phi"]->Fill(CRMEs.front()->phi(),aWeight);
        } else {}

        histo1d_["2M_CRME_all_eta"]->Fill(CRMEs.front()->eta(),aWeight);

        histo1d_["2M_CRME_M1_pt"]->Fill(lvecMa.pt(),aWeight);
        histo1d_["2M_CRME_M1_eta"]->Fill(lvecMa.eta(),aWeight);
        histo1d_["2M_CRME_M1_phi"]->Fill(lvecMa.phi(),aWeight);

        histo1d_["2M_CRME_M2_pt"]->Fill(lvecMb.pt(),aWeight);
        histo1d_["2M_CRME_M2_eta"]->Fill(lvecMb.eta(),aWeight);
        histo1d_["2M_CRME_M2_phi"]->Fill(lvecMb.phi(),aWeight);

        histo1d_["2M_CRME_lll_pt"]->Fill(lvecCRllME.pt(),aWeight);
        histo1d_["2M_CRME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
        histo1d_["2M_CRME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
        histo1d_["2M_CRME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

        if ( CRMEs.front()->isEB() )
          histo1d_["2M_Et_CR_EB_CRME"]->Fill( u5x5Et, aWeight );
      }
    } // CRMEs.size()==1

    if ( CRMEs.size()==0 && antiMEs.size()==1 ) {
      const auto lvecAME = antiMEs.front()->polarP4()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy();
      const auto lvecAllME = lvecAME + lvecMa + lvecMb;
      const double dr2l1ME = reco::deltaR2(antiMEs.front()->eta(),antiMEs.front()->phi(),lvecMa.eta(),lvecMa.phi());
      const double dr2l2ME = reco::deltaR2(antiMEs.front()->eta(),antiMEs.front()->phi(),lvecMb.eta(),lvecMb.phi());
      const double dr2l1l2 = reco::deltaR2(lvecMa.eta(),lvecMa.phi(),lvecMb.eta(),lvecMb.phi());
      const double mlll = lvecAllME.M();

      const float u5x5Et = estimateU5x5Et(antiMEs.front());

      if (mlll > 50.) {
        if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
          histo1d_["2M_antiME_Et"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight);
          histo1d_["2M_antiME_eta"]->Fill(antiMEs.front()->eta(),aWeight);
          histo1d_["2M_antiME_phi"]->Fill(antiMEs.front()->phi(),aWeight);
        } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
          histo1d_["2M_antiME_noGsf_Et"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight);
          histo1d_["2M_antiME_noGsf_eta"]->Fill(antiMEs.front()->eta(),aWeight);
          histo1d_["2M_antiME_noGsf_phi"]->Fill(antiMEs.front()->phi(),aWeight);
        } else {}

        histo1d_["2M_antiME_M1_pt"]->Fill(lvecMa.pt(),aWeight);
        histo1d_["2M_antiME_M1_eta"]->Fill(lvecMa.eta(),aWeight);
        histo1d_["2M_antiME_M1_phi"]->Fill(lvecMa.phi(),aWeight);

        histo1d_["2M_antiME_M2_pt"]->Fill(lvecMb.pt(),aWeight);
        histo1d_["2M_antiME_M2_eta"]->Fill(lvecMb.eta(),aWeight);
        histo1d_["2M_antiME_M2_phi"]->Fill(lvecMb.phi(),aWeight);

        histo1d_["2M_antiME_lll_pt"]->Fill(lvecAllME.pt(),aWeight);
        histo1d_["2M_antiME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
        histo1d_["2M_antiME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
        histo1d_["2M_antiME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

        TF1* osfunc = osboth_;
        const TF1* ssfunc = ssboth_;
        int ibin = osPar_->FindBin(antiMEs.front()->superCluster()->eta());
        osfunc->SetParameter(2,osPar_->GetBinContent(ibin));
        osfunc->SetParError(2,osPar_->GetBinError(ibin));
        const double ffOS = osfunc->Eval(u5x5Et);
        const double ffSS = ssfunc->Eval(u5x5Et);
        const double xval[1] = {u5x5Et};
        double ciOS[1], ciSS[1];
        const TFitResultPtr osFit = osFit_;
        osFit->GetConfidenceIntervals(1,1,0,xval,ciOS,0.95,false);
        ssFit_->GetConfidenceIntervals(1,1,0,xval,ciSS,0.95,false);
        histo1d_["2M_antiME_lll_invM_CR"]->Fill(mlll,aWeight);
        histo1d_["2M_antiME_lll_invM_CR_xOSFF"]->Fill(mlll,aWeight*ffOS);
        histo1d_["2M_antiME_lll_invM_CR_xSSFF"]->Fill(mlll,aWeight*ffSS);
        histo1d_["2M_antiME_lll_invM_CR_xOSFF_up"]->Fill(mlll,aWeight*(ffOS+ciOS[0]));
        histo1d_["2M_antiME_lll_invM_CR_xSSFF_up"]->Fill(mlll,aWeight*(ffSS+ciSS[0]));
        histo1d_["2M_antiME_lll_invM_CR_xOSFF_dn"]->Fill(mlll,aWeight*(ffOS-ciOS[0]));
        histo1d_["2M_antiME_lll_invM_CR_xSSFF_dn"]->Fill(mlll,aWeight*(ffSS-ciSS[0]));

        if ( mlll < 200. ) {
          histo1d_["2M_antiME_eta_xSSFF"]->Fill(antiMEs.front()->eta(),aWeight*ffSS);
          histo1d_["2M_antiME_eta_xOSFF"]->Fill(antiMEs.front()->eta(),aWeight*ffOS);
          histo1d_["2M_antiME_eta_xSSFF_up"]->Fill(antiMEs.front()->eta(),aWeight*(ffSS+ciSS[0]));
          histo1d_["2M_antiME_eta_xOSFF_up"]->Fill(antiMEs.front()->eta(),aWeight*(ffOS+ciOS[0]));
          histo1d_["2M_antiME_eta_xSSFF_dn"]->Fill(antiMEs.front()->eta(),aWeight*(ffSS-ciSS[0]));
          histo1d_["2M_antiME_eta_xOSFF_dn"]->Fill(antiMEs.front()->eta(),aWeight*(ffOS-ciOS[0]));

          if ( antiMEs.front()->isEB() ) {
            histo1d_["2M_Et_CR_EB_antiME"]->Fill( u5x5Et, aWeight );
            histo1d_["2M_Et_CR_EB_antiME_xSSFF"]->Fill( u5x5Et, aWeight*ffSS);
            histo1d_["2M_Et_CR_EB_antiME_xSSFF_up"]->Fill( u5x5Et, aWeight*(ffSS+ciSS[0]));
            histo1d_["2M_Et_CR_EB_antiME_xSSFF_dn"]->Fill( u5x5Et, aWeight*(ffSS-ciSS[0]));
            histo1d_["2M_Et_CR_EB_antiME_xOSFF"]->Fill( u5x5Et, aWeight*ffOS);
            histo1d_["2M_Et_CR_EB_antiME_xOSFF_up"]->Fill( u5x5Et, aWeight*(ffOS+ciOS[0]));
            histo1d_["2M_Et_CR_EB_antiME_xOSFF_dn"]->Fill( u5x5Et, aWeight*(ffOS-ciOS[0]));
          }
        } // CR
      } // mlll > 50.
    } // antiMEs.size()==1
  } // nHighPtMuons==2 && acceptEles.size()==1

  if ( nHighPtMuons==1 && acceptEles.size()==1 ) {
    histo1d_["cutflow_1M"]->Fill( 5.5, aWeight );

    const pat::MuonRef& aMu = allHighPtMuons.front();

    double momcorr1 = 1.;

    if (isMC_) {
      edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
      iEvent.getByToken(genptcToken_, genptcHandle);

      momcorr1 = mucorrHelper_.rochesterMC(aMu,genptcHandle);
    } else {
      momcorr1 = mucorrHelper_.rochesterData(aMu);
    }

    const auto lvecMa = math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                                     aMu->tunePMuonBestTrack()->eta(),
                                                     aMu->tunePMuonBestTrack()->phi(),
                                                     mumass_)*momcorr1;

    const auto& aMET = metHandle->at(0);

    if ( CRMEs.size()==1 && antiMEs.size()==0 ) {
      histo1d_["cutflow_1M"]->Fill( 6.5, aWeight );

      const auto lvecCRME = CRMEs.front()->polarP4()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy();
      const auto lvecCRll = lvecCRME + lvecMa;

      const auto tpMET = aMET.p4() - lvecMa - lvecCRME
                                   + aMu->p4() + CRMEs.front()->p4();
      const auto lvecCRllMET = lvecCRME + lvecMa + tpMET;
      const double mtllMET = lvecCRllMET.mt();

      double ffMM = ffFunc_->Eval(mtllMET);
      ffMM = ffMM > 0. ? ffMM : 0.;
      const double xvalMM[1] = {mtllMET};
      double ciMM[1];
      fitResult_->GetConfidenceIntervals(1,1,0,xvalMM,ciMM,0.95,false);

      // const float u5x5Et = estimateU5x5Et(CRMEs.front());

      if ( tpMET.pt() > ptThres_ && lvecCRll.M() > 50. ) {
        histo1d_["cutflow_1M"]->Fill( 7.5, aWeight );

        double dphi = reco::deltaPhi(lvecMa.phi(),tpMET.phi());
        bool passDPhi = std::abs(dphi) < drThres_;
        bool passDPhiCR = (std::abs(dphi) < drThresCR_) && !passDPhi;
        bool passRatioPt = false;

        double ratioPt = ( lvecCRME + lvecMa ).pt() / tpMET.pt(); // ( tpMET.pt() + lvecMa.pt() ) / lvecCRME.pt();

        if ( ratioThresLo_ < ratioPt && ratioPt < ratioThresHi_ )
          passRatioPt = true;

        histo1d_["1M_ABCD_MET_dphi"]->Fill(dphi, aWeight);
        histo1d_["1M_ABCD_MET_ratioPt"]->Fill(ratioPt, aWeight);

        histo2d_["1M_mt_dphi"]->Fill(mtllMET,dphi,aWeight);
        histo2d_["1M_mt_ratioPt"]->Fill(mtllMET,ratioPt,aWeight);
        histo2d_["1M_dphi_ratioPt"]->Fill(dphi,ratioPt,aWeight);

        if (passDPhi)
          histo1d_["cutflow_1M"]->Fill( 8.5, aWeight );

        if ( passDPhi && passRatioPt ) {
          histo1d_["cutflow_1M"]->Fill( 9.5, aWeight );

          if (mtllMET < 250. || isMC_) { // blinded
            histo1d_["1M_CRME_MM_pt"]->Fill( lvecMa.pt(), aWeight );
            histo1d_["1M_CRME_MM_eta"]->Fill( lvecMa.eta(), aWeight );
            histo1d_["1M_CRME_MM_phi"]->Fill( lvecMa.phi(), aWeight );

            histo1d_["1M_CRME_MET_pt"]->Fill( tpMET.pt(), aWeight );
            histo1d_["1M_CRME_MET_phi"]->Fill( tpMET.phi(), aWeight );

            if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
              histo1d_["1M_CRME_Et"]->Fill(lvecCRME.pt(),aWeight);
              histo1d_["1M_CRME_eta"]->Fill(lvecCRME.eta(),aWeight);
              histo1d_["1M_CRME_phi"]->Fill(lvecCRME.phi(),aWeight);
            } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
              histo1d_["1M_CRME_noGsf_Et"]->Fill(lvecCRME.pt(),aWeight);
              histo1d_["1M_CRME_noGsf_eta"]->Fill(lvecCRME.eta(),aWeight);
              histo1d_["1M_CRME_noGsf_phi"]->Fill(lvecCRME.phi(),aWeight);
            } else {}

            histo1d_["1M_CRME_all_Et"]->Fill(lvecCRME.pt(),aWeight);
            histo1d_["1M_CRME_all_eta"]->Fill(lvecCRME.eta(),aWeight);

            histo1d_["1M_CRME_MET_dphi"]->Fill( dphi, aWeight );
            histo1d_["1M_CRME_ratioPt"]->Fill( ratioPt, aWeight );
            histo1d_["1M_CRME_mt"]->Fill( mtllMET, aWeight );

            histo1d_["1M_FF_CRME_MMMET_pt"]->Fill( (tpMET+lvecMa).pt(), aWeight );
            histo1d_["1M_FF_CRME_MET_pt"]->Fill( tpMET.pt(), aWeight );
          }
        } else if ( passDPhi && !passRatioPt ) {
          histo1d_["1M_CRME_antiRpt_MM_pt"]->Fill( lvecMa.pt(), aWeight );
          histo1d_["1M_CRME_antiRpt_MM_pt_xFF"]->Fill( lvecMa.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_antiRpt_MM_eta"]->Fill( lvecMa.eta(), aWeight );
          histo1d_["1M_CRME_antiRpt_all_Et"]->Fill( lvecCRME.pt(), aWeight );
          histo1d_["1M_CRME_antiRpt_all_Et_xFF"]->Fill( lvecCRME.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_antiRpt_all_eta"]->Fill( lvecCRME.eta(), aWeight );
          histo1d_["1M_CRME_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["1M_CRME_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_antiRpt_MET_phi"]->Fill( tpMET.phi(), aWeight );
          histo1d_["1M_CRME_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["1M_CRME_antiRpt_ratioPt"]->Fill( ratioPt, aWeight );
          histo1d_["1M_CRME_antiRpt_mt"]->Fill( mtllMET, aWeight );
          histo1d_["1M_CRME_antiRpt_mt_xFF"]->Fill( mtllMET, aWeight*ffMM );
          histo1d_["1M_CRME_antiRpt_mt_xFF_up"]->Fill( mtllMET, aWeight*(ffMM+ciMM[0]) );
          histo1d_["1M_CRME_antiRpt_mt_xFF_dn"]->Fill( mtllMET, aWeight*(ffMM-ciMM[0]) );
        } else if ( passDPhiCR && passRatioPt ) {
          histo1d_["1M_CRME_CRdphi_MM_pt"]->Fill( lvecMa.pt(), aWeight );
          histo1d_["1M_CRME_CRdphi_MM_eta"]->Fill( lvecMa.eta(), aWeight );
          histo1d_["1M_CRME_CRdphi_MM_phi"]->Fill( lvecMa.phi(), aWeight );

          histo1d_["1M_CRME_CRdphi_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["1M_CRME_CRdphi_MET_phi"]->Fill( tpMET.phi(), aWeight );

          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["1M_CRME_CRdphi_Et"]->Fill(lvecCRME.pt(),aWeight);
            histo1d_["1M_CRME_CRdphi_eta"]->Fill(lvecCRME.eta(),aWeight);
            histo1d_["1M_CRME_CRdphi_phi"]->Fill(lvecCRME.phi(),aWeight);
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["1M_CRME_CRdphi_noGsf_Et"]->Fill(lvecCRME.pt(),aWeight);
            histo1d_["1M_CRME_CRdphi_noGsf_eta"]->Fill(lvecCRME.eta(),aWeight);
            histo1d_["1M_CRME_CRdphi_noGsf_phi"]->Fill(lvecCRME.phi(),aWeight);
          } else {}

          histo1d_["1M_CRME_CRdphi_all_Et"]->Fill(lvecCRME.pt(),aWeight);
          histo1d_["1M_CRME_CRdphi_all_eta"]->Fill(lvecCRME.eta(),aWeight);

          histo1d_["1M_CRME_CRdphi_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["1M_CRME_CRdphi_ratioPt"]->Fill( ratioPt, aWeight );
          histo1d_["1M_CRME_CRdphi_mt"]->Fill( mtllMET, aWeight );
        } else if ( passDPhiCR && !passRatioPt ) {
          histo1d_["1M_CRME_CRdphi_antiRpt_MM_pt"]->Fill( lvecMa.pt(), aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_MM_pt_xFF"]->Fill( lvecMa.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_CRdphi_antiRpt_MM_eta"]->Fill( lvecMa.eta(), aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_all_Et"]->Fill( lvecCRME.pt(), aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_all_Et_xFF"]->Fill( lvecCRME.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_CRdphi_antiRpt_all_eta"]->Fill( lvecCRME.eta(), aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_CRdphi_antiRpt_MET_phi"]->Fill( tpMET.phi(), aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_ratioPt"]->Fill( ratioPt, aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_mt"]->Fill( mtllMET, aWeight );
          histo1d_["1M_CRME_CRdphi_antiRpt_mt_xFF"]->Fill( mtllMET, aWeight*ffMM );
          histo1d_["1M_CRME_CRdphi_antiRpt_mt_xFF_up"]->Fill( mtllMET, aWeight*(ffMM+ciMM[0]) );
          histo1d_["1M_CRME_CRdphi_antiRpt_mt_xFF_dn"]->Fill( mtllMET, aWeight*(ffMM-ciMM[0]) );
        } else if ( !passDPhi && !passDPhiCR && !passRatioPt ) {
          histo1d_["1M_CRME_antiDphi_antiRpt_MM_pt"]->Fill( lvecMa.pt(), aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_MM_pt_xFF"]->Fill( lvecMa.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_antiDphi_antiRpt_MM_eta"]->Fill( lvecMa.eta(), aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_all_Et"]->Fill( lvecCRME.pt(), aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_all_Et_xFF"]->Fill( lvecCRME.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_antiDphi_antiRpt_all_eta"]->Fill( lvecCRME.eta(), aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["1M_CRME_antiDphi_antiRpt_MET_phi"]->Fill( tpMET.phi(), aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_ratioPt"]->Fill( ratioPt, aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_mt"]->Fill( mtllMET, aWeight );
          histo1d_["1M_CRME_antiDphi_antiRpt_mt_xFF"]->Fill( mtllMET, aWeight*ffMM );
          histo1d_["1M_CRME_antiDphi_antiRpt_mt_xFF_up"]->Fill( mtllMET, aWeight*(ffMM+ciMM[0]) );
          histo1d_["1M_CRME_antiDphi_antiRpt_mt_xFF_dn"]->Fill( mtllMET, aWeight*(ffMM-ciMM[0]) );
        } else {
          histo1d_["1M_CRME_antiDphi_MM_pt"]->Fill( lvecMa.pt(), aWeight );
          histo1d_["1M_CRME_antiDphi_MM_eta"]->Fill( lvecMa.eta(), aWeight );
          histo1d_["1M_CRME_antiDphi_MM_phi"]->Fill( lvecMa.phi(), aWeight );

          histo1d_["1M_CRME_antiDphi_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["1M_CRME_antiDphi_MET_phi"]->Fill( tpMET.phi(), aWeight );

          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["1M_CRME_antiDphi_Et"]->Fill(lvecCRME.pt(),aWeight);
            histo1d_["1M_CRME_antiDphi_eta"]->Fill(lvecCRME.eta(),aWeight);
            histo1d_["1M_CRME_antiDphi_phi"]->Fill(lvecCRME.phi(),aWeight);
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["1M_CRME_antiDphi_noGsf_Et"]->Fill(lvecCRME.pt(),aWeight);
            histo1d_["1M_CRME_antiDphi_noGsf_eta"]->Fill(lvecCRME.eta(),aWeight);
            histo1d_["1M_CRME_antiDphi_noGsf_phi"]->Fill(lvecCRME.phi(),aWeight);
          } else {}

          histo1d_["1M_CRME_antiDphi_all_Et"]->Fill(lvecCRME.pt(),aWeight);
          histo1d_["1M_CRME_antiDphi_all_eta"]->Fill(lvecCRME.eta(),aWeight);

          histo1d_["1M_CRME_antiDphi_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["1M_CRME_antiDphi_ratioPt"]->Fill( ratioPt, aWeight );
          histo1d_["1M_CRME_antiDphi_mt"]->Fill( mtllMET, aWeight );

          if ( mtllMET < 250. ) {
            histo1d_["1M_FF_CRME_antiDphi_MMMET_pt"]->Fill( (tpMET+lvecMa).pt(), aWeight );
            histo1d_["1M_FF_CRME_antiDphi_MET_pt"]->Fill( tpMET.pt(), aWeight );
          }
        } // MET dPhi ratioPt CR
      } // MET ptThres
    } // CRMEs.size()==1 && antiMEs.size()==0
  } // nHighPtMuons==1 && acceptEles.size()==1

  return;
}

DEFINE_FWK_MODULE(MergedEMuCRanalyzer);
