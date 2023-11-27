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

#include "TF1.h"
#include "TFitResult.h"

class MergedMuCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedMuCRanalyzer(const edm::ParameterSet&);
  ~MergedMuCRanalyzer() override {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  bool isMC_;

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

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const std::string trigHistName_;
  const std::string idHistName_;
  const std::string idHistNameTrkHighPt_;
  const std::string isoHistName_;
  const std::string isoHistNameTrkHighPt_;
  const edm::FileInPath MMFFpath_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath muonIdSFpath_;
  const edm::FileInPath muonIsoSFpath_;
  const edm::FileInPath purwgtPath_;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  const double ptThres_;
  const double ptMuThres_;
  const double drThres_;
  const double drThresCR_;
  const double ratioThresLo_;
  const double ratioThresHi_;
  const double mumass_ = 0.1056583745;

  MuonCorrectionHelper mucorrHelper_;

  std::unique_ptr<TFile> MMFFfile_;
  TF1* ffFunc_;
  TFitResultPtr fitResult_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;
};

MergedMuCRanalyzer::MergedMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
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
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
trigHistName_(iConfig.getParameter<std::string>("trigHistName")),
idHistName_(iConfig.getParameter<std::string>("idHistName")),
idHistNameTrkHighPt_(iConfig.getParameter<std::string>("idHistNameTrkHighPt")),
isoHistName_(iConfig.getParameter<std::string>("isoHistName")),
isoHistNameTrkHighPt_(iConfig.getParameter<std::string>("isoHistNameTrkHighPt")),
MMFFpath_(iConfig.getParameter<edm::FileInPath>("MMFFpath")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
muonIdSFpath_(iConfig.getParameter<edm::FileInPath>("muonIdSFpath")),
muonIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonIsoSFpath")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptMuThres_(iConfig.getParameter<double>("ptMuThres")),
drThres_(iConfig.getParameter<double>("drThres")),
drThresCR_(iConfig.getParameter<double>("drThresCR")),
ratioThresLo_(iConfig.getParameter<double>("ratioThresLo")),
ratioThresHi_(iConfig.getParameter<double>("ratioThresHi")),
mucorrHelper_(rochesterPath_,triggerSFpath_,muonIdSFpath_,muonIsoSFpath_,trigHistName_) {
  usesResource("TFileService");
}

void MergedMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  MMFFfile_ = std::make_unique<TFile>(MMFFpath_.fullPath().c_str(),"READ");
  ffFunc_ = static_cast<TF1*>(MMFFfile_->FindObjectAny("MMFF"));
  fitResult_ = (static_cast<TH1D*>(MMFFfile_->Get("1M_MMFF_numer_rebin")))->Fit(ffFunc_,"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["totWeightedSum_4M"] = fs->make<TH1D>("totWeightedSum_4M","totWeightedSum_4M",1,0.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);

  histo1d_["cutflow_4M"] = fs->make<TH1D>("cutflow_4M","cutflow",10,0.,10.);
  histo1d_["cutflow_3M"] = fs->make<TH1D>("cutflow_3M","cutflow",15,0.,15.);

  histo1d_["3M_MM_pt"] = fs->make<TH1D>("3M_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_MM_eta"] = fs->make<TH1D>("3M_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_MM_phi"] = fs->make<TH1D>("3M_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_M1_pt"] = fs->make<TH1D>("3M_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_M1_eta"] = fs->make<TH1D>("3M_M1_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_M1_phi"] = fs->make<TH1D>("3M_M1_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_M2_pt"] = fs->make<TH1D>("3M_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_M2_eta"] = fs->make<TH1D>("3M_M2_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_M2_phi"] = fs->make<TH1D>("3M_M2_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_MET_pt"] = fs->make<TH1D>("3M_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_MET_phi"] = fs->make<TH1D>("3M_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["3M_MET_dphi"] = fs->make<TH1D>("3M_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_mt"] = fs->make<TH1D>("3M_mt","m_{T};m_{T};",500,0.,2500.);

  histo1d_["3M_antiRpt_MM_pt"] = fs->make<TH1D>("3M_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("3M_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_M1_pt"] = fs->make<TH1D>("3M_antiRpt_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_M2_pt"] = fs->make<TH1D>("3M_antiRpt_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MET_pt"] = fs->make<TH1D>("3M_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("3M_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiRpt_MET_dphi"] = fs->make<TH1D>("3M_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_antiRpt_mt"] = fs->make<TH1D>("3M_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF"] = fs->make<TH1D>("3M_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF_up"] = fs->make<TH1D>("3M_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("3M_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_MET_ratioPt"] = fs->make<TH1D>("3M_antiDphi_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);

  histo1d_["3M_CRdphi_MM_pt"] = fs->make<TH1D>("3M_CRdphi_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_MM_eta"] = fs->make<TH1D>("3M_CRdphi_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_CRdphi_MM_phi"] = fs->make<TH1D>("3M_CRdphi_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_CRdphi_M1_pt"] = fs->make<TH1D>("3M_CRdphi_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_M1_eta"] = fs->make<TH1D>("3M_CRdphi_M1_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_CRdphi_M1_phi"] = fs->make<TH1D>("3M_CRdphi_M1_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_CRdphi_M2_pt"] = fs->make<TH1D>("3M_CRdphi_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_M2_eta"] = fs->make<TH1D>("3M_CRdphi_M2_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_CRdphi_M2_phi"] = fs->make<TH1D>("3M_CRdphi_M2_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_CRdphi_MET_pt"] = fs->make<TH1D>("3M_CRdphi_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_MET_phi"] = fs->make<TH1D>("3M_CRdphi_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["3M_CRdphi_MET_dphi"] = fs->make<TH1D>("3M_CRdphi_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_CRdphi_mt"] = fs->make<TH1D>("3M_CRdphi_mt","m_{T};m_{T};",500,0.,2500.);

  histo1d_["3M_CRdphi_antiRpt_MM_pt"] = fs->make<TH1D>("3M_CRdphi_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("3M_CRdphi_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_antiRpt_M1_pt"] = fs->make<TH1D>("3M_CRdphi_antiRpt_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_antiRpt_M2_pt"] = fs->make<TH1D>("3M_CRdphi_antiRpt_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_antiRpt_MET_pt"] = fs->make<TH1D>("3M_CRdphi_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("3M_CRdphi_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_CRdphi_antiRpt_MET_dphi"] = fs->make<TH1D>("3M_CRdphi_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_CRdphi_antiRpt_mt"] = fs->make<TH1D>("3M_CRdphi_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_CRdphi_antiRpt_mt_xFF"] = fs->make<TH1D>("3M_CRdphi_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_CRdphi_antiRpt_mt_xFF_up"] = fs->make<TH1D>("3M_CRdphi_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_CRdphi_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("3M_CRdphi_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_CRdphi_antiDphi_MET_ratioPt"] = fs->make<TH1D>("3M_CRdphi_antiDphi_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);

  histo1d_["3M_antiDphi_MM_pt"] = fs->make<TH1D>("3M_antiDphi_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_MM_eta"] = fs->make<TH1D>("3M_antiDphi_MM_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_antiDphi_MM_phi"] = fs->make<TH1D>("3M_antiDphi_MM_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_antiDphi_M1_pt"] = fs->make<TH1D>("3M_antiDphi_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_M1_eta"] = fs->make<TH1D>("3M_antiDphi_M1_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_antiDphi_M1_phi"] = fs->make<TH1D>("3M_antiDphi_M1_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_antiDphi_M2_pt"] = fs->make<TH1D>("3M_antiDphi_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_M2_eta"] = fs->make<TH1D>("3M_antiDphi_M2_eta","3Eta;#eta;",200,-2.5,2.5);
  histo1d_["3M_antiDphi_M2_phi"] = fs->make<TH1D>("3M_antiDphi_M2_phi","Phi;#phi;",128,-3.2,3.2);

  histo1d_["3M_antiDphi_MET_pt"] = fs->make<TH1D>("3M_antiDphi_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_MET_phi"] = fs->make<TH1D>("3M_antiDphi_MET_phi","Phi;#phi;",128,-3.2,3.2);
  histo1d_["3M_antiDphi_MET_dphi"] = fs->make<TH1D>("3M_antiDphi_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_antiDphi_MET_ratioPt"] = fs->make<TH1D>("3M_antiDphi_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);
  histo1d_["3M_antiDphi_mt"] = fs->make<TH1D>("3M_antiDphi_mt","m_{T};m_{T};",500,0.,2500.);

  histo1d_["3M_antiDphi_M1M2_dR"] = fs->make<TH1D>("3M_antiDphi_M1M2_dR","#Delta R(M1,M2);#Delta R;",100,0.,0.5);
  histo1d_["3M_antiDphi_M1MM_dR"] = fs->make<TH1D>("3M_antiDphi_M1MM_dR","#Delta R(M1,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_antiDphi_M2MM_dR"] = fs->make<TH1D>("3M_antiDphi_M2MM_dR","#Delta R(M2,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_antiDphi_M1M2_invM"] = fs->make<TH1D>("3M_antiDphi_M1M2_invM","M(M1,M2);M(M1,M2);",400,0.,200.);
  histo1d_["3M_antiDphi_M1M2_invM_zoomed"] = fs->make<TH1D>("3M_antiDphi_M1M2_invM_zoomed","M(M1,M2);M(M1,M2);",400,0.,20.);
  histo1d_["3M_antiDphi_M1M2MM_invM"] = fs->make<TH1D>("3M_antiDphi_M1M2MM_invM","M;GeV;",400,0.,400);

  histo1d_["3M_antiDphi_antiRpt_MM_pt"] = fs->make<TH1D>("3M_antiDphi_antiRpt_MM_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_antiRpt_MM_pt_xFF"] = fs->make<TH1D>("3M_antiDphi_antiRpt_MM_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_antiRpt_M1_pt"] = fs->make<TH1D>("3M_antiDphi_antiRpt_M1_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_antiRpt_M2_pt"] = fs->make<TH1D>("3M_antiDphi_antiRpt_M2_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_antiRpt_MET_pt"] = fs->make<TH1D>("3M_antiDphi_antiRpt_MET_pt","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_antiRpt_MET_pt_xFF"] = fs->make<TH1D>("3M_antiDphi_antiRpt_MET_pt_xFF","Pt;p_{T};",200,0.,500.);
  histo1d_["3M_antiDphi_antiRpt_MET_dphi"] = fs->make<TH1D>("3M_antiDphi_antiRpt_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_antiDphi_antiRpt_mt"] = fs->make<TH1D>("3M_antiDphi_antiRpt_mt","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_antiRpt_mt_xFF"] = fs->make<TH1D>("3M_antiDphi_antiRpt_mt_xFF","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_antiRpt_mt_xFF_up"] = fs->make<TH1D>("3M_antiDphi_antiRpt_mt_xFF_up","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_antiRpt_mt_xFF_dn"] = fs->make<TH1D>("3M_antiDphi_antiRpt_mt_xFF_dn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_antiDphi_MET_ratioPt"] = fs->make<TH1D>("3M_antiDphi_antiDphi_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);

  histo1d_["3M_ABCD_MET_dphi"] = fs->make<TH1D>("3M_ABCD_MET_dphi","dPhi;#Delta#phi;",128,-3.2,3.2);
  histo1d_["3M_ABCD_MET_ratioPt"] = fs->make<TH1D>("3M_ABCD_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);

  histo1d_["3M_check_resolved_M1M2_dR"] = fs->make<TH1D>("3M_check_resolved_M1M2_dR","#Delta R(M1,M2);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_resolved_M1MM_dR"] = fs->make<TH1D>("3M_check_resolved_M1MM_dR","#Delta R(M1,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_resolved_M2MM_dR"] = fs->make<TH1D>("3M_check_resolved_M2MM_dR","#Delta R(M2,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_resolved_M1M2_invM"] = fs->make<TH1D>("3M_check_resolved_M1M2_invM","M(M1,M2);M(M1,M2);",400,0.,400.);
  histo1d_["3M_check_resolved_M1M2_invM_zoomed"] = fs->make<TH1D>("3M_check_resolved_M1M2_invM_zoomed","M(M1,M2);M(M1,M2);",400,0.,20.);
  histo1d_["3M_check_resolved_M1M2MM_invM"] = fs->make<TH1D>("3M_check_resolved_M1M2MM_invM","M;GeV;",400,0.,400);

  histo1d_["3M_check_M1M2_dR"] = fs->make<TH1D>("3M_check_M1M2_dR","#Delta R(M1,M2);#Delta R;",100,0.,0.5);
  histo1d_["3M_check_M1MM_dR"] = fs->make<TH1D>("3M_check_M1MM_dR","#Delta R(M1,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_M2MM_dR"] = fs->make<TH1D>("3M_check_M2MM_dR","#Delta R(M2,MM);#Delta R;",128,0.,6.4);
  histo1d_["3M_check_M1M2_invM"] = fs->make<TH1D>("3M_check_M1M2_invM","M(M1,M2);M(M1,M2);",400,0.,200.);
  histo1d_["3M_check_M1M2_invM_zoomed"] = fs->make<TH1D>("3M_check_M1M2_invM_zoomed","M(M1,M2);M(M1,M2);",400,0.,20.);
  histo1d_["3M_check_M1M2MM_invM"] = fs->make<TH1D>("3M_check_M1M2MM_invM","M;GeV;",400,0.,400);
  histo1d_["3M_check_M1_iso"] = fs->make<TH1D>("3M_check_M1_iso","M1 Track iso;#Sigma p_{T iso};",400,0.,100.);
  histo1d_["3M_check_M2_iso"] = fs->make<TH1D>("3M_check_M2_iso","M2 Track iso;#Sigma p_{T iso};",400,0.,100.);
  histo1d_["3M_check_MM_iso"] = fs->make<TH1D>("3M_check_MM_iso","MM Track iso;#Sigma p_{T iso};",400,0.,100.);
  histo1d_["3M_check_M1_iso_musubtract"] = fs->make<TH1D>("3M_check_M1_iso_musubtract","M1 Track iso - M2 p_{T};#Sigma p_{T iso} - p_{T #mu};",400,0.,100.);
  histo1d_["3M_check_M2_iso_musubtract"] = fs->make<TH1D>("3M_check_M2_iso_musubtract","M2 Track iso - M1 p_{T};#Sigma p_{T iso} - p_{T #mu};",400,0.,100.);

  histo2d_["3M_mt_dphi"] = fs->make<TH2D>("3M_mt_dphi","m_{T} vs #Delta#phi;m_{T};#Delta#phi",100,0.,500.,128,-3.2,3.2);
  histo2d_["3M_mt_ratioPt"] = fs->make<TH2D>("3M_mt_ratioPt","m_{T} vs p_{T} ratio;m_{T};p_{T} ratio",100,0.,500.,100,0.,5.);
  histo2d_["3M_dphi_ratioPt"] = fs->make<TH2D>("3M_dphi_ratioPt","#Delta#phi vs p_{T} ratio;#Delta#phi;p_{T} ratio",128,-3.2,3.2,100,0.,5.);
}

void MergedMuCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( ( std::abs(genptc->pdgId())==13 ) && genptc->fromHardProcessFinalState() )
        promptMuons.push_back(genptc.castTo<reco::GenParticleRef>());
    }

    if (promptMuons.size()==4)
      histo1d_["totWeightedSum_4M"]->Fill(0.5,aWeight);
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

  histo1d_["cutflow_4M"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 0.5, aWeight );

  if (!isFired)
    return;

  histo1d_["cutflow_4M"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 1.5, aWeight );

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

  histo1d_["cutflow_4M"]->Fill( 2.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 2.5, aWeight );

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

    if ( firstIso/firstMuon->pt() < 0.1 )
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

    if ( firstIso/firstMuon->pt() < 0.1 )
      isolatedHighPtTrackerMuons.push_back(firstMuon);
  }

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  if ( isolatedHighPtMuons.size() < 2 )
    return;

  if ( isolatedHighPtMuons.front()->pt() < 52. )
    return;

  histo1d_["cutflow_4M"]->Fill( 3.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 3.5, aWeight );

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

  histo1d_["cutflow_4M"]->Fill( 4.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 4.5, aWeight );

  if (isMC_) {
    for (const auto& aMu : isolatedHighPtMuons) {
      aWeight *= mucorrHelper_.idSF(aMu,idHistName_);
      aWeight *= mucorrHelper_.isoSF(aMu,isoHistName_);
    }

    for (const auto& aMu : isolatedHighPtTrackerMuons) {
      aWeight *= mucorrHelper_.idSF(aMu,idHistNameTrkHighPt_);
      aWeight *= mucorrHelper_.isoSF(aMu,isoHistNameTrkHighPt_);
    }
  }

  // concatenate(-ish) highPtMuons & highPtTrackerMuons - order is important!
  std::vector<pat::MuonRef> allHighPtMuons(isolatedHighPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), isolatedHighPtTrackerMuons.begin(), isolatedHighPtTrackerMuons.end() );

  unsigned int nHighPtMuons = allHighPtMuons.size();

  if ( nHighPtMuons==4 )
    histo1d_["cutflow_4M"]->Fill( 5.5, aWeight );

  if ( nHighPtMuons==3 ) {
    histo1d_["cutflow_3M"]->Fill( 5.5, aWeight );

    pat::MuonRef aMu = allHighPtMuons.at(0);
    pat::MuonRef bMu = allHighPtMuons.at(1);
    pat::MuonRef cMu = allHighPtMuons.at(2);
    std::vector<std::pair<pat::MuonRef,pat::MuonRef>> muPairs;
    muPairs.push_back(std::make_pair(aMu,bMu));
    muPairs.push_back(std::make_pair(bMu,cMu));
    muPairs.push_back(std::make_pair(cMu,aMu));
    const auto lvecMa = math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                                     aMu->tunePMuonBestTrack()->eta(),
                                                     aMu->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecMb = math::PtEtaPhiMLorentzVector(bMu->tunePMuonBestTrack()->pt(),
                                                     bMu->tunePMuonBestTrack()->eta(),
                                                     bMu->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecMc = math::PtEtaPhiMLorentzVector(cMu->tunePMuonBestTrack()->pt(),
                                                     cMu->tunePMuonBestTrack()->eta(),
                                                     cMu->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const double mass3M = (lvecMa+lvecMb+lvecMc).M();

    if ( mass3M < 50. )
      return;

    histo1d_["cutflow_3M"]->Fill( 6.5, aWeight );

    auto sortByDR = [] (const std::pair<pat::MuonRef,pat::MuonRef>& apair, const std::pair<pat::MuonRef,pat::MuonRef>& bpair) {
      const auto& a1 = apair.first->tunePMuonBestTrack();
      const auto& a2 = apair.second->tunePMuonBestTrack();
      const auto& b1 = bpair.first->tunePMuonBestTrack();
      const auto& b2 = bpair.second->tunePMuonBestTrack();
      return reco::deltaR2(a1->eta(),a1->phi(),a2->eta(),a2->phi()) < reco::deltaR2(b1->eta(),b1->phi(),b2->eta(),b2->phi());
    };

    std::sort(muPairs.begin(),muPairs.end(),sortByDR);
    const auto& M1M2 = muPairs.at(0);
    const auto& M2MM = muPairs.at(1);
    const auto& M1MM = muPairs.at(2);
    const auto lvecM1 = math::PtEtaPhiMLorentzVector(M1M2.first->tunePMuonBestTrack()->pt(),
                                                     M1M2.first->tunePMuonBestTrack()->eta(),
                                                     M1M2.first->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecM2 = math::PtEtaPhiMLorentzVector(M1M2.second->tunePMuonBestTrack()->pt(),
                                                     M1M2.second->tunePMuonBestTrack()->eta(),
                                                     M1M2.second->tunePMuonBestTrack()->phi(),
                                                     mumass_);
    const auto lvecM1M2 = lvecM1 + lvecM2;

    histo1d_["3M_check_resolved_M1M2_dR"]->Fill( reco::deltaR(M1M2.first->tunePMuonBestTrack()->eta(),
                                                              M1M2.first->tunePMuonBestTrack()->phi(),
                                                              M1M2.second->tunePMuonBestTrack()->eta(),
                                                              M1M2.second->tunePMuonBestTrack()->phi()) , aWeight );
    histo1d_["3M_check_resolved_M2MM_dR"]->Fill( reco::deltaR(M2MM.first->tunePMuonBestTrack()->eta(),
                                                              M2MM.first->tunePMuonBestTrack()->phi(),
                                                              M2MM.second->tunePMuonBestTrack()->eta(),
                                                              M2MM.second->tunePMuonBestTrack()->phi()) , aWeight );
    histo1d_["3M_check_resolved_M1MM_dR"]->Fill( reco::deltaR(M1MM.first->tunePMuonBestTrack()->eta(),
                                                              M1MM.first->tunePMuonBestTrack()->phi(),
                                                              M1MM.second->tunePMuonBestTrack()->eta(),
                                                              M1MM.second->tunePMuonBestTrack()->phi()) , aWeight );
    histo1d_["3M_check_resolved_M1M2_invM"]->Fill( lvecM1M2.M() , aWeight );
    histo1d_["3M_check_resolved_M1M2_invM_zoomed"]->Fill( lvecM1M2.M() , aWeight );
    histo1d_["3M_check_resolved_M1M2MM_invM"]->Fill( mass3M , aWeight );

    double drThres2 = drThres_*drThres_;
    std::vector<std::pair<pat::MuonRef,pat::MuonRef>> cands;

    // find a collimated pair of muons
    for (unsigned idx = 0; idx < isolatedHighPtMuons.size(); idx++) { // one of them should be global
      const auto& cand1 = allHighPtMuons.at(idx);

      for (unsigned jdx = idx+1; jdx < allHighPtMuons.size(); jdx++) {
        const auto& cand2 = allHighPtMuons.at(jdx);

        // they shouldn't be the same, throw exception otherwise
        if ( cand1==cand2 )
          throw cms::Exception("LogicError") << "Error: MergedMuCRanalyzer::analyze - attempting to compare the same muons " << std::endl;

        double dr2 = reco::deltaR2(cand1->tunePMuonBestTrack()->eta(),cand1->tunePMuonBestTrack()->phi(),
                                   cand2->tunePMuonBestTrack()->eta(),cand2->tunePMuonBestTrack()->phi());

        if ( dr2 < drThres2 )
          cands.push_back(std::make_pair(cand1,cand2));
      }
    }

    if ( cands.empty() )
      return;

    std::sort(cands.begin(),cands.end(),sortByDR);

    // third muon is the merged muon - should be highPtMuon
    const pat::MuonRef& firstMuon = cands.front().first;
    const pat::MuonRef& secondMuon = cands.front().second;
    pat::MuonRef mergedMuon;

    for (unsigned idx = 0; idx < isolatedHighPtMuons.size(); idx++) {
      const auto& cand = isolatedHighPtMuons.at(idx);

      if ( cand!=firstMuon && cand!=secondMuon ) {
        mergedMuon = cand;
        break;
      }
    }

    if ( !firstMuon.isNull() && !secondMuon.isNull() && !mergedMuon.isNull() ) {
      histo1d_["cutflow_3M"]->Fill( 7.5, aWeight );

      double momcorr1 = 1.;
      double momcorr2 = 1.;
      double momcorrMM = 1.;

      if (isMC_) {
        edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
        iEvent.getByToken(genptcToken_, genptcHandle);

        momcorr1 = mucorrHelper_.rochesterMC(firstMuon,genptcHandle);
        momcorr2 = mucorrHelper_.rochesterMC(secondMuon,genptcHandle);
        momcorrMM = mucorrHelper_.rochesterMC(mergedMuon,genptcHandle);
      } else {
        momcorr1 = mucorrHelper_.rochesterData(firstMuon);
        momcorr2 = mucorrHelper_.rochesterData(secondMuon);
        momcorrMM = mucorrHelper_.rochesterData(mergedMuon);
      }

      const auto lvecFirst = math::PtEtaPhiMLorentzVector(firstMuon->tunePMuonBestTrack()->pt(),
                                                          firstMuon->tunePMuonBestTrack()->eta(),
                                                          firstMuon->tunePMuonBestTrack()->phi(),
                                                          mumass_);
      const auto lvecSecond = math::PtEtaPhiMLorentzVector(secondMuon->tunePMuonBestTrack()->pt(),
                                                           secondMuon->tunePMuonBestTrack()->eta(),
                                                           secondMuon->tunePMuonBestTrack()->phi(),
                                                           mumass_);
      const auto lvecMerged = math::PtEtaPhiMLorentzVector(mergedMuon->tunePMuonBestTrack()->pt(),
                                                           mergedMuon->tunePMuonBestTrack()->eta(),
                                                           mergedMuon->tunePMuonBestTrack()->phi(),
                                                           mumass_);

      auto firstP4 = lvecFirst*momcorr1;
      auto secondP4 = lvecSecond*momcorr2;
      auto mergedP4 = lvecMerged*momcorrMM;

      histo1d_["3M_check_M1M2_invM"]->Fill( (firstP4+secondP4).M() , aWeight );
      histo1d_["3M_check_M1M2_invM_zoomed"]->Fill( (firstP4+secondP4).M() , aWeight );
      histo1d_["3M_check_M1M2MM_invM"]->Fill( (firstP4+secondP4+mergedP4).M() , aWeight );

      histo1d_["3M_check_M1M2_dR"]->Fill( reco::deltaR(firstP4.eta(),firstP4.phi(),secondP4.eta(),secondP4.phi()) , aWeight );
      histo1d_["3M_check_M1MM_dR"]->Fill( reco::deltaR(firstP4.eta(),firstP4.phi(),mergedP4.eta(),mergedP4.phi()) , aWeight );
      histo1d_["3M_check_M2MM_dR"]->Fill( reco::deltaR(secondP4.eta(),secondP4.phi(),mergedP4.eta(),mergedP4.phi()) , aWeight );
      histo1d_["3M_check_M1_iso"]->Fill( firstMuon->trackIso() , aWeight );
      histo1d_["3M_check_M2_iso"]->Fill( secondMuon->trackIso() , aWeight );
      histo1d_["3M_check_MM_iso"]->Fill( mergedMuon->trackIso() , aWeight );
      bool insideIsoVeto = reco::deltaR2(firstMuon->innerTrack()->eta(),firstMuon->innerTrack()->phi(),
                                         secondMuon->innerTrack()->eta(),secondMuon->innerTrack()->phi()) < 0.0001;
      histo1d_["3M_check_M1_iso_musubtract"]->Fill( insideIsoVeto ? firstMuon->trackIso() :
                                                                    firstMuon->trackIso() - secondMuon->innerTrack()->pt() , aWeight );
      histo1d_["3M_check_M2_iso_musubtract"]->Fill( insideIsoVeto ? secondMuon->trackIso() :
                                                                    secondMuon->trackIso() - firstMuon->innerTrack()->pt() , aWeight );

      const auto& aMET = metHandle->at(0);
      const auto tpMET = aMET.p4() - firstP4 - secondP4 - mergedP4
                                   + firstMuon->p4() + secondMuon->p4() + mergedMuon->p4();

      if ( tpMET.pt() > ptThres_ ) {
        histo1d_["cutflow_3M"]->Fill( 8.5, aWeight );

        double dphi = reco::deltaPhi(mergedP4.phi(),tpMET.phi());
        auto p4 = firstP4 + secondP4 + mergedP4 + tpMET;
        double mt = p4.mt();

        bool passDPhi = std::abs(dphi) < drThres_;
        bool passDPhiCR = (std::abs(dphi) < drThresCR_) && !passDPhi;
        bool passRatioPt = false;

        double ratioPt = (firstP4+secondP4+mergedP4).pt()/tpMET.pt();

        if ( ratioThresLo_ < ratioPt && ratioPt < ratioThresHi_ )
          passRatioPt = true;

        const double ffMM = std::max(ffFunc_->Eval(mt),0.);
        const double xvalMET[1] = {mt};
        double ciMM[1];
        fitResult_->GetConfidenceIntervals(1,1,0,xvalMET,ciMM,0.95,false);

        histo1d_["nPV"]->Fill( static_cast<float>(pvHandle->size())+0.5, aWeight );
        histo1d_["3M_ABCD_MET_dphi"]->Fill(dphi, aWeight);
        histo1d_["3M_ABCD_MET_ratioPt"]->Fill(ratioPt, aWeight);

        histo2d_["3M_mt_dphi"]->Fill(mt,dphi,aWeight);
        histo2d_["3M_mt_ratioPt"]->Fill(mt,ratioPt,aWeight);
        histo2d_["3M_dphi_ratioPt"]->Fill(dphi,ratioPt,aWeight);

        if (passDPhi)
          histo1d_["cutflow_3M"]->Fill( 9.5, aWeight );

        if ( passDPhi && passRatioPt ) {
          histo1d_["cutflow_3M"]->Fill( 10.5, aWeight );

          if (mt < 250. || isMC_) { // blinded
            histo1d_["3M_MM_pt"]->Fill( mergedP4.pt(), aWeight );
            histo1d_["3M_MM_eta"]->Fill( mergedP4.eta(), aWeight );
            histo1d_["3M_MM_phi"]->Fill( mergedP4.phi(), aWeight );

            histo1d_["3M_MET_pt"]->Fill( tpMET.pt(), aWeight );
            histo1d_["3M_MET_phi"]->Fill( tpMET.phi(), aWeight );

            histo1d_["3M_M1_pt"]->Fill( firstP4.pt(), aWeight );
            histo1d_["3M_M1_eta"]->Fill( firstP4.eta(), aWeight );
            histo1d_["3M_M1_phi"]->Fill( firstP4.phi(), aWeight );

            histo1d_["3M_M2_pt"]->Fill( secondP4.pt(), aWeight );
            histo1d_["3M_M2_eta"]->Fill( secondP4.eta(), aWeight );
            histo1d_["3M_M2_phi"]->Fill( secondP4.phi(), aWeight );

            histo1d_["3M_MET_dphi"]->Fill( dphi, aWeight );
            histo1d_["3M_mt"]->Fill( mt, aWeight );
          }
        } else if ( passDPhi && !passRatioPt ) {
          histo1d_["3M_antiRpt_MM_pt"]->Fill( mergedP4.pt(), aWeight );
          histo1d_["3M_antiRpt_MM_pt_xFF"]->Fill( mergedP4.pt(), aWeight*ffMM );
          histo1d_["3M_antiRpt_M1_pt"]->Fill( firstP4.pt(), aWeight );
          histo1d_["3M_antiRpt_M2_pt"]->Fill( secondP4.pt(), aWeight );
          histo1d_["3M_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["3M_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["3M_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["3M_antiRpt_mt"]->Fill( mt, aWeight );
          histo1d_["3M_antiRpt_mt_xFF"]->Fill( mt, aWeight*ffMM );
          histo1d_["3M_antiRpt_mt_xFF_up"]->Fill( mt, aWeight*(ffMM+ciMM[0]) );
          histo1d_["3M_antiRpt_mt_xFF_dn"]->Fill( mt, aWeight*(ffMM-std::min(ciMM[0],ffMM)) );
          histo1d_["3M_antiDphi_MET_ratioPt"]->Fill( ratioPt, aWeight );
        } else if ( passDPhiCR && passRatioPt ) {
          histo1d_["3M_CRdphi_MM_pt"]->Fill( mergedP4.pt(), aWeight );
          histo1d_["3M_CRdphi_MM_eta"]->Fill( mergedP4.eta(), aWeight );
          histo1d_["3M_CRdphi_MM_phi"]->Fill( mergedP4.phi(), aWeight );

          histo1d_["3M_CRdphi_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["3M_CRdphi_MET_phi"]->Fill( tpMET.phi(), aWeight );

          histo1d_["3M_CRdphi_M1_pt"]->Fill( firstP4.pt(), aWeight );
          histo1d_["3M_CRdphi_M1_eta"]->Fill( firstP4.eta(), aWeight );
          histo1d_["3M_CRdphi_M1_phi"]->Fill( firstP4.phi(), aWeight );

          histo1d_["3M_CRdphi_M2_pt"]->Fill( secondP4.pt(), aWeight );
          histo1d_["3M_CRdphi_M2_eta"]->Fill( secondP4.eta(), aWeight );
          histo1d_["3M_CRdphi_M2_phi"]->Fill( secondP4.phi(), aWeight );

          histo1d_["3M_CRdphi_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["3M_CRdphi_mt"]->Fill( mt, aWeight );
        } else if ( passDPhiCR && !passRatioPt ) {
          histo1d_["3M_CRdphi_antiRpt_MM_pt"]->Fill( mergedP4.pt(), aWeight );
          histo1d_["3M_CRdphi_antiRpt_MM_pt_xFF"]->Fill( mergedP4.pt(), aWeight*ffMM );
          histo1d_["3M_CRdphi_antiRpt_M1_pt"]->Fill( firstP4.pt(), aWeight );
          histo1d_["3M_CRdphi_antiRpt_M2_pt"]->Fill( secondP4.pt(), aWeight );
          histo1d_["3M_CRdphi_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["3M_CRdphi_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["3M_CRdphi_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["3M_CRdphi_antiRpt_mt"]->Fill( mt, aWeight );
          histo1d_["3M_CRdphi_antiRpt_mt_xFF"]->Fill( mt, aWeight*ffMM );
          histo1d_["3M_CRdphi_antiRpt_mt_xFF_up"]->Fill( mt, aWeight*(ffMM+ciMM[0]) );
          histo1d_["3M_CRdphi_antiRpt_mt_xFF_dn"]->Fill( mt, aWeight*(ffMM-std::min(ciMM[0],ffMM)) );
          histo1d_["3M_CRdphi_antiDphi_MET_ratioPt"]->Fill( ratioPt, aWeight );
        } else if ( !passDPhi && !passDPhiCR && passRatioPt ) {
          histo1d_["3M_antiDphi_MM_pt"]->Fill( mergedP4.pt(), aWeight );
          histo1d_["3M_antiDphi_MM_eta"]->Fill( mergedP4.eta(), aWeight );
          histo1d_["3M_antiDphi_MM_phi"]->Fill( mergedP4.phi(), aWeight );

          histo1d_["3M_antiDphi_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["3M_antiDphi_MET_phi"]->Fill( tpMET.phi(), aWeight );

          histo1d_["3M_antiDphi_M1_pt"]->Fill( firstP4.pt(), aWeight );
          histo1d_["3M_antiDphi_M1_eta"]->Fill( firstP4.eta(), aWeight );
          histo1d_["3M_antiDphi_M1_phi"]->Fill( firstP4.phi(), aWeight );

          histo1d_["3M_antiDphi_M2_pt"]->Fill( secondP4.pt(), aWeight );
          histo1d_["3M_antiDphi_M2_eta"]->Fill( secondP4.eta(), aWeight );
          histo1d_["3M_antiDphi_M2_phi"]->Fill( secondP4.phi(), aWeight );

          histo1d_["3M_antiDphi_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["3M_antiDphi_MET_ratioPt"]->Fill( ratioPt, aWeight );
          histo1d_["3M_antiDphi_mt"]->Fill( mt, aWeight );

          histo1d_["3M_antiDphi_M1M2_invM"]->Fill( (firstP4+secondP4).M() , aWeight );
          histo1d_["3M_antiDphi_M1M2_invM_zoomed"]->Fill( (firstP4+secondP4).M() , aWeight );
          histo1d_["3M_antiDphi_M1M2MM_invM"]->Fill( (firstP4+secondP4+mergedP4).M() , aWeight );

          histo1d_["3M_antiDphi_M1M2_dR"]->Fill( reco::deltaR(firstP4.eta(),firstP4.phi(),secondP4.eta(),secondP4.phi()) , aWeight );
          histo1d_["3M_antiDphi_M1MM_dR"]->Fill( reco::deltaR(firstP4.eta(),firstP4.phi(),mergedP4.eta(),mergedP4.phi()) , aWeight );
          histo1d_["3M_antiDphi_M2MM_dR"]->Fill( reco::deltaR(secondP4.eta(),secondP4.phi(),mergedP4.eta(),mergedP4.phi()) , aWeight );
        } else if ( !passDPhi && !passDPhiCR && !passRatioPt ) {
          histo1d_["3M_antiDphi_antiRpt_MM_pt"]->Fill( mergedP4.pt(), aWeight );
          histo1d_["3M_antiDphi_antiRpt_MM_pt_xFF"]->Fill( mergedP4.pt(), aWeight*ffMM );
          histo1d_["3M_antiDphi_antiRpt_M1_pt"]->Fill( firstP4.pt(), aWeight );
          histo1d_["3M_antiDphi_antiRpt_M2_pt"]->Fill( secondP4.pt(), aWeight );
          histo1d_["3M_antiDphi_antiRpt_MET_pt"]->Fill( tpMET.pt(), aWeight );
          histo1d_["3M_antiDphi_antiRpt_MET_pt_xFF"]->Fill( tpMET.pt(), aWeight*ffMM );
          histo1d_["3M_antiDphi_antiRpt_MET_dphi"]->Fill( dphi, aWeight );
          histo1d_["3M_antiDphi_antiRpt_mt"]->Fill( mt, aWeight );
          histo1d_["3M_antiDphi_antiRpt_mt_xFF"]->Fill( mt, aWeight*ffMM );
          histo1d_["3M_antiDphi_antiRpt_mt_xFF_up"]->Fill( mt, aWeight*(ffMM+ciMM[0]) );
          histo1d_["3M_antiDphi_antiRpt_mt_xFF_dn"]->Fill( mt, aWeight*(ffMM-std::min(ciMM[0],ffMM)) );
          histo1d_["3M_antiDphi_antiDphi_MET_ratioPt"]->Fill( ratioPt, aWeight );
        } // MET dPhi ratioPt CR
      } // MET ptThres
    } // find a pair of collimated muons
  } // nHighPtMuons==3

  return;
}

DEFINE_FWK_MODULE(MergedMuCRanalyzer);
