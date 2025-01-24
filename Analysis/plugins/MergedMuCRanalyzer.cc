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
#include "DataFormats/PatCandidates/interface/Electron.h"
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
#include "TH2D.h"
#include "TTree.h"
#include "TFitResult.h"

#include "correction.h"

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
  const edm::EDGetTokenT<float> prefweight_token;
  const edm::EDGetTokenT<float> prefweightUp_token_;
  const edm::EDGetTokenT<float> prefweightDn_token_;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;

  const edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;

  const edm::FileInPath MMFFpath_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath triggerSFpath_;
  const edm::FileInPath muonIdIsoSFpath_;
  const edm::FileInPath muonBoostIsoSFpath_;
  const edm::FileInPath muonRecoSFpath_;

  std::unique_ptr<correction::CorrectionSet> purwgt_;
  const std::string puname_;

  const std::string year_;

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

  TTree* interestEvtTree_ = nullptr;
  unsigned int runNo_ = 0;
  unsigned int lumiNo_ = 0;
  unsigned long long evtNo_ = 0;
};

MergedMuCRanalyzer::MergedMuCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
prefweight_token(consumes<float>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
prefweightUp_token_(consumes<float>(edm::InputTag("prefiringweight:nonPrefiringProbUp"))),
prefweightDn_token_(consumes<float>(edm::InputTag("prefiringweight:nonPrefiringProbDown"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
MMFFpath_(iConfig.getParameter<edm::FileInPath>("MMFFpath")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
triggerSFpath_(iConfig.getParameter<edm::FileInPath>("triggerSF")),
muonIdIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonIdIsoSFpath")),
muonBoostIsoSFpath_(iConfig.getParameter<edm::FileInPath>("muonBoostIsoSFpath")),
muonRecoSFpath_(iConfig.getParameter<edm::FileInPath>("muonRecoSFpath")),
purwgt_(std::move(correction::CorrectionSet::from_file((iConfig.getParameter<edm::FileInPath>("PUrwgt")).fullPath()))),
puname_(iConfig.getParameter<std::string>("PUname")),
year_(iConfig.getParameter<std::string>("year")),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptMuThres_(iConfig.getParameter<double>("ptMuThres")),
drThres_(iConfig.getParameter<double>("drThres")),
drThresCR_(iConfig.getParameter<double>("drThresCR")),
ratioThresLo_(iConfig.getParameter<double>("ratioThresLo")),
ratioThresHi_(iConfig.getParameter<double>("ratioThresHi")),
mucorrHelper_(rochesterPath_,triggerSFpath_,muonIdIsoSFpath_,muonBoostIsoSFpath_,muonRecoSFpath_) {
  usesResource("TFileService");
}

void MergedMuCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  MMFFfile_ = std::make_unique<TFile>(MMFFpath_.fullPath().c_str(),"READ");
  ffFunc_ = static_cast<TF1*>(MMFFfile_->FindObjectAny("MMFF"));
  fitResult_ = (static_cast<TH1D*>(MMFFfile_->Get("1M_MMFF_numer_rebin")))->Fit(ffFunc_,"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["totWeightedSum_4M"] = fs->make<TH1D>("totWeightedSum_4M","totWeightedSum_4M",1,0.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);

  histo1d_["GEN_pt_1st"] = fs->make<TH1D>("GEN_pt_1st",";p_{T};",1000,0.,2000.);
  histo1d_["GEN_pt_2nd"] = fs->make<TH1D>("GEN_pt_2nd",";p_{T};",1000,0.,2000.);
  histo1d_["GEN_pt_3rd"] = fs->make<TH1D>("GEN_pt_3rd",";p_{T};",1000,0.,2000.);
  histo1d_["GEN_pt_4th"] = fs->make<TH1D>("GEN_pt_4th",";p_{T};",1000,0.,2000.);

  histo1d_["resolved_trig_denom"] = fs->make<TH1D>("resolved_trig_denom",";p_{T};",1000,0.,2000.);
  histo1d_["resolved_trig_numer"] = fs->make<TH1D>("resolved_trig_numer",";p_{T};",1000,0.,2000.);
  histo1d_["merged_trig_denom"] = fs->make<TH1D>("merged_trig_denom",";p_{T};",1000,0.,2000.);
  histo1d_["merged_trig_numer"] = fs->make<TH1D>("merged_trig_numer",";p_{T};",1000,0.,2000.);

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
  histo1d_["3M_mt_idUp"] = fs->make<TH1D>("3M_mt_idUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_idDn"] = fs->make<TH1D>("3M_mt_idDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_isoUp"] = fs->make<TH1D>("3M_mt_isoUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_isoDn"] = fs->make<TH1D>("3M_mt_isoDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_muBoostIsoUp"] = fs->make<TH1D>("3M_mt_muBoostIsoUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_muBoostIsoDn"] = fs->make<TH1D>("3M_mt_muBoostIsoDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_trigUp"] = fs->make<TH1D>("3M_mt_trigUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_trigDn"] = fs->make<TH1D>("3M_mt_trigDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_recoUp"] = fs->make<TH1D>("3M_mt_recoUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_recoDn"] = fs->make<TH1D>("3M_mt_recoDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_JESup"] = fs->make<TH1D>("3M_mt_JESup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_JESdn"] = fs->make<TH1D>("3M_mt_JESdn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_JERup"] = fs->make<TH1D>("3M_mt_JERup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_JERdn"] = fs->make<TH1D>("3M_mt_JERdn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_prefireUp"] = fs->make<TH1D>("3M_mt_prefireUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_prefireDn"] = fs->make<TH1D>("3M_mt_prefireDn","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_PUrwgtUp"] = fs->make<TH1D>("3M_mt_PUrwgtUp","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_mt_PUrwgtDn"] = fs->make<TH1D>("3M_mt_PUrwgtDn","m_{T};m_{T};",500,0.,2500.);

  histo1d_["3M_GEN_ratioPt"] = fs->make<TH1D>("3M_GEN_ratioPt",";R(p_{T}^{reco}/p_{T}^{GEN});",250,0.,5.);
  histo1d_["3M_GEN_ratioMET"] = fs->make<TH1D>("3M_GEN_ratioMET",";R(MET/p_{T}^{GEN});",250,0.,5.);
  histo1d_["3M_GEN_ratioPtSum"] = fs->make<TH1D>("3M_GEN_ratioPtSum",";R(p_{T}^{reco}+MET / #Sigma p_{T}^{GEN});",250,0.,5.);
  histo1d_["3M_GEN_ratioPtResolved1"] = fs->make<TH1D>("3M_GEN_ratioPtResolved1",";R(p_{T}^{reco}/p_{T}^{GEN});",250,0.,5.);
  histo1d_["3M_GEN_ratioPtResolved2"] = fs->make<TH1D>("3M_GEN_ratioPtResolved2",";R(p_{T}^{reco}/p_{T}^{GEN});",250,0.,5.);

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
  histo1d_["3M_antiRpt_MET_ratioPt"] = fs->make<TH1D>("3M_antiRpt_MET_ratioPt","p_{T} ratio;R(p_{T});",200,0.,5.);
  histo1d_["3M_antiRpt_mt_xFF_JESup"] = fs->make<TH1D>("3M_antiRpt_mt_xFF_JESup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiRpt_mt_xFF_JESdn"] = fs->make<TH1D>("3M_antiRpt_mt_xFF_JESdn","m_{T};m_{T};",500,0.,2500.);

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
  histo1d_["3M_CRdphi_mt_JESup"] = fs->make<TH1D>("3M_CRdphi_mt_JESup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_CRdphi_mt_JESdn"] = fs->make<TH1D>("3M_CRdphi_mt_JESdn","m_{T};m_{T};",500,0.,2500.);

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
  histo1d_["3M_CRdphi_antiRpt_mt_xFF_JESup"] = fs->make<TH1D>("3M_CRdphi_antiRpt_mt_xFF_JESup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_CRdphi_antiRpt_mt_xFF_JESdn"] = fs->make<TH1D>("3M_CRdphi_antiRpt_mt_xFF_JESdn","m_{T};m_{T};",500,0.,2500.);

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
  histo1d_["3M_antiDphi_mt_JESup"] = fs->make<TH1D>("3M_antiDphi_mt_JESup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_mt_JESdn"] = fs->make<TH1D>("3M_antiDphi_mt_JESdn","m_{T};m_{T};",500,0.,2500.);

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
  histo1d_["3M_antiDphi_antiRpt_mt_xFF_JESup"] = fs->make<TH1D>("3M_antiDphi_antiRpt_mt_xFF_JESup","m_{T};m_{T};",500,0.,2500.);
  histo1d_["3M_antiDphi_antiRpt_mt_xFF_JESdn"] = fs->make<TH1D>("3M_antiDphi_antiRpt_mt_xFF_JESdn","m_{T};m_{T};",500,0.,2500.);

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

  interestEvtTree_ = fs->make<TTree>("evtTree","evtTree");
  interestEvtTree_->Branch("runNo",&runNo_,"runNo/i");
  interestEvtTree_->Branch("lumiNo",&lumiNo_,"lumiNo/i");
  interestEvtTree_->Branch("evtNo",&evtNo_,"evtNo/l");
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

  std::vector<reco::GenParticleRef> promptMuons;

  double aWeight = 1.;
  double purwgtNo = 1.;
  double purwgtUp = 1.;
  double purwgtDn = 1.;
  float prefireNo = 1.;
  float prefireUp = 1.;
  float prefireDn = 1.;

  if (isMC_) {
    edm::Handle<float> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    prefireNo = *theprefweight;

    edm::Handle<float> theprefweightUp;
    iEvent.getByToken(prefweightUp_token_, theprefweightUp);
    prefireUp = *theprefweightUp;

    edm::Handle<float> theprefweightDn;
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

    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(genptcToken_, genptcHandle);

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( ( std::abs(genptc->pdgId())==13 ) && genptc->fromHardProcessFinalState() )
        promptMuons.push_back(genptc.castTo<reco::GenParticleRef>());
    }

    auto sortByPt = [](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
      return a->pt() > b->pt();
    };

    if (promptMuons.size()==4) {
      histo1d_["totWeightedSum_4M"]->Fill(0.5,aWeight);

      std::sort(promptMuons.begin(),promptMuons.end(),sortByPt);

      histo1d_["GEN_pt_1st"]->Fill(promptMuons.at(0)->pt(),aWeight);
      histo1d_["GEN_pt_2nd"]->Fill(promptMuons.at(1)->pt(),aWeight);
      histo1d_["GEN_pt_3rd"]->Fill(promptMuons.at(2)->pt(),aWeight);
      histo1d_["GEN_pt_4th"]->Fill(promptMuons.at(3)->pt(),aWeight);
    }
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
  std::vector<pat::MuonRef> boostedMuons;

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

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 ) {
      isolatedHighPtMuons.push_back(firstMuon);

      if (!isoCands.empty())
        boostedMuons.push_back(firstMuon);
    }
  }

  for (unsigned idx = 0; idx < highPtTrackerMuons.size(); idx++) {
    const auto& firstMuon = highPtTrackerMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 ) {
      isolatedHighPtTrackerMuons.push_back(firstMuon);

      if (!isoCands.empty())
        boostedMuons.push_back(firstMuon);
    }
  }

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  if ( isolatedHighPtMuons.size() < 2 )
    return;

  // concatenate(-ish) highPtMuons & highPtTrackerMuons - order is important!
  std::vector<pat::MuonRef> allHighPtMuons(isolatedHighPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), isolatedHighPtTrackerMuons.begin(), isolatedHighPtTrackerMuons.end() );
  std::vector<pat::MuonRef> sortedHighPtMuons(allHighPtMuons);
  std::sort(sortedHighPtMuons.begin(),sortedHighPtMuons.end(),sortByTuneP);

  unsigned int nHighPtMuons = allHighPtMuons.size();

  if ( sortedHighPtMuons.front()->tunePMuonBestTrack()->pt() < 52. )
    return;

  histo1d_["cutflow_4M"]->Fill( 3.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 3.5, aWeight );

  bool trigMatched = false;
  std::vector<double> trigSyst;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    const auto& leadMu = sortedHighPtMuons.front();
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

  std::vector<pat::MuonRef> nonHighPtMuons;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( !aMuon->isTrackerMuon() || aMuon->track().isNull() )
      continue;

    if ( aMuon->tunePMuonBestTrack()->pt() < ptMuThres_ || std::abs(aMuon->eta()) > 2.4 )
      continue;

    const auto castMu = aMuon.castTo<pat::MuonRef>();

    if ( std::find( allHighPtMuons.begin(), allHighPtMuons.end(), castMu ) != allHighPtMuons.end() )
      continue;

    bool hits = aMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
                && aMuon->innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

    bool ip = std::abs(aMuon->innerTrack()->dxy(primaryVertex.position())) < 0.2
              && std::abs(aMuon->innerTrack()->dz(primaryVertex.position())) < 0.5;

    if ( aMuon->numberOfMatchedStations() >= 1 && hits && ip )
      nonHighPtMuons.push_back(castMu);
  }

  if (!nonHighPtMuons.empty())
    return;

  histo1d_["cutflow_4M"]->Fill( 4.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 4.5, aWeight );

  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

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
      auto apair = mucorrHelper_.boostIsoSFupdn(aMu,primaryVertex);
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
        isoSyst.push_back( mucorrHelper_.looseIsoSFsyst(aMu)/mucorrHelper_.looseIsoSF(aMu) );
    }
  }

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

        momcorr1 = mucorrHelper_.nominalMC(firstMuon,genptcHandle);
        momcorr2 = mucorrHelper_.nominalMC(secondMuon,genptcHandle);
        momcorrMM = mucorrHelper_.nominalMC(mergedMuon,genptcHandle);
      } else {
        momcorr1 = mucorrHelper_.nominalData(firstMuon);
        momcorr2 = mucorrHelper_.nominalData(secondMuon);
        momcorrMM = mucorrHelper_.nominalData(mergedMuon);
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
      // const auto corrMET = METXYCorr_Met_MetPhi(aMET.pt(),aMET.phi(),
      //                                           iEvent.id().run(),TString(year_),isMC_,
      //                                           pvHandle->size(), true, false);
      const auto metp4 = aMET.p4(); // math::PtEtaPhiMLorentzVector(corrMET.first,aMET.p4().eta(),corrMET.second,aMET.p4().M());
      const auto tpMET = metp4 - firstP4 - secondP4 - mergedP4
                               + firstMuon->p4() + secondMuon->p4() + mergedMuon->p4();
      const auto jesUp = math::PtEtaPhiMLorentzVector(aMET.shiftedPt(pat::MET::JetEnUp),0.,aMET.shiftedPhi(pat::MET::JetEnUp),0.);
      const auto jesDn = math::PtEtaPhiMLorentzVector(aMET.shiftedPt(pat::MET::JetEnDown),0.,aMET.shiftedPhi(pat::MET::JetEnDown),0.);
      const auto deltaJESup = jesUp - aMET.p4();
      const auto deltaJESdn = jesDn - aMET.p4();
      const auto jerUp = math::PtEtaPhiMLorentzVector(aMET.shiftedPt(pat::MET::JetResUp),0.,aMET.shiftedPhi(pat::MET::JetResUp),0.);
      const auto jerDn = math::PtEtaPhiMLorentzVector(aMET.shiftedPt(pat::MET::JetResDown),0.,aMET.shiftedPhi(pat::MET::JetResDown),0.);
      const auto deltaJERup = jerUp - aMET.p4();
      const auto deltaJERdn = jerDn - aMET.p4();

      if ( tpMET.pt() > ptThres_ ) {
        histo1d_["cutflow_3M"]->Fill( 8.5, aWeight );

        double dphi = reco::deltaPhi(mergedP4.phi(),tpMET.phi());
        auto p4 = firstP4 + secondP4 + mergedP4 + tpMET;
        double mt = std::min(p4.mt(),2499.9);
        const double mtJESup = std::min((p4+deltaJESup).mt(),2499.9);
        const double mtJESdn = std::min((p4+deltaJESdn).mt(),2499.9);
        const double mtJERup = std::min((p4+deltaJERup).mt(),2499.9);
        const double mtJERdn = std::min((p4+deltaJERdn).mt(),2499.9);

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

          if ( true /*mt < 500. || isMC_ (unblinded)*/ ) { // blinded
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
            histo1d_["3M_mt_idUp"]->Fill( mt, aWeight*systSFratio(idSyst).first );
            histo1d_["3M_mt_idDn"]->Fill( mt, aWeight*systSFratio(idSyst).second );
            histo1d_["3M_mt_isoUp"]->Fill( mt, aWeight*systSFratio(isoSyst).first );
            histo1d_["3M_mt_isoDn"]->Fill( mt, aWeight*systSFratio(isoSyst).second );
            histo1d_["3M_mt_muBoostIsoUp"]->Fill( mt, boostIsoSFupdn(aWeight).first );
            histo1d_["3M_mt_muBoostIsoDn"]->Fill( mt, boostIsoSFupdn(aWeight).second );
            histo1d_["3M_mt_trigUp"]->Fill( mt, aWeight*systSFratio(trigSyst).first );
            histo1d_["3M_mt_trigDn"]->Fill( mt, aWeight*systSFratio(trigSyst).second );
            histo1d_["3M_mt_recoUp"]->Fill( mt, aWeight*systSFratio(recoSyst).first );
            histo1d_["3M_mt_recoDn"]->Fill( mt, aWeight*systSFratio(recoSyst).second );
            histo1d_["3M_mt_JESup"]->Fill( mtJESup, aWeight );
            histo1d_["3M_mt_JESdn"]->Fill( mtJESdn, aWeight );
            histo1d_["3M_mt_JERup"]->Fill( mtJERup, aWeight );
            histo1d_["3M_mt_JERdn"]->Fill( mtJERdn, aWeight );
            histo1d_["3M_mt_PUrwgtUp"]->Fill( mt, aWeight*purwgtUp/purwgtNo );
            histo1d_["3M_mt_PUrwgtDn"]->Fill( mt, aWeight*purwgtDn/purwgtNo );
            histo1d_["3M_mt_prefireUp"]->Fill( mt, aWeight*prefireUp/prefireNo );
            histo1d_["3M_mt_prefireDn"]->Fill( mt, aWeight*prefireDn/prefireNo );

            if (mt > 800.) {
              runNo_ = iEvent.id().run();
              lumiNo_ = iEvent.id().luminosityBlock();
              evtNo_ = iEvent.id().event();
              interestEvtTree_->Fill();
            }

            if (isMC_) {
              float ptGen = 1e-3, ptsumGen = 1e-3, diff = std::numeric_limits<float>::max();
              float ptGen1 = 1e-3, diff1 = std::numeric_limits<float>::max();
              float ptGen2 = 1e-3, diff2 = std::numeric_limits<float>::max();
              reco::GenParticleRef mergedMuonRef, cleanedMuRef;
              std::vector<reco::GenParticleRef> mergedPair;

              for (const auto muGen : promptMuons) {
                if (reco::deltaR(mergedP4.eta(),mergedP4.phi(),muGen->eta(),muGen->phi()) < 0.3) {
                  ptsumGen += muGen->pt();
                  mergedPair.push_back(muGen);

                  if ( std::abs(muGen->pt() - mergedP4.pt()) < diff ) {
                    ptGen = muGen->pt();
                    diff = std::abs(muGen->pt() - mergedP4.pt());
                    mergedMuonRef = muGen;
                  }
                }

                if (reco::deltaR(firstP4.eta(),firstP4.phi(),muGen->eta(),muGen->phi()) < 0.3) {
                  if ( std::abs(muGen->pt() - firstP4.pt()) < diff1 ) {
                    ptGen1 = muGen->pt();
                    diff1 = std::abs(muGen->pt() - firstP4.pt());
                  }
                }

                if (reco::deltaR(secondP4.eta(),secondP4.phi(),muGen->eta(),muGen->phi()) < 0.3) {
                  if ( std::abs(muGen->pt() - secondP4.pt()) < diff2 ) {
                    ptGen2 = muGen->pt();
                    diff2 = std::abs(muGen->pt() - secondP4.pt());
                  }
                }
              }

              for (const auto& muGen : mergedPair) {
                if (muGen!=mergedMuonRef) {
                  cleanedMuRef = muGen;
                  break;
                }
              }

              if (ptGen > 200.)
                histo1d_["3M_GEN_ratioPt"]->Fill(mergedP4.pt()/ptGen, aWeight);

              if (!cleanedMuRef.isNull() && cleanedMuRef->pt() > 200.)
                histo1d_["3M_GEN_ratioMET"]->Fill(tpMET.pt()/cleanedMuRef->pt(), aWeight);

              if (ptGen1 > 200.)
                histo1d_["3M_GEN_ratioPtResolved1"]->Fill(firstP4.pt()/ptGen1, aWeight);

              if (ptGen2 > 200.)
                histo1d_["3M_GEN_ratioPtResolved2"]->Fill(secondP4.pt()/ptGen2, aWeight);

              histo1d_["3M_GEN_ratioPtSum"]->Fill((mergedP4.pt()+tpMET.pt())/ptsumGen, aWeight);
            }

            auto checkTrig = [&trigObjs] (const math::PtEtaPhiMLorentzVector& aLvec) -> bool {
              for (const auto& trigObj : trigObjs) {
                if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),aLvec.eta(),aLvec.phi()) < 0.01 )
                  return true;
              }

              return false;
            };

            histo1d_["resolved_trig_denom"]->Fill(firstP4.pt(), aWeight);
            histo1d_["resolved_trig_denom"]->Fill(secondP4.pt(), aWeight);
            histo1d_["merged_trig_denom"]->Fill(mergedP4.pt(), aWeight);

            if (checkTrig(firstP4))
              histo1d_["resolved_trig_numer"]->Fill(firstP4.pt(), aWeight);

            if (checkTrig(secondP4))
              histo1d_["resolved_trig_numer"]->Fill(secondP4.pt(), aWeight);

            if (checkTrig(mergedP4))
              histo1d_["merged_trig_numer"]->Fill(mergedP4.pt(), aWeight);
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
          histo1d_["3M_antiRpt_mt_xFF_JESup"]->Fill( mtJESup, aWeight*ffMM );
          histo1d_["3M_antiRpt_mt_xFF_JESdn"]->Fill( mtJESdn, aWeight*ffMM );
          histo1d_["3M_antiRpt_MET_ratioPt"]->Fill( ratioPt, aWeight );
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
          histo1d_["3M_CRdphi_mt_JESup"]->Fill( mtJESup, aWeight );
          histo1d_["3M_CRdphi_mt_JESdn"]->Fill( mtJESdn, aWeight );
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
          histo1d_["3M_CRdphi_antiRpt_mt_xFF_JESup"]->Fill( mtJESup, aWeight*ffMM );
          histo1d_["3M_CRdphi_antiRpt_mt_xFF_JESdn"]->Fill( mtJESdn, aWeight*ffMM );
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
          histo1d_["3M_antiDphi_mt_JESup"]->Fill( mtJESup, aWeight );
          histo1d_["3M_antiDphi_mt_JESdn"]->Fill( mtJESdn, aWeight );

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
          histo1d_["3M_antiDphi_antiRpt_mt_xFF_JESup"]->Fill( mtJESup, aWeight*ffMM );
          histo1d_["3M_antiDphi_antiRpt_mt_xFF_JESdn"]->Fill( mtJESdn, aWeight*ffMM );
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
