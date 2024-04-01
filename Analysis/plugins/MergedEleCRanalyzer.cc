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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"
#include "ZprimeTo4l/Analysis/interface/ElectronSystematicsHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TF2.h"
#include "TFile.h"
#include "TFitResult.h"

class MergedEleCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleCRanalyzer(const edm::ParameterSet&);
  virtual ~MergedEleCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  const bool isMC_;
  const bool selectHT_;
  const bool selectZpT_;
  const double maxHT_;
  const double maxVpT_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;

  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dEtaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dPhiInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5EnergyToken_;

  const edm::FileInPath FFpath_;
  const edm::FileInPath purwgtPath_;
  const std::vector<std::string> trigFilters_;
  const std::vector<std::string> trigUnseededFilters_;
  const std::vector<std::string> METfilterList_;
  const double etThres1_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::unique_ptr<TFile> FFfile_;
  TF1* ssboth_;
  TF2* osboth_;
  TFitResultPtr ssFit_;
  TFitResultPtr osFit_;

  ElectronSystematicsHelper systHelperEle_;

  TTree* SRtree_ = nullptr;
  float SRinvM_ = -1.;
  float SRwgt_ = 0.;
  int SRisOS_ = -1;
  float SRpt1_ = -1.;
  float SReta1_ = std::numeric_limits<float>::max();
  float SRphi1_ = std::numeric_limits<float>::max();
  float SRpt2_ = -1.;
  float SReta2_ = std::numeric_limits<float>::max();
  float SRphi2_ = std::numeric_limits<float>::max();

  TTree* mixedCRtree_ = nullptr;
  float mixedCRinvM_ = -1.;
  float mixedCRwgt_ = 0.;
  int mixedCRisOS_ = -1;
  float mixedCRpt1_ = -1.;
  float mixedCReta1_ = std::numeric_limits<float>::max();
  float mixedCRphi1_ = std::numeric_limits<float>::max();
  float mixedCRpt2_ = -1.;
  float mixedCReta2_ = std::numeric_limits<float>::max();
  float mixedCRphi2_ = std::numeric_limits<float>::max();

  TTree* antiCRtree_ = nullptr;
  float antiCRinvM_ = -1.;
  float antiCRwgt_ = 0.;
  int antiCRisOS_ = -1;
  float antiCRpt1_ = -1.;
  float antiCReta1_ = std::numeric_limits<float>::max();
  float antiCRphi1_ = std::numeric_limits<float>::max();
  float antiCRpt2_ = -1.;
  float antiCReta2_ = std::numeric_limits<float>::max();
  float antiCRphi2_ = std::numeric_limits<float>::max();
};

MergedEleCRanalyzer::MergedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
selectHT_(iConfig.getParameter<bool>("selectHT")),
selectZpT_(iConfig.getParameter<bool>("selectZpT")),
maxHT_(iConfig.getParameter<double>("maxHT")),
maxVpT_(iConfig.getParameter<double>("maxVpT")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
union5x5dEtaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dEtaIn"))),
union5x5dPhiInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dPhiIn"))),
union5x5EnergyToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5Energy"))),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
trigFilters_(iConfig.getParameter<std::vector<std::string>>("trigFilters")),
trigUnseededFilters_(iConfig.getParameter<std::vector<std::string>>("trigUnseededFilters")),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
etThres1_(iConfig.getParameter<double>("etThres1")),
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

  systHelperEle_.SetMergedEleSF(iConfig.getParameter<double>("mergedEleSFmuHasTrk"),
                                iConfig.getParameter<double>("mergedEleSFmuNoTrk"));
  systHelperEle_.SetMergedEleSFcl95(iConfig.getParameter<double>("mergedEleSFcl95HasTrk"),
                                    iConfig.getParameter<double>("mergedEleSFcl95NoTrk"));
  systHelperEle_.SetMergedEleSFupper(iConfig.getParameter<double>("mergedEleSFupperHasTrk"),
                                     iConfig.getParameter<double>("mergedEleSFupperNoTrk"));
  systHelperEle_.SetMergedElePol1(iConfig.getParameter<std::string>("mergedElePolHasTrkStr"),
                                  iConfig.getParameter<std::string>("mergedElePolNoTrkStr"));

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

  systHelperEle_.SetTrigMergedSFandErr(iConfig.getParameter<double>("mergedTrigSF"),
                                       iConfig.getParameter<double>("mergedTrigSFerr"));

  systHelperEle_.SetAbcdScaleSmear(iConfig.getParameter<double>("abcdScaleAbove"),
                                   iConfig.getParameter<double>("abcdScaleBelow"),
                                   iConfig.getParameter<double>("abcdSmear"));
}

void MergedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  ssboth_ = static_cast<TF1*>(FFfile_->FindObjectAny("ssboth"));
  osboth_ = static_cast<TF2*>(FFfile_->FindObjectAny("osboth"));
  ssFit_ = (static_cast<TH1D*>(FFfile_->Get("2E_Et_SSCR_EB_mixedME")))->Fit(ssboth_,"RS");
  osFit_ = (static_cast<TH2D*>(FFfile_->Get("2E_Et_eta_OSCR_EB_mixedME")))->Fit(osboth_,"RS");

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["lheHT"] = fs->make<TH1D>("lheHT","lheHT",600,0.,3000.);
  histo1d_["lheHT_cut"] = fs->make<TH1D>("lheHT_cut","lheHT_cut",600,0.,3000.);
  histo1d_["lheVpT"] = fs->make<TH1D>("lheVpT","lheVpT",500,0.,1000.);
  histo1d_["lheVpT_cut"] = fs->make<TH1D>("lheVpT_cut","lheVpT",500,0.,1000.);
  histo1d_["N-1_modifiedHEEP_EB"] = fs->make<TH1D>("N-1_modifiedHEEP_EB","N-1 efficiency (modified HEEP EB)",16,-1.,15.);
  histo1d_["N-1_modifiedHEEP_EE"] = fs->make<TH1D>("N-1_modifiedHEEP_EE","N-1 efficiency (modified HEEP EE)",16,-1.,15.);
  histo1d_["N-1_HEEP_EB"] = fs->make<TH1D>("N-1_HEEP_EB","N-1 efficiency (HEEP EB)",16,-1.,15.);
  histo1d_["N-1_HEEP_EE"] = fs->make<TH1D>("N-1_HEEP_EE","N-1 efficiency (HEEP EE)",16,-1.,15.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",100,0.,1.);
  histo1d_["mva_bkgEt2EB"] = fs->make<TH1D>("mva_bkgEt2EB","MVA score",100,0.,1.);
  histo1d_["cutflow_4E"] = fs->make<TH1D>("cutflow_4E","cutflow (4E)",10,0.,10.);
  histo1d_["cutflow_3E"] = fs->make<TH1D>("cutflow_3E","cutflow (3E)",10,0.,10.);
  histo1d_["cutflow_2E"] = fs->make<TH1D>("cutflow_2E","cutflow (2E)",10,0.,10.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);
  histo1d_["PUsummary"] = fs->make<TH1D>("PUsummary","PUsummary",99,0.,99.);

  histo1d_["invM_ll_CRME_EB"] = fs->make<TH1D>("invM_ll_CRME_EB","ME M(ll);[GeV];",1000,0.,10.);

  histo1d_["2E_CRME1_Et"] = fs->make<TH1D>("2E_CRME1_Et","CR ME1 E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME1_eta"] = fs->make<TH1D>("2E_CRME1_eta","CR ME1 #eta",50,-2.5,2.5);
  histo1d_["2E_CRME1_phi"] = fs->make<TH1D>("2E_CRME1_phi","CR ME1 #phi",64,-3.2,3.2);

  histo1d_["2E_CRME1_noGsf_Et"] = fs->make<TH1D>("2E_CRME1_noGsf_Et","CR ME1 (w/o GSF) E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME1_noGsf_eta"] = fs->make<TH1D>("2E_CRME1_noGsf_eta","CR ME1 (w/o GSF) #eta",50,-2.5,2.5);
  histo1d_["2E_CRME1_noGsf_phi"] = fs->make<TH1D>("2E_CRME1_noGsf_phi","CR ME1 (w/o GSF) #phi",64,-3.2,3.2);

  histo1d_["2E_CRME2_Et"] = fs->make<TH1D>("2E_CRME2_Et","CR ME2 E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME2_eta"] = fs->make<TH1D>("2E_CRME2_eta","CR ME2 #eta",50,-2.5,2.5);
  histo1d_["2E_CRME2_phi"] = fs->make<TH1D>("2E_CRME2_phi","CR ME2 #phi",64,-3.2,3.2);

  histo1d_["2E_CRME2_noGsf_Et"] = fs->make<TH1D>("2E_CRME2_noGsf_Et","CR ME2 (w/o GSF) E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME2_noGsf_eta"] = fs->make<TH1D>("2E_CRME2_noGsf_eta","CR ME2 (w/o GSF) #eta",50,-2.5,2.5);
  histo1d_["2E_CRME2_noGsf_phi"] = fs->make<TH1D>("2E_CRME2_noGsf_phi","CR ME2 (w/o GSF) #phi",64,-3.2,3.2);

  histo1d_["2E_CRME_ll_invM"] = fs->make<TH1D>("2E_CRME_ll_invM","CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_ll_pt"] = fs->make<TH1D>("2E_CRME_ll_pt","CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_ll_dr"] = fs->make<TH1D>("2E_CRME_ll_dr","CR dR(ll)",64,0.,6.4);

  histo1d_["2E_CRME_SSll_invM"] = fs->make<TH1D>("2E_CRME_SSll_invM","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_mergedEleScale"] = fs->make<TH1D>("2E_CRME_SSll_invM_mergedEleScale","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_mergedEleSmear"] = fs->make<TH1D>("2E_CRME_SSll_invM_mergedEleSmear","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_heepIdUp"] = fs->make<TH1D>("2E_CRME_SSll_invM_heepIdUp","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_heepIdDn"] = fs->make<TH1D>("2E_CRME_SSll_invM_heepIdDn","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_mergedEleIdUp"] = fs->make<TH1D>("2E_CRME_SSll_invM_mergedEleIdUp","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_mergedEleIdDn"] = fs->make<TH1D>("2E_CRME_SSll_invM_mergedEleIdDn","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_elRecoUp"] = fs->make<TH1D>("2E_CRME_SSll_invM_elRecoUp","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_elRecoDn"] = fs->make<TH1D>("2E_CRME_SSll_invM_elRecoDn","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_elTrigUp"] = fs->make<TH1D>("2E_CRME_SSll_invM_elTrigUp","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_invM_elTrigDn"] = fs->make<TH1D>("2E_CRME_SSll_invM_elTrigDn","SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_SSll_pt"] = fs->make<TH1D>("2E_CRME_SSll_pt","SS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_SSll_dr"] = fs->make<TH1D>("2E_CRME_SSll_dr","SS CR dR(ll)",64,0.,6.4);
  histo1d_["2E_CRME_SSll_Et"] = fs->make<TH1D>("2E_CRME_SSll_Et","SS CR ME E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME_SSll_eta"] = fs->make<TH1D>("2E_CRME_SSll_eta","SS CR ME #eta",50,-2.5,2.5);

  histo1d_["2E_CRME_OSll_invM"] = fs->make<TH1D>("2E_CRME_OSll_invM","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_mergedEleScale"] = fs->make<TH1D>("2E_CRME_OSll_invM_mergedEleScale","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_mergedEleSmear"] = fs->make<TH1D>("2E_CRME_OSll_invM_mergedEleSmear","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_heepIdUp"] = fs->make<TH1D>("2E_CRME_OSll_invM_heepIdUp","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_heepIdDn"] = fs->make<TH1D>("2E_CRME_OSll_invM_heepIdDn","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_elRecoUp"] = fs->make<TH1D>("2E_CRME_OSll_invM_elRecoUp","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_elRecoDn"] = fs->make<TH1D>("2E_CRME_OSll_invM_elRecoDn","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_elTrigUp"] = fs->make<TH1D>("2E_CRME_OSll_invM_elTrigUp","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_elTrigDn"] = fs->make<TH1D>("2E_CRME_OSll_invM_elTrigDn","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_mergedEleIdUp"] = fs->make<TH1D>("2E_CRME_OSll_invM_mergedEleIdUp","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_invM_mergedEleIdDn"] = fs->make<TH1D>("2E_CRME_OSll_invM_mergedEleIdDn","OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_CRME_OSll_pt"] = fs->make<TH1D>("2E_CRME_OSll_pt","OS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_OSll_dr"] = fs->make<TH1D>("2E_CRME_OSll_dr","OS CR dR(ll)",64,0.,6.4);
  histo1d_["2E_CRME_OSll_Et"] = fs->make<TH1D>("2E_CRME_OSll_Et","OS CR ME E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME_OSll_eta"] = fs->make<TH1D>("2E_CRME_OSll_eta","OS CR ME #eta",50,-2.5,2.5);

  histo1d_["2E_mixed_ME_Et"] = fs->make<TH1D>("2E_mixed_ME_Et","mixed CR ME E_{T};GeV;",100,0.,500.);
  histo1d_["2E_mixed_ME_eta"] = fs->make<TH1D>("2E_mixed_ME_eta","mixed CR ME #eta",100,-2.5,2.5);
  histo1d_["2E_mixed_ME_phi"] = fs->make<TH1D>("2E_mixed_ME_phi","mixed CR ME #phi",64,-3.2,3.2);

  histo1d_["2E_mixed_ME_noGsf_Et"] = fs->make<TH1D>("2E_mixed_ME_noGsf_Et","mixed CR ME (w/o GSF) E_{T};GeV;",100,0.,500.);
  histo1d_["2E_mixed_ME_noGsf_eta"] = fs->make<TH1D>("2E_mixed_ME_noGsf_eta","mixed CR ME (w/o GSF) #eta",100,-2.5,2.5);
  histo1d_["2E_mixed_ME_noGsf_phi"] = fs->make<TH1D>("2E_mixed_ME_noGsf_phi","mixed CR ME (w/o GSF) #phi",64,-3.2,3.2);

  histo1d_["2E_mixed_antiME_Et"] = fs->make<TH1D>("2E_mixed_antiME_Et","mixed CR anti ME E_{T};GeV;",100,0.,500.);
  histo1d_["2E_mixed_antiME_eta"] = fs->make<TH1D>("2E_mixed_antiME_eta","mixed CR anti ME #eta",100,-2.5,2.5);
  histo1d_["2E_mixed_antiME_phi"] = fs->make<TH1D>("2E_mixed_antiME_phi","mixed CR anti ME #phi",64,-3.2,3.2);

  histo1d_["2E_mixed_antiME_noGsf_Et"] = fs->make<TH1D>("2E_mixed_antiME_noGsf_Et","mixed CR anti ME (w/o GSF) E_{T};GeV;",100,0.,500.);
  histo1d_["2E_mixed_antiME_noGsf_eta"] = fs->make<TH1D>("2E_mixed_antiME_noGsf_eta","mixed CR anti ME (w/o GSF) #eta",100,-2.5,2.5);
  histo1d_["2E_mixed_antiME_noGsf_phi"] = fs->make<TH1D>("2E_mixed_antiME_noGsf_phi","mixed CR anti ME (w/o GSF) #phi",64,-3.2,3.2);

  histo1d_["2E_mixedME_ll_invM"] = fs->make<TH1D>("2E_mixedME_ll_invM","mixed CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_ll_pt"] = fs->make<TH1D>("2E_mixedME_ll_pt","mixed CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_ll_dr"] = fs->make<TH1D>("2E_mixedME_ll_dr","mixed CR dR(ll)",64,0.,6.4);

  histo1d_["2E_mixedME_SSll_invM"] = fs->make<TH1D>("2E_mixedME_SSll_invM","mixed SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_xFF"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF","mixed SS CR M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_xFF_preCorr"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF_preCorr","mixed SS CR M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_xFF_up"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF_up","mixed SS CR M(ll) x Fake factor (up);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_xFF_dn"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF_dn","mixed SS CR M(ll) x Fake factor (down);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_heepIdUp"] = fs->make<TH1D>("2E_mixedME_SSll_invM_heepIdUp","mixed SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_heepIdDn"] = fs->make<TH1D>("2E_mixedME_SSll_invM_heepIdDn","mixed SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_mergedEleIdUp"] = fs->make<TH1D>("2E_mixedME_SSll_invM_mergedEleIdUp","mixed SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_invM_mergedEleIdDn"] = fs->make<TH1D>("2E_mixedME_SSll_invM_mergedEleIdDn","mixed SS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_SSll_pt"] = fs->make<TH1D>("2E_mixedME_SSll_pt","mixed SS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_SSll_dr"] = fs->make<TH1D>("2E_mixedME_SSll_dr","mixed SS CR dR(ll)",64,0.,6.4);
  histo1d_["2E_mixedME_SSll_Et"] = fs->make<TH1D>("2E_mixedME_SSll_Et","mixed SS CR ME E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_mixedME_SSll_Et_noCorr"] = fs->make<TH1D>("2E_mixedME_SSll_Et_noCorr","mixed SS CR ME E_{T} (no S&S corr);GeV;",200,0.,1000.);
  histo1d_["2E_mixedME_SSll_eta"] = fs->make<TH1D>("2E_mixedME_SSll_eta","mixed SS CR ME #eta",100,-2.5,2.5);
  histo1d_["2E_mixedAntiME_SSll_Et_xFF"] = fs->make<TH1D>("2E_mixedAntiME_SSll_Et_xFF","mixed SS CR anti ME E_{T} x Fake factor;GeV;",100,0.,500.);
  histo1d_["2E_mixedAntiME_SSll_eta_xFF"] = fs->make<TH1D>("2E_mixedAntiME_SSll_eta_xFF","mixed SS CR anti ME #eta x Fake factor",100,-2.5,2.5);
  histo1d_["2E_mixedAntiME_SSll_Et_xFF_up"] = fs->make<TH1D>("2E_mixedAntiME_SSll_Et_xFF_up","mixed SS CR anti ME E_{T} x Fake factor (up);GeV;",100,0.,500.);
  histo1d_["2E_mixedAntiME_SSll_eta_xFF_up"] = fs->make<TH1D>("2E_mixedAntiME_SSll_eta_xFF_up","mixed SS CR anti ME #eta x Fake factor (up)",100,-2.5,2.5);
  histo1d_["2E_mixedAntiME_SSll_Et_xFF_dn"] = fs->make<TH1D>("2E_mixedAntiME_SSll_Et_xFF_dn","mixed SS CR anti ME E_{T} x Fake factor (down);GeV;",100,0.,500.);
  histo1d_["2E_mixedAntiME_SSll_eta_xFF_dn"] = fs->make<TH1D>("2E_mixedAntiME_SSll_eta_xFF_dn","mixed SS CR anti ME #eta x Fake factor (down)",100,-2.5,2.5);

  histo1d_["2E_mixedME_OSll_invM"] = fs->make<TH1D>("2E_mixedME_OSll_invM","mixed OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_xFF"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF","mixed OS CR M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_xFF_preCorr"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF_preCorr","mixed OS CR M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_xFF_up"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF_up","mixed OS CR M(ll) x Fake factor (up);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_xFF_dn"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF_dn","mixed OS CR M(ll) x Fake factor (down);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_heepIdUp"] = fs->make<TH1D>("2E_mixedME_OSll_invM_heepIdUp","mixed OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_heepIdDn"] = fs->make<TH1D>("2E_mixedME_OSll_invM_heepIdDn","mixed OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_mergedEleIdUp"] = fs->make<TH1D>("2E_mixedME_OSll_invM_mergedEleIdUp","mixed OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_invM_mergedEleIdDn"] = fs->make<TH1D>("2E_mixedME_OSll_invM_mergedEleIdDn","mixed OS CR M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_mixedME_OSll_pt"] = fs->make<TH1D>("2E_mixedME_OSll_pt","mixed OS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_OSll_dr"] = fs->make<TH1D>("2E_mixedME_OSll_dr","mixed OS CR dR(ll)",64,0.,6.4);
  histo1d_["2E_mixedME_OSll_Et"] = fs->make<TH1D>("2E_mixedME_OSll_Et","mixed OS CR ME E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_mixedME_OSll_Et_noCorr"] = fs->make<TH1D>("2E_mixedME_OSll_Et_noCorr","mixed OS CR ME E_{T} (no S&S corr);GeV;",200,0.,1000.);
  histo1d_["2E_mixedME_OSll_eta"] = fs->make<TH1D>("2E_mixedME_OSll_eta","mixed OS CR ME #eta",100,-2.5,2.5);
  histo1d_["2E_mixedAntiME_OSll_Et_xFF"] = fs->make<TH1D>("2E_mixedAntiME_OSll_Et_xFF","mixed OS CR anti ME E_{T} x Fake factor;GeV;",100,0.,500.);
  histo1d_["2E_mixedAntiME_OSll_eta_xFF"] = fs->make<TH1D>("2E_mixedAntiME_OSll_eta_xFF","mixed OS CR anti ME #eta x Fake factor",100,-2.5,2.5);
  histo1d_["2E_mixedAntiME_OSll_Et_xFF_up"] = fs->make<TH1D>("2E_mixedAntiME_OSll_Et_xFF_up","mixed OS CR anti ME E_{T} x Fake factor (up);GeV;",100,0.,500.);
  histo1d_["2E_mixedAntiME_OSll_eta_xFF_up"] = fs->make<TH1D>("2E_mixedAntiME_OSll_eta_xFF_up","mixed OS CR anti ME #eta x Fake factor (up)",100,-2.5,2.5);
  histo1d_["2E_mixedAntiME_OSll_Et_xFF_dn"] = fs->make<TH1D>("2E_mixedAntiME_OSll_Et_xFF_dn","mixed OS CR anti ME E_{T} x Fake factor (down);GeV;",100,0.,500.);
  histo1d_["2E_mixedAntiME_OSll_eta_xFF_dn"] = fs->make<TH1D>("2E_mixedAntiME_OSll_eta_xFF_dn","mixed OS CR anti ME #eta x Fake factor (down)",100,-2.5,2.5);

  histo1d_["2E_antiME1_Et"] = fs->make<TH1D>("2E_antiME1_Et","anti ME1 E_{T};GeV;",200,0.,500.);
  histo1d_["2E_antiME1_eta"] = fs->make<TH1D>("2E_antiME1_eta","anti ME1 #eta",200,-2.5,2.5);
  histo1d_["2E_antiME1_phi"] = fs->make<TH1D>("2E_antiME1_phi","anti ME1 #phi",128,-3.2,3.2);

  histo1d_["2E_antiME1_noGsf_Et"] = fs->make<TH1D>("2E_antiME1_noGsf_Et","anti ME1 (w/o GSF) E_{T};GeV;",200,0.,500.);
  histo1d_["2E_antiME1_noGsf_eta"] = fs->make<TH1D>("2E_antiME1_noGsf_eta","anti ME1 (w/o GSF) #eta",200,-2.5,2.5);
  histo1d_["2E_antiME1_noGsf_phi"] = fs->make<TH1D>("2E_antiME1_noGsf_phi","anti ME1 (w/o GSF) #phi",128,-3.2,3.2);

  histo1d_["2E_antiME2_Et"] = fs->make<TH1D>("2E_antiME2_Et","anti ME2 E_{T};GeV;",200,0.,500.);
  histo1d_["2E_antiME2_eta"] = fs->make<TH1D>("2E_antiME2_eta","anti ME2 #eta",200,-2.5,2.5);
  histo1d_["2E_antiME2_phi"] = fs->make<TH1D>("2E_antiME2_phi","anti ME2 #phi",128,-3.2,3.2);

  histo1d_["2E_antiME2_noGsf_Et"] = fs->make<TH1D>("2E_antiME2_noGsf_Et","anti ME2 (w/o GSF) E_{T};GeV;",200,0.,500.);
  histo1d_["2E_antiME2_noGsf_eta"] = fs->make<TH1D>("2E_antiME2_noGsf_eta","anti ME2 (w/o GSF) #eta",200,-2.5,2.5);
  histo1d_["2E_antiME2_noGsf_phi"] = fs->make<TH1D>("2E_antiME2_noGsf_phi","anti ME2 (w/o GSF) #phi",128,-3.2,3.2);

  histo1d_["2E_antiME_ll_invM"] = fs->make<TH1D>("2E_antiME_ll_invM","anti M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_ll_pt"] = fs->make<TH1D>("2E_antiME_ll_pt","anti p_{T}(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_ll_dr"] = fs->make<TH1D>("2E_antiME_ll_dr","anti dR(ll)",128,0.,6.4);

  histo1d_["2E_antiME_SSll_invM_CR"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR","SS anti M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF","SS anti M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF_preCorr"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF_preCorr","SS anti M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF_up"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF_up","SS anti M(ll) x Fake factor (up);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF_dn"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF_dn","SS anti M(ll) x Fake factor (down);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_SSll_pt"] = fs->make<TH1D>("2E_antiME_SSll_pt","SS anti p_{T}(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_dr"] = fs->make<TH1D>("2E_antiME_SSll_dr","SS anti dR(ll)",128,0.,6.4);
  histo1d_["2E_antiME_SSll_Et_noCorr"] = fs->make<TH1D>("2E_antiME_SSll_Et_noCorr","SS CR anti ME E_{T} (no S&S corr);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_eta"] = fs->make<TH1D>("2E_antiME_SSll_eta","SS CR anti ME #eta",200,-2.5,2.5);
  histo1d_["2E_antiME_SSll_Et_xFF"] = fs->make<TH1D>("2E_antiME_SSll_Et_xFF","SS CR anti ME E_{T} x Fake factor;GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_eta_xFF"] = fs->make<TH1D>("2E_antiME_SSll_eta_xFF","SS CR anti ME #eta x Fake factor",200,-2.5,2.5);
  histo1d_["2E_antiME_SSll_Et_xFF_up"] = fs->make<TH1D>("2E_antiME_SSll_Et_xFF_up","SS CR anti ME E_{T} x Fake factor (up);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_eta_xFF_up"] = fs->make<TH1D>("2E_antiME_SSll_eta_xFF_up","SS CR anti ME #eta x Fake factor (up)",200,-2.5,2.5);
  histo1d_["2E_antiME_SSll_Et_xFF_dn"] = fs->make<TH1D>("2E_antiME_SSll_Et_xFF_dn","SS CR anti ME E_{T} x Fake factor (down);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_eta_xFF_dn"] = fs->make<TH1D>("2E_antiME_SSll_eta_xFF_dn","SS CR anti ME #eta x Fake factor (down)",200,-2.5,2.5);

  histo1d_["2E_antiME_OSll_invM_CR"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR","OS anti M(ll);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF","OS anti M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF_preCorr"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF_preCorr","OS anti M(ll) x Fake factor;GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF_up"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF_up","OS anti M(ll) x Fake factor (up);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF_dn"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF_dn","OS anti M(ll) x Fake factor (down);GeV;",1000,0.,2500.);
  histo1d_["2E_antiME_OSll_pt"] = fs->make<TH1D>("2E_antiME_OSll_pt","OS anti p_{T}(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_dr"] = fs->make<TH1D>("2E_antiME_OSll_dr","OS anti dR(ll)",128,0.,6.4);
  histo1d_["2E_antiME_OSll_Et_noCorr"] = fs->make<TH1D>("2E_antiME_OSll_Et_noCorr","OS CR anti ME E_{T} (no S&S corr);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_eta"] = fs->make<TH1D>("2E_antiME_OSll_eta","OS CR anti ME #eta",200,-2.5,2.5);
  histo1d_["2E_antiME_OSll_Et_xFF"] = fs->make<TH1D>("2E_antiME_OSll_Et_xFF","OS CR anti ME E_{T} x Fake factor;GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_eta_xFF"] = fs->make<TH1D>("2E_antiME_OSll_eta_xFF","OS CR anti ME #eta x Fake factor",200,-2.5,2.5);
  histo1d_["2E_antiME_OSll_Et_xFF_up"] = fs->make<TH1D>("2E_antiME_OSll_Et_xFF_up","OS CR anti ME E_{T} x Fake factor (up);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_eta_xFF_up"] = fs->make<TH1D>("2E_antiME_OSll_eta_xFF_up","OS CR anti ME #eta x Fake factor (up)",200,-2.5,2.5);
  histo1d_["2E_antiME_OSll_Et_xFF_dn"] = fs->make<TH1D>("2E_antiME_OSll_Et_xFF_dn","OS CR anti ME E_{T} x Fake factor (down);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_eta_xFF_dn"] = fs->make<TH1D>("2E_antiME_OSll_eta_xFF_dn","OS CR anti ME #eta x Fake factor (down)",200,-2.5,2.5);

  histo1d_["2E_Et_SSCR_EB_antiME"] = fs->make<TH1D>("2E_Et_SSCR_EB_antiME","SSCR EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_antiME_xFF"] = fs->make<TH1D>("2E_Et_SSCR_EB_antiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_antiME_xFF_up"] = fs->make<TH1D>("2E_Et_SSCR_EB_antiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_antiME_xFF_dn"] = fs->make<TH1D>("2E_Et_SSCR_EB_antiME_xFF_dn",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_CRME"] = fs->make<TH1D>("2E_Et_SSCR_EB_CRME","SSCR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedME"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedME","SSCR (mixed) EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedAntiME"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedAntiME","SSCR (mixed) EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedAntiME_xFF"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedAntiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedAntiME_xFF_up"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedAntiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedAntiME_xFF_dn"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedAntiME_xFF_dn",";GeV;",200,0.,1000.);

  histo1d_["2E_Et_OSCR_EB_antiME"] = fs->make<TH1D>("2E_Et_OSCR_EB_antiME","OSCR EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_antiME_xFF"] = fs->make<TH1D>("2E_Et_OSCR_EB_antiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_antiME_xFF_up"] = fs->make<TH1D>("2E_Et_OSCR_EB_antiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_antiME_xFF_dn"] = fs->make<TH1D>("2E_Et_OSCR_EB_antiME_xFF_dn",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_CRME"] = fs->make<TH1D>("2E_Et_OSCR_EB_CRME","OSCR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedME"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedME","OSCR (mixed) EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedAntiME"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedAntiME","OSCR (mixed) EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedAntiME_xFF"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedAntiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedAntiME_xFF_up"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedAntiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedAntiME_xFF_dn"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedAntiME_xFF_dn",";GeV;",200,0.,1000.);

  histo1d_["2E_Et_OSCR_EB0p8_antiME"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_antiME","OSCR EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_antiME_xFF"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_antiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_antiME_xFF_up"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_antiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_antiME_xFF_dn"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_antiME_xFF_dn",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_CRME"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_CRME","OSCR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_mixedME"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_mixedME","OSCR (mixed) EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_mixedAntiME","OSCR (mixed) EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME_xFF"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_mixedAntiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME_xFF_up"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_mixedAntiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME_xFF_dn"] = fs->make<TH1D>("2E_Et_OSCR_EB0p8_mixedAntiME_xFF_dn",";GeV;",200,0.,1000.);

  histo1d_["2E_Et_OSCR_EB1p5_antiME"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_antiME","OSCR EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_antiME_xFF"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_antiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_antiME_xFF_up"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_antiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_antiME_xFF_dn"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_antiME_xFF_dn",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_CRME"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_CRME","OSCR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_mixedME"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_mixedME","OSCR (mixed) EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_mixedAntiME","OSCR (mixed) EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME_xFF"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_mixedAntiME_xFF",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME_xFF_up"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_mixedAntiME_xFF_up",";GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME_xFF_dn"] = fs->make<TH1D>("2E_Et_OSCR_EB1p5_mixedAntiME_xFF_dn",";GeV;",200,0.,1000.);

  histo1d_["2E_eta_OSCR_EB_antiME"] = fs->make<TH1D>("2E_eta_OSCR_EB_antiME",";#eta;",40,-2.,2.);
  histo1d_["2E_eta_OSCR_EB_antiME_xFF"] = fs->make<TH1D>("2E_eta_OSCR_EB_antiME_xFF",";#eta;",40,-2.,2.);
  histo1d_["2E_eta_OSCR_EB_CRME"] = fs->make<TH1D>("2E_eta_OSCR_EB_CRME",";#eta;",40,-2.,2.);
  histo1d_["2E_eta_OSCR_EB_mixedME"] = fs->make<TH1D>("2E_eta_OSCR_EB_mixedME",";#eta;",40,-2.,2.);
  histo1d_["2E_eta_OSCR_EB_mixedAntiME"] = fs->make<TH1D>("2E_eta_OSCR_EB_mixedAntiME",";#eta;",40,-2.,2.);
  histo1d_["2E_eta_OSCR_EB_mixedAntiME_xFF"] = fs->make<TH1D>("2E_eta_OSCR_EB_mixedAntiME_xFF",";#eta;",40,-2.,2.);

  std::vector<double> xbinsOS = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                                 225, 250, 275, 300, 400, 500, 700, 1000};
  const int nbinsOS = xbinsOS.size()-1;

  histo2d_["2E_Et_eta_OSCR_EB_antiME"] = fs->make<TH2D>("2E_Et_eta_OSCR_EB_antiME",";#eta;E_{T}",40,-2.,2.,nbinsOS,&(xbinsOS[0]));
  histo2d_["2E_Et_eta_OSCR_EB_antiME_xFF"] = fs->make<TH2D>("2E_Et_eta_OSCR_EB_antiME_xFF",";#eta;E_{T}",40,-2.,2.,nbinsOS,&(xbinsOS[0]));
  histo2d_["2E_Et_eta_OSCR_EB_CRME"] = fs->make<TH2D>("2E_Et_eta_OSCR_EB_CRME",";#eta;E_{T}",40,-2.,2.,nbinsOS,&(xbinsOS[0]));
  histo2d_["2E_Et_eta_OSCR_EB_mixedME"] = fs->make<TH2D>("2E_Et_eta_OSCR_EB_mixedME",";#eta;E_{T}",40,-2.,2.,nbinsOS,&(xbinsOS[0]));
  histo2d_["2E_Et_eta_OSCR_EB_mixedAntiME"] = fs->make<TH2D>("2E_Et_eta_OSCR_EB_mixedAntiME",";#eta;E_{T}",40,-2.,2.,nbinsOS,&(xbinsOS[0]));
  histo2d_["2E_Et_eta_OSCR_EB_mixedAntiME_xFF"] = fs->make<TH2D>("2E_Et_eta_OSCR_EB_mixedAntiME_xFF",";#eta;E_{T}",40,-2.,2.,nbinsOS,&(xbinsOS[0]));

  histo1d_["3E_CRME_Et"] = fs->make<TH1D>("3E_CRME_Et","CR ME E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CRME_eta"] = fs->make<TH1D>("3E_CRME_eta","CR ME #eta",50,-2.5,2.5);
  histo1d_["3E_CRME_phi"] = fs->make<TH1D>("3E_CRME_phi","CR ME #phi",64,-3.2,3.2);

  histo1d_["3E_CRME_noGsf_Et"] = fs->make<TH1D>("3E_CRME_noGsf_Et","CR ME (w/o GSF) E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CRME_noGsf_eta"] = fs->make<TH1D>("3E_CRME_noGsf_eta","CR ME (w/o GSF) #eta",50,-2.5,2.5);
  histo1d_["3E_CRME_noGsf_phi"] = fs->make<TH1D>("3E_CRME_noGsf_phi","CR ME (w/o GSF) #phi",64,-3.2,3.2);

  histo1d_["3E_CRME_all_Et"] = fs->make<TH1D>("3E_CRME_all_Et","CR ME (w/ & w/o GSF) E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CRME_all_eta"] = fs->make<TH1D>("3E_CRME_all_eta","CR ME (w/ & w/o GSF) #eta",50,-2.5,2.5);
  histo1d_["3E_CRME_all_phi"] = fs->make<TH1D>("3E_CRME_all_phi","CR ME (w/ & w/o GSF) #phi",64,-3.2,3.2);

  histo1d_["3E_CR_NME1_Et"] = fs->make<TH1D>("3E_CR_NME1_Et","CR nonME1 E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CR_NME1_eta"] = fs->make<TH1D>("3E_CR_NME1_eta","CR nonME1 #eta",50,-2.5,2.5);
  histo1d_["3E_CR_NME1_phi"] = fs->make<TH1D>("3E_CR_NME1_phi","CR nonME1 #phi",64,-3.2,3.2);

  histo1d_["3E_CR_NME2_Et"] = fs->make<TH1D>("3E_CR_NME2_Et","CR nonME2 E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CR_NME2_eta"] = fs->make<TH1D>("3E_CR_NME2_eta","CR nonME2 #eta",50,-2.5,2.5);
  histo1d_["3E_CR_NME2_phi"] = fs->make<TH1D>("3E_CR_NME2_phi","CR nonME2 #phi",64,-3.2,3.2);

  histo1d_["3E_CRME_lll_invM"] = fs->make<TH1D>("3E_CRME_lll_invM","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_mergedEleScale"] = fs->make<TH1D>("3E_CRME_lll_invM_mergedEleScale","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_mergedEleSmear"] = fs->make<TH1D>("3E_CRME_lll_invM_mergedEleSmear","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_scaleUp"] = fs->make<TH1D>("3E_CRME_lll_invM_scaleUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_scaleDn"] = fs->make<TH1D>("3E_CRME_lll_invM_scaleDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_sigmaUp"] = fs->make<TH1D>("3E_CRME_lll_invM_sigmaUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_sigmaDn"] = fs->make<TH1D>("3E_CRME_lll_invM_sigmaDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_heepIdUp"] = fs->make<TH1D>("3E_CRME_lll_invM_heepIdUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_heepIdDn"] = fs->make<TH1D>("3E_CRME_lll_invM_heepIdDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_elRecoUp"] = fs->make<TH1D>("3E_CRME_lll_invM_elRecoUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_elRecoDn"] = fs->make<TH1D>("3E_CRME_lll_invM_elRecoDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_elTrigUp"] = fs->make<TH1D>("3E_CRME_lll_invM_elTrigUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_elTrigDn"] = fs->make<TH1D>("3E_CRME_lll_invM_elTrigDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_mergedEleIdUp"] = fs->make<TH1D>("3E_CRME_lll_invM_mergedEleIdUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_invM_mergedEleIdDn"] = fs->make<TH1D>("3E_CRME_lll_invM_mergedEleIdDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_CRME_lll_pt"] = fs->make<TH1D>("3E_CRME_lll_pt","p_{T}(lll);GeV;",100,0.,500.);
  histo1d_["3E_CRME_lll_dr_l1ME"] = fs->make<TH1D>("3E_CRME_lll_dr_l1ME","dR(l1ME)",64,0.,6.4);
  histo1d_["3E_CRME_lll_dr_l2ME"] = fs->make<TH1D>("3E_CRME_lll_dr_l2ME","dR(l2ME)",64,0.,6.4);
  histo1d_["3E_CRME_lll_dr_l1l2"] = fs->make<TH1D>("3E_CRME_lll_dr_l1l2","dR(l1l2)",64,0.,6.4);

  histo1d_["3E_antiME_Et"] = fs->make<TH1D>("3E_antiME_Et","anti ME E_{T};GeV;",200,0.,500.);
  histo1d_["3E_antiME_eta"] = fs->make<TH1D>("3E_antiME_eta","anti ME #eta",200,-2.5,2.5);
  histo1d_["3E_antiME_phi"] = fs->make<TH1D>("3E_antiME_phi","anti ME #phi",128,-3.2,3.2);

  histo1d_["3E_antiME_noGsf_Et"] = fs->make<TH1D>("3E_antiME_noGsf_Et","anti ME (w/o GSF) E_{T};GeV;",200,0.,500.);
  histo1d_["3E_antiME_noGsf_eta"] = fs->make<TH1D>("3E_antiME_noGsf_eta","anti ME (w/o GSF) #eta",200,-2.5,2.5);
  histo1d_["3E_antiME_noGsf_phi"] = fs->make<TH1D>("3E_antiME_noGsf_phi","anti ME (w/o GSF) #phi",128,-3.2,3.2);

  histo1d_["3E_antiME_NME1_Et"] = fs->make<TH1D>("3E_antiME_NME1_Et","anti nonME1 E_{T};GeV;",200,0.,500.);
  histo1d_["3E_antiME_NME1_eta"] = fs->make<TH1D>("3E_antiME_NME1_eta","anti nonME1 #eta",200,-2.5,2.5);
  histo1d_["3E_antiME_NME1_phi"] = fs->make<TH1D>("3E_antiME_NME1_phi","anti nonME1 #phi",128,-3.2,3.2);

  histo1d_["3E_antiME_NME2_Et"] = fs->make<TH1D>("3E_antiME_NME2_Et","anti nonME2 E_{T};GeV;",200,0.,500.);
  histo1d_["3E_antiME_NME2_eta"] = fs->make<TH1D>("3E_antiME_NME2_eta","anti nonME2 #eta",200,-2.5,2.5);
  histo1d_["3E_antiME_NME2_phi"] = fs->make<TH1D>("3E_antiME_NME2_phi","anti nonME2 #phi",128,-3.2,3.2);

  histo1d_["3E_antiME_lll_invM_CR"] = fs->make<TH1D>("3E_antiME_lll_invM_CR","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_heepIdUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_heepIdUp","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_heepIdDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_heepIdDn","CR M(lll);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_preCorr"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_preCorr","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_heepIdUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_heepIdUp","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_heepIdDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_heepIdDn","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_elRecoUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_elRecoUp","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_elRecoDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_elRecoDn","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_elTrigUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_elTrigUp","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_elTrigDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_elTrigDn","CR M(lll) x SS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_preCorr"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_preCorr","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_heepIdUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_heepIdUp","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_heepIdDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_heepIdDn","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_elRecoUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_elRecoUp","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_elRecoDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_elRecoDn","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_elTrigUp"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_elTrigUp","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_elTrigDn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_elTrigDn","CR M(lll) x OS Fake factor;GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_up"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_up","CR M(lll) x SS Fake factor (up);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_up"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_up","CR M(lll) x OS Fake factor (up);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_dn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_dn","CR M(lll) x SS Fake factor (down);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_dn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_dn","CR M(lll) x OS Fake factor (down);GeV;",1000,0.,2500.);
  histo1d_["3E_antiME_lll_pt"] = fs->make<TH1D>("3E_antiME_lll_pt","p_{T}(lll);GeV;",200,0.,1000.);
  histo1d_["3E_antiME_lll_dr_l1ME"] = fs->make<TH1D>("3E_antiME_lll_dr_l1ME","dR(l1ME)",128,0.,6.4);
  histo1d_["3E_antiME_lll_dr_l2ME"] = fs->make<TH1D>("3E_antiME_lll_dr_l2ME","dR(l2ME)",128,0.,6.4);
  histo1d_["3E_antiME_lll_dr_l1l2"] = fs->make<TH1D>("3E_antiME_lll_dr_l1l2","dR(l1l2)",128,0.,6.4);
  histo1d_["3E_antiME_Et_xSSFF"] = fs->make<TH1D>("3E_antiME_Et_xSSFF","CR anti ME E_{T} x SS Fake factor;GeV;",80,0.,200.);
  histo1d_["3E_antiME_eta_xSSFF"] = fs->make<TH1D>("3E_antiME_eta_xSSFF","CR anti ME #eta x SS Fake factor",100,-2.5,2.5);
  histo1d_["3E_antiME_Et_xOSFF"] = fs->make<TH1D>("3E_antiME_Et_xOSFF","CR anti ME E_{T} x OS Fake factor;GeV;",80,0.,200.);
  histo1d_["3E_antiME_eta_xOSFF"] = fs->make<TH1D>("3E_antiME_eta_xOSFF","CR anti ME #eta x OS Fake factor",100,-2.5,2.5);
  histo1d_["3E_antiME_Et_xSSFF_up"] = fs->make<TH1D>("3E_antiME_Et_xSSFF_up","CR anti ME E_{T} x SS Fake factor (up);GeV;",80,0.,200.);
  histo1d_["3E_antiME_eta_xSSFF_up"] = fs->make<TH1D>("3E_antiME_eta_xSSFF_up","CR anti ME #eta x SS Fake factor (up)",100,-2.5,2.5);
  histo1d_["3E_antiME_Et_xOSFF_up"] = fs->make<TH1D>("3E_antiME_Et_xOSFF_up","CR anti ME E_{T} x OS Fake factor (up);GeV;",80,0.,200.);
  histo1d_["3E_antiME_eta_xOSFF_up"] = fs->make<TH1D>("3E_antiME_eta_xOSFF_up","CR anti ME #eta x OS Fake factor (up)",100,-2.5,2.5);
  histo1d_["3E_antiME_Et_xSSFF_dn"] = fs->make<TH1D>("3E_antiME_Et_xSSFF_dn","CR anti ME E_{T} x SS Fake factor (down);GeV;",80,0.,200.);
  histo1d_["3E_antiME_eta_xSSFF_dn"] = fs->make<TH1D>("3E_antiME_eta_xSSFF_dn","CR anti ME #eta x SS Fake factor (down)",100,-2.5,2.5);
  histo1d_["3E_antiME_Et_xOSFF_dn"] = fs->make<TH1D>("3E_antiME_Et_xOSFF_dn","CR anti ME E_{T} x OS Fake factor (down);GeV;",80,0.,200.);
  histo1d_["3E_antiME_eta_xOSFF_dn"] = fs->make<TH1D>("3E_antiME_eta_xOSFF_dn","CR anti ME #eta x OS Fake factor (down)",100,-2.5,2.5);

  histo1d_["3E_Et_CR_EB_antiME"] = fs->make<TH1D>("3E_Et_CR_EB_antiME",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_antiME_xSSFF"] = fs->make<TH1D>("3E_Et_CR_EB_antiME_xSSFF",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_antiME_xSSFF_up"] = fs->make<TH1D>("3E_Et_CR_EB_antiME_xSSFF_up",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_antiME_xSSFF_dn"] = fs->make<TH1D>("3E_Et_CR_EB_antiME_xSSFF_dn",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_antiME_xOSFF"] = fs->make<TH1D>("3E_Et_CR_EB_antiME_xOSFF",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_antiME_xOSFF_up"] = fs->make<TH1D>("3E_Et_CR_EB_antiME_xOSFF_up",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_antiME_xOSFF_dn"] = fs->make<TH1D>("3E_Et_CR_EB_antiME_xOSFF_dn",";GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_CRME"] = fs->make<TH1D>("3E_Et_CR_EB_CRME","3EVR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["3E_eta_CR_EB_CRME"] = fs->make<TH1D>("3E_eta_CR_EB_CRME",";#eta;",100,-2.5,2.5);

  SRtree_ = fs->make<TTree>("SRtree","SRtree");
  SRtree_->Branch("invM",&SRinvM_,"invM/F");
  SRtree_->Branch("wgt",&SRwgt_,"wgt/F");
  SRtree_->Branch("isOS",&SRisOS_,"isOS/I");
  SRtree_->Branch("pt1",&SRpt1_,"pt1/F");
  SRtree_->Branch("eta1",&SReta1_,"eta1/F");
  SRtree_->Branch("phi1",&SRphi1_,"phi1/F");
  SRtree_->Branch("pt2",&SRpt2_,"pt2/F");
  SRtree_->Branch("eta2",&SReta2_,"eta2/F");
  SRtree_->Branch("phi2",&SRphi2_,"phi2/F");

  mixedCRtree_ = fs->make<TTree>("mixedCRtree","mixedCRtree");
  mixedCRtree_->Branch("invM",&mixedCRinvM_,"invM/F");
  mixedCRtree_->Branch("wgt",&mixedCRwgt_,"wgt/F");
  mixedCRtree_->Branch("isOS",&mixedCRisOS_,"isOS/I");
  mixedCRtree_->Branch("pt1",&mixedCRpt1_,"pt1/F");
  mixedCRtree_->Branch("eta1",&mixedCReta1_,"eta1/F");
  mixedCRtree_->Branch("phi1",&mixedCRphi1_,"phi1/F");
  mixedCRtree_->Branch("pt2",&mixedCRpt2_,"pt2/F");
  mixedCRtree_->Branch("eta2",&mixedCReta2_,"eta2/F");
  mixedCRtree_->Branch("phi2",&mixedCRphi2_,"phi2/F");

  antiCRtree_ = fs->make<TTree>("antiCRtree","antiCRtree");
  antiCRtree_->Branch("invM",&antiCRinvM_,"invM/F");
  antiCRtree_->Branch("wgt",&antiCRwgt_,"wgt/F");
  antiCRtree_->Branch("isOS",&antiCRisOS_,"isOS/I");
  antiCRtree_->Branch("pt1",&antiCRpt1_,"pt1/F");
  antiCRtree_->Branch("eta1",&antiCReta1_,"eta1/F");
  antiCRtree_->Branch("phi1",&antiCRphi1_,"phi1/F");
  antiCRtree_->Branch("pt2",&antiCRpt2_,"pt2/F");
  antiCRtree_->Branch("eta2",&antiCReta2_,"eta2/F");
  antiCRtree_->Branch("phi2",&antiCRphi2_,"phi2/F");
}

void MergedEleCRanalyzer::endJob() {
  FFfile_->Close();
}

void MergedEleCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
        histo1d_["PUsummary"]->Fill( static_cast<float>(apu->getTrueNumInteractions())+0.5, aWeight );

        break;
      }
    }
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  if (selectZpT_ || selectHT_) {
    edm::Handle<LHEEventProduct> lheEventHandle;
    iEvent.getByToken(lheToken_, lheEventHandle);

    const auto& hepeup = lheEventHandle->hepeup();
    const auto& pup = hepeup.PUP;
    int lep = -1, lepBar = -1, nu = -1, nuBar = -1;
    double lheVpt = -1.;
    double lheHT = -1.;

    for (unsigned int i = 0, n = pup.size(); i < n; ++i) {
      int status = hepeup.ISTUP[i];
      int idabs = std::abs(hepeup.IDUP[i]);

      if ((status == 1) && ((idabs == 21) || (idabs > 0 && idabs < 7))) { //# gluons and quarks
        double pt = std::hypot(pup[i][0], pup[i][1]);  // first entry is px, second py
        lheHT += pt;
      } else if (idabs == 12 || idabs == 14 || idabs == 16) {
        (hepeup.IDUP[i] > 0 ? nu : nuBar) = i;
      } else if (idabs == 11 || idabs == 13 || idabs == 15) {
        (hepeup.IDUP[i] > 0 ? lep : lepBar) = i;
      }
    }

    std::pair<int, int> v(0, 0);
    if (lep != -1 && lepBar != -1)
      v = std::make_pair(lep, lepBar);
    else if (lep != -1 && nuBar != -1)
      v = std::make_pair(lep, nuBar);
    else if (nu != -1 && lepBar != -1)
      v = std::make_pair(nu, lepBar);
    else if (nu != -1 && nuBar != -1)
      v = std::make_pair(nu, nuBar);
    if (v.first != -1 && v.second != -1) {
      lheVpt = std::hypot(pup[v.first][0] + pup[v.second][0], pup[v.first][1] + pup[v.second][1]);
    }

    histo1d_["lheHT"]->Fill( lheHT );
    histo1d_["lheVpT"]->Fill( lheVpt );

    if (selectHT_) {
      if (lheHT > maxHT_)
        return;
    }

    if (selectZpT_) {
      if (lheVpt > maxVpT_)
        return;
    }

    histo1d_["lheHT_cut"]->Fill( lheHT );
    histo1d_["lheVpT_cut"]->Fill( lheVpt );
  }

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  histo1d_["cutflow_4E"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_3E"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_2E"]->Fill( 0.5, aWeight );

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

  histo1d_["cutflow_4E"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_3E"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_2E"]->Fill( 1.5, aWeight );

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dPhiInHandle;
  iEvent.getByToken(union5x5dPhiInToken_, union5x5dPhiInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  // first two (modified)HEEP electrons should be Et > 35 GeV
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

    std::string strModHEEP = "N-1_modifiedHEEP";
    std::string strHEEP = "N-1_HEEP";

    if (aEle->isEB()) {
      strModHEEP += "_EB";
      strHEEP += "_EB";
    } else {
      strModHEEP += "_EE";
      strHEEP += "_EE";
    }

    for (unsigned jdx = 0; jdx < 12; jdx++) {
      int32_t mask = 0x00000001 << jdx;
      int32_t bitmap = aEle->userInt("modifiedHeepElectronID");
      int32_t nm1 = mask | bitmap;

      if ( nm1 == 0x00000FFF )
        histo1d_[strModHEEP]->Fill(jdx, aWeight);
    }

    for (unsigned jdx = 0; jdx < 12; jdx++) {
      int32_t mask = 0x00000001 << jdx;
      int32_t bitmap = aEle->userInt("heepElectronID-HEEPV70");
      int32_t nm1 = mask | bitmap;

      if ( nm1 == 0x00000FFF )
        histo1d_[strHEEP]->Fill(jdx, aWeight);
    }

    bool selected = false;

    if ( aEle->electronID("modifiedHeepElectronID") ) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      acceptEles.push_back(castEle);
      selected = true;
    }

    if (!selected) {
      int32_t bitmap = aEle->userInt("modifiedHeepElectronID");
      int32_t mask = 0x000007B0; // = 0111 1011 0000
      int32_t pass = bitmap | mask;
      bool passMaskedId = pass==0x00000FFF; // HEEP ID has 12 cuts

      if (passMaskedId) {
        auto castEle = aEle.castTo<pat::ElectronRef>();
        nonHeepEles.push_back(castEle);
      }
    }
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByEt);

  bool accepted = false;

  if ( acceptEles.size() > 1 ) {
    if ( acceptEles.at(1)->et() > etThres1_ )
      accepted = true;
  }

  if (!accepted)
    return;

  if (!nonHeepEles.empty())
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

  std::vector<pat::ElectronRef> CRMEs;
  std::vector<pat::ElectronRef> antiMEs;

  histo1d_["nPV"]->Fill( static_cast<float>(pvHandle->size())+0.5, aWeight );

  auto estimateU5x5Eta = [&union5x5dEtaInHandle] (const pat::ElectronRef& aEle) -> float {
    const float eta1stGSF = -(aEle->deltaEtaSeedClusterTrackAtVtx() - aEle->superCluster()->seed()->eta());
    const float u5x5Eta = (*union5x5dEtaInHandle)[aEle] + eta1stGSF;

    return u5x5Eta;
  };

  auto estimateU5x5Et = [&union5x5dEtaInHandle,&union5x5EnergyHandle,&estimateU5x5Eta] (const pat::ElectronRef& aEle) -> float {
    const float u5x5Eta = estimateU5x5Eta(aEle);
    const float u5x5Et = (*union5x5EnergyHandle)[aEle]/std::cosh(u5x5Eta);

    return u5x5Et;
  };

  auto lvecME = [&estimateU5x5Et,&estimateU5x5Eta,&union5x5dEtaInHandle,&union5x5dPhiInHandle] (const pat::ElectronRef& aEle) {
    if ( aEle->userInt("mvaMergedElectronCategories")==1 )
      return math::PtEtaPhiMLorentzVector(aEle->p4());

    const float u5x5Eta = estimateU5x5Eta(aEle);

    const float phi1stGSF = reco::reduceRange( -( aEle->deltaPhiSuperClusterTrackAtVtx() - aEle->superCluster()->phi() ) );
    const float u5x5Phi = reco::reduceRange( (*union5x5dPhiInHandle)[aEle] + phi1stGSF );

    return math::PtEtaPhiMLorentzVector(estimateU5x5Et(aEle),u5x5Eta,u5x5Phi,0.);
  };

  for (unsigned int idx = 0; idx < acceptEles.size(); ++idx) {
    const auto aEle = acceptEles.at(idx);

    const auto& orgGsfTrk = aEle->gsfTrack();
    const auto& addGsfTrk = (*addGsfTrkMap)[aEle];
    const auto& addPackedCand = (*addPackedCandHandle)[aEle];

    const reco::Track* addTrk = addGsfTrk.get();

    if ( addGsfTrk==orgGsfTrk && addPackedCand.isNonnull() )
      addTrk = addPackedCand->bestTrack();

    const float u5x5Et = estimateU5x5Et(aEle);

    if ( u5x5Et < 50. )
      continue;

    bool isNotMerged = MergedLeptonHelperFct::isNotMerged(aEle,eleHandle,addGsfTrk);

    if ( aEle->electronID("mvaMergedElectron") && !isNotMerged )
      CRMEs.push_back(aEle);

    // electrons with mva score = -1 are out of our interest
    if ( !static_cast<bool>(aEle->electronID("mvaMergedElectron")) && aEle->userFloat("mvaMergedElectronValues")!=-1. && !isNotMerged )
      antiMEs.push_back(aEle);

    // for filling histograms
    if ( acceptEles.size()==2 || acceptEles.size()==3 ) {
      if ( aEle->userInt("mvaMergedElectronCategories")==1 ) { // no 2nd GSF
        histo1d_["mva_bkgEt2EB"]->Fill( aEle->userFloat("mvaMergedElectronValues") , aWeight );
      } else { // has 2nd GSF
        if (!isNotMerged)
          histo1d_["mva_HasTrkEB"]->Fill( aEle->userFloat("mvaMergedElectronValues") , aWeight );

        if ( aEle->electronID("mvaMergedElectron") ) {
          const auto lvecE = math::PtEtaPhiMLorentzVector(aEle->gsfTrack()->pt(),aEle->gsfTrack()->eta(),aEle->gsfTrack()->phi(),0.);
          const auto lvecGSF = math::PtEtaPhiMLorentzVector(addTrk->pt(),addTrk->eta(),addTrk->phi(),0.);
          const auto lvecMEll = lvecE + lvecGSF;

          if ( std::abs(aEle->superCluster()->eta()) < 1.5 )
            histo1d_["invM_ll_CRME_EB"]->Fill( lvecMEll.M() , aWeight );
        }
      }
    }
  } // electron loop

  std::sort(CRMEs.begin(),CRMEs.end(),sortByEt);
  std::sort(antiMEs.begin(),antiMEs.end(),sortByEt);

  // SF
  if (isMC_) {
    aWeight *= systHelperEle_.GetTrigSF(acceptEles.front());
    aWeight *= systHelperEle_.GetTrigUnseededSF(acceptEles.at(1));

    for (const auto& aEle : acceptEles) {
      aWeight *= systHelperEle_.GetRecoSF(aEle);
      aWeight *= systHelperEle_.GetModifiedHeepSF(aEle);
    }

    for (const auto& aEle : CRMEs)
      aWeight *= systHelperEle_.GetMergedEleSF(aEle);
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

  auto mergedEleSFcl95 = [this,&estimateU5x5Et,&CRMEs] (const double wgt) -> std::pair<double,double> {
    if (!isMC_)
      return std::make_pair(wgt,wgt);

    double wgtUp = wgt;
    double wgtDn = wgt;

    for (const auto& aEle : CRMEs) {
      auto apair = systHelperEle_.GetMergedEleSFcl95UpDn(aEle,estimateU5x5Et(aEle));
      wgtUp *= apair.first;
      wgtDn *= apair.second;
    } // acceptEles

    return std::make_pair(wgtUp,wgtDn);
  };

  auto recoSFcl95 = [this,&acceptEles] (const double wgt) -> std::pair<double,double> {
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

    return std::make_pair(wgtUp,wgtDn);
  };

  auto trigSFcl95 = [this,&acceptEles,&CRMEs] (const double wgt) -> std::pair<double,double> {
    if (!isMC_)
      return std::make_pair(wgt,wgt);

    double wgtUp = wgt;
    double wgtDn = wgt;

    std::pair<double,double> pair1, pair2;

    if (CRMEs.size()==1) {
      if (CRMEs.front()==acceptEles.front()) {
        pair1 = systHelperEle_.GetTrigSFMergedUpDn(acceptEles.front());
        pair2 = systHelperEle_.GetTrigUnseededSFcl95UpDn(acceptEles.at(1));
      } else if (CRMEs.front()==acceptEles.at(1)) {
        pair1 = systHelperEle_.GetTrigSFcl95UpDn(acceptEles.front());
        pair2 = systHelperEle_.GetTrigSFMergedUpDn(acceptEles.at(1));
      }
    } else if (CRMEs.size()==2) {
      pair1 = systHelperEle_.GetTrigSFMergedUpDn(acceptEles.front());
      pair2 = systHelperEle_.GetTrigSFMergedUpDn(acceptEles.at(1));
    } else {
      pair1 = systHelperEle_.GetTrigSFcl95UpDn(acceptEles.front());
      pair2 = systHelperEle_.GetTrigUnseededSFcl95UpDn(acceptEles.at(1));
    }

    wgtUp *= pair1.first;
    wgtUp *= pair2.first;
    wgtDn *= pair1.second;
    wgtDn *= pair2.second;

    return std::make_pair(wgtUp,wgtDn);
  };

  auto ciSSOS = [this,&estimateU5x5Et,&estimateU5x5Eta] (const pat::ElectronRef& aEle, double& ffSS, double& ffOS) -> std::pair<double,double> {
    const double et = estimateU5x5Et(aEle);
    const double eta = estimateU5x5Eta(aEle);
    ffOS = osboth_->Eval(eta,et);
    ffSS = ssboth_->Eval(et);
    const double xval[1] = {et};
    const double xyval[2] = {eta,et};

    double ciOS[1], ciSS[1];
    osFit_->GetConfidenceIntervals(1,2,1,xyval,ciOS,0.95,false);
    ssFit_->GetConfidenceIntervals(1,1,0,xval,ciSS,0.95,false);

    return std::make_pair( ciSS[0], ciOS[0] );
  };

  switch (acceptEles.size()) {
    case 4:
      histo1d_["cutflow_4E"]->Fill( 2.5, aWeight );
      if ( CRMEs.size()==0 )
        histo1d_["cutflow_4E"]->Fill( 3.5, aWeight );
      break;
    case 3:
      histo1d_["cutflow_3E"]->Fill( 2.5, aWeight );

      if ( CRMEs.size()==1 ) {
        std::vector<pat::ElectronRef> nonMEs;

        for (const auto& aEle : acceptEles) {
          if ( aEle!=CRMEs.front() )
            nonMEs.push_back(aEle);
        }

        auto lvecCorr = [] (const pat::ElectronRef& aEle, const std::string& opt) {
          return aEle->polarP4()*aEle->userFloat(opt)/aEle->energy();
        };

        const auto lvecCRME = lvecME(CRMEs.front());
        const auto lvecCRME_scale = systHelperEle_.mergedEleScale(CRMEs.front())*lvecCRME;
        const auto lvecCRME_smear = systHelperEle_.mergedEleSmear(CRMEs.front(),(*union5x5EnergyHandle)[CRMEs.front()])*lvecCRME;
        const auto lvecNME1 = lvecCorr(nonMEs.at(0),"ecalTrkEnergyPostCorr");
        const auto lvecNME2 = lvecCorr(nonMEs.at(1),"ecalTrkEnergyPostCorr");
        const auto lvecNME1_scaleUp = lvecCorr(nonMEs.at(0),"energyScaleUp");
        const auto lvecNME2_scaleUp = lvecCorr(nonMEs.at(1),"energyScaleUp");
        const auto lvecNME1_scaleDn = lvecCorr(nonMEs.at(0),"energyScaleDown");
        const auto lvecNME2_scaleDn = lvecCorr(nonMEs.at(1),"energyScaleDown");
        const auto lvecNME1_sigmaUp = lvecCorr(nonMEs.at(0),"energySigmaUp");
        const auto lvecNME2_sigmaUp = lvecCorr(nonMEs.at(1),"energySigmaUp");
        const auto lvecNME1_sigmaDn = lvecCorr(nonMEs.at(0),"energySigmaDown");
        const auto lvecNME2_sigmaDn = lvecCorr(nonMEs.at(1),"energySigmaDown");
        const auto lvecCRllME = lvecCRME + lvecNME1 + lvecNME2;
        const double dr2l1ME = reco::deltaR2(lvecCRME.eta(),lvecCRME.phi(),nonMEs.front()->eta(),nonMEs.front()->phi());
        const double dr2l2ME = reco::deltaR2(lvecCRME.eta(),lvecCRME.phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double dr2l1l2 = reco::deltaR2(nonMEs.front()->eta(),nonMEs.front()->phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double mlll = lvecCRllME.M();
        const double mlll_scaleUp = (lvecCRME+lvecNME1_scaleUp+lvecNME2_scaleUp).M();
        const double mlll_scaleDn = (lvecCRME+lvecNME1_scaleDn+lvecNME2_scaleDn).M();
        const double mlll_sigmaUp = (lvecCRME+lvecNME1_sigmaUp+lvecNME2_sigmaUp).M();
        const double mlll_sigmaDn = (lvecCRME+lvecNME1_sigmaDn+lvecNME2_sigmaDn).M();
        const double mlll_mergedEleScale = (lvecCRME_scale+lvecNME1+lvecNME2).M();
        const double mlll_mergedEleSmear = (lvecCRME_smear+lvecNME1+lvecNME2).M();

        if ( mlll > 50. && (mlll < 200. || isMC_) ) {
          histo1d_["3E_CRME_lll_invM"]->Fill(mlll,aWeight);
          histo1d_["3E_CRME_lll_invM_mergedEleScale"]->Fill(mlll_mergedEleScale,aWeight);
          histo1d_["3E_CRME_lll_invM_mergedEleSmear"]->Fill(mlll_mergedEleSmear,aWeight);
          histo1d_["3E_CRME_lll_invM_scaleUp"]->Fill(mlll_scaleUp,aWeight);
          histo1d_["3E_CRME_lll_invM_scaleDn"]->Fill(mlll_scaleDn,aWeight);
          histo1d_["3E_CRME_lll_invM_sigmaUp"]->Fill(mlll_sigmaUp,aWeight);
          histo1d_["3E_CRME_lll_invM_sigmaDn"]->Fill(mlll_sigmaDn,aWeight);
          histo1d_["3E_CRME_lll_invM_heepIdUp"]->Fill(mlll,modHeepSFcl95(aWeight).first);
          histo1d_["3E_CRME_lll_invM_heepIdDn"]->Fill(mlll,modHeepSFcl95(aWeight).second);
          histo1d_["3E_CRME_lll_invM_elRecoUp"]->Fill(mlll,recoSFcl95(aWeight).first);
          histo1d_["3E_CRME_lll_invM_elRecoDn"]->Fill(mlll,recoSFcl95(aWeight).second);
          histo1d_["3E_CRME_lll_invM_elTrigUp"]->Fill(mlll,trigSFcl95(aWeight).first);
          histo1d_["3E_CRME_lll_invM_elTrigDn"]->Fill(mlll,trigSFcl95(aWeight).second);
          histo1d_["3E_CRME_lll_invM_mergedEleIdUp"]->Fill(mlll,mergedEleSFcl95(aWeight).first);
          histo1d_["3E_CRME_lll_invM_mergedEleIdDn"]->Fill(mlll,mergedEleSFcl95(aWeight).second);
        }

        if ( mlll > 50. && mlll < 200. ) {
          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["3E_CRME_Et"]->Fill(lvecCRME.Et(),aWeight);
            histo1d_["3E_CRME_eta"]->Fill(lvecCRME.eta(),aWeight);
            histo1d_["3E_CRME_phi"]->Fill(lvecCRME.phi(),aWeight);
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["3E_CRME_noGsf_Et"]->Fill(lvecCRME.Et(),aWeight);
            histo1d_["3E_CRME_noGsf_eta"]->Fill(lvecCRME.eta(),aWeight);
            histo1d_["3E_CRME_noGsf_phi"]->Fill(lvecCRME.phi(),aWeight);
          } else {}

          histo1d_["3E_CRME_all_Et"]->Fill(lvecCRME.Et(),aWeight);
          histo1d_["3E_CRME_all_eta"]->Fill(lvecCRME.eta(),aWeight);
          histo1d_["3E_CRME_all_phi"]->Fill(lvecCRME.phi(),aWeight);

          histo1d_["3E_CR_NME1_Et"]->Fill(nonMEs.front()->et()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy(),aWeight);
          histo1d_["3E_CR_NME1_eta"]->Fill(nonMEs.front()->eta(),aWeight);
          histo1d_["3E_CR_NME1_phi"]->Fill(nonMEs.front()->phi(),aWeight);

          histo1d_["3E_CR_NME2_Et"]->Fill(nonMEs.at(1)->et()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy(),aWeight);
          histo1d_["3E_CR_NME2_eta"]->Fill(nonMEs.at(1)->eta(),aWeight);
          histo1d_["3E_CR_NME2_phi"]->Fill(nonMEs.at(1)->phi(),aWeight);

          histo1d_["3E_CRME_lll_pt"]->Fill(lvecCRllME.pt(),aWeight);
          histo1d_["3E_CRME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
          histo1d_["3E_CRME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
          histo1d_["3E_CRME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

          histo1d_["3E_eta_CR_EB_CRME"]->Fill( estimateU5x5Eta(CRMEs.front()), aWeight );

          if ( CRMEs.front()->isEB() )
            histo1d_["3E_Et_CR_EB_CRME"]->Fill( estimateU5x5Et(CRMEs.front()), aWeight );
        }
      } // CRMEs.size()==1

      if ( CRMEs.size()==0 && antiMEs.size()>0 ) {
        std::vector<pat::ElectronRef> nonMEs;

        for (const auto& aEle : acceptEles) {
          if ( aEle!=antiMEs.front() )
            nonMEs.push_back(aEle);
        }

        const auto lvecAMEpreCorr = lvecME(antiMEs.front());
        const auto lvecAME = systHelperEle_.GetSingleAbcdScaleSmear(lvecAMEpreCorr)*lvecAMEpreCorr;
        const auto lvecNME1 = nonMEs.front()->polarP4()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy();
        const auto lvecNME2 = nonMEs.at(1)->polarP4()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy();
        const auto lvecAllME = lvecAME + lvecNME1 + lvecNME2;
        const double dr2l1ME = reco::deltaR2(lvecAME.eta(),lvecAME.phi(),nonMEs.front()->eta(),nonMEs.front()->phi());
        const double dr2l2ME = reco::deltaR2(lvecAME.eta(),lvecAME.phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double dr2l1l2 = reco::deltaR2(nonMEs.front()->eta(),nonMEs.front()->phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double mlll = lvecAllME.M();

        const auto lvecAllMEpreCorr = lvecAMEpreCorr + lvecNME1 + lvecNME2;
        const auto mlllPreCorr = lvecAllMEpreCorr.M();

        if (mlll > 50.) {
          if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["3E_antiME_Et"]->Fill(lvecAME.Et(),aWeight);
            histo1d_["3E_antiME_eta"]->Fill(lvecAME.eta(),aWeight);
            histo1d_["3E_antiME_phi"]->Fill(lvecAME.phi(),aWeight);
          } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["3E_antiME_noGsf_Et"]->Fill(lvecAME.Et(),aWeight);
            histo1d_["3E_antiME_noGsf_eta"]->Fill(lvecAME.eta(),aWeight);
            histo1d_["3E_antiME_noGsf_phi"]->Fill(lvecAME.phi(),aWeight);
          } else {}

          histo1d_["3E_antiME_NME1_Et"]->Fill(nonMEs.front()->et()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy(),aWeight);
          histo1d_["3E_antiME_NME1_eta"]->Fill(nonMEs.front()->eta(),aWeight);
          histo1d_["3E_antiME_NME1_phi"]->Fill(nonMEs.front()->phi(),aWeight);

          histo1d_["3E_antiME_NME2_Et"]->Fill(nonMEs.at(1)->et()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy(),aWeight);
          histo1d_["3E_antiME_NME2_eta"]->Fill(nonMEs.at(1)->eta(),aWeight);
          histo1d_["3E_antiME_NME2_phi"]->Fill(nonMEs.at(1)->phi(),aWeight);

          histo1d_["3E_antiME_lll_pt"]->Fill(lvecAllME.pt(),aWeight);
          histo1d_["3E_antiME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
          histo1d_["3E_antiME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
          histo1d_["3E_antiME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

          double ffSS = 0., ffOS = 0.;
          auto ci = ciSSOS(antiMEs.front(),ffSS,ffOS);
          histo1d_["3E_antiME_lll_invM_CR"]->Fill(mlllPreCorr,aWeight);
          histo1d_["3E_antiME_lll_invM_CR_heepIdUp"]->Fill(mlllPreCorr,modHeepSFcl95(aWeight).first);
          histo1d_["3E_antiME_lll_invM_CR_heepIdDn"]->Fill(mlllPreCorr,modHeepSFcl95(aWeight).second);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF"]->Fill(mlll,aWeight*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_preCorr"]->Fill(mlllPreCorr,aWeight*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_heepIdUp"]->Fill(mlll,modHeepSFcl95(aWeight).first*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_heepIdDn"]->Fill(mlll,modHeepSFcl95(aWeight).second*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_elRecoUp"]->Fill(mlll,recoSFcl95(aWeight).first*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_elRecoDn"]->Fill(mlll,recoSFcl95(aWeight).second*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_elTrigUp"]->Fill(mlll,trigSFcl95(aWeight).first*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_elTrigDn"]->Fill(mlll,trigSFcl95(aWeight).second*ffOS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF"]->Fill(mlll,aWeight*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_preCorr"]->Fill(mlllPreCorr,aWeight*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_heepIdUp"]->Fill(mlll,modHeepSFcl95(aWeight).first*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_heepIdDn"]->Fill(mlll,modHeepSFcl95(aWeight).second*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_elRecoUp"]->Fill(mlll,recoSFcl95(aWeight).first*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_elRecoDn"]->Fill(mlll,recoSFcl95(aWeight).second*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_elTrigUp"]->Fill(mlll,trigSFcl95(aWeight).first*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_elTrigDn"]->Fill(mlll,trigSFcl95(aWeight).second*ffSS);
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_up"]->Fill(mlll,aWeight*(ffOS+ci.second));
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_up"]->Fill(mlll,aWeight*(ffSS+ci.first));
          histo1d_["3E_antiME_lll_invM_CR_xOSFF_dn"]->Fill(mlll,aWeight*(ffOS-ci.second));
          histo1d_["3E_antiME_lll_invM_CR_xSSFF_dn"]->Fill(mlll,aWeight*(ffSS-ci.first));

          if ( mlll < 200. ) {
            const double u5x5Et = estimateU5x5Et(antiMEs.front());
            const double u5x5Eta = estimateU5x5Eta(antiMEs.front());

            histo1d_["3E_antiME_Et_xSSFF"]->Fill(u5x5Et,aWeight*ffSS);
            histo1d_["3E_antiME_eta_xSSFF"]->Fill(u5x5Eta,aWeight*ffSS);
            histo1d_["3E_antiME_Et_xOSFF"]->Fill(u5x5Et,aWeight*ffOS);
            histo1d_["3E_antiME_eta_xOSFF"]->Fill(u5x5Eta,aWeight*ffOS);
            histo1d_["3E_antiME_Et_xSSFF_up"]->Fill(u5x5Et,aWeight*(ffSS+ci.first));
            histo1d_["3E_antiME_eta_xSSFF_up"]->Fill(u5x5Eta,aWeight*(ffSS+ci.first));
            histo1d_["3E_antiME_Et_xOSFF_up"]->Fill(u5x5Et,aWeight*(ffOS+ci.second));
            histo1d_["3E_antiME_eta_xOSFF_up"]->Fill(u5x5Eta,aWeight*(ffOS+ci.second));
            histo1d_["3E_antiME_Et_xSSFF_dn"]->Fill(u5x5Et,aWeight*(ffSS-ci.first));
            histo1d_["3E_antiME_eta_xSSFF_dn"]->Fill(u5x5Eta,aWeight*(ffSS-ci.first));
            histo1d_["3E_antiME_Et_xOSFF_dn"]->Fill(u5x5Et,aWeight*(ffOS-ci.second));
            histo1d_["3E_antiME_eta_xOSFF_dn"]->Fill(u5x5Eta,aWeight*(ffOS-ci.second));

            if ( antiMEs.front()->isEB() ) {
              histo1d_["3E_Et_CR_EB_antiME"]->Fill( u5x5Et, aWeight );
              histo1d_["3E_Et_CR_EB_antiME_xSSFF"]->Fill( u5x5Et, aWeight*ffSS);
              histo1d_["3E_Et_CR_EB_antiME_xSSFF_up"]->Fill( u5x5Et, aWeight*(ffSS+ci.first));
              histo1d_["3E_Et_CR_EB_antiME_xSSFF_dn"]->Fill( u5x5Et, aWeight*(ffSS-ci.first));
              histo1d_["3E_Et_CR_EB_antiME_xOSFF"]->Fill( u5x5Et, aWeight*ffOS);
              histo1d_["3E_Et_CR_EB_antiME_xOSFF_up"]->Fill( u5x5Et, aWeight*(ffOS+ci.second));
              histo1d_["3E_Et_CR_EB_antiME_xOSFF_dn"]->Fill( u5x5Et, aWeight*(ffOS-ci.second));
            }
          } // CR
        } // mlll > 50.
      } // antiMEs.size()==1

      break;
    case 2:
      histo1d_["cutflow_2E"]->Fill( 2.5, aWeight );
      if ( CRMEs.size() > 0 )
        histo1d_["cutflow_2E"]->Fill( 3.5, aWeight );
      if ( CRMEs.size()==2 )
        histo1d_["cutflow_2E"]->Fill( 4.5, aWeight );

      if (CRMEs.size()==2) {
        const auto lvecCRME1 = lvecME(CRMEs.at(0));
        const auto lvecCRME2 = lvecME(CRMEs.at(1));
        const auto lvecCRME1_scale = systHelperEle_.mergedEleScale(CRMEs.at(0))*lvecCRME1;
        const auto lvecCRME1_smear = systHelperEle_.mergedEleSmear(CRMEs.at(0),(*union5x5EnergyHandle)[CRMEs.at(0)])*lvecCRME1;
        const auto lvecCRME2_scale = systHelperEle_.mergedEleScale(CRMEs.at(1))*lvecCRME2;
        const auto lvecCRME2_smear = systHelperEle_.mergedEleSmear(CRMEs.at(1),(*union5x5EnergyHandle)[CRMEs.at(1)])*lvecCRME2;
        const auto lvecCRll = lvecCRME1 + lvecCRME2;
        const double dr2CRll = reco::deltaR2(lvecCRME1.eta(),lvecCRME1.phi(),lvecCRME2.eta(),lvecCRME2.phi());
        const double mll = lvecCRll.M();
        const double mll_scale = (lvecCRME1_scale+lvecCRME2_scale).M();
        const double mll_smear = (lvecCRME1_smear+lvecCRME2_smear).M();

        if ( mll > 50. && (mll < 200. || isMC_) ) {
          histo1d_["2E_CRME_ll_invM"]->Fill( mll, aWeight );

          SRinvM_ = mll;
          SRwgt_ = aWeight;
          SRpt1_ = lvecCRME1.pt();
          SReta1_ = lvecCRME1.eta();
          SRphi1_ = lvecCRME1.phi();
          SRpt2_ = lvecCRME2.pt();
          SReta2_ = lvecCRME2.eta();
          SRphi2_ = lvecCRME2.phi();

          if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() < 0. ) { // OS
            histo1d_["2E_CRME_OSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_CRME_OSll_invM_mergedEleScale"]->Fill( mll_scale, aWeight );
            histo1d_["2E_CRME_OSll_invM_mergedEleSmear"]->Fill( mll_smear, aWeight );
            histo1d_["2E_CRME_OSll_invM_heepIdUp"]->Fill( mll, modHeepSFcl95(aWeight).first );
            histo1d_["2E_CRME_OSll_invM_heepIdDn"]->Fill( mll, modHeepSFcl95(aWeight).second );
            histo1d_["2E_CRME_OSll_invM_elRecoUp"]->Fill( mll, recoSFcl95(aWeight).first );
            histo1d_["2E_CRME_OSll_invM_elRecoDn"]->Fill( mll, recoSFcl95(aWeight).second );
            histo1d_["2E_CRME_OSll_invM_elTrigUp"]->Fill( mll, trigSFcl95(aWeight).first );
            histo1d_["2E_CRME_OSll_invM_elTrigDn"]->Fill( mll, trigSFcl95(aWeight).second );
            histo1d_["2E_CRME_OSll_invM_mergedEleIdUp"]->Fill( mll, mergedEleSFcl95(aWeight).first );
            histo1d_["2E_CRME_OSll_invM_mergedEleIdDn"]->Fill( mll, mergedEleSFcl95(aWeight).second );

            SRisOS_ = 1;
          } else if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_CRME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_CRME_SSll_invM_mergedEleScale"]->Fill( mll_scale, aWeight );
            histo1d_["2E_CRME_SSll_invM_mergedEleSmear"]->Fill( mll_smear, aWeight );
            histo1d_["2E_CRME_SSll_invM_heepIdUp"]->Fill( mll, modHeepSFcl95(aWeight).first );
            histo1d_["2E_CRME_SSll_invM_heepIdDn"]->Fill( mll, modHeepSFcl95(aWeight).second );
            histo1d_["2E_CRME_SSll_invM_elRecoUp"]->Fill( mll, recoSFcl95(aWeight).first );
            histo1d_["2E_CRME_SSll_invM_elRecoDn"]->Fill( mll, recoSFcl95(aWeight).second );
            histo1d_["2E_CRME_SSll_invM_elTrigUp"]->Fill( mll, trigSFcl95(aWeight).first );
            histo1d_["2E_CRME_SSll_invM_elTrigDn"]->Fill( mll, trigSFcl95(aWeight).second );
            histo1d_["2E_CRME_SSll_invM_mergedEleIdUp"]->Fill( mll, mergedEleSFcl95(aWeight).first );
            histo1d_["2E_CRME_SSll_invM_mergedEleIdDn"]->Fill( mll, mergedEleSFcl95(aWeight).second );

            SRisOS_ = 0;
          }

          SRtree_->Fill();
        }

        if ( mll > 50. && mll < 200. ) {
          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["2E_CRME1_Et"]->Fill( lvecCRME1.Et(), aWeight );
            histo1d_["2E_CRME1_eta"]->Fill( lvecCRME1.eta(), aWeight );
            histo1d_["2E_CRME1_phi"]->Fill( lvecCRME1.phi(), aWeight );
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["2E_CRME1_noGsf_Et"]->Fill( lvecCRME1.Et(), aWeight );
            histo1d_["2E_CRME1_noGsf_eta"]->Fill( lvecCRME1.eta(), aWeight );
            histo1d_["2E_CRME1_noGsf_phi"]->Fill( lvecCRME1.phi(), aWeight );
          } else {}

          if ( CRMEs.at(1)->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["2E_CRME2_Et"]->Fill( lvecCRME2.Et(), aWeight );
            histo1d_["2E_CRME2_eta"]->Fill( lvecCRME2.eta(), aWeight );
            histo1d_["2E_CRME2_phi"]->Fill( lvecCRME2.phi(), aWeight );
          } else if ( CRMEs.at(1)->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["2E_CRME2_noGsf_Et"]->Fill( lvecCRME2.Et(), aWeight );
            histo1d_["2E_CRME2_noGsf_eta"]->Fill( lvecCRME2.eta(), aWeight );
            histo1d_["2E_CRME2_noGsf_phi"]->Fill( lvecCRME2.phi(), aWeight );
          } else {}

          histo1d_["2E_CRME_ll_pt"]->Fill( lvecCRll.pt(), aWeight );
          histo1d_["2E_CRME_ll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

          if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() < 0. ) { // OS
            histo1d_["2E_CRME_OSll_pt"]->Fill( lvecCRll.pt(), aWeight );
            histo1d_["2E_CRME_OSll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

            for (const auto& aME : CRMEs) {
              const double u5x5Et = estimateU5x5Et(aME);
              const double u5x5Eta = estimateU5x5Eta(aME);

              histo1d_["2E_CRME_OSll_Et"]->Fill( u5x5Et, aWeight );
              histo1d_["2E_CRME_OSll_eta"]->Fill( u5x5Eta, aWeight );

              if ( aME->isEB() ) {
                histo1d_["2E_Et_OSCR_EB_CRME"]->Fill( u5x5Et, aWeight );
                histo1d_["2E_eta_OSCR_EB_CRME"]->Fill( u5x5Eta, aWeight );
                histo2d_["2E_Et_eta_OSCR_EB_CRME"]->Fill( u5x5Eta, u5x5Et, aWeight );

                if ( std::abs(aME->superCluster()->eta()) < 0.8 )
                  histo1d_["2E_Et_OSCR_EB0p8_CRME"]->Fill( u5x5Et, aWeight );
                else
                  histo1d_["2E_Et_OSCR_EB1p5_CRME"]->Fill( u5x5Et, aWeight );
              }
            }
          } else if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_CRME_SSll_pt"]->Fill( lvecCRll.pt(), aWeight );
            histo1d_["2E_CRME_SSll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

            for (const auto& aME : CRMEs) {
              const double u5x5Et = estimateU5x5Et(aME);
              const double u5x5Eta = estimateU5x5Eta(aME);

              histo1d_["2E_CRME_SSll_Et"]->Fill( u5x5Et, aWeight );
              histo1d_["2E_CRME_SSll_eta"]->Fill( u5x5Eta, aWeight );

              if ( aME->isEB() )
                histo1d_["2E_Et_SSCR_EB_CRME"]->Fill( u5x5Et, aWeight );
            }
          }
        } // CR
      } // CRMEs.size()==2

      if (CRMEs.size()==1 && antiMEs.size()==1) {
        const auto lvecCRME = lvecME(CRMEs.front());
        const auto lvecAnMEpreCorr = lvecME(antiMEs.front());
        const auto lvecAnME = systHelperEle_.GetAbcdScaleSmear(lvecCRME,lvecAnMEpreCorr).second*lvecAnMEpreCorr;
        const auto lvecMixedll = lvecCRME + lvecAnME;
        const double dr2mixedll = reco::deltaR2(lvecCRME.eta(),lvecCRME.phi(),lvecAnME.eta(),lvecAnME.phi());
        const double mll = lvecMixedll.M();

        const auto lvecMixedllpreCorr = lvecCRME + lvecAnMEpreCorr;
        const auto mllPreCorr = lvecMixedllpreCorr.M();

        double ffSS = 0., ffOS = 0.;
        auto ci = ciSSOS(antiMEs.front(),ffSS,ffOS);

        if ( mll > 50. ) {
          histo1d_["2E_mixedME_ll_invM"]->Fill( mll, aWeight );

          mixedCRinvM_ = mllPreCorr;
          mixedCRwgt_ = aWeight;
          mixedCRpt1_ = lvecCRME.pt();
          mixedCReta1_ = lvecCRME.eta();
          mixedCRphi1_ = lvecCRME.phi();
          mixedCRpt2_ = lvecAnMEpreCorr.pt();
          mixedCReta2_ = lvecAnMEpreCorr.eta();
          mixedCRphi2_ = lvecAnMEpreCorr.phi();

          if ( CRMEs.front()->charge()*antiMEs.front()->charge() < 0. ) { // OS
            histo1d_["2E_mixedME_OSll_invM"]->Fill( mllPreCorr, aWeight );
            histo1d_["2E_mixedME_OSll_invM_heepIdUp"]->Fill( mllPreCorr, modHeepSFcl95(aWeight).first );
            histo1d_["2E_mixedME_OSll_invM_heepIdDn"]->Fill( mllPreCorr, modHeepSFcl95(aWeight).second );
            histo1d_["2E_mixedME_OSll_invM_mergedEleIdUp"]->Fill( mllPreCorr, mergedEleSFcl95(aWeight).first );
            histo1d_["2E_mixedME_OSll_invM_mergedEleIdDn"]->Fill( mllPreCorr, mergedEleSFcl95(aWeight).second );
            histo1d_["2E_mixedME_OSll_invM_xFF"]->Fill( mll, 0.5*aWeight*ffOS ); // double-counting
            histo1d_["2E_mixedME_OSll_invM_xFF_preCorr"]->Fill( mllPreCorr, 0.5*aWeight*ffOS );
            histo1d_["2E_mixedME_OSll_invM_xFF_up"]->Fill( mll, 0.5*aWeight*(ffOS+ci.second) );
            histo1d_["2E_mixedME_OSll_invM_xFF_dn"]->Fill( mll, 0.5*aWeight*(ffOS-ci.second) );

            mixedCRwgt_ *= (0.5*ffOS);
            mixedCRisOS_ = 1;
          } else if ( CRMEs.front()->charge()*antiMEs.front()->charge() > 0. ) { // SS
            histo1d_["2E_mixedME_SSll_invM"]->Fill( mllPreCorr, aWeight );
            histo1d_["2E_mixedME_SSll_invM_heepIdUp"]->Fill( mllPreCorr, modHeepSFcl95(aWeight).first );
            histo1d_["2E_mixedME_SSll_invM_heepIdDn"]->Fill( mllPreCorr, modHeepSFcl95(aWeight).second );
            histo1d_["2E_mixedME_SSll_invM_mergedEleIdUp"]->Fill( mllPreCorr, mergedEleSFcl95(aWeight).first );
            histo1d_["2E_mixedME_SSll_invM_mergedEleIdDn"]->Fill( mllPreCorr, mergedEleSFcl95(aWeight).second );
            histo1d_["2E_mixedME_SSll_invM_xFF"]->Fill( mll, 0.5*aWeight*ffSS ); // double-counting
            histo1d_["2E_mixedME_SSll_invM_xFF_preCorr"]->Fill( mllPreCorr, 0.5*aWeight*ffSS );
            histo1d_["2E_mixedME_SSll_invM_xFF_up"]->Fill( mll, 0.5*aWeight*(ffSS+ci.first) );
            histo1d_["2E_mixedME_SSll_invM_xFF_dn"]->Fill( mll, 0.5*aWeight*(ffSS-ci.first) );

            mixedCRwgt_ *= (0.5*ffSS);
            mixedCRisOS_ = 0;
          }

          mixedCRtree_->Fill();
        }

        if ( mll > 50. && mll < 200. ) {
          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["2E_mixed_ME_Et"]->Fill( lvecCRME.Et(), aWeight );
            histo1d_["2E_mixed_ME_eta"]->Fill( lvecCRME.eta(), aWeight );
            histo1d_["2E_mixed_ME_phi"]->Fill( lvecCRME.phi(), aWeight );
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["2E_mixed_ME_noGsf_Et"]->Fill( lvecCRME.Et(), aWeight );
            histo1d_["2E_mixed_ME_noGsf_eta"]->Fill( lvecCRME.eta(), aWeight );
            histo1d_["2E_mixed_ME_noGsf_phi"]->Fill( lvecCRME.phi(), aWeight );
          } else {}

          if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["2E_mixed_antiME_Et"]->Fill( lvecAnME.Et(), aWeight );
            histo1d_["2E_mixed_antiME_eta"]->Fill( lvecAnME.eta(), aWeight );
            histo1d_["2E_mixed_antiME_phi"]->Fill( lvecAnME.phi(), aWeight );
          } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["2E_mixed_antiME_noGsf_Et"]->Fill( lvecAnME.Et(), aWeight );
            histo1d_["2E_mixed_antiME_noGsf_eta"]->Fill( lvecAnME.eta(), aWeight );
            histo1d_["2E_mixed_antiME_noGsf_phi"]->Fill( lvecAnME.phi(), aWeight );
          } else {}

          histo1d_["2E_mixedME_ll_pt"]->Fill( lvecMixedll.pt(), aWeight );
          histo1d_["2E_mixedME_ll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

          if ( CRMEs.front()->charge()*antiMEs.front()->charge() < 0. ) { // OS
            histo1d_["2E_mixedME_OSll_pt"]->Fill( lvecMixedll.pt(), aWeight );
            histo1d_["2E_mixedME_OSll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

            const double u5x5EtPass = estimateU5x5Et(CRMEs.front());
            const double u5x5EtaPass = estimateU5x5Eta(CRMEs.front());
            const double u5x5EtFail = estimateU5x5Et(antiMEs.front());
            const double u5x5EtaFail = estimateU5x5Eta(antiMEs.front());

            histo1d_["2E_mixedME_OSll_Et"]->Fill( u5x5EtPass, aWeight );
            histo1d_["2E_mixedME_OSll_Et_noCorr"]->Fill( CRMEs.front()->et(), aWeight );
            histo1d_["2E_mixedME_OSll_eta"]->Fill( u5x5EtaPass, aWeight );
            histo1d_["2E_mixedAntiME_OSll_Et_xFF"]->Fill( u5x5EtFail, aWeight*ffOS );
            histo1d_["2E_mixedAntiME_OSll_eta_xFF"]->Fill( u5x5EtaFail, aWeight*ffOS );
            histo1d_["2E_mixedAntiME_OSll_Et_xFF_up"]->Fill( u5x5EtFail, aWeight*(ffOS+ci.second) );
            histo1d_["2E_mixedAntiME_OSll_eta_xFF_up"]->Fill( u5x5EtaFail, aWeight*(ffOS+ci.second) );
            histo1d_["2E_mixedAntiME_OSll_Et_xFF_dn"]->Fill( u5x5EtFail, aWeight*(ffOS-ci.second) );
            histo1d_["2E_mixedAntiME_OSll_eta_xFF_dn"]->Fill( u5x5EtaFail, aWeight*(ffOS-ci.second) );

            if ( CRMEs.front()->isEB() ) {
              histo1d_["2E_Et_OSCR_EB_mixedME"]->Fill( u5x5EtPass, aWeight );
              histo1d_["2E_eta_OSCR_EB_mixedME"]->Fill( u5x5EtaPass, aWeight );
              histo2d_["2E_Et_eta_OSCR_EB_mixedME"]->Fill( u5x5EtaPass, u5x5EtPass, aWeight );

              if ( std::abs(CRMEs.front()->superCluster()->eta()) < 0.8 )
                histo1d_["2E_Et_OSCR_EB0p8_mixedME"]->Fill( u5x5EtPass, aWeight );
              else
                histo1d_["2E_Et_OSCR_EB1p5_mixedME"]->Fill( u5x5EtPass, aWeight );
            }

            if ( antiMEs.front()->isEB() ) {
              histo1d_["2E_Et_OSCR_EB_mixedAntiME"]->Fill( u5x5EtFail, aWeight );
              histo1d_["2E_Et_OSCR_EB_mixedAntiME_xFF"]->Fill(u5x5EtFail, aWeight*ffOS );
              histo1d_["2E_Et_OSCR_EB_mixedAntiME_xFF_up"]->Fill(u5x5EtFail, aWeight*(ffOS+ci.second) );
              histo1d_["2E_Et_OSCR_EB_mixedAntiME_xFF_dn"]->Fill(u5x5EtFail, aWeight*(ffOS-ci.second) );
              histo1d_["2E_eta_OSCR_EB_mixedAntiME"]->Fill( u5x5EtaFail, aWeight );
              histo1d_["2E_eta_OSCR_EB_mixedAntiME_xFF"]->Fill( u5x5EtaFail, aWeight*ffOS );
              histo2d_["2E_Et_eta_OSCR_EB_mixedAntiME"]->Fill( u5x5EtaFail, u5x5EtFail, aWeight );
              histo2d_["2E_Et_eta_OSCR_EB_mixedAntiME_xFF"]->Fill( u5x5EtaFail, u5x5EtFail, aWeight*ffOS );

              if ( std::abs(antiMEs.front()->superCluster()->eta()) < 0.8 ) {
                histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME"]->Fill( u5x5EtFail, aWeight );
                histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME_xFF"]->Fill( u5x5EtFail, aWeight*ffOS );
                histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME_xFF_up"]->Fill( u5x5EtFail, aWeight*(ffOS+ci.second) );
                histo1d_["2E_Et_OSCR_EB0p8_mixedAntiME_xFF_dn"]->Fill( u5x5EtFail, aWeight*(ffOS-ci.second) );
              } else {
                histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME"]->Fill( u5x5EtFail, aWeight );
                histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME_xFF"]->Fill( u5x5EtFail, aWeight*ffOS );
                histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME_xFF_up"]->Fill( u5x5EtFail, aWeight*(ffOS+ci.second) );
                histo1d_["2E_Et_OSCR_EB1p5_mixedAntiME_xFF_dn"]->Fill( u5x5EtFail, aWeight*(ffOS-ci.second) );
              }
            }
          } else if ( CRMEs.front()->charge()*antiMEs.front()->charge() > 0. ) { // SS
            histo1d_["2E_mixedME_SSll_pt"]->Fill( lvecMixedll.pt(), aWeight );
            histo1d_["2E_mixedME_SSll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

            const double u5x5EtPass = estimateU5x5Et(CRMEs.front());
            const double u5x5EtaPass = estimateU5x5Eta(CRMEs.front());
            const double u5x5EtFail = estimateU5x5Et(antiMEs.front());
            const double u5x5EtaFail = estimateU5x5Eta(antiMEs.front());

            histo1d_["2E_mixedME_SSll_Et"]->Fill( u5x5EtPass, aWeight );
            histo1d_["2E_mixedME_SSll_Et_noCorr"]->Fill( CRMEs.front()->et(), aWeight );
            histo1d_["2E_mixedME_SSll_eta"]->Fill( u5x5EtaPass, aWeight );
            histo1d_["2E_mixedAntiME_SSll_Et_xFF"]->Fill( u5x5EtFail, aWeight*ffSS );
            histo1d_["2E_mixedAntiME_SSll_eta_xFF"]->Fill( u5x5EtaFail, aWeight*ffSS );
            histo1d_["2E_mixedAntiME_SSll_Et_xFF_up"]->Fill( u5x5EtFail, aWeight*(ffSS+ci.first) );
            histo1d_["2E_mixedAntiME_SSll_eta_xFF_up"]->Fill( u5x5EtaFail, aWeight*(ffSS+ci.first) );
            histo1d_["2E_mixedAntiME_SSll_Et_xFF_dn"]->Fill( u5x5EtFail, aWeight*(ffSS-ci.first) );
            histo1d_["2E_mixedAntiME_SSll_eta_xFF_dn"]->Fill( u5x5EtaFail, aWeight*(ffSS-ci.first) );

            if ( CRMEs.front()->isEB() )
              histo1d_["2E_Et_SSCR_EB_mixedME"]->Fill( u5x5EtPass, aWeight );

            if ( antiMEs.front()->isEB() ) {
              histo1d_["2E_Et_SSCR_EB_mixedAntiME"]->Fill( u5x5EtFail, aWeight );
              histo1d_["2E_Et_SSCR_EB_mixedAntiME_xFF"]->Fill( u5x5EtFail, aWeight*ffSS );
              histo1d_["2E_Et_SSCR_EB_mixedAntiME_xFF_up"]->Fill( u5x5EtFail, aWeight*(ffSS+ci.first) );
              histo1d_["2E_Et_SSCR_EB_mixedAntiME_xFF_dn"]->Fill( u5x5EtFail, aWeight*(ffSS-ci.first) );
            }
          }
        } // CR
      } // CRMEs.size()==1 && antiMEs.size()==1

      if (antiMEs.size()==2) {
        const auto lvecAME1preCorr = lvecME(antiMEs.at(0));
        const auto lvecAME2preCorr = lvecME(antiMEs.at(1));
        const auto lvecAME1 = systHelperEle_.GetAbcdScaleSmear(lvecAME1preCorr,lvecAME2preCorr).first*lvecAME1preCorr;
        const auto lvecAME2 = systHelperEle_.GetAbcdScaleSmear(lvecAME1preCorr,lvecAME2preCorr).second*lvecAME2preCorr;
        const auto lvecAll1 = lvecAME1 + lvecAME2preCorr;
        const auto lvecAll2 = lvecAME1preCorr + lvecAME2;
        const auto lvecAllPreCorr = lvecAME1preCorr + lvecAME2preCorr;
        const double dr2All = reco::deltaR2(lvecAME1.eta(),lvecAME1.phi(),lvecAME2.eta(),lvecAME2.phi());
        const double mll1 = lvecAll1.M();
        const double mll2 = lvecAll2.M();
        const double mllPreCorr = lvecAllPreCorr.M();

        if ( mllPreCorr > 50. ) {
          if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["2E_antiME1_Et"]->Fill( lvecAME1.Et(), aWeight );
            histo1d_["2E_antiME1_eta"]->Fill( lvecAME1.eta(), aWeight );
            histo1d_["2E_antiME1_phi"]->Fill( lvecAME1.phi(), aWeight );
          } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["2E_antiME1_noGsf_Et"]->Fill( lvecAME1.Et(), aWeight );
            histo1d_["2E_antiME1_noGsf_eta"]->Fill( lvecAME1.eta(), aWeight );
            histo1d_["2E_antiME1_noGsf_phi"]->Fill( lvecAME1.phi(), aWeight );
          } else {}

          if ( antiMEs.at(1)->userInt("mvaMergedElectronCategories")!=1 ) {
            histo1d_["2E_antiME2_Et"]->Fill( lvecAME2.Et(), aWeight );
            histo1d_["2E_antiME2_eta"]->Fill( lvecAME2.eta(), aWeight );
            histo1d_["2E_antiME2_phi"]->Fill( lvecAME2.phi(), aWeight );
          } else if ( antiMEs.at(1)->userInt("mvaMergedElectronCategories")==1 ) {
            histo1d_["2E_antiME2_noGsf_Et"]->Fill( lvecAME2.Et(), aWeight );
            histo1d_["2E_antiME2_noGsf_eta"]->Fill( lvecAME2.eta(), aWeight );
            histo1d_["2E_antiME2_noGsf_phi"]->Fill( lvecAME2.phi(), aWeight );
          } else {}

          histo1d_["2E_antiME_ll_invM"]->Fill( mllPreCorr, aWeight );
          histo1d_["2E_antiME_ll_pt"]->Fill( lvecAllPreCorr.pt(), aWeight );
          histo1d_["2E_antiME_ll_dr"]->Fill( std::sqrt(dr2All), aWeight );

          double ffSS1 = 0., ffOS1 = 0., ffSS2 = 0., ffOS2 = 0.;
          auto ci1 = ciSSOS(antiMEs.at(0),ffSS1,ffOS1);
          auto ci2 = ciSSOS(antiMEs.at(1),ffSS2,ffOS2);

          antiCRinvM_ = mllPreCorr;
          antiCRwgt_ = aWeight;
          antiCRpt1_ = lvecAME1preCorr.pt();
          antiCReta1_ = lvecAME1preCorr.eta();
          antiCRphi1_ = lvecAME1preCorr.phi();
          antiCRpt2_ = lvecAME2preCorr.pt();
          antiCReta2_ = lvecAME2preCorr.eta();
          antiCRphi2_ = lvecAME2preCorr.phi();

          if ( antiMEs.front()->charge()*antiMEs.at(1)->charge() < 0. ) { // OS
            histo1d_["2E_antiME_OSll_invM_CR"]->Fill( mllPreCorr, aWeight );
            histo1d_["2E_antiME_OSll_invM_CR_xFF"]->Fill( mll1, aWeight*ffOS1 );
            histo1d_["2E_antiME_OSll_invM_CR_xFF"]->Fill( mll2, aWeight*ffOS2 );
            histo1d_["2E_antiME_OSll_invM_CR_xFF_preCorr"]->Fill( mllPreCorr, aWeight*ffOS1 );
            histo1d_["2E_antiME_OSll_invM_CR_xFF_preCorr"]->Fill( mllPreCorr, aWeight*ffOS2 );
            histo1d_["2E_antiME_OSll_invM_CR_xFF_up"]->Fill( mll1, aWeight*(ffOS1+ci1.second) );
            histo1d_["2E_antiME_OSll_invM_CR_xFF_up"]->Fill( mll2, aWeight*(ffOS2+ci2.second) );
            histo1d_["2E_antiME_OSll_invM_CR_xFF_dn"]->Fill( mll1, aWeight*(ffOS1-ci1.second) );
            histo1d_["2E_antiME_OSll_invM_CR_xFF_dn"]->Fill( mll2, aWeight*(ffOS2-ci2.second) );

            histo1d_["2E_antiME_OSll_pt"]->Fill( lvecAllPreCorr.pt(), aWeight );
            histo1d_["2E_antiME_OSll_dr"]->Fill( std::sqrt(dr2All), aWeight );

            antiCRwgt_ *= (ffOS1+ffOS2);
            antiCRisOS_ = 1;

            if (mllPreCorr < 200.) {
              for (const auto& aME : antiMEs) {
                double ffSS = 0., ffOS = 0.;
                auto aci = ciSSOS(aME,ffSS,ffOS);

                const double u5x5Et = estimateU5x5Et(aME);
                const double u5x5Eta = estimateU5x5Eta(aME);

                histo1d_["2E_antiME_OSll_Et_noCorr"]->Fill( aME->et(), aWeight );
                histo1d_["2E_antiME_OSll_eta"]->Fill( u5x5Eta, aWeight );

                histo1d_["2E_antiME_OSll_Et_xFF"]->Fill( u5x5Et, aWeight*ffOS );
                histo1d_["2E_antiME_OSll_eta_xFF"]->Fill( u5x5Eta, aWeight*ffOS );
                histo1d_["2E_antiME_OSll_Et_xFF_up"]->Fill( u5x5Et, aWeight*(ffOS+aci.second) );
                histo1d_["2E_antiME_OSll_eta_xFF_up"]->Fill( u5x5Eta, aWeight*(ffOS+aci.second) );
                histo1d_["2E_antiME_OSll_Et_xFF_dn"]->Fill( u5x5Et, aWeight*(ffOS-aci.second) );
                histo1d_["2E_antiME_OSll_eta_xFF_dn"]->Fill( u5x5Eta, aWeight*(ffOS-aci.second) );

                if ( aME->isEB() ) {
                  histo1d_["2E_Et_OSCR_EB_antiME"]->Fill( u5x5Et, aWeight );
                  histo1d_["2E_Et_OSCR_EB_antiME_xFF"]->Fill( u5x5Et, aWeight*ffOS );
                  histo1d_["2E_Et_OSCR_EB_antiME_xFF_up"]->Fill( u5x5Et, aWeight*(ffOS+aci.second) );
                  histo1d_["2E_Et_OSCR_EB_antiME_xFF_dn"]->Fill( u5x5Et, aWeight*(ffOS-aci.second) );
                  histo1d_["2E_eta_OSCR_EB_antiME"]->Fill( u5x5Eta, aWeight );
                  histo1d_["2E_eta_OSCR_EB_antiME_xFF"]->Fill( u5x5Eta, aWeight*ffOS );
                  histo2d_["2E_Et_eta_OSCR_EB_antiME"]->Fill( u5x5Eta, u5x5Et, aWeight );
                  histo2d_["2E_Et_eta_OSCR_EB_antiME_xFF"]->Fill( u5x5Eta, u5x5Et, aWeight*ffOS );

                  if ( std::abs(aME->superCluster()->eta()) < 0.8 ) {
                    histo1d_["2E_Et_OSCR_EB0p8_antiME"]->Fill( u5x5Et, aWeight );
                    histo1d_["2E_Et_OSCR_EB0p8_antiME_xFF"]->Fill( u5x5Et, aWeight*ffOS );
                    histo1d_["2E_Et_OSCR_EB0p8_antiME_xFF_up"]->Fill( u5x5Et, aWeight*(ffOS+aci.second) );
                    histo1d_["2E_Et_OSCR_EB0p8_antiME_xFF_dn"]->Fill( u5x5Et, aWeight*(ffOS-aci.second) );
                  } else {
                    histo1d_["2E_Et_OSCR_EB1p5_antiME"]->Fill( u5x5Et, aWeight );
                    histo1d_["2E_Et_OSCR_EB1p5_antiME_xFF"]->Fill( u5x5Et, aWeight*ffOS );
                    histo1d_["2E_Et_OSCR_EB1p5_antiME_xFF_up"]->Fill( u5x5Et, aWeight*(ffOS+aci.second) );
                    histo1d_["2E_Et_OSCR_EB1p5_antiME_xFF_dn"]->Fill( u5x5Et, aWeight*(ffOS-aci.second) );
                  }
                }
              }
            }
          } else if ( antiMEs.front()->charge()*antiMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_antiME_SSll_invM_CR"]->Fill( mllPreCorr, aWeight );
            histo1d_["2E_antiME_SSll_invM_CR_xFF"]->Fill( mll1, aWeight*ffSS1 );
            histo1d_["2E_antiME_SSll_invM_CR_xFF"]->Fill( mll2, aWeight*ffSS2 );
            histo1d_["2E_antiME_SSll_invM_CR_xFF_preCorr"]->Fill( mllPreCorr, aWeight*ffSS1 );
            histo1d_["2E_antiME_SSll_invM_CR_xFF_preCorr"]->Fill( mllPreCorr, aWeight*ffSS2 );
            histo1d_["2E_antiME_SSll_invM_CR_xFF_up"]->Fill( mll1, aWeight*(ffSS1+ci1.first) );
            histo1d_["2E_antiME_SSll_invM_CR_xFF_up"]->Fill( mll2, aWeight*(ffSS2+ci2.first) );
            histo1d_["2E_antiME_SSll_invM_CR_xFF_dn"]->Fill( mll1, aWeight*(ffSS1-ci1.first) );
            histo1d_["2E_antiME_SSll_invM_CR_xFF_dn"]->Fill( mll2, aWeight*(ffSS2-ci2.first) );

            histo1d_["2E_antiME_SSll_pt"]->Fill( lvecAllPreCorr.pt(), aWeight );
            histo1d_["2E_antiME_SSll_dr"]->Fill( std::sqrt(dr2All), aWeight );

            antiCRwgt_ *= (ffSS1+ffSS2);
            antiCRisOS_ = 0;

            if (mllPreCorr < 200.) {
              for (const auto& aME : antiMEs) {
                double ffSS = 0., ffOS = 0.;
                auto aci = ciSSOS(aME,ffSS,ffOS);

                const double u5x5Et = estimateU5x5Et(aME);
                const double u5x5Eta = estimateU5x5Eta(aME);

                histo1d_["2E_antiME_SSll_Et_noCorr"]->Fill( aME->et(), aWeight );
                histo1d_["2E_antiME_SSll_eta"]->Fill( u5x5Eta, aWeight );

                histo1d_["2E_antiME_SSll_Et_xFF"]->Fill( u5x5Et, aWeight*ffSS );
                histo1d_["2E_antiME_SSll_eta_xFF"]->Fill( u5x5Eta, aWeight*ffSS );
                histo1d_["2E_antiME_SSll_Et_xFF_up"]->Fill( u5x5Et, aWeight*(ffSS+aci.first) );
                histo1d_["2E_antiME_SSll_eta_xFF_up"]->Fill( u5x5Eta, aWeight*(ffSS+aci.first) );
                histo1d_["2E_antiME_SSll_Et_xFF_dn"]->Fill( u5x5Et, aWeight*(ffSS-aci.first) );
                histo1d_["2E_antiME_SSll_eta_xFF_dn"]->Fill( u5x5Eta, aWeight*(ffSS-aci.first) );

                if ( aME->isEB() ) {
                  histo1d_["2E_Et_SSCR_EB_antiME"]->Fill( u5x5Et, aWeight );
                  histo1d_["2E_Et_SSCR_EB_antiME_xFF"]->Fill( u5x5Et, aWeight*ffSS );
                  histo1d_["2E_Et_SSCR_EB_antiME_xFF_up"]->Fill( u5x5Et, aWeight*(ffSS+aci.first) );
                  histo1d_["2E_Et_SSCR_EB_antiME_xFF_dn"]->Fill( u5x5Et, aWeight*(ffSS-aci.first) );
                }
              }
            } // CR
          } // OSSS

          antiCRtree_->Fill();
        } // mll > 50.
      } // antiME.size()==2

      break;
  } // switch acceptEles.size()
}

DEFINE_FWK_MODULE(MergedEleCRanalyzer);
