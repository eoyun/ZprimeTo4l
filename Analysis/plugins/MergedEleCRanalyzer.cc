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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "TH2D.h"
#include "TF1.h"
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
  const bool isDY_;
  const bool selectHT_;
  const bool selectZpT_;
  const double maxHT_;
  const double maxVpT_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;

  const edm::FileInPath recoSFpath_;
  const edm::FileInPath recoSFlowPtPath_;
  const edm::FileInPath FFpath_;
  const std::vector<std::string> trigList_;
  const double etThres1_;
  const double ssBoundary_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  std::unique_ptr<TFile> recoSFfile_;
  std::unique_ptr<TFile> recoSFlowPtFile_;
  TH2D* recoSF_;
  TH2D* recoSFlowPt_;

  std::unique_ptr<TFile> FFfile_;
  TF1* ssboth_;
  TF1* osboth_;
  TFitResultPtr ssFit_;
  TFitResultPtr osFit_;
};

MergedEleCRanalyzer::MergedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
isDY_(iConfig.getUntrackedParameter<bool>("isDY")),
selectHT_(iConfig.getUntrackedParameter<bool>("selectHT")),
selectZpT_(iConfig.getUntrackedParameter<bool>("selectZpT")),
maxHT_(iConfig.getParameter<double>("maxHT")),
maxVpT_(iConfig.getParameter<double>("maxVpT")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
recoSFpath_(iConfig.getParameter<edm::FileInPath>("recoSFpath")),
recoSFlowPtPath_(iConfig.getParameter<edm::FileInPath>("recoSFlowPtPath")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
etThres1_(iConfig.getParameter<double>("etThres1")),
ssBoundary_(iConfig.getParameter<double>("ssBoundary")) {
  usesResource("TFileService");
}

void MergedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  recoSFfile_ = std::make_unique<TFile>(recoSFpath_.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));

  recoSFlowPtFile_ = std::make_unique<TFile>(recoSFlowPtPath_.fullPath().c_str(),"READ");
  recoSFlowPt_ = static_cast<TH2D*>(recoSFlowPtFile_->Get("EGamma_SF2D"));

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  ssboth_ = static_cast<TF1*>(FFfile_->FindObjectAny("ssboth"));
  osboth_ = static_cast<TF1*>(FFfile_->FindObjectAny("osboth"));
  ssFit_ = (static_cast<TH1D*>(FFfile_->Get("2E_mixedME_SSll_Et_noCorr_rebin")))->Fit(ssboth_,"RS");
  osFit_ = (static_cast<TH1D*>(FFfile_->Get("2E_mixedME_OSll_Et_noCorr_rebin")))->Fit(osboth_,"RS");

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["lheHT"] = fs->make<TH1D>("lheHT","lheHT",600,0.,3000.);
  histo1d_["lheHT_cut"] = fs->make<TH1D>("lheHT_cut","lheHT_cut",600,0.,3000.);
  histo1d_["lheVpT"] = fs->make<TH1D>("lheVpT","lheVpT",500,0.,1000.);
  histo1d_["lheVpT_cut"] = fs->make<TH1D>("lheVpT_cut","lheVpT",500,0.,1000.);
  histo1d_["N-1_modifiedHEEP_EB"] = fs->make<TH1D>("N-1_modifiedHEEP_EB","N-1 efficiency (modified HEEP EB)",16,-1.,15.);
  histo1d_["N-1_modifiedHEEP_EE"] = fs->make<TH1D>("N-1_modifiedHEEP_EE","N-1 efficiency (modified HEEP EE)",16,-1.,15.);
  histo1d_["N-1_HEEP_EB"] = fs->make<TH1D>("N-1_HEEP_EB","N-1 efficiency (HEEP EB)",16,-1.,15.);
  histo1d_["N-1_HEEP_EE"] = fs->make<TH1D>("N-1_HEEP_EE","N-1 efficiency (HEEP EE)",16,-1.,15.);
  histo1d_["mva_DR1Et2EB"] = fs->make<TH1D>("mva_DR1Et2EB","MVA score",100,0.,1.);
  histo1d_["mva_DR2Et1EB"] = fs->make<TH1D>("mva_DR2Et1EB","MVA score",100,0.,1.);
  histo1d_["mva_DR2Et2EB"] = fs->make<TH1D>("mva_DR2Et2EB","MVA score",100,0.,1.);
  histo1d_["mva_bkgEt2EB"] = fs->make<TH1D>("mva_bkgEt2EB","MVA score",100,0.,1.);
  histo1d_["mva_DR1Et1EB_extended"] = fs->make<TH1D>("mva_DR1Et1EB_extended","MVA score",100,0.,1.);
  histo1d_["mva_bkgEt1EB_extended"] = fs->make<TH1D>("mva_bkgEt1EB_extended","MVA score",100,0.,1.);
  histo1d_["cutflow_4E"] = fs->make<TH1D>("cutflow_4E","cutflow (4E)",10,0.,10.);
  histo1d_["cutflow_3E"] = fs->make<TH1D>("cutflow_3E","cutflow (3E)",10,0.,10.);
  histo1d_["cutflow_2E"] = fs->make<TH1D>("cutflow_2E","cutflow (2E)",10,0.,10.);

  histo2d_["invM_dr"] = fs->make<TH2D>("invM_dr","M(ll) vs dR(ll)",200,0.,1000.,128,0.,6.4);
  histo2d_["invM_dr_noGsf"] = fs->make<TH2D>("invM_dr_noGsf","M(ll) vs dR(ll)",200,0.,1000.,128,0.,6.4);
  histo2d_["invM_dr_mixed"] = fs->make<TH2D>("invM_dr_mixed","M(ll) vs dR(ll)",200,0.,1000.,128,0.,6.4);

  histo1d_["invM_ll_CRME_EB"] = fs->make<TH1D>("invM_ll_CRME_EB","ME M(ll);[GeV];",500,0.,5.);
  histo1d_["invM_ll_CRME_Et200_EB"] = fs->make<TH1D>("invM_ll_CRME_Et200_EB","ME M(ll);[GeV];",500,0.,5.);

  histo2d_["et_dr_final"] = fs->make<TH2D>("et_dr_final","Et vs dR(l,GSF)",200,0.,2000.,100,0.,0.5);
  histo2d_["et_dr_HEEP"] = fs->make<TH2D>("et_dr_HEEP","Et vs dR(l,GSF)",200,0.,2000.,100,0.,0.5);

  histo1d_["et_noGsf_final"] = fs->make<TH1D>("et_noGsf_final","Et (CR);E_{T} [GeV];",200,0.,2000.);
  histo1d_["et_noGsf_HEEP"] = fs->make<TH1D>("et_noGsf_HEEP","Et (pass HEEP);E_{T} [GeV];",200,0.,2000.);

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

  histo1d_["2E_CRME_ll_invM"] = fs->make<TH1D>("2E_CRME_ll_invM","CR M(ll);GeV;",200,0.,200.);
  histo1d_["2E_CRME_ll_pt"] = fs->make<TH1D>("2E_CRME_ll_pt","CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_ll_dr"] = fs->make<TH1D>("2E_CRME_ll_dr","CR dR(ll)",64,0.,6.4);

  histo1d_["2E_CRME_SSll_invM"] = fs->make<TH1D>("2E_CRME_SSll_invM","SS CR M(ll);GeV;",200,0.,200.);
  histo1d_["2E_CRME_SSll_pt"] = fs->make<TH1D>("2E_CRME_SSll_pt","SS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_SSll_dr"] = fs->make<TH1D>("2E_CRME_SSll_dr","SS CR dR(ll)",64,0.,6.4);
  histo1d_["2E_CRME_SSll_Et"] = fs->make<TH1D>("2E_CRME_SSll_Et","SS CR ME E_{T};GeV;",100,0.,500.);
  histo1d_["2E_CRME_SSll_eta"] = fs->make<TH1D>("2E_CRME_SSll_eta","SS CR ME #eta",50,-2.5,2.5);

  histo1d_["2E_CRME_OSll_invM"] = fs->make<TH1D>("2E_CRME_OSll_invM","OS CR M(ll);GeV;",200,0.,200.);
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

  histo1d_["2E_mixedME_ll_invM"] = fs->make<TH1D>("2E_mixedME_ll_invM","mixed CR M(ll);GeV;",200,0.,200.);
  histo1d_["2E_mixedME_ll_pt"] = fs->make<TH1D>("2E_mixedME_ll_pt","mixed CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_ll_dr"] = fs->make<TH1D>("2E_mixedME_ll_dr","mixed CR dR(ll)",64,0.,6.4);

  histo1d_["2E_mixedME_SSll_invM"] = fs->make<TH1D>("2E_mixedME_SSll_invM","mixed SS CR M(ll);GeV;",200,0.,200.);
  histo1d_["2E_mixedME_SSll_invM_xFF"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF","mixed SS CR M(ll) x Fake factor;GeV;",200,0.,200.);
  histo1d_["2E_mixedME_SSll_invM_xFF_up"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF_up","mixed SS CR M(ll) x Fake factor (up);GeV;",200,0.,200.);
  histo1d_["2E_mixedME_SSll_invM_xFF_dn"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF_dn","mixed SS CR M(ll) x Fake factor (down);GeV;",200,0.,200.);
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

  histo1d_["2E_mixedME_OSll_invM"] = fs->make<TH1D>("2E_mixedME_OSll_invM","mixed OS CR M(ll);GeV;",200,0.,200.);
  histo1d_["2E_mixedME_OSll_invM_xFF"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF","mixed OS CR M(ll) x Fake factor;GeV;",200,0.,200.);
  histo1d_["2E_mixedME_OSll_invM_xFF_up"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF_up","mixed OS CR M(ll) x Fake factor (up);GeV;",200,0.,200.);
  histo1d_["2E_mixedME_OSll_invM_xFF_dn"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF_dn","mixed OS CR M(ll) x Fake factor (down);GeV;",200,0.,200.);
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

  histo1d_["2E_antiME_ll_invM"] = fs->make<TH1D>("2E_antiME_ll_invM","anti M(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_ll_invM_CR"] = fs->make<TH1D>("2E_antiME_ll_invM_CR","anti M(ll);GeV;",200,0.,200.);
  histo1d_["2E_antiME_ll_pt"] = fs->make<TH1D>("2E_antiME_ll_pt","anti p_{T}(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_ll_dr"] = fs->make<TH1D>("2E_antiME_ll_dr","anti dR(ll)",128,0.,6.4);

  histo1d_["2E_antiME_SSll_invM"] = fs->make<TH1D>("2E_antiME_SSll_invM","SS anti M(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_invM_CR"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR","SS anti M(ll);GeV;",200,0.,200.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF","SS anti M(ll) x Fake factor;GeV;",200,0.,200.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF_up"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF_up","SS anti M(ll) x Fake factor (up);GeV;",200,0.,200.);
  histo1d_["2E_antiME_SSll_invM_CR_xFF_dn"] = fs->make<TH1D>("2E_antiME_SSll_invM_CR_xFF_dn","SS anti M(ll) x Fake factor (down);GeV;",200,0.,200.);
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

  histo1d_["2E_antiME_OSll_invM"] = fs->make<TH1D>("2E_antiME_OSll_invM","OS anti M(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_invM_CR"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR","OS anti M(ll);GeV;",200,0.,200.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF","OS anti M(ll) x Fake factor;GeV;",200,0.,200.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF_up"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF_up","OS anti M(ll) x Fake factor (up);GeV;",200,0.,200.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF_dn"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF_dn","OS anti M(ll) x Fake factor (down);GeV;",200,0.,200.);
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
  histo1d_["2E_Et_SSCR_EB_CRME"] = fs->make<TH1D>("2E_Et_SSCR_EB_CRME","SSCR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedME"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedME","SSCR (mixed) EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_EB_mixedAntiME"] = fs->make<TH1D>("2E_Et_SSCR_EB_mixedAntiME","SSCR (mixed) EB denominator E_{T};GeV;",200,0.,1000.);

  histo1d_["2E_Et_OSCR_EB_antiME"] = fs->make<TH1D>("2E_Et_OSCR_EB_antiME","OSCR EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_CRME"] = fs->make<TH1D>("2E_Et_OSCR_EB_CRME","OSCR EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedME"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedME","OSCR (mixed) EB numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_EB_mixedAntiME"] = fs->make<TH1D>("2E_Et_OSCR_EB_mixedAntiME","OSCR (mixed) EB denominator E_{T};GeV;",200,0.,1000.);

  histo1d_["2E_SFCR_OSll_invM_pass"] = fs->make<TH1D>("2E_SFCR_OSll_invM_pass","OS passing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_fail"] = fs->make<TH1D>("2E_SFCR_OSll_invM_fail","OS failing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_pass_Et-50to100"] = fs->make<TH1D>("2E_SFCR_OSll_invM_pass_Et-50to100","OS passing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_fail_Et-50to100"] = fs->make<TH1D>("2E_SFCR_OSll_invM_fail_Et-50to100","OS failing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_pass_Et-100to300"] = fs->make<TH1D>("2E_SFCR_OSll_invM_pass_Et-100to300","OS passing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_fail_Et-100to300"] = fs->make<TH1D>("2E_SFCR_OSll_invM_fail_Et-100to300","OS failing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_pass_Et-300to1000"] = fs->make<TH1D>("2E_SFCR_OSll_invM_pass_Et-300to1000","OS passing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_fail_Et-300to1000"] = fs->make<TH1D>("2E_SFCR_OSll_invM_fail_Et-300to1000","OS failing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_Et_numer"] = fs->make<TH1D>("2E_SFCR_OSll_Et_numer","OS numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_SFCR_OSll_Et_denom"] = fs->make<TH1D>("2E_SFCR_OSll_Et_denom","OS denominator E_{T};GeV;",200,0.,1000.);

  histo1d_["2E_SFCR_OSll_invM_pass_EB"] = fs->make<TH1D>("2E_SFCR_OSll_invM_pass_EB","OS passing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_invM_fail_EB"] = fs->make<TH1D>("2E_SFCR_OSll_invM_fail_EB","OS failing M(ll);GeV;",200,0.,200.);
  histo1d_["2E_SFCR_OSll_Et_numer_EB"] = fs->make<TH1D>("2E_SFCR_OSll_Et_numer_EB","OS numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_SFCR_OSll_Et_denom_EB"] = fs->make<TH1D>("2E_SFCR_OSll_Et_denom_EB","OS denominator E_{T};GeV;",200,0.,1000.);

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

  histo1d_["3E_CRME_lll_invM"] = fs->make<TH1D>("3E_CRME_lll_invM","CR M(lll);GeV;",40,0.,200.);
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

  histo1d_["3E_antiME_lll_invM"] = fs->make<TH1D>("3E_antiME_lll_invM","CR M(lll);GeV;",200,0.,1000.);
  histo1d_["3E_antiME_lll_invM_CR"] = fs->make<TH1D>("3E_antiME_lll_invM_CR","CR M(lll);GeV;",200,0,200.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF","CR M(lll) x SS Fake factor;GeV;",200,0,200.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF","CR M(lll) x OS Fake factor;GeV;",200,0,200.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_up"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_up","CR M(lll) x SS Fake factor (up);GeV;",200,0,200.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_up"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_up","CR M(lll) x OS Fake factor (up);GeV;",200,0,200.);
  histo1d_["3E_antiME_lll_invM_CR_xSSFF_dn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xSSFF_dn","CR M(lll) x SS Fake factor (down);GeV;",200,0,200.);
  histo1d_["3E_antiME_lll_invM_CR_xOSFF_dn"] = fs->make<TH1D>("3E_antiME_lll_invM_CR_xOSFF_dn","CR M(lll) x OS Fake factor (down);GeV;",200,0,200.);
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

  histo1d_["3E_Et_CR_EB_antiME"] = fs->make<TH1D>("3E_Et_CR_EB_antiME","3EVR EB denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_EB_CRME"] = fs->make<TH1D>("3E_Et_CR_EB_CRME","3EVR EB numerator E_{T};GeV;",200,0.,1000.);
}

void MergedEleCRanalyzer::endJob() {
  recoSFfile_->Close();
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
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  // prevent DY-ZG double-counting
  if (isDY_) {
    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(srcGenPtc_, genptcHandle);

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( ( std::abs(genptc->pdgId())==22 ) && genptc->fromHardProcessFinalState() )
        return;
    }
  }

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

  histo1d_["cutflow_4E"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_3E"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_2E"]->Fill( 0.5, aWeight );

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

  histo1d_["cutflow_4E"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_3E"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_2E"]->Fill( 1.5, aWeight );

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->ptrAt(0);
  else
    return;

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  // first two (modified)HEEP electrons should be Et > 35 GeV
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

    if ( aEle->electronID("modifiedHeepElectronID") ) {
      auto castEle = aEle.castTo<pat::ElectronRef>();
      acceptEles.push_back(castEle);
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
  std::vector<pat::ElectronRef> CRMEs;
  std::vector<pat::ElectronRef> antiMEs;

  for (unsigned int idx = 0; idx < acceptEles.size(); ++idx) {
    auto aEle = acceptEles.at(idx);

    double drGSF = std::sqrt(reco::deltaR2(aEle->gsfTrack()->eta(),
                                           aEle->gsfTrack()->phi(),
                                           (*addGsfTrkMap)[aEle]->eta(),
                                           (*addGsfTrkMap)[aEle]->phi()));

    std::string postfix = "";

    if ( aEle->et() < 50. )
      continue;

    if ( aEle->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
      // no 2nd GSF
      postfix = "_bkgEt2EB";

      if ( aEle->et() < ssBoundary_ )
        postfix = "_bkgEt1EB_extended";

      histo1d_["et_noGsf_HEEP"]->Fill( aEle->et(), aWeight );
      histo1d_[static_cast<std::string>("mva")+postfix]
        ->Fill( aEle->userFloat("mvaMergedElectronValues") , aWeight );

    } else {
      // has 2nd GSF
      histo2d_["et_dr_HEEP"]->Fill( aEle->et(), drGSF, aWeight );

      switch ( static_cast<MergedLeptonIDs::GSFtype>(aEle->userInt("mvaMergedElectronCategories")) ) { // better way to do this?
        case MergedLeptonIDs::GSFtype::DR1Et2EB:
          postfix = "_DR1Et2EB";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et1EB:
          postfix = "_DR2Et1EB";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et2EB:
          postfix = "_DR2Et2EB";
          break;
        case MergedLeptonIDs::GSFtype::extendedCR:
          postfix = "_DR1Et1EB_extended"; // TODO empty histo for the moment
          break;
        default:
          throw cms::Exception("LogicError") << "opening angle between the original and additional GSF track does not fit into any of MergedLeptonIDs::openingAngle categories" << std::endl;
          break;
      }

      histo1d_[static_cast<std::string>("mva")+postfix]
        ->Fill( aEle->userFloat("mvaMergedElectronValues") , aWeight );
    }

    if ( aEle->electronID("mvaMergedElectron") && aEle->passConversionVeto() ) {
      if ( aEle->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::DR1Et2EB) ||
           aEle->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
        if ( aEle->et() > ssBoundary_ )
          mergedEls.push_back(aEle);
      } else {
        mergedEls.push_back(aEle);
      }
    }

    if ( aEle->electronID("mvaMergedElectron") && aEle->passConversionVeto() && aEle->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
      histo2d_["et_dr_final"]->Fill( aEle->et(), drGSF, aWeight );
      const auto lvecE = math::PtEtaPhiMLorentzVector(aEle->gsfTrack()->pt(),aEle->gsfTrack()->eta(),aEle->gsfTrack()->phi(),0.);
      const auto lvecGSF = math::PtEtaPhiMLorentzVector((*addGsfTrkMap)[aEle]->pt(),(*addGsfTrkMap)[aEle]->eta(),(*addGsfTrkMap)[aEle]->phi(),0.);
      const auto lvecMEll = lvecE + lvecGSF;

      if ( std::abs(aEle->superCluster()->eta()) < 1.5 ) {
        histo1d_["invM_ll_CRME_EB"]->Fill( lvecMEll.M() , aWeight );

        if (aEle->et() > 200.)
          histo1d_["invM_ll_CRME_Et200_EB"]->Fill( lvecMEll.M() , aWeight );
      }
    }

    if ( aEle->electronID("mvaMergedElectron") && aEle->passConversionVeto() && aEle->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) )
      histo1d_["et_noGsf_final"]->Fill( aEle->et(), aWeight );

    if ( aEle->electronID("mvaMergedElectron") && aEle->passConversionVeto() )
      CRMEs.push_back(aEle);

    // electrons with mva score = -1 are out of our interest
    if ( !( static_cast<bool>(aEle->electronID("mvaMergedElectron")) && aEle->passConversionVeto() ) && aEle->userFloat("mvaMergedElectronValues")!=-1. )
      antiMEs.push_back(aEle);
  } // electron loop

  std::sort(mergedEls.begin(),mergedEls.end(),sortByEt);
  std::sort(CRMEs.begin(),CRMEs.end(),sortByEt);
  std::sort(antiMEs.begin(),antiMEs.end(),sortByEt);

  // recoSF
  if (isMC_) {
    for (const auto& aEle : acceptEles) {
      TH2D* aSF = aEle->pt() >= 20. ? recoSF_ : recoSFlowPt_;
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      aWeight *= aSF->GetBinContent( aSF->FindBin(aEle->superCluster()->eta(), apt ) );
    }
  }

  switch (acceptEles.size()) {
    case 4:
      histo1d_["cutflow_4E"]->Fill( 2.5, aWeight );
      if ( mergedEls.size()==0 )
        histo1d_["cutflow_4E"]->Fill( 3.5, aWeight );
      break;
    case 3:
      histo1d_["cutflow_3E"]->Fill( 2.5, aWeight );

      if ( mergedEls.size()==1 )
        histo1d_["cutflow_3E"]->Fill( 3.5, aWeight );

      if ( CRMEs.size()==1 ) {
        std::vector<pat::ElectronRef> nonMEs;

        for (const auto& aEle : acceptEles) {
          if ( aEle!=CRMEs.front() )
            nonMEs.push_back(aEle);
        }

        const auto lvecCRME = CRMEs.front()->polarP4()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy();
        const auto lvecNME1 = nonMEs.front()->polarP4()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy();
        const auto lvecNME2 = nonMEs.at(1)->polarP4()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy();
        const auto lvecCRllME = lvecCRME + lvecNME1 + lvecNME2;
        const double dr2l1ME = reco::deltaR2(CRMEs.front()->eta(),CRMEs.front()->phi(),nonMEs.front()->eta(),nonMEs.front()->phi());
        const double dr2l2ME = reco::deltaR2(CRMEs.front()->eta(),CRMEs.front()->phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double dr2l1l2 = reco::deltaR2(nonMEs.front()->eta(),nonMEs.front()->phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double mlll = lvecCRllME.M();

        if ( mlll > 50. && mlll < 200. ) {
          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["3E_CRME_Et"]->Fill(CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(),aWeight);
            histo1d_["3E_CRME_eta"]->Fill(CRMEs.front()->eta(),aWeight);
            histo1d_["3E_CRME_phi"]->Fill(CRMEs.front()->phi(),aWeight);
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["3E_CRME_noGsf_Et"]->Fill(CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(),aWeight);
            histo1d_["3E_CRME_noGsf_eta"]->Fill(CRMEs.front()->eta(),aWeight);
            histo1d_["3E_CRME_noGsf_phi"]->Fill(CRMEs.front()->phi(),aWeight);
          } else {}

          histo1d_["3E_CRME_all_Et"]->Fill(CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(),aWeight);
          histo1d_["3E_CRME_all_eta"]->Fill(CRMEs.front()->eta(),aWeight);
          histo1d_["3E_CRME_all_phi"]->Fill(CRMEs.front()->phi(),aWeight);

          histo1d_["3E_CR_NME1_Et"]->Fill(nonMEs.front()->et()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy(),aWeight);
          histo1d_["3E_CR_NME1_eta"]->Fill(nonMEs.front()->eta(),aWeight);
          histo1d_["3E_CR_NME1_phi"]->Fill(nonMEs.front()->phi(),aWeight);

          histo1d_["3E_CR_NME2_Et"]->Fill(nonMEs.at(1)->et()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy(),aWeight);
          histo1d_["3E_CR_NME2_eta"]->Fill(nonMEs.at(1)->eta(),aWeight);
          histo1d_["3E_CR_NME2_phi"]->Fill(nonMEs.at(1)->phi(),aWeight);

          histo1d_["3E_CRME_lll_invM"]->Fill(mlll,aWeight);
          histo1d_["3E_CRME_lll_pt"]->Fill(lvecCRllME.pt(),aWeight);
          histo1d_["3E_CRME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
          histo1d_["3E_CRME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
          histo1d_["3E_CRME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

          if ( CRMEs.front()->isEB() )
            histo1d_["3E_Et_CR_EB_CRME"]->Fill( CRMEs.front()->et(), aWeight );
        }
      } // CRMEs.size()==1

      if ( CRMEs.size()==0 && antiMEs.size()>0 ) {
        std::vector<pat::ElectronRef> nonMEs;

        for (const auto& aEle : acceptEles) {
          if ( aEle!=antiMEs.front() )
            nonMEs.push_back(aEle);
        }

        const auto lvecAME = antiMEs.front()->polarP4()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy();
        const auto lvecNME1 = nonMEs.front()->polarP4()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy();
        const auto lvecNME2 = nonMEs.at(1)->polarP4()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy();
        const auto lvecAllME = lvecAME + lvecNME1 + lvecNME2;
        const double dr2l1ME = reco::deltaR2(antiMEs.front()->eta(),antiMEs.front()->phi(),nonMEs.front()->eta(),nonMEs.front()->phi());
        const double dr2l2ME = reco::deltaR2(antiMEs.front()->eta(),antiMEs.front()->phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double dr2l1l2 = reco::deltaR2(nonMEs.front()->eta(),nonMEs.front()->phi(),nonMEs.at(1)->eta(),nonMEs.at(1)->phi());
        const double mlll = lvecAllME.M();

        if (mlll > 50.) {
          if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["3E_antiME_Et"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight);
            histo1d_["3E_antiME_eta"]->Fill(antiMEs.front()->eta(),aWeight);
            histo1d_["3E_antiME_phi"]->Fill(antiMEs.front()->phi(),aWeight);
          } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["3E_antiME_noGsf_Et"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight);
            histo1d_["3E_antiME_noGsf_eta"]->Fill(antiMEs.front()->eta(),aWeight);
            histo1d_["3E_antiME_noGsf_phi"]->Fill(antiMEs.front()->phi(),aWeight);
          } else {}

          histo1d_["3E_antiME_NME1_Et"]->Fill(nonMEs.front()->et()*nonMEs.front()->userFloat("ecalTrkEnergyPostCorr")/nonMEs.front()->energy(),aWeight);
          histo1d_["3E_antiME_NME1_eta"]->Fill(nonMEs.front()->eta(),aWeight);
          histo1d_["3E_antiME_NME1_phi"]->Fill(nonMEs.front()->phi(),aWeight);

          histo1d_["3E_antiME_NME2_Et"]->Fill(nonMEs.at(1)->et()*nonMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/nonMEs.at(1)->energy(),aWeight);
          histo1d_["3E_antiME_NME2_eta"]->Fill(nonMEs.at(1)->eta(),aWeight);
          histo1d_["3E_antiME_NME2_phi"]->Fill(nonMEs.at(1)->phi(),aWeight);

          histo1d_["3E_antiME_lll_invM"]->Fill(mlll,aWeight);
          histo1d_["3E_antiME_lll_pt"]->Fill(lvecAllME.pt(),aWeight);
          histo1d_["3E_antiME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
          histo1d_["3E_antiME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
          histo1d_["3E_antiME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

          if ( mlll < 200. ) {
            const double aEt = antiMEs.front()->et();
            const TF1* osfunc = osboth_;
            const TF1* ssfunc = ssboth_;
            const double ffOS = osfunc->Eval(aEt);
            const double ffSS = ssfunc->Eval(aEt);
            const double xval[1] = {aEt};
            double ciOS[1], ciSS[1];
            osFit_->GetConfidenceIntervals(1,1,0,xval,ciOS,0.6827,false);
            ssFit_->GetConfidenceIntervals(1,1,0,xval,ciSS,0.6827,false);
            histo1d_["3E_antiME_lll_invM_CR"]->Fill(mlll,aWeight);
            histo1d_["3E_antiME_lll_invM_CR_xOSFF"]->Fill(mlll,aWeight*ffOS);
            histo1d_["3E_antiME_lll_invM_CR_xSSFF"]->Fill(mlll,aWeight*ffSS);
            histo1d_["3E_antiME_lll_invM_CR_xOSFF_up"]->Fill(mlll,aWeight*(ffOS+ciOS[0]));
            histo1d_["3E_antiME_lll_invM_CR_xSSFF_up"]->Fill(mlll,aWeight*(ffSS+ciSS[0]));
            histo1d_["3E_antiME_lll_invM_CR_xOSFF_dn"]->Fill(mlll,aWeight*(ffOS-ciOS[0]));
            histo1d_["3E_antiME_lll_invM_CR_xSSFF_dn"]->Fill(mlll,aWeight*(ffSS-ciSS[0]));

            histo1d_["3E_antiME_Et_xSSFF"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight*ffSS);
            histo1d_["3E_antiME_eta_xSSFF"]->Fill(antiMEs.front()->eta(),aWeight*ffSS);
            histo1d_["3E_antiME_Et_xOSFF"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight*ffOS);
            histo1d_["3E_antiME_eta_xOSFF"]->Fill(antiMEs.front()->eta(),aWeight*ffOS);
            histo1d_["3E_antiME_Et_xSSFF_up"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight*(ffSS+ciSS[0]));
            histo1d_["3E_antiME_eta_xSSFF_up"]->Fill(antiMEs.front()->eta(),aWeight*(ffSS+ciSS[0]));
            histo1d_["3E_antiME_Et_xOSFF_up"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight*(ffOS+ciOS[0]));
            histo1d_["3E_antiME_eta_xOSFF_up"]->Fill(antiMEs.front()->eta(),aWeight*(ffOS+ciOS[0]));
            histo1d_["3E_antiME_Et_xSSFF_dn"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight*(ffSS-ciSS[0]));
            histo1d_["3E_antiME_eta_xSSFF_dn"]->Fill(antiMEs.front()->eta(),aWeight*(ffSS-ciSS[0]));
            histo1d_["3E_antiME_Et_xOSFF_dn"]->Fill(antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(),aWeight*(ffOS-ciOS[0]));
            histo1d_["3E_antiME_eta_xOSFF_dn"]->Fill(antiMEs.front()->eta(),aWeight*(ffOS-ciOS[0]));

            if ( antiMEs.front()->isEB() )
              histo1d_["3E_Et_CR_EB_antiME"]->Fill( antiMEs.front()->et(), aWeight );
          } // CR
        } // mlll > 50.
      } // antiMEs.size()==1

      break;
    case 2:
      histo1d_["cutflow_2E"]->Fill( 2.5, aWeight );
      if ( mergedEls.size() > 0 )
        histo1d_["cutflow_2E"]->Fill( 3.5, aWeight );
      if ( mergedEls.size()==2 )
        histo1d_["cutflow_2E"]->Fill( 4.5, aWeight );

      const auto lvecE1 = acceptEles.front()->polarP4()*acceptEles.front()->userFloat("ecalTrkEnergyPostCorr")/acceptEles.front()->energy();
      const auto lvecE2 = acceptEles.at(1)->polarP4()*acceptEles.at(1)->userFloat("ecalTrkEnergyPostCorr")/acceptEles.at(1)->energy();
      const auto lvecll = lvecE1 + lvecE2;
      const double mll_SF = lvecll.M();

      // special CR for pseudo-SF estimation
      if ( mll_SF > 50. && mll_SF < 200. && acceptEles.front()->charge()*acceptEles.at(1)->charge() < 0. ) {
        for (const auto& aEle : acceptEles) {
          if ( std::find(CRMEs.begin(),CRMEs.end(),aEle)!=CRMEs.end() ) {
            histo1d_["2E_SFCR_OSll_invM_pass"]->Fill( mll_SF, aWeight );

            if ( std::abs(aEle->superCluster()->eta()) < 1.5 )
              histo1d_["2E_SFCR_OSll_invM_pass_EB"]->Fill( mll_SF, aWeight );

            if ( aEle->et() > 50. && aEle->et() < 100. )
              histo1d_["2E_SFCR_OSll_invM_pass_Et-50to100"]->Fill( mll_SF, aWeight );
            else if ( aEle->et() > 100. && aEle->et() < 300. )
              histo1d_["2E_SFCR_OSll_invM_pass_Et-100to300"]->Fill( mll_SF, aWeight );
            else if ( aEle->et() > 300. && aEle->et() < 1000. )
              histo1d_["2E_SFCR_OSll_invM_pass_Et-300to1000"]->Fill( mll_SF, aWeight );
            else {}

          } else {
            histo1d_["2E_SFCR_OSll_invM_fail"]->Fill( mll_SF, aWeight );

            if ( std::abs(aEle->superCluster()->eta()) < 1.5 )
              histo1d_["2E_SFCR_OSll_invM_fail_EB"]->Fill( mll_SF, aWeight );

            if ( aEle->et() > 50. && aEle->et() < 100. )
              histo1d_["2E_SFCR_OSll_invM_fail_Et-50to100"]->Fill( mll_SF, aWeight );
            else if ( aEle->et() > 100. && aEle->et() < 300. )
              histo1d_["2E_SFCR_OSll_invM_fail_Et-100to300"]->Fill( mll_SF, aWeight );
            else if ( aEle->et() > 300. && aEle->et() < 1000. )
              histo1d_["2E_SFCR_OSll_invM_fail_Et-300to1000"]->Fill( mll_SF, aWeight );
            else {}

          }
        }

        if ( mll_SF > 80. && mll_SF < 100. ) {
          for (const auto& aEle : acceptEles) {
            histo1d_["2E_SFCR_OSll_Et_denom"]->Fill( aEle->et(), aWeight );

            if ( std::abs(aEle->superCluster()->eta()) < 1.5 )
              histo1d_["2E_SFCR_OSll_Et_denom_EB"]->Fill( aEle->et(), aWeight );

            if ( std::find(CRMEs.begin(),CRMEs.end(),aEle)!=CRMEs.end() ) {
              histo1d_["2E_SFCR_OSll_Et_numer"]->Fill( aEle->et(), aWeight );

              if ( std::abs(aEle->superCluster()->eta()) < 1.5 )
                histo1d_["2E_SFCR_OSll_Et_numer_EB"]->Fill( aEle->et(), aWeight );
            }
          }
        }
      } // SFCR

      if (CRMEs.size()==2) {
        const auto lvecCRME1 = CRMEs.front()->polarP4()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy();
        const auto lvecCRME2 = CRMEs.at(1)->polarP4()*CRMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/CRMEs.at(1)->energy();
        const auto lvecCRll = lvecCRME1 + lvecCRME2;
        const double dr2CRll = reco::deltaR2(lvecCRME1.eta(),lvecCRME1.phi(),lvecCRME2.eta(),lvecCRME2.phi());
        const double mll = lvecCRll.M();

        if ( mll > 50. && mll < 200. ) {
          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_CRME1_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_CRME1_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_CRME1_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_CRME1_noGsf_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_CRME1_noGsf_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_CRME1_noGsf_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else {}

          if ( CRMEs.at(1)->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_CRME2_Et"]->Fill( CRMEs.at(1)->et()*CRMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/CRMEs.at(1)->energy(), aWeight );
            histo1d_["2E_CRME2_eta"]->Fill( CRMEs.at(1)->eta(), aWeight );
            histo1d_["2E_CRME2_phi"]->Fill( CRMEs.at(1)->phi(), aWeight );
          } else if ( CRMEs.at(1)->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_CRME2_noGsf_Et"]->Fill( CRMEs.at(1)->et()*CRMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/CRMEs.at(1)->energy(), aWeight );
            histo1d_["2E_CRME2_noGsf_eta"]->Fill( CRMEs.at(1)->eta(), aWeight );
            histo1d_["2E_CRME2_noGsf_phi"]->Fill( CRMEs.at(1)->phi(), aWeight );
          } else {}

          histo1d_["2E_CRME_ll_invM"]->Fill( mll, aWeight );
          histo1d_["2E_CRME_ll_pt"]->Fill( lvecCRll.pt(), aWeight );
          histo1d_["2E_CRME_ll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

          if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() < 0. ) { // OS
            histo1d_["2E_CRME_OSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_CRME_OSll_pt"]->Fill( lvecCRll.pt(), aWeight );
            histo1d_["2E_CRME_OSll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

            for (const auto& aME : CRMEs) {
              histo1d_["2E_CRME_OSll_Et"]->Fill( aME->et()*aME->userFloat("ecalTrkEnergyPostCorr")/aME->energy(), aWeight );
              histo1d_["2E_CRME_OSll_eta"]->Fill( aME->eta(), aWeight );

              if ( aME->isEB() )
                histo1d_["2E_Et_OSCR_EB_CRME"]->Fill( aME->et(), aWeight );
            }
          } else if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_CRME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_CRME_SSll_pt"]->Fill( lvecCRll.pt(), aWeight );
            histo1d_["2E_CRME_SSll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

            for (const auto& aME : CRMEs) {
              histo1d_["2E_CRME_SSll_Et"]->Fill( aME->et()*aME->userFloat("ecalTrkEnergyPostCorr")/aME->energy(), aWeight );
              histo1d_["2E_CRME_SSll_eta"]->Fill( aME->eta(), aWeight );

              if ( aME->isEB() )
                histo1d_["2E_Et_SSCR_EB_CRME"]->Fill( aME->et(), aWeight );
            }
          }
        } // CR
      } // CRMEs.size()==2

      if (CRMEs.size()==1 && antiMEs.size()==1) {
        const auto lvecCRME = CRMEs.front()->polarP4()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy();
        const auto lvecAnME = antiMEs.front()->polarP4()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy();
        const auto lvecMixedll = lvecCRME + lvecAnME;
        const double dr2mixedll = reco::deltaR2(lvecCRME.eta(),lvecCRME.phi(),lvecAnME.eta(),lvecAnME.phi());
        const double mll = lvecMixedll.M();

        if ( mll > 50. && mll < 200. ) {
          if ( CRMEs.front()->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_mixed_ME_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_ME_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_ME_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else if ( CRMEs.front()->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_mixed_ME_noGsf_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_ME_noGsf_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_ME_noGsf_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else {}

          if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_mixed_antiME_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_antiME_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_antiME_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_mixed_antiME_noGsf_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_antiME_noGsf_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_antiME_noGsf_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else {}

          histo1d_["2E_mixedME_ll_invM"]->Fill( mll, aWeight );
          histo1d_["2E_mixedME_ll_pt"]->Fill( lvecMixedll.pt(), aWeight );
          histo1d_["2E_mixedME_ll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

          const double aEt = antiMEs.front()->et();
          const TF1* osfunc = osboth_;
          const TF1* ssfunc = ssboth_;
          const double ffOS = osfunc->Eval(aEt);
          const double ffSS = ssfunc->Eval(aEt);
          const double xval[1] = {aEt};
          double ciOS[1], ciSS[1];
          osFit_->GetConfidenceIntervals(1,1,0,xval,ciOS,0.6827,false);
          ssFit_->GetConfidenceIntervals(1,1,0,xval,ciSS,0.6827,false);

          if ( CRMEs.front()->charge()*antiMEs.front()->charge() < 0. ) { // OS
            histo1d_["2E_mixedME_OSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_mixedME_OSll_invM_xFF"]->Fill( mll, 0.5*aWeight*ffOS ); // double-counting
            histo1d_["2E_mixedME_OSll_invM_xFF_up"]->Fill( mll, 0.5*aWeight*(ffOS+ciOS[0]) ); // double-counting
            histo1d_["2E_mixedME_OSll_invM_xFF_dn"]->Fill( mll, 0.5*aWeight*(ffOS-ciOS[0]) ); // double-counting
            histo1d_["2E_mixedME_OSll_pt"]->Fill( lvecMixedll.pt(), aWeight );
            histo1d_["2E_mixedME_OSll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

            histo1d_["2E_mixedME_OSll_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_mixedME_OSll_Et_noCorr"]->Fill( CRMEs.front()->et(), aWeight );
            histo1d_["2E_mixedME_OSll_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_mixedAntiME_OSll_Et_xFF"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*ffOS );
            histo1d_["2E_mixedAntiME_OSll_eta_xFF"]->Fill( antiMEs.front()->eta(), aWeight*ffOS );
            histo1d_["2E_mixedAntiME_OSll_Et_xFF_up"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffOS+ciOS[0]) );
            histo1d_["2E_mixedAntiME_OSll_eta_xFF_up"]->Fill( antiMEs.front()->eta(), aWeight*(ffOS+ciOS[0]) );
            histo1d_["2E_mixedAntiME_OSll_Et_xFF_dn"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffOS-ciOS[0]) );
            histo1d_["2E_mixedAntiME_OSll_eta_xFF_dn"]->Fill( antiMEs.front()->eta(), aWeight*(ffOS-ciOS[0]) );

            if ( CRMEs.front()->isEB() )
              histo1d_["2E_Et_OSCR_EB_mixedME"]->Fill( CRMEs.front()->et(), aWeight );

            if ( antiMEs.front()->isEB() )
              histo1d_["2E_Et_OSCR_EB_mixedAntiME"]->Fill( antiMEs.front()->et(), aWeight );
          } else if ( CRMEs.front()->charge()*antiMEs.front()->charge() > 0. ) { // SS
            histo1d_["2E_mixedME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_mixedME_SSll_invM_xFF"]->Fill( mll, 0.5*aWeight*ffSS ); // double-counting
            histo1d_["2E_mixedME_SSll_invM_xFF_up"]->Fill( mll, 0.5*aWeight*(ffSS+ciSS[0]) ); // double-counting
            histo1d_["2E_mixedME_SSll_invM_xFF_dn"]->Fill( mll, 0.5*aWeight*(ffSS-ciSS[0]) ); // double-counting
            histo1d_["2E_mixedME_SSll_pt"]->Fill( lvecMixedll.pt(), aWeight );
            histo1d_["2E_mixedME_SSll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

            histo1d_["2E_mixedME_SSll_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_mixedME_SSll_Et_noCorr"]->Fill( CRMEs.front()->et(), aWeight );
            histo1d_["2E_mixedME_SSll_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_mixedAntiME_SSll_Et_xFF"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*ffSS );
            histo1d_["2E_mixedAntiME_SSll_eta_xFF"]->Fill( antiMEs.front()->eta(), aWeight*ffSS );
            histo1d_["2E_mixedAntiME_SSll_Et_xFF_up"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffSS+ciSS[0]) );
            histo1d_["2E_mixedAntiME_SSll_eta_xFF_up"]->Fill( antiMEs.front()->eta(), aWeight*(ffSS+ciSS[0]) );
            histo1d_["2E_mixedAntiME_SSll_Et_xFF_dn"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffSS-ciSS[0]) );
            histo1d_["2E_mixedAntiME_SSll_eta_xFF_dn"]->Fill( antiMEs.front()->eta(), aWeight*(ffSS-ciSS[0]) );

            if ( CRMEs.front()->isEB() )
              histo1d_["2E_Et_SSCR_EB_mixedME"]->Fill( CRMEs.front()->et(), aWeight );

            if ( antiMEs.front()->isEB() )
              histo1d_["2E_Et_SSCR_EB_mixedAntiME"]->Fill( antiMEs.front()->et(), aWeight );
          }
        } // CR
      } // CRMEs.size()==1 && antiMEs.size()==1

      if (antiMEs.size()==2) {
        const auto lvecAME1 = antiMEs.front()->polarP4()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy();
        const auto lvecAME2 = antiMEs.at(1)->polarP4()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy();
        const auto lvecAll = lvecAME1 + lvecAME2;
        const double dr2All = reco::deltaR2(lvecAME1.eta(),lvecAME1.phi(),lvecAME2.eta(),lvecAME2.phi());
        const double mll = lvecAll.M();

        if ( mll > 50. ) {
          if ( antiMEs.front()->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_antiME1_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_antiME1_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_antiME1_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else if ( antiMEs.front()->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_antiME1_noGsf_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_antiME1_noGsf_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_antiME1_noGsf_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else {}

          if ( antiMEs.at(1)->userInt("mvaMergedElectronCategories")!=static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_antiME2_Et"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight );
            histo1d_["2E_antiME2_eta"]->Fill( antiMEs.at(1)->eta(), aWeight );
            histo1d_["2E_antiME2_phi"]->Fill( antiMEs.at(1)->phi(), aWeight );
          } else if ( antiMEs.at(1)->userInt("mvaMergedElectronCategories")==static_cast<int>(MergedLeptonIDs::GSFtype::bkgEt2EB) ) {
            histo1d_["2E_antiME2_noGsf_Et"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight );
            histo1d_["2E_antiME2_noGsf_eta"]->Fill( antiMEs.at(1)->eta(), aWeight );
            histo1d_["2E_antiME2_noGsf_phi"]->Fill( antiMEs.at(1)->phi(), aWeight );
          } else {}

          histo1d_["2E_antiME_ll_invM"]->Fill( mll, aWeight );
          histo1d_["2E_antiME_ll_pt"]->Fill( lvecAll.pt(), aWeight );
          histo1d_["2E_antiME_ll_dr"]->Fill( std::sqrt(dr2All), aWeight );

          const double aEt1 = antiMEs.front()->et();
          const TF1* osfunc1 = osboth_;
          const TF1* ssfunc1 = ssboth_;
          const double ffOS1 = osfunc1->Eval(aEt1);
          const double ffSS1 = ssfunc1->Eval(aEt1);
          const double aEt2 = antiMEs.at(1)->et();
          const TF1* osfunc2 = osboth_;
          const TF1* ssfunc2 = ssboth_;
          const double ffOS2 = osfunc2->Eval(aEt2);
          const double ffSS2 = ssfunc2->Eval(aEt2);
          const double xval[2] = {aEt1,aEt2};
          double ciOS[2], ciSS[2];
          osFit_->GetConfidenceIntervals(2,1,0,xval,ciOS,0.6827,false);
          ssFit_->GetConfidenceIntervals(2,1,0,xval,ciSS,0.6827,false);

          if (mll < 200.)
            histo1d_["2E_antiME_ll_invM_CR"]->Fill( lvecAll.M(), aWeight );

          if ( antiMEs.front()->charge()*antiMEs.at(1)->charge() < 0. ) { // OS
            histo1d_["2E_antiME_OSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_antiME_OSll_pt"]->Fill( lvecAll.pt(), aWeight );
            histo1d_["2E_antiME_OSll_dr"]->Fill( std::sqrt(dr2All), aWeight );

            if (mll < 200.) {
              histo1d_["2E_antiME_OSll_invM_CR"]->Fill( mll, aWeight );
              histo1d_["2E_antiME_OSll_invM_CR_xFF"]->Fill( mll, aWeight*ffOS1 );
              histo1d_["2E_antiME_OSll_invM_CR_xFF"]->Fill( mll, aWeight*ffOS2 );
              histo1d_["2E_antiME_OSll_invM_CR_xFF_up"]->Fill( mll, aWeight*(ffOS1+ciOS[0]) );
              histo1d_["2E_antiME_OSll_invM_CR_xFF_up"]->Fill( mll, aWeight*(ffOS2+ciOS[1]) );
              histo1d_["2E_antiME_OSll_invM_CR_xFF_dn"]->Fill( mll, aWeight*(ffOS1-ciOS[0]) );
              histo1d_["2E_antiME_OSll_invM_CR_xFF_dn"]->Fill( mll, aWeight*(ffOS2-ciOS[1]) );

              histo1d_["2E_antiME_OSll_Et_xFF"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*ffOS1 );
              histo1d_["2E_antiME_OSll_eta_xFF"]->Fill( antiMEs.front()->eta(), aWeight*ffOS1 );
              histo1d_["2E_antiME_OSll_Et_xFF"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight*ffOS2 );
              histo1d_["2E_antiME_OSll_eta_xFF"]->Fill( antiMEs.at(1)->eta(), aWeight*ffOS2 );
              histo1d_["2E_antiME_OSll_Et_xFF_up"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffOS1+ciOS[0]) );
              histo1d_["2E_antiME_OSll_eta_xFF_up"]->Fill( antiMEs.front()->eta(), aWeight*(ffOS1+ciOS[0]) );
              histo1d_["2E_antiME_OSll_Et_xFF_up"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight*(ffOS2+ciOS[1]) );
              histo1d_["2E_antiME_OSll_eta_xFF_up"]->Fill( antiMEs.at(1)->eta(), aWeight*(ffOS2+ciOS[1]) );
              histo1d_["2E_antiME_OSll_Et_xFF_dn"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffOS1-ciOS[0]) );
              histo1d_["2E_antiME_OSll_eta_xFF_dn"]->Fill( antiMEs.front()->eta(), aWeight*(ffOS1-ciOS[0]) );
              histo1d_["2E_antiME_OSll_Et_xFF_dn"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight*(ffOS2-ciOS[1]) );
              histo1d_["2E_antiME_OSll_eta_xFF_dn"]->Fill( antiMEs.at(1)->eta(), aWeight*(ffOS2-ciOS[1]) );

              for (const auto& aME : antiMEs) {
                histo1d_["2E_antiME_OSll_Et_noCorr"]->Fill( aME->et(), aWeight );
                histo1d_["2E_antiME_OSll_eta"]->Fill( aME->eta(), aWeight );

                if ( aME->isEB() )
                  histo1d_["2E_Et_OSCR_EB_antiME"]->Fill( aME->et(), aWeight );
              }
            }
          } else if ( antiMEs.front()->charge()*antiMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_antiME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_antiME_SSll_pt"]->Fill( lvecAll.pt(), aWeight );
            histo1d_["2E_antiME_SSll_dr"]->Fill( std::sqrt(dr2All), aWeight );

            if (mll < 200.) {
              histo1d_["2E_antiME_SSll_invM_CR"]->Fill( mll, aWeight );
              histo1d_["2E_antiME_SSll_invM_CR_xFF"]->Fill( mll, aWeight*ffSS1 );
              histo1d_["2E_antiME_SSll_invM_CR_xFF"]->Fill( mll, aWeight*ffSS2 );
              histo1d_["2E_antiME_SSll_invM_CR_xFF_up"]->Fill( mll, aWeight*(ffSS1+ciSS[0]) );
              histo1d_["2E_antiME_SSll_invM_CR_xFF_up"]->Fill( mll, aWeight*(ffSS2+ciSS[1]) );
              histo1d_["2E_antiME_SSll_invM_CR_xFF_dn"]->Fill( mll, aWeight*(ffSS1-ciSS[0]) );
              histo1d_["2E_antiME_SSll_invM_CR_xFF_dn"]->Fill( mll, aWeight*(ffSS2-ciSS[1]) );

              histo1d_["2E_antiME_SSll_Et_xFF"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*ffSS1 );
              histo1d_["2E_antiME_SSll_eta_xFF"]->Fill( antiMEs.front()->eta(), aWeight*ffSS1 );
              histo1d_["2E_antiME_SSll_Et_xFF"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight*ffSS2 );
              histo1d_["2E_antiME_SSll_eta_xFF"]->Fill( antiMEs.at(1)->eta(), aWeight*ffSS2 );
              histo1d_["2E_antiME_SSll_Et_xFF_up"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffSS1+ciSS[0]) );
              histo1d_["2E_antiME_SSll_eta_xFF_up"]->Fill( antiMEs.front()->eta(), aWeight*(ffSS1+ciSS[0]) );
              histo1d_["2E_antiME_SSll_Et_xFF_up"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight*(ffSS2+ciSS[1]) );
              histo1d_["2E_antiME_SSll_eta_xFF_up"]->Fill( antiMEs.at(1)->eta(), aWeight*(ffSS2+ciSS[1]) );
              histo1d_["2E_antiME_SSll_Et_xFF_dn"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight*(ffSS1-ciSS[0]) );
              histo1d_["2E_antiME_SSll_eta_xFF_dn"]->Fill( antiMEs.front()->eta(), aWeight*(ffSS1-ciSS[0]) );
              histo1d_["2E_antiME_SSll_Et_xFF_dn"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight*(ffSS2-ciSS[1]) );
              histo1d_["2E_antiME_SSll_eta_xFF_dn"]->Fill( antiMEs.at(1)->eta(), aWeight*(ffSS2-ciSS[1]) );

              for (const auto& aME : antiMEs) {
                histo1d_["2E_antiME_SSll_Et_noCorr"]->Fill( aME->et(), aWeight );
                histo1d_["2E_antiME_SSll_eta"]->Fill( aME->eta(), aWeight );

                if ( aME->isEB() )
                  histo1d_["2E_Et_SSCR_EB_antiME"]->Fill( aME->et(), aWeight );
              }
            } // CR
          } // OSSS
        } // mll > 50.
      } // antiME.size()==2

      break;
  } // switch acceptEles.size()
}

DEFINE_FWK_MODULE(MergedEleCRanalyzer);
