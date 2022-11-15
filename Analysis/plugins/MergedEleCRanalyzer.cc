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

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "TH2D.h"
#include "TF1.h"
#include "TFile.h"

class MergedEleCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleCRanalyzer(const edm::ParameterSet&);
  virtual ~MergedEleCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  bool isMC_;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> cutflow_modifiedHEEPToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> cutflow_HEEPToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> status_mergedElectronToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> mva_mergedElectronToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> GSFtype_mergedElectronToken_;

  const edm::FileInPath recoSFpath_;
  const edm::FileInPath FFpath_;
  const std::vector<std::string> trigList_;
  const double etThres1_;
  const double ssBoundary_;
  const double osBoundary_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  std::unique_ptr<TFile> recoSFfile_;
  TH2D* recoSF_;

  std::unique_ptr<TFile> FFfile_;
  TF1* sslow_;
  TF1* sshigh_;
  TF1* oslow_;
  TF1* oshigh_;
};

MergedEleCRanalyzer::MergedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
cutflow_modifiedHEEPToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("cutflow_modifiedHEEP"))),
cutflow_HEEPToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("cutflow_HEEP"))),
status_mergedElectronToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("status_mergedElectron"))),
mva_mergedElectronToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mva_mergedElectron"))),
GSFtype_mergedElectronToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("GSFtype_mergedElectron"))),
recoSFpath_(iConfig.getParameter<edm::FileInPath>("recoSFpath")),
FFpath_(iConfig.getParameter<edm::FileInPath>("FFpath")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
etThres1_(iConfig.getParameter<double>("etThres1")),
ssBoundary_(iConfig.getParameter<double>("ssBoundary")),
osBoundary_(iConfig.getParameter<double>("osBoundary")) {
  usesResource("TFileService");
}

void MergedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  recoSFfile_ = std::make_unique<TFile>(recoSFpath_.fullPath().c_str(),"READ");
  recoSF_ = static_cast<TH2D*>(recoSFfile_->Get("EGamma_SF2D"));

  FFfile_ = std::make_unique<TFile>(FFpath_.fullPath().c_str(),"READ");
  sslow_ = static_cast<TF1*>(FFfile_->FindObjectAny("sslow"));
  sshigh_ = static_cast<TF1*>(FFfile_->FindObjectAny("sshigh"));
  oslow_ = static_cast<TF1*>(FFfile_->FindObjectAny("oslow"));
  oshigh_ = static_cast<TF1*>(FFfile_->FindObjectAny("oshigh"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow_modifiedHEEP_EB"] = fs->make<TH1D>("cutflow_modifiedHEEP_EB","cutflow (modified HEEP EB)",11,-1.,10.);
  histo1d_["cutflow_modifiedHEEP_EE"] = fs->make<TH1D>("cutflow_modifiedHEEP_EE","cutflow (modified HEEP EE)",11,-1.,10.);
  histo1d_["cutflow_HEEP_EB"] = fs->make<TH1D>("cutflow_HEEP_EB","cutflow (HEEP EB)",16,-1.,15.);
  histo1d_["cutflow_HEEP_EE"] = fs->make<TH1D>("cutflow_HEEP_EE","cutflow (HEEP EE)",16,-1.,15.);
  histo1d_["status_mergedElectron_EB"] = fs->make<TH1D>("status_mergedElectron_EB","ME status (EB)",40,-1.,39.);
  histo1d_["status_mergedElectron_EE"] = fs->make<TH1D>("status_mergedElectron_EE","ME status (EE)",40,-1.,39.);
  histo1d_["mva_DR1Et2EB"] = fs->make<TH1D>("mva_DR1Et2EB","MVA score",100,0.,1.);
  histo1d_["mva_DR2Et1EB"] = fs->make<TH1D>("mva_DR2Et1EB","MVA score",100,0.,1.);
  histo1d_["mva_DR2Et2EB"] = fs->make<TH1D>("mva_DR2Et2EB","MVA score",100,0.,1.);
  histo1d_["mva_DR2Et1EE"] = fs->make<TH1D>("mva_DR2Et1EE","MVA score",100,0.,1.);
  histo1d_["mva_DR2Et2EE"] = fs->make<TH1D>("mva_DR2Et2EE","MVA score",100,0.,1.);
  histo1d_["mva_bkgEt2EB"] = fs->make<TH1D>("mva_bkgEt2EB","MVA score",100,0.,1.);
  histo1d_["mva_DR1Et1EB_extended"] = fs->make<TH1D>("mva_DR1Et1EB_extended","MVA score",100,0.,1.);
  histo1d_["mva_bkgEt1EB_extended"] = fs->make<TH1D>("mva_bkgEt1EB_extended","MVA score",100,0.,1.);
  histo1d_["cutflow_4E"] = fs->make<TH1D>("cutflow_4E","cutflow (4E)",10,0.,10.);
  histo1d_["cutflow_3E"] = fs->make<TH1D>("cutflow_3E","cutflow (3E)",10,0.,10.);
  histo1d_["cutflow_2E"] = fs->make<TH1D>("cutflow_2E","cutflow (2E)",10,0.,10.);

  histo2d_["invM_dr"] = fs->make<TH2D>("invM_dr","M(ll) vs dR(ll)",200,0.,1000.,128,0.,6.4);
  histo2d_["invM_dr_noGsf"] = fs->make<TH2D>("invM_dr_noGsf","M(ll) vs dR(ll)",200,0.,1000.,128,0.,6.4);
  histo2d_["invM_dr_mixed"] = fs->make<TH2D>("invM_dr_mixed","M(ll) vs dR(ll)",200,0.,1000.,128,0.,6.4);

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

  histo1d_["2E_CRME_ll_invM"] = fs->make<TH1D>("2E_CRME_ll_invM","CR M(ll);GeV;",40,0.,200.);
  histo1d_["2E_CRME_ll_pt"] = fs->make<TH1D>("2E_CRME_ll_pt","CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_ll_dr"] = fs->make<TH1D>("2E_CRME_ll_dr","CR dR(ll)",64,0.,6.4);

  histo1d_["2E_CRME_SSll_invM"] = fs->make<TH1D>("2E_CRME_SSll_invM","SS CR M(ll);GeV;",40,0.,200.);
  histo1d_["2E_CRME_SSll_pt"] = fs->make<TH1D>("2E_CRME_SSll_pt","SS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_SSll_dr"] = fs->make<TH1D>("2E_CRME_SSll_dr","SS CR dR(ll)",64,0.,6.4);

  histo1d_["2E_CRME_OSll_invM"] = fs->make<TH1D>("2E_CRME_OSll_invM","OS CR M(ll);GeV;",40,0.,200.);
  histo1d_["2E_CRME_OSll_pt"] = fs->make<TH1D>("2E_CRME_OSll_pt","OS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_CRME_OSll_dr"] = fs->make<TH1D>("2E_CRME_OSll_dr","OS CR dR(ll)",64,0.,6.4);

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

  histo1d_["2E_mixedME_ll_invM"] = fs->make<TH1D>("2E_mixedME_ll_invM","mixed CR M(ll);GeV;",80,0.,200.);
  histo1d_["2E_mixedME_ll_pt"] = fs->make<TH1D>("2E_mixedME_ll_pt","mixed CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_ll_dr"] = fs->make<TH1D>("2E_mixedME_ll_dr","mixed CR dR(ll)",64,0.,6.4);

  histo1d_["2E_mixedME_SSll_invM"] = fs->make<TH1D>("2E_mixedME_SSll_invM","mixed SS CR M(ll);GeV;",80,0.,200.);
  histo1d_["2E_mixedME_SSll_invM_xFF"] = fs->make<TH1D>("2E_mixedME_SSll_invM_xFF","mixed SS CR M(ll) x Fake factor;GeV;",80,0.,200.);
  histo1d_["2E_mixedME_SSll_pt"] = fs->make<TH1D>("2E_mixedME_SSll_pt","mixed SS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_SSll_dr"] = fs->make<TH1D>("2E_mixedME_SSll_dr","mixed SS CR dR(ll)",64,0.,6.4);

  histo1d_["2E_mixedME_OSll_invM"] = fs->make<TH1D>("2E_mixedME_OSll_invM","mixed OS CR M(ll);GeV;",80,0.,200.);
  histo1d_["2E_mixedME_OSll_invM_xFF"] = fs->make<TH1D>("2E_mixedME_OSll_invM_xFF","mixed OS CR M(ll) x Fake factor;GeV;",80,0.,200.);
  histo1d_["2E_mixedME_OSll_pt"] = fs->make<TH1D>("2E_mixedME_OSll_pt","mixed OS CR p_{T}(ll);GeV;",100,0.,500.);
  histo1d_["2E_mixedME_OSll_dr"] = fs->make<TH1D>("2E_mixedME_OSll_dr","mixed OS CR dR(ll)",64,0.,6.4);

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
  histo1d_["2E_antiME_SSll_pt"] = fs->make<TH1D>("2E_antiME_SSll_pt","SS anti p_{T}(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_SSll_dr"] = fs->make<TH1D>("2E_antiME_SSll_dr","SS anti dR(ll)",128,0.,6.4);

  histo1d_["2E_antiME_OSll_invM"] = fs->make<TH1D>("2E_antiME_OSll_invM","OS anti M(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_invM_CR"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR","OS anti M(ll);GeV;",200,0.,200.);
  histo1d_["2E_antiME_OSll_invM_CR_xFF"] = fs->make<TH1D>("2E_antiME_OSll_invM_CR_xFF","OS anti M(ll) x Fake factor;GeV;",200,0.,200.);
  histo1d_["2E_antiME_OSll_pt"] = fs->make<TH1D>("2E_antiME_OSll_pt","OS anti p_{T}(ll);GeV;",200,0.,1000.);
  histo1d_["2E_antiME_OSll_dr"] = fs->make<TH1D>("2E_antiME_OSll_dr","OS anti dR(ll)",128,0.,6.4);

  histo1d_["2E_Et_SSCR_antiME"] = fs->make<TH1D>("2E_Et_SSCR_antiME","SSCR denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_CRME"] = fs->make<TH1D>("2E_Et_SSCR_CRME","SSCR numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_mixedME"] = fs->make<TH1D>("2E_Et_SSCR_mixedME","SSCR (mixed) numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_SSCR_mixedAntiME"] = fs->make<TH1D>("2E_Et_SSCR_mixedAntiME","SSCR (mixed) denominator E_{T};GeV;",200,0.,1000.);

  histo1d_["2E_Et_OSCR_antiME"] = fs->make<TH1D>("2E_Et_OSCR_antiME","OSCR denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_CRME"] = fs->make<TH1D>("2E_Et_OSCR_CRME","OSCR numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_mixedME"] = fs->make<TH1D>("2E_Et_OSCR_mixedME","OSCR (mixed) numerator E_{T};GeV;",200,0.,1000.);
  histo1d_["2E_Et_OSCR_mixedAntiME"] = fs->make<TH1D>("2E_Et_OSCR_mixedAntiME","OSCR (mixed) denominator E_{T};GeV;",200,0.,1000.);

  histo1d_["3E_CRME_Et"] = fs->make<TH1D>("3E_CRME_Et","CR ME E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CRME_eta"] = fs->make<TH1D>("3E_CRME_eta","CR ME #eta",50,-2.5,2.5);
  histo1d_["3E_CRME_phi"] = fs->make<TH1D>("3E_CRME_phi","CR ME #phi",64,-3.2,3.2);

  histo1d_["3E_CRME_noGsf_Et"] = fs->make<TH1D>("3E_CRME_noGsf_Et","CR ME (w/o GSF) E_{T};GeV;",40,0.,200.);
  histo1d_["3E_CRME_noGsf_eta"] = fs->make<TH1D>("3E_CRME_noGsf_eta","CR ME (w/o GSF) #eta",50,-2.5,2.5);
  histo1d_["3E_CRME_noGsf_phi"] = fs->make<TH1D>("3E_CRME_noGsf_phi","CR ME (w/o GSF) #phi",64,-3.2,3.2);

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
  histo1d_["3E_antiME_lll_pt"] = fs->make<TH1D>("3E_antiME_lll_pt","p_{T}(lll);GeV;",200,0.,1000.);
  histo1d_["3E_antiME_lll_dr_l1ME"] = fs->make<TH1D>("3E_antiME_lll_dr_l1ME","dR(l1ME)",128,0.,6.4);
  histo1d_["3E_antiME_lll_dr_l2ME"] = fs->make<TH1D>("3E_antiME_lll_dr_l2ME","dR(l2ME)",128,0.,6.4);
  histo1d_["3E_antiME_lll_dr_l1l2"] = fs->make<TH1D>("3E_antiME_lll_dr_l1l2","dR(l1l2)",128,0.,6.4);

  histo1d_["3E_Et_CR_antiME"] = fs->make<TH1D>("3E_Et_CR_antiME","3EVR denominator E_{T};GeV;",200,0.,1000.);
  histo1d_["3E_Et_CR_CRME"] = fs->make<TH1D>("3E_Et_CR_CRME","3EVR numerator E_{T};GeV;",200,0.,1000.);
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

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

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

  edm::Handle<edm::ValueMap<int>> cutflow_modifiedHEEPHandle;
  iEvent.getByToken(cutflow_modifiedHEEPToken_, cutflow_modifiedHEEPHandle);

  edm::Handle<edm::ValueMap<int>> cutflow_HEEPHandle;
  iEvent.getByToken(cutflow_HEEPToken_, cutflow_HEEPHandle);

  edm::Handle<edm::ValueMap<int>> status_mergedElectronHandle;
  iEvent.getByToken(status_mergedElectronToken_, status_mergedElectronHandle);

  edm::Handle<edm::ValueMap<float>> mva_mergedElectronHandle;
  iEvent.getByToken(mva_mergedElectronToken_, mva_mergedElectronHandle);

  edm::Handle<edm::ValueMap<int>> GSFtype_mergedElectronHandle;
  iEvent.getByToken(GSFtype_mergedElectronToken_, GSFtype_mergedElectronHandle);

  // first two (modified)HEEP electrons should be Et > 35 GeV
  std::vector<pat::ElectronRef> acceptEles;

  auto sortByEt = [](const pat::ElectronRef& a, const pat::ElectronRef& b) {
    return a->et() > b->et();
  };

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    const auto& aEle = eleHandle->refAt(idx);
    const auto& orgGsfTrk = aEle->gsfTrack();
    const auto& addGsfTrk = (*addGsfTrkMap)[aEle];

    if ( std::abs(aEle->eta()) > 2.5 )
      continue;

    // veto EBEE gap
    if ( std::abs(aEle->eta()) > 1.44 && std::abs(aEle->eta()) < 1.57 )
      continue;

    std::string strModHEEP = "cutflow_modifiedHEEP";
    std::string strHEEP = "cutflow_HEEP";

    if (aEle->isEB()) {
      strModHEEP += "_EB";
      strHEEP += "_EB";
    } else {
      strModHEEP += "_EE";
      strHEEP += "_EE";
    }

    MergedLeptonIDs::fillCutflow( histo1d_[strModHEEP], (*cutflow_modifiedHEEPHandle)[aEle], aWeight );
    MergedLeptonIDs::fillCutflow( histo1d_[strHEEP], (*cutflow_HEEPHandle)[aEle] , aWeight );

    if ( MergedLeptonIDs::isSameGsfTrack(addGsfTrk,orgGsfTrk) ) {
      if ( (*cutflow_HEEPHandle)[aEle] == static_cast<int>(MergedLeptonIDs::cutflowElectron::showerShape) ) {
        auto castEle = aEle.castTo<pat::ElectronRef>();
        acceptEles.push_back(castEle);
      }
    } else {
      if ( (*cutflow_modifiedHEEPHandle)[aEle] == static_cast<int>(MergedLeptonIDs::cutflowElectron::dxy) ) {
        auto castEle = aEle.castTo<pat::ElectronRef>();
        acceptEles.push_back(castEle);
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

  std::vector<pat::ElectronRef> mergedEls;
  std::vector<pat::ElectronRef> CRMEs;
  std::vector<pat::ElectronRef> antiMEs;

  for (unsigned int idx = 0; idx < acceptEles.size(); ++idx) {
    auto aEle = acceptEles.at(idx);

    auto aGSFtype = static_cast<MergedLeptonIDs::GSFtype>((*GSFtype_mergedElectronHandle)[aEle]);
    auto aStatus = static_cast<MergedLeptonIDs::cutflowElectron>((*status_mergedElectronHandle)[aEle]);

    double drGSF = std::sqrt(reco::deltaR2(aEle->gsfTrack()->eta(),
                                           aEle->gsfTrack()->phi(),
                                           (*addGsfTrkMap)[aEle]->eta(),
                                           (*addGsfTrkMap)[aEle]->phi()));

    std::string postfix = "";
    std::string strEBEE = "";

    if (aEle->isEB())
      strEBEE = "_EB";
    else
      strEBEE = "_EE";

    if ( aEle->et() < 50. )
      continue;

    if ( aStatus==MergedLeptonIDs::cutflowElectron::no2ndGsf ||
         aStatus==MergedLeptonIDs::cutflowElectron::passedMVA2 ||
         aStatus==MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail ||
         aStatus==MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass ) {
      // no 2nd GSF
      postfix = "_bkgEt2EB";

      if ( aStatus==MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail ||
           aStatus==MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass )
        postfix = "_bkgEt1EB_extended";

      histo1d_["et_noGsf_HEEP"]->Fill( aEle->et(), aWeight );
      histo1d_[static_cast<std::string>("mva")+postfix]
        ->Fill( (*mva_mergedElectronHandle)[aEle] , aWeight );

    } else if ( aStatus==MergedLeptonIDs::cutflowElectron::has2ndGsf ||
                aStatus==MergedLeptonIDs::cutflowElectron::passedMVA1 ||
                aStatus==MergedLeptonIDs::cutflowElectron::has2ndGsfCRfail ||
                aStatus==MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass ) {
      // has 2nd GSF
      histo2d_["et_dr_HEEP"]->Fill( aEle->et(), drGSF, aWeight );

      switch (aGSFtype) { // better way to do this?
        case MergedLeptonIDs::GSFtype::DR1Et2EB:
          postfix = "_DR1Et2EB";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et1EB:
          postfix = "_DR2Et1EB";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et2EB:
          postfix = "_DR2Et2EB";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et1EE:
          postfix = "_DR2Et1EE";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et2EE:
          postfix = "_DR2Et2EE";
          break;
        case MergedLeptonIDs::GSFtype::extendedCR:
          postfix = "_DR1Et1EB_extended";
          break;
        default:
          throw cms::Exception("LogicError") << "opening angle between the original and additional GSF track does not fit into any of MergedLeptonIDs::openingAngle categories" << std::endl;
          break;
      }

      histo1d_[static_cast<std::string>("mva")+postfix]
        ->Fill( (*mva_mergedElectronHandle)[aEle] , aWeight );
    } else {}

    histo1d_[static_cast<std::string>("status_mergedElectron")+strEBEE]
      ->Fill( static_cast<float>( (*status_mergedElectronHandle)[aEle] ) + 0.5 , aWeight );

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) )
      mergedEls.push_back(aEle);

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass) )
      histo2d_["et_dr_final"]->Fill( aEle->et(), drGSF, aWeight );

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass) )
      histo1d_["et_noGsf_final"]->Fill( aEle->et(), aWeight );

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass) )
      CRMEs.push_back(aEle);

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsf) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsf) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRfail) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail) )
      antiMEs.push_back(aEle);

  } // electron loop

  std::sort(mergedEls.begin(),mergedEls.end(),sortByEt);
  std::sort(CRMEs.begin(),CRMEs.end(),sortByEt);
  std::sort(antiMEs.begin(),antiMEs.end(),sortByEt);

  double objwgt = 1.;

  if (isMC_) {
    for (const auto& aEle : acceptEles) {
      double apt = aEle->pt() >= 500. ? 499.999 : aEle->pt();
      objwgt *= recoSF_->GetBinContent( recoSF_->FindBin(aEle->eta(), apt ) );
    }
  }

  aWeight *= objwgt;

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
          if ( (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
               (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass) ) {
            histo1d_["3E_CRME_Et"]->Fill(CRMEs.front()->et(),aWeight);
            histo1d_["3E_CRME_eta"]->Fill(CRMEs.front()->eta(),aWeight);
            histo1d_["3E_CRME_phi"]->Fill(CRMEs.front()->phi(),aWeight);
          } else if ( (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ||
                      (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass) ) {
            histo1d_["3E_CRME_noGsf_Et"]->Fill(CRMEs.front()->et(),aWeight);
            histo1d_["3E_CRME_noGsf_eta"]->Fill(CRMEs.front()->eta(),aWeight);
            histo1d_["3E_CRME_noGsf_phi"]->Fill(CRMEs.front()->phi(),aWeight);
          } else {}

          histo1d_["3E_CR_NME1_Et"]->Fill(nonMEs.front()->et(),aWeight);
          histo1d_["3E_CR_NME1_eta"]->Fill(nonMEs.front()->eta(),aWeight);
          histo1d_["3E_CR_NME1_phi"]->Fill(nonMEs.front()->phi(),aWeight);

          histo1d_["3E_CR_NME2_Et"]->Fill(nonMEs.at(1)->et(),aWeight);
          histo1d_["3E_CR_NME2_eta"]->Fill(nonMEs.at(1)->eta(),aWeight);
          histo1d_["3E_CR_NME2_phi"]->Fill(nonMEs.at(1)->phi(),aWeight);

          histo1d_["3E_CRME_lll_invM"]->Fill(mlll,aWeight);
          histo1d_["3E_CRME_lll_pt"]->Fill(lvecCRllME.pt(),aWeight);
          histo1d_["3E_CRME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
          histo1d_["3E_CRME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
          histo1d_["3E_CRME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

          histo1d_["3E_Et_CR_CRME"]->Fill( CRMEs.front()->et(), aWeight );
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
          if ( (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsf) ||
               (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRfail) ) {
            histo1d_["3E_antiME_Et"]->Fill(antiMEs.front()->et(),aWeight);
            histo1d_["3E_antiME_eta"]->Fill(antiMEs.front()->eta(),aWeight);
            histo1d_["3E_antiME_phi"]->Fill(antiMEs.front()->phi(),aWeight);
          } else if ( (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsf) ||
                      (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail) ) {
            histo1d_["3E_antiME_noGsf_Et"]->Fill(antiMEs.front()->et(),aWeight);
            histo1d_["3E_antiME_noGsf_eta"]->Fill(antiMEs.front()->eta(),aWeight);
            histo1d_["3E_antiME_noGsf_phi"]->Fill(antiMEs.front()->phi(),aWeight);
          } else {}

          histo1d_["3E_antiME_NME1_Et"]->Fill(nonMEs.front()->et(),aWeight);
          histo1d_["3E_antiME_NME1_eta"]->Fill(nonMEs.front()->eta(),aWeight);
          histo1d_["3E_antiME_NME1_phi"]->Fill(nonMEs.front()->phi(),aWeight);

          histo1d_["3E_antiME_NME2_Et"]->Fill(nonMEs.at(1)->et(),aWeight);
          histo1d_["3E_antiME_NME2_eta"]->Fill(nonMEs.at(1)->eta(),aWeight);
          histo1d_["3E_antiME_NME2_phi"]->Fill(nonMEs.at(1)->phi(),aWeight);

          histo1d_["3E_antiME_lll_invM"]->Fill(mlll,aWeight);
          histo1d_["3E_antiME_lll_pt"]->Fill(lvecAllME.pt(),aWeight);
          histo1d_["3E_antiME_lll_dr_l1ME"]->Fill(std::sqrt(dr2l1ME),aWeight);
          histo1d_["3E_antiME_lll_dr_l2ME"]->Fill(std::sqrt(dr2l2ME),aWeight);
          histo1d_["3E_antiME_lll_dr_l1l2"]->Fill(std::sqrt(dr2l1l2),aWeight);

          if ( mlll < 200. ) {
            const double aEt = antiMEs.front()->et();
            const double ffOS = ( aEt > osBoundary_ ) ? oshigh_->Eval(aEt) : oslow_->Eval(aEt);
            const double ffSS = ( aEt > ssBoundary_ ) ? sshigh_->Eval(aEt) : sslow_->Eval(aEt);
            histo1d_["3E_antiME_lll_invM_CR"]->Fill(mlll,aWeight);
            histo1d_["3E_antiME_lll_invM_CR_xOSFF"]->Fill(mlll,aWeight*ffOS);
            histo1d_["3E_antiME_lll_invM_CR_xSSFF"]->Fill(mlll,aWeight*ffSS);
            histo1d_["3E_Et_CR_antiME"]->Fill( antiMEs.front()->et(), aWeight );
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

      if (CRMEs.size()==2) {
        const auto lvecCRME1 = CRMEs.front()->polarP4()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy();
        const auto lvecCRME2 = CRMEs.at(1)->polarP4()*CRMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/CRMEs.at(1)->energy();
        const auto lvecCRll = lvecCRME1 + lvecCRME2;
        const double dr2CRll = reco::deltaR2(lvecCRME1.eta(),lvecCRME1.phi(),lvecCRME2.eta(),lvecCRME2.phi());
        const double mll = lvecCRll.M();

        if ( mll > 50. && mll < 200. ) {
          if ( (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
               (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass) ) {
            histo1d_["2E_CRME1_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_CRME1_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_CRME1_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else if ( (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ||
                      (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass) ) {
            histo1d_["2E_CRME1_noGsf_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_CRME1_noGsf_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_CRME1_noGsf_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else {}

          if ( (*status_mergedElectronHandle)[CRMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
               (*status_mergedElectronHandle)[CRMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass) ) {
            histo1d_["2E_CRME2_Et"]->Fill( CRMEs.at(1)->et()*CRMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/CRMEs.at(1)->energy(), aWeight );
            histo1d_["2E_CRME2_eta"]->Fill( CRMEs.at(1)->eta(), aWeight );
            histo1d_["2E_CRME2_phi"]->Fill( CRMEs.at(1)->phi(), aWeight );
          } else if ( (*status_mergedElectronHandle)[CRMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ||
                      (*status_mergedElectronHandle)[CRMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass) ) {
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

            for (const auto& aME : CRMEs)
              histo1d_["2E_Et_OSCR_CRME"]->Fill( aME->et(), aWeight );
          } else if ( CRMEs.front()->charge()*CRMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_CRME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_CRME_SSll_pt"]->Fill( lvecCRll.pt(), aWeight );
            histo1d_["2E_CRME_SSll_dr"]->Fill( std::sqrt(dr2CRll), aWeight );

            for (const auto& aME : CRMEs)
              histo1d_["2E_Et_SSCR_CRME"]->Fill( aME->et(), aWeight );
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
          if ( (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
               (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRpass) ) {
            histo1d_["2E_mixed_ME_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_ME_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_ME_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else if ( (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ||
                      (*status_mergedElectronHandle)[CRMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRpass) ) {
            histo1d_["2E_mixed_ME_noGsf_Et"]->Fill( CRMEs.front()->et()*CRMEs.front()->userFloat("ecalTrkEnergyPostCorr")/CRMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_ME_noGsf_eta"]->Fill( CRMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_ME_noGsf_phi"]->Fill( CRMEs.front()->phi(), aWeight );
          } else {}

          if ( (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsf) ||
               (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRfail) ) {
            histo1d_["2E_mixed_antiME_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_antiME_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_antiME_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else if ( (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsf) ||
                      (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail) ) {
            histo1d_["2E_mixed_antiME_noGsf_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_mixed_antiME_noGsf_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_mixed_antiME_noGsf_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else {}

          histo1d_["2E_mixedME_ll_invM"]->Fill( mll, aWeight );
          histo1d_["2E_mixedME_ll_pt"]->Fill( lvecMixedll.pt(), aWeight );
          histo1d_["2E_mixedME_ll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );

          const double aEt = antiMEs.front()->et();
          const double ffOS = ( aEt > osBoundary_ ) ? oshigh_->Eval(aEt) : oslow_->Eval(aEt);
          const double ffSS = ( aEt > ssBoundary_ ) ? sshigh_->Eval(aEt) : sslow_->Eval(aEt);

          if ( CRMEs.front()->charge()*antiMEs.front()->charge() < 0. ) { // OS
            histo1d_["2E_mixedME_OSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_mixedME_OSll_invM_xFF"]->Fill( mll, 0.5*aWeight*ffOS ); // double-counting
            histo1d_["2E_mixedME_OSll_pt"]->Fill( lvecMixedll.pt(), aWeight );
            histo1d_["2E_mixedME_OSll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );
            histo1d_["2E_Et_OSCR_mixedME"]->Fill( CRMEs.front()->et(), aWeight );
            histo1d_["2E_Et_OSCR_mixedAntiME"]->Fill( antiMEs.front()->et(), aWeight );
          } else if ( CRMEs.front()->charge()*antiMEs.front()->charge() > 0. ) { // SS
            histo1d_["2E_mixedME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_mixedME_SSll_invM_xFF"]->Fill( mll, 0.5*aWeight*ffSS ); // double-counting
            histo1d_["2E_mixedME_SSll_pt"]->Fill( lvecMixedll.pt(), aWeight );
            histo1d_["2E_mixedME_SSll_dr"]->Fill( std::sqrt(dr2mixedll), aWeight );
            histo1d_["2E_Et_SSCR_mixedME"]->Fill( CRMEs.front()->et(), aWeight );
            histo1d_["2E_Et_SSCR_mixedAntiME"]->Fill( antiMEs.front()->et(), aWeight );
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
          if ( (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsf) ||
               (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRfail) ) {
            histo1d_["2E_antiME1_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_antiME1_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_antiME1_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else if ( (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsf) ||
                      (*status_mergedElectronHandle)[antiMEs.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail) ) {
            histo1d_["2E_antiME1_noGsf_Et"]->Fill( antiMEs.front()->et()*antiMEs.front()->userFloat("ecalTrkEnergyPostCorr")/antiMEs.front()->energy(), aWeight );
            histo1d_["2E_antiME1_noGsf_eta"]->Fill( antiMEs.front()->eta(), aWeight );
            histo1d_["2E_antiME1_noGsf_phi"]->Fill( antiMEs.front()->phi(), aWeight );
          } else {}

          if ( (*status_mergedElectronHandle)[antiMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsf) ||
               (*status_mergedElectronHandle)[antiMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::has2ndGsfCRfail) ) {
            histo1d_["2E_antiME2_Et"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight );
            histo1d_["2E_antiME2_eta"]->Fill( antiMEs.at(1)->eta(), aWeight );
            histo1d_["2E_antiME2_phi"]->Fill( antiMEs.at(1)->phi(), aWeight );
          } else if ( (*status_mergedElectronHandle)[antiMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsf) ||
                      (*status_mergedElectronHandle)[antiMEs.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::no2ndGsfCRfail) ) {
            histo1d_["2E_antiME2_noGsf_Et"]->Fill( antiMEs.at(1)->et()*antiMEs.at(1)->userFloat("ecalTrkEnergyPostCorr")/antiMEs.at(1)->energy(), aWeight );
            histo1d_["2E_antiME2_noGsf_eta"]->Fill( antiMEs.at(1)->eta(), aWeight );
            histo1d_["2E_antiME2_noGsf_phi"]->Fill( antiMEs.at(1)->phi(), aWeight );
          } else {}

          histo1d_["2E_antiME_ll_invM"]->Fill( mll, aWeight );
          histo1d_["2E_antiME_ll_pt"]->Fill( lvecAll.pt(), aWeight );
          histo1d_["2E_antiME_ll_dr"]->Fill( std::sqrt(dr2All), aWeight );

          const double aEt1 = antiMEs.front()->et();
          const double ffOS1 = ( aEt1 > osBoundary_ ) ? oshigh_->Eval(aEt1) : oslow_->Eval(aEt1);
          const double ffSS1 = ( aEt1 > ssBoundary_ ) ? sshigh_->Eval(aEt1) : sslow_->Eval(aEt1);
          const double aEt2 = antiMEs.at(1)->et();
          const double ffOS2 = ( aEt2 > osBoundary_ ) ? oshigh_->Eval(aEt2) : oslow_->Eval(aEt2);
          const double ffSS2 = ( aEt2 > ssBoundary_ ) ? sshigh_->Eval(aEt2) : sslow_->Eval(aEt2);

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

              for (const auto& aME : antiMEs)
                histo1d_["2E_Et_OSCR_antiME"]->Fill( aME->et(), aWeight );
            }
          } else if ( antiMEs.front()->charge()*antiMEs.at(1)->charge() > 0. ) { // SS
            histo1d_["2E_antiME_SSll_invM"]->Fill( mll, aWeight );
            histo1d_["2E_antiME_SSll_pt"]->Fill( lvecAll.pt(), aWeight );
            histo1d_["2E_antiME_SSll_dr"]->Fill( std::sqrt(dr2All), aWeight );

            if (mll < 200.) {
              histo1d_["2E_antiME_SSll_invM_CR"]->Fill( mll, aWeight );
              histo1d_["2E_antiME_SSll_invM_CR_xFF"]->Fill( mll, aWeight*ffSS1 );
              histo1d_["2E_antiME_SSll_invM_CR_xFF"]->Fill( mll, aWeight*ffSS2 );

              for (const auto& aME : antiMEs)
                histo1d_["2E_Et_SSCR_antiME"]->Fill( aME->et(), aWeight );
            } // CR
          } // OSSS
        } // mll > 50.
      } // antiME.size()==2

      break;
  } // switch acceptEles.size()
}

DEFINE_FWK_MODULE(MergedEleCRanalyzer);
