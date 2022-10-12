#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
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

class MergedEleCRanalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleCRanalyzer(const edm::ParameterSet&);
  virtual ~MergedEleCRanalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  bool isMC_;

  const edm::EDGetTokenT<edm::View<reco::GsfElectron>> srcEle_;
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

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;
};

MergedEleCRanalyzer::MergedEleCRanalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
srcEle_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
cutflow_modifiedHEEPToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("cutflow_modifiedHEEP"))),
cutflow_HEEPToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("cutflow_HEEP"))),
status_mergedElectronToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("status_mergedElectron"))),
mva_mergedElectronToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mva_mergedElectron"))),
GSFtype_mergedElectronToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("GSFtype_mergedElectron"))) {
  usesResource("TFileService");
}

void MergedEleCRanalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow_modifiedHEEP_EB"] = fs->make<TH1D>("cutflow_modifiedHEEP_EB","cutflow",11,-1.,10.);
  histo1d_["cutflow_modifiedHEEP_EE"] = fs->make<TH1D>("cutflow_modifiedHEEP_EE","cutflow",11,-1.,10.);
  histo1d_["cutflow_HEEP_EB"] = fs->make<TH1D>("cutflow_HEEP_EB","cutflow",16,-1.,15.);
  histo1d_["cutflow_HEEP_EE"] = fs->make<TH1D>("cutflow_HEEP_EE","cutflow",16,-1.,15.);
  histo1d_["status_mergedElectron_DR1Et2EB"] = fs->make<TH1D>("status_mergedElectron_DR1Et2EB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR2Et1EB"] = fs->make<TH1D>("status_mergedElectron_DR2Et1EB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR2Et2EB"] = fs->make<TH1D>("status_mergedElectron_DR2Et2EB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR1Et2EE"] = fs->make<TH1D>("status_mergedElectron_DR1Et2EE","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR2Et1EE"] = fs->make<TH1D>("status_mergedElectron_DR2Et1EE","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR2Et2EE"] = fs->make<TH1D>("status_mergedElectron_DR2Et2EE","status",30,-1.,29.);
  histo1d_["status_mergedElectron_bkgEt2EB"] = fs->make<TH1D>("status_mergedElectron_bkgEt2EB","status",30,-1.,29.);
  histo1d_["mva_DR1Et2EB"] = fs->make<TH1D>("mva_DR1Et2EB","MVA score",200,-1.,1.);
  histo1d_["mva_DR2Et1EB"] = fs->make<TH1D>("mva_DR2Et1EB","MVA score",200,-1.,1.);
  histo1d_["mva_DR2Et2EB"] = fs->make<TH1D>("mva_DR2Et2EB","MVA score",200,-1.,1.);
  histo1d_["mva_DR1Et2EE"] = fs->make<TH1D>("mva_DR1Et2EE","MVA score",200,-1.,1.);
  histo1d_["mva_DR2Et1EE"] = fs->make<TH1D>("mva_DR2Et1EE","MVA score",200,-1.,1.);
  histo1d_["mva_DR2Et2EE"] = fs->make<TH1D>("mva_DR2Et2EE","MVA score",200,-1.,1.);
  histo1d_["mva_bkgEt2EB"] = fs->make<TH1D>("mva_bkgEt2EB","MVA score",200,-1.,1.);
  histo1d_["cutflow_4E"] = fs->make<TH1D>("cutflow_4E","cutflow",10,0.,10.);
  histo1d_["cutflow_3E"] = fs->make<TH1D>("cutflow_3E","cutflow",10,0.,10.);
  histo1d_["cutflow_2E"] = fs->make<TH1D>("cutflow_2E","cutflow",10,0.,10.);

  histo1d_["pt_3E_ME"] = fs->make<TH1D>("pt_3E_ME","Pt",200,0.,500.);
  histo1d_["eta_3E_ME"] = fs->make<TH1D>("eta_3E_ME","3Eta",200,-2.5,2.5);
  histo1d_["phi_3E_ME"] = fs->make<TH1D>("phi_3E_ME","Phi",128,-3.2,3.2);

  histo1d_["pt_2E_ME1"] = fs->make<TH1D>("pt_2E_ME1","Pt",200,0.,500.);
  histo1d_["eta_2E_ME1"] = fs->make<TH1D>("eta_2E_ME1","Eta",200,-2.5,2.5);
  histo1d_["phi_2E_ME1"] = fs->make<TH1D>("phi_2E_ME1","Phi",128,-3.2,3.2);

  histo1d_["pt_2E_ME2"] = fs->make<TH1D>("pt_2E_ME2","Pt",200,0.,500.);
  histo1d_["eta_2E_ME2"] = fs->make<TH1D>("eta_2E_ME2","Eta",200,-2.5,2.5);
  histo1d_["phi_2E_ME2"] = fs->make<TH1D>("phi_2E_ME2","Phi",128,-3.2,3.2);

  histo1d_["invM_2E_ll"] = fs->make<TH1D>("invM_2E_ll","M(ll)",200,0.,500.);
  histo1d_["rap_2E_ll"] = fs->make<TH1D>("rap_2E_ll","rapidity",200,-2.5,2.5);
  histo1d_["pt_2E_ll"] = fs->make<TH1D>("pt_2E_ll","Pt(ll)",200,0.,500.);
  histo1d_["dr_2E_ll"] = fs->make<TH1D>("dr_2E_ll","dR(ll)",128,0.,6.4);
  histo1d_["deta_2E_ll"] = fs->make<TH1D>("deta_2E_ll","dEta(ll)",200,-2.5,2.5);
  histo1d_["dphi_2E_ll"] = fs->make<TH1D>("dphi_2E_ll","dPhi(ll)",128,-3.2,3.2);

  histo1d_["pt_3E_ME_noGsf"] = fs->make<TH1D>("pt_3E_ME_noGsf","Pt",200,0.,500.);
  histo1d_["eta_3E_ME_noGsf"] = fs->make<TH1D>("eta_3E_ME_noGsf","3Eta",200,-2.5,2.5);
  histo1d_["phi_3E_ME_noGsf"] = fs->make<TH1D>("phi_3E_ME_noGsf","Phi",128,-3.2,3.2);

  histo1d_["pt_2E_ME1_noGsf"] = fs->make<TH1D>("pt_2E_ME1_noGsf","Pt",200,0.,500.);
  histo1d_["eta_2E_ME1_noGsf"] = fs->make<TH1D>("eta_2E_ME1_noGsf","Eta",200,-2.5,2.5);
  histo1d_["phi_2E_ME1_noGsf"] = fs->make<TH1D>("phi_2E_ME1_noGsf","Phi",128,-3.2,3.2);

  histo1d_["pt_2E_ME2_noGsf"] = fs->make<TH1D>("pt_2E_ME2_noGsf","Pt",200,0.,500.);
  histo1d_["eta_2E_ME2_noGsf"] = fs->make<TH1D>("eta_2E_ME2_noGsf","Eta",200,-2.5,2.5);
  histo1d_["phi_2E_ME2_noGsf"] = fs->make<TH1D>("phi_2E_ME2_noGsf","Phi",128,-3.2,3.2);

  histo1d_["invM_2E_ll_noGsf"] = fs->make<TH1D>("invM_2E_ll_noGsf","M(ll)",200,0.,500.);
  histo1d_["rap_2E_ll_noGsf"] = fs->make<TH1D>("rap_2E_ll_noGsf","rapidity",200,-2.5,2.5);
  histo1d_["pt_2E_ll_noGsf"] = fs->make<TH1D>("pt_2E_ll_noGsf","Pt(ll)",200,0.,500.);
  histo1d_["dr_2E_ll_noGsf"] = fs->make<TH1D>("dr_2E_ll_noGsf","dR(ll)",128,0.,6.4);
  histo1d_["deta_2E_ll_noGsf"] = fs->make<TH1D>("deta_2E_ll_noGsf","dEta(ll)",200,-2.5,2.5);
  histo1d_["dphi_2E_ll_noGsf"] = fs->make<TH1D>("dphi_2E_ll_noGsf","dPhi(ll)",128,-3.2,3.2);

  histo1d_["invM_2E_ll_mixed"] = fs->make<TH1D>("invM_2E_ll_mixed","M(ll)",200,0.,500.);
  histo1d_["rap_2E_ll_mixed"] = fs->make<TH1D>("rap_2E_ll_mixed","rapidity",200,-2.5,2.5);
  histo1d_["pt_2E_ll_mixed"] = fs->make<TH1D>("pt_2E_ll_mixed","Pt(ll)",200,0.,500.);
  histo1d_["dr_2E_ll_mixed"] = fs->make<TH1D>("dr_2E_ll_mixed","dR(ll)",128,0.,6.4);
  histo1d_["deta_2E_ll_mixed"] = fs->make<TH1D>("deta_2E_ll_mixed","dEta(ll)",200,-2.5,2.5);
  histo1d_["dphi_2E_ll_mixed"] = fs->make<TH1D>("dphi_2E_ll_mixed","dPhi(ll)",128,-3.2,3.2);

  histo2d_["invM_dr"] = fs->make<TH2D>("invM_dr","M(ll) vs dR(ll)",200,0.,500.,128,0.,6.4);
  histo2d_["invM_dr_noGsf"] = fs->make<TH2D>("invM_dr_noGsf","M(ll) vs dR(ll)",200,0.,500.,128,0.,6.4);
  histo2d_["invM_dr_mixed"] = fs->make<TH2D>("invM_dr_mixed","M(ll) vs dR(ll)",200,0.,500.,128,0.,6.4);

  histo2d_["et_dr_final"] = fs->make<TH2D>("et_dr_final","Et vs dR(l,GSF)",200,0.,1000.,128,0.,6.4);
  histo2d_["et_dr_HEEP"] = fs->make<TH2D>("et_dr_HEEP","Et vs dR(l,GSF)",200,0.,1000.,128,0.,6.4);

  histo1d_["et_noGsf_final"] = fs->make<TH1D>("et_noGsf_final","Et",200,0.,1000.);
  histo1d_["et_noGsf_HEEP"] = fs->make<TH1D>("et_noGsf_HEEP","Et",200,0.,1000.);
}

void MergedEleCRanalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::GsfElectron>> eleHandle;
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

  std::string trigs[] = {
    // https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  };

  const unsigned int nTrigs = sizeof(trigs)/sizeof(*trigs);
  const unsigned int nTrig = trigResultHandle.product()->size();
  std::vector<std::pair<std::string, int>> indices;
  edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    std::string trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != nTrigs; jTrig++) {
      if (trigName.find(trigs[jTrig].substr(0, trigs[jTrig].find("*"))) != std::string::npos) {
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
  std::vector<edm::Ptr<reco::GsfElectron>> acceptEles;

  auto sortByPt = [](const edm::Ptr<reco::GsfElectron> a, const edm::Ptr<reco::GsfElectron> b) {
    return a->pt() > b->pt();
  };

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    auto aEle = eleHandle->ptrAt(idx);
    auto orgGsfTrk = aEle->gsfTrack();
    auto addGsfTrk = (*addGsfTrkMap)[aEle];

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
      if ( (*cutflow_HEEPHandle)[aEle] == static_cast<int>(MergedLeptonIDs::cutflowElectron::showerShape) )
        acceptEles.push_back(aEle);
    } else {
      if ( (*cutflow_modifiedHEEPHandle)[aEle] == static_cast<int>(MergedLeptonIDs::cutflowElectron::dxy) )
        acceptEles.push_back(aEle);
    }
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByPt);

  bool accepted = false;

  if ( acceptEles.size() > 1 ) {
    if ( acceptEles.at(1)->pt() > 35. )
      accepted = true;
  }

  if (!accepted)
    return;

  std::vector<edm::Ptr<reco::GsfElectron>> mergedEls;

  for (unsigned int idx = 0; idx < acceptEles.size(); ++idx) {
    auto aEle = acceptEles.at(idx);

    auto aGSFtype = static_cast<MergedLeptonIDs::GSFtype>((*GSFtype_mergedElectronHandle)[aEle]);
    auto aStatus = static_cast<MergedLeptonIDs::cutflowElectron>((*status_mergedElectronHandle)[aEle]);

    double et = MergedLeptonIDs::Et(aEle);
    double drGSF = std::sqrt(reco::deltaR2(aEle->gsfTrack()->eta(),
                                           aEle->gsfTrack()->phi(),
                                           (*addGsfTrkMap)[aEle]->eta(),
                                           (*addGsfTrkMap)[aEle]->phi()));

    std::string postfix = "";

    if ( aStatus==MergedLeptonIDs::cutflowElectron::no2ndGsf ||
         aStatus==MergedLeptonIDs::cutflowElectron::passedMVA2 ) {
      // no 2nd GSF
      histo1d_["et_noGsf_HEEP"]->Fill( et, aWeight );

      if ( std::abs(aEle->eta()) < 1.5 ) {
        postfix = "_bkgEt2EB";
        histo1d_[static_cast<std::string>("mva")+postfix]
          ->Fill( (*mva_mergedElectronNoGsfHandle)[aEle] , aWeight );
      }
    } else if ( aStatus==MergedLeptonIDs::cutflowElectron::has2ndGsf ||
                aStatus==MergedLeptonIDs::cutflowElectron::passedMVA1 ) {
      // has 2nd GSF
      histo2d_["et_dr_HEEP"]->Fill( et, drGSF, aWeight );

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
        case MergedLeptonIDs::GSFtype::DR1Et2EE:
          postfix = "_DR1Et2EE";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et1EE:
          postfix = "_DR2Et1EE";
          break;
        case MergedLeptonIDs::GSFtype::DR2Et2EE:
          postfix = "_DR2Et2EE";
          break;
        default:
          throw cms::Exception("LogicError") << "opening angle between the original and additional GSF track does not fit into any of MergedLeptonIDs::openingAngle categories" << std::endl;
          break;
      }

      histo1d_[static_cast<std::string>("mva")+postfix]
        ->Fill( (*mva_mergedElectronHandle)[aEle] , aWeight );
    } else {}

    histo1d_[static_cast<std::string>("status_mergedElectron")+postfix]
      ->Fill( static_cast<float>( (*status_mergedElectronHandle)[aEle] ) + 0.5 , aWeight );

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) {
      mergedEls.push_back(aEle);
      histo2d_["et_dr_final"]->Fill( et, drGSF, aWeight );
    } else if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ) {
      mergedEls.push_back(aEle);
      histo1d_["et_noGsf_final"]->Fill( et, aWeight );
    } else {}

  } // electron loop

  std::sort(mergedEls.begin(),mergedEls.end(),sortByPt);

  switch (acceptEles.size()) {
    case 4:
      histo1d_["cutflow_4E"]->Fill( 2.5, aWeight );
      if ( mergedEls.size()==0 )
        histo1d_["cutflow_4E"]->Fill( 3.5, aWeight );
      break;
    case 3:
      histo1d_["cutflow_3E"]->Fill( 2.5, aWeight );

      if ( mergedEls.size()==1 ) {
        histo1d_["cutflow_3E"]->Fill( 3.5, aWeight );

        if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) {
          histo1d_["pt_3E_ME"]->Fill( mergedEls.front()->pt(), aWeight );
          histo1d_["eta_3E_ME"]->Fill( mergedEls.front()->eta(), aWeight );
          histo1d_["phi_3E_ME"]->Fill( mergedEls.front()->phi(), aWeight );
        } else if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ) {
          histo1d_["pt_3E_ME_noGsf"]->Fill( mergedEls.front()->pt(), aWeight );
          histo1d_["eta_3E_ME_noGsf"]->Fill( mergedEls.front()->eta(), aWeight );
          histo1d_["phi_3E_ME_noGsf"]->Fill( mergedEls.front()->phi(), aWeight );
        } else {}
      }

      break;
    case 2:
      histo1d_["cutflow_2E"]->Fill( 2.5, aWeight );
      if ( mergedEls.size() > 0 )
        histo1d_["cutflow_2E"]->Fill( 3.5, aWeight );
      if ( mergedEls.size()==2 ) {
        histo1d_["cutflow_2E"]->Fill( 4.5, aWeight );

        if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) {
          histo1d_["pt_2E_ME1"]->Fill( mergedEls.front()->pt(), aWeight );
          histo1d_["eta_2E_ME1"]->Fill( mergedEls.front()->eta(), aWeight );
          histo1d_["phi_2E_ME1"]->Fill( mergedEls.front()->phi(), aWeight );
        } else if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ) {
          histo1d_["pt_2E_ME1_noGsf"]->Fill( mergedEls.front()->pt(), aWeight );
          histo1d_["eta_2E_ME1_noGsf"]->Fill( mergedEls.front()->eta(), aWeight );
          histo1d_["phi_2E_ME1_noGsf"]->Fill( mergedEls.front()->phi(), aWeight );
        } else {}

        if ( (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) {
          histo1d_["pt_2E_ME2"]->Fill( mergedEls.at(1)->pt(), aWeight );
          histo1d_["eta_2E_ME2"]->Fill( mergedEls.at(1)->eta(), aWeight );
          histo1d_["phi_2E_ME2"]->Fill( mergedEls.at(1)->phi(), aWeight );
        } else if ( (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ) {
          histo1d_["pt_2E_ME2_noGsf"]->Fill( mergedEls.at(1)->pt(), aWeight );
          histo1d_["eta_2E_ME2_noGsf"]->Fill( mergedEls.at(1)->eta(), aWeight );
          histo1d_["phi_2E_ME2_noGsf"]->Fill( mergedEls.at(1)->phi(), aWeight );
        } else {}

        const auto lvecME1 = mergedEls.front()->polarP4();
        const auto lvecME2 = mergedEls.at(1)->polarP4();
        const auto lvecll = lvecME1 + lvecME2;
        const double dr2ll = reco::deltaR2(lvecME1.eta(),lvecME1.phi(),lvecME2.eta(),lvecME2.phi());

        if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) &&
             (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) {
          histo1d_["invM_2E_ll"]->Fill( lvecll.M(), aWeight );
          histo1d_["rap_2E_ll"]->Fill( lvecll.Rapidity(), aWeight );
          histo1d_["pt_2E_ll"]->Fill( lvecll.pt(), aWeight );
          histo1d_["dr_2E_ll"]->Fill( std::sqrt(dr2ll), aWeight );
          histo1d_["deta_2E_ll"]->Fill( lvecME1.eta()-lvecME2.eta(), aWeight );
          histo1d_["dphi_2E_ll"]->Fill( reco::deltaPhi(lvecME1.phi(),lvecME2.phi()), aWeight );
          histo2d_["invM_dr"]->Fill( lvecll.M(), std::sqrt(dr2ll), aWeight );
        }

        if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) &&
             (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA2) ) {
          histo1d_["invM_2E_ll_noGsf"]->Fill( lvecll.M(), aWeight );
          histo1d_["rap_2E_ll_noGsf"]->Fill( lvecll.Rapidity(), aWeight );
          histo1d_["pt_2E_ll_noGsf"]->Fill( lvecll.pt(), aWeight );
          histo1d_["dr_2E_ll_noGsf"]->Fill( std::sqrt(dr2ll), aWeight );
          histo1d_["deta_2E_ll_noGsf"]->Fill( lvecME1.eta()-lvecME2.eta(), aWeight );
          histo1d_["dphi_2E_ll_noGsf"]->Fill( reco::deltaPhi(lvecME1.phi(),lvecME2.phi()), aWeight );
          histo2d_["invM_dr_noGsf"]->Fill( lvecll.M(), std::sqrt(dr2ll), aWeight );
        }

        if ( ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedAllMVA) &&
               (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) ||
             ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) &&
               (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedAllMVA) ) ) {
          histo1d_["invM_2E_ll_mixed"]->Fill( lvecll.M(), aWeight );
          histo1d_["rap_2E_ll_mixed"]->Fill( lvecll.Rapidity(), aWeight );
          histo1d_["pt_2E_ll_mixed"]->Fill( lvecll.pt(), aWeight );
          histo1d_["dr_2E_ll_mixed"]->Fill( std::sqrt(dr2ll), aWeight );
          histo1d_["deta_2E_ll_mixed"]->Fill( lvecME1.eta()-lvecME2.eta(), aWeight );
          histo1d_["dphi_2E_ll_mixed"]->Fill( reco::deltaPhi(lvecME1.phi(),lvecME2.phi()), aWeight );
          histo2d_["invM_dr_mixed"]->Fill( lvecll.M(), std::sqrt(dr2ll), aWeight );
        }
      }

      break;
  }
}

DEFINE_FWK_MODULE(MergedEleCRanalyzer);
