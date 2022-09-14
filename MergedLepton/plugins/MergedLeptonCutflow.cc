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

class MergedLeptonCutflow : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonCutflow(const edm::ParameterSet&);
  virtual ~MergedLeptonCutflow() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  bool isMC_;

  edm::EDGetTokenT<edm::View<reco::GsfElectron>> srcEle_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::EDGetTokenT<double> prefweight_token;

  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

  edm::EDGetTokenT<edm::ValueMap<int>> cutflow_modifiedHEEPToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> cutflow_HEEPToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> status_mergedElectronToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> mva_mergedElectronToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> mva_mergedElectronNoGsfToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> openingAngle_mergedElectronToken_;
  edm::EDGetTokenT<int> nresolvedElectronToken_;
  edm::EDGetTokenT<int> nmergedElectronToken_;

  std::map<std::string,TH1*> histo1d_;
};

MergedLeptonCutflow::MergedLeptonCutflow(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
srcEle_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
cutflow_modifiedHEEPToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("cutflow_modifiedHEEP"))),
cutflow_HEEPToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("cutflow_HEEP"))),
status_mergedElectronToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("status_mergedElectron"))),
mva_mergedElectronToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mva_mergedElectron"))),
mva_mergedElectronNoGsfToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("mva_mergedElectronNoGsf"))),
openingAngle_mergedElectronToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("openingAngle_mergedElectron"))),
nresolvedElectronToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nresolvedElectron"))),
nmergedElectronToken_(consumes<int>(iConfig.getParameter<edm::InputTag>("nmergedElectron"))) {
  usesResource("TFileService");
}

void MergedLeptonCutflow::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow_modifiedHEEP_EB"] = fs->make<TH1D>("cutflow_modifiedHEEP_EB","cutflow",11,-1.,10.);
  histo1d_["cutflow_modifiedHEEP_EE"] = fs->make<TH1D>("cutflow_modifiedHEEP_EE","cutflow",11,-1.,10.);
  histo1d_["cutflow_HEEP_EB"] = fs->make<TH1D>("cutflow_HEEP_EB","cutflow",16,-1.,15.);
  histo1d_["cutflow_HEEP_EE"] = fs->make<TH1D>("cutflow_HEEP_EE","cutflow",16,-1.,15.);
  histo1d_["status_mergedElectron_DR1EB"] = fs->make<TH1D>("status_mergedElectron_DR1EB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR2EB"] = fs->make<TH1D>("status_mergedElectron_DR2EB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR3EB"] = fs->make<TH1D>("status_mergedElectron_DR3EB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR1EE"] = fs->make<TH1D>("status_mergedElectron_DR1EE","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR2EE"] = fs->make<TH1D>("status_mergedElectron_DR2EE","status",30,-1.,29.);
  histo1d_["status_mergedElectron_DR3EE"] = fs->make<TH1D>("status_mergedElectron_DR3EE","status",30,-1.,29.);
  histo1d_["status_mergedElectron_bkgEB"] = fs->make<TH1D>("status_mergedElectron_bkgEB","status",30,-1.,29.);
  histo1d_["status_mergedElectron_bkgEE"] = fs->make<TH1D>("status_mergedElectron_bkgEE","status",30,-1.,29.);
  histo1d_["mva1_DR1EB"] = fs->make<TH1D>("mva1_DR1EB","MVA score",200,-1.,1.);
  histo1d_["mva1_DR2EB"] = fs->make<TH1D>("mva1_DR2EB","MVA score",200,-1.,1.);
  histo1d_["mva1_DR3EB"] = fs->make<TH1D>("mva1_DR3EB","MVA score",200,-1.,1.);
  histo1d_["mva1_DR1EE"] = fs->make<TH1D>("mva1_DR1EE","MVA score",200,-1.,1.);
  histo1d_["mva1_DR2EE"] = fs->make<TH1D>("mva1_DR2EE","MVA score",200,-1.,1.);
  histo1d_["mva1_DR3EE"] = fs->make<TH1D>("mva1_DR3EE","MVA score",200,-1.,1.);
  histo1d_["mva1_bkgEB"] = fs->make<TH1D>("mva1_bkgEB","MVA score",200,-1.,1.);
  histo1d_["mva1_bkgEE"] = fs->make<TH1D>("mva1_bkgEE","MVA score",200,-1.,1.);
  histo1d_["mva2_bkgEB"] = fs->make<TH1D>("mva2_bkgEB","MVA score",200,-1.,1.);
  histo1d_["mva2_bkgEE"] = fs->make<TH1D>("mva2_bkgEE","MVA score",200,-1.,1.);
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

  histo1d_["pt_3E_ME_noGsf"] = fs->make<TH1D>("pt_3E_ME_noGsf","Pt",200,0.,500.);
  histo1d_["eta_3E_ME_noGsf"] = fs->make<TH1D>("eta_3E_ME_noGsf","3Eta",200,-2.5,2.5);
  histo1d_["phi_3E_ME_noGsf"] = fs->make<TH1D>("phi_3E_ME_noGsf","Phi",128,-3.2,3.2);

  histo1d_["pt_2E_ME1_noGsf"] = fs->make<TH1D>("pt_2E_ME1_noGsf","Pt",200,0.,500.);
  histo1d_["eta_2E_ME1_noGsf"] = fs->make<TH1D>("eta_2E_ME1_noGsf","Eta",200,-2.5,2.5);
  histo1d_["phi_2E_ME1_noGsf"] = fs->make<TH1D>("phi_2E_ME1_noGsf","Phi",128,-3.2,3.2);

  histo1d_["pt_2E_ME2_noGsf"] = fs->make<TH1D>("pt_2E_ME2_noGsf","Pt",200,0.,500.);
  histo1d_["eta_2E_ME2_noGsf"] = fs->make<TH1D>("eta_2E_ME2_noGsf","Eta",200,-2.5,2.5);
  histo1d_["phi_2E_ME2_noGsf"] = fs->make<TH1D>("phi_2E_ME2_noGsf","Phi",128,-3.2,3.2);
}

void MergedLeptonCutflow::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  edm::Handle<edm::ValueMap<int>> cutflow_modifiedHEEPHandle;
  iEvent.getByToken(cutflow_modifiedHEEPToken_, cutflow_modifiedHEEPHandle);

  edm::Handle<edm::ValueMap<int>> cutflow_HEEPHandle;
  iEvent.getByToken(cutflow_HEEPToken_, cutflow_HEEPHandle);

  edm::Handle<edm::ValueMap<int>> status_mergedElectronHandle;
  iEvent.getByToken(status_mergedElectronToken_, status_mergedElectronHandle);

  edm::Handle<edm::ValueMap<float>> mva_mergedElectronHandle;
  iEvent.getByToken(mva_mergedElectronToken_, mva_mergedElectronHandle);

  edm::Handle<edm::ValueMap<float>> mva_mergedElectronNoGsfHandle;
  iEvent.getByToken(mva_mergedElectronNoGsfToken_, mva_mergedElectronNoGsfHandle);

  edm::Handle<edm::ValueMap<int>> openingAngle_mergedElectronHandle;
  iEvent.getByToken(openingAngle_mergedElectronToken_, openingAngle_mergedElectronHandle);

  edm::Handle<int> nresolvedElectronHandle;
  iEvent.getByToken(nresolvedElectronToken_, nresolvedElectronHandle);

  edm::Handle<int> nmergedElectronHandle;
  iEvent.getByToken(nmergedElectronToken_, nmergedElectronHandle);

  // first two (modified)HEEP electrons should be Et > 35 GeV
  std::vector<edm::Ptr<reco::GsfElectron>> acceptEles;

  auto sortByPt = [](const edm::Ptr<reco::GsfElectron> a, const edm::Ptr<reco::GsfElectron> b) {
    return a->pt() > b->pt();
  };

  for (unsigned int idx = 0; idx < eleHandle->size(); idx++) {
    auto aEle = eleHandle->ptrAt(idx);

    if ( (*cutflow_modifiedHEEPHandle)[aEle] == static_cast<int>(MergedLeptonIDs::cutflowElectron::dxy) )
      acceptEles.push_back(aEle);
  }

  std::sort(acceptEles.begin(),acceptEles.end(),sortByPt);

  bool accepted = false;

  if ( acceptEles.size() > 1 ) {
    if ( acceptEles.at(1)->pt() > 35. )
      accepted = true;
  }

  if (!accepted)
    return;

  int npassed_modHEEP = 0;
  std::vector<edm::Ptr<reco::GsfElectron>> mergedEls;

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    auto aEle = eleHandle->ptrAt(idx);

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

    if ( (*cutflow_modifiedHEEPHandle)[aEle] == static_cast<int>(MergedLeptonIDs::cutflowElectron::dxy) )
      npassed_modHEEP++;

    auto aopenangle = static_cast<MergedLeptonIDs::openingAngle>((*openingAngle_mergedElectronHandle)[aEle]);

    std::string postfix = "";

    if (aopenangle==MergedLeptonIDs::openingAngle::nullAngle) {
      // no 2nd GSF
      if ( std::abs(aEle->eta()) < 1.5 )
        postfix = "_bkgEB";
      else
        postfix = "_bkgEE";

      histo1d_[static_cast<std::string>("mva2")+postfix]
        ->Fill( (*mva_mergedElectronNoGsfHandle)[aEle] , aWeight );

    } else {
      // has 2nd GSF
      switch (aopenangle) { // better way to do this?
        case MergedLeptonIDs::openingAngle::DR1EB:
          postfix = "_DR1EB";
          break;
        case MergedLeptonIDs::openingAngle::DR2EB:
          postfix = "_DR2EB";
          break;
        case MergedLeptonIDs::openingAngle::DR3EB:
          postfix = "_DR3EB";
          break;
        case MergedLeptonIDs::openingAngle::DR1EE:
          postfix = "_DR1EE";
          break;
        case MergedLeptonIDs::openingAngle::DR2EE:
          postfix = "_DR2EE";
          break;
        case MergedLeptonIDs::openingAngle::DR3EE:
          postfix = "_DR3EE";
          break;
        default:
          throw cms::Exception("LogicError") << "opening angle between the original and additional GSF track does not fit into any of MergedLeptonIDs::openingAngle categories" << std::endl;
          break;
      }
    }

    histo1d_[static_cast<std::string>("status_mergedElectron")+postfix]
      ->Fill( static_cast<float>( (*status_mergedElectronHandle)[aEle] ) + 0.5 , aWeight );

    histo1d_[static_cast<std::string>("mva1")+postfix]
      ->Fill( (*mva_mergedElectronHandle)[aEle] , aWeight );

    if ( (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ||
         (*status_mergedElectronHandle)[aEle]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedAllMVA) )
      mergedEls.push_back(aEle);

  } // electron loop

  std::sort(mergedEls.begin(),mergedEls.end(),sortByPt);

  switch (npassed_modHEEP) {
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
        } else if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedAllMVA) ) {
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
        } else if ( (*status_mergedElectronHandle)[mergedEls.front()]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedAllMVA) ) {
          histo1d_["pt_2E_ME1_noGsf"]->Fill( mergedEls.front()->pt(), aWeight );
          histo1d_["eta_2E_ME1_noGsf"]->Fill( mergedEls.front()->eta(), aWeight );
          histo1d_["phi_2E_ME1_noGsf"]->Fill( mergedEls.front()->phi(), aWeight );
        } else {}

        if ( (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedMVA1) ) {
          histo1d_["pt_2E_ME2"]->Fill( mergedEls.at(1)->pt(), aWeight );
          histo1d_["eta_2E_ME2"]->Fill( mergedEls.at(1)->eta(), aWeight );
          histo1d_["phi_2E_ME2"]->Fill( mergedEls.at(1)->phi(), aWeight );
        } else if ( (*status_mergedElectronHandle)[mergedEls.at(1)]==static_cast<int>(MergedLeptonIDs::cutflowElectron::passedAllMVA) ) {
          histo1d_["pt_2E_ME2_noGsf"]->Fill( mergedEls.at(1)->pt(), aWeight );
          histo1d_["eta_2E_ME2_noGsf"]->Fill( mergedEls.at(1)->eta(), aWeight );
          histo1d_["phi_2E_ME2_noGsf"]->Fill( mergedEls.at(1)->phi(), aWeight );
        } else {}
      }

      break;
  }
}

DEFINE_FWK_MODULE(MergedLeptonCutflow);
