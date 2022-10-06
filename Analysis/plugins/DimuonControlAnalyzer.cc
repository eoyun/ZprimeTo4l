#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

class DimuonControlAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DimuonControlAnalyzer(const edm::ParameterSet&);
  virtual ~DimuonControlAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  bool isMC_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;

  const double ptThres_;

  std::map<std::string,TH1*> histo1d_;
};

DimuonControlAnalyzer::DimuonControlAnalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
ptThres_(iConfig.getParameter<double>("ptThres")) {
  usesResource("TFileService");
}

void DimuonControlAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  histo1d_["cutflow_2M"] = fs->make<TH1D>("cutflow_2M","cutflow",10,0.,10.);

  histo1d_["pt_M1"] = fs->make<TH1D>("pt_M1","Pt",500,0.,500.);
  histo1d_["eta_M1"] = fs->make<TH1D>("eta_M1","Eta",200,-2.5,2.5);
  histo1d_["phi_M1"] = fs->make<TH1D>("phi_M1","Phi",128,-3.2,3.2);

  histo1d_["pt_M2"] = fs->make<TH1D>("pt_M2","Pt",500,0.,500.);
  histo1d_["eta_M2"] = fs->make<TH1D>("eta_M2","Eta",200,-2.5,2.5);
  histo1d_["phi_M2"] = fs->make<TH1D>("phi_M2","Phi",128,-3.2,3.2);

  histo1d_["invM_ll"] = fs->make<TH1D>("invM_ll","M(ll)",500,0.,500.);
  histo1d_["invM_window_ll"] = fs->make<TH1D>("invM_window_ll","M(ll)",200,70,120.);
  histo1d_["rap_ll"] = fs->make<TH1D>("rap_ll","rapidity(ll)",200,-2.5,2.5);
  histo1d_["pt_ll"] = fs->make<TH1D>("pt_ll","Pt(ll)",500,0.,500.);
  histo1d_["dr_ll"] = fs->make<TH1D>("dr_ll","dR(ll)",128,0.,6.4);
}

void DimuonControlAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

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
  iEvent.getByToken(triggerToken_, trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  std::string trigs[] = {
    // https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
  };

  const unsigned int nTrigs = sizeof(trigs)/sizeof(*trigs);
  const unsigned int nTrig = trigResultHandle.product()->size();
  const edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    const std::string& trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != nTrigs; jTrig++) {
      if (trigName.find(trigs[jTrig].substr(0, trigs[jTrig].find("*"))) != std::string::npos) {
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
      for (unsigned int jTrig = 0; jTrig < nTrigs; jTrig++) {
        if ( name.find(trigs[jTrig].substr(0, trigs[jTrig].find("*"))) != std::string::npos &&
             trigObjInst.hasPathName(name,true,true) ) {
          trigObjs.push_back(trigObj);
        }
      } // wanted triggers
    } // fired triggers
  } // trigger objs

  histo1d_["cutflow_2M"]->Fill( 0.5, aWeight );

  if (!isFired)
    return;

  histo1d_["cutflow_2M"]->Fill( 1.5, aWeight );

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  std::vector<edm::RefToBase<pat::Muon>> highPtMuons;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->pt() < ptThres_ )
      continue;

    if ( aMuon->isHighPtMuon(primaryVertex) )
      highPtMuons.push_back(aMuon);
  }

  std::vector<edm::RefToBase<pat::Muon>> isolatedHighPtMuons;

  for (unsigned idx = 0; idx < highPtMuons.size(); idx++) {
    const auto& firstMuon = highPtMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    if ( firstIso/firstMuon->pt() < 0.05 )
      isolatedHighPtMuons.push_back(firstMuon);
  }

  auto sortByPt = [](const edm::RefToBase<pat::Muon> a, const edm::RefToBase<pat::Muon> b) {
    return a->pt() > b->pt();
  };

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByPt);

  if ( isolatedHighPtMuons.size() < 2 || isolatedHighPtMuons.front()->pt() < 52. )
    return;

  histo1d_["cutflow_2M"]->Fill( 2.5, aWeight );

  bool trigMatched = false;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),isolatedHighPtMuons.front()->eta(),isolatedHighPtMuons.front()->phi()) < 0.09 )
      trigMatched = true;
  }

  if ( !trigMatched )
    return;

  histo1d_["cutflow_2M"]->Fill( 3.5, aWeight );

  histo1d_["pt_M1"]->Fill( isolatedHighPtMuons.front()->pt(), aWeight );
  histo1d_["eta_M1"]->Fill( isolatedHighPtMuons.front()->eta(), aWeight );
  histo1d_["phi_M1"]->Fill( isolatedHighPtMuons.front()->phi(), aWeight );

  histo1d_["pt_M2"]->Fill( isolatedHighPtMuons.at(1)->pt(), aWeight );
  histo1d_["eta_M2"]->Fill( isolatedHighPtMuons.at(1)->eta(), aWeight );
  histo1d_["phi_M2"]->Fill( isolatedHighPtMuons.at(1)->phi(), aWeight );

  const auto& lvecM1 = isolatedHighPtMuons.front()->polarP4();
  const auto& lvecM2 = isolatedHighPtMuons.at(1)->polarP4();
  const auto lvecll = lvecM1 + lvecM2;
  const double dr2ll = reco::deltaR2(lvecM1.eta(),lvecM1.phi(),lvecM2.eta(),lvecM2.phi());

  histo1d_["invM_ll"]->Fill( lvecll.M(), aWeight );
  histo1d_["invM_window_ll"]->Fill( lvecll.M(), aWeight );
  histo1d_["rap_ll"]->Fill( lvecll.Rapidity(), aWeight );
  histo1d_["pt_ll"]->Fill( lvecll.pt(), aWeight );
  histo1d_["dr_ll"]->Fill( std::sqrt(dr2ll), aWeight );

  return;
}

DEFINE_FWK_MODULE(DimuonControlAnalyzer);
