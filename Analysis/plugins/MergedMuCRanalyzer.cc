#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"

class MergedMuonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedMuonAnalyzer(const edm::ParameterSet&);
  ~MergedMuonAnalyzer() override {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  bool isMC_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;

  const edm::EDGetTokenT<edm::View<reco::Muon>> muonToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  const double ptThres_;
  const double drThres_;
  const double ratioBarrel_;
  const double ratioEndcap_;

  std::map<std::string,TH1*> histo1d_;
};

MergedMuonAnalyzer::MergedMuonAnalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getUntrackedParameter<bool>("isMC")),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
muonToken_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")),
ratioBarrel_(iConfig.getParameter<double>("ratioBarrel")),
ratioEndcap_(iConfig.getParameter<double>("ratioEndcap")) {
  usesResource("TFileService");
}

void MergedMuonAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  histo1d_["cutflow_4M"] = fs->make<TH1D>("cutflow_4M","cutflow",10,0.,10.);
  histo1d_["cutflow_3M"] = fs->make<TH1D>("cutflow_3M","cutflow",10,0.,10.);

  histo1d_["pt_3M_MM"] = fs->make<TH1D>("pt_3M_MM","Pt",200,0.,500.);
  histo1d_["eta_3M_MM"] = fs->make<TH1D>("eta_3M_MM","3Eta",200,-2.5,2.5);
  histo1d_["phi_3M_MM"] = fs->make<TH1D>("phi_3M_MM","Phi",128,-3.2,3.2);

  histo1d_["pt_3M_M1"] = fs->make<TH1D>("pt_3M_M1","Pt",200,0.,500.);
  histo1d_["eta_3M_M1"] = fs->make<TH1D>("eta_3M_M1","3Eta",200,-2.5,2.5);
  histo1d_["phi_3M_M1"] = fs->make<TH1D>("phi_3M_M1","Phi",128,-3.2,3.2);

  histo1d_["pt_3M_M2"] = fs->make<TH1D>("pt_3M_M2","Pt",200,0.,500.);
  histo1d_["eta_3M_M2"] = fs->make<TH1D>("eta_3M_M2","3Eta",200,-2.5,2.5);
  histo1d_["phi_3M_M2"] = fs->make<TH1D>("phi_3M_M2","Phi",128,-3.2,3.2);

  histo1d_["pt_3M_MET"] = fs->make<TH1D>("pt_3M_MET","Pt",200,0.,500.);
  histo1d_["phi_3M_MET"] = fs->make<TH1D>("phi_3M_MET","Phi",128,-3.2,3.2);

  histo1d_["dphi_3M_MET"] = fs->make<TH1D>("dphi_3M_MET","dPhi",128,-3.2,3.2);
}

void MergedMuonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Muon>> muonHandle;
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
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  std::string trigs[] = {
    // https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
    // https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    // https://twiki.cern.ch/twiki/bin/view/CMS/MuonHLT2016
    // "HLT_DoubleEle33_CaloIdL_MW_v*",
    // "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
    "HLT_Mu50_v*",
    "HLT_TkMu50_v*"
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

  histo1d_["cutflow_4M"]->Fill( 0.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 0.5, aWeight );

  if (!isFired)
    return;

  histo1d_["cutflow_4M"]->Fill( 1.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 1.5, aWeight );

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  std::vector<edm::Ptr<reco::Muon>> highPtMuons;
  std::vector<edm::Ptr<reco::Muon>> highPtTrackerMuons; // but not highPtMuon
  std::vector<edm::Ptr<reco::Muon>> looseMuons; // but neither of above

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    auto aMuon = muonHandle->ptrAt(idx);

    if ( aMuon->pt() < 25. )
      continue;

    if ( muon::isHighPtMuon(*aMuon,primaryVertex) )
      highPtMuons.push_back(aMuon);
    // else if ( MergedLeptonIDs::isHighPtTrackerMuon(*aMuon,primaryVertex) )
    else if ( muon::isTrackerHighPtMuon(*aMuon,primaryVertex) )
      highPtTrackerMuons.push_back(aMuon);
    else if ( muon::isLooseMuon(*aMuon) )
      looseMuons.push_back(aMuon);
    else {}
  }

  auto sortByPt = [](const edm::Ptr<reco::Muon> a, const edm::Ptr<reco::Muon> b) {
    return a->pt() > b->pt();
  };

  std::sort(highPtMuons.begin(),highPtMuons.end(),sortByPt);
  std::sort(highPtTrackerMuons.begin(),highPtTrackerMuons.end(),sortByPt);

  if ( highPtMuons.size() < 2 )
    return;

  if ( highPtMuons.front()->pt() < 50. )
    return;

  histo1d_["cutflow_4M"]->Fill( 2.5, aWeight );
  histo1d_["cutflow_3M"]->Fill( 2.5, aWeight );

  // concatenate(-ish) highPtMuons & highPtTrackerMuons - order is important!
  std::vector<edm::Ptr<reco::Muon>> allHighPtMuons(highPtMuons);
  allHighPtMuons.insert( allHighPtMuons.end(), highPtTrackerMuons.begin(), highPtTrackerMuons.end() );

  unsigned int nHighPtMuons = allHighPtMuons.size();

  if ( nHighPtMuons==4 ) {
    histo1d_["cutflow_4M"]->Fill( 3.5, aWeight );

    if ( looseMuons.size()==0 )
      histo1d_["cutflow_4M"]->Fill( 4.5, aWeight );
  }

  if ( nHighPtMuons==3 ) {
    histo1d_["cutflow_3M"]->Fill( 3.5, aWeight );

    if ( looseMuons.size()==0 ) {
      histo1d_["cutflow_3M"]->Fill( 4.5, aWeight );

      double drThres2 = drThres_*drThres_;

      edm::Ptr<reco::Muon> firstMuon, secondMuon;

      // find a collimated pair of muons
      for (unsigned idx = 0; idx < highPtMuons.size(); idx++) { // one of them should be global
        const auto& cand1 = allHighPtMuons.at(idx);

        bool found = false;

        for (unsigned jdx = idx+1; jdx < allHighPtMuons.size(); jdx++) {
          const auto& cand2 = allHighPtMuons.at(jdx);

          // they shouldn't be the same, throw exception otherwise
          if ( cand1==cand2 )
            throw cms::Exception("LogicError") << "Error: MergedMuonAnalyzer::analyze - attempting to compare the same muons " << std::endl;

          double dr2 = reco::deltaR2(cand1->eta(),cand1->phi(),cand2->eta(),cand2->phi());

          if ( dr2 < drThres2 ) {
            firstMuon = cand1;
            secondMuon = cand2;
            found = true;

            break;
          }
        }

        if (found)
          break;
      }

      // third muon is the merged muon - should be highPtMuon
      edm::Ptr<reco::Muon> mergedMuon;

      for (unsigned idx = 0; idx < highPtMuons.size(); idx++) {
        const auto& cand = highPtMuons.at(idx);

        if ( cand!=firstMuon && cand!=secondMuon ) {
          mergedMuon = cand;
          break;
        }
      }

      if ( !firstMuon.isNull() && !secondMuon.isNull() && !mergedMuon.isNull() ) {
        histo1d_["cutflow_3M"]->Fill( 5.5, aWeight );

        const auto& aMET = metHandle->at(0);

        if ( aMET.pt() > ptThres_ ) {
          histo1d_["cutflow_3M"]->Fill( 6.5, aWeight );

          double dphi = reco::deltaPhi(mergedMuon->phi(),aMET.phi());

          if ( dphi < drThres_ ) {
            histo1d_["cutflow_3M"]->Fill( 7.5, aWeight );

            double ptsum1 = firstMuon->pt()+secondMuon->pt();
            double ptsum2 = mergedMuon->pt()+aMET.pt();

            if ( std::abs(mergedMuon->eta()) < 1.2 ) {
              if ( ratioBarrel_*ptsum1 > ptsum2 )
                return;
            } else {
              if ( ratioEndcap_*ptsum1 > ptsum2 )
                return;
            }

            histo1d_["cutflow_3M"]->Fill( 8.5, aWeight );

            histo1d_["pt_3M_MM"]->Fill( mergedMuon->pt(), aWeight );
            histo1d_["eta_3M_MM"]->Fill( mergedMuon->eta(), aWeight );
            histo1d_["phi_3M_MM"]->Fill( mergedMuon->phi(), aWeight );

            histo1d_["pt_3M_MET"]->Fill( aMET.pt(), aWeight );
            histo1d_["phi_3M_MET"]->Fill( aMET.phi(), aWeight );

            histo1d_["pt_3M_M1"]->Fill( firstMuon->pt(), aWeight );
            histo1d_["eta_3M_M1"]->Fill( firstMuon->eta(), aWeight );
            histo1d_["phi_3M_M1"]->Fill( firstMuon->phi(), aWeight );

            histo1d_["pt_3M_M2"]->Fill( secondMuon->pt(), aWeight );
            histo1d_["eta_3M_M2"]->Fill( secondMuon->eta(), aWeight );
            histo1d_["phi_3M_M2"]->Fill( secondMuon->phi(), aWeight );

            histo1d_["dphi_3M_MET"]->Fill( dphi, aWeight );
          } // MET dPhi
        } // MET ptThres
      } // find a pair of collimated muons
    } // loose muon veto
  } // nHighPtMuons==3

  return;
}

DEFINE_FWK_MODULE(MergedMuonAnalyzer);
