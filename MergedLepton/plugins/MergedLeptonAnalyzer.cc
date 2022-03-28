#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

class MergedLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonAnalyzer(const edm::ParameterSet&);
  virtual ~MergedLeptonAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  static bool sortByTuneP(const edm::Ptr<reco::Muon> a, const edm::Ptr<reco::Muon> b);
  void initTObj(edm::Service<TFileService>&, const std::string&);

  template <typename ObjStruct>
  void initTree(std::map<std::string,ObjStruct>& values,
                edm::Service<TFileService>& fs,
                const std::string& name,
                const TString& typestr,
                const std::string& prefix,
                const std::string& postfix);

  int categorize(std::vector<edm::Ptr<reco::Muon>>&, const std::vector<edm::Ptr<reco::GenParticle>>&);
  bool isHighPtTrackerMuon(const reco::Muon& muon, const reco::Vertex& vtx);

  void fillByCategory(std::vector<edm::Ptr<reco::Muon>>& muons,
                      std::vector<edm::Ptr<reco::Muon>>& candidates,
                      const std::vector<edm::Ptr<reco::GenParticle>>& promptLeptons,
                      std::vector<edm::Ptr<reco::GenParticle>>::const_iterator muonItr,
                      const edm::Handle<edm::View<reco::Muon>>& muonHandle,
                      const reco::Vertex& vtx);

  void fillByGsfTrack(std::vector<edm::Ptr<reco::GsfElectron>>& eles,
                      std::vector<edm::Ptr<reco::GenParticle>>::const_iterator eleItr,
                      const edm::Handle<edm::ValueMap<float>>& trkIsoMapHandle,
                      const edm::Handle<edm::ValueMap<float>>& ecalIsoMapHandle,
                      const edm::Handle<edm::ValueMap<int>>& nrSatCrysHandle,
                      const edm::Handle<std::vector<reco::Conversion>>& conversions,
                      const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                      const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkHandle,
                      const reco::Vertex& vtx);

  void fillMuons(const edm::Ptr<reco::Muon>&, const edm::Ptr<reco::GenParticle>&, const edm::Ptr<reco::GenParticle>&, const reco::Vertex&, std::string);
  void fillMets(const edm::Ptr<reco::Muon>&, const edm::Ptr<reco::GenParticle>&, const pat::MET&, const std::vector<edm::Ptr<reco::Muon>>&, std::string);
  void fillMets(const edm::Ptr<reco::Candidate>&, const pat::MET&, const std::vector<edm::Ptr<reco::Candidate>>&, std::string);

  void fillElectrons(const edm::Ptr<reco::GsfElectron>& el,
                     const reco::Vertex& vtx,
                     const float& trkIso,
                     const float& ecalIso,
                     const int& nrSatCrys,
                     const edm::Handle<reco::ConversionCollection>& conversions,
                     const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                     const double& genPt,
                     const std::string& prefix);

  void fillGsfTracks(const reco::GsfTrackRef& addGsfTrk,
                     const reco::GsfTrackRef& orgGsfTrk,
                     const reco::Vertex& vtx,
                     const std::string& prefix);

  bool isModifiedHEEP(const reco::GsfElectron& el,
                      const reco::Vertex& primaryVertex,
                      const float& trkIso,
                      const float& ecalIso,
                      const int& nrSatCrys,
                      const double& rho,
                      int& cutflow);
  bool hasPassedHEEP(const reco::GsfElectron& el,
                     const reco::Vertex& primaryVertex,
                     const int& nrSatCrys,
                     const double& rho,
                     int& cutflow);

  edm::EDGetTokenT<edm::View<reco::Muon>> srcMuon_;
  edm::EDGetTokenT<edm::View<reco::GsfElectron>> srcEle_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<edm::View<pat::MET>> pfMETToken_;
  edm::EDGetTokenT<edm::View<pat::MET>> puppiMETToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> nrSatCrysMapToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  edm::ESHandle<TransientTrackBuilder> TTBuilder_;
  edm::ESHandle<MagneticField> magField_;
  edm::ParameterSet vtxFitterPset_;

  const double ptThres_;
  const double drThres_;

  typedef struct {
    int numberOfValidTrackerHits, numberOfValidPixelHits, numberOfValidStripHits,
    trackerLayersWithMeasurement, pixelLayersWithMeasurement, stripLayersWithMeasurement,
    trackerLayersWithoutMeasurement, pixelLayersWithoutMeasurement, stripLayersWithoutMeasurement;
    int isHighPt, isTrackerHighPt, isPFloose, isPFmedium, isPFtight;
    int isGlobal, isTracker, isStdAlone, bestType, tunePtype, nShower, nMatchedStations, nExpectedStations;
    float trackerVoM, pixelVoM, stripVoM, pfPt, tunepPt, genPt, genEta, genPhi, PFoGen, TPoGen, trackerPt, TRKoGen;
    float pairPt, pairEta, pairPhi;
    float globalChi2, globalNormChi2, trackerChi2, trackerNormChi2, tunepChi2, tunepNormChi2;
  } MuonStruct;

  typedef struct {
    float PFphi, PFdPhi, PFpt, PFoGen, PFoMu, MuEta, PFSumEt, TPphi, TPdPhi, TPpt, TPoGen, TPoMu, TPSumEt, dSumEt, sumEtRatio, sumEtRatioTP;
    int nLepton;
  } METstruct;

  typedef struct {
    float pt, eta, phi, en, genPt;
    int charge;
    float enSC, etSC, etaSC, phiSC, etaSCWidth, phiSCWidth;
    float full5x5_sigmaIetaIeta, full5x5_sigmaIphiIphi;
    float full5x5_E1x5, full5x5_E2x5, full5x5_E5x5, full5x5_hOverE, full5x5_r9;
    float dEtaIn, dPhiIn, dPhiSeed, dEtaEle, dPhiEle, dEtaSeed;
    float EseedOverP, EOverP;
    float ecalEn, ecalErr, trkErr, combErr, PFcombErr;
    float dr03EcalRecHitSumEt, dr03HcalDepth1TowerSumEt;
    int nrSatCrys;
    float modTrkIso, modEcalIso;
    int lostHits, nValidHits, nValidPixelHits, GsfHits;
    float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz;
    float Gsfpt, Gsfeta, Gsfphi;
    int expectedMissingInnerHits;
    float convVtxFitProb, convVtxChi2, convDist, convDcot, convRadius;
    int passConversionVeto, nbrem;
    float fbrem, fbremSC;
  } ElectronStruct;

  typedef struct {
    float Gsfpt, Gsfeta, Gsfphi;
    int lostHits, nValidHits, nValidPixelHits;
    float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz;
    int vtxValid;
    float vtx_dx, vtx_dy, vtx_dz, vtx_chi2, vtx_xErr, vtx_yErr, vtx_zErr;
    float vtx_pt, vtx_rapidity, vtx_phi, vtx_M;
  } AddGsfStruct;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TTree*> tree_;
  std::map<std::string,MuonStruct> values_;
  std::map<std::string,METstruct> metvalues_;
  std::map<std::string,ElectronStruct> elvalues_;
  std::map<std::string,AddGsfStruct> gsfvalues_;

  TString mustr_ = TString("numberOfValidTrackerHits/I:numberOfValidPixelHits:numberOfValidStripHits:")
  + "trackerLayersWithMeasurement:pixelLayersWithMeasurement:stripLayersWithMeasurement:"
  + "trackerLayersWithoutMeasurement:pixelLayersWithoutMeasurement:stripLayersWithoutMeasurement:"
  + "isHighPt:isTrackerHighPt:isPFloose:isPFmedium:isPFtight:"
  + "isGlobal:isTracker:isStdAlone:bestType:tunePtype:nShower:nMatchedStations:nExpectedStations:"
  + "trackerVoM/F:pixelVoM:stripVoM:pfPt:tunepPt:genPt:genEta:genPhi:PFoGen:TPoGen:trackerPt:TRKoGen:"
  + "pairPt:pairEta:pairPhi:"
  + "globalChi2:globalNormChi2:trackerChi2:trackerNormChi2:tunepChi2:tunepNormChi2";

  TString metstr_ = TString("PFphi/F:PFdPhi:PFpt:PFoGen:PFoMu:MuEta:PFSumEt:")
  + "TPphi:TPdPhi:TPpt:TPoGen:TPoMu:TPSumEt:dSumEt:sumEtRatio:sumEtRatioTP:nLepton/I";

  TString elstr_ = TString("pt/F:eta:phi:en:genPt:charge/I:")
  + "enSC/F:etSC:etaSC:phiSC:etaSCWidth:phiSCWidth:"
  + "full5x5_sigmaIetaIeta:full5x5_sigmaIphiIphi:"
  + "full5x5_E1x5:full5x5_E2x5:full5x5_E5x5:full5x5_hOverE:full5x5_r9:"
  + "dEtaIn:dPhiIn:dPhiSeed:dEtaEle:dPhiEle:dEtaSeed:"
  + "EseedOverP:EOverP:"
  + "ecalEn:ecalErr:trkErr:combErr:PFcombErr:"
  + "dr03EcalRecHitSumEt:dr03HcalDepth1TowerSumEt:"
  + "nrSatCrys/I:"
  + "modTrkIso/F:modEcalIso:"
  + "lostHits/I:nValidHits:nValidPixelHits:GsfHits:"
  + "chi2/F:d0:d0Err:dxyErr:vz:dzErr:dxy:dz:"
  + "Gsfpt/F:Gsfeta:Gsfphi:"
  + "expectedMissingInnerHits/I:"
  + "convVtxFitProb/F:convVtxChi2:convDist:convDcot:convRadius:"
  + "passConversionVeto/I:nbrem:"
  + "fbrem/F:fbremSC";

  TString addgsfstr_ = TString("Gsfpt/F:Gsfeta:Gsfphi:")
  + "lostHits/I:nValidHits:nValidPixelHits:"
  + "chi2/F:d0:d0Err:dxyErr:vz:dzErr:dxy:dz:"
  + "vtxValid/I:"
  + "vtx_dx/F:vtx_dy:vtx_dz:vtx_chi2:vtx_xErr:vtx_yErr:vtx_zErr:"
  + "vtx_pt:vtx_rapidity:vtx_phi:vtx_M";
};

MergedLeptonAnalyzer::MergedLeptonAnalyzer(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcEle_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pfMETToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcPFMET"))),
puppiMETToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcPuppiMET"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
nrSatCrysMapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("nrSatCrysMap"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
vtxFitterPset_(iConfig.getParameter<edm::ParameterSet>("KFParameters")),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")) {
  usesResource("TFileService");
}

void MergedLeptonAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",10,0.,10.);
  histo1d_["cutflow_HEEP"] = fs->make<TH1D>("cutflow_HEEP","cutflow_HEEP",15,0.,15.);
  histo1d_["pt_H"] = fs->make<TH1D>("pt_H","p_{T}(H)",200,0.,200.);
  histo1d_["isEcalDriven"] = fs->make<TH1D>("isEcalDriven","isEcalDriven",2,0.,2.);
  histo1d_["isTrackerDriven"] = fs->make<TH1D>("isTrackerDriven","isTrackerDriven",2,0.,2.);

  initTObj(fs,"solo");
  initTObj(fs,"tag");
  initTObj(fs,"probe");
  initTObj(fs,"highPt1");
  initTObj(fs,"highPt2");

  initTree<METstruct>(metvalues_,fs,"METstruct",metstr_,"pf","MET");
  initTree<METstruct>(metvalues_,fs,"METstruct",metstr_,"puppi","MET");
  initTree<METstruct>(metvalues_,fs,"METstruct",metstr_,"stdPf","MET");
  initTree<METstruct>(metvalues_,fs,"METstruct",metstr_,"stdPuppi","MET");

  initTree<ElectronStruct>(elvalues_,fs,"ElectronStruct",elstr_,"heep1","el");
  initTree<ElectronStruct>(elvalues_,fs,"ElectronStruct",elstr_,"heep2","el");
  initTree<ElectronStruct>(elvalues_,fs,"ElectronStruct",elstr_,"mergedEl1","el");
  initTree<ElectronStruct>(elvalues_,fs,"ElectronStruct",elstr_,"mergedEl2","el");

  initTree<AddGsfStruct>(gsfvalues_,fs,"AddGsfStruct",addgsfstr_,"heep1Gsf","addGsf");
  initTree<AddGsfStruct>(gsfvalues_,fs,"AddGsfStruct",addgsfstr_,"heep2Gsf","addGsf");
  initTree<AddGsfStruct>(gsfvalues_,fs,"AddGsfStruct",addgsfstr_,"mergedEl1Gsf","addGsf");
}

void MergedLeptonAnalyzer::initTObj(edm::Service<TFileService>& fs, const std::string& prefix) {
  histo1d_[prefix+"_muType"] = fs->make<TH1D>(TString(prefix)+"_muType","muon type",3,0.,3.);
  histo1d_[prefix+"_algo"] = fs->make<TH1D>(TString(prefix)+"_algo","muon algo",46,0.,46.);
  histo1d_[prefix+"_orgAlgo"] = fs->make<TH1D>(TString(prefix)+"_orgAlgo","muon original algo",46,0.,46.);
  histo1d_[prefix+"_algoMask"] = fs->make<TH1D>(TString(prefix)+"_algoMask","muon algo mask",46,0.,46.);
  histo1d_[prefix+"_quality"] = fs->make<TH1D>(TString(prefix)+"_quality","muon quality",8,0.,8.);
  histo1d_[prefix+"_passIDs"] = fs->make<TH1D>(TString(prefix)+"_passIDs","muon passed IDs",5,0.,5.);

  initTree<MuonStruct>(values_,fs,"MuonStruct",mustr_,prefix,"muon");
}

template <typename ObjStruct>
void MergedLeptonAnalyzer::initTree(std::map<std::string,ObjStruct>& values,
                                    edm::Service<TFileService>& fs,
                                    const std::string& name,
                                    const TString& typestr,
                                    const std::string& prefix,
                                    const std::string& postfix) {
  values[prefix+"_"+postfix] = ObjStruct();
  tree_[prefix+"_"+postfix+"Tree"] = fs->make<TTree>(TString(prefix)+"_"+postfix+"Tree",TString(postfix)+"Tree");
  tree_[prefix+"_"+postfix+"Tree"]->Branch(TString(name),&(values[prefix+"_"+postfix]),typestr);
}

void MergedLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<reco::GsfElectron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::View<pat::MET>> pfMEThandle;
  iEvent.getByToken(pfMETToken_, pfMEThandle);

  edm::Handle<edm::View<pat::MET>> puppiMEThandle;
  iEvent.getByToken(puppiMETToken_, puppiMEThandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> ecalIsoMapHandle;
  iEvent.getByToken(ecalIsoToken_, ecalIsoMapHandle);

  edm::Handle<edm::ValueMap<int>> nrSatCrysHandle;
  iEvent.getByToken(nrSatCrysMapToken_, nrSatCrysHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_,rhoHandle);

  edm::Handle<reco::ConversionCollection> conversionsHandle;
  iEvent.getByToken(conversionsToken_, conversionsHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder_);
  iSetup.get<IdealMagneticFieldRecord>().get(magField_);

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  std::vector<edm::Ptr<reco::GenParticle>> promptLeptons;
  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    auto genptc = genptcHandle->ptrAt(idx);

    if ( genptc->isHardProcess() && genptc->pdgId()==35 )
      histo1d_["pt_H"]->Fill(genptc->pt());

    if ( ( std::abs(genptc->pdgId())==11 || std::abs(genptc->pdgId())==13 ) && genptc->isPromptFinalState() && genptc->pt() > ptThres_ )
      promptLeptons.push_back(genptc);
  }

  if (promptLeptons.size()!=4)
    return;

  auto sisterLambda = [](const edm::Ptr<reco::GenParticle> a, const edm::Ptr<reco::GenParticle> b) {
    return (a->mother() == b->mother()) ? (a->pt() > b->pt()) : (a->mother() < b->mother());
  };

  std::sort(promptLeptons.begin(),promptLeptons.end(),sisterLambda);

  std::vector<edm::Ptr<reco::Muon>> highptMuons1, highptMuons2;
  std::vector<edm::Ptr<reco::GsfElectron>> heeps1, heeps2;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    auto aMuon = muonHandle->ptrAt(idx);

    if ( !(aMuon->pt() > ptThres_) || !muon::isHighPtMuon(*aMuon,primaryVertex) )
      continue;

    bool matched = false;
    size_t igen = 0;

    for (; igen < promptLeptons.size(); ++igen) {
      auto genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if (matched)
      (igen < 2) ? highptMuons1.push_back(aMuon) : highptMuons2.push_back(aMuon);
  }

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    auto aEle = eleHandle->ptrAt(idx);
    int cutflow = 0;

    if ( !(aEle->pt() > ptThres_) )
      continue;

    bool matched = false;
    size_t igen = 0;

    for (; igen < promptLeptons.size(); ++igen) {
      auto genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aEle->eta(),aEle->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if ( !matched )
      continue;

    histo1d_["isEcalDriven"]->Fill(static_cast<float>(aEle->ecalDrivenSeed())+0.5);
    histo1d_["isTrackerDriven"]->Fill(static_cast<float>(aEle->trackerDrivenSeed())+0.5);

    bool passModifiedHEEP = isModifiedHEEP(*aEle,
                                           primaryVertex,
                                           (*trkIsoMapHandle)[aEle],
                                           (*ecalIsoMapHandle)[aEle],
                                           (*nrSatCrysHandle)[aEle],
                                           *rhoHandle,
                                           cutflow);
    histo1d_["cutflow"]->Fill(static_cast<float>(cutflow)+0.5);

    int cutflow_HEEP = 0;
    hasPassedHEEP(*aEle,
                  primaryVertex,
                  (*nrSatCrysHandle)[aEle],
                  *rhoHandle,
                  cutflow_HEEP);
    histo1d_["cutflow_HEEP"]->Fill(static_cast<float>(cutflow_HEEP)+0.5);

    if ( !passModifiedHEEP )
      return; // WARNING veto events

    if (matched)
      (igen < 2) ? heeps1.push_back(aEle) : heeps2.push_back(aEle);
  }

  bool noEle = true;
  bool noMu = true;

  // assume 4e or 4mu only
  if ( heeps1.size() > 0 && heeps2.size() > 0 )
    noEle = false;

  if ( highptMuons1.size() > 0 && highptMuons2.size() > 0 )
    noMu = false;

  if ( noEle && noMu )
    return;

  std::vector<edm::Ptr<reco::Muon>> candidates1, candidates2;
  std::vector<edm::Ptr<reco::Muon>> finalMuons;

  auto pushMuon = [](std::vector<edm::Ptr<reco::Muon>>& fnal,
                     const std::vector<edm::Ptr<reco::Muon>>& highpt,
                     const std::vector<edm::Ptr<reco::Muon>>& cands) {
    int ncand = fnal.size();
    fnal.push_back(highpt.front());

    if ( highpt.size() > 1 )
      fnal.push_back(highpt.at(1));
    else {
      if ( cands.size() > 0 )
        fnal.push_back(cands.front());
    }

    return static_cast<int>(fnal.size()) - ncand;
  };

  int nmuon = 0, nCand1 = 0, nCand2 = 0;

  if ( !noMu ) {
    fillByCategory(highptMuons1,candidates1,promptLeptons,promptLeptons.begin(),muonHandle,primaryVertex);
    fillByCategory(highptMuons2,candidates2,promptLeptons,std::next(promptLeptons.begin(),2),muonHandle,primaryVertex);
    nCand1 = pushMuon(finalMuons,highptMuons1, candidates1);
    nCand2 = pushMuon(finalMuons,highptMuons2, candidates2);
    nmuon = nCand1 + nCand2;
  }

  if ( !noEle ) {
    fillByGsfTrack(heeps1,
                   promptLeptons.begin(),
                   trkIsoMapHandle,
                   ecalIsoMapHandle,
                   nrSatCrysHandle,
                   conversionsHandle,
                   beamSpotHandle,
                   addGsfTrkHandle,
                   primaryVertex);
    fillByGsfTrack(heeps2,
                   std::next(promptLeptons.begin(),2),
                   trkIsoMapHandle,
                   ecalIsoMapHandle,
                   nrSatCrysHandle,
                   conversionsHandle,
                   beamSpotHandle,
                   addGsfTrkHandle,
                   primaryVertex);
  }

  if ( nmuon < 4 && !noMu ) {
    edm::Ptr<reco::Muon> mergedMuon;

    if ( nmuon < 3 )
      mergedMuon = ( highptMuons1.front()->tunePMuonBestTrack()->pt() < highptMuons2.front()->tunePMuonBestTrack()->pt() )
                   ? highptMuons1.front() : highptMuons2.front();
    else
      mergedMuon = ( nCand1 < 2 ) ? highptMuons1.front() : highptMuons2.front();

    float tunepPt = mergedMuon->tunePMuonBestTrack()->pt();
    auto matched = ( mergedMuon==highptMuons1.front() ) ? promptLeptons.begin() : std::next(promptLeptons.begin(),2);
    auto probe = std::next(matched,1);

    if ( std::abs(tunepPt-(*matched)->pt()) > std::abs(tunepPt-(*probe)->pt()) ) {
      probe = matched;
      matched = std::next(matched,1);
    }

    fillMets(mergedMuon, *probe, pfMEThandle->at(0), finalMuons, "pf");
    fillMets(mergedMuon, *probe, puppiMEThandle->at(0), finalMuons, "puppi");
  } else if ( nmuon==4 || noMu ) {
    std::vector<edm::Ptr<reco::Candidate>> finalCands;
    auto fillCandidate = [&finalCands](const edm::Ptr<reco::Candidate>& cand) { finalCands.push_back(cand); };
    std::for_each(finalMuons.begin(),finalMuons.end(),fillCandidate);
    std::for_each(heeps1.begin(),heeps1.end(),fillCandidate);
    std::for_each(heeps2.begin(),heeps2.end(),fillCandidate);

    if ( finalCands.empty() )
      return;

    auto findNN = [&finalCands](const pat::MET& met) {
      edm::Ptr<reco::Candidate> candidate;
      double dphicand = DBL_MAX;

      for (auto aLepton : finalCands) {
        double dphi = reco::deltaPhi(aLepton->phi(),met.phi());

        if (dphi < dphicand)
          candidate = aLepton;
      }

      return candidate;
    };

    fillMets(findNN(pfMEThandle->at(0)),pfMEThandle->at(0), finalCands, "stdPf");
    fillMets(findNN(puppiMEThandle->at(0)),puppiMEThandle->at(0), finalCands, "stdPuppi");
  }
}

int MergedLeptonAnalyzer::categorize(std::vector<edm::Ptr<reco::Muon>>& muons, const std::vector<edm::Ptr<reco::GenParticle>>& promptLeptons) {
  if ( muons.size() > 1 )
    std::sort(muons.begin(),muons.end(),sortByTuneP);
  else if ( muons.size() == 1 ) {
    int nMatched = 0;
    const auto highptmu = muons.front();

    for (const auto genptc : promptLeptons) {
      double dr2 = reco::deltaR2(highptmu->eta(),highptmu->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres_*drThres_)
        nMatched++;
    }

    if ( nMatched < 2 )
      return -1; // not merged
  }

  return static_cast<int>(muons.size());
}

void MergedLeptonAnalyzer::fillByCategory(std::vector<edm::Ptr<reco::Muon>>& muons,
  std::vector<edm::Ptr<reco::Muon>>& candidates,
  const std::vector<edm::Ptr<reco::GenParticle>>& promptLeptons,
  std::vector<edm::Ptr<reco::GenParticle>>::const_iterator muonItr,
  const edm::Handle<edm::View<reco::Muon>>& muonHandle,
  const reco::Vertex& vtx) {
  switch (categorize(muons,promptLeptons)) {
    case -1:
      return;
    case 0:
      return;
    case 1: {
      auto aMuon = muons.front();
      float tunepPt = aMuon->tunePMuonBestTrack()->pt();

      auto matched = muonItr;
      auto probe = std::next(muonItr,1);

      if ( std::abs(tunepPt-(*matched)->pt()) > std::abs(tunepPt-(*probe)->pt()) ) {
        probe = matched;
        matched = std::next(matched,1);
      }

      for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
        auto otherMuon = muonHandle->ptrAt(idx);

        if ( aMuon==otherMuon )
          continue;

        if ( otherMuon->pt() > ptThres_ ) {
          double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),otherMuon->eta(),otherMuon->phi());

          if ( dr2 < drThres_*drThres_ )
            candidates.push_back(otherMuon);
        }
      }

      if (candidates.size()==0) {
        fillMuons(aMuon,*matched,*probe,vtx,"solo");
        return;
      }

      std::sort(candidates.begin(),candidates.end(),sortByTuneP);
      fillMuons(aMuon,*matched,*probe,vtx,"tag");
      fillMuons(candidates.front(),*probe,*matched,vtx,"probe");

      return;
    } default:
      fillMuons(muons.front(),*muonItr,*std::next(muonItr,1),vtx,"highPt1");
      fillMuons(muons.at(1),*std::next(muonItr,1),*muonItr,vtx,"highPt2");
      break;
  }

  return;
}

void MergedLeptonAnalyzer::fillByGsfTrack(std::vector<edm::Ptr<reco::GsfElectron>>& eles,
                                          std::vector<edm::Ptr<reco::GenParticle>>::const_iterator eleItr,
                                          const edm::Handle<edm::ValueMap<float>>& trkIsoMapHandle,
                                          const edm::Handle<edm::ValueMap<float>>& ecalIsoMapHandle,
                                          const edm::Handle<edm::ValueMap<int>>& nrSatCrysHandle,
                                          const edm::Handle<std::vector<reco::Conversion>>& conversions,
                                          const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                                          const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkHandle,
                                          const reco::Vertex& vtx) {
  if ( eles.size()==1 ) {
    auto addGsfTrk = (*addGsfTrkHandle)[eles.front()];
    auto orgGsfTrk = eles.front()->gsfTrack();

    if ( addGsfTrk->pt()==orgGsfTrk->pt() && addGsfTrk->eta()==orgGsfTrk->eta() && addGsfTrk->phi()==orgGsfTrk->phi() ) {
      fillElectrons(eles.front(),
                    vtx,
                    (*trkIsoMapHandle)[eles.front()],
                    (*ecalIsoMapHandle)[eles.front()],
                    (*nrSatCrysHandle)[eles.front()],
                    conversions,
                    beamSpotHandle,
                    (*eleItr)->pt()+(*std::next(eleItr,1))->pt(),
                    "mergedEl2");
    } else {
      fillElectrons(eles.front(),
                    vtx,
                    (*trkIsoMapHandle)[eles.front()],
                    (*ecalIsoMapHandle)[eles.front()],
                    (*nrSatCrysHandle)[eles.front()],
                    conversions,
                    beamSpotHandle,
                    (*eleItr)->pt()+(*std::next(eleItr,1))->pt(),
                    "mergedEl1");
      fillGsfTracks((*addGsfTrkHandle)[eles.front()],eles.front()->gsfTrack(),vtx,"mergedEl1Gsf");
    }

    return;
  }

  // default case: eles.size() > 1
  std::sort(eles.begin(),eles.end(),[](const edm::Ptr<reco::GsfElectron>& a, const edm::Ptr<reco::GsfElectron>& b) { return a->pt() > b->pt(); });
  fillElectrons(eles.at(0),
                vtx,
                (*trkIsoMapHandle)[eles.at(0)],
                (*ecalIsoMapHandle)[eles.at(0)],
                (*nrSatCrysHandle)[eles.at(0)],
                conversions,
                beamSpotHandle,
                (*eleItr)->pt(),
                "heep1");
  fillElectrons(eles.at(1),
                vtx,
                (*trkIsoMapHandle)[eles.at(1)],
                (*ecalIsoMapHandle)[eles.at(1)],
                (*nrSatCrysHandle)[eles.at(1)],
                conversions,
                beamSpotHandle,
                (*std::next(eleItr,1))->pt(),
                "heep2");
  fillGsfTracks((*addGsfTrkHandle)[eles.at(0)],eles.at(0)->gsfTrack(),vtx,"heep1Gsf");
  fillGsfTracks((*addGsfTrkHandle)[eles.at(1)],eles.at(1)->gsfTrack(),vtx,"heep2Gsf");

  return;
}

bool MergedLeptonAnalyzer::sortByTuneP(const edm::Ptr<reco::Muon> a, const edm::Ptr<reco::Muon> b) {
  return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
}

void MergedLeptonAnalyzer::fillMuons(const edm::Ptr<reco::Muon>& aMuon,
                                     const edm::Ptr<reco::GenParticle>& matched,
                                     const edm::Ptr<reco::GenParticle>& probe,
                                     const reco::Vertex& vtx,
                                     std::string prefix) {
  if (aMuon->isGlobalMuon())
    histo1d_[prefix+"_muType"]->Fill(0.5);
  if (aMuon->isTrackerMuon())
    histo1d_[prefix+"_muType"]->Fill(1.5);
  if (aMuon->isStandAloneMuon())
    histo1d_[prefix+"_muType"]->Fill(2.5);

  if (muon::isHighPtMuon(*aMuon,vtx))
    histo1d_[prefix+"_passIDs"]->Fill(0.5);
  if (isHighPtTrackerMuon(*aMuon,vtx))
    histo1d_[prefix+"_passIDs"]->Fill(1.5);
  if (muon::isLooseMuon(*aMuon))
    histo1d_[prefix+"_passIDs"]->Fill(2.5);
  if (muon::isMediumMuon(*aMuon))
    histo1d_[prefix+"_passIDs"]->Fill(3.5);
  if (muon::isTightMuon(*aMuon,vtx))
    histo1d_[prefix+"_passIDs"]->Fill(4.5);

  values_[prefix+"_muon"].isHighPt = muon::isHighPtMuon(*aMuon,vtx);
  values_[prefix+"_muon"].isTrackerHighPt = isHighPtTrackerMuon(*aMuon,vtx);
  values_[prefix+"_muon"].isPFloose = muon::isLooseMuon(*aMuon);
  values_[prefix+"_muon"].isPFmedium = muon::isMediumMuon(*aMuon);
  values_[prefix+"_muon"].isPFtight = muon::isTightMuon(*aMuon,vtx);

  values_[prefix+"_muon"].isGlobal = aMuon->isGlobalMuon();
  values_[prefix+"_muon"].isTracker = aMuon->isTrackerMuon();
  values_[prefix+"_muon"].isStdAlone = aMuon->isStandAloneMuon();
  values_[prefix+"_muon"].bestType = static_cast<int>(aMuon->muonBestTrackType());
  values_[prefix+"_muon"].tunePtype = static_cast<int>(aMuon->tunePMuonBestTrackType());
  values_[prefix+"_muon"].nShower = aMuon->numberOfShowers();
  values_[prefix+"_muon"].nMatchedStations = aMuon->numberOfMatchedStations();
  values_[prefix+"_muon"].nExpectedStations = aMuon->expectedNnumberOfMatchedStations();

  values_[prefix+"_muon"].pfPt = aMuon->pt();
  values_[prefix+"_muon"].tunepPt = aMuon->tunePMuonBestTrack()->pt();
  values_[prefix+"_muon"].genPt = matched->pt();
  values_[prefix+"_muon"].genEta = matched->eta();
  values_[prefix+"_muon"].genPhi = matched->phi();
  values_[prefix+"_muon"].PFoGen = aMuon->pt()/matched->pt();
  values_[prefix+"_muon"].TPoGen = aMuon->tunePMuonBestTrack()->pt()/matched->pt();
  values_[prefix+"_muon"].trackerPt = aMuon->isTrackerMuon() ? aMuon->innerTrack()->pt() : -1.;
  values_[prefix+"_muon"].TRKoGen = aMuon->isTrackerMuon() ? aMuon->innerTrack()->pt()/matched->pt() : -1.;
  values_[prefix+"_muon"].pairPt = probe->pt();
  values_[prefix+"_muon"].pairEta = probe->eta();
  values_[prefix+"_muon"].pairPhi = probe->phi();
  values_[prefix+"_muon"].globalChi2 = aMuon->isGlobalMuon() ? aMuon->globalTrack()->chi2() : -1.;
  values_[prefix+"_muon"].globalNormChi2 = aMuon->isGlobalMuon() ? aMuon->globalTrack()->normalizedChi2() : -1.;
  values_[prefix+"_muon"].trackerChi2 = aMuon->isTrackerMuon() ? aMuon->innerTrack()->chi2() : -1.;
  values_[prefix+"_muon"].trackerNormChi2 = aMuon->isTrackerMuon() ? aMuon->innerTrack()->normalizedChi2() : -1.;
  values_[prefix+"_muon"].tunepChi2 = aMuon->isGlobalMuon() ? aMuon->tunePMuonBestTrack()->chi2() : -1.;
  values_[prefix+"_muon"].tunepNormChi2 = aMuon->isGlobalMuon() ? aMuon->tunePMuonBestTrack()->normalizedChi2() : -1.;

  if (!aMuon->isTrackerMuon()) {
    values_[prefix+"_muon"].numberOfValidTrackerHits = -1;
    values_[prefix+"_muon"].numberOfValidPixelHits = -1;
    values_[prefix+"_muon"].numberOfValidStripHits = -1;
    values_[prefix+"_muon"].trackerLayersWithMeasurement = -1;
    values_[prefix+"_muon"].pixelLayersWithMeasurement = -1;
    values_[prefix+"_muon"].stripLayersWithMeasurement = -1;
    values_[prefix+"_muon"].trackerLayersWithoutMeasurement = -1;
    values_[prefix+"_muon"].pixelLayersWithoutMeasurement = -1;
    values_[prefix+"_muon"].stripLayersWithoutMeasurement = -1;
    values_[prefix+"_muon"].trackerVoM = -1.;
    values_[prefix+"_muon"].pixelVoM = -1.;
    values_[prefix+"_muon"].stripVoM = -1.;
    tree_[prefix+"_muonTree"]->Fill();

    return;
  }

  auto innerTrack = aMuon->innerTrack();
  values_[prefix+"_muon"].numberOfValidTrackerHits = innerTrack->hitPattern().numberOfValidTrackerHits();
  values_[prefix+"_muon"].numberOfValidPixelHits = innerTrack->hitPattern().numberOfValidPixelHits();
  values_[prefix+"_muon"].numberOfValidStripHits = innerTrack->hitPattern().numberOfValidStripHits();
  values_[prefix+"_muon"].trackerLayersWithMeasurement = innerTrack->hitPattern().trackerLayersWithMeasurement();
  values_[prefix+"_muon"].pixelLayersWithMeasurement = innerTrack->hitPattern().pixelLayersWithMeasurement();
  values_[prefix+"_muon"].stripLayersWithMeasurement = innerTrack->hitPattern().stripLayersWithMeasurement();
  values_[prefix+"_muon"].trackerLayersWithoutMeasurement = innerTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  values_[prefix+"_muon"].pixelLayersWithoutMeasurement = innerTrack->hitPattern().pixelLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  values_[prefix+"_muon"].stripLayersWithoutMeasurement = innerTrack->hitPattern().stripLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  values_[prefix+"_muon"].trackerVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidTrackerHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
  values_[prefix+"_muon"].pixelVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidPixelHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());
  values_[prefix+"_muon"].stripVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidStripHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().stripLayersWithMeasurement());
  tree_[prefix+"_muonTree"]->Fill();

  histo1d_[prefix+"_algo"]->Fill(static_cast<float>(innerTrack->algo())+0.5);
  histo1d_[prefix+"_orgAlgo"]->Fill(static_cast<float>(innerTrack->originalAlgo())+0.5);

  for (int en = reco::TrackBase::undefAlgorithm; en != reco::TrackBase::algoSize; ++en) {
    if ( (innerTrack->algoMask())[en] )
      histo1d_[prefix+"_algoMask"]->Fill(static_cast<float>(en)+0.5);
  }

  for (int en = reco::TrackBase::loose; en != reco::TrackBase::qualitySize; ++en) {
    if ( innerTrack->quality( static_cast<reco::TrackBase::TrackQuality>(en) ) )
      histo1d_[prefix+"_quality"]->Fill(static_cast<float>(en)+0.5);
  }
}

void MergedLeptonAnalyzer::fillElectrons(const edm::Ptr<reco::GsfElectron>& el,
                                         const reco::Vertex& vtx,
                                         const float& trkIso,
                                         const float& ecalIso,
                                         const int& nrSatCrys,
                                         const edm::Handle<reco::ConversionCollection>& conversions,
                                         const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                                         const double& genPt,
                                         const std::string& prefix) {
  auto sq = [](const double& val) { return val*val; };
  double R = std::sqrt(sq(el->superCluster()->x()) + sq(el->superCluster()->y()) + sq(el->superCluster()->z()));
  double Rt = std::sqrt(sq(el->superCluster()->x()) + sq(el->superCluster()->y()));

  elvalues_[prefix+"_el"].pt = el->pt();
  elvalues_[prefix+"_el"].eta = el->eta();
  elvalues_[prefix+"_el"].phi = el->phi();
  elvalues_[prefix+"_el"].en = el->energy();
  elvalues_[prefix+"_el"].genPt = genPt;
  elvalues_[prefix+"_el"].charge = el->charge();

  elvalues_[prefix+"_el"].enSC = el->superCluster()->energy();
  elvalues_[prefix+"_el"].etSC = (el->superCluster()->energy())*(Rt/R);
  elvalues_[prefix+"_el"].etaSC = el->superCluster()->eta();
  elvalues_[prefix+"_el"].phiSC = el->superCluster()->phi();
  elvalues_[prefix+"_el"].etaSCWidth = el->superCluster()->etaWidth();
  elvalues_[prefix+"_el"].phiSCWidth = el->superCluster()->phiWidth();

  elvalues_[prefix+"_el"].full5x5_sigmaIetaIeta = el->full5x5_sigmaIetaIeta();
  elvalues_[prefix+"_el"].full5x5_sigmaIphiIphi = el->full5x5_sigmaIphiIphi();
  elvalues_[prefix+"_el"].full5x5_E1x5 = el->full5x5_e1x5();
  elvalues_[prefix+"_el"].full5x5_E2x5 = el->full5x5_e2x5Max();
  elvalues_[prefix+"_el"].full5x5_E5x5 = el->full5x5_e5x5();
  elvalues_[prefix+"_el"].full5x5_hOverE = el->full5x5_hcalOverEcal();
  elvalues_[prefix+"_el"].full5x5_r9 = el->full5x5_r9();

  elvalues_[prefix+"_el"].dEtaIn = el->deltaEtaSuperClusterTrackAtVtx();
  elvalues_[prefix+"_el"].dPhiIn = el->deltaPhiSuperClusterTrackAtVtx();
  elvalues_[prefix+"_el"].dPhiSeed = el->deltaPhiSeedClusterTrackAtCalo();
  elvalues_[prefix+"_el"].dEtaEle = el->deltaEtaEleClusterTrackAtCalo();
  elvalues_[prefix+"_el"].dPhiEle = el->deltaPhiEleClusterTrackAtCalo();
  elvalues_[prefix+"_el"].dEtaSeed = el->deltaEtaSeedClusterTrackAtVtx();
  elvalues_[prefix+"_el"].EseedOverP = el->eSeedClusterOverP();
  elvalues_[prefix+"_el"].EOverP = el->eSuperClusterOverP();

  elvalues_[prefix+"_el"].ecalEn = el->ecalEnergy();
  elvalues_[prefix+"_el"].ecalErr = el->ecalEnergyError();
  elvalues_[prefix+"_el"].trkErr = el->trackMomentumError();
  elvalues_[prefix+"_el"].combErr = el->p4Error(reco::GsfElectron::P4_COMBINATION);
  elvalues_[prefix+"_el"].PFcombErr = el->p4Error(reco::GsfElectron::P4_PFLOW_COMBINATION);

  elvalues_[prefix+"_el"].dr03EcalRecHitSumEt = el->dr03EcalRecHitSumEt();
  elvalues_[prefix+"_el"].dr03HcalDepth1TowerSumEt = el->dr03HcalDepth1TowerSumEt();
  elvalues_[prefix+"_el"].nrSatCrys = nrSatCrys;
  elvalues_[prefix+"_el"].modTrkIso = trkIso;
  elvalues_[prefix+"_el"].modEcalIso = ecalIso;

  reco::GsfTrackRef TrackRef = el->gsfTrack();
  elvalues_[prefix+"_el"].lostHits = TrackRef->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  elvalues_[prefix+"_el"].nValidHits = TrackRef->hitPattern().numberOfValidHits();
  elvalues_[prefix+"_el"].nValidPixelHits = TrackRef->hitPattern().numberOfValidPixelHits();
  elvalues_[prefix+"_el"].chi2 = TrackRef->normalizedChi2();
  elvalues_[prefix+"_el"].GsfHits = TrackRef->hitPattern().trackerLayersWithMeasurement();

  elvalues_[prefix+"_el"].d0 = TrackRef->d0();
  elvalues_[prefix+"_el"].d0Err = TrackRef->d0Error();
  elvalues_[prefix+"_el"].dxyErr = TrackRef->dxyError();
  elvalues_[prefix+"_el"].vz = TrackRef->vz();
  elvalues_[prefix+"_el"].dzErr = TrackRef->dzError();

  elvalues_[prefix+"_el"].dxy = TrackRef->dxy(vtx.position());
  elvalues_[prefix+"_el"].dz = TrackRef->dz(vtx.position());

  elvalues_[prefix+"_el"].Gsfpt = TrackRef->pt();
  elvalues_[prefix+"_el"].Gsfeta = TrackRef->eta();
  elvalues_[prefix+"_el"].Gsfphi = TrackRef->phi();

  elvalues_[prefix+"_el"].expectedMissingInnerHits = el->gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS); // 94X

  elvalues_[prefix+"_el"].convVtxFitProb = -1.;
  elvalues_[prefix+"_el"].convVtxChi2 = -1.;
  elvalues_[prefix+"_el"].passConversionVeto = !ConversionTools::hasMatchedConversion(*el, *conversions, beamSpotHandle->position());
  const reco::Conversion* convRef = ConversionTools::matchedConversion(*el, *conversions, beamSpotHandle->position());

  if ( convRef!=nullptr ) {
    reco::Vertex convVtx = convRef->conversionVertex();
    if (convVtx.isValid()) {
      elvalues_[prefix+"_el"].convVtxFitProb = TMath::Prob(convVtx.chi2(), convVtx.ndof());
      elvalues_[prefix+"_el"].convVtxChi2 = convVtx.normalizedChi2();
    }
  }

  elvalues_[prefix+"_el"].convDist = el->convDist();
  elvalues_[prefix+"_el"].convDcot = el->convDcot();
  elvalues_[prefix+"_el"].convRadius = el->convRadius();

  elvalues_[prefix+"_el"].fbrem = el->fbrem();
  elvalues_[prefix+"_el"].nbrem = el->numberOfBrems();
  elvalues_[prefix+"_el"].fbremSC = el->superClusterFbrem();

  tree_[prefix+"_elTree"]->Fill();
}

void MergedLeptonAnalyzer::fillGsfTracks(const reco::GsfTrackRef& addGsfTrk,
                                         const reco::GsfTrackRef& orgGsfTrk,
                                         const reco::Vertex& vtx,
                                         const std::string& prefix) {
  gsfvalues_[prefix+"_addGsf"].Gsfpt = addGsfTrk->pt();
  gsfvalues_[prefix+"_addGsf"].Gsfeta = addGsfTrk->eta();
  gsfvalues_[prefix+"_addGsf"].Gsfphi = addGsfTrk->phi();
  gsfvalues_[prefix+"_addGsf"].lostHits = addGsfTrk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS);
  gsfvalues_[prefix+"_addGsf"].nValidHits = addGsfTrk->hitPattern().numberOfValidHits();
  gsfvalues_[prefix+"_addGsf"].nValidPixelHits = addGsfTrk->hitPattern().numberOfValidPixelHits();
  gsfvalues_[prefix+"_addGsf"].chi2 = addGsfTrk->normalizedChi2();
  gsfvalues_[prefix+"_addGsf"].d0 = addGsfTrk->d0();
  gsfvalues_[prefix+"_addGsf"].d0Err = addGsfTrk->d0Error();
  gsfvalues_[prefix+"_addGsf"].dxyErr = addGsfTrk->dxyError();
  gsfvalues_[prefix+"_addGsf"].vz = addGsfTrk->vz();
  gsfvalues_[prefix+"_addGsf"].dzErr = addGsfTrk->dzError();

  gsfvalues_[prefix+"_addGsf"].dxy = addGsfTrk->dxy(vtx.position());
  gsfvalues_[prefix+"_addGsf"].dz = addGsfTrk->dz(vtx.position());

  KalmanVertexFitter fitter(vtxFitterPset_);
  TransientVertex theVertex;
  std::vector<reco::TransientTrack> Gsf12TT;

  if ( orgGsfTrk->pt()!=addGsfTrk->pt() || orgGsfTrk->eta()!=addGsfTrk->eta() ) {
    Gsf12TT.push_back(TTBuilder_->build(orgGsfTrk));
    Gsf12TT.push_back(TTBuilder_->build(addGsfTrk));
    theVertex = fitter.vertex(Gsf12TT);
  }

  gsfvalues_[prefix+"_addGsf"].vtxValid = theVertex.isValid();
  gsfvalues_[prefix+"_addGsf"].vtx_dx = 0.;
  gsfvalues_[prefix+"_addGsf"].vtx_dy = 0.;
  gsfvalues_[prefix+"_addGsf"].vtx_dz = 0.;

  gsfvalues_[prefix+"_addGsf"].vtx_chi2 = -1.;
  gsfvalues_[prefix+"_addGsf"].vtx_xErr = -1.;
  gsfvalues_[prefix+"_addGsf"].vtx_yErr = -1.;
  gsfvalues_[prefix+"_addGsf"].vtx_zErr = -1.;

  gsfvalues_[prefix+"_addGsf"].vtx_pt = -1.;
  gsfvalues_[prefix+"_addGsf"].vtx_rapidity = 0.;
  gsfvalues_[prefix+"_addGsf"].vtx_phi = 0.;
  gsfvalues_[prefix+"_addGsf"].vtx_M = -1.;

  if (theVertex.isValid()) {
    reco::Vertex theRecoVtx = theVertex;

    gsfvalues_[prefix+"_addGsf"].vtx_dx = theRecoVtx.x() - vtx.position().x();
    gsfvalues_[prefix+"_addGsf"].vtx_dy = theRecoVtx.y() - vtx.position().y();
    gsfvalues_[prefix+"_addGsf"].vtx_dz = theRecoVtx.z() - vtx.position().z();

    gsfvalues_[prefix+"_addGsf"].vtx_chi2 = theRecoVtx.normalizedChi2();
    gsfvalues_[prefix+"_addGsf"].vtx_xErr = theRecoVtx.xError();
    gsfvalues_[prefix+"_addGsf"].vtx_yErr = theRecoVtx.yError();
    gsfvalues_[prefix+"_addGsf"].vtx_zErr = theRecoVtx.zError();

    auto momentumSum = theRecoVtx.p4(0.0005109990615);
    gsfvalues_[prefix+"_addGsf"].vtx_pt = momentumSum.Pt();
    gsfvalues_[prefix+"_addGsf"].vtx_rapidity = momentumSum.Rapidity();
    gsfvalues_[prefix+"_addGsf"].vtx_phi = momentumSum.Phi();
    gsfvalues_[prefix+"_addGsf"].vtx_M = momentumSum.M();
  }

  tree_[prefix+"_addGsfTree"]->Fill();
}

void MergedLeptonAnalyzer::fillMets(const edm::Ptr<reco::Muon>& aMuon,
                                    const edm::Ptr<reco::GenParticle>& matched,
                                    const pat::MET& met,
                                    const std::vector<edm::Ptr<reco::Muon>>& finalMuons,
                                    std::string prefix) {
  metvalues_[prefix+"_MET"].PFphi = met.phi();
  metvalues_[prefix+"_MET"].PFdPhi = reco::deltaPhi(met.phi(),matched->phi());
  metvalues_[prefix+"_MET"].PFpt = met.pt();
  metvalues_[prefix+"_MET"].PFoGen = met.pt()/matched->pt();
  metvalues_[prefix+"_MET"].PFoMu = met.pt()/aMuon->pt();
  metvalues_[prefix+"_MET"].MuEta = aMuon->eta();
  metvalues_[prefix+"_MET"].PFSumEt = met.sumEt();

  pat::MET::Vector2 tpvec;
  tpvec.px = met.px() + aMuon->px() - aMuon->tunePMuonBestTrack()->px();
  tpvec.py = met.py() + aMuon->py() - aMuon->tunePMuonBestTrack()->py();
  metvalues_[prefix+"_MET"].TPphi = tpvec.phi();
  metvalues_[prefix+"_MET"].TPdPhi = reco::deltaPhi(tpvec.phi(),matched->phi());
  metvalues_[prefix+"_MET"].TPpt = tpvec.pt();
  metvalues_[prefix+"_MET"].TPoGen = tpvec.pt()/matched->pt();
  metvalues_[prefix+"_MET"].TPoMu = tpvec.pt()/aMuon->tunePMuonBestTrack()->pt();

  float sumpt = 0.;
  float sumTPpt = 0.;

  for (auto mu : finalMuons) {
    sumpt += mu->pt();
    sumTPpt += mu->tunePMuonBestTrack()->pt();
  }

  metvalues_[prefix+"_MET"].TPSumEt = met.sumEt() + sumTPpt - sumpt;
  metvalues_[prefix+"_MET"].dSumEt = met.sumEt() - sumpt;
  metvalues_[prefix+"_MET"].sumEtRatio = met.sumEt()/sumpt;
  metvalues_[prefix+"_MET"].sumEtRatioTP = metvalues_[prefix+"_MET"].TPSumEt/sumTPpt;
  metvalues_[prefix+"_MET"].nLepton = static_cast<int>(finalMuons.size());

  tree_[prefix+"_METTree"]->Fill();
}

void MergedLeptonAnalyzer::fillMets(const edm::Ptr<reco::Candidate>& aLepton,
                                    const pat::MET& met,
                                    const std::vector<edm::Ptr<reco::Candidate>>& finalLeptons,
                                    std::string prefix) {
  metvalues_[prefix+"_MET"].PFphi = met.phi();
  metvalues_[prefix+"_MET"].PFdPhi = reco::deltaPhi(met.phi(),aLepton->phi());
  metvalues_[prefix+"_MET"].PFpt = met.pt();
  metvalues_[prefix+"_MET"].PFoGen = 1.;
  metvalues_[prefix+"_MET"].PFoMu = met.pt()/aLepton->pt();
  metvalues_[prefix+"_MET"].MuEta = aLepton->eta();
  metvalues_[prefix+"_MET"].PFSumEt = met.sumEt();

  pat::MET::Vector2 tpvec;
  tpvec.px = met.px();
  tpvec.py = met.py();
  metvalues_[prefix+"_MET"].TPphi = tpvec.phi();
  metvalues_[prefix+"_MET"].TPdPhi = reco::deltaPhi(met.phi(),aLepton->phi());
  metvalues_[prefix+"_MET"].TPpt = tpvec.pt();
  metvalues_[prefix+"_MET"].TPoGen = 1.;
  metvalues_[prefix+"_MET"].TPoMu = tpvec.pt()/aLepton->pt();

  float sumpt = 0.;
  float sumTPpt = 0.;

  for (auto lep : finalLeptons) {
    sumpt += lep->pt();
    sumTPpt += lep->pt();
  }

  metvalues_[prefix+"_MET"].TPSumEt = met.sumEt() + sumTPpt - sumpt;
  metvalues_[prefix+"_MET"].dSumEt = met.sumEt() - sumpt;
  metvalues_[prefix+"_MET"].sumEtRatio = met.sumEt()/sumpt;
  metvalues_[prefix+"_MET"].sumEtRatioTP = metvalues_[prefix+"_MET"].TPSumEt/sumTPpt;
  metvalues_[prefix+"_MET"].nLepton = static_cast<int>(finalLeptons.size());

  tree_[prefix+"_METTree"]->Fill();
}

bool MergedLeptonAnalyzer::isHighPtTrackerMuon(const reco::Muon& muon, const reco::Vertex& vtx) {
  if (!muon.isTrackerMuon())
    return false;

  bool muMatchedSt = muon.numberOfMatchedStations() > 1;
  if (!muMatchedSt) {
    if (muon.isTrackerMuon() && muon.numberOfMatchedStations() == 1) {
      if (muon.expectedNnumberOfMatchedStations() < 2 || !(muon.stationMask() == 1 || muon.stationMask() == 16) ||
          muon.numberOfMatchedRPCLayers() > 2)
        muMatchedSt = true;
    }
  }

  bool muID = muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
              muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError() / muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip =
      std::abs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && std::abs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;
}

bool MergedLeptonAnalyzer::isModifiedHEEP(const reco::GsfElectron& el,
                                          const reco::Vertex& primaryVertex,
                                          const float& trkIso,
                                          const float& ecalIso,
                                          const int& nrSatCrys,
                                          const double& rho,
                                          int& cutflow) {
  auto sq = [](const double& val) { return val*val; };
  double R = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()) + sq(el.superCluster()->z()));
  double Rt = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()));
  double etSC = (el.superCluster()->energy())*(Rt/R);
  double etaSC = el.superCluster()->eta();
  double caloIsoVal = ecalIso + el.dr03HcalDepth1TowerSumEt();
  if ( etSC < 20. )
    return false;
  cutflow++;
  if ( !el.ecalDriven() )
    return false;
  cutflow++;
  if ( !(std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.06) )
    return false;
  cutflow++;
  if ( nrSatCrys > 0 )
    return false;
  cutflow++;
  if ( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1 )
    return false;
  cutflow++;
  if ( trkIso > 5. )
    return false;
  cutflow++;
  if ( std::abs(etaSC) < 1.4442 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (1./el.energy() + 0.05) );
    bool caloIso = ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.02;

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut );
  } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (5./el.energy() + 0.05) );
    bool caloIso = ( etSC < 50. ) ? ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho )
                                  : ( caloIsoVal < 2 + 0.03*(etSC-50.) + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.05;

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut );
  }

  return false;
}

bool MergedLeptonAnalyzer::hasPassedHEEP(const reco::GsfElectron& el,
                                         const reco::Vertex& primaryVertex,
                                         const int& nrSatCrys,
                                         const double& rho,
                                         int& cutflow) {
  auto sq = [](const double& val) { return val*val; };
  double R = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()) + sq(el.superCluster()->z()));
  double Rt = std::sqrt(sq(el.superCluster()->x()) + sq(el.superCluster()->y()));
  double etSC = (el.superCluster()->energy())*(Rt/R);
  double etaSC = el.superCluster()->eta();
  double caloIsoVal = el.dr03EcalRecHitSumEt() + el.dr03HcalDepth1TowerSumEt();
  if ( etSC < 20. )
    return false;
  cutflow++;
  if ( !el.ecalDriven() )
    return false;
  cutflow++;
  if ( !(std::abs(el.deltaPhiSuperClusterTrackAtVtx()) < 0.06) )
    return false;
  cutflow++;
  if ( nrSatCrys > 0 )
    return false;
  cutflow++;
  if ( el.gsfTrack()->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) > 1 )
    return false;
  cutflow++;
  if ( el.dr03TkSumPtHEEP() > 5. )
    return false;
  cutflow++;
  if ( std::abs(etaSC) < 1.4442 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (1./el.energy() + 0.05) );
    bool caloIso = ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.02;
    bool dEtaInCut = std::abs( el.deltaEtaSeedClusterTrackAtVtx() ) < 0.004;
    bool ssCut = ( el.full5x5_e2x5Max()/el.full5x5_e5x5() > 0.94 ||
                   el.full5x5_e1x5()/el.full5x5_e5x5() > 0.83 );

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;
    if (dEtaInCut)
      cutflow++;
    if (ssCut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut && dEtaInCut && ssCut );
  } else if ( std::abs(etaSC) > 1.566 && std::abs(etaSC) < 2.5 ) {
    bool HoEcut = ( el.full5x5_hcalOverEcal() < (5./el.energy() + 0.05) );
    bool caloIso = ( etSC < 50. ) ? ( caloIsoVal < 2 + 0.03*etSC + 0.28*rho )
                                  : ( caloIsoVal < 2 + 0.03*(etSC-50.) + 0.28*rho );
    bool dxycut = std::abs( el.gsfTrack()->dxy(primaryVertex.position()) ) < 0.05;
    bool dEtaInCut = std::abs( el.deltaEtaSeedClusterTrackAtVtx() ) < 0.006;
    bool ssCut = ( el.full5x5_sigmaIetaIeta() < 0.03 );

    if (HoEcut)
      cutflow++;
    if (caloIso)
      cutflow++;
    if (dxycut)
      cutflow++;
    if (dEtaInCut)
      cutflow++;
    if (ssCut)
      cutflow++;

    return ( caloIso && HoEcut && dxycut && dEtaInCut && ssCut );
  }

  return false;
}

DEFINE_FWK_MODULE(MergedLeptonAnalyzer);
