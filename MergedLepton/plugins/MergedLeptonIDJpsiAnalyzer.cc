#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"

#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include "RecoVertex/KinematicFitPrimitives/interface/KinematicState.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "TrackingTools/GeomPropagators/interface/StateOnTrackerBound.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateOnSurface.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "RecoEgamma/EgammaElectronAlgos/interface/ElectronUtilities.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

// produce TTree for merged electron training with H->AA->4e events

class MergedLeptonIDJpsiAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonIDJpsiAnalyzer(const edm::ParameterSet&);
  virtual ~MergedLeptonIDJpsiAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  double extrapolateToSC(const pat::Electron& aEle,
                         const reco::TrackBase& addTrk,
                         const reco::BeamSpot& beamSpot,
                         const edm::EventSetup& iSetup) const;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dEtaInSeed2ndToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPhiInSC2ndToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaTrackToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaCaloToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> normDParaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5covIeIeToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5covIeIpToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dEtaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dPhiInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5EnergyToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedPFcandToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> ESrecHitToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EErecHitToken_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  edm::ConsumesCollector collector_ = consumesCollector();
  
  const std::vector<std::string> trigList_;

  const edm::FileInPath purwgtPath_;

  const bool isMC_;

  const double ptThresTag_;
  const double IPthresTag_;
  const double dzThres_;
  const double probThres_;
  const double ptThresK_;
  const double d0Thres_;
  const double cosAlpha2dThres_;
  const edm::ESGetToken<MagneticField, IdealMagneticFieldRecord> magneticToken_;
  const edm::ESGetToken<GeometricSearchTracker, TrackerRecoGeometryRecord> geotrkToken_;
  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken_;

  // PDG mass & error
  const double elmass_ = 0.0005109989461;
  const double elmassErr_ = 0.0000000000031;
  const double pionMass_ = 0.13957039;
  const double pionMassErr_ = 0.00000018;
  const double kaonMass_ = 0.493677;
  const double kaonMassErr_ = 0.000016;
  const double jpsiMass_ = 3.0969;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  TTree* EStree_ = nullptr;
  std::vector<float> ESenergy;
  std::vector<float> EStime;
  std::vector<int> ESsix;
  std::vector<int> ESsiy;
  std::vector<int> ESzside;
  std::vector<int> ESstrip;
  std::vector<int> ESplane;
  std::vector<float> ESx;
  std::vector<float> ESy;
  std::vector<float> ESz;

  TTree* EEtree_ = nullptr;
  std::vector<float> EEenergy;
  std::vector<float> EEtime;
  std::vector<int> EEzside;
  std::vector<int> EEix;
  std::vector<int> EEiy;
  std::vector<float> EEx;
  std::vector<float> EEy;
  std::vector<float> EEz;

  TTree* tree_ = nullptr;
  float invM_ = -1.;
  float u5x5Et_ = -1.;
  float et_ = -1.;
  float wgt_ = 0.;
  int passME_ = -1;
  int passModHeep_ = -1;
  int passTrig25_ = -1;
  int passTrig25unseeded_ = -1;
  int passTrig25HE_ = -1.;
  int passTrig25caloIdL_ = -1;
  int passTrig33_ = -1;
  int passTrig33unseeded_ = -1;

  TTree* Btree_ = nullptr;
  float BinvMtrk_ = -1.;
  float BinvMcalo_ = -1.;
  float Bu5x5Et_ = -1.;
  float Bwgt_ = 0.;
  int BpassME_ = -1;
  int BpassModHeep_ = -1;
  float ptK_ = -1.;
  float etaK_ = std::numeric_limits<float>::max();
  float phiK_ = std::numeric_limits<float>::max();
  float u5x5En_ = -1.;
  float etaJpsi_ = std::numeric_limits<float>::max();
  float phiJpsi_ = std::numeric_limits<float>::max();

public:
  struct dielectron {
    dielectron(const KinematicState& refitEle1,
               const KinematicState& refitEle2,
               const pat::ElectronRef& aEle,
               const reco::TrackBase& bEle,
               const pat::MuonRef& trigMu)
    : refitFirstEle(refitEle1),
      refitSecondEle(refitEle2),
      firstEle(aEle),
      secondEle(bEle),
      trigMuon(trigMu) {}

    KinematicState refitFirstEle;
    KinematicState refitSecondEle;
    pat::ElectronRef firstEle;
    reco::TrackBase secondEle;
    pat::MuonRef trigMuon;
  };

  struct dielectronFit {
    dielectronFit(const KinematicState& dielState,
                  const dielectron& diel,
                  const double achi2,
                  const double andof)
    : dielectronState(dielState),
      dielec(diel),
      chi2(achi2),
      ndof(andof) {}

    KinematicState dielectronState;
    dielectron dielec;
    double chi2;
    double ndof;

    bool operator < (const dielectronFit& other) const {
      return TMath::Prob(chi2,static_cast<int>(std::rint(ndof))) > TMath::Prob(other.chi2,static_cast<int>(std::rint(other.ndof)));
    };
  };

  struct decaychain {
    decaychain(const dielectron& diel,
               const KinematicState& refitCand,
               const pat::PackedCandidateRef& cand1,
               const KinematicState& bmeson,
               const float bmesonChi2,
               const float bmesonNdof,
               const float d0,
               const float cos)
    : dielec(diel),
      refitThirdTrk(refitCand),
      thirdTrk(cand1),
      Bmeson(bmeson),
      BmesonChi2(bmesonChi2),
      BmesonNdof(bmesonNdof),
      d0thirdTrk(d0),
      cosAlpha2d(cos) {}

    dielectron dielec;
    KinematicState refitThirdTrk;
    pat::PackedCandidateRef thirdTrk;
    KinematicState Bmeson;
    float BmesonChi2;
    float BmesonNdof;
    float d0thirdTrk;
    float cosAlpha2d;

    bool operator < (const decaychain& other) const {
      return TMath::Prob(BmesonChi2,static_cast<int>(std::rint(BmesonNdof))) > TMath::Prob(other.BmesonChi2,static_cast<int>(std::rint(other.BmesonNdof)));
    };
  };
};

MergedLeptonIDJpsiAnalyzer::MergedLeptonIDJpsiAnalyzer(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
dEtaInSeed2ndToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dEtaInSeed2nd"))),
dPhiInSC2ndToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPhiInSC2nd"))),
alphaTrackToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaTrack"))),
alphaCaloToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaCalo"))),
normDParaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("normalizedDParaIn"))),
union5x5covIeIeToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5covIeIe"))),
union5x5covIeIpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5covIeIp"))),
union5x5dEtaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dEtaIn"))),
union5x5dPhiInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dPhiIn"))),
union5x5EnergyToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5Energy"))),
packedPFcandToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedPFcand"))),
ESrecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ESrecHits"))),
EErecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EErecHits"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
isMC_(iConfig.getParameter<bool>("isMC")),
ptThresTag_(iConfig.getParameter<double>("ptThresTag")),
IPthresTag_(iConfig.getParameter<double>("IPthresTag")),
dzThres_(iConfig.getParameter<double>("dzThres")),
probThres_(iConfig.getParameter<double>("probThres")),
ptThresK_(iConfig.getParameter<double>("ptThresK")),
d0Thres_(iConfig.getParameter<double>("d0Thres")),
cosAlpha2dThres_(iConfig.getParameter<double>("cosAlpha2dThres")),
magneticToken_(esConsumes()),
geotrkToken_(esConsumes()),
ttbToken_(esConsumes(edm::ESInputTag("","TransientTrackBuilder"))),
geometryToken_(esConsumes())
{
  std::cout<<"hello"<<std::endl;
  std::cout<<"hello2"<<std::endl;
  
  usesResource("TFileService");
}

void MergedLeptonIDJpsiAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",30,0.,30.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",200,-1.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);
  histo1d_["PUsummary"] = fs->make<TH1D>("PUsummary","PUsummary",99,0.,99.);

  EStree_ = fs->make<TTree>("ESTree","ESTree");
  EStree_->Branch("ESenergy",&ESenergy,32000,0);
  EStree_->Branch("EStime",&EStime,32000,0);
  EStree_->Branch("ESsix",&ESsix,32000,0);
  EStree_->Branch("ESsiy",&ESsiy,32000,0);
  EStree_->Branch("ESzside",&ESzside,32000,0);
  EStree_->Branch("ESstrip",&ESstrip,32000,0);
  EStree_->Branch("ESplane",&ESplane,32000,0);
  EStree_->Branch("ESx",&ESx,32000,0);
  EStree_->Branch("ESy",&ESy,32000,0);
  EStree_->Branch("ESz",&ESz,32000,0);

  EEtree_ = fs->make<TTree>("EETree","EETree");
  EEtree_->Branch("EEenergy",&EEenergy,32000,0);
  EEtree_->Branch("EEtime",&EEtime,32000,0);
  EEtree_->Branch("EEix",&EEix,32000,0);
  EEtree_->Branch("EEiy",&EEiy,32000,0);
  EEtree_->Branch("EEzside",&EEzside,32000,0);
  EEtree_->Branch("EEx",&EEx,32000,0);
  EEtree_->Branch("EEy",&EEy,32000,0);
  EEtree_->Branch("EEz",&EEz,32000,0);

  tree_ = fs->make<TTree>("dielTree","dielTree");
  tree_->Branch("invM",&invM_,"invM/F");
  tree_->Branch("u5x5Et",&u5x5Et_,"u5x5Et/F");
  tree_->Branch("Et",&et_,"Et/F");
  tree_->Branch("wgt",&wgt_,"wgt/F");
  tree_->Branch("passME",&passME_,"passME/I");
  tree_->Branch("passModHeep",&passModHeep_,"passModHeep/I");
  tree_->Branch("passTrig25",&passTrig25_,"passTrig25/I");
  tree_->Branch("passTrig25unseeded",&passTrig25unseeded_,"passTrig25unseeded/I");
  tree_->Branch("passTrig25HE",&passTrig25HE_,"passTrig25HE/I");
  tree_->Branch("passTrig25caloIdL",&passTrig25caloIdL_,"passTrig25caloIdL/I");
  tree_->Branch("passTrig33",&passTrig33_,"passTrig33/I");
  tree_->Branch("passTrig33unseeded",&passTrig33unseeded_,"passTrig33unseeded/I");

  Btree_ = fs->make<TTree>("BmesonTree","BmesonTree");
  Btree_->Branch("invMtrk",&BinvMtrk_,"invMtrk/F");
  Btree_->Branch("invMcalo",&BinvMcalo_,"invMcalo/F");
  Btree_->Branch("u5x5Et",&Bu5x5Et_,"u5x5Et/F");
  Btree_->Branch("wgt",&Bwgt_,"wgt/F");
  Btree_->Branch("passME",&BpassME_,"passME/I");
  Btree_->Branch("passModHeep",&BpassModHeep_,"passModHeep/I");
  Btree_->Branch("ptK",&ptK_,"ptK/F");
  Btree_->Branch("etaK",&etaK_,"etaK/F");
  Btree_->Branch("phiK",&phiK_,"phiK/F");
  Btree_->Branch("u5x5En",&u5x5En_,"u5x5En/F");
  Btree_->Branch("etaJpsi",&etaJpsi_,"etaJpsi/F");
  Btree_->Branch("phiJpsi",&phiJpsi_,"phiJpsi/F");

  auto create2trkHisto = [this,&fs] (const std::string& prefix) {
    histo1d_[prefix+"2trk_diel_invM"] = fs->make<TH1D>((prefix+"2trk_diel_invM").c_str(),";Mass [GeV];",1000,0.,5.);
    histo1d_[prefix+"2trk_diel_pt"] = fs->make<TH1D>((prefix+"2trk_diel_pt").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"2trk_diel_pt_preFit"] = fs->make<TH1D>((prefix+"2trk_diel_pt_preFit").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"2trk_diel_etSC"] = fs->make<TH1D>((prefix+"2trk_diel_etSC").c_str(),";E_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"2trk_diel_etU5x5"] = fs->make<TH1D>((prefix+"2trk_diel_etU5x5").c_str(),";E_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"2trk_diel_normChi2"] = fs->make<TH1D>((prefix+"2trk_diel_normChi2").c_str(),";#chi^{2}/ndof;",100,0.,10.);
    histo1d_[prefix+"2trk_diel_rap"] = fs->make<TH1D>((prefix+"2trk_diel_rap").c_str(),";rapidity;",200,-2.5,2.5);
    histo1d_[prefix+"2trk_diel_dR"] = fs->make<TH1D>((prefix+"2trk_diel_dR").c_str(),";#Delta R;",500,0.,0.5);
    histo1d_[prefix+"2trk_diel_dEta"] = fs->make<TH1D>((prefix+"2trk_diel_dEta").c_str(),";#Delta#eta(e1,e2);",400,-0.4,0.4);
    histo1d_[prefix+"2trk_diel_dPhi"] = fs->make<TH1D>((prefix+"2trk_diel_dPhi").c_str(),";#Delta#phi(e1,e2);",400,-0.4,0.4);
    histo1d_[prefix+"2trk_diel_prob"] = fs->make<TH1D>((prefix+"2trk_diel_prob").c_str(),";Prob;",200,0.,1.);
    histo1d_[prefix+"2trk_diel_alphaCalo"] = fs->make<TH1D>((prefix+"2trk_diel_alphaCalo").c_str(),";#alpha_{calo};",360,-1.8,1.8);
    histo1d_[prefix+"2trk_diel_alphaTrk"] = fs->make<TH1D>((prefix+"2trk_diel_alphaTrk").c_str(),";#alpha_{trk};",360,-1.8,1.8);
    histo1d_[prefix+"2trk_diel_dPerpIn"] = fs->make<TH1D>((prefix+"2trk_diel_dPerpIn").c_str(),";#Delta u_{in};",400,-0.05,0.05);
    histo1d_[prefix+"2trk_diel_normDParaIn"] = fs->make<TH1D>((prefix+"2trk_diel_normDParaIn").c_str(),";#Delta v_{in}/#Delta R;",200,-0.5,1.5);
    histo1d_[prefix+"2trk_diel_u5x5dEtaIn"] = fs->make<TH1D>((prefix+"2trk_diel_u5x5dEtaIn").c_str(),";#Delta #eta_{in}(u5x5);",500,-0.1,0.1);
    histo1d_[prefix+"2trk_diel_u5x5dPhiIn"] = fs->make<TH1D>((prefix+"2trk_diel_u5x5dPhiIn").c_str(),";#Delta #phi_{in}(u5x5);",500,-0.1,0.1);
    histo1d_[prefix+"2trk_diel_covIeIe"] = fs->make<TH1D>((prefix+"2trk_diel_covIeIe").c_str(),";#sigma_{i#eta i#eta};",400,0.,2.);
    histo1d_[prefix+"2trk_diel_covIeIp"] = fs->make<TH1D>((prefix+"2trk_diel_covIeIp").c_str(),";#sigma_{i#eta i#phi};",400,-2.,2.);
    histo1d_[prefix+"2trk_diel_HoE"] = fs->make<TH1D>((prefix+"2trk_diel_HoE").c_str(),";H/E;",200,0.,0.2);
    histo1d_[prefix+"2trk_diel_sigIeIe"] = fs->make<TH1D>((prefix+"2trk_diel_sigIeIe").c_str(),";#sigma_{i#eta i#eta};",300,0.,0.03);
    histo1d_[prefix+"2trk_diel_trackIso"] = fs->make<TH1D>((prefix+"2trk_diel_trackIso").c_str(),";#Sigma p_{T iso};",300,0.,30.);
    histo1d_[prefix+"2trk_pt_1st"] = fs->make<TH1D>((prefix+"2trk_pt_1st").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"2trk_eta_1st"] = fs->make<TH1D>((prefix+"2trk_eta_1st").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"2trk_dzTrig_1st"] = fs->make<TH1D>((prefix+"2trk_dzTrig_1st").c_str(),";#Delta z(trk, #mu_{trig});",200,-1.,1.);
    histo1d_[prefix+"2trk_dxy_1st"] = fs->make<TH1D>((prefix+"2trk_dxy_1st").c_str(),";d_{xy};",200,-1.,1.);
    histo1d_[prefix+"2trk_IPsig_1st"] = fs->make<TH1D>((prefix+"2trk_IPsig_1st").c_str(),";|d_{xy}/#sigma|;",150,0.,15.);
    histo1d_[prefix+"2trk_pt_2nd"] = fs->make<TH1D>((prefix+"2trk_pt_2nd").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"2trk_eta_2nd"] = fs->make<TH1D>((prefix+"2trk_eta_2nd").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"2trk_dzTrig_2nd"] = fs->make<TH1D>((prefix+"2trk_dzTrig_2nd").c_str(),";#Delta z(trk, #mu_{trig});",200,-1.,1.);
    histo1d_[prefix+"2trk_dxy_2nd"] = fs->make<TH1D>((prefix+"2trk_dxy_2nd").c_str(),";d_{xy};",200,-1.,1.);
    histo1d_[prefix+"2trk_IPsig_2nd"] = fs->make<TH1D>((prefix+"2trk_IPsig_2nd").c_str(),";|d_{xy}/#sigma|;",150,0.,15.);
    histo1d_[prefix+"2trk_mva_HasTrkEB"] = fs->make<TH1D>((prefix+"2trk_mva_HasTrkEB").c_str(),"MVA score",200,-1.,1.);
  };

  create2trkHisto("");
  create2trkHisto("PassModifiedHeep_");
  create2trkHisto("FailModifiedHeep_");
  create2trkHisto("PassMergedEle_");
  create2trkHisto("FailMergedEle_");

  create2trkHisto("Jpsi");
  create2trkHisto("JpsiPassModifiedHeep_");
  create2trkHisto("JpsiFailModifiedHeep_");
  create2trkHisto("JpsiPassMergedEle_");
  create2trkHisto("JpsiFailMergedEle_");

  create2trkHisto("JpsiPt20to25");
  create2trkHisto("JpsiPt20to25PassModifiedHeep_");
  create2trkHisto("JpsiPt20to25FailModifiedHeep_");
  create2trkHisto("JpsiPt20to25PassMergedEle_");
  create2trkHisto("JpsiPt20to25FailMergedEle_");

  create2trkHisto("JpsiPt25to30");
  create2trkHisto("JpsiPt25to30PassModifiedHeep_");
  create2trkHisto("JpsiPt25to30FailModifiedHeep_");
  create2trkHisto("JpsiPt25to30PassMergedEle_");
  create2trkHisto("JpsiPt25to30FailMergedEle_");

  create2trkHisto("JpsiPt30to35");
  create2trkHisto("JpsiPt30to35PassModifiedHeep_");
  create2trkHisto("JpsiPt30to35FailModifiedHeep_");
  create2trkHisto("JpsiPt30to35PassMergedEle_");
  create2trkHisto("JpsiPt30to35FailMergedEle_");

  create2trkHisto("JpsiPt35to40");
  create2trkHisto("JpsiPt35to40PassModifiedHeep_");
  create2trkHisto("JpsiPt35to40FailModifiedHeep_");
  create2trkHisto("JpsiPt35to40PassMergedEle_");
  create2trkHisto("JpsiPt35to40FailMergedEle_");

  create2trkHisto("JpsiPt40to50");
  create2trkHisto("JpsiPt40to50PassModifiedHeep_");
  create2trkHisto("JpsiPt40to50FailModifiedHeep_");
  create2trkHisto("JpsiPt40to50PassMergedEle_");
  create2trkHisto("JpsiPt40to50FailMergedEle_");

  create2trkHisto("JpsiPt50to60");
  create2trkHisto("JpsiPt50to60PassModifiedHeep_");
  create2trkHisto("JpsiPt50to60FailModifiedHeep_");
  create2trkHisto("JpsiPt50to60PassMergedEle_");
  create2trkHisto("JpsiPt50to60FailMergedEle_");

  create2trkHisto("JpsiPt60toInf");
  create2trkHisto("JpsiPt60toInfPassModifiedHeep_");
  create2trkHisto("JpsiPt60toInfFailModifiedHeep_");
  create2trkHisto("JpsiPt60toInfPassMergedEle_");
  create2trkHisto("JpsiPt60toInfFailMergedEle_");

  auto create3trkHisto = [this,&fs] (const std::string& prefix) {
    histo1d_[prefix+"3trk_diel_invM"] = fs->make<TH1D>((prefix+"3trk_diel_invM").c_str(),";Mass [GeV];",1000,0.,5.);
    histo1d_[prefix+"3trk_diel_invM_wide"] = fs->make<TH1D>((prefix+"3trk_diel_invM_wide").c_str(),";Mass [GeV];",500,0.,50.);
    histo1d_[prefix+"3trk_diel_pt"] = fs->make<TH1D>((prefix+"3trk_diel_pt").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"3trk_diel_etSC"] = fs->make<TH1D>((prefix+"3trk_diel_etSC").c_str(),";E_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"3trk_diel_etU5x5"] = fs->make<TH1D>((prefix+"3trk_diel_etU5x5").c_str(),";E_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"3trk_diel_rap"] = fs->make<TH1D>((prefix+"3trk_diel_rap").c_str(),";rapidity;",200,-2.5,2.5);
    histo1d_[prefix+"3trk_diel_dR"] = fs->make<TH1D>((prefix+"3trk_diel_dR").c_str(),";#Delta R(e1,e2);",500,0.,0.5);
    histo1d_[prefix+"3trk_diel_dEta"] = fs->make<TH1D>((prefix+"3trk_diel_dEta").c_str(),";#Delta#eta(e1,e2);",400,-0.4,0.4);
    histo1d_[prefix+"3trk_diel_dPhi"] = fs->make<TH1D>((prefix+"3trk_diel_dPhi").c_str(),";#Delta#phi(e1,e2);",400,-0.4,0.4);
    histo1d_[prefix+"3trk_diel_alphaCalo"] = fs->make<TH1D>((prefix+"3trk_diel_alphaCalo").c_str(),";#alpha_{calo};",360,-1.8,1.8);
    histo1d_[prefix+"3trk_diel_alphaTrk"] = fs->make<TH1D>((prefix+"3trk_diel_alphaTrk").c_str(),";#alpha_{trk};",360,-1.8,1.8);

    histo1d_[prefix+"3trk_diel_dPerpIn"] = fs->make<TH1D>((prefix+"3trk_diel_dPerpIn").c_str(),";#Delta u_{in};",500,-0.1,0.1);
    histo1d_[prefix+"3trk_diel_normDParaIn"] = fs->make<TH1D>((prefix+"3trk_diel_normDParaIn").c_str(),";#Delta v_{in}/#Delta R;",200,-0.5,1.5);
    histo1d_[prefix+"3trk_diel_u5x5dEtaIn"] = fs->make<TH1D>((prefix+"3trk_diel_u5x5dEtaIn").c_str(),";#Delta #eta_{in}(u5x5);",500,-0.1,0.1);
    histo1d_[prefix+"3trk_diel_u5x5dPhiIn"] = fs->make<TH1D>((prefix+"3trk_diel_u5x5dPhiIn").c_str(),";#Delta #phi_{in}(u5x5);",500,-0.1,0.1);
    histo1d_[prefix+"3trk_diel_covIeIe"] = fs->make<TH1D>((prefix+"3trk_diel_covIeIe").c_str(),";#sigma_{i#eta i#eta};",400,0.,2.);
    histo1d_[prefix+"3trk_diel_covIeIp"] = fs->make<TH1D>((prefix+"3trk_diel_covIeIp").c_str(),";#sigma_{i#eta i#phi};",400,-2.,2.);

    histo1d_[prefix+"3trk_diel_HoE"] = fs->make<TH1D>((prefix+"3trk_diel_HoE").c_str(),";H/E;",200,0.,0.2);
    histo1d_[prefix+"3trk_diel_sigIeIe"] = fs->make<TH1D>((prefix+"3trk_diel_sigIeIe").c_str(),";#sigma_{i#eta i#eta};",300,0.,0.03);
    histo1d_[prefix+"3trk_pt_1st"] = fs->make<TH1D>((prefix+"3trk_pt_1st").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"3trk_eta_1st"] = fs->make<TH1D>((prefix+"3trk_eta_1st").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"3trk_dzTrig_1st"] = fs->make<TH1D>((prefix+"3trk_dzTrig_1st").c_str(),";#Delta z(trk, #mu_{trig});",200,-0.5,0.5);
    histo1d_[prefix+"3trk_IPsig_1st"] = fs->make<TH1D>((prefix+"3trk_IPsig_1st").c_str(),";|d_{xy}/#sigma|;",150,0.,15.);
    histo1d_[prefix+"3trk_pt_2nd"] = fs->make<TH1D>((prefix+"3trk_pt_2nd").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"3trk_eta_2nd"] = fs->make<TH1D>((prefix+"3trk_eta_2nd").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"3trk_dzTrig_2nd"] = fs->make<TH1D>((prefix+"3trk_dzTrig_2nd").c_str(),";#Delta z(trk, #mu_{trig});",200,-0.5,0.5);
    histo1d_[prefix+"3trk_IPsig_2nd"] = fs->make<TH1D>((prefix+"3trk_IPsig_2nd").c_str(),";|d_{xy}/#sigma|;",150,0.,15.);

    histo1d_[prefix+"3trk_pt_3rd"] = fs->make<TH1D>((prefix+"3trk_pt_3rd").c_str(),";p_{T} [GeV];",200,0.,100.);
    histo1d_[prefix+"3trk_eta_3rd"] = fs->make<TH1D>((prefix+"3trk_eta_3rd").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"3trk_drAtSC_3rd"] = fs->make<TH1D>((prefix+"3trk_drAtSC_3rd").c_str(),";#Delta R_{calo}(SC,h);",640,0.,3.2);
    histo1d_[prefix+"3trk_dzTrig_3rd"] = fs->make<TH1D>((prefix+"3trk_dzTrig_3rd").c_str(),";#Delta z(trk, #mu_{trig});",200,-0.5,0.5);
    histo1d_[prefix+"3trk_IPsig_3rd"] = fs->make<TH1D>((prefix+"3trk_IPsig_3rd").c_str(),";|d_{xy}/#sigma|;",150,0.,15.);
    histo1d_[prefix+"3trk_d0_3rd"] = fs->make<TH1D>((prefix+"3trk_d0_3rd").c_str(),";d_{3D}(ee,h);",200,-0.1,0.1);

    histo1d_[prefix+"3trk_B_invM"] = fs->make<TH1D>((prefix+"3trk_B_invM").c_str(),";Mass [GeV];",1000,0.,10.);
    histo1d_[prefix+"3trk_B_invM_wide"] = fs->make<TH1D>((prefix+"3trk_B_invM_wide").c_str(),";Mass [GeV];",500,0.,50.);
    histo1d_[prefix+"3trk_B_invM_u5x5"] = fs->make<TH1D>((prefix+"3trk_B_invM_u5x5").c_str(),";Mass [GeV];",1000,0.,10.);
    histo1d_[prefix+"3trk_B_pt"] = fs->make<TH1D>((prefix+"3trk_B_pt").c_str(),";p_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"3trk_B_normChi2"] = fs->make<TH1D>((prefix+"3trk_B_normChi2").c_str(),";#chi^{2}/ndof;",100,0.,10.);
    histo1d_[prefix+"3trk_B_prob"] = fs->make<TH1D>((prefix+"3trk_B_prob").c_str(),";Prob;",200,0.,1.);
    histo1d_[prefix+"3trk_B_rap"] = fs->make<TH1D>((prefix+"3trk_B_rap").c_str(),";rapidity;",200,-2.5,2.5);
    histo1d_[prefix+"3trk_B_dR"] = fs->make<TH1D>((prefix+"3trk_B_dR").c_str(),";#Delta R(ee,h);",500,0.,1.0);
    histo1d_[prefix+"3trk_B_dR_wide"] = fs->make<TH1D>((prefix+"3trk_B_dR_wide").c_str(),";#Delta R(ee,h);",640,0.,6.4);
    histo1d_[prefix+"3trk_B_cosAlpha2d"] = fs->make<TH1D>((prefix+"3trk_B_cosAlpha2d").c_str(),";cos_{2D}(#alpha);",100,0.95,1.);

    histo1d_[prefix+"3trk_Pi1Pi2_invM"] = fs->make<TH1D>((prefix+"3trk_Pi1Pi2_invM").c_str(),";Mass [GeV];",800,0.,8.);
    histo1d_[prefix+"3trk_Pi1K2_invM"] = fs->make<TH1D>((prefix+"3trk_Pi1K2_invM").c_str(),";Mass [GeV];",800,0.,8.);
    histo1d_[prefix+"3trk_K1Pi2_invM"] = fs->make<TH1D>((prefix+"3trk_K1Pi2_invM").c_str(),";Mass [GeV];",800,0.,8.);

    histo1d_[prefix+"3trk_e1K_invM"] = fs->make<TH1D>((prefix+"3trk_e1K_invM").c_str(),";Mass [GeV];",800,0.,8.);
    histo1d_[prefix+"3trk_e2K_invM"] = fs->make<TH1D>((prefix+"3trk_e2K_invM").c_str(),";Mass [GeV];",800,0.,8.);

    histo2d_[prefix+"3trk_dalitz_E1X_E2X"] = fs->make<TH2F>((prefix+"3trk_dalitz_E1X_E2X").c_str(),";m^{2}(e1,X);m^{2}(e2,X);",200,0.,50.,200,0.,50.);
    histo2d_[prefix+"3trk_dalitz_E1E2_E2X"] = fs->make<TH2F>((prefix+"3trk_dalitz_E1E2_E2X").c_str(),";m^{2}(e1,e2);m^{2}(e2,X);",200,0.,15.,200,0.,50.);
    histo2d_[prefix+"3trk_dalitz_E1E2_E1X"] = fs->make<TH2F>((prefix+"3trk_dalitz_E1E2_E1X").c_str(),";m^{2}(e1,e2);m^{2}(e1,X);",200,0.,15.,200,0.,50.);
  };

  create3trkHisto("");
  create3trkHisto("PassModifiedHeep_");
  create3trkHisto("FailModifiedHeep_");
  create3trkHisto("PassMergedEle_");
  create3trkHisto("FailMergedEle_");
  create3trkHisto("PassModifiedHeepPassMergedEle_");
}

void MergedLeptonIDJpsiAnalyzer::endJob() {
  purwgtFile_->Close();
}

double MergedLeptonIDJpsiAnalyzer::extrapolateToSC(const pat::Electron& aEle,
                                               const reco::TrackBase& addTrk,
                                               const reco::BeamSpot& beamSpot,
                                               const edm::EventSetup& iSetup) const {
  // Get magField & tracker
  //edm::ESHandle<MagneticField> magFieldHandle;
  //iSetup.get<IdealMagneticFieldRecord>().get(magFieldHandle);
  const MagneticField* magFieldHandle = &iSetup.getData(magneticToken_);

  //edm::ESHandle<GeometricSearchTracker> trackerSearchHandle;
  //iSetup.get<TrackerRecoGeometryRecord>().get(trackerSearchHandle);
  const GeometricSearchTracker* trackerSearchHandle = &iSetup.getData(geotrkToken_);

  const auto& pixelBarrelLayers = trackerSearchHandle->pixelBarrelLayers();
  BarrelDetLayer* innermostLayer = nullptr;
  float innermostRadius = std::numeric_limits<float>::max();

  for (const auto* alayer : pixelBarrelLayers) {
    float aradius = alayer->surface().rSpan().first;
    if ( aradius < innermostRadius ) {
      innermostRadius = aradius;
      innermostLayer = const_cast<BarrelDetLayer*>(alayer);
    }
  }

  // GlobalTag matters here (tracker alignment)
  FreeTrajectoryState freestate(GlobalTrajectoryParameters(GlobalPoint(addTrk.referencePoint().x(),
                                                                       addTrk.referencePoint().y(),
                                                                       addTrk.referencePoint().z()),
                                                           GlobalVector(addTrk.momentum().x(),
                                                                        addTrk.momentum().y(),
                                                                        addTrk.momentum().z()),
                                                           addTrk.charge(),
                                                           magFieldHandle),
                                CurvilinearTrajectoryError(addTrk.covariance()));
  auto propagator = std::make_unique<AnalyticalPropagator>(magFieldHandle);
  auto extrapolator = std::make_unique<TransverseImpactPointExtrapolator>(*propagator);
  TrajectoryStateOnSurface innTSOS = propagator->propagate(freestate,innermostLayer->surface());
  StateOnTrackerBound stateOnBound(propagator.get());
  TrajectoryStateOnSurface outTSOS = stateOnBound(freestate);
  double dEtaSC = std::numeric_limits<float>::max();
  double dPhiSC = std::numeric_limits<float>::max();

  if ( innTSOS.isValid() && outTSOS.isValid() ) {
    TrajectoryStateOnSurface sclTSOS = extrapolator->extrapolate(*(innTSOS.freeState()), // with TSOS assert fails (not a real measurement)
                                                                 GlobalPoint(aEle.superCluster()->x(),
                                                                             aEle.superCluster()->y(),
                                                                             aEle.superCluster()->z()));
    if (!sclTSOS.isValid())
      sclTSOS = outTSOS;

    // single state (KF track)
    GlobalPoint sclPos = sclTSOS.globalPosition();

    EleRelPointPair scAtVtx(aEle.superCluster()->position(),sclPos,beamSpot.position());
    dPhiSC = scAtVtx.dPhi();
    dEtaSC = scAtVtx.dEta();
  }

  return std::sqrt( dEtaSC*dEtaSC + dPhiSC*dPhiSC );
}

void MergedLeptonIDJpsiAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
        auto npu = apu->getTrueNumInteractions();
        aWeight *= purwgt_->GetBinContent( purwgt_->FindBin(apu->getTrueNumInteractions()) );
        histo1d_["PUsummary"]->Fill( static_cast<float>(npu)+0.5, aWeight );

        break;
      }
    }
  }
  
  const CaloGeometry* caloGeom = &iSetup.getData(geometryToken_);

  edm::Handle<EcalRecHitCollection> ESrecHitHandle;
  iEvent.getByToken(ESrecHitToken_, ESrecHitHandle);
  edm::Handle<EcalRecHitCollection> EErecHitHandle;
  iEvent.getByToken(EErecHitToken_, EErecHitHandle);
  ESenergy.clear();
  EStime.clear();
  ESsix.clear();
  ESsiy.clear();
  ESstrip.clear();
  ESplane.clear();
  ESzside.clear();
  ESx.clear();
  ESy.clear();
  ESz.clear();
  for (auto& silicon : *ESrecHitHandle){
    ESenergy.push_back(silicon.energy());
    EStime.push_back(silicon.time());
    auto idES = ESDetId(silicon.detid());
    const auto& siliconGeo = caloGeom->getGeometry(silicon.detid());

    ESsix.push_back(idES.six());
    ESsiy.push_back(idES.siy());
    ESstrip.push_back(idES.strip());
    ESzside.push_back(idES.zside());
    ESplane.push_back(idES.plane());
    ESx.push_back((float)siliconGeo->getPosition().x());
    ESy.push_back((float)siliconGeo->getPosition().y());
    ESz.push_back((float)siliconGeo->getPosition().z());
  }
  EStree_->Fill();

  EEenergy.clear();
  EEtime.clear();
  EEix.clear();
  EEiy.clear();
  EEzside.clear();
  EEx.clear();
  EEy.clear();
  EEz.clear();

  for (auto& crystal : *EErecHitHandle){
    EEenergy.push_back(crystal.energy());
    EEtime.push_back(crystal.time());
    auto idEE = EEDetId(crystal.detid());
    const auto& crystalGeo = caloGeom->getGeometry(crystal.detid());
    EEix.push_back(idEE.ix());
    EEiy.push_back(idEE.iy());
    EEzside.push_back(idEE.zside());
    EEx.push_back((float)crystalGeo->getPosition().x());
    EEy.push_back((float)crystalGeo->getPosition().y());
    EEz.push_back((float)crystalGeo->getPosition().z());
    
  }
  EEtree_->Fill();

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);
  histo1d_["cutflow"]->Fill(0.5,aWeight);

  reco::VertexRef primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->refAt(0).castTo<reco::VertexRef>();
  else
    return;

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  std::vector<std::pair<std::string, int>> indices;
  edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;
  bool isFired_mu9 = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    std::string trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
        if (trigResultHandle.product()->accept(iTrig)) {
          isFired = true;

          if (trigList_.at(jTrig)=="HLT_Mu9_IP6_part*")
            isFired_mu9 = true;
        }
      }
    } // wanted triggers
  } // fired triggers

  if (!isFired)
    return;

  histo1d_["cutflow"]->Fill(1.5,aWeight);
  histo1d_["nPV"]->Fill( static_cast<float>(pvHandle->size())+0.5, aWeight );

  auto retrieveTrigObj = [&trigObjHandle,&trigResultHandle,&trigList,&iEvent] (const std::vector<std::string>& alist, bool useFilters=false)
                         -> std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> {
    std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjs;

    for (unsigned iTrig = 0; iTrig < trigObjHandle->size(); iTrig++) {
      const auto& trigObj = trigObjHandle->refAt(iTrig);
      auto trigObjInst = trigObjHandle->at(iTrig); // workaround for copy
      trigObjInst.unpackPathNames(trigList);
      trigObjInst.unpackFilterLabels(iEvent, *trigResultHandle);
      const auto& pathNames = trigObjInst.pathNames();

      if (useFilters) {
        for (const auto& aname : alist) {
          if (trigObjInst.hasFilterLabel(aname))
            trigObjs.push_back(trigObj);
        }
      } else {
        for (const auto name : pathNames) {
          for (unsigned int jTrig = 0; jTrig < alist.size(); jTrig++) {
            if ( name.find(alist.at(jTrig).substr(0, alist.at(jTrig).find("*"))) != std::string::npos &&
                 trigObjInst.hasPathName(name,true,true) ) {
              trigObjs.push_back(trigObj);
            }
          } // wanted triggers
        } // fired triggers
      }
    } // trigger objs

    return trigObjs;
  };

  auto matchTrigObj = [] (const math::XYZTLorentzVectorD& p4,
                          std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>>& trigObjs) {
    for (const auto& trigObj : trigObjs) {
      if ( reco::deltaR2(p4.eta(),p4.phi(),trigObj->eta(),trigObj->phi()) < 0.01 )
        return true;
    }

    return false;
  };

  auto trigObjs = retrieveTrigObj(trigList_);
  auto trigObjs_doubleEle25 = retrieveTrigObj({"hltEle25CaloIdLMWPMS2Filter"},true);
  auto trigObjs_doubleEle25unseeded = retrieveTrigObj({"hltDiEle25CaloIdLMWPMS2UnseededFilter"},true);
  auto trigObjs_doubleEle25HE = retrieveTrigObj({"hltEG25HEFilter"},true);
  auto trigObjs_doubleEle25caloIdL = retrieveTrigObj({"hltEG25CaloIdLClusterShapeFilter"},true);
  auto trigObjs_doubleEle33 = retrieveTrigObj({"hltEle33CaloIdLMWPMS2Filter"},true);
  auto trigObjs_doubleEle33unseeded = retrieveTrigObj({"hltDiEle33CaloIdLMWPMS2UnseededFilter"},true);

  // preselections motivated from R(K) analysis i.e. BPH-22-005
  // first select trigger muon
  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  if (muonHandle->size()==0)
    return;

  std::vector<pat::MuonRef> mediumMuons;

  for (unsigned iMu = 0; iMu < muonHandle->size(); iMu++) {
    const auto& aMu = muonHandle->refAt(iMu);

    if (aMu->innerTrack().isNull())
      continue;

    if (!muon::isMediumMuon(*aMu))
      continue;

    const auto& innerTrack = aMu->innerTrack();
    const double dxy = innerTrack->dxy(primaryVertex->position());
    const double dxyErr = innerTrack->dxyError(primaryVertex->position(),primaryVertex->covariance());
    const double ipSig = dxy/dxyErr;

    if ( aMu->pt() > ptThresTag_ && std::abs(aMu->eta()) < 1.5 && std::abs(ipSig) > IPthresTag_ )
      if ( isFired_mu9 || aMu->pt() > 12. )
        mediumMuons.push_back(aMu.castTo<pat::MuonRef>());
  }

  if (mediumMuons.empty())
    return;

  bool matched = false;
  std::vector<pat::MuonRef> trigMuons;

  for (const auto& aMu : mediumMuons) {
    if ( matchTrigObj(aMu->p4(),trigObjs) ) {
      trigMuons.push_back(aMu);
      matched = true;
    }
  }

  if (!matched)
    return;

  histo1d_["cutflow"]->Fill(2.5,aWeight);

  std::sort( trigMuons.begin(), trigMuons.end(), [](const pat::MuonRef& a, const pat::MuonRef& b) { return a->pt() > b->pt(); } );

  // next select one electron and one track (first GSF and then KF)
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  if ( eleHandle->size()==0 )
    return;

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<edm::ValueMap<float>> dPerpInHandle;
  iEvent.getByToken(dPerpInToken_, dPerpInHandle);

  edm::Handle<edm::ValueMap<float>> dEtaInSeed2ndHandle;
  iEvent.getByToken(dEtaInSeed2ndToken_, dEtaInSeed2ndHandle);

  edm::Handle<edm::ValueMap<float>> dPhiInSC2ndHandle;
  iEvent.getByToken(dPhiInSC2ndToken_, dPhiInSC2ndHandle);

  edm::Handle<edm::ValueMap<float>> alphaTrackHandle;
  iEvent.getByToken(alphaTrackToken_, alphaTrackHandle);

  edm::Handle<edm::ValueMap<float>> alphaCaloHandle;
  iEvent.getByToken(alphaCaloToken_, alphaCaloHandle);

  edm::Handle<edm::ValueMap<float>> normDParaInHandle;
  iEvent.getByToken(normDParaInToken_, normDParaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5covIeIeHandle;
  iEvent.getByToken(union5x5covIeIeToken_, union5x5covIeIeHandle);

  edm::Handle<edm::ValueMap<float>> union5x5covIeIpHandle;
  iEvent.getByToken(union5x5covIeIpToken_, union5x5covIeIpHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dPhiInHandle;
  iEvent.getByToken(union5x5dPhiInToken_, union5x5dPhiInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::View<pat::PackedCandidate>> packedPFcandHandle;
  iEvent.getByToken(packedPFcandToken_, packedPFcandHandle);

  
  //edm::ESHandle<TransientTrackBuilder> TTbuilder;
  //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTbuilder);
  const TransientTrackBuilder& TTbuilder = iSetup.getData(ttbToken_);

  for (unsigned iEle = 0; iEle < eleHandle->size(); iEle++) {
    const auto& aEle = eleHandle->refAt(iEle);
    const auto& orgGsfTrk = aEle->gsfTrack();

    const auto& addGsfTrk = (*addGsfTrkHandle)[aEle];
    const auto& addPackedCand = (*addPackedCandHandle)[aEle];
    const float alphaTrack = (*alphaTrackHandle)[aEle];
    const float alphaCalo = (*alphaCaloHandle)[aEle];
    const float trackIso = (*trkIsoMapHandle)[aEle];

    const float eta1stGSF = -(aEle->deltaEtaSeedClusterTrackAtVtx() - aEle->superCluster()->seed()->eta());
    const float u5x5Eta = (*union5x5dEtaInHandle)[aEle] + eta1stGSF;
    const float u5x5Et = (*union5x5EnergyHandle)[aEle]/std::cosh(u5x5Eta);

    if ( u5x5Et < 10. || std::abs(aEle->superCluster()->eta()) > 1.4442 )
      continue;

    int32_t bitmap = aEle->userInt("modifiedHeepElectronID");
    int32_t mask = 0x00000381; // = 0011 1000 0001 - 1st for min Et, 7th for trk iso, 8th for EM+HadD1 iso, 9th for dxy
    int32_t pass = bitmap | mask;
    bool passMaskedId = pass==0x00000FFF; // HEEP ID has 12 cuts

    if (!aEle->ecalDriven())
      continue;

    if (addGsfTrk==orgGsfTrk && addPackedCand.isNull())
      continue;

    // first look for GSF and then packedPFcand
    const reco::Track* addTrk = addGsfTrk.get();

    if ( addGsfTrk==orgGsfTrk && addPackedCand.isNonnull() )
      addTrk = addPackedCand->bestTrack();

    if ( orgGsfTrk->charge()*addTrk->charge() > 0 )
      continue;

    std::vector<dielectronFit> diels;
    std::vector<dielectronFit> Jpsis;

    for (const auto& trigMuon : trigMuons) {
      KinematicParticleFactoryFromTransientTrack kinFactory;
      KinematicParticleVertexFitter kinFitter;
      float elmassErr = elmassErr_; // why this thing requires lvalue? :(
      const reco::TransientTrack firstEle = TTbuilder.build(orgGsfTrk);
      const reco::TransientTrack secondEle = TTbuilder.build(addTrk);
      double chi = 0.;
      double ndf = 0.;

      std::vector<RefCountedKinematicParticle> dielvec;
      dielvec.push_back(kinFactory.particle(firstEle,elmass_,chi,ndf,elmassErr));
      dielvec.push_back(kinFactory.particle(secondEle,elmass_,chi,ndf,elmassErr));
      RefCountedKinematicTree dielectronTree = kinFitter.fit(dielvec);

      if ( std::abs( orgGsfTrk->vz() - trigMuon->innerTrack()->vz() ) > dzThres_ ||
           std::abs( addTrk->vz() - trigMuon->innerTrack()->vz() ) > dzThres_ )
        continue;

      if ( reco::deltaR2(addTrk->eta(),addTrk->phi(),trigMuon->innerTrack()->eta(),trigMuon->innerTrack()->phi()) < 0.01 )
        continue;

      if ( !dielectronTree->isValid() || !dielectronTree->currentDecayVertex()->vertexIsValid() )
        continue;

      dielectronTree->movePointerToTheTop();
      RefCountedKinematicParticle aPtc = dielectronTree->currentParticle();
      RefCountedKinematicVertex aVtx = dielectronTree->currentDecayVertex();

      const double prob = TMath::Prob(aVtx->chiSquared(),static_cast<int>(std::rint(aVtx->degreesOfFreedom())));

      if ( prob < probThres_ )
        continue;

      const auto dielState = aPtc->currentState();
      dielectronTree->movePointerToTheFirstChild();
      const auto refitEle1 = dielectronTree->currentParticle()->currentState();
      dielectronTree->movePointerToTheNextChild();
      const auto refitEle2 = dielectronTree->currentParticle()->currentState();

      diels.push_back(dielectronFit(dielState,
                                    dielectron(refitEle1,
                                               refitEle2,
                                               aEle.castTo<pat::ElectronRef>(),
                                               *addTrk,
                                               trigMuon),
                                    aVtx->chiSquared(),
                                    aVtx->degreesOfFreedom()));

      if ( dielState.mass() < 2.5 || dielState.mass() > 4. )
        continue;

      Jpsis.push_back(dielectronFit(dielState,
                                    dielectron(refitEle1,
                                               refitEle2,
                                               aEle.castTo<pat::ElectronRef>(),
                                               *addTrk,
                                               trigMuon),
                                    aVtx->chiSquared(),
                                    aVtx->degreesOfFreedom()));
    }

    const auto lvecE1preFit = math::PtEtaPhiMLorentzVector( orgGsfTrk->pt(), orgGsfTrk->eta(), orgGsfTrk->phi(), elmass_ );
    const auto lvecE2preFit = math::PtEtaPhiMLorentzVector( addTrk->pt(), addTrk->eta(), addTrk->phi(), elmass_ );

    bool passMergedElectronID = aEle->electronID("mvaMergedElectron");
    const float mvascore = aEle->userFloat("mvaMergedElectronValues");

    if ( aEle->userInt("mvaMergedElectronCategories")==0 )
      histo1d_["mva_HasTrkEB"]->Fill(mvascore, aWeight);

    auto fill2trkHisto = [&,this] (const std::string& prefix, const dielectronFit& diel) {
      const auto astate = diel.dielectronState;
      const auto amomentum = astate.globalMomentum();
      const double dxyOrgTrk = orgGsfTrk->dxy(primaryVertex->position());
      const double dxyAddTrk = addTrk->dxy(primaryVertex->position());
      const double IPorgTrk = dxyOrgTrk / orgGsfTrk->dxyError(primaryVertex->position(),primaryVertex->covariance());
      const double IPAddTrk = dxyAddTrk / addTrk->dxyError(primaryVertex->position(),primaryVertex->covariance());

      const KinematicState refitPtc1 = diel.dielec.refitFirstEle;
      const auto& vec3P1 = refitPtc1.globalMomentum();
      const auto lvecE1 = math::PtEtaPhiMLorentzVector( vec3P1.perp(), vec3P1.eta(), vec3P1.phi(), elmass_ );
      const KinematicState refitPtc2 = diel.dielec.refitSecondEle;
      const auto& vec3P2 = refitPtc2.globalMomentum();
      const auto lvecE2 = math::PtEtaPhiMLorentzVector( vec3P2.perp(), vec3P2.eta(), vec3P2.phi(), elmass_ );

      histo1d_[prefix+"2trk_diel_invM"]->Fill(astate.mass(),aWeight);
      histo1d_[prefix+"2trk_diel_pt"]->Fill(amomentum.perp(),aWeight);
      histo1d_[prefix+"2trk_diel_pt_preFit"]->Fill( (lvecE1preFit+lvecE2preFit).pt() , aWeight);
      histo1d_[prefix+"2trk_diel_etSC"]->Fill( aEle->superCluster()->energy()*
                                               aEle->superCluster()->position().rho()/
                                               aEle->superCluster()->position().r(), aWeight );
      histo1d_[prefix+"2trk_diel_etU5x5"]->Fill(u5x5Et, aWeight);
      histo1d_[prefix+"2trk_diel_normChi2"]->Fill(diel.chi2/diel.ndof,aWeight);
      histo1d_[prefix+"2trk_diel_rap"]->Fill(math::PtEtaPhiMLorentzVector(amomentum.perp(),amomentum.eta(),amomentum.phi(),astate.mass()).Rapidity(),aWeight);
      histo1d_[prefix+"2trk_diel_dR"]->Fill(reco::deltaR(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi()),aWeight);
      histo1d_[prefix+"2trk_diel_dEta"]->Fill(lvecE2.eta()-lvecE1.eta(),aWeight);
      histo1d_[prefix+"2trk_diel_dPhi"]->Fill( reco::deltaPhi(lvecE2.phi(),lvecE1.phi()), aWeight);
      histo1d_[prefix+"2trk_diel_prob"]->Fill( TMath::Prob(diel.chi2,static_cast<int>(std::rint(diel.ndof))) ,aWeight);

      histo1d_[prefix+"2trk_diel_alphaCalo"]->Fill(alphaCalo,aWeight);
      histo1d_[prefix+"2trk_diel_alphaTrk"]->Fill(alphaTrack,aWeight);
      histo1d_[prefix+"2trk_diel_normDParaIn"]->Fill((*normDParaInHandle)[aEle],aWeight);
      histo1d_[prefix+"2trk_diel_dPerpIn"]->Fill((*dPerpInHandle)[aEle],aWeight);
      histo1d_[prefix+"2trk_diel_u5x5dEtaIn"]->Fill((*union5x5dEtaInHandle)[aEle],aWeight);
      histo1d_[prefix+"2trk_diel_u5x5dPhiIn"]->Fill((*union5x5dPhiInHandle)[aEle],aWeight);
      histo1d_[prefix+"2trk_diel_covIeIe"]->Fill((*union5x5covIeIeHandle)[aEle],aWeight);
      histo1d_[prefix+"2trk_diel_covIeIp"]->Fill((*union5x5covIeIpHandle)[aEle],aWeight);

      histo1d_[prefix+"2trk_diel_trackIso"]->Fill(trackIso,aWeight);
      histo1d_[prefix+"2trk_diel_HoE"]->Fill( aEle->hcalOverEcal(), aWeight );
      histo1d_[prefix+"2trk_diel_sigIeIe"]->Fill( aEle->full5x5_sigmaIetaIeta(), aWeight );
      histo1d_[prefix+"2trk_pt_1st"]->Fill(lvecE1.pt(),aWeight);
      histo1d_[prefix+"2trk_eta_1st"]->Fill(lvecE1.eta(),aWeight);
      histo1d_[prefix+"2trk_dzTrig_1st"]->Fill(orgGsfTrk->vz()-diel.dielec.trigMuon->innerTrack()->vz(),aWeight);
      histo1d_[prefix+"2trk_dxy_1st"]->Fill( dxyOrgTrk, aWeight );
      histo1d_[prefix+"2trk_IPsig_1st"]->Fill( std::abs(IPorgTrk), aWeight );
      histo1d_[prefix+"2trk_pt_2nd"]->Fill(lvecE2.pt(),aWeight);
      histo1d_[prefix+"2trk_eta_2nd"]->Fill(lvecE2.eta(),aWeight);
      histo1d_[prefix+"2trk_dzTrig_2nd"]->Fill(addTrk->vz()-diel.dielec.trigMuon->innerTrack()->vz(),aWeight);
      histo1d_[prefix+"2trk_dxy_2nd"]->Fill( dxyAddTrk, aWeight );
      histo1d_[prefix+"2trk_IPsig_2nd"]->Fill( std::abs(IPAddTrk), aWeight );
      histo1d_[prefix+"2trk_mva_HasTrkEB"]->Fill(mvascore, aWeight);
    };

    auto fill2trkHistos = [&,this] (const std::string& prefix, const dielectronFit& diel) {
      fill2trkHisto(prefix,diel);

      if (passMaskedId)
        fill2trkHisto(prefix+"PassModifiedHeep_",diel);
      else
        fill2trkHisto(prefix+"FailModifiedHeep_",diel);

      if (passMergedElectronID)
        fill2trkHisto(prefix+"PassMergedEle_",diel);
      else
        fill2trkHisto(prefix+"FailMergedEle_",diel);
    };

    if (!diels.empty()) {
      std::sort(diels.begin(),diels.end());

      fill2trkHistos("",diels.front());
    }

    if (!Jpsis.empty()) {
      std::sort(Jpsis.begin(),Jpsis.end());

      fill2trkHistos("Jpsi",Jpsis.front());

      invM_ = Jpsis.front().dielectronState.mass();
      u5x5Et_ = u5x5Et;
      et_ = Jpsis.front().dielec.firstEle->et();
      wgt_ = aWeight;
      passME_ = static_cast<int>(passMergedElectronID);
      passModHeep_ = static_cast<int>(Jpsis.front().dielec.firstEle->electronID("modifiedHeepElectronID"));
      passTrig25_ = static_cast<int>(matchTrigObj(Jpsis.front().dielec.firstEle->p4(),trigObjs_doubleEle25));
      passTrig25unseeded_ = static_cast<int>(matchTrigObj(Jpsis.front().dielec.firstEle->p4(),trigObjs_doubleEle25unseeded));
      passTrig25HE_ = static_cast<int>(matchTrigObj(Jpsis.front().dielec.firstEle->p4(),trigObjs_doubleEle25HE));
      passTrig25caloIdL_ = static_cast<int>(matchTrigObj(Jpsis.front().dielec.firstEle->p4(),trigObjs_doubleEle25caloIdL));
      passTrig33_ = static_cast<int>(matchTrigObj(Jpsis.front().dielec.firstEle->p4(),trigObjs_doubleEle33));
      passTrig33unseeded_ = static_cast<int>(matchTrigObj(Jpsis.front().dielec.firstEle->p4(),trigObjs_doubleEle33unseeded));
      tree_->Fill();

      if (u5x5Et >= 20. && u5x5Et < 25.)
        fill2trkHistos("JpsiPt20to25",Jpsis.front());
      else if (u5x5Et >= 25. && u5x5Et < 30.)
        fill2trkHistos("JpsiPt25to30",Jpsis.front());
      else if (u5x5Et >= 30. && u5x5Et < 35.)
        fill2trkHistos("JpsiPt30to35",Jpsis.front());
      else if (u5x5Et >= 35. && u5x5Et < 40.)
        fill2trkHistos("JpsiPt35to40",Jpsis.front());
      else if (u5x5Et >= 40. && u5x5Et < 50.)
        fill2trkHistos("JpsiPt40to50",Jpsis.front());
      else if (u5x5Et >= 50. && u5x5Et < 60.)
        fill2trkHistos("JpsiPt50to60",Jpsis.front());
      else if (u5x5Et >= 60.)
        fill2trkHistos("JpsiPt60toInf",Jpsis.front());
    }

    std::vector<decaychain> Bdecays;

    for (const auto& trigMuon : trigMuons) {
      // require dZ threshold
      if ( std::abs( orgGsfTrk->vz() - trigMuon->innerTrack()->vz() ) > dzThres_ ||
           std::abs( addTrk->vz() - trigMuon->innerTrack()->vz() ) > dzThres_ )
        continue;

      if ( reco::deltaR2(addTrk->eta(),addTrk->phi(),trigMuon->innerTrack()->eta(),trigMuon->innerTrack()->phi()) < 0.01 )
        continue;

      if (packedPFcandHandle->empty())
        continue;

      for (unsigned iCand = 0; iCand < packedPFcandHandle->size(); iCand++) {
        const auto& aCand = packedPFcandHandle->refAt(iCand).castTo<pat::PackedCandidateRef>();

        if ( std::abs(aCand->pdgId())==13 || std::abs(aCand->pdgId())==11 )
          continue;

        if ( !aCand->hasTrackDetails() || !aCand->trackHighPurity() )
          continue;

        if ( addPackedCand.isNonnull() && addPackedCand==aCand )
          continue;

        const auto* aBestTrack = aCand->bestTrack();

        if ( aBestTrack->pt() < ptThresK_ || std::abs(aBestTrack->eta()) > 2.4 )
          continue;

        if ( std::abs( aBestTrack->vz() - trigMuon->innerTrack()->vz() ) > dzThres_ )
          continue;

        if ( reco::deltaR2(aBestTrack->eta(),aBestTrack->phi(),trigMuon->innerTrack()->eta(),trigMuon->innerTrack()->phi()) < 0.01 )
          continue;

        KinematicParticleFactoryFromTransientTrack kinFactory;
        KinematicParticleVertexFitter kinFitter;
        float elmassErr = elmassErr_; // why this thing requires lvalue? :(
        const reco::TransientTrack firstEle = TTbuilder.build(orgGsfTrk);
        const reco::TransientTrack secondEle = TTbuilder.build(addTrk);
        double chi = 0.;
        double ndf = 0.;

        float kaonMassErr = kaonMassErr_;
        const reco::TransientTrack aTransTrack = TTbuilder.build(aBestTrack);

        std::vector<RefCountedKinematicParticle> Bcandidate;
        Bcandidate.push_back(kinFactory.particle(firstEle,elmass_,chi,ndf,elmassErr));
        Bcandidate.push_back(kinFactory.particle(secondEle,elmass_,chi,ndf,elmassErr));
        Bcandidate.push_back(kinFactory.particle(aTransTrack,kaonMass_,chi,ndf,kaonMassErr));
        RefCountedKinematicTree atree = kinFitter.fit(Bcandidate);

        if ( !atree->isValid() || !atree->currentDecayVertex()->vertexIsValid() )
          continue;

        atree->movePointerToTheTop();
        RefCountedKinematicParticle BmesonPtc = atree->currentParticle();
        RefCountedKinematicVertex BmesonVtx = atree->currentDecayVertex();

        const double prob = TMath::Prob(BmesonVtx->chiSquared(),static_cast<int>(std::rint(BmesonVtx->degreesOfFreedom())));

        if ( prob < probThres_ )
          continue;

        atree->movePointerToTheFirstChild();
        const KinematicState refitPtc1 = atree->currentParticle()->currentState();
        const auto& vec3P1 = refitPtc1.globalMomentum();
        const auto lvecE1 = math::PtEtaPhiMLorentzVector( vec3P1.perp(), vec3P1.eta(), vec3P1.phi(), elmass_ );
        atree->movePointerToTheNextChild();
        const KinematicState refitPtc2 = atree->currentParticle()->currentState();
        const auto& vec3P2 = refitPtc2.globalMomentum();
        const auto lvecE2 = math::PtEtaPhiMLorentzVector( vec3P2.perp(), vec3P2.eta(), vec3P2.phi(), elmass_ );
        atree->movePointerToTheNextChild();
        const KinematicState refitPtc3 = atree->currentParticle()->currentState();
        const auto& vec3P3 = refitPtc3.globalMomentum();
        const auto lvecK = math::PtEtaPhiMLorentzVector( vec3P3.perp(), vec3P3.eta(), vec3P3.phi(), kaonMass_ );

        const reco::Vertex aVtx = *BmesonVtx;
        const auto lvecEE = lvecE1+lvecE2;
        std::pair<bool, Measurement1D> ip3dKaon = IPTools::signedImpactParameter3D(aTransTrack, GlobalVector(lvecEE.x(),lvecEE.y(),lvecEE.z()), aVtx);
        const double d0Kaon = ip3dKaon.second.value();

        if ( std::abs(d0Kaon) > d0Thres_ )
          continue;

        const double dvx = BmesonVtx->position().x() - beamSpotHandle->position().x();
        const double dvy = BmesonVtx->position().y() - beamSpotHandle->position().y();
        const double px = BmesonPtc->currentState().globalMomentum().x();
        const double py = BmesonPtc->currentState().globalMomentum().y();
        const double cosAlpha2d = (px*dvx+py*dvy) / ( std::sqrt(dvx*dvx+dvy*dvy)*std::sqrt(px*px+py*py) );

        if ( cosAlpha2d < cosAlpha2dThres_ )
          continue;

        if ( lvecEE.M() < 2.5 || lvecEE.M() > 4. )
          continue;

        if ( reco::deltaR2(lvecEE.eta(),lvecEE.phi(),lvecK.eta(),lvecK.phi()) > 0.64 )
          continue;

        if ( (lvecE1+lvecK).M() < 2. || (lvecE2+lvecK).M() < 2. )
          continue;

        const double u5x5En = (*union5x5EnergyHandle)[aEle];
        const auto lvecJpsi = math::PtEtaPhiELorentzVector( std::sqrt(u5x5En*u5x5En - jpsiMass_*jpsiMass_)/std::cosh(lvecEE.eta()),
                                                            lvecEE.eta(),
                                                            lvecEE.phi(),
                                                            u5x5En );
        const auto lvecBu = lvecJpsi + lvecK;

        if ( lvecBu.M() < 4.5 || lvecBu.M() > 6. )
          continue;

        const auto aDielectron = dielectron(refitPtc1,
                                            refitPtc2,
                                            aEle.castTo<pat::ElectronRef>(),
                                            *addTrk,
                                            trigMuon);

        auto achain = decaychain( aDielectron,
                                  refitPtc3,
                                  aCand,
                                  BmesonPtc->currentState(),
                                  BmesonVtx->chiSquared(),
                                  BmesonVtx->degreesOfFreedom(),
                                  d0Kaon,
                                  cosAlpha2d );

        Bdecays.push_back( std::move(achain) );
      } // aCand
    }

    std::sort(Bdecays.begin(),Bdecays.end());

    auto estimateJpsiVec = [&,this] (const decaychain& achain) {
      const auto lvecE1 = math::PtEtaPhiMLorentzVector( achain.dielec.refitFirstEle.globalMomentum().perp(),
                                                        achain.dielec.refitFirstEle.globalMomentum().eta(),
                                                        achain.dielec.refitFirstEle.globalMomentum().phi(), elmass_ );
      const auto lvecE2 = math::PtEtaPhiMLorentzVector( achain.dielec.refitSecondEle.globalMomentum().perp(),
                                                        achain.dielec.refitSecondEle.globalMomentum().eta(),
                                                        achain.dielec.refitSecondEle.globalMomentum().phi(), elmass_ );
      const double u5x5En = (*union5x5EnergyHandle)[achain.dielec.firstEle];
      const auto lvecEE = lvecE1+lvecE2;
      const auto lvecJpsi = math::PtEtaPhiELorentzVector( std::sqrt(u5x5En*u5x5En - jpsiMass_*jpsiMass_)/std::cosh(lvecEE.eta()),
                                                          lvecEE.eta(),
                                                          lvecEE.phi(),
                                                          u5x5En );

      return lvecJpsi;
    };

    auto BmesonMassU5x5 = [&,this] (const decaychain& achain) -> double {
      const auto lvecK = math::PtEtaPhiMLorentzVector( achain.refitThirdTrk.globalMomentum().perp(),
                                                       achain.refitThirdTrk.globalMomentum().eta(),
                                                       achain.refitThirdTrk.globalMomentum().phi(), kaonMass_ );
      const auto lvecJpsi = estimateJpsiVec(achain);

      return (lvecJpsi+lvecK).M();
    };

    auto fill3trkHisto = [&,this] (const std::string& prefix, const decaychain& achain) {
      const auto lvecP1 = math::PtEtaPhiMLorentzVector( achain.dielec.refitFirstEle.globalMomentum().perp(),
                                                        achain.dielec.refitFirstEle.globalMomentum().eta(),
                                                        achain.dielec.refitFirstEle.globalMomentum().phi(), pionMass_ );
      const auto lvecP2 = math::PtEtaPhiMLorentzVector( achain.dielec.refitSecondEle.globalMomentum().perp(),
                                                        achain.dielec.refitSecondEle.globalMomentum().eta(),
                                                        achain.dielec.refitSecondEle.globalMomentum().phi(), pionMass_ );
      const auto lvecK1 = math::PtEtaPhiMLorentzVector( achain.dielec.refitFirstEle.globalMomentum().perp(),
                                                        achain.dielec.refitFirstEle.globalMomentum().eta(),
                                                        achain.dielec.refitFirstEle.globalMomentum().phi(), kaonMass_ );
      const auto lvecK2 = math::PtEtaPhiMLorentzVector( achain.dielec.refitSecondEle.globalMomentum().perp(),
                                                        achain.dielec.refitSecondEle.globalMomentum().eta(),
                                                        achain.dielec.refitSecondEle.globalMomentum().phi(), kaonMass_ );
      const auto lvecE1 = math::PtEtaPhiMLorentzVector( achain.dielec.refitFirstEle.globalMomentum().perp(),
                                                        achain.dielec.refitFirstEle.globalMomentum().eta(),
                                                        achain.dielec.refitFirstEle.globalMomentum().phi(), elmass_ );
      const auto lvecE2 = math::PtEtaPhiMLorentzVector( achain.dielec.refitSecondEle.globalMomentum().perp(),
                                                        achain.dielec.refitSecondEle.globalMomentum().eta(),
                                                        achain.dielec.refitSecondEle.globalMomentum().phi(), elmass_ );
      const auto lvecK = math::PtEtaPhiMLorentzVector( achain.refitThirdTrk.globalMomentum().perp(),
                                                       achain.refitThirdTrk.globalMomentum().eta(),
                                                       achain.refitThirdTrk.globalMomentum().phi(), kaonMass_ );
      const auto lvecEE = lvecE1+lvecE2;
      const auto bMomentum = achain.Bmeson.globalMomentum();

      histo1d_[prefix+"3trk_Pi1Pi2_invM"]->Fill((lvecP1+lvecP2).M(),aWeight);
      histo1d_[prefix+"3trk_Pi1K2_invM"]->Fill((lvecP1+lvecK2).M(),aWeight);
      histo1d_[prefix+"3trk_K1Pi2_invM"]->Fill((lvecK1+lvecP2).M(),aWeight);

      histo1d_[prefix+"3trk_e1K_invM"]->Fill((lvecE1+lvecK).M(),aWeight);
      histo1d_[prefix+"3trk_e2K_invM"]->Fill((lvecE2+lvecK).M(),aWeight);

      histo1d_[prefix+"3trk_diel_invM"]->Fill(lvecEE.M(),aWeight);
      histo1d_[prefix+"3trk_diel_invM_wide"]->Fill(lvecEE.M(),aWeight);
      histo1d_[prefix+"3trk_diel_pt"]->Fill(lvecEE.pt(),aWeight);
      histo1d_[prefix+"3trk_diel_etSC"]->Fill( achain.dielec.firstEle->superCluster()->energy()*
                                               achain.dielec.firstEle->superCluster()->position().rho()/
                                               achain.dielec.firstEle->superCluster()->position().r(), aWeight );
      histo1d_[prefix+"3trk_diel_etU5x5"]->Fill(u5x5Et,aWeight);
      histo1d_[prefix+"3trk_diel_rap"]->Fill(lvecEE.Rapidity(), aWeight);
      histo1d_[prefix+"3trk_diel_dR"]->Fill(reco::deltaR(lvecE1.eta(),lvecE1.phi(),lvecE2.eta(),lvecE2.phi()) , aWeight);
      histo1d_[prefix+"3trk_diel_dEta"]->Fill( lvecE2.eta() - lvecE1.eta(), aWeight );
      histo1d_[prefix+"3trk_diel_dPhi"]->Fill( reco::deltaPhi(lvecE2.phi(),lvecE1.phi()), aWeight );

      histo1d_[prefix+"3trk_diel_alphaCalo"]->Fill(alphaCalo,aWeight);
      histo1d_[prefix+"3trk_diel_alphaTrk"]->Fill(alphaTrack,aWeight);
      histo1d_[prefix+"3trk_diel_dPerpIn"]->Fill((*normDParaInHandle)[aEle],aWeight);
      histo1d_[prefix+"3trk_diel_normDParaIn"]->Fill((*normDParaInHandle)[aEle],aWeight);
      histo1d_[prefix+"3trk_diel_u5x5dEtaIn"]->Fill((*union5x5dEtaInHandle)[aEle],aWeight);
      histo1d_[prefix+"3trk_diel_u5x5dPhiIn"]->Fill((*union5x5dPhiInHandle)[aEle],aWeight);
      histo1d_[prefix+"3trk_diel_covIeIe"]->Fill((*union5x5covIeIeHandle)[aEle],aWeight);
      histo1d_[prefix+"3trk_diel_covIeIp"]->Fill((*union5x5covIeIpHandle)[aEle],aWeight);

      histo1d_[prefix+"3trk_diel_HoE"]->Fill( achain.dielec.firstEle->hcalOverEcal(), aWeight );
      histo1d_[prefix+"3trk_diel_sigIeIe"]->Fill( achain.dielec.firstEle->full5x5_sigmaIetaIeta(), aWeight );
      histo1d_[prefix+"3trk_pt_1st"]->Fill(achain.dielec.refitFirstEle.globalMomentum().perp(), aWeight);
      histo1d_[prefix+"3trk_eta_1st"]->Fill(achain.dielec.refitFirstEle.globalMomentum().eta(), aWeight);
      histo1d_[prefix+"3trk_pt_2nd"]->Fill(achain.dielec.refitSecondEle.globalMomentum().perp(), aWeight);
      histo1d_[prefix+"3trk_eta_2nd"]->Fill(achain.dielec.refitSecondEle.globalMomentum().eta(), aWeight);
      histo1d_[prefix+"3trk_dzTrig_1st"]->Fill( achain.dielec.firstEle->gsfTrack()->vz() - achain.dielec.trigMuon->vz(), aWeight );
      histo1d_[prefix+"3trk_IPsig_1st"]->Fill( achain.dielec.firstEle->gsfTrack()->dxy(primaryVertex->position()) /
                                               achain.dielec.firstEle->gsfTrack()->dxyError(primaryVertex->position(),primaryVertex->covariance()) , aWeight );
      histo1d_[prefix+"3trk_dzTrig_2nd"]->Fill( achain.dielec.secondEle.vz() - achain.dielec.trigMuon->vz(), aWeight );
      histo1d_[prefix+"3trk_IPsig_2nd"]->Fill( achain.dielec.secondEle.dxy(primaryVertex->position()) /
                                               achain.dielec.secondEle.dxyError(primaryVertex->position(),primaryVertex->covariance()) , aWeight );

      histo1d_[prefix+"3trk_pt_3rd"]->Fill( achain.refitThirdTrk.globalMomentum().perp(), aWeight );
      histo1d_[prefix+"3trk_eta_3rd"]->Fill( achain.refitThirdTrk.globalMomentum().eta(), aWeight );
      histo1d_[prefix+"3trk_drAtSC_3rd"]->Fill( extrapolateToSC(*(achain.dielec.firstEle),*(achain.thirdTrk->bestTrack()),*beamSpotHandle,iSetup) ,aWeight);
      histo1d_[prefix+"3trk_dzTrig_3rd"]->Fill( achain.thirdTrk->bestTrack()->vz() - achain.dielec.trigMuon->vz(), aWeight );
      histo1d_[prefix+"3trk_IPsig_3rd"]->Fill( achain.thirdTrk->bestTrack()->dxy(primaryVertex->position()) /
                                               achain.thirdTrk->bestTrack()->dxyError(primaryVertex->position(),primaryVertex->covariance()) , aWeight );
      histo1d_[prefix+"3trk_d0_3rd"]->Fill(achain.d0thirdTrk, aWeight);

      histo1d_[prefix+"3trk_B_invM"]->Fill(achain.Bmeson.mass(),aWeight);
      histo1d_[prefix+"3trk_B_invM_wide"]->Fill(achain.Bmeson.mass(),aWeight);
      histo1d_[prefix+"3trk_B_invM_u5x5"]->Fill(BmesonMassU5x5(achain),aWeight);
      histo1d_[prefix+"3trk_B_pt"]->Fill(bMomentum.perp(),aWeight);
      histo1d_[prefix+"3trk_B_normChi2"]->Fill(achain.BmesonChi2/achain.BmesonNdof,aWeight);
      histo1d_[prefix+"3trk_B_prob"]->Fill( TMath::Prob(achain.BmesonChi2,static_cast<int>(std::rint(achain.BmesonNdof))) , aWeight );
      histo1d_[prefix+"3trk_B_rap"]->Fill(math::PtEtaPhiMLorentzVector(bMomentum.perp(),bMomentum.eta(),bMomentum.phi(),achain.Bmeson.mass()).Rapidity(),aWeight);
      histo1d_[prefix+"3trk_B_dR"]->Fill( reco::deltaR(lvecK.eta(),lvecK.phi(),lvecEE.eta(),lvecEE.phi()) ,aWeight);
      histo1d_[prefix+"3trk_B_dR_wide"]->Fill( reco::deltaR(lvecK.eta(),lvecK.phi(),lvecEE.eta(),lvecEE.phi()) ,aWeight);
      histo1d_[prefix+"3trk_B_cosAlpha2d"]->Fill( achain.cosAlpha2d, aWeight );

      const double massE1X = (lvecE1+lvecK).M();
      const double massE2X = (lvecE2+lvecK).M();
      const double massE1E2 = lvecEE.M();

      histo2d_[prefix+"3trk_dalitz_E1X_E2X"]->Fill( massE1X*massE1X, massE2X*massE2X ,aWeight );
      histo2d_[prefix+"3trk_dalitz_E1E2_E2X"]->Fill( massE1E2*massE1E2, massE2X*massE2X, aWeight );
      histo2d_[prefix+"3trk_dalitz_E1E2_E1X"]->Fill( massE1E2*massE1E2, massE1X*massE1X, aWeight );
    };

    if (!Bdecays.empty()) {
      auto fill3trkHistos = [&,this] (const std::string& prefix) {
        fill3trkHisto(prefix,Bdecays.front());

        if (passMaskedId)
          fill3trkHisto(prefix+"PassModifiedHeep_",Bdecays.front());
        else
          fill3trkHisto(prefix+"FailModifiedHeep_",Bdecays.front());

        if (passMergedElectronID)
          fill3trkHisto(prefix+"PassMergedEle_",Bdecays.front());
        else
          fill3trkHisto(prefix+"FailMergedEle_",Bdecays.front());

        if (passMaskedId && passMergedElectronID)
          fill3trkHisto(prefix+"PassModifiedHeepPassMergedEle_",Bdecays.front());
      };

      const auto lvecJpsi = estimateJpsiVec(Bdecays.front());

      BinvMtrk_ = Bdecays.front().Bmeson.mass();
      BinvMcalo_ = BmesonMassU5x5(Bdecays.front());
      Bu5x5Et_ = u5x5Et;
      Bwgt_ = aWeight;
      BpassME_ = static_cast<int>(passMergedElectronID);
      BpassModHeep_ = static_cast<int>(passMaskedId);
      ptK_ = Bdecays.front().refitThirdTrk.globalMomentum().perp();
      etaK_ = Bdecays.front().refitThirdTrk.globalMomentum().eta();
      phiK_ = Bdecays.front().refitThirdTrk.globalMomentum().phi();
      u5x5En_ = (*union5x5EnergyHandle)[aEle];
      etaJpsi_ = lvecJpsi.eta();
      phiJpsi_ = lvecJpsi.phi();

      Btree_->Fill();

      fill3trkHistos("");
    } // if Bdecays
  } // slimmedElectrons
}

DEFINE_FWK_MODULE(MergedLeptonIDJpsiAnalyzer);
