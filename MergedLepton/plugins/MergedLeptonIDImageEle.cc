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

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedDEtaInSeed.h"

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
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"

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
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/PreshowerCluster.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"


#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

// produce TTree for merged electron training with H->AA->4e events

class MergedLeptonIDImageEle : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonIDImageEle(const edm::ParameterSet&);
  virtual ~MergedLeptonIDImageEle() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  edm::ConsumesCollector iC = consumesCollector();

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
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
  const edm::EDGetTokenT<edm::View<reco::GsfElectron>> emObjectToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweight_token;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

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

  TTree* ImageTree_ = nullptr;
  std::vector<std::vector<float>> EEImage_branch; 
  std::vector<std::vector<float>> ES1Image_branch; 
  std::vector<std::vector<float>> ES2Image_branch; 
  int label_num_ele;
  int label_num_ele_hard;
  float pT_gsfele;
  int isAddTrk;
  std::vector<int> EleXY;

  int imageSize_;
  int ESimageSize_;

  PositionCalc posCalcLog_;

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

MergedLeptonIDImageEle::MergedLeptonIDImageEle(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
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
emObjectToken_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("emObjectProducer"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
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
geometryToken_(esConsumes()),
imageSize_(iConfig.getParameter<int>("imageSize")),
ESimageSize_(iConfig.getParameter<int>("ESimageSize"))
{
  std::cout<<"hello"<<std::endl;
  std::cout<<"hello2"<<std::endl;
  
  usesResource("TFileService");
}

void MergedLeptonIDImageEle::beginJob() {
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

  ImageTree_ = fs->make<TTree>("ImageTree","ImageTree");
  ImageTree_->Branch("EEImage",&EEImage_branch,32000,0);
  ImageTree_->Branch("ES1Image",&ES1Image_branch,32000,0);
  ImageTree_->Branch("ES2Image",&ES2Image_branch,32000,0);
  ImageTree_->Branch("NumEle",&label_num_ele,"NumEle/I");
  ImageTree_->Branch("NumEleHard",&label_num_ele_hard,"NumEleHard/I");
  ImageTree_->Branch("pT",&pT_gsfele,"pT/F");
  ImageTree_->Branch("IsAddTrk",&isAddTrk,"IsAddTrk/I");
  ImageTree_->Branch("EleXY",&EleXY,32000,0);
}

void MergedLeptonIDImageEle::endJob() {
  purwgtFile_->Close();
}


void MergedLeptonIDImageEle::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  edm::Handle< edm::View<reco::GsfElectron> > emObjectHandle;
  iEvent.getByToken(emObjectToken_,emObjectHandle);

  edm::Handle<EcalRecHitCollection> ESrecHitHandle;
  iEvent.getByToken(ESrecHitToken_, ESrecHitHandle);
  
  edm::Handle<EcalRecHitCollection> EErecHitHandle;
  iEvent.getByToken(EErecHitToken_, EErecHitHandle);
  
  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(genptcToken_, genptcHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  std::vector<reco::GenParticleRef> promptEles;
  std::vector<reco::GenParticleRef> Eles;
  for (unsigned int idx =0; idx<genptcHandle->size();++idx){
	  const auto& genPtc = genptcHandle->refAt(idx);
	  //std::cout <<"flag : "<< genPtc->statusFlags().flags_<<" | status : "<<genPtc->status()<<" | pdg : " <<genPtc->pdgId()<<" | hard "<<genPtc->isHardProcess()  <<std::endl;    
	  if ( ( std::abs(genPtc->pdgId())==11 ) && genPtc->fromHardProcessFinalState() ) promptEles.push_back(genPtc.castTo<reco::GenParticleRef>());
	  if ( ( std::abs(genPtc->pdgId())==11 ) && genPtc->isPromptFinalState() ) Eles.push_back(genPtc.castTo<reco::GenParticleRef>());
  }

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
  int idx =0;
  //std::cout<<"wtf"<<std::endl;

  //auto dEtaInSeedCalc = ModifiedDEtaInSeed(posCalcLog_,iC);
  for (const auto& electron : *emObjectHandle){
    if (electron.eta()<1.653 || electron.eta()>2.5){
      continue;
    }
    const auto& aEle = emObjectHandle->refAt(idx);
    idx ++;
    const auto& orgGsfTrk = electron.gsfTrack();
    const auto& addGsfTrk = (*addGsfTrkMap)[aEle];
    const auto& addPackedCand = (*addPackedCandHandle)[aEle];

    const reco::Track* addTrk = addGsfTrk.get();

    if ( addGsfTrk==orgGsfTrk && addPackedCand.isNonnull() )
      addTrk = addPackedCand->bestTrack();
    if (orgGsfTrk.get()->pt() < 20 || addTrk->pt() < 10) continue;
    if(addTrk == orgGsfTrk.get()) isAddTrk = 0;
    else isAddTrk = 1;

    //const bool isPackedCand = addGsfTrk==orgGsfTrk && addPackedCand.isNonnull();
    //const reco::TrackBase* addTrk = isPackedCand ? addPackedCand->bestTrack() : addGsfTrk.get();
    //std::cout<<"hello"<<std::endl;
    const auto& seedPosition = electron.superCluster()->seed()->position(); // supercluster seed position
    
    //calculate 1st and 2nd ele position
    //
    //float eta_1st = seedPosition.eta() - electron.deltaEtaSeedClusterTrackAtCalo();
    //float phi_1st = reco::reduceRange( seedPosition.phi() - electron.deltaPhiSeedClusterTrackAtCalo());
    
    //auto beamSpot = beamSpotHandle.product();

    //double dEtaInSeed2nd = std::numeric_limits<float>::max();
    //double dPhiInSC2nd = std::numeric_limits<float>::max();

    //auto scAtVtx = EleRelPointPair(math::XYZPoint(),math::XYZPoint(),beamSpot->position());
    //auto seedAtCalo = EleRelPointPair(math::XYZPoint(),math::XYZPoint(),beamSpot->position());

    //if ( dEtaInSeedCalc.extrapolate(electron,*addTrk,beamSpot->position(),iSetup,scAtVtx,seedAtCalo) ) {
      //dPhiInSC2nd = scAtVtx.dPhi();
      //dEtaInSeed2nd = scAtVtx.dEta() - electron.superCluster()->eta() + electron.superCluster()->seed()->eta();
    //}

    //if ( dEtaInSeed2nd==std::numeric_limits<float>::max() || dPhiInSC2nd==std::numeric_limits<float>::max() )
      //continue;

    //float eta_2nd = -( dEtaInSeed2nd - electron.superCluster()->seed()->eta() );
    //float phi_2nd = reco::reduceRange( -( dPhiInSC2nd - electron.superCluster()->phi() ) );

    // eta phi to x y
   // float radius_1st = 295./sinh(eta_1st);
   // float x_1st = radius_1st * cos(phi_1st);
   // float y_1st = radius_1st * sin(phi_1st);
   // 
   // float radius_2nd = 295./sinh(eta_2nd);
   // float x_2nd = radius_2nd * cos(phi_2nd);
   // float y_2nd = radius_2nd * sin(phi_2nd);
    std::vector<float> genelexy;
    //std::cout<<"pT : "<<electron.pt()<<std::endl;
    //
    //
    //std::cout<<"5x5 energy : "<<electron.e5x5()<<std::endl;
    //std::cout<<electron.superCluster()->seed()->eta()<<" | "<<seedPosition.eta()<<std::endl;
    //std::cout<<orgGsfTrk.get()->pt()<<" | add : "<<addTrk->pt()<<" | "<<isAddTrk<<std::endl;
    int halfSize = imageSize_ / 2;
    int EShalfSize = ESimageSize_ / 2;
    std::vector<std::vector<float>> EEImage(imageSize_,std::vector<float>(imageSize_,0.0));
    std::vector<std::vector<float>> ESImage_plane1(ESimageSize_*32,std::vector<float>(ESimageSize_,0.0));
    std::vector<std::vector<float>> ESImage_plane2(ESimageSize_,std::vector<float>(ESimageSize_*32,0.0));
    const EEDetId* matchedCrystal = nullptr;
    const ESDetId* matchedSilicon = nullptr;
    int matched_ix=electron.superCluster()->seedCrysIEtaOrIx();
    int matched_iy=electron.superCluster()->seedCrysIPhiOrIy();
    float minDistance = std::numeric_limits<float>::max();
    int matched_gen_ele = 0;
    int matched_gen_prompt_ele = 0;
    pT_gsfele = electron.pt();
    std::vector<float> geneleetaphi; 
    for (const auto& ele : promptEles){
       float dEta = ele->eta()-seedPosition.eta();
       float dPhi = ele->phi()-seedPosition.phi();
       float dR = std::sqrt(dEta * dEta + dPhi * dPhi);
       if (dR < 0.1) {
	  matched_gen_prompt_ele ++;
	  float radius = 306./sinh(ele->eta());
	  genelexy.push_back(radius * cos(ele->phi()));
	  genelexy.push_back(radius * sin(ele->phi()));
	  geneleetaphi.push_back(ele->eta());
	  geneleetaphi.push_back(ele->phi());
       }
    }
    if (matched_gen_prompt_ele < 2) continue;
    for (const auto& ele : Eles){
       float dEta = ele->eta()-seedPosition.eta();
       float dPhi = ele->phi()-seedPosition.phi();
       float dR = std::sqrt(dEta * dEta + dPhi * dPhi);
       if (dR < 0.1) matched_gen_ele ++;
    }
    //for (const auto& hit : *EErecHitHandle){
    //  const auto& detID = hit.id();
    //  EEDetId id_xtal(hit.detid());
    //  const auto& hitPosition = caloGeom->getGeometry(detID);
    //  if (hitPosition->getPosition().z() * seedPosition.z() < 0) continue;
    //  float deltaX = hitPosition->getPosition().x()-seedPosition.x();
    //  float deltaY = hitPosition->getPosition().y()-seedPosition.y();

    //  float distance = std::sqrt(deltaX * deltaX + deltaY * deltaY);
    //  if (distance < minDistance){
    //    minDistance = distance;
    //    matchedCrystal = &id_xtal;
    //    matched_ix = id_xtal.ix();
    //    matched_iy = id_xtal.iy();
    //  }
    //}
    //if (matchedCrystal){
    for (const auto& hit : *EErecHitHandle){
      const auto& detID = hit.id();
      auto id_xtal =EEDetId(hit.detid());
      const auto& hitPosition = caloGeom->getGeometry(detID);
    	if (hitPosition->getPosition().z() * seedPosition.z() < 0) continue;

      int dX = matched_ix-id_xtal.ix();
      int dY = matched_iy-id_xtal.iy();
      //std::cout<<dX<<" | "<<abs(dX)<<" | "<<id_xtal.ix()<<" | "<<matchedCrystal->ix()<<std::endl;
      if (abs(dX) <= halfSize && abs(dY) <= halfSize){
        int iX = dX + halfSize;
        int iY = dY + halfSize;
        EEImage[iX][iY] = hit.energy();
      }
    }
    //}
   // for (const auto& row : EEImage) {
   //   for (const auto& pixel : row) {
   //     std::cout << pixel << " ";
   //   }
   //   std::cout << std::endl;
   // }
    //matched_ix=ESDetId(electron.superCluster()->preshowerClusters().seed()).six();
    //matched_iy=ESDetId(electron.superCluster()->preshowerClusters().seed()).siy();
    matched_ix=std::numeric_limits<int>::min();
    matched_iy=std::numeric_limits<int>::min();
    minDistance = std::numeric_limits<float>::max();
    float mindX_1st = std::numeric_limits<float>::max();
    float mindY_1st = std::numeric_limits<float>::max();
    float mindX_2nd = std::numeric_limits<float>::max();
    float mindY_2nd = std::numeric_limits<float>::max();
    int six_1st= std::numeric_limits<int>::min();
    int siy_1st= std::numeric_limits<int>::min();
    int six_2nd= std::numeric_limits<int>::min();
    int siy_2nd= std::numeric_limits<int>::min();
    int stripx_1st = std::numeric_limits<int>::min();
    int stripy_1st = std::numeric_limits<int>::min();
    int stripx_2nd = std::numeric_limits<int>::min();
    int stripy_2nd = std::numeric_limits<int>::min();
    
    for (const auto& hit : *ESrecHitHandle){
      const auto& detID = hit.id();
      ESDetId id_silicon(hit.detid());
      const auto& hitPosition = caloGeom->getGeometry(detID);
      if (hitPosition->getPosition().z() * seedPosition.z() < 0) continue;

      float deltaEta = hitPosition->getPosition().eta()-electron.superCluster()->seed()->eta();
      float deltaPhi = hitPosition->getPosition().phi()-electron.superCluster()->seed()->phi();

      if(id_silicon.plane() == 1){
	float dX_1st_iter = pow(hitPosition->getPosition().x() - genelexy.at(0),2);
        float dX_2nd_iter = pow(hitPosition->getPosition().x() - genelexy.at(2),2);
	if (dX_1st_iter < mindX_1st ){
		stripx_1st = id_silicon.strip();
		six_1st = id_silicon.six();
		mindX_1st = dX_1st_iter;
	}
	if (dX_2nd_iter < mindX_2nd ){
		stripx_2nd = id_silicon.strip();
		six_2nd = id_silicon.six();
		mindX_2nd = dX_2nd_iter;
	}
      
      }
      if(id_silicon.plane() == 2){
	float dY_1st_iter = pow(hitPosition->getPosition().y() - genelexy.at(1),2);
        float dY_2nd_iter = pow(hitPosition->getPosition().y() - genelexy.at(3),2);
	if (dY_1st_iter < mindY_1st ){
		stripy_1st = id_silicon.strip();
		siy_1st = id_silicon.siy();
		mindY_1st = dY_1st_iter;
	}
	if (dY_2nd_iter < mindY_2nd ){
		stripy_2nd = id_silicon.strip();
		siy_2nd = id_silicon.siy();
		mindY_2nd = dY_2nd_iter;
	}
      
      }
      float distance = std::sqrt(deltaEta * deltaEta + deltaPhi * deltaPhi);
      if (distance < minDistance){
        minDistance = distance;
        matchedSilicon = &id_silicon;
        matched_ix = id_silicon.six();
        matched_iy = id_silicon.siy();
      }
      
    }
    //std::cout<<matched_ix <<" | "<<matched_iy<<std::endl;
    if (matchedSilicon){
      
      for (const auto& hit : *ESrecHitHandle){
        const auto& detID = hit.id();
        ESDetId id_silicon(hit.detid());
        const auto& hitPosition = caloGeom->getGeometry(detID);
	//std::cout<<hitPosition->getPosition().z()<<std::endl;
        if (hitPosition->getPosition().z() * seedPosition.z() < 0) continue;
        int dX = matched_ix-id_silicon.six();
        int dY = matched_iy-id_silicon.siy();
        if (abs(dX) <= EShalfSize && abs(dY) <= EShalfSize){
          int iX = dX + EShalfSize;
          int iY = dY + EShalfSize;
          if (id_silicon.plane() == 1)ESImage_plane1[iX*32+(int)id_silicon.strip()-1][iY] = hit.energy();
	  if (id_silicon.plane() == 2)ESImage_plane2[iX][iY*32+id_silicon.strip()-1] = hit.energy();
          //if (id_silicon.plane() == 1)ESImage_plane1[iX*32+(int)id_silicon.strip()-1][iY] = id_silicon.siy();
          //if (id_silicon.plane() == 2)ESImage_plane2[iX][iY*32+id_silicon.strip()-1] = id_silicon.siy();
        }
      
      
      }
    }
    float test_radius = seedPosition.z()/sinh(seedPosition.eta());
    //std::cout<< seedPosition.z()<<" | "<<test_radius<<" | "<< test_radius*sin(seedPosition.phi())<<" | "<<test_radius*cos(seedPosition.phi())<<std::endl;
    //std::cout << geneleetaphi.at(0) << " | " <<geneleetaphi.at(1)<<" | | "<<seedPosition.y()<<" | "<<seedPosition.x()<<std::endl; 
    //for (const auto& row : ESImage_plane1) {
    //  for (const auto& pixel : row) {
    //    std::cout << pixel << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << "----" << std::endl;
    //for (const auto& row : ESImage_plane2) {
    //  for (const auto& pixel : row) {
    //    std::cout << pixel << " ";
    //  }
    //  std::cout << std::endl;
    //}
    //std::cout << "----" << std::endl;
    EleXY.push_back((matched_ix - six_1st + EShalfSize) * 32 + stripx_1st - 1);
    EleXY.push_back((matched_iy - siy_1st + EShalfSize) * 32 + stripy_1st - 1);
    EleXY.push_back((matched_ix - six_2nd + EShalfSize) * 32 + stripx_2nd - 1);
    EleXY.push_back((matched_iy - siy_2nd + EShalfSize) * 32 + stripy_2nd - 1);
    EEImage_branch = EEImage;
    ES1Image_branch = ESImage_plane1;
    ES2Image_branch = ESImage_plane2;
    //std::cout << matched_gen_prompt_ele <<" | "<<matched_gen_ele<<std::endl;
    //std::cout << "----" << std::endl;
    label_num_ele_hard = matched_gen_prompt_ele;
    label_num_ele = matched_gen_ele;
    ImageTree_->Fill();
    EleXY.clear();
    //Ele2nd_XY.clear();
  }



  return;
}

DEFINE_FWK_MODULE(MergedLeptonIDImageEle);
