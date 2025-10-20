#include <memory>
#include <iostream>


#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedDEtaInSeed.h"


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
#include "DataFormats/MuonReco/interface/Muon.h"
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

#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "TrackingTools/GsfTools/interface/GsfPropagatorAdapter.h"
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
#include "DataFormats/MuonReco/interface/MuonEnergy.h"


#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

// produce TTree for merged electron training with H->AA->4e events

class MergedMuon : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedMuon(const edm::ParameterSet&);
  virtual ~MergedMuon() {}
  bool extrapolate(const reco::GsfElectron& aEle, const reco::TrackBase& addTrk,
                   const math::XYZPoint& beamSpotPos, const edm::EventSetup& iSetup,
                   EleRelPointPair& scAtVtx, EleRelPointPair& seedAtCalo);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  const edm::EDGetTokenT<edm::View<reco::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaTrackToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaCaloToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> normDParaInToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedPFcandToken_;
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
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> topologyToken_;

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

  TTree* muon_ = nullptr;

  float pT_muon;
  float eta_muon;
  float phi_muon;
  float em_muon;
  float emS9_muon;
  float emS25_muon;
  float emMax_muon;
  float had_muon;
  float hadS9_muon;
  float hadMax_muon;
  float ho_muon;
  float hoS9_muon;

  TTree* mergedMuon_ = nullptr;

  float pT_mergedMuon;
  float eta_mergedMuon;
  float phi_mergedMuon;
  float em_mergedMuon;
  float emS9_mergedMuon;
  float emS25_mergedMuon;
  float emMax_mergedMuon;
  float had_mergedMuon;
  float hadS9_mergedMuon;
  float hadMax_mergedMuon;
  float ho_mergedMuon;
  float hoS9_mergedMuon;

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

MergedMuon::MergedMuon(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
alphaTrackToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaTrack"))),
alphaCaloToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaCalo"))),
normDParaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("normalizedDParaIn"))),
packedPFcandToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedPFcand"))),
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
topologyToken_(esConsumes())
{
  std::cout<<"hello"<<std::endl;
  std::cout<<"hello2"<<std::endl;
  
  usesResource("TFileService");
}

bool MergedMuon::extrapolate(const reco::GsfElectron& aEle, const reco::TrackBase& addTrk,
                                     const math::XYZPoint& beamSpotPos, const edm::EventSetup& iSetup,
                                     EleRelPointPair& scAtVtx, EleRelPointPair& seedAtCalo) {
  // track-cluster matching (see RecoEgamma/EgammaElectronAlgos/src/GsfElectronAlgo.cc)
  //edm::ESHandle<MagneticField> magFieldHandle;
  //iSetup.get<IdealMagneticFieldRecord>().get(magFieldHandle);
  const MagneticField* magFieldHandle = &iSetup.getData(magneticToken_);
  // at innermost/outermost point
  // edm::ESHandle<TrackerGeometry> trackerHandle;
  // iSetup.get<TrackerDigiGeometryRecord>().get(trackerHandle);
  // auto mtsTransform = std::make_unique<MultiTrajectoryStateTransform>(trackerHandle.product(),magFieldHandle.product());
  // TrajectoryStateOnSurface innTSOS = mtsTransform->innerStateOnSurface(*addGsfTrk);
  // TrajectoryStateOnSurface outTSOS = mtsTransform->outerStateOnSurface(*addGsfTrk);

  // unfortunately above requires trackExtra - which is only available in RECO
  // instead we start from track vtx then propagate free state to inner/outer surface
  // no recHit hence make a reasonable approximation that
  // inner surface is pixel barrel & outer surface is tracker envelope

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
  auto gsfPropagator = std::make_unique<GsfPropagatorAdapter>(AnalyticalPropagator(magFieldHandle));
  auto extrapolator = std::make_unique<TransverseImpactPointExtrapolator>(*gsfPropagator);
  TrajectoryStateOnSurface innTSOS = gsfPropagator->propagate(freestate,innermostLayer->surface());
  StateOnTrackerBound stateOnBound(gsfPropagator.get());
  TrajectoryStateOnSurface outTSOS = stateOnBound(freestate);

  if ( innTSOS.isValid() && outTSOS.isValid() ) {
    // at seed
    TrajectoryStateOnSurface seedTSOS = extrapolator->extrapolate(*(outTSOS.freeState()), // with TSOS assert fails (not a real measurement)
                                                                  GlobalPoint(aEle.superCluster()->seed()->position().x(),
                                                                              aEle.superCluster()->seed()->position().y(),
                                                                              aEle.superCluster()->seed()->position().z()));
    //std::cout<<aEle.superCluster()->seed()->position().x()<<" | "<<aEle.superCluster()->seed()->position().y()<<" | "<<aEle.superCluster()->seed()->position().z()<<" | "<<std::endl;
    if (!seedTSOS.isValid()){
      seedTSOS = outTSOS;
    }
    TrajectoryStateOnSurface sclTSOS = extrapolator->extrapolate(*(innTSOS.freeState()), // with TSOS assert fails (not a real measurement)
                                                                 GlobalPoint(aEle.superCluster()->x(),
                                                                             aEle.superCluster()->y(),
                                                                             aEle.superCluster()->z()));
    if (!sclTSOS.isValid())
      sclTSOS = outTSOS;

    GlobalPoint seedPos, sclPos;
    multiTrajectoryStateMode::positionFromModeCartesian(seedTSOS,seedPos);
    multiTrajectoryStateMode::positionFromModeCartesian(sclTSOS,sclPos);

    scAtVtx = EleRelPointPair(aEle.superCluster()->position(),sclPos,beamSpotPos);
    seedAtCalo = EleRelPointPair(aEle.superCluster()->seed()->position(),seedPos,beamSpotPos);

    return true;
  }

  return false;
}

void MergedMuon::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",30,0.,30.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",200,-1.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);
  histo1d_["PUsummary"] = fs->make<TH1D>("PUsummary","PUsummary",99,0.,99.);

  
  muon_ = fs->make<TTree>("muonTree","muonTree");
  muon_->Branch("pT",&pT_muon,"pT/F"); 
  muon_->Branch("eta",&eta_muon,"eta/F"); 
  muon_->Branch("phi",&phi_muon,"phi/F"); 
  muon_->Branch("em",&em_muon,"em/F"); 
  muon_->Branch("emS9",&emS9_muon,"ems9/F"); 
  muon_->Branch("emS25",&emS25_muon,"emS25/F"); 
  muon_->Branch("emMax",&emMax_muon,"emMax/F"); 
  muon_->Branch("had",&had_muon,"had/F"); 
  muon_->Branch("hadS9",&hadS9_muon,"hadS9/F"); 
  muon_->Branch("hadMax",&hadMax_muon,"hadMax/F"); 
  muon_->Branch("ho",&ho_muon,"ho/F"); 
  muon_->Branch("hoS9",&hoS9_muon,"hoS9/F"); 

  mergedMuon_ = fs->make<TTree>("mergedMuonTree","mergedMuonTree");
  mergedMuon_->Branch("pT",&pT_mergedMuon,"pT/F"); 
  mergedMuon_->Branch("eta",&eta_mergedMuon,"eta/F"); 
  mergedMuon_->Branch("phi",&phi_mergedMuon,"phi/F"); 
  mergedMuon_->Branch("em",&em_mergedMuon,"em/F"); 
  mergedMuon_->Branch("emS9",&emS9_mergedMuon,"ems9/F"); 
  mergedMuon_->Branch("emS25",&emS25_mergedMuon,"emS25/F"); 
  mergedMuon_->Branch("emMax",&emMax_mergedMuon,"emMax/F"); 
  mergedMuon_->Branch("had",&had_mergedMuon,"had/F"); 
  mergedMuon_->Branch("hadS9",&hadS9_mergedMuon,"hadS9/F"); 
  mergedMuon_->Branch("hadMax",&hadMax_mergedMuon,"hadMax/F"); 
  mergedMuon_->Branch("ho",&ho_mergedMuon,"ho/F"); 
  mergedMuon_->Branch("hoS9",&hoS9_mergedMuon,"hoS9/F"); 

}


void MergedMuon::endJob() {
  purwgtFile_->Close();
}


void MergedMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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
  
  
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);
  
  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(genptcToken_, genptcHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  std::vector<reco::GenParticleRef> promptMuons;
  std::vector<reco::GenParticleRef> Muons;
  for (unsigned int idx =0; idx<genptcHandle->size();++idx){
	  const auto& genPtc = genptcHandle->refAt(idx);
	  //std::cout <<"flag : "<< genPtc->statusFlags().flags_<<" | status : "<<genPtc->status()<<" | pdg : " <<genPtc->pdgId()<<" | hard "<<genPtc->isHardProcess()  <<std::endl;    
	  if ( ( std::abs(genPtc->pdgId())==13 ) && genPtc->fromHardProcessFinalState() ) promptMuons.push_back(genPtc.castTo<reco::GenParticleRef>());
	  if ( ( std::abs(genPtc->pdgId())==13 ) && genPtc->isPromptFinalState() ) Muons.push_back(genPtc.castTo<reco::GenParticleRef>());
  }

  const reco::Vertex* pv = nullptr;
  for (auto const& v : *pvHandle) {
    if (!v.isFake() && v.ndof() > 4 && fabs(v.z()) <= 24. && v.position().Rho() <= 2.) {
      pv = &v;
      break; // 첫 번째 good vertex를 PV로 사용
    }
  }
  if (pv){
    for (const auto& muon : *muonHandle){
      int close_muon = 0;
      for (const auto& muon_sub : *muonHandle){
        double dR = sqrt(pow(muon.eta() - muon_sub.eta(),2)+pow(muon.phi() - muon_sub.phi(),2));
        if (dR < 0.1) close_muon ++;
      }
      int gen_muon = 0;
      for (auto gen : promptMuons){
        double dR = sqrt(pow(gen->eta()-muon.eta(),2)+pow(gen->phi()-muon.phi(),2));
        if (dR < 0.1){
          gen_muon ++;
        }      
      }
      if (!muon::isTightMuon(muon,*pv)) continue;
      //std::cout<<muon.isEnergyValid()<<" | "<<muon.calEnergy().towerS9<<" | "<<muon.calEnergy().emS25<<" | "<<muon.calEnergy().hadS9<<" | "<<muon.calEnergy().hoS9<<" | "<<muon.pt()<<" | "<<gen_muon<<" | "<<muon.eta()<<" | "<<muon.phi()<<" | "<<close_muon<<std::endl;
      if (gen_muon == 1 && close_muon == 1){
      	pT_muon = muon.pt();
          eta_muon = muon.eta();
          phi_muon = muon.phi();
          em_muon = muon.calEnergy().em;
          emS9_muon = muon.calEnergy().emS9;
          emS25_muon = muon.calEnergy().emS25;
          emMax_muon = muon.calEnergy().emMax;
          had_muon = muon.calEnergy().had;
          hadS9_muon = muon.calEnergy().hadS9;
          hadMax_muon = muon.calEnergy().hadMax;
          ho_muon = muon.calEnergy().ho;
          hoS9_muon = muon.calEnergy().hoS9;
          muon_->Fill();
      } 
      if (gen_muon == 2 && close_muon == 1){
      	pT_mergedMuon = muon.pt();
          eta_mergedMuon = muon.eta();
          phi_mergedMuon = muon.phi();
          em_mergedMuon = muon.calEnergy().em;
          emS9_mergedMuon = muon.calEnergy().emS9;
          emS25_mergedMuon = muon.calEnergy().emS25;
          emMax_mergedMuon = muon.calEnergy().emMax;
          had_mergedMuon = muon.calEnergy().had;
          hadS9_mergedMuon = muon.calEnergy().hadS9;
          hadMax_mergedMuon = muon.calEnergy().hadMax;
          ho_mergedMuon = muon.calEnergy().ho;
          hoS9_mergedMuon = muon.calEnergy().hoS9;
          mergedMuon_->Fill();
      } 
    
    }
  }
  return;
}

DEFINE_FWK_MODULE(MergedMuon);
