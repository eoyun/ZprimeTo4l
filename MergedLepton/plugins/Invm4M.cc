#include <memory>
#include <iostream>


#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedDEtaInSeed.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedMuonTkIsolFromCands.h"

#include "DataFormats/PatCandidates/interface/MET.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

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
#include "DataFormats/MuonReco/interface/MuonSelectors.h" // muon::segmentCompatibility


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
#include "DataFormats/MuonReco/interface/MuonPFIsolation.h"


#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

// produce TTree for merged electron training with H->AA->4e events

class Invm4M : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit Invm4M(const edm::ParameterSet&);
  virtual ~Invm4M() {}
  bool extrapolate(const reco::GsfElectron& aEle, const reco::TrackBase& addTrk,
                   const math::XYZPoint& beamSpotPos, const edm::EventSetup& iSetup,
                   EleRelPointPair& scAtVtx, EleRelPointPair& seedAtCalo);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  const edm::EDGetTokenT<edm::View<reco::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const std::vector<std::string> METfilterList_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const std::vector<edm::EDGetTokenT<edm::View<pat::PackedCandidate>>> trackCandsTokens_;
  const std::vector<MergedMuonTkIsolFromCands::PIDVeto> trackCandsVetos_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> packedPFcandToken_;
  MergedMuonTkIsolFromCands muonTkIsoCalc_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<double> prefweightToken_;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

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

  TTree* InvMFLTree_ = nullptr;

  std::vector<float> dR_gen;
  std::vector<float> gen_pt;
  std::vector<float> gen_eta;
  std::vector<float> gen_phi;
  float weight_evt;
  std::vector<float> reco_pt;
  std::vector<float> reco_eta;
  std::vector<float> reco_phi;
  std::vector<int> reco_idx;
  std::vector<bool> reco_highptid;
  std::vector<bool> reco_trackerhighptid;
  float MET;
  float MET_phi;
  float MET_phi_cor;
  float MET_cor;


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

Invm4M::Invm4M(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
trackCandsTokens_([&iConfig, this]() {
  std::vector<edm::EDGetTokenT<edm::View<pat::PackedCandidate>>> tokens;
  const auto trackCands = iConfig.getParameter<std::vector<edm::InputTag>>("trackCands");
  tokens.reserve(trackCands.size());
  for (const auto& tag : trackCands)
    tokens.push_back(consumes<edm::View<pat::PackedCandidate>>(tag));
  return tokens;
}()),
trackCandsVetos_([&iConfig]() {
  std::vector<MergedMuonTkIsolFromCands::PIDVeto> vetos;
  const auto vetoNames = iConfig.getParameter<std::vector<std::string>>("trackCandsVetos");
  vetos.reserve(vetoNames.size());
  for (const auto& name : vetoNames)
    vetos.push_back(MergedMuonTkIsolFromCands::pidVetoFromStr(name));
  return vetos;
}()),
packedPFcandToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("packedPFcand"))),
muonTkIsoCalc_(iConfig.getParameter<edm::ParameterSet>("muonTkIsoCalc")),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
prefweightToken_([&iConfig, this]() {
  const auto prefTag = iConfig.getParameter<edm::InputTag>("prefiringWeight");
  if (prefTag.label().empty())
    return edm::EDGetTokenT<double>();
  return consumes<double>(prefTag);
}()),
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

  if (trackCandsTokens_.size() != trackCandsVetos_.size())
    throw cms::Exception("ConfigError") << "trackCands and trackCandsVetos must have same size";

  usesResource("TFileService");
}

bool Invm4M::extrapolate(const reco::GsfElectron& aEle, const reco::TrackBase& addTrk,
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

void Invm4M::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",30,0.,30.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",200,-1.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);
  histo1d_["PUsummary"] = fs->make<TH1D>("PUsummary","PUsummary",99,0.,99.);

  InvMFLTree_ = fs->make<TTree>("invMFLTree","invMFLTree");
  InvMFLTree_->Branch("dR_gen",&dR_gen,32000,0); 
  InvMFLTree_->Branch("gen_pt",&gen_pt,32000,0); 
  InvMFLTree_->Branch("gen_eta",&gen_eta,32000,0); 
  InvMFLTree_->Branch("gen_phi",&gen_phi,32000,0); 
  InvMFLTree_->Branch("reco_pt",&reco_pt,32000,0); 
  InvMFLTree_->Branch("reco_eta",&reco_eta,32000,0); 
  InvMFLTree_->Branch("reco_phi",&reco_phi,32000,0); 
  InvMFLTree_->Branch("reco_idx",&reco_idx,32000,0); 
  InvMFLTree_->Branch("reco_highptid",&reco_highptid,32000,0); 
  InvMFLTree_->Branch("reco_trackerhightptid",&reco_trackerhighptid,32000,0); 
  InvMFLTree_->Branch("weight_evt",&weight_evt,"weight_evt/F");
  InvMFLTree_->Branch("MET",&MET,"MET/F");
  InvMFLTree_->Branch("MET_phi",&MET_phi,"MET_phi/F");
  InvMFLTree_->Branch("MET_cor",&MET_cor,"MET_cor/F");
  InvMFLTree_->Branch("MET_phi_cor",&MET_phi_cor,"MET_phi_cor/F");

}


void Invm4M::endJob() {
  purwgtFile_->Close();
}


void Invm4M::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  double aWeight = 1.;
  double mcweight = 0.;
  if (isMC_) {
    double prefiringweight = 1.;
    if (!prefweightToken_.isUninitialized()) {
      edm::Handle<double> theprefweight;
      if (iEvent.getByToken(prefweightToken_, theprefweight) && theprefweight.isValid())
        prefiringweight = *theprefweight;
    }

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    mcweight = genInfo->weight();
    //std::cout<<mcweight<<" : weight"<<std::endl;

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

  edm::Handle<edm::View<pat::MET>> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  edm::Handle<edm::TriggerResults> METfilterHandle;
  iEvent.getByToken(METfilterToken_,METfilterHandle);
  edm::TriggerNames METfilters = iEvent.triggerNames(*METfilterHandle);
  
  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(genptcToken_, genptcHandle);

  std::vector<edm::Handle<edm::View<pat::PackedCandidate>>> trackCandsHandles(trackCandsTokens_.size());
  for (size_t iCand = 0; iCand < trackCandsTokens_.size(); ++iCand)
    iEvent.getByToken(trackCandsTokens_[iCand], trackCandsHandles[iCand]);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  std::vector<reco::GenParticleRef> promptMuons;
  std::vector<reco::GenParticleRef> Muons;

  unsigned int nPassedFilters = 0;

  for (unsigned int iTrig = 0; iTrig < METfilterHandle.product()->size(); iTrig++) {
    const std::string trigname = METfilters.triggerName(iTrig);
    //std::cout<<trigname <<"filter name of MET"<<std::endl;
    if (METfilterHandle.product()->accept(iTrig)) {
      for (const auto& filterName : METfilterList_) {
	//std::cout<< filterName<<" : name "<<std::endl;
	if (trigname.find(filterName) != std::string::npos)
	  {nPassedFilters++;
	  std::cout<<"passed!!"<<std::endl;}
	//else {std::cout<<"fail!!"<<std::endl;}
	
      }
    }
  }

  //if (nPassedFilters!=METfilterList_.size())
  //  return;

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
    if (promptMuons.size() != 4) return;
  
    // ---------- pair 찾기 (원본 로직 그대로) ----------
    int p1a = -1, p1b = -1;
    for (size_t i = 0; i < 4 && p1a < 0; ++i) {
      for (size_t j = i + 1; j < 4; ++j) {
        if (reco::deltaR(*(promptMuons.at(i)), *(promptMuons.at(j))) < 0.1f) {
          p1a = (int)i; p1b = (int)j;
          break;
        }
      }
    }
    if (p1a < 0) return; // 가까운 pair 없음
  
    // 나머지 둘이 자동으로 둘째 pair
    int p2a = -1, p2b = -1;
    for (size_t k = 0; k < 4; ++k) {
      if ((int)k == p1a || (int)k == p1b) continue;
      if (p2a < 0) p2a = (int)k;
      else         p2b = (int)k;
    }
  
    // (선택) 둘째 pair도 가까운지 확인
    // if (reco::deltaR(*(promptMuons.at(p2a)), *(promptMuons.at(p2b))) >= 0.1f) return;
  
    std::array<std::pair<int,int>, 2> pairs = {{ {p1a, p1b}, {p2a, p2b} }};
  
    // ---------- pair별 reco matching ----------
    auto matchPair = [&](const reco::Candidate& genA,
                         const reco::Candidate& genB,
                         int& matchA, int& matchB) {
      matchA = -1; matchB = -1;
      float bestA = 999.f, secA = 999.f;
      float bestB = 999.f, secB = 999.f;
      int   bestAi = -1, secAi = -1;
      int   bestBi = -1, secBi = -1;
  
      const float dR_cut = 0.03f;
  
      for (size_t i = 0; i < muonHandle->size(); ++i) {
        const reco::Muon& m = (*muonHandle)[i];
  
        float dA = reco::deltaR(m, genA);
        if (dA < dR_cut) {
          if (dA < bestA)     { secA = bestA; secAi = bestAi; bestA = dA; bestAi = (int)i; }
          else if (dA < secA) { secA = dA;    secAi = (int)i; }
        }
        float dB = reco::deltaR(m, genB);
        if (dB < dR_cut) {
          if (dB < bestB)     { secB = bestB; secBi = bestBi; bestB = dB; bestBi = (int)i; }
          else if (dB < secB) { secB = dB;    secBi = (int)i; }
        }
      }
  
      matchA = bestAi;
      matchB = bestBi;
  
      if (matchA >= 0 && matchA == matchB) {
        if (bestA <= bestB) {
          matchB = (secBi != matchA) ? secBi : -1;
        } else {
          matchA = (secAi != matchB) ? secAi : -1;
        }
      }
      if (matchA >= 0 && matchA == matchB) matchB = -1;
    };
  
    std::vector<int> match_idx(4, -1);
    for (auto& p : pairs) {
      int mA = -1, mB = -1;
      matchPair(*(promptMuons.at(p.first)),
                *(promptMuons.at(p.second)),
                mA, mB);
      match_idx[p.first]  = mA;
      match_idx[p.second] = mB;
    }
  
    // ---------- tree 채우기 ----------
    std::array<int, 4> order = {{
      pairs[0].first, pairs[0].second,
      pairs[1].first, pairs[1].second
    }};
  
    gen_pt.clear();   gen_eta.clear();  gen_phi.clear();
    reco_pt.clear();  reco_eta.clear(); reco_phi.clear();
    reco_idx.clear();
    reco_highptid.clear();
    reco_trackerhighptid.clear();
    dR_gen.clear();
  
    for (int idx : order) {
      const auto& gen = *(promptMuons.at(idx));
      gen_pt .push_back(gen.pt());
      gen_eta.push_back(gen.eta());
      gen_phi.push_back(gen.phi());
  
      if (match_idx[idx] >= 0) {
        const reco::Muon& m = (*muonHandle)[match_idx[idx]];
        reco_pt .push_back(m.pt());
        reco_eta.push_back(m.eta());
        reco_phi.push_back(m.phi());
        reco_idx.push_back(match_idx[idx]);
        reco_highptid       .push_back( muon::isHighPtMuon(m, *pv)        ? 1 : 0 );
        reco_trackerhighptid.push_back( muon::isTrackerHighPtMuon(m, *pv) ? 1 : 0 );
      } else {
        reco_pt .push_back(-999.f);
        reco_eta.push_back(-999.f);
        reco_phi.push_back(-999.f);
        reco_idx.push_back(-1);
        reco_highptid       .push_back(-1);
        reco_trackerhighptid.push_back(-1);
      }
    }
  
    for (auto& p : pairs)
      dR_gen.push_back( reco::deltaR(*(promptMuons.at(p.first)),
                                     *(promptMuons.at(p.second))) );
  
    weight_evt = mcweight;
    const auto& aMET = metHandle->at(0);
    MET         = aMET.pt();
    MET_phi     = aMET.phi();
    MET_cor     = aMET.corPt(pat::MET::Type1);
    MET_phi_cor = aMET.corPhi(pat::MET::TypeXY);
  
    InvMFLTree_->Fill();
  }
  return;
}

DEFINE_FWK_MODULE(Invm4M);
