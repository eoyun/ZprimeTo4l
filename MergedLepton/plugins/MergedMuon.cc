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
  const edm::EDGetTokenT<edm::View<pat::MET>> metToken_;
  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;
  const std::vector<std::string> METfilterList_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
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
  bool isGlobal_muon;
  bool isTracker_muon;
  int numOfMatchedStations_muon;
  int numOfChambers_muon;
  double segCompatibility_muon;
  double caloCompatibility_muon;
  float chi2LocalPosition_muon;
  float trkKink_muon;
  float glbKink_muon;
  float pfrelIso04_muon;
  float gen_eta_muon;
  float gen_phi_muon;
  float gen_pt_muon;
  float weight_muon;

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
  bool isGlobal_mergedMuon;
  bool isTracker_mergedMuon;
  int numOfMatchedStations_mergedMuon;
  int numOfChambers_mergedMuon;
  double segCompatibility_mergedMuon;
  double caloCompatibility_mergedMuon;
  float chi2LocalPosition_mergedMuon;
  float trkKink_mergedMuon;
  float glbKink_mergedMuon;
  float gen_eta_mergedMuon;
  float gen_phi_mergedMuon;
  float gen_pt_mergedMuon;
  float gen_sub_eta_mergedMuon;
  float gen_sub_phi_mergedMuon;
  float gen_sub_pt_mergedMuon;
  float weight_mergedMuon;

  TTree* muonEfficiencyTree_ = nullptr;

  float dR_gen;
  int num_reco_muon;
  bool flag_Id;
  bool flag_Id_woID;
  bool flag_Id_any;
  float eff_weight;
  float eff_gen_1_pt;
  float eff_gen_1_eta;
  float eff_gen_1_phi;
  float eff_gen_2_pt;
  float eff_gen_2_eta;
  float eff_gen_2_phi;
  float eff_reco_1_pt;
  float eff_reco_1_eta;
  float eff_reco_1_phi;
  float eff_reco_1_dr;
  int eff_reco_1_idx;
  bool eff_reco_1_highptid;
  bool eff_reco_1_trackerhighptid;
  float eff_reco_1_em;
  float eff_reco_1_emS9;
  float eff_reco_1_emS25;
  float eff_reco_1_emMax;
  float eff_reco_1_had;
  float eff_reco_1_hadS9;
  float eff_reco_1_hadMax;
  float eff_reco_1_ho;
  float eff_reco_1_hoS9;
  bool eff_reco_1_isGlobal;
  bool eff_reco_1_isTracker;
  bool eff_reco_1_isPFMuon;
  bool eff_reco_1_isPFIsolationValid;
  int eff_reco_1_numOfMatchedStations;
  int eff_reco_1_numOfChambers;
  double eff_reco_1_segCompatibility;
  double eff_reco_1_caloCompatibility;
  float eff_reco_1_chi2LocalPosition;
  float eff_reco_1_trkKink;
  float eff_reco_1_glbKink;
  float eff_reco_1_pfrelIso04;
  bool eff_reco_1_isLoose;
  bool eff_reco_1_isMedium;
  bool eff_reco_1_isTight;
  int eff_reco_1_nValidMuonHits;
  float eff_reco_1_normChi2_global;
  int eff_reco_1_nValidPixelHits;
  int eff_reco_1_trkLayers;
  int eff_reco_1_pixelLayers;
  float eff_reco_1_normChi2_inner;
  float eff_reco_1_validFraction;
  float eff_reco_1_tunePtErrorOverPt;
  float eff_reco_1_dxy_PV;
  float eff_reco_1_dz_PV;
  int eff_reco_1_nMatches;
  float eff_reco_1_innerTrk_pt;
  float eff_reco_1_innerTrk_eta;
  float eff_reco_1_innerTrk_phi;
  float eff_reco_1_outerTrk_pt;
  float eff_reco_1_outerTrk_eta;
  float eff_reco_1_outerTrk_phi;
  float eff_reco_1_globalTrk_pt;
  float eff_reco_1_globalTrk_eta;
  float eff_reco_1_globalTrk_phi;
  float eff_reco_1_bestTrk_pt;
  float eff_reco_1_bestTrk_eta;
  float eff_reco_1_bestTrk_phi;
  float eff_reco_1_tunePTrk_pt;
  float eff_reco_1_tunePTrk_eta;
  float eff_reco_1_tunePTrk_phi;
  float eff_reco_2_pt;
  float eff_reco_2_eta;
  float eff_reco_2_phi;
  float eff_reco_2_dr;
  int eff_reco_2_idx;
  bool eff_reco_2_highptid;
  bool eff_reco_2_trackerhighptid;
  float eff_reco_2_em;
  float eff_reco_2_emS9;
  float eff_reco_2_emS25;
  float eff_reco_2_emMax;
  float eff_reco_2_had;
  float eff_reco_2_hadS9;
  float eff_reco_2_hadMax;
  float eff_reco_2_ho;
  float eff_reco_2_hoS9;
  bool eff_reco_2_isGlobal;
  bool eff_reco_2_isTracker;
  bool eff_reco_2_isPFMuon;
  bool eff_reco_2_isPFIsolationValid;
  int eff_reco_2_numOfMatchedStations;
  int eff_reco_2_numOfChambers;
  double eff_reco_2_segCompatibility;
  double eff_reco_2_caloCompatibility;
  float eff_reco_2_chi2LocalPosition;
  float eff_reco_2_trkKink;
  float eff_reco_2_glbKink;
  float eff_reco_2_pfrelIso04;
  bool eff_reco_2_isLoose;
  bool eff_reco_2_isMedium;
  bool eff_reco_2_isTight;
  int eff_reco_2_nValidMuonHits;
  float eff_reco_2_normChi2_global;
  int eff_reco_2_nValidPixelHits;
  int eff_reco_2_trkLayers;
  int eff_reco_2_pixelLayers;
  float eff_reco_2_normChi2_inner;
  float eff_reco_2_validFraction;
  float eff_reco_2_tunePtErrorOverPt;
  float eff_reco_2_dxy_PV;
  float eff_reco_2_dz_PV;
  int eff_reco_2_nMatches;
  float eff_reco_2_innerTrk_pt;
  float eff_reco_2_innerTrk_eta;
  float eff_reco_2_innerTrk_phi;
  float eff_reco_2_outerTrk_pt;
  float eff_reco_2_outerTrk_eta;
  float eff_reco_2_outerTrk_phi;
  float eff_reco_2_globalTrk_pt;
  float eff_reco_2_globalTrk_eta;
  float eff_reco_2_globalTrk_phi;
  float eff_reco_2_bestTrk_pt;
  float eff_reco_2_bestTrk_eta;
  float eff_reco_2_bestTrk_phi;
  float eff_reco_2_tunePTrk_pt;
  float eff_reco_2_tunePTrk_eta;
  float eff_reco_2_tunePTrk_phi;
  bool add1_track_flag;
  float add1_track_pt;
  float add1_track_eta;
  float add1_track_phi;
  bool add2_track_flag;
  float add2_track_pt;
  float add2_track_eta;
  float add2_track_phi;
  float eff_MET;
  float eff_MET_phi;
  float eff_MET_cor;
  float eff_MET_cor_XY_phi;

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
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
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

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",30,0.,30.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",200,-1.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);

  
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
  muon_->Branch("isGlobal",&isGlobal_muon,"isGlobal/B"); 
  muon_->Branch("isTracker",&isTracker_muon,"isTracker/B"); 
  muon_->Branch("numOfMatchedStations",&numOfMatchedStations_muon,"numOfMatchedStations/I"); 
  muon_->Branch("numOfChambers",&numOfChambers_muon,"numOfChambers/I"); 
  muon_->Branch("segCompatibility",&segCompatibility_muon,"segCompatibility/D"); 
  muon_->Branch("caloCompatibility",&caloCompatibility_muon,"caloCompatibility/D"); 
  muon_->Branch("chi2LocalPosition",&chi2LocalPosition_muon,"chi2LocalPosition/F"); 
  muon_->Branch("trkKink",&trkKink_muon,"trkKink/F"); 
  muon_->Branch("glbKink",&glbKink_muon,"glbKink/F"); 
  muon_->Branch("pfrelIso04",&pfrelIso04_muon,"pfrelIso04/F"); 
  muon_->Branch("gen_phi",&gen_phi_muon,"gen_phi/F");
  muon_->Branch("gen_eta",&gen_eta_muon,"gen_eta/F");
  muon_->Branch("gen_pt",&gen_pt_muon,"gen_pt/F");
  muon_->Branch("weight",&weight_muon,"weight/F");

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
  mergedMuon_->Branch("isGlobal",&isGlobal_mergedMuon,"isGlobal/B"); 
  mergedMuon_->Branch("isTracker",&isTracker_mergedMuon,"isTracker/B"); 
  mergedMuon_->Branch("numOfMatchedStations",&numOfMatchedStations_mergedMuon,"numOfMatchedStations/I"); 
  mergedMuon_->Branch("numOfChambers",&numOfChambers_mergedMuon,"numOfChambers/I"); 
  mergedMuon_->Branch("segCompatibility",&segCompatibility_mergedMuon,"segCompatibility/D"); 
  mergedMuon_->Branch("caloCompatibility",&caloCompatibility_mergedMuon,"caloCompatibility/D"); 
  mergedMuon_->Branch("chi2LocalPosition",&chi2LocalPosition_mergedMuon,"chi2LocalPosition/F"); 
  mergedMuon_->Branch("trkKink",&trkKink_mergedMuon,"trkKink/F"); 
  mergedMuon_->Branch("glbKink",&glbKink_mergedMuon,"glbKink/F"); 
  mergedMuon_->Branch("gen_phi",&gen_phi_mergedMuon,"gen_phi/F");
  mergedMuon_->Branch("gen_eta",&gen_eta_mergedMuon,"gen_eta/F");
  mergedMuon_->Branch("gen_pt",&gen_pt_mergedMuon,"gen_pt/F");
  mergedMuon_->Branch("gen_sub_phi",&gen_sub_phi_mergedMuon,"gen_sub_phi/F");
  mergedMuon_->Branch("gen_sub_eta",&gen_sub_eta_mergedMuon,"gen_sub_eta/F");
  mergedMuon_->Branch("gen_sub_pt",&gen_sub_pt_mergedMuon,"gen_sub_pt/F");
  mergedMuon_->Branch("weight",&weight_mergedMuon,"weight/F");

  muonEfficiencyTree_ = fs->make<TTree>("muonEfficiencyTree","muonEfficiencyTree");
  muonEfficiencyTree_->Branch("dR_gen",&dR_gen,"dR_gen/F");
  muonEfficiencyTree_->Branch("num_reco_muon",&num_reco_muon,"num_reco_muon/I");
  muonEfficiencyTree_->Branch("flag_Id",&flag_Id,"flag_Id/B");
  muonEfficiencyTree_->Branch("flag_Id_woID",&flag_Id_woID,"flag_Id_woID/B");
  muonEfficiencyTree_->Branch("flag_Id_any",&flag_Id_any,"flag_Id_any/B");
  muonEfficiencyTree_->Branch("eff_gen_1_pt",&eff_gen_1_pt,"eff_gen_1_pt/F");
  muonEfficiencyTree_->Branch("eff_gen_1_phi",&eff_gen_1_phi,"eff_gen_1_phi/F");
  muonEfficiencyTree_->Branch("eff_gen_1_eta",&eff_gen_1_eta,"eff_gen_1_eta/F");
  muonEfficiencyTree_->Branch("eff_gen_2_pt",&eff_gen_2_pt,"eff_gen_2_pt/F");
  muonEfficiencyTree_->Branch("eff_gen_2_phi",&eff_gen_2_phi,"eff_gen_2_phi/F");
  muonEfficiencyTree_->Branch("eff_gen_2_eta",&eff_gen_2_eta,"eff_gen_2_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_pt",&eff_reco_1_pt,"eff_reco_1_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_phi",&eff_reco_1_phi,"eff_reco_1_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_1_dr",&eff_reco_1_dr,"eff_reco_1_dr/F");
  muonEfficiencyTree_->Branch("eff_reco_1_eta",&eff_reco_1_eta,"eff_reco_1_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_idx",&eff_reco_1_idx,"eff_reco_1_idx/I");
  muonEfficiencyTree_->Branch("eff_reco_1_highptid",&eff_reco_1_highptid,"eff_reco_1_highptid/B");
  muonEfficiencyTree_->Branch("eff_reco_1_trackerhighptid",&eff_reco_1_trackerhighptid,"eff_reco_1_trackerhighptid/B");
  muonEfficiencyTree_->Branch("eff_reco_1_em",&eff_reco_1_em,"eff_reco_1_em/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_emS9",&eff_reco_1_emS9,"eff_reco_1_ems9/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_emS25",&eff_reco_1_emS25,"eff_reco_1_emS25/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_emMax",&eff_reco_1_emMax,"eff_reco_1_emMax/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_had",&eff_reco_1_had,"eff_reco_1_had/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_hadS9",&eff_reco_1_hadS9,"eff_reco_1_hadS9/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_hadMax",&eff_reco_1_hadMax,"eff_reco_1_hadMax/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_ho",&eff_reco_1_ho,"eff_reco_1_ho/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_hoS9",&eff_reco_1_hoS9,"eff_reco_1_hoS9/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_isGlobal",&eff_reco_1_isGlobal,"eff_reco_1_isGlobal/B"); 
  muonEfficiencyTree_->Branch("eff_reco_1_isTracker",&eff_reco_1_isTracker,"eff_reco_1_isTracker/B"); 
  muonEfficiencyTree_->Branch("eff_reco_1_isPFMuon",&eff_reco_1_isPFMuon,"eff_reco_1_isPFMuon/B"); 
  muonEfficiencyTree_->Branch("eff_reco_1_isPFIsolationValid",&eff_reco_1_isPFIsolationValid,"eff_reco_1_isPFIsolationValid/B"); 
  muonEfficiencyTree_->Branch("eff_reco_1_numOfMatchedStations",&eff_reco_1_numOfMatchedStations,"eff_reco_1_numOfMatchedStations/I"); 
  muonEfficiencyTree_->Branch("eff_reco_1_numOfChambers",&eff_reco_1_numOfChambers,"eff_reco_1_numOfChambers/I"); 
  muonEfficiencyTree_->Branch("eff_reco_1_segCompatibility",&eff_reco_1_segCompatibility,"eff_reco_1_segCompatibility/D"); 
  muonEfficiencyTree_->Branch("eff_reco_1_caloCompatibility",&eff_reco_1_caloCompatibility,"eff_reco_1_caloCompatibility/D"); 
  muonEfficiencyTree_->Branch("eff_reco_1_chi2LocalPosition",&eff_reco_1_chi2LocalPosition,"eff_reco_1_chi2LocalPosition/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_trkKink",&eff_reco_1_trkKink,"eff_reco_1_trkKink/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_glbKink",&eff_reco_1_glbKink,"eff_reco_1_glbKink/F"); 
  muonEfficiencyTree_->Branch("eff_reco_1_pfrelIso04",&eff_reco_1_pfrelIso04,"eff_reco_1_pfrelIso04/F");
  muonEfficiencyTree_->Branch("eff_reco_1_isLoose",&eff_reco_1_isLoose,"eff_reco_1_isLoose/O");
  muonEfficiencyTree_->Branch("eff_reco_1_isMedium",&eff_reco_1_isMedium,"eff_reco_1_isMedium/O");
  muonEfficiencyTree_->Branch("eff_reco_1_isTight",&eff_reco_1_isTight,"eff_reco_1_isTight/O");
  muonEfficiencyTree_->Branch("eff_reco_1_nValidMuonHits",&eff_reco_1_nValidMuonHits,"eff_reco_1_nValidMuonHits/I");
  muonEfficiencyTree_->Branch("eff_reco_1_normChi2_global",&eff_reco_1_normChi2_global,"eff_reco_1_normChi2_global/F");
  muonEfficiencyTree_->Branch("eff_reco_1_nValidPixelHits",&eff_reco_1_nValidPixelHits,"eff_reco_1_nValidPixelHits/I");
  muonEfficiencyTree_->Branch("eff_reco_1_trkLayers",&eff_reco_1_trkLayers,"eff_reco_1_trkLayers/I");
  muonEfficiencyTree_->Branch("eff_reco_1_pixelLayers",&eff_reco_1_pixelLayers,"eff_reco_1_pixelLayers/I");
  muonEfficiencyTree_->Branch("eff_reco_1_normChi2_inner",&eff_reco_1_normChi2_inner,"eff_reco_1_normChi2_inner/F");
  muonEfficiencyTree_->Branch("eff_reco_1_validFraction",&eff_reco_1_validFraction,"eff_reco_1_validFraction/F");
  muonEfficiencyTree_->Branch("eff_reco_1_tunePtErrorOverPt",&eff_reco_1_tunePtErrorOverPt,"eff_reco_1_tunePtErrorOverPt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_dxy_PV",&eff_reco_1_dxy_PV,"eff_reco_1_dxy_PV/F");
  muonEfficiencyTree_->Branch("eff_reco_1_dz_PV",&eff_reco_1_dz_PV,"eff_reco_1_dz_PV/F");
  muonEfficiencyTree_->Branch("eff_reco_1_nMatches",&eff_reco_1_nMatches,"eff_reco_1_nMatches/I");
  muonEfficiencyTree_->Branch("eff_reco_1_innerTrk_pt",&eff_reco_1_innerTrk_pt,"eff_reco_1_innerTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_innerTrk_eta",&eff_reco_1_innerTrk_eta,"eff_reco_1_innerTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_innerTrk_phi",&eff_reco_1_innerTrk_phi,"eff_reco_1_innerTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_1_outerTrk_pt",&eff_reco_1_outerTrk_pt,"eff_reco_1_outerTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_outerTrk_eta",&eff_reco_1_outerTrk_eta,"eff_reco_1_outerTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_outerTrk_phi",&eff_reco_1_outerTrk_phi,"eff_reco_1_outerTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_1_globalTrk_pt",&eff_reco_1_globalTrk_pt,"eff_reco_1_globalTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_globalTrk_eta",&eff_reco_1_globalTrk_eta,"eff_reco_1_globalTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_globalTrk_phi",&eff_reco_1_globalTrk_phi,"eff_reco_1_globalTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_1_bestTrk_pt",&eff_reco_1_bestTrk_pt,"eff_reco_1_bestTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_bestTrk_eta",&eff_reco_1_bestTrk_eta,"eff_reco_1_bestTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_bestTrk_phi",&eff_reco_1_bestTrk_phi,"eff_reco_1_bestTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_1_tunePTrk_pt",&eff_reco_1_tunePTrk_pt,"eff_reco_1_tunePTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_1_tunePTrk_eta",&eff_reco_1_tunePTrk_eta,"eff_reco_1_tunePTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_1_tunePTrk_phi",&eff_reco_1_tunePTrk_phi,"eff_reco_1_tunePTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_2_pt",&eff_reco_2_pt,"eff_reco_2_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_phi",&eff_reco_2_phi,"eff_reco_2_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_2_eta",&eff_reco_2_eta,"eff_reco_2_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_2_dr",&eff_reco_2_dr,"eff_reco_2_dr/F");
  muonEfficiencyTree_->Branch("eff_reco_2_idx",&eff_reco_2_idx,"eff_reco_2_idx/I");
  muonEfficiencyTree_->Branch("eff_reco_2_highptid",&eff_reco_2_highptid,"eff_reco_2_highptid/B");
  muonEfficiencyTree_->Branch("eff_reco_2_trackerhighptid",&eff_reco_2_trackerhighptid,"eff_reco_2_trackerhighptid/B");
  muonEfficiencyTree_->Branch("eff_reco_2_em",&eff_reco_2_em,"eff_reco_2_em/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_emS9",&eff_reco_2_emS9,"eff_reco_2_ems9/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_emS25",&eff_reco_2_emS25,"eff_reco_2_emS25/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_emMax",&eff_reco_2_emMax,"eff_reco_2_emMax/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_had",&eff_reco_2_had,"eff_reco_2_had/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_hadS9",&eff_reco_2_hadS9,"eff_reco_2_hadS9/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_hadMax",&eff_reco_2_hadMax,"eff_reco_2_hadMax/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_ho",&eff_reco_2_ho,"eff_reco_2_ho/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_hoS9",&eff_reco_2_hoS9,"eff_reco_2_hoS9/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_isGlobal",&eff_reco_2_isGlobal,"eff_reco_2_isGlobal/B"); 
  muonEfficiencyTree_->Branch("eff_reco_2_isTracker",&eff_reco_2_isTracker,"eff_reco_2_isTracker/B"); 
  muonEfficiencyTree_->Branch("eff_reco_2_isPFMuon",&eff_reco_2_isPFMuon,"eff_reco_2_isPFMuon/B"); 
  muonEfficiencyTree_->Branch("eff_reco_2_isPFIsolationValid",&eff_reco_2_isPFIsolationValid,"eff_reco_2_isPFIsolationValid/B"); 
  muonEfficiencyTree_->Branch("eff_reco_2_numOfMatchedStations",&eff_reco_2_numOfMatchedStations,"eff_reco_2_numOfMatchedStations/I"); 
  muonEfficiencyTree_->Branch("eff_reco_2_numOfChambers",&eff_reco_2_numOfChambers,"eff_reco_2_numOfChambers/I"); 
  muonEfficiencyTree_->Branch("eff_reco_2_segCompatibility",&eff_reco_2_segCompatibility,"eff_reco_2_segCompatibility/D"); 
  muonEfficiencyTree_->Branch("eff_reco_2_caloCompatibility",&eff_reco_2_caloCompatibility,"eff_reco_2_caloCompatibility/D"); 
  muonEfficiencyTree_->Branch("eff_reco_2_chi2LocalPosition",&eff_reco_2_chi2LocalPosition,"eff_reco_2_chi2LocalPosition/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_trkKink",&eff_reco_2_trkKink,"eff_reco_2_trkKink/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_glbKink",&eff_reco_2_glbKink,"eff_reco_2_glbKink/F"); 
  muonEfficiencyTree_->Branch("eff_reco_2_pfrelIso04",&eff_reco_2_pfrelIso04,"eff_reco_2_pfrelIso04/F");
  muonEfficiencyTree_->Branch("eff_reco_2_isLoose",&eff_reco_2_isLoose,"eff_reco_2_isLoose/O");
  muonEfficiencyTree_->Branch("eff_reco_2_isMedium",&eff_reco_2_isMedium,"eff_reco_2_isMedium/O");
  muonEfficiencyTree_->Branch("eff_reco_2_isTight",&eff_reco_2_isTight,"eff_reco_2_isTight/O");
  muonEfficiencyTree_->Branch("eff_reco_2_nValidMuonHits",&eff_reco_2_nValidMuonHits,"eff_reco_2_nValidMuonHits/I");
  muonEfficiencyTree_->Branch("eff_reco_2_normChi2_global",&eff_reco_2_normChi2_global,"eff_reco_2_normChi2_global/F");
  muonEfficiencyTree_->Branch("eff_reco_2_nValidPixelHits",&eff_reco_2_nValidPixelHits,"eff_reco_2_nValidPixelHits/I");
  muonEfficiencyTree_->Branch("eff_reco_2_trkLayers",&eff_reco_2_trkLayers,"eff_reco_2_trkLayers/I");
  muonEfficiencyTree_->Branch("eff_reco_2_pixelLayers",&eff_reco_2_pixelLayers,"eff_reco_2_pixelLayers/I");
  muonEfficiencyTree_->Branch("eff_reco_2_normChi2_inner",&eff_reco_2_normChi2_inner,"eff_reco_2_normChi2_inner/F");
  muonEfficiencyTree_->Branch("eff_reco_2_validFraction",&eff_reco_2_validFraction,"eff_reco_2_validFraction/F");
  muonEfficiencyTree_->Branch("eff_reco_2_tunePtErrorOverPt",&eff_reco_2_tunePtErrorOverPt,"eff_reco_2_tunePtErrorOverPt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_dxy_PV",&eff_reco_2_dxy_PV,"eff_reco_2_dxy_PV/F");
  muonEfficiencyTree_->Branch("eff_reco_2_dz_PV",&eff_reco_2_dz_PV,"eff_reco_2_dz_PV/F");
  muonEfficiencyTree_->Branch("eff_reco_2_nMatches",&eff_reco_2_nMatches,"eff_reco_2_nMatches/I");
  muonEfficiencyTree_->Branch("eff_reco_2_innerTrk_pt",&eff_reco_2_innerTrk_pt,"eff_reco_2_innerTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_innerTrk_eta",&eff_reco_2_innerTrk_eta,"eff_reco_2_innerTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_2_innerTrk_phi",&eff_reco_2_innerTrk_phi,"eff_reco_2_innerTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_2_outerTrk_pt",&eff_reco_2_outerTrk_pt,"eff_reco_2_outerTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_outerTrk_eta",&eff_reco_2_outerTrk_eta,"eff_reco_2_outerTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_2_outerTrk_phi",&eff_reco_2_outerTrk_phi,"eff_reco_2_outerTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_2_globalTrk_pt",&eff_reco_2_globalTrk_pt,"eff_reco_2_globalTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_globalTrk_eta",&eff_reco_2_globalTrk_eta,"eff_reco_2_globalTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_2_globalTrk_phi",&eff_reco_2_globalTrk_phi,"eff_reco_2_globalTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_2_bestTrk_pt",&eff_reco_2_bestTrk_pt,"eff_reco_2_bestTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_bestTrk_eta",&eff_reco_2_bestTrk_eta,"eff_reco_2_bestTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_2_bestTrk_phi",&eff_reco_2_bestTrk_phi,"eff_reco_2_bestTrk_phi/F");
  muonEfficiencyTree_->Branch("eff_reco_2_tunePTrk_pt",&eff_reco_2_tunePTrk_pt,"eff_reco_2_tunePTrk_pt/F");
  muonEfficiencyTree_->Branch("eff_reco_2_tunePTrk_eta",&eff_reco_2_tunePTrk_eta,"eff_reco_2_tunePTrk_eta/F");
  muonEfficiencyTree_->Branch("eff_reco_2_tunePTrk_phi",&eff_reco_2_tunePTrk_phi,"eff_reco_2_tunePTrk_phi/F");
  muonEfficiencyTree_->Branch("add1_trk_flag",&add1_track_flag,"add1_trk_flag/B");
  muonEfficiencyTree_->Branch("add1_trk_pt",&add1_track_pt,"add1_trk_pt/F");
  muonEfficiencyTree_->Branch("add1_trk_eta",&add1_track_eta,"add1_trk_eta/F");
  muonEfficiencyTree_->Branch("add1_trk_phi",&add1_track_phi,"add1_trk_phi/F");
  muonEfficiencyTree_->Branch("add2_trk_flag",&add2_track_flag,"add2_trk_flag/B");
  muonEfficiencyTree_->Branch("add2_trk_pt",&add2_track_pt,"add2_trk_pt/F");
  muonEfficiencyTree_->Branch("add2_trk_eta",&add2_track_eta,"add2_trk_eta/F");
  muonEfficiencyTree_->Branch("add2_trk_phi",&add2_track_phi,"add2_trk_phi/F");
  muonEfficiencyTree_->Branch("eff_weight",&eff_weight,"eff_weight/F");
  muonEfficiencyTree_->Branch("eff_MET",&eff_MET,"eff_MET/F");
  muonEfficiencyTree_->Branch("eff_MET_phi",&eff_MET_phi,"eff_MET_phi/F");
  muonEfficiencyTree_->Branch("eff_MET_cor",&eff_MET_cor,"eff_MET_cor/F");
  muonEfficiencyTree_->Branch("eff_MET_cor_XY_phi",&eff_MET_cor_XY_phi,"eff_MET_cor_XY_phi/F");
}


void MergedMuon::endJob() {
}


void MergedMuon::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  double mcweight = 1.;
  if (isMC_) {
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    mcweight = genInfo->weight();
  }
  
  
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<pat::MET>> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  // edm::Handle<edm::TriggerResults> METfilterHandle;
  // iEvent.getByToken(METfilterToken_,METfilterHandle);
  // edm::TriggerNames METfilters = iEvent.triggerNames(*METfilterHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(genptcToken_, genptcHandle);

  std::vector<edm::Handle<edm::View<pat::PackedCandidate>>> trackCandsHandles(trackCandsTokens_.size());
  for (size_t iCand = 0; iCand < trackCandsTokens_.size(); ++iCand)
    iEvent.getByToken(trackCandsTokens_[iCand], trackCandsHandles[iCand]);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  std::vector<reco::GenParticleRef> promptMuons;
  std::vector<reco::GenParticleRef> Muons;

  // unsigned int nPassedFilters = 0;
  //
  // for (unsigned int iTrig = 0; iTrig < METfilterHandle.product()->size(); iTrig++) {
  //   const std::string trigname = METfilters.triggerName(iTrig);
  //   if (METfilterHandle.product()->accept(iTrig)) {
  //     for (const auto& filterName : METfilterList_) {
  //       if (trigname.find(filterName) != std::string::npos)
  //         nPassedFilters++;
  //     }
  //   }
  // }
  //
  // if (nPassedFilters != METfilterList_.size())
  //   return;
  //
  // edm::Handle<edm::TriggerResults> trigResultHandle;
  // iEvent.getByToken(triggerToken_, trigResultHandle);
  // edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);
  // bool isFired = false;
  // for (unsigned int iTrig = 0; iTrig != trigResultHandle.product()->size(); iTrig++) {
  //   const std::string trigName = trigList.triggerName(iTrig);
  //   for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
  //     const std::string& wanted = trigList_.at(jTrig);
  //     const std::string stem = wanted.substr(0, wanted.find("*"));
  //     if (trigName.find(stem) != std::string::npos) {
  //       if (trigResultHandle.product()->accept(iTrig))
  //         isFired = true;
  //     }
  //   }
  // }
  // if (!isFired)
  //   return;

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
    const auto& ttBuilder = iSetup.getData(ttbToken_);
    std::vector<bool> usedPrompt(promptMuons.size(), false);
    for (size_t i = 0; i < promptMuons.size(); ++i) {
      if (usedPrompt[i]) continue;
      int index = -1;

      // gen pair 찾기 (이미 페어로 잡힌 gen muon은 재사용하지 않음)
      for (size_t j = i + 1; j < promptMuons.size(); ++j) {
        if (usedPrompt[j]) continue;
        const auto& gen1 = *(promptMuons.at(i));
        const auto& gen2 = *(promptMuons.at(j));

        if (reco::deltaR(gen1, gen2) >= 0.1) continue;

        // leading / subleading pt cut
        float pt1 = gen1.pt();
        float pt2 = gen2.pt();
        float leadPt = std::max(pt1, pt2);
        float subleadPt = std::min(pt1, pt2);

        if (leadPt < 50.0) continue;
        if (subleadPt < 20.0) continue;

        index = (int)j;
        break;
      }

      if (index < 0) continue;

      usedPrompt[i] = true;
      usedPrompt[(size_t)index] = true;

      const auto& gen1 = *(promptMuons.at(i));
      const auto& gen2 = *(promptMuons.at(index));
    
      int idx_any1 = -1;
      int idx_any2 = -1;
      int idx_any1_woID = -1;
      int idx_any2_woID = -1;
      int idx_G = -1;
    
      
      int sub_muon_high_pt_flag = 0;
      bool add1TrkFlag = false;
      float add1TrkPt = -999.f;
      float add1TrkEta = 999.f;
      float add1TrkPhi = 999.f;
      bool add2TrkFlag = false;
      float add2TrkPt = -999.f;
      float add2TrkEta = 999.f;
      float add2TrkPhi = 999.f;
      
    
      // gen1/gen2에 대한 best / second-best reco match 저장
      int best_idx1 = -1, second_idx1 = -1;
      int best_idx2 = -1, second_idx2 = -1;
      float best_dr1 = 999.f, second_dr1 = 999.f;
      float best_dr2 = 999.f, second_dr2 = 999.f;
    
      for (size_t iMuon = 0; iMuon < muonHandle->size(); ++iMuon) {
        const reco::Muon& muon = (*muonHandle)[iMuon];
        // 기존 counting/flag용 후보는 gen1 주변 reco만 보던 구조를 유지
        if (reco::deltaR(muon, gen1) > 0.1) continue;
    
        if (muon::isHighPtMuon(muon, *pv)) sub_muon_high_pt_flag++;
    
        bool T = muon::isTrackerHighPtMuon(muon, *pv);
        bool G = muon::isHighPtMuon(muon, *pv);
    
        // woID reco count
        if (idx_any1_woID < 0) idx_any1_woID = (int)iMuon;
        else if ((int)iMuon != idx_any1_woID && idx_any2_woID < 0) idx_any2_woID = (int)iMuon;
    
        // ID reco count
        if (T || G) {
          if (G && idx_G < 0) idx_G = (int)iMuon;
    
          if (idx_any1 < 0) idx_any1 = (int)iMuon;
          else if ((int)iMuon != idx_any1 && idx_any2 < 0) idx_any2 = (int)iMuon;
        }
    
        // ---------- gen1 matching candidate ----------
        float dr1  = reco::deltaR(muon, gen1);
        float dpt1 = std::abs(muon.pt() - gen1.pt()) / gen1.pt();
        if (dr1 < 0.03 && dpt1 < 0.30 && muon.charge() == gen1.charge()) {
          if (dr1 < best_dr1) {
            second_dr1 = best_dr1;
            second_idx1 = best_idx1;
            best_dr1 = dr1;
            best_idx1 = (int)iMuon;
          } else if (dr1 < second_dr1) {
            second_dr1 = dr1;
            second_idx1 = (int)iMuon;
          }
        }

        // ---------- gen2 matching candidate ----------
        float dr2  = reco::deltaR(muon, gen2);
        float dpt2 = std::abs(muon.pt() - gen2.pt()) / gen2.pt();
        if (dr2 < 0.03 && dpt2 < 0.30 && muon.charge() == gen2.charge()) {
          if (dr2 < best_dr2) {
            second_dr2 = best_dr2;
            second_idx2 = best_idx2;
            best_dr2 = dr2;
            best_idx2 = (int)iMuon;
          } else if (dr2 < second_dr2) {
            second_dr2 = dr2;
            second_idx2 = (int)iMuon;
          }
        }
      }
    
      // 최종 unique matching
      int match_idx1 = best_idx1;
      int match_idx2 = best_idx2;
    
      // 둘이 같은 reco muon을 잡은 경우 처리
      if (match_idx1 >= 0 && match_idx2 >= 0 && match_idx1 == match_idx2) {
        if (best_dr1 <= best_dr2) {
          // gen1이 그 reco를 유지, gen2는 second-best 시도
          match_idx2 = second_idx2;
          if (match_idx2 == match_idx1) match_idx2 = -1;
        } else {
          // gen2가 그 reco를 유지, gen1은 second-best 시도
          match_idx1 = second_idx1;
          if (match_idx1 == match_idx2) match_idx1 = -1;
        }
      }
    
      // 혹시 second-best도 같은 경우 방지
      if (match_idx1 >= 0 && match_idx2 >= 0 && match_idx1 == match_idx2) {
        match_idx2 = -1;
      }
    
      // tree 변수 채우기 전 초기화
      eff_gen_1_pt  = gen1.pt();
      eff_gen_1_phi = gen1.phi();
      eff_gen_1_eta = gen1.eta();
    
      eff_gen_2_pt  = gen2.pt();
      eff_gen_2_phi = gen2.phi();
      eff_gen_2_eta = gen2.eta();
    
      eff_reco_1_idx = -1;
      eff_reco_2_idx = -1;
    
      eff_reco_1_pt  = -999.f;
      eff_reco_1_phi = -999.f;
      eff_reco_1_eta = -999.f;
      eff_reco_1_dr  = -999.f;
      eff_reco_1_highptid = 0;
      eff_reco_1_trackerhighptid = 0;
      eff_reco_1_em = -999.f;
      eff_reco_1_emS9 = -999.f;
      eff_reco_1_emS25 = -999.f;
      eff_reco_1_emMax = -999.f;
      eff_reco_1_had = -999.f;
      eff_reco_1_hadS9 = -999.f;
      eff_reco_1_hadMax = -999.f;
      eff_reco_1_ho = -999.f;
      eff_reco_1_hoS9 = -999.f;
      eff_reco_1_isGlobal = false;
      eff_reco_1_isTracker = false;
      eff_reco_1_isPFMuon = false;
      eff_reco_1_isPFIsolationValid = false;
      eff_reco_1_numOfMatchedStations = -999;
      eff_reco_1_numOfChambers = -999;
      eff_reco_1_segCompatibility = -999.;
      eff_reco_1_caloCompatibility = -999.;
      eff_reco_1_chi2LocalPosition = -999.f;
      eff_reco_1_trkKink = -999.f;
      eff_reco_1_glbKink = -999.f;
      eff_reco_1_pfrelIso04 = -999.f;
      eff_reco_1_isLoose = false;
      eff_reco_1_isMedium = false;
      eff_reco_1_isTight = false;
      eff_reco_1_nValidMuonHits = -999;
      eff_reco_1_normChi2_global = -999.f;
      eff_reco_1_nValidPixelHits = -999;
      eff_reco_1_trkLayers = -999;
      eff_reco_1_pixelLayers = -999;
      eff_reco_1_normChi2_inner = -999.f;
      eff_reco_1_validFraction = -999.f;
      eff_reco_1_tunePtErrorOverPt = -999.f;
      eff_reco_1_dxy_PV = -999.f;
      eff_reco_1_dz_PV = -999.f;
      eff_reco_1_nMatches = -999;
      eff_reco_1_innerTrk_pt = -999.f;
      eff_reco_1_innerTrk_eta = -999.f;
      eff_reco_1_innerTrk_phi = -999.f;
      eff_reco_1_outerTrk_pt = -999.f;
      eff_reco_1_outerTrk_eta = -999.f;
      eff_reco_1_outerTrk_phi = -999.f;
      eff_reco_1_globalTrk_pt = -999.f;
      eff_reco_1_globalTrk_eta = -999.f;
      eff_reco_1_globalTrk_phi = -999.f;
      eff_reco_1_bestTrk_pt = -999.f;
      eff_reco_1_bestTrk_eta = -999.f;
      eff_reco_1_bestTrk_phi = -999.f;
      eff_reco_1_tunePTrk_pt = -999.f;
      eff_reco_1_tunePTrk_eta = -999.f;
      eff_reco_1_tunePTrk_phi = -999.f;

      eff_reco_2_pt  = -999.f;
      eff_reco_2_phi = -999.f;
      eff_reco_2_eta = -999.f;
      eff_reco_2_dr  = -999.f;
      eff_reco_2_highptid = 0;
      eff_reco_2_trackerhighptid = 0;
      eff_reco_2_em = -999.f;
      eff_reco_2_emS9 = -999.f;
      eff_reco_2_emS25 = -999.f;
      eff_reco_2_emMax = -999.f;
      eff_reco_2_had = -999.f;
      eff_reco_2_hadS9 = -999.f;
      eff_reco_2_hadMax = -999.f;
      eff_reco_2_ho = -999.f;
      eff_reco_2_hoS9 = -999.f;
      eff_reco_2_isGlobal = false;
      eff_reco_2_isTracker = false;
      eff_reco_2_isPFMuon = false;
      eff_reco_2_isPFIsolationValid = false;
      eff_reco_2_numOfMatchedStations = -999;
      eff_reco_2_numOfChambers = -999;
      eff_reco_2_segCompatibility = -999.;
      eff_reco_2_caloCompatibility = -999.;
      eff_reco_2_chi2LocalPosition = -999.f;
      eff_reco_2_trkKink = -999.f;
      eff_reco_2_glbKink = -999.f;
      eff_reco_2_pfrelIso04 = -999.f;
      eff_reco_2_isLoose = false;
      eff_reco_2_isMedium = false;
      eff_reco_2_isTight = false;
      eff_reco_2_nValidMuonHits = -999;
      eff_reco_2_normChi2_global = -999.f;
      eff_reco_2_nValidPixelHits = -999;
      eff_reco_2_trkLayers = -999;
      eff_reco_2_pixelLayers = -999;
      eff_reco_2_normChi2_inner = -999.f;
      eff_reco_2_validFraction = -999.f;
      eff_reco_2_tunePtErrorOverPt = -999.f;
      eff_reco_2_dxy_PV = -999.f;
      eff_reco_2_dz_PV = -999.f;
      eff_reco_2_nMatches = -999;
      eff_reco_2_innerTrk_pt = -999.f;
      eff_reco_2_innerTrk_eta = -999.f;
      eff_reco_2_innerTrk_phi = -999.f;
      eff_reco_2_outerTrk_pt = -999.f;
      eff_reco_2_outerTrk_eta = -999.f;
      eff_reco_2_outerTrk_phi = -999.f;
      eff_reco_2_globalTrk_pt = -999.f;
      eff_reco_2_globalTrk_eta = -999.f;
      eff_reco_2_globalTrk_phi = -999.f;
      eff_reco_2_bestTrk_pt = -999.f;
      eff_reco_2_bestTrk_eta = -999.f;
      eff_reco_2_bestTrk_phi = -999.f;
      eff_reco_2_tunePTrk_pt = -999.f;
      eff_reco_2_tunePTrk_eta = -999.f;
      eff_reco_2_tunePTrk_phi = -999.f;

      if (match_idx1 >= 0) {
        eff_reco_1_pt  = (*muonHandle)[match_idx1].pt();
        eff_reco_1_phi = (*muonHandle)[match_idx1].phi();
        eff_reco_1_eta = (*muonHandle)[match_idx1].eta();
        eff_reco_1_dr  = reco::deltaR((*muonHandle)[match_idx1], gen1);
        eff_reco_1_idx = match_idx1;
        eff_reco_1_highptid =
            muon::isHighPtMuon((*muonHandle)[match_idx1], *pv);
        eff_reco_1_trackerhighptid =
        muon::isTrackerHighPtMuon((*muonHandle)[match_idx1], *pv);
        eff_reco_1_em = (*muonHandle)[match_idx1].calEnergy().em;
        eff_reco_1_emS9 = (*muonHandle)[match_idx1].calEnergy().emS9;
        eff_reco_1_emS25 = (*muonHandle)[match_idx1].calEnergy().emS25;
        eff_reco_1_emMax = (*muonHandle)[match_idx1].calEnergy().emMax;
        eff_reco_1_had = (*muonHandle)[match_idx1].calEnergy().had;
        eff_reco_1_hadS9 = (*muonHandle)[match_idx1].calEnergy().hadS9;
        eff_reco_1_hadMax = (*muonHandle)[match_idx1].calEnergy().hadMax;
        eff_reco_1_ho = (*muonHandle)[match_idx1].calEnergy().ho;
        eff_reco_1_hoS9 = (*muonHandle)[match_idx1].calEnergy().hoS9;
        eff_reco_1_isGlobal = (*muonHandle)[match_idx1].isGlobalMuon();
        eff_reco_1_isTracker = (*muonHandle)[match_idx1].isTrackerMuon();
        eff_reco_1_isPFMuon = (*muonHandle)[match_idx1].isPFMuon();
        eff_reco_1_isPFIsolationValid = (*muonHandle)[match_idx1].isPFIsolationValid();
        eff_reco_1_numOfMatchedStations = (*muonHandle)[match_idx1].numberOfMatchedStations();
        eff_reco_1_numOfChambers = (*muonHandle)[match_idx1].numberOfChambers();
        eff_reco_1_segCompatibility = muon::segmentCompatibility((*muonHandle)[match_idx1]);
        eff_reco_1_caloCompatibility = (*muonHandle)[match_idx1].caloCompatibility();
        eff_reco_1_chi2LocalPosition = (*muonHandle)[match_idx1].combinedQuality().chi2LocalPosition;
        eff_reco_1_trkKink = (*muonHandle)[match_idx1].combinedQuality().trkKink;
        eff_reco_1_glbKink = (*muonHandle)[match_idx1].combinedQuality().glbKink;
	if (eff_reco_1_isPFMuon && eff_reco_1_isPFIsolationValid){
      	  auto iso04 = (*muonHandle)[match_idx1].pfIsolationR04();
	  eff_reco_1_pfrelIso04 =
            ( iso04.sumChargedHadronPt
            + std::max(0.f,
                  iso04.sumNeutralHadronEt
                + iso04.sumPhotonEt
                - 0.5f * iso04.sumPUPt ) )
            / (*muonHandle)[match_idx1].pt();
	}
        eff_reco_1_isLoose  = muon::isLooseMuon((*muonHandle)[match_idx1]);
        eff_reco_1_isMedium = muon::isMediumMuon((*muonHandle)[match_idx1]);
        eff_reco_1_isTight  = muon::isTightMuon((*muonHandle)[match_idx1], *pv);
        const auto& glbTrk1 = (*muonHandle)[match_idx1].globalTrack();
        if (glbTrk1.isNonnull()) {
          eff_reco_1_nValidMuonHits  = glbTrk1->hitPattern().numberOfValidMuonHits();
          eff_reco_1_normChi2_global = glbTrk1->normalizedChi2();
        }
        const auto& innTrk1 = (*muonHandle)[match_idx1].innerTrack();
        if (innTrk1.isNonnull()) {
          eff_reco_1_nValidPixelHits = innTrk1->hitPattern().numberOfValidPixelHits();
          eff_reco_1_trkLayers       = innTrk1->hitPattern().trackerLayersWithMeasurement();
          eff_reco_1_pixelLayers     = innTrk1->hitPattern().pixelLayersWithMeasurement();
          eff_reco_1_normChi2_inner  = innTrk1->normalizedChi2();
          eff_reco_1_validFraction   = innTrk1->validFraction();
          eff_reco_1_dxy_PV          = innTrk1->dxy(pv->position());
          eff_reco_1_dz_PV           = innTrk1->dz(pv->position());
        }
        const auto& tunePTrk1 = (*muonHandle)[match_idx1].tunePMuonBestTrack();
        if (tunePTrk1.isNonnull() && tunePTrk1->pt() > 0.f)
          eff_reco_1_tunePtErrorOverPt = tunePTrk1->ptError() / tunePTrk1->pt();
        eff_reco_1_nMatches = (*muonHandle)[match_idx1].numberOfMatches(reco::Muon::SegmentAndTrackArbitration);
        if (innTrk1.isNonnull()) {
          eff_reco_1_innerTrk_pt  = innTrk1->pt();
          eff_reco_1_innerTrk_eta = innTrk1->eta();
          eff_reco_1_innerTrk_phi = innTrk1->phi();
        }
        const auto& outTrk1 = (*muonHandle)[match_idx1].outerTrack();
        if (outTrk1.isNonnull()) {
          eff_reco_1_outerTrk_pt  = outTrk1->pt();
          eff_reco_1_outerTrk_eta = outTrk1->eta();
          eff_reco_1_outerTrk_phi = outTrk1->phi();
        }
        if (glbTrk1.isNonnull()) {
          eff_reco_1_globalTrk_pt  = glbTrk1->pt();
          eff_reco_1_globalTrk_eta = glbTrk1->eta();
          eff_reco_1_globalTrk_phi = glbTrk1->phi();
        }
        const auto& bestTrk1 = (*muonHandle)[match_idx1].muonBestTrack();
        if (bestTrk1.isNonnull()) {
          eff_reco_1_bestTrk_pt  = bestTrk1->pt();
          eff_reco_1_bestTrk_eta = bestTrk1->eta();
          eff_reco_1_bestTrk_phi = bestTrk1->phi();
        }
        if (tunePTrk1.isNonnull()) {
          eff_reco_1_tunePTrk_pt  = tunePTrk1->pt();
          eff_reco_1_tunePTrk_eta = tunePTrk1->eta();
          eff_reco_1_tunePTrk_phi = tunePTrk1->phi();
        }
        const reco::Muon& muon = (*muonHandle)[match_idx1];
        const auto addPackedCand = muonTkIsoCalc_.additionalPackedCandSelector(muon, trackCandsHandles, trackCandsVetos_, ttBuilder);
        const reco::Track* addPackedBestTrk = addPackedCand.isNonnull() ? addPackedCand->bestTrack() : nullptr;
        if (addPackedCand.isNonnull()) {
	  add1TrkFlag = true;
	  add1TrkPt = addPackedBestTrk->pt();
	  add1TrkEta = addPackedBestTrk->eta();
	  add1TrkPhi = addPackedBestTrk->phi();
	}
      }
    
      if (match_idx2 >= 0) {
        eff_reco_2_pt  = (*muonHandle)[match_idx2].pt();
        eff_reco_2_phi = (*muonHandle)[match_idx2].phi();
        eff_reco_2_eta = (*muonHandle)[match_idx2].eta();
        eff_reco_2_dr  = reco::deltaR((*muonHandle)[match_idx2], gen2);
        eff_reco_2_idx = match_idx2;
        eff_reco_2_highptid =
            muon::isHighPtMuon((*muonHandle)[match_idx2], *pv);
        eff_reco_2_trackerhighptid =
            muon::isTrackerHighPtMuon((*muonHandle)[match_idx2], *pv);
        eff_reco_2_em = (*muonHandle)[match_idx2].calEnergy().em;
        eff_reco_2_emS9 = (*muonHandle)[match_idx2].calEnergy().emS9;
        eff_reco_2_emS25 = (*muonHandle)[match_idx2].calEnergy().emS25;
        eff_reco_2_emMax = (*muonHandle)[match_idx2].calEnergy().emMax;
        eff_reco_2_had = (*muonHandle)[match_idx2].calEnergy().had;
        eff_reco_2_hadS9 = (*muonHandle)[match_idx2].calEnergy().hadS9;
        eff_reco_2_hadMax = (*muonHandle)[match_idx2].calEnergy().hadMax;
        eff_reco_2_ho = (*muonHandle)[match_idx2].calEnergy().ho;
        eff_reco_2_hoS9 = (*muonHandle)[match_idx2].calEnergy().hoS9;
        eff_reco_2_isGlobal = (*muonHandle)[match_idx2].isGlobalMuon();
        eff_reco_2_isTracker = (*muonHandle)[match_idx2].isTrackerMuon();
        eff_reco_2_isPFMuon = (*muonHandle)[match_idx2].isPFMuon();
        eff_reco_2_isPFIsolationValid = (*muonHandle)[match_idx2].isPFIsolationValid();
        eff_reco_2_numOfMatchedStations = (*muonHandle)[match_idx2].numberOfMatchedStations();
        eff_reco_2_numOfChambers = (*muonHandle)[match_idx2].numberOfChambers();
        eff_reco_2_segCompatibility = muon::segmentCompatibility((*muonHandle)[match_idx2]);
        eff_reco_2_caloCompatibility = (*muonHandle)[match_idx2].caloCompatibility();
        eff_reco_2_chi2LocalPosition = (*muonHandle)[match_idx2].combinedQuality().chi2LocalPosition;
        eff_reco_2_trkKink = (*muonHandle)[match_idx2].combinedQuality().trkKink;
        eff_reco_2_glbKink = (*muonHandle)[match_idx2].combinedQuality().glbKink;
	if (eff_reco_2_isPFMuon && eff_reco_2_isPFIsolationValid){
      	  auto iso04 = (*muonHandle)[match_idx2].pfIsolationR04();
	  eff_reco_2_pfrelIso04 =
            ( iso04.sumChargedHadronPt
            + std::max(0.f,
                  iso04.sumNeutralHadronEt
                + iso04.sumPhotonEt
                - 0.5f * iso04.sumPUPt ) )
            / (*muonHandle)[match_idx2].pt();
	}
        eff_reco_2_isLoose  = muon::isLooseMuon((*muonHandle)[match_idx2]);
        eff_reco_2_isMedium = muon::isMediumMuon((*muonHandle)[match_idx2]);
        eff_reco_2_isTight  = muon::isTightMuon((*muonHandle)[match_idx2], *pv);
        const auto& glbTrk2 = (*muonHandle)[match_idx2].globalTrack();
        if (glbTrk2.isNonnull()) {
          eff_reco_2_nValidMuonHits  = glbTrk2->hitPattern().numberOfValidMuonHits();
          eff_reco_2_normChi2_global = glbTrk2->normalizedChi2();
        }
        const auto& innTrk2 = (*muonHandle)[match_idx2].innerTrack();
        if (innTrk2.isNonnull()) {
          eff_reco_2_nValidPixelHits = innTrk2->hitPattern().numberOfValidPixelHits();
          eff_reco_2_trkLayers       = innTrk2->hitPattern().trackerLayersWithMeasurement();
          eff_reco_2_pixelLayers     = innTrk2->hitPattern().pixelLayersWithMeasurement();
          eff_reco_2_normChi2_inner  = innTrk2->normalizedChi2();
          eff_reco_2_validFraction   = innTrk2->validFraction();
          eff_reco_2_dxy_PV          = innTrk2->dxy(pv->position());
          eff_reco_2_dz_PV           = innTrk2->dz(pv->position());
        }
        const auto& tunePTrk2 = (*muonHandle)[match_idx2].tunePMuonBestTrack();
        if (tunePTrk2.isNonnull() && tunePTrk2->pt() > 0.f)
          eff_reco_2_tunePtErrorOverPt = tunePTrk2->ptError() / tunePTrk2->pt();
        eff_reco_2_nMatches = (*muonHandle)[match_idx2].numberOfMatches(reco::Muon::SegmentAndTrackArbitration);
        if (innTrk2.isNonnull()) {
          eff_reco_2_innerTrk_pt  = innTrk2->pt();
          eff_reco_2_innerTrk_eta = innTrk2->eta();
          eff_reco_2_innerTrk_phi = innTrk2->phi();
        }
        const auto& outTrk2 = (*muonHandle)[match_idx2].outerTrack();
        if (outTrk2.isNonnull()) {
          eff_reco_2_outerTrk_pt  = outTrk2->pt();
          eff_reco_2_outerTrk_eta = outTrk2->eta();
          eff_reco_2_outerTrk_phi = outTrk2->phi();
        }
        if (glbTrk2.isNonnull()) {
          eff_reco_2_globalTrk_pt  = glbTrk2->pt();
          eff_reco_2_globalTrk_eta = glbTrk2->eta();
          eff_reco_2_globalTrk_phi = glbTrk2->phi();
        }
        const auto& bestTrk2 = (*muonHandle)[match_idx2].muonBestTrack();
        if (bestTrk2.isNonnull()) {
          eff_reco_2_bestTrk_pt  = bestTrk2->pt();
          eff_reco_2_bestTrk_eta = bestTrk2->eta();
          eff_reco_2_bestTrk_phi = bestTrk2->phi();
        }
        if (tunePTrk2.isNonnull()) {
          eff_reco_2_tunePTrk_pt  = tunePTrk2->pt();
          eff_reco_2_tunePTrk_eta = tunePTrk2->eta();
          eff_reco_2_tunePTrk_phi = tunePTrk2->phi();
        }
        const reco::Muon& muon = (*muonHandle)[match_idx2];
        const auto addPackedCand = muonTkIsoCalc_.additionalPackedCandSelector(muon, trackCandsHandles, trackCandsVetos_, ttBuilder);
        const reco::Track* addPackedBestTrk = addPackedCand.isNonnull() ? addPackedCand->bestTrack() : nullptr;
        if (addPackedCand.isNonnull()) {
	  add2TrkFlag = true;
	  add2TrkPt = addPackedBestTrk->pt();
	  add2TrkEta = addPackedBestTrk->eta();
	  add2TrkPhi = addPackedBestTrk->phi();
	}
      }

      add1_track_flag = add1TrkFlag;
      add1_track_pt = add1TrkPt;
      add1_track_eta = add1TrkEta;
      add1_track_phi = add1TrkPhi;
      add2_track_flag = add2TrkFlag;
      add2_track_pt = add2TrkPt;
      add2_track_eta = add2TrkEta;
      add2_track_phi = add2TrkPhi;
    
      dR_gen = reco::deltaR(gen1, gen2);
      num_reco_muon = sub_muon_high_pt_flag;
      flag_Id = (idx_G >= 0) && (idx_any2 >= 0);
      flag_Id_woID = (idx_any2_woID >= 0);
      flag_Id_any = (idx_any2 >= 0);
      eff_weight = mcweight;
      const auto& aMET = metHandle->at(0);
      eff_MET_cor_XY_phi = aMET.corPhi(pat::MET::TypeXY);
      eff_MET_phi = aMET.phi();
      eff_MET = aMET.pt();
      eff_MET_cor = aMET.corPt(pat::MET::Type1);
      muonEfficiencyTree_->Fill();
    }
    for (size_t i = 0; i < muonHandle->size(); ++i) {
      const reco::Muon& muon = (*muonHandle)[i];

      const reco::TrackRef muTrkRef = muon.muonBestTrack();

      //std::cout<<"hello "<<i<<<" | "<<addPackedCand.isNonnull()<std::endl;

      int close_muon = 0;
      for (size_t j = 0; j < muonHandle->size(); ++j) {
        const reco::Muon& muon_sub = (*muonHandle)[j];
        double dR = reco::deltaR(muon, muon_sub);
	//double dR = sqrt(pow(muon_sub.eta()-muon.eta(),2)+pow(muon_sub.phi()-muon.phi(),2));
	//if (muon::isHighPtMuon(muon_sub,*pv) && (&muon != &muon_sub) && dR < 0.1) sub_muon_high_pt_flag ++;
	if (!muon::isLooseMuon(muon_sub)) continue;
        if (dR < 0.1 && i != j ) close_muon ++;
      }
      int gen_muon = 0;
      std::vector<reco::GenParticleRef> target_muon;
      for (auto gen : promptMuons){
        double dR = reco::deltaR(*gen, muon);
        if (dR < 0.1){
          gen_muon ++;
	  target_muon.push_back(gen);
        }      
      }
      if (!muon::isHighPtMuon(muon,*pv) ) continue;
      // if (gen_muon == 2){
      //    dR_gen = reco::deltaR(*(target_muon.at(0)), *(target_muon.at(1)));
      //    num_reco_muon = sub_muon_high_pt_flag;
      //    muonEfficiencyTree_->Fill();
      // }
      //std::cout<<muon.isEnergyValid()<<" | "<<muon.calEnergy().towerS9<<" | "<<muon.calEnergy().emS25<<" | "<<muon.calEnergy().hadS9<<" | "<<muon.calEnergy().hoS9<<" | "<<muon.pt()<<" | "<<gen_muon<<" | "<<muon.eta()<<" | "<<muon.phi()<<" | "<<close_muon<<std::endl;
      if (gen_muon == 1 && muon.isPFMuon() && muon.isPFIsolationValid()){
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
          isGlobal_muon = muon.isGlobalMuon();
          isTracker_muon = muon.isTrackerMuon();
          numOfMatchedStations_muon = muon.numberOfMatchedStations();
          numOfChambers_muon = muon.numberOfChambers();
          segCompatibility_muon = muon::segmentCompatibility(muon);
          caloCompatibility_muon = muon.caloCompatibility();
          chi2LocalPosition_muon = muon.combinedQuality().chi2LocalPosition;
          trkKink_muon = muon.combinedQuality().trkKink;
          glbKink_muon = muon.combinedQuality().glbKink;
	  auto iso04 = muon.pfIsolationR04();
	  pfrelIso04_muon =
            ( iso04.sumChargedHadronPt
            + std::max(0.f,
                  iso04.sumNeutralHadronEt
                + iso04.sumPhotonEt
                - 0.5f * iso04.sumPUPt ) )
            / muon.pt();
	  gen_eta_muon = target_muon.at(0)->eta();
	  gen_phi_muon = target_muon.at(0)->phi();
	  gen_pt_muon = target_muon.at(0)->pt();
	  weight_muon = mcweight;
          muon_->Fill();
      } 
      if (gen_muon == 2 && close_muon == 0){
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
          isGlobal_mergedMuon = muon.isGlobalMuon();
          isTracker_mergedMuon = muon.isTrackerMuon();
          numOfMatchedStations_mergedMuon = muon.numberOfMatchedStations();
          numOfChambers_mergedMuon = muon.numberOfChambers();
          segCompatibility_mergedMuon = muon::segmentCompatibility(muon);
          caloCompatibility_mergedMuon = muon.caloCompatibility();
          chi2LocalPosition_mergedMuon = muon.combinedQuality().chi2LocalPosition;
          trkKink_mergedMuon = muon.combinedQuality().trkKink;
          glbKink_mergedMuon = muon.combinedQuality().glbKink;
	  gen_eta_mergedMuon = target_muon.at(0)->eta();
	  gen_phi_mergedMuon = target_muon.at(0)->phi();
	  gen_pt_mergedMuon = target_muon.at(0)->pt();
	  gen_sub_eta_mergedMuon = target_muon.at(1)->eta();
	  gen_sub_phi_mergedMuon = target_muon.at(1)->phi();
	  gen_sub_pt_mergedMuon = target_muon.at(1)->pt();
	  weight_mergedMuon = mcweight;
          mergedMuon_->Fill();
      } 
      target_muon.clear();
    
    }
  }
  return;
}

DEFINE_FWK_MODULE(MergedMuon);
