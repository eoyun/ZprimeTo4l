#include <memory>
#include <iostream>


#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedDEtaInSeed.h"

#include "DataFormats/PatCandidates/interface/MET.h"

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

class MuonJpsi : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MuonJpsi(const edm::ParameterSet&);
  virtual ~MuonJpsi() {}
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
  const double mumass_ = 0.1056583745;
  const double mumassErr_ = 0.0000000024;
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

  TTree* jpsi_ = nullptr;

  int runNo_jpsi;
  int lumiNo_jpsi;
  int evtNo_jpsi;
  int weight_mu_jpsi;
  std::vector<float> invM_jpsi;
  std::vector<float> pT_mu_sub_jpsi;
  std::vector<float> eta_mu_sub_jpsi;
  std::vector<float> phi_mu_sub_jpsi;
  std::vector<int> idx_mu_sub_jpsi;
  std::vector<bool> looseID_mu_sub_jpsi;
  std::vector<bool> highptID_mu_sub_jpsi;
  std::vector<bool> trackerhighptID_mu_sub_jpsi;
  std::vector<bool> mediumID_mu_sub_jpsi;
  std::vector<bool> tightID_mu_sub_jpsi;
  std::vector<bool> isGlobal_mu_sub_jpsi;
  std::vector<bool> isTracker_mu_sub_jpsi;
  std::vector<bool> isPFMuon_mu_sub_jpsi;
  std::vector<bool> isPFIsolationValid_mu_sub_jpsi;
  std::vector<float> em_mu_sub_jpsi;
  std::vector<float> emS9_mu_sub_jpsi;
  std::vector<float> emS25_mu_sub_jpsi;
  std::vector<float> emMax_mu_sub_jpsi;
  std::vector<float> had_mu_sub_jpsi;
  std::vector<float> hadS9_mu_sub_jpsi;
  std::vector<float> hadMax_mu_sub_jpsi;
  std::vector<float> ho_mu_sub_jpsi;
  std::vector<float> hoS9_mu_sub_jpsi;
  std::vector<int> numOfMatchedStations_mu_sub_jpsi;
  std::vector<int> numOfChambers_mu_sub_jpsi;
  std::vector<double> segCompatibility_mu_sub_jpsi;
  std::vector<double> caloCompatibility_mu_sub_jpsi;
  std::vector<float> chi2LocalPosition_mu_sub_jpsi;
  std::vector<float> trkKink_mu_sub_jpsi;
  std::vector<float> glbKink_mu_sub_jpsi;
  std::vector<float> pfrelIso04_mu_sub_jpsi;
  std::vector<int> nValidMuonHits_mu_sub_jpsi;
  std::vector<float> normChi2_global_mu_sub_jpsi;
  std::vector<int> nValidPixelHits_mu_sub_jpsi;
  std::vector<int> trkLayers_mu_sub_jpsi;
  std::vector<int> pixelLayers_mu_sub_jpsi;
  std::vector<float> normChi2_inner_mu_sub_jpsi;
  std::vector<float> validFraction_mu_sub_jpsi;
  std::vector<float> tunePtErrorOverPt_mu_sub_jpsi;
  std::vector<float> dxy_PV_mu_sub_jpsi;
  std::vector<float> dz_PV_mu_sub_jpsi;
  std::vector<int> nMatches_mu_sub_jpsi;
  std::vector<float> innerTrk_pt_mu_sub_jpsi;
  std::vector<float> innerTrk_eta_mu_sub_jpsi;
  std::vector<float> innerTrk_phi_mu_sub_jpsi;
  std::vector<float> outerTrk_pt_mu_sub_jpsi;
  std::vector<float> outerTrk_eta_mu_sub_jpsi;
  std::vector<float> outerTrk_phi_mu_sub_jpsi;
  std::vector<float> globalTrk_pt_mu_sub_jpsi;
  std::vector<float> globalTrk_eta_mu_sub_jpsi;
  std::vector<float> globalTrk_phi_mu_sub_jpsi;
  std::vector<float> bestTrk_pt_mu_sub_jpsi;
  std::vector<float> bestTrk_eta_mu_sub_jpsi;
  std::vector<float> bestTrk_phi_mu_sub_jpsi;
  std::vector<float> tunePTrk_pt_mu_sub_jpsi;
  std::vector<float> tunePTrk_eta_mu_sub_jpsi;
  std::vector<float> tunePTrk_phi_mu_sub_jpsi;
  std::vector<int>   muon_charge_mu_sub_jpsi;
  std::vector<int>   innerTrk_charge_mu_sub_jpsi;
  std::vector<int>   outerTrk_charge_mu_sub_jpsi;
  std::vector<int>   globalTrk_charge_mu_sub_jpsi;
  std::vector<int>   bestTrk_charge_mu_sub_jpsi;
  std::vector<int>   tunePTrk_charge_mu_sub_jpsi;
  std::vector<float> innerTrk_ptError_mu_sub_jpsi;
  std::vector<float> outerTrk_ptError_mu_sub_jpsi;
  std::vector<float> globalTrk_ptError_mu_sub_jpsi;
  std::vector<float> bestTrk_ptError_mu_sub_jpsi;
  std::vector<float> tunePTrk_ptError_mu_sub_jpsi;
  std::vector<float> dxy_bestTrk_PV_mu_sub_jpsi;
  std::vector<float> dz_bestTrk_PV_mu_sub_jpsi;
  std::vector<float> pT_mu_jpsi;
  std::vector<float> eta_mu_jpsi;
  std::vector<float> phi_mu_jpsi;
  std::vector<int> idx_mu_jpsi;
  std::vector<bool> looseID_mu_jpsi;
  std::vector<bool> highptID_mu_jpsi;
  std::vector<bool> trackerhighptID_mu_jpsi;
  std::vector<bool> mediumID_mu_jpsi;
  std::vector<bool> tightID_mu_jpsi;
  std::vector<bool> isGlobal_mu_jpsi;
  std::vector<bool> isTracker_mu_jpsi;
  std::vector<bool> isPFMuon_mu_jpsi;
  std::vector<bool> isPFIsolationValid_mu_jpsi;
  std::vector<float> em_mu_jpsi;
  std::vector<float> emS9_mu_jpsi;
  std::vector<float> emS25_mu_jpsi;
  std::vector<float> emMax_mu_jpsi;
  std::vector<float> had_mu_jpsi;
  std::vector<float> hadS9_mu_jpsi;
  std::vector<float> hadMax_mu_jpsi;
  std::vector<float> ho_mu_jpsi;
  std::vector<float> hoS9_mu_jpsi;
  std::vector<int> numOfMatchedStations_mu_jpsi;
  std::vector<int> numOfChambers_mu_jpsi;
  std::vector<double> segCompatibility_mu_jpsi;
  std::vector<double> caloCompatibility_mu_jpsi;
  std::vector<float> chi2LocalPosition_mu_jpsi;
  std::vector<float> trkKink_mu_jpsi;
  std::vector<float> glbKink_mu_jpsi;
  std::vector<float> pfrelIso04_mu_jpsi;
  std::vector<int> nValidMuonHits_mu_jpsi;
  std::vector<float> normChi2_global_mu_jpsi;
  std::vector<int> nValidPixelHits_mu_jpsi;
  std::vector<int> trkLayers_mu_jpsi;
  std::vector<int> pixelLayers_mu_jpsi;
  std::vector<float> normChi2_inner_mu_jpsi;
  std::vector<float> validFraction_mu_jpsi;
  std::vector<float> tunePtErrorOverPt_mu_jpsi;
  std::vector<float> dxy_PV_mu_jpsi;
  std::vector<float> dz_PV_mu_jpsi;
  std::vector<int> nMatches_mu_jpsi;
  std::vector<float> innerTrk_pt_mu_jpsi;
  std::vector<float> innerTrk_eta_mu_jpsi;
  std::vector<float> innerTrk_phi_mu_jpsi;
  std::vector<float> outerTrk_pt_mu_jpsi;
  std::vector<float> outerTrk_eta_mu_jpsi;
  std::vector<float> outerTrk_phi_mu_jpsi;
  std::vector<float> globalTrk_pt_mu_jpsi;
  std::vector<float> globalTrk_eta_mu_jpsi;
  std::vector<float> globalTrk_phi_mu_jpsi;
  std::vector<float> bestTrk_pt_mu_jpsi;
  std::vector<float> bestTrk_eta_mu_jpsi;
  std::vector<float> bestTrk_phi_mu_jpsi;
  std::vector<float> tunePTrk_pt_mu_jpsi;
  std::vector<float> tunePTrk_eta_mu_jpsi;
  std::vector<float> tunePTrk_phi_mu_jpsi;
  std::vector<int>   muon_charge_mu_jpsi;
  std::vector<int>   innerTrk_charge_mu_jpsi;
  std::vector<int>   outerTrk_charge_mu_jpsi;
  std::vector<int>   globalTrk_charge_mu_jpsi;
  std::vector<int>   bestTrk_charge_mu_jpsi;
  std::vector<int>   tunePTrk_charge_mu_jpsi;
  std::vector<float> innerTrk_ptError_mu_jpsi;
  std::vector<float> outerTrk_ptError_mu_jpsi;
  std::vector<float> globalTrk_ptError_mu_jpsi;
  std::vector<float> bestTrk_ptError_mu_jpsi;
  std::vector<float> tunePTrk_ptError_mu_jpsi;
  std::vector<float> dxy_bestTrk_PV_mu_jpsi;
  std::vector<float> dz_bestTrk_PV_mu_jpsi;

  TTree* muon_ = nullptr;

  int runNo;
  int lumiNo;
  int evtNo;
  float pT_muoniso;
  float eta_muoniso;
  float phi_muoniso;
  bool mediumID_muoniso;
  float iso04_muoniso;
  float MET;
  float MET_phi;
  float MT;
  float dphi;
  int weight_mu;
  std::vector<float> invM;
  std::vector<float> pT_mu_sub;
  std::vector<float> eta_mu_sub;
  std::vector<float> phi_mu_sub;
  std::vector<int> idx_mu_sub;
  std::vector<bool> looseID_mu_sub;
  std::vector<bool> highptID_mu_sub;
  std::vector<bool> trackerhighptID_mu_sub;
  std::vector<bool> mediumID_mu_sub;
  std::vector<bool> tightID_mu_sub;
  std::vector<bool> isGlobal_mu_sub;
  std::vector<bool> isTracker_mu_sub;
  std::vector<bool> isPFMuon_mu_sub;
  std::vector<bool> isPFIsolationValid_mu_sub;
  std::vector<float> em_mu_sub;
  std::vector<float> emS9_mu_sub;
  std::vector<float> emS25_mu_sub;
  std::vector<float> emMax_mu_sub;
  std::vector<float> had_mu_sub;
  std::vector<float> hadS9_mu_sub;
  std::vector<float> hadMax_mu_sub;
  std::vector<float> ho_mu_sub;
  std::vector<float> hoS9_mu_sub;
  std::vector<int> numOfMatchedStations_mu_sub;
  std::vector<int> numOfChambers_mu_sub;
  std::vector<double> segCompatibility_mu_sub;
  std::vector<double> caloCompatibility_mu_sub;
  std::vector<float> chi2LocalPosition_mu_sub;
  std::vector<float> trkKink_mu_sub;
  std::vector<float> glbKink_mu_sub;
  std::vector<float> pfrelIso04_mu_sub;
  std::vector<int> nValidMuonHits_mu_sub;
  std::vector<float> normChi2_global_mu_sub;
  std::vector<int> nValidPixelHits_mu_sub;
  std::vector<int> trkLayers_mu_sub;
  std::vector<int> pixelLayers_mu_sub;
  std::vector<float> normChi2_inner_mu_sub;
  std::vector<float> validFraction_mu_sub;
  std::vector<float> tunePtErrorOverPt_mu_sub;
  std::vector<float> dxy_PV_mu_sub;
  std::vector<float> dz_PV_mu_sub;
  std::vector<int> nMatches_mu_sub;
  std::vector<float> innerTrk_pt_mu_sub;
  std::vector<float> innerTrk_eta_mu_sub;
  std::vector<float> innerTrk_phi_mu_sub;
  std::vector<float> outerTrk_pt_mu_sub;
  std::vector<float> outerTrk_eta_mu_sub;
  std::vector<float> outerTrk_phi_mu_sub;
  std::vector<float> globalTrk_pt_mu_sub;
  std::vector<float> globalTrk_eta_mu_sub;
  std::vector<float> globalTrk_phi_mu_sub;
  std::vector<float> bestTrk_pt_mu_sub;
  std::vector<float> bestTrk_eta_mu_sub;
  std::vector<float> bestTrk_phi_mu_sub;
  std::vector<float> tunePTrk_pt_mu_sub;
  std::vector<float> tunePTrk_eta_mu_sub;
  std::vector<float> tunePTrk_phi_mu_sub;
  std::vector<int>   muon_charge_mu_sub;
  std::vector<int>   innerTrk_charge_mu_sub;
  std::vector<int>   outerTrk_charge_mu_sub;
  std::vector<int>   globalTrk_charge_mu_sub;
  std::vector<int>   bestTrk_charge_mu_sub;
  std::vector<int>   tunePTrk_charge_mu_sub;
  std::vector<float> innerTrk_ptError_mu_sub;
  std::vector<float> outerTrk_ptError_mu_sub;
  std::vector<float> globalTrk_ptError_mu_sub;
  std::vector<float> bestTrk_ptError_mu_sub;
  std::vector<float> tunePTrk_ptError_mu_sub;
  std::vector<float> dxy_bestTrk_PV_mu_sub;
  std::vector<float> dz_bestTrk_PV_mu_sub;
  std::vector<float> pT_mu;
  std::vector<float> eta_mu;
  std::vector<float> phi_mu;
  std::vector<int> idx_mu;
  std::vector<bool> looseID_mu;
  std::vector<bool> highptID_mu;
  std::vector<bool> trackerhighptID_mu;
  std::vector<bool> mediumID_mu;
  std::vector<bool> tightID_mu;
  std::vector<bool> isGlobal_mu;
  std::vector<bool> isTracker_mu;
  std::vector<bool> isPFMuon_mu;
  std::vector<bool> isPFIsolationValid_mu;
  std::vector<float> em_mu;
  std::vector<float> emS9_mu;
  std::vector<float> emS25_mu;
  std::vector<float> emMax_mu;
  std::vector<float> had_mu;
  std::vector<float> hadS9_mu;
  std::vector<float> hadMax_mu;
  std::vector<float> ho_mu;
  std::vector<float> hoS9_mu;
  std::vector<int> numOfMatchedStations_mu;
  std::vector<int> numOfChambers_mu;
  std::vector<double> segCompatibility_mu;
  std::vector<double> caloCompatibility_mu;
  std::vector<float> chi2LocalPosition_mu;
  std::vector<float> trkKink_mu;
  std::vector<float> glbKink_mu;
  std::vector<float> pfrelIso04_mu;
  std::vector<int> nValidMuonHits_mu;
  std::vector<float> normChi2_global_mu;
  std::vector<int> nValidPixelHits_mu;
  std::vector<int> trkLayers_mu;
  std::vector<int> pixelLayers_mu;
  std::vector<float> normChi2_inner_mu;
  std::vector<float> validFraction_mu;
  std::vector<float> tunePtErrorOverPt_mu;
  std::vector<float> dxy_PV_mu;
  std::vector<float> dz_PV_mu;
  std::vector<int> nMatches_mu;
  std::vector<float> innerTrk_pt_mu;
  std::vector<float> innerTrk_eta_mu;
  std::vector<float> innerTrk_phi_mu;
  std::vector<float> outerTrk_pt_mu;
  std::vector<float> outerTrk_eta_mu;
  std::vector<float> outerTrk_phi_mu;
  std::vector<float> globalTrk_pt_mu;
  std::vector<float> globalTrk_eta_mu;
  std::vector<float> globalTrk_phi_mu;
  std::vector<float> bestTrk_pt_mu;
  std::vector<float> bestTrk_eta_mu;
  std::vector<float> bestTrk_phi_mu;
  std::vector<float> tunePTrk_pt_mu;
  std::vector<float> tunePTrk_eta_mu;
  std::vector<float> tunePTrk_phi_mu;
  std::vector<int>   muon_charge_mu;
  std::vector<int>   innerTrk_charge_mu;
  std::vector<int>   outerTrk_charge_mu;
  std::vector<int>   globalTrk_charge_mu;
  std::vector<int>   bestTrk_charge_mu;
  std::vector<int>   tunePTrk_charge_mu;
  std::vector<float> innerTrk_ptError_mu;
  std::vector<float> outerTrk_ptError_mu;
  std::vector<float> globalTrk_ptError_mu;
  std::vector<float> bestTrk_ptError_mu;
  std::vector<float> tunePTrk_ptError_mu;
  std::vector<float> dxy_bestTrk_PV_mu;
  std::vector<float> dz_bestTrk_PV_mu;




  TTree* jpsiGen_ = nullptr;

  int runNo_gen;
  int lumiNo_gen;
  int evtNo_gen;
  float weight_gen;
  std::vector<float> jpsi_pt_gen;
  std::vector<float> jpsi_eta_gen;
  std::vector<float> jpsi_phi_gen;
  std::vector<float> jpsi_mass_gen;
  std::vector<float> gen_1_pt;
  std::vector<float> gen_1_eta;
  std::vector<float> gen_1_phi;
  std::vector<int>   gen_1_charge;
  std::vector<float> gen_2_pt;
  std::vector<float> gen_2_eta;
  std::vector<float> gen_2_phi;
  std::vector<int>   gen_2_charge;
  std::vector<float> gen_dR;
  std::vector<float> gen_invM;

  PositionCalc posCalcLog_;

};

MuonJpsi::MuonJpsi(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
metToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcMET"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
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

bool MuonJpsi::extrapolate(const reco::GsfElectron& aEle, const reco::TrackBase& addTrk,
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

void MuonJpsi::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",30,0.,30.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",200,-1.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);
  histo1d_["PUsummary"] = fs->make<TH1D>("PUsummary","PUsummary",99,0.,99.);

  jpsi_ = fs->make<TTree>("jpsiTree","jpsiTree");
  jpsi_->Branch("runNo_jpsi",&runNo_jpsi,"runNo_jpsi/I");
  jpsi_->Branch("lumiNo_jpsi",&lumiNo_jpsi,"lumiNo_jpsi/I");
  jpsi_->Branch("evtNo_jpsi",&evtNo_jpsi,"evtNo_jpsi/I");
  jpsi_->Branch("weight_mu_jpsi",&weight_mu_jpsi,"weight_mu_jpsi/I");
  jpsi_->Branch("invM_jpsi",&invM_jpsi,32000,0);
  jpsi_->Branch("pT_mu_sub_jpsi",&pT_mu_sub_jpsi,32000,0);
  jpsi_->Branch("eta_mu_sub_jpsi",&eta_mu_sub_jpsi,32000,0);
  jpsi_->Branch("phi_mu_sub_jpsi",&phi_mu_sub_jpsi,32000,0);
  jpsi_->Branch("idx_mu_sub_jpsi",&idx_mu_sub_jpsi,32000,0);
  jpsi_->Branch("looseID_mu_sub_jpsi",&looseID_mu_sub_jpsi,32000,0);
  jpsi_->Branch("highptID_mu_sub_jpsi",&highptID_mu_sub_jpsi,32000,0);
  jpsi_->Branch("trackerhighptID_mu_sub_jpsi",&trackerhighptID_mu_sub_jpsi,32000,0);
  jpsi_->Branch("mediumID_mu_sub_jpsi",&mediumID_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tightID_mu_sub_jpsi",&tightID_mu_sub_jpsi,32000,0);
  jpsi_->Branch("isGlobal_mu_sub_jpsi",&isGlobal_mu_sub_jpsi,32000,0);
  jpsi_->Branch("isTracker_mu_sub_jpsi",&isTracker_mu_sub_jpsi,32000,0);
  jpsi_->Branch("isPFMuon_mu_sub_jpsi",&isPFMuon_mu_sub_jpsi,32000,0);
  jpsi_->Branch("isPFIsolationValid_mu_sub_jpsi",&isPFIsolationValid_mu_sub_jpsi,32000,0);
  jpsi_->Branch("em_mu_sub_jpsi",&em_mu_sub_jpsi,32000,0);
  jpsi_->Branch("emS9_mu_sub_jpsi",&emS9_mu_sub_jpsi,32000,0);
  jpsi_->Branch("emS25_mu_sub_jpsi",&emS25_mu_sub_jpsi,32000,0);
  jpsi_->Branch("emMax_mu_sub_jpsi",&emMax_mu_sub_jpsi,32000,0);
  jpsi_->Branch("had_mu_sub_jpsi",&had_mu_sub_jpsi,32000,0);
  jpsi_->Branch("hadS9_mu_sub_jpsi",&hadS9_mu_sub_jpsi,32000,0);
  jpsi_->Branch("hadMax_mu_sub_jpsi",&hadMax_mu_sub_jpsi,32000,0);
  jpsi_->Branch("ho_mu_sub_jpsi",&ho_mu_sub_jpsi,32000,0);
  jpsi_->Branch("hoS9_mu_sub_jpsi",&hoS9_mu_sub_jpsi,32000,0);
  jpsi_->Branch("numOfMatchedStations_mu_sub_jpsi",&numOfMatchedStations_mu_sub_jpsi,32000,0);
  jpsi_->Branch("numOfChambers_mu_sub_jpsi",&numOfChambers_mu_sub_jpsi,32000,0);
  jpsi_->Branch("segCompatibility_mu_sub_jpsi",&segCompatibility_mu_sub_jpsi,32000,0);
  jpsi_->Branch("caloCompatibility_mu_sub_jpsi",&caloCompatibility_mu_sub_jpsi,32000,0);
  jpsi_->Branch("chi2LocalPosition_mu_sub_jpsi",&chi2LocalPosition_mu_sub_jpsi,32000,0);
  jpsi_->Branch("trkKink_mu_sub_jpsi",&trkKink_mu_sub_jpsi,32000,0);
  jpsi_->Branch("glbKink_mu_sub_jpsi",&glbKink_mu_sub_jpsi,32000,0);
  jpsi_->Branch("pfrelIso04_mu_sub_jpsi",&pfrelIso04_mu_sub_jpsi,32000,0);
  jpsi_->Branch("nValidMuonHits_mu_sub_jpsi",&nValidMuonHits_mu_sub_jpsi,32000,0);
  jpsi_->Branch("normChi2_global_mu_sub_jpsi",&normChi2_global_mu_sub_jpsi,32000,0);
  jpsi_->Branch("nValidPixelHits_mu_sub_jpsi",&nValidPixelHits_mu_sub_jpsi,32000,0);
  jpsi_->Branch("trkLayers_mu_sub_jpsi",&trkLayers_mu_sub_jpsi,32000,0);
  jpsi_->Branch("pixelLayers_mu_sub_jpsi",&pixelLayers_mu_sub_jpsi,32000,0);
  jpsi_->Branch("normChi2_inner_mu_sub_jpsi",&normChi2_inner_mu_sub_jpsi,32000,0);
  jpsi_->Branch("validFraction_mu_sub_jpsi",&validFraction_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tunePtErrorOverPt_mu_sub_jpsi",&tunePtErrorOverPt_mu_sub_jpsi,32000,0);
  jpsi_->Branch("dxy_PV_mu_sub_jpsi",&dxy_PV_mu_sub_jpsi,32000,0);
  jpsi_->Branch("dz_PV_mu_sub_jpsi",&dz_PV_mu_sub_jpsi,32000,0);
  jpsi_->Branch("nMatches_mu_sub_jpsi",&nMatches_mu_sub_jpsi,32000,0);
  jpsi_->Branch("innerTrk_pt_mu_sub_jpsi",&innerTrk_pt_mu_sub_jpsi,32000,0);
  jpsi_->Branch("innerTrk_eta_mu_sub_jpsi",&innerTrk_eta_mu_sub_jpsi,32000,0);
  jpsi_->Branch("innerTrk_phi_mu_sub_jpsi",&innerTrk_phi_mu_sub_jpsi,32000,0);
  jpsi_->Branch("outerTrk_pt_mu_sub_jpsi",&outerTrk_pt_mu_sub_jpsi,32000,0);
  jpsi_->Branch("outerTrk_eta_mu_sub_jpsi",&outerTrk_eta_mu_sub_jpsi,32000,0);
  jpsi_->Branch("outerTrk_phi_mu_sub_jpsi",&outerTrk_phi_mu_sub_jpsi,32000,0);
  jpsi_->Branch("globalTrk_pt_mu_sub_jpsi",&globalTrk_pt_mu_sub_jpsi,32000,0);
  jpsi_->Branch("globalTrk_eta_mu_sub_jpsi",&globalTrk_eta_mu_sub_jpsi,32000,0);
  jpsi_->Branch("globalTrk_phi_mu_sub_jpsi",&globalTrk_phi_mu_sub_jpsi,32000,0);
  jpsi_->Branch("bestTrk_pt_mu_sub_jpsi",&bestTrk_pt_mu_sub_jpsi,32000,0);
  jpsi_->Branch("bestTrk_eta_mu_sub_jpsi",&bestTrk_eta_mu_sub_jpsi,32000,0);
  jpsi_->Branch("bestTrk_phi_mu_sub_jpsi",&bestTrk_phi_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_pt_mu_sub_jpsi",&tunePTrk_pt_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_eta_mu_sub_jpsi",&tunePTrk_eta_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_phi_mu_sub_jpsi",&tunePTrk_phi_mu_sub_jpsi,32000,0);
  jpsi_->Branch("muon_charge_mu_sub_jpsi",&muon_charge_mu_sub_jpsi,32000,0);
  jpsi_->Branch("innerTrk_charge_mu_sub_jpsi",&innerTrk_charge_mu_sub_jpsi,32000,0);
  jpsi_->Branch("outerTrk_charge_mu_sub_jpsi",&outerTrk_charge_mu_sub_jpsi,32000,0);
  jpsi_->Branch("globalTrk_charge_mu_sub_jpsi",&globalTrk_charge_mu_sub_jpsi,32000,0);
  jpsi_->Branch("bestTrk_charge_mu_sub_jpsi",&bestTrk_charge_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_charge_mu_sub_jpsi",&tunePTrk_charge_mu_sub_jpsi,32000,0);
  jpsi_->Branch("innerTrk_ptError_mu_sub_jpsi",&innerTrk_ptError_mu_sub_jpsi,32000,0);
  jpsi_->Branch("outerTrk_ptError_mu_sub_jpsi",&outerTrk_ptError_mu_sub_jpsi,32000,0);
  jpsi_->Branch("globalTrk_ptError_mu_sub_jpsi",&globalTrk_ptError_mu_sub_jpsi,32000,0);
  jpsi_->Branch("bestTrk_ptError_mu_sub_jpsi",&bestTrk_ptError_mu_sub_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_ptError_mu_sub_jpsi",&tunePTrk_ptError_mu_sub_jpsi,32000,0);
  jpsi_->Branch("dxy_bestTrk_PV_mu_sub_jpsi",&dxy_bestTrk_PV_mu_sub_jpsi,32000,0);
  jpsi_->Branch("dz_bestTrk_PV_mu_sub_jpsi",&dz_bestTrk_PV_mu_sub_jpsi,32000,0);
  jpsi_->Branch("pT_mu_jpsi",&pT_mu_jpsi,32000,0);
  jpsi_->Branch("eta_mu_jpsi",&eta_mu_jpsi,32000,0);
  jpsi_->Branch("phi_mu_jpsi",&phi_mu_jpsi,32000,0);
  jpsi_->Branch("idx_mu_jpsi",&idx_mu_jpsi,32000,0);
  jpsi_->Branch("looseID_mu_jpsi",&looseID_mu_jpsi,32000,0);
  jpsi_->Branch("highptID_mu_jpsi",&highptID_mu_jpsi,32000,0);
  jpsi_->Branch("trackerhighptID_mu_jpsi",&trackerhighptID_mu_jpsi,32000,0);
  jpsi_->Branch("mediumID_mu_jpsi",&mediumID_mu_jpsi,32000,0);
  jpsi_->Branch("tightID_mu_jpsi",&tightID_mu_jpsi,32000,0);
  jpsi_->Branch("isGlobal_mu_jpsi",&isGlobal_mu_jpsi,32000,0);
  jpsi_->Branch("isTracker_mu_jpsi",&isTracker_mu_jpsi,32000,0);
  jpsi_->Branch("isPFMuon_mu_jpsi",&isPFMuon_mu_jpsi,32000,0);
  jpsi_->Branch("isPFIsolationValid_mu_jpsi",&isPFIsolationValid_mu_jpsi,32000,0);
  jpsi_->Branch("em_mu_jpsi",&em_mu_jpsi,32000,0);
  jpsi_->Branch("emS9_mu_jpsi",&emS9_mu_jpsi,32000,0);
  jpsi_->Branch("emS25_mu_jpsi",&emS25_mu_jpsi,32000,0);
  jpsi_->Branch("emMax_mu_jpsi",&emMax_mu_jpsi,32000,0);
  jpsi_->Branch("had_mu_jpsi",&had_mu_jpsi,32000,0);
  jpsi_->Branch("hadS9_mu_jpsi",&hadS9_mu_jpsi,32000,0);
  jpsi_->Branch("hadMax_mu_jpsi",&hadMax_mu_jpsi,32000,0);
  jpsi_->Branch("ho_mu_jpsi",&ho_mu_jpsi,32000,0);
  jpsi_->Branch("hoS9_mu_jpsi",&hoS9_mu_jpsi,32000,0);
  jpsi_->Branch("numOfMatchedStations_mu_jpsi",&numOfMatchedStations_mu_jpsi,32000,0);
  jpsi_->Branch("numOfChambers_mu_jpsi",&numOfChambers_mu_jpsi,32000,0);
  jpsi_->Branch("segCompatibility_mu_jpsi",&segCompatibility_mu_jpsi,32000,0);
  jpsi_->Branch("caloCompatibility_mu_jpsi",&caloCompatibility_mu_jpsi,32000,0);
  jpsi_->Branch("chi2LocalPosition_mu_jpsi",&chi2LocalPosition_mu_jpsi,32000,0);
  jpsi_->Branch("trkKink_mu_jpsi",&trkKink_mu_jpsi,32000,0);
  jpsi_->Branch("glbKink_mu_jpsi",&glbKink_mu_jpsi,32000,0);
  jpsi_->Branch("pfrelIso04_mu_jpsi",&pfrelIso04_mu_jpsi,32000,0);
  jpsi_->Branch("nValidMuonHits_mu_jpsi",&nValidMuonHits_mu_jpsi,32000,0);
  jpsi_->Branch("normChi2_global_mu_jpsi",&normChi2_global_mu_jpsi,32000,0);
  jpsi_->Branch("nValidPixelHits_mu_jpsi",&nValidPixelHits_mu_jpsi,32000,0);
  jpsi_->Branch("trkLayers_mu_jpsi",&trkLayers_mu_jpsi,32000,0);
  jpsi_->Branch("pixelLayers_mu_jpsi",&pixelLayers_mu_jpsi,32000,0);
  jpsi_->Branch("normChi2_inner_mu_jpsi",&normChi2_inner_mu_jpsi,32000,0);
  jpsi_->Branch("validFraction_mu_jpsi",&validFraction_mu_jpsi,32000,0);
  jpsi_->Branch("tunePtErrorOverPt_mu_jpsi",&tunePtErrorOverPt_mu_jpsi,32000,0);
  jpsi_->Branch("dxy_PV_mu_jpsi",&dxy_PV_mu_jpsi,32000,0);
  jpsi_->Branch("dz_PV_mu_jpsi",&dz_PV_mu_jpsi,32000,0);
  jpsi_->Branch("nMatches_mu_jpsi",&nMatches_mu_jpsi,32000,0);
  jpsi_->Branch("innerTrk_pt_mu_jpsi",&innerTrk_pt_mu_jpsi,32000,0);
  jpsi_->Branch("innerTrk_eta_mu_jpsi",&innerTrk_eta_mu_jpsi,32000,0);
  jpsi_->Branch("innerTrk_phi_mu_jpsi",&innerTrk_phi_mu_jpsi,32000,0);
  jpsi_->Branch("outerTrk_pt_mu_jpsi",&outerTrk_pt_mu_jpsi,32000,0);
  jpsi_->Branch("outerTrk_eta_mu_jpsi",&outerTrk_eta_mu_jpsi,32000,0);
  jpsi_->Branch("outerTrk_phi_mu_jpsi",&outerTrk_phi_mu_jpsi,32000,0);
  jpsi_->Branch("globalTrk_pt_mu_jpsi",&globalTrk_pt_mu_jpsi,32000,0);
  jpsi_->Branch("globalTrk_eta_mu_jpsi",&globalTrk_eta_mu_jpsi,32000,0);
  jpsi_->Branch("globalTrk_phi_mu_jpsi",&globalTrk_phi_mu_jpsi,32000,0);
  jpsi_->Branch("bestTrk_pt_mu_jpsi",&bestTrk_pt_mu_jpsi,32000,0);
  jpsi_->Branch("bestTrk_eta_mu_jpsi",&bestTrk_eta_mu_jpsi,32000,0);
  jpsi_->Branch("bestTrk_phi_mu_jpsi",&bestTrk_phi_mu_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_pt_mu_jpsi",&tunePTrk_pt_mu_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_eta_mu_jpsi",&tunePTrk_eta_mu_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_phi_mu_jpsi",&tunePTrk_phi_mu_jpsi,32000,0);
  jpsi_->Branch("muon_charge_mu_jpsi",&muon_charge_mu_jpsi,32000,0);
  jpsi_->Branch("innerTrk_charge_mu_jpsi",&innerTrk_charge_mu_jpsi,32000,0);
  jpsi_->Branch("outerTrk_charge_mu_jpsi",&outerTrk_charge_mu_jpsi,32000,0);
  jpsi_->Branch("globalTrk_charge_mu_jpsi",&globalTrk_charge_mu_jpsi,32000,0);
  jpsi_->Branch("bestTrk_charge_mu_jpsi",&bestTrk_charge_mu_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_charge_mu_jpsi",&tunePTrk_charge_mu_jpsi,32000,0);
  jpsi_->Branch("innerTrk_ptError_mu_jpsi",&innerTrk_ptError_mu_jpsi,32000,0);
  jpsi_->Branch("outerTrk_ptError_mu_jpsi",&outerTrk_ptError_mu_jpsi,32000,0);
  jpsi_->Branch("globalTrk_ptError_mu_jpsi",&globalTrk_ptError_mu_jpsi,32000,0);
  jpsi_->Branch("bestTrk_ptError_mu_jpsi",&bestTrk_ptError_mu_jpsi,32000,0);
  jpsi_->Branch("tunePTrk_ptError_mu_jpsi",&tunePTrk_ptError_mu_jpsi,32000,0);
  jpsi_->Branch("dxy_bestTrk_PV_mu_jpsi",&dxy_bestTrk_PV_mu_jpsi,32000,0);
  jpsi_->Branch("dz_bestTrk_PV_mu_jpsi",&dz_bestTrk_PV_mu_jpsi,32000,0);

  muon_ = fs->make<TTree>("muonTree","muonTree");
  muon_->Branch("runNo",&runNo,"runNo/I");
  muon_->Branch("lumiNo",&lumiNo,"lumiNo/I");
  muon_->Branch("evtNo",&evtNo,"evtNo/I");
  muon_->Branch("pT_isomu",&pT_muoniso,"pT_isomu/F");
  muon_->Branch("eta_isomu",&eta_muoniso,"eta_isomu/F");
  muon_->Branch("phi_isomu",&phi_muoniso,"phi_isomu/F");
  muon_->Branch("mediumID_isomu",&mediumID_muoniso,"mediumID_isomu/O");
  muon_->Branch("iso04_isomu",&iso04_muoniso,"iso04_isomu/F");
  muon_->Branch("weight_mu",&weight_mu,"weight_mu/I");
  muon_->Branch("MET",&MET,"MET/F");
  muon_->Branch("MET_phi",&MET_phi,"MET_phi/F"); 
  muon_->Branch("MT",&MT,"MT/F"); 
  muon_->Branch("dphi",&dphi,"dphi/F"); 
  muon_->Branch("invM",&invM,32000,0);
  muon_->Branch("pT_mu_sub",&pT_mu_sub,32000,0);
  muon_->Branch("eta_mu_sub",&eta_mu_sub,32000,0);
  muon_->Branch("phi_mu_sub",&phi_mu_sub,32000,0);
  muon_->Branch("idx_mu_sub",&idx_mu_sub,32000,0);
  muon_->Branch("looseID_mu_sub",&looseID_mu_sub,32000,0);
  muon_->Branch("highptID_mu_sub",&highptID_mu_sub,32000,0);
  muon_->Branch("trackerhighptID_mu_sub",&trackerhighptID_mu_sub,32000,0);
  muon_->Branch("mediumID_mu_sub",&mediumID_mu_sub,32000,0);
  muon_->Branch("tightID_mu_sub",&tightID_mu_sub,32000,0);
  muon_->Branch("isGlobal_mu_sub",&isGlobal_mu_sub,32000,0);
  muon_->Branch("isTracker_mu_sub",&isTracker_mu_sub,32000,0);
  muon_->Branch("isPFMuon_mu_sub",&isPFMuon_mu_sub,32000,0);
  muon_->Branch("isPFIsolationValid_mu_sub",&isPFIsolationValid_mu_sub,32000,0);
  muon_->Branch("em_mu_sub",&em_mu_sub,32000,0);
  muon_->Branch("emS9_mu_sub",&emS9_mu_sub,32000,0);
  muon_->Branch("emS25_mu_sub",&emS25_mu_sub,32000,0);
  muon_->Branch("emMax_mu_sub",&emMax_mu_sub,32000,0);
  muon_->Branch("had_mu_sub",&had_mu_sub,32000,0);
  muon_->Branch("hadS9_mu_sub",&hadS9_mu_sub,32000,0);
  muon_->Branch("hadMax_mu_sub",&hadMax_mu_sub,32000,0);
  muon_->Branch("ho_mu_sub",&ho_mu_sub,32000,0);
  muon_->Branch("hoS9_mu_sub",&hoS9_mu_sub,32000,0);
  muon_->Branch("numOfMatchedStations_mu_sub",&numOfMatchedStations_mu_sub,32000,0);
  muon_->Branch("numOfChambers_mu_sub",&numOfChambers_mu_sub,32000,0);
  muon_->Branch("segCompatibility_mu_sub",&segCompatibility_mu_sub,32000,0);
  muon_->Branch("caloCompatibility_mu_sub",&caloCompatibility_mu_sub,32000,0);
  muon_->Branch("chi2LocalPosition_mu_sub",&chi2LocalPosition_mu_sub,32000,0);
  muon_->Branch("trkKink_mu_sub",&trkKink_mu_sub,32000,0);
  muon_->Branch("glbKink_mu_sub",&glbKink_mu_sub,32000,0);
  muon_->Branch("pfrelIso04_mu_sub",&pfrelIso04_mu_sub,32000,0);
  muon_->Branch("nValidMuonHits_mu_sub",&nValidMuonHits_mu_sub,32000,0);
  muon_->Branch("normChi2_global_mu_sub",&normChi2_global_mu_sub,32000,0);
  muon_->Branch("nValidPixelHits_mu_sub",&nValidPixelHits_mu_sub,32000,0);
  muon_->Branch("trkLayers_mu_sub",&trkLayers_mu_sub,32000,0);
  muon_->Branch("pixelLayers_mu_sub",&pixelLayers_mu_sub,32000,0);
  muon_->Branch("normChi2_inner_mu_sub",&normChi2_inner_mu_sub,32000,0);
  muon_->Branch("validFraction_mu_sub",&validFraction_mu_sub,32000,0);
  muon_->Branch("tunePtErrorOverPt_mu_sub",&tunePtErrorOverPt_mu_sub,32000,0);
  muon_->Branch("dxy_PV_mu_sub",&dxy_PV_mu_sub,32000,0);
  muon_->Branch("dz_PV_mu_sub",&dz_PV_mu_sub,32000,0);
  muon_->Branch("nMatches_mu_sub",&nMatches_mu_sub,32000,0);
  muon_->Branch("innerTrk_pt_mu_sub",&innerTrk_pt_mu_sub,32000,0);
  muon_->Branch("innerTrk_eta_mu_sub",&innerTrk_eta_mu_sub,32000,0);
  muon_->Branch("innerTrk_phi_mu_sub",&innerTrk_phi_mu_sub,32000,0);
  muon_->Branch("outerTrk_pt_mu_sub",&outerTrk_pt_mu_sub,32000,0);
  muon_->Branch("outerTrk_eta_mu_sub",&outerTrk_eta_mu_sub,32000,0);
  muon_->Branch("outerTrk_phi_mu_sub",&outerTrk_phi_mu_sub,32000,0);
  muon_->Branch("globalTrk_pt_mu_sub",&globalTrk_pt_mu_sub,32000,0);
  muon_->Branch("globalTrk_eta_mu_sub",&globalTrk_eta_mu_sub,32000,0);
  muon_->Branch("globalTrk_phi_mu_sub",&globalTrk_phi_mu_sub,32000,0);
  muon_->Branch("bestTrk_pt_mu_sub",&bestTrk_pt_mu_sub,32000,0);
  muon_->Branch("bestTrk_eta_mu_sub",&bestTrk_eta_mu_sub,32000,0);
  muon_->Branch("bestTrk_phi_mu_sub",&bestTrk_phi_mu_sub,32000,0);
  muon_->Branch("tunePTrk_pt_mu_sub",&tunePTrk_pt_mu_sub,32000,0);
  muon_->Branch("tunePTrk_eta_mu_sub",&tunePTrk_eta_mu_sub,32000,0);
  muon_->Branch("tunePTrk_phi_mu_sub",&tunePTrk_phi_mu_sub,32000,0);
  muon_->Branch("muon_charge_mu_sub",&muon_charge_mu_sub,32000,0);
  muon_->Branch("innerTrk_charge_mu_sub",&innerTrk_charge_mu_sub,32000,0);
  muon_->Branch("outerTrk_charge_mu_sub",&outerTrk_charge_mu_sub,32000,0);
  muon_->Branch("globalTrk_charge_mu_sub",&globalTrk_charge_mu_sub,32000,0);
  muon_->Branch("bestTrk_charge_mu_sub",&bestTrk_charge_mu_sub,32000,0);
  muon_->Branch("tunePTrk_charge_mu_sub",&tunePTrk_charge_mu_sub,32000,0);
  muon_->Branch("innerTrk_ptError_mu_sub",&innerTrk_ptError_mu_sub,32000,0);
  muon_->Branch("outerTrk_ptError_mu_sub",&outerTrk_ptError_mu_sub,32000,0);
  muon_->Branch("globalTrk_ptError_mu_sub",&globalTrk_ptError_mu_sub,32000,0);
  muon_->Branch("bestTrk_ptError_mu_sub",&bestTrk_ptError_mu_sub,32000,0);
  muon_->Branch("tunePTrk_ptError_mu_sub",&tunePTrk_ptError_mu_sub,32000,0);
  muon_->Branch("dxy_bestTrk_PV_mu_sub",&dxy_bestTrk_PV_mu_sub,32000,0);
  muon_->Branch("dz_bestTrk_PV_mu_sub",&dz_bestTrk_PV_mu_sub,32000,0);
  muon_->Branch("pT_mu",&pT_mu,32000,0);
  muon_->Branch("eta_mu",&eta_mu,32000,0);
  muon_->Branch("phi_mu",&phi_mu,32000,0);
  muon_->Branch("idx_mu",&idx_mu,32000,0);
  muon_->Branch("looseID_mu",&looseID_mu,32000,0);
  muon_->Branch("highptID_mu",&highptID_mu,32000,0);
  muon_->Branch("trackerhighptID_mu",&trackerhighptID_mu,32000,0);
  muon_->Branch("mediumID_mu",&mediumID_mu,32000,0);
  muon_->Branch("tightID_mu",&tightID_mu,32000,0);
  muon_->Branch("isGlobal_mu",&isGlobal_mu,32000,0);
  muon_->Branch("isTracker_mu",&isTracker_mu,32000,0);
  muon_->Branch("isPFMuon_mu",&isPFMuon_mu,32000,0);
  muon_->Branch("isPFIsolationValid_mu",&isPFIsolationValid_mu,32000,0);
  muon_->Branch("em_mu",&em_mu,32000,0);
  muon_->Branch("emS9_mu",&emS9_mu,32000,0);
  muon_->Branch("emS25_mu",&emS25_mu,32000,0);
  muon_->Branch("emMax_mu",&emMax_mu,32000,0);
  muon_->Branch("had_mu",&had_mu,32000,0);
  muon_->Branch("hadS9_mu",&hadS9_mu,32000,0);
  muon_->Branch("hadMax_mu",&hadMax_mu,32000,0);
  muon_->Branch("ho_mu",&ho_mu,32000,0);
  muon_->Branch("hoS9_mu",&hoS9_mu,32000,0);
  muon_->Branch("numOfMatchedStations_mu",&numOfMatchedStations_mu,32000,0);
  muon_->Branch("numOfChambers_mu",&numOfChambers_mu,32000,0);
  muon_->Branch("segCompatibility_mu",&segCompatibility_mu,32000,0);
  muon_->Branch("caloCompatibility_mu",&caloCompatibility_mu,32000,0);
  muon_->Branch("chi2LocalPosition_mu",&chi2LocalPosition_mu,32000,0);
  muon_->Branch("trkKink_mu",&trkKink_mu,32000,0);
  muon_->Branch("glbKink_mu",&glbKink_mu,32000,0);
  muon_->Branch("pfrelIso04_mu",&pfrelIso04_mu,32000,0);
  muon_->Branch("nValidMuonHits_mu",&nValidMuonHits_mu,32000,0);
  muon_->Branch("normChi2_global_mu",&normChi2_global_mu,32000,0);
  muon_->Branch("nValidPixelHits_mu",&nValidPixelHits_mu,32000,0);
  muon_->Branch("trkLayers_mu",&trkLayers_mu,32000,0);
  muon_->Branch("pixelLayers_mu",&pixelLayers_mu,32000,0);
  muon_->Branch("normChi2_inner_mu",&normChi2_inner_mu,32000,0);
  muon_->Branch("validFraction_mu",&validFraction_mu,32000,0);
  muon_->Branch("tunePtErrorOverPt_mu",&tunePtErrorOverPt_mu,32000,0);
  muon_->Branch("dxy_PV_mu",&dxy_PV_mu,32000,0);
  muon_->Branch("dz_PV_mu",&dz_PV_mu,32000,0);
  muon_->Branch("nMatches_mu",&nMatches_mu,32000,0);
  muon_->Branch("innerTrk_pt_mu",&innerTrk_pt_mu,32000,0);
  muon_->Branch("innerTrk_eta_mu",&innerTrk_eta_mu,32000,0);
  muon_->Branch("innerTrk_phi_mu",&innerTrk_phi_mu,32000,0);
  muon_->Branch("outerTrk_pt_mu",&outerTrk_pt_mu,32000,0);
  muon_->Branch("outerTrk_eta_mu",&outerTrk_eta_mu,32000,0);
  muon_->Branch("outerTrk_phi_mu",&outerTrk_phi_mu,32000,0);
  muon_->Branch("globalTrk_pt_mu",&globalTrk_pt_mu,32000,0);
  muon_->Branch("globalTrk_eta_mu",&globalTrk_eta_mu,32000,0);
  muon_->Branch("globalTrk_phi_mu",&globalTrk_phi_mu,32000,0);
  muon_->Branch("bestTrk_pt_mu",&bestTrk_pt_mu,32000,0);
  muon_->Branch("bestTrk_eta_mu",&bestTrk_eta_mu,32000,0);
  muon_->Branch("bestTrk_phi_mu",&bestTrk_phi_mu,32000,0);
  muon_->Branch("tunePTrk_pt_mu",&tunePTrk_pt_mu,32000,0);
  muon_->Branch("tunePTrk_eta_mu",&tunePTrk_eta_mu,32000,0);
  muon_->Branch("tunePTrk_phi_mu",&tunePTrk_phi_mu,32000,0);
  muon_->Branch("muon_charge_mu",&muon_charge_mu,32000,0);
  muon_->Branch("innerTrk_charge_mu",&innerTrk_charge_mu,32000,0);
  muon_->Branch("outerTrk_charge_mu",&outerTrk_charge_mu,32000,0);
  muon_->Branch("globalTrk_charge_mu",&globalTrk_charge_mu,32000,0);
  muon_->Branch("bestTrk_charge_mu",&bestTrk_charge_mu,32000,0);
  muon_->Branch("tunePTrk_charge_mu",&tunePTrk_charge_mu,32000,0);
  muon_->Branch("innerTrk_ptError_mu",&innerTrk_ptError_mu,32000,0);
  muon_->Branch("outerTrk_ptError_mu",&outerTrk_ptError_mu,32000,0);
  muon_->Branch("globalTrk_ptError_mu",&globalTrk_ptError_mu,32000,0);
  muon_->Branch("bestTrk_ptError_mu",&bestTrk_ptError_mu,32000,0);
  muon_->Branch("tunePTrk_ptError_mu",&tunePTrk_ptError_mu,32000,0);
  muon_->Branch("dxy_bestTrk_PV_mu",&dxy_bestTrk_PV_mu,32000,0);
  muon_->Branch("dz_bestTrk_PV_mu",&dz_bestTrk_PV_mu,32000,0);

  jpsiGen_ = fs->make<TTree>("jpsiGenTree","jpsiGenTree");
  jpsiGen_->Branch("runNo_gen",&runNo_gen,"runNo_gen/I");
  jpsiGen_->Branch("lumiNo_gen",&lumiNo_gen,"lumiNo_gen/I");
  jpsiGen_->Branch("evtNo_gen",&evtNo_gen,"evtNo_gen/I");
  jpsiGen_->Branch("weight_gen",&weight_gen,"weight_gen/F");
  jpsiGen_->Branch("jpsi_pt_gen",&jpsi_pt_gen,32000,0);
  jpsiGen_->Branch("jpsi_eta_gen",&jpsi_eta_gen,32000,0);
  jpsiGen_->Branch("jpsi_phi_gen",&jpsi_phi_gen,32000,0);
  jpsiGen_->Branch("jpsi_mass_gen",&jpsi_mass_gen,32000,0);
  jpsiGen_->Branch("gen_1_pt",&gen_1_pt,32000,0);
  jpsiGen_->Branch("gen_1_eta",&gen_1_eta,32000,0);
  jpsiGen_->Branch("gen_1_phi",&gen_1_phi,32000,0);
  jpsiGen_->Branch("gen_1_charge",&gen_1_charge,32000,0);
  jpsiGen_->Branch("gen_2_pt",&gen_2_pt,32000,0);
  jpsiGen_->Branch("gen_2_eta",&gen_2_eta,32000,0);
  jpsiGen_->Branch("gen_2_phi",&gen_2_phi,32000,0);
  jpsiGen_->Branch("gen_2_charge",&gen_2_charge,32000,0);
  jpsiGen_->Branch("gen_dR",&gen_dR,32000,0);
  jpsiGen_->Branch("gen_invM",&gen_invM,32000,0);
}


void MuonJpsi::endJob() {
  purwgtFile_->Close();
}


void MuonJpsi::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  runNo = iEvent.id().run();
  lumiNo = iEvent.id().luminosityBlock();
  evtNo = iEvent.id().event();
  runNo_jpsi = runNo;
  lumiNo_jpsi = lumiNo;
  evtNo_jpsi = evtNo;

  invM_jpsi.clear();
  pT_mu_sub_jpsi.clear();
  eta_mu_sub_jpsi.clear();
  phi_mu_sub_jpsi.clear();
  idx_mu_sub_jpsi.clear();
  looseID_mu_sub_jpsi.clear();
  highptID_mu_sub_jpsi.clear();
  trackerhighptID_mu_sub_jpsi.clear();
  mediumID_mu_sub_jpsi.clear();
  tightID_mu_sub_jpsi.clear();
  isGlobal_mu_sub_jpsi.clear();
  isTracker_mu_sub_jpsi.clear();
  isPFMuon_mu_sub_jpsi.clear();
  isPFIsolationValid_mu_sub_jpsi.clear();
  em_mu_sub_jpsi.clear();
  emS9_mu_sub_jpsi.clear();
  emS25_mu_sub_jpsi.clear();
  emMax_mu_sub_jpsi.clear();
  had_mu_sub_jpsi.clear();
  hadS9_mu_sub_jpsi.clear();
  hadMax_mu_sub_jpsi.clear();
  ho_mu_sub_jpsi.clear();
  hoS9_mu_sub_jpsi.clear();
  numOfMatchedStations_mu_sub_jpsi.clear();
  numOfChambers_mu_sub_jpsi.clear();
  segCompatibility_mu_sub_jpsi.clear();
  caloCompatibility_mu_sub_jpsi.clear();
  chi2LocalPosition_mu_sub_jpsi.clear();
  trkKink_mu_sub_jpsi.clear();
  glbKink_mu_sub_jpsi.clear();
  pfrelIso04_mu_sub_jpsi.clear();
  nValidMuonHits_mu_sub_jpsi.clear();
  normChi2_global_mu_sub_jpsi.clear();
  nValidPixelHits_mu_sub_jpsi.clear();
  trkLayers_mu_sub_jpsi.clear();
  pixelLayers_mu_sub_jpsi.clear();
  normChi2_inner_mu_sub_jpsi.clear();
  validFraction_mu_sub_jpsi.clear();
  tunePtErrorOverPt_mu_sub_jpsi.clear();
  dxy_PV_mu_sub_jpsi.clear();
  dz_PV_mu_sub_jpsi.clear();
  nMatches_mu_sub_jpsi.clear();
  innerTrk_pt_mu_sub_jpsi.clear();
  innerTrk_eta_mu_sub_jpsi.clear();
  innerTrk_phi_mu_sub_jpsi.clear();
  outerTrk_pt_mu_sub_jpsi.clear();
  outerTrk_eta_mu_sub_jpsi.clear();
  outerTrk_phi_mu_sub_jpsi.clear();
  globalTrk_pt_mu_sub_jpsi.clear();
  globalTrk_eta_mu_sub_jpsi.clear();
  globalTrk_phi_mu_sub_jpsi.clear();
  bestTrk_pt_mu_sub_jpsi.clear();
  bestTrk_eta_mu_sub_jpsi.clear();
  bestTrk_phi_mu_sub_jpsi.clear();
  tunePTrk_pt_mu_sub_jpsi.clear();
  tunePTrk_eta_mu_sub_jpsi.clear();
  tunePTrk_phi_mu_sub_jpsi.clear();
  muon_charge_mu_sub_jpsi.clear();
  innerTrk_charge_mu_sub_jpsi.clear();
  outerTrk_charge_mu_sub_jpsi.clear();
  globalTrk_charge_mu_sub_jpsi.clear();
  bestTrk_charge_mu_sub_jpsi.clear();
  tunePTrk_charge_mu_sub_jpsi.clear();
  innerTrk_ptError_mu_sub_jpsi.clear();
  outerTrk_ptError_mu_sub_jpsi.clear();
  globalTrk_ptError_mu_sub_jpsi.clear();
  bestTrk_ptError_mu_sub_jpsi.clear();
  tunePTrk_ptError_mu_sub_jpsi.clear();
  dxy_bestTrk_PV_mu_sub_jpsi.clear();
  dz_bestTrk_PV_mu_sub_jpsi.clear();
  pT_mu_jpsi.clear();
  eta_mu_jpsi.clear();
  phi_mu_jpsi.clear();
  idx_mu_jpsi.clear();
  looseID_mu_jpsi.clear();
  highptID_mu_jpsi.clear();
  trackerhighptID_mu_jpsi.clear();
  mediumID_mu_jpsi.clear();
  tightID_mu_jpsi.clear();
  isGlobal_mu_jpsi.clear();
  isTracker_mu_jpsi.clear();
  isPFMuon_mu_jpsi.clear();
  isPFIsolationValid_mu_jpsi.clear();
  em_mu_jpsi.clear();
  emS9_mu_jpsi.clear();
  emS25_mu_jpsi.clear();
  emMax_mu_jpsi.clear();
  had_mu_jpsi.clear();
  hadS9_mu_jpsi.clear();
  hadMax_mu_jpsi.clear();
  ho_mu_jpsi.clear();
  hoS9_mu_jpsi.clear();
  numOfMatchedStations_mu_jpsi.clear();
  numOfChambers_mu_jpsi.clear();
  segCompatibility_mu_jpsi.clear();
  caloCompatibility_mu_jpsi.clear();
  chi2LocalPosition_mu_jpsi.clear();
  trkKink_mu_jpsi.clear();
  glbKink_mu_jpsi.clear();
  pfrelIso04_mu_jpsi.clear();
  nValidMuonHits_mu_jpsi.clear();
  normChi2_global_mu_jpsi.clear();
  nValidPixelHits_mu_jpsi.clear();
  trkLayers_mu_jpsi.clear();
  pixelLayers_mu_jpsi.clear();
  normChi2_inner_mu_jpsi.clear();
  validFraction_mu_jpsi.clear();
  tunePtErrorOverPt_mu_jpsi.clear();
  dxy_PV_mu_jpsi.clear();
  dz_PV_mu_jpsi.clear();
  nMatches_mu_jpsi.clear();
  innerTrk_pt_mu_jpsi.clear();
  innerTrk_eta_mu_jpsi.clear();
  innerTrk_phi_mu_jpsi.clear();
  outerTrk_pt_mu_jpsi.clear();
  outerTrk_eta_mu_jpsi.clear();
  outerTrk_phi_mu_jpsi.clear();
  globalTrk_pt_mu_jpsi.clear();
  globalTrk_eta_mu_jpsi.clear();
  globalTrk_phi_mu_jpsi.clear();
  bestTrk_pt_mu_jpsi.clear();
  bestTrk_eta_mu_jpsi.clear();
  bestTrk_phi_mu_jpsi.clear();
  tunePTrk_pt_mu_jpsi.clear();
  tunePTrk_eta_mu_jpsi.clear();
  tunePTrk_phi_mu_jpsi.clear();
  muon_charge_mu_jpsi.clear();
  innerTrk_charge_mu_jpsi.clear();
  outerTrk_charge_mu_jpsi.clear();
  globalTrk_charge_mu_jpsi.clear();
  bestTrk_charge_mu_jpsi.clear();
  tunePTrk_charge_mu_jpsi.clear();
  innerTrk_ptError_mu_jpsi.clear();
  outerTrk_ptError_mu_jpsi.clear();
  globalTrk_ptError_mu_jpsi.clear();
  bestTrk_ptError_mu_jpsi.clear();
  tunePTrk_ptError_mu_jpsi.clear();
  dxy_bestTrk_PV_mu_jpsi.clear();
  dz_bestTrk_PV_mu_jpsi.clear();

  jpsi_pt_gen.clear();
  jpsi_eta_gen.clear();
  jpsi_phi_gen.clear();
  jpsi_mass_gen.clear();
  gen_1_pt.clear();
  gen_1_eta.clear();
  gen_1_phi.clear();
  gen_1_charge.clear();
  gen_2_pt.clear();
  gen_2_eta.clear();
  gen_2_phi.clear();
  gen_2_charge.clear();
  gen_dR.clear();
  gen_invM.clear();

  invM.clear();
  pT_mu_sub.clear();
  eta_mu_sub.clear();
  phi_mu_sub.clear();
  idx_mu_sub.clear();
  looseID_mu_sub.clear();
  highptID_mu_sub.clear();
  trackerhighptID_mu_sub.clear();
  mediumID_mu_sub.clear();
  tightID_mu_sub.clear();
  isGlobal_mu_sub.clear();
  isTracker_mu_sub.clear();
  isPFMuon_mu_sub.clear();
  isPFIsolationValid_mu_sub.clear();
  em_mu_sub.clear();
  emS9_mu_sub.clear();
  emS25_mu_sub.clear();
  emMax_mu_sub.clear();
  had_mu_sub.clear();
  hadS9_mu_sub.clear();
  hadMax_mu_sub.clear();
  ho_mu_sub.clear();
  hoS9_mu_sub.clear();
  numOfMatchedStations_mu_sub.clear();
  numOfChambers_mu_sub.clear();
  segCompatibility_mu_sub.clear();
  caloCompatibility_mu_sub.clear();
  chi2LocalPosition_mu_sub.clear();
  trkKink_mu_sub.clear();
  glbKink_mu_sub.clear();
  pfrelIso04_mu_sub.clear();
  nValidMuonHits_mu_sub.clear();
  normChi2_global_mu_sub.clear();
  nValidPixelHits_mu_sub.clear();
  trkLayers_mu_sub.clear();
  pixelLayers_mu_sub.clear();
  normChi2_inner_mu_sub.clear();
  validFraction_mu_sub.clear();
  tunePtErrorOverPt_mu_sub.clear();
  dxy_PV_mu_sub.clear();
  dz_PV_mu_sub.clear();
  nMatches_mu_sub.clear();
  innerTrk_pt_mu_sub.clear();
  innerTrk_eta_mu_sub.clear();
  innerTrk_phi_mu_sub.clear();
  outerTrk_pt_mu_sub.clear();
  outerTrk_eta_mu_sub.clear();
  outerTrk_phi_mu_sub.clear();
  globalTrk_pt_mu_sub.clear();
  globalTrk_eta_mu_sub.clear();
  globalTrk_phi_mu_sub.clear();
  bestTrk_pt_mu_sub.clear();
  bestTrk_eta_mu_sub.clear();
  bestTrk_phi_mu_sub.clear();
  tunePTrk_pt_mu_sub.clear();
  tunePTrk_eta_mu_sub.clear();
  tunePTrk_phi_mu_sub.clear();
  muon_charge_mu_sub.clear();
  innerTrk_charge_mu_sub.clear();
  outerTrk_charge_mu_sub.clear();
  globalTrk_charge_mu_sub.clear();
  bestTrk_charge_mu_sub.clear();
  tunePTrk_charge_mu_sub.clear();
  innerTrk_ptError_mu_sub.clear();
  outerTrk_ptError_mu_sub.clear();
  globalTrk_ptError_mu_sub.clear();
  bestTrk_ptError_mu_sub.clear();
  tunePTrk_ptError_mu_sub.clear();
  dxy_bestTrk_PV_mu_sub.clear();
  dz_bestTrk_PV_mu_sub.clear();
  pT_mu.clear();
  eta_mu.clear();
  phi_mu.clear();
  idx_mu.clear();
  looseID_mu.clear();
  highptID_mu.clear();
  trackerhighptID_mu.clear();
  mediumID_mu.clear();
  tightID_mu.clear();
  isGlobal_mu.clear();
  isTracker_mu.clear();
  isPFMuon_mu.clear();
  isPFIsolationValid_mu.clear();
  em_mu.clear();
  emS9_mu.clear();
  emS25_mu.clear();
  emMax_mu.clear();
  had_mu.clear();
  hadS9_mu.clear();
  hadMax_mu.clear();
  ho_mu.clear();
  hoS9_mu.clear();
  numOfMatchedStations_mu.clear();
  numOfChambers_mu.clear();
  segCompatibility_mu.clear();
  caloCompatibility_mu.clear();
  chi2LocalPosition_mu.clear();
  trkKink_mu.clear();
  glbKink_mu.clear();
  pfrelIso04_mu.clear();
  nValidMuonHits_mu.clear();
  normChi2_global_mu.clear();
  nValidPixelHits_mu.clear();
  trkLayers_mu.clear();
  pixelLayers_mu.clear();
  normChi2_inner_mu.clear();
  validFraction_mu.clear();
  tunePtErrorOverPt_mu.clear();
  dxy_PV_mu.clear();
  dz_PV_mu.clear();
  nMatches_mu.clear();
  innerTrk_pt_mu.clear();
  innerTrk_eta_mu.clear();
  innerTrk_phi_mu.clear();
  outerTrk_pt_mu.clear();
  outerTrk_eta_mu.clear();
  outerTrk_phi_mu.clear();
  globalTrk_pt_mu.clear();
  globalTrk_eta_mu.clear();
  globalTrk_phi_mu.clear();
  bestTrk_pt_mu.clear();
  bestTrk_eta_mu.clear();
  bestTrk_phi_mu.clear();
  tunePTrk_pt_mu.clear();
  tunePTrk_eta_mu.clear();
  tunePTrk_phi_mu.clear();
  muon_charge_mu.clear();
  innerTrk_charge_mu.clear();
  outerTrk_charge_mu.clear();
  globalTrk_charge_mu.clear();
  bestTrk_charge_mu.clear();
  tunePTrk_charge_mu.clear();
  innerTrk_ptError_mu.clear();
  outerTrk_ptError_mu.clear();
  globalTrk_ptError_mu.clear();
  bestTrk_ptError_mu.clear();
  tunePTrk_ptError_mu.clear();
  dxy_bestTrk_PV_mu.clear();
  dz_bestTrk_PV_mu.clear();

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);
  double aWeight = 1.;

  if (isMC_) {
    //edm::Handle<double> theprefweight;
    //iEvent.getByToken(prefweight_token, theprefweight);
    //double prefiringweight = *theprefweight;

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    //aWeight = prefiringweight*mcweight/std::abs(mcweight);
    aWeight = mcweight/std::abs(mcweight);
    //std::cout<<"weight : "<<aWeight<<std::endl;

    edm::Handle<edm::View<PileupSummaryInfo>> pusummary;
    iEvent.getByToken(pileupToken_, pusummary);

    //for (unsigned int idx = 0; idx < pusummary->size(); ++idx) {
    //  const auto& apu = pusummary->refAt(idx);

    //  int bx = apu->getBunchCrossing();

    //  if (bx==0) { // in-time PU only
    //    auto npu = apu->getTrueNumInteractions();
    //    aWeight *= purwgt_->GetBinContent( purwgt_->FindBin(apu->getTrueNumInteractions()) );
    //    histo1d_["PUsummary"]->Fill( static_cast<float>(npu)+0.5, aWeight );

    //    break;
    //  }
    //}

    // Gen J/psi -> mu+ mu- 수집 (mother=J/psi인 gen muon만)
    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(genptcToken_, genptcHandle);

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& gp = genptcHandle->refAt(idx);
      if (gp->pdgId() != 443) continue;                // J/psi(1S)
      if (!gp->isLastCopy()) continue;                 // 최종 상태 J/psi만
      const reco::Candidate* mu1 = nullptr;
      const reco::Candidate* mu2 = nullptr;
      for (size_t k = 0; k < gp->numberOfDaughters(); ++k) {
        const reco::Candidate* dau = gp->daughter(k);
        if (!dau) continue;
        if (std::abs(dau->pdgId()) != 13) continue;
        if (!mu1)      mu1 = dau;
        else if (!mu2) mu2 = dau;
      }
      if (!mu1 || !mu2) continue;
      if (mu1->pt() < mu2->pt()) std::swap(mu1, mu2);   // leading이 앞으로

      jpsi_pt_gen.push_back(gp->pt());
      jpsi_eta_gen.push_back(gp->eta());
      jpsi_phi_gen.push_back(gp->phi());
      jpsi_mass_gen.push_back(gp->mass());
      gen_1_pt.push_back(mu1->pt());
      gen_1_eta.push_back(mu1->eta());
      gen_1_phi.push_back(mu1->phi());
      gen_1_charge.push_back(mu1->charge());
      gen_2_pt.push_back(mu2->pt());
      gen_2_eta.push_back(mu2->eta());
      gen_2_phi.push_back(mu2->phi());
      gen_2_charge.push_back(mu2->charge());
      gen_dR.push_back(reco::deltaR(*mu1, *mu2));
      const auto v1 = math::PtEtaPhiMLorentzVector(mu1->pt(), mu1->eta(), mu1->phi(), mumass_);
      const auto v2 = math::PtEtaPhiMLorentzVector(mu2->pt(), mu2->eta(), mu2->phi(), mumass_);
      gen_invM.push_back((float)(v1 + v2).M());
    }

    runNo_gen  = iEvent.id().run();
    lumiNo_gen = iEvent.id().luminosityBlock();
    evtNo_gen  = iEvent.id().event();
    weight_gen = (float)aWeight;
    jpsiGen_->Fill();
  }


  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);
  
  edm::Handle<edm::View<pat::MET>> metHandle;
  iEvent.getByToken(metToken_, metHandle);

  edm::Handle<edm::TriggerResults> METfilterHandle;
  iEvent.getByToken(METfilterToken_,METfilterHandle);
  edm::TriggerNames METfilters = iEvent.triggerNames(*METfilterHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  unsigned int nPassedFilters = 0;

  for (unsigned int iTrig = 0; iTrig < METfilterHandle.product()->size(); iTrig++) {
    const std::string trigname = METfilters.triggerName(iTrig);
 //   std::cout<<trigname<<std::endl;
    if (METfilterHandle.product()->accept(iTrig)) {
      for (const auto& filterName : METfilterList_) {
	if (trigname.find(filterName) != std::string::npos)
	  nPassedFilters++;
	
      }
    }
  }
//  std::cout<<nPassedFilters<<std::endl;
  if (nPassedFilters!=METfilterList_.size())
    return;


  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  std::vector<std::pair<std::string, int>> indices;
  edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);
  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    std::string trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
	//std::cout<<trigList_.at(jTrig)<<" | "<<trigName<<std::endl;
        if (trigResultHandle.product()->accept(iTrig)) {
          isFired = true;
	  //std::cout << "Trigger " << trigList.triggerName(iTrig) << " fired!" << std::endl;

	}
      }
    
    }
  } 
  //if (isFired) std::cout<<"hello "<<std::endl;
  const reco::Vertex* pv = nullptr;
  if (pvHandle.isValid()) {
    for (auto const& v : *pvHandle) {
      if (!v.isFake() && v.ndof() > 4 && std::abs(v.z()) < 24 && v.position().Rho() < 2.0) {
        pv = &v;
        break; // 1st good PV
      }
    }
  }
  if (isFired && pv){
    const auto& aMET = metHandle->at(0);
    int index_isomu = -1;
    bool noIsoMuflag = true;
    for (size_t i = 0; i < muonHandle->size(); ++i) {
      const auto& muon = muonHandle->at(i);
      auto isPFmuon = muon.isPFMuon();
      auto isPFIsoValid = muon.isPFIsolationValid();
      auto isMedium = muon::isMediumMuon(muon);
      float pfrelIso04 = -1.;
      if (isPFmuon && isPFIsoValid){
      	  auto iso04 = muon.pfIsolationR04();
	  pfrelIso04 =
            ( iso04.sumChargedHadronPt
            + std::max(0.f,
                  iso04.sumNeutralHadronEt
                + iso04.sumPhotonEt
                - 0.5f * iso04.sumPUPt ) )
            / muon.pt();
	  if (muon.pt()> 27 && isMedium && pfrelIso04 < 0.15){
	    index_isomu = (int) i;
	    pT_muoniso = muon.pt();
	    eta_muoniso = muon.eta();
	    phi_muoniso = muon.phi();
	    iso04_muoniso = pfrelIso04;
	    mediumID_muoniso = isMedium;
	    dphi = reco::deltaPhi(aMET,muon);
	    MET = aMET.corPt(pat::MET::Type1);
            const auto muon_vec     = math::PtEtaPhiMLorentzVector(muon.pt(), muon.eta(), muon.phi(), mumass_);
            const auto MET_vec     = math::PtEtaPhiMLorentzVector(aMET.corPt(pat::MET::Type1), aMET.p4().eta(), aMET.phi(), aMET.p4().M());
	    const auto sum = muon_vec + MET_vec;
	    MT = (float)sum.Mt();
	    noIsoMuflag = false;


	    break;
	  }
          
	 

        
      }

    } 
    //std::cout<<"hello !!"<<std::endl;
    for (size_t i = 0; i < muonHandle->size(); ++i) {
      const auto& muon = muonHandle->at(i);
      
      for (size_t j = i + 1; j < muonHandle->size(); ++j) {
        const auto& muon_sub = muonHandle->at(j);
        float dR = reco::deltaR(muon, muon_sub);
        if (dR < 0.15 && muon.pt() > 50) {
          const auto muon_vec     = math::PtEtaPhiMLorentzVector(muon.pt(), muon.eta(), muon.phi(), mumass_);
          const auto muon_sub_vec     = math::PtEtaPhiMLorentzVector(muon_sub.pt(), muon_sub.eta(), muon_sub.phi(), mumass_);
          const auto sum = muon_vec + muon_sub_vec;
           
          invM_jpsi.push_back(sum.M());
          pT_mu_jpsi.push_back(muon.pt());
          eta_mu_jpsi.push_back(muon.eta());
          phi_mu_jpsi.push_back(muon.phi());
          idx_mu_jpsi.push_back(i);
          looseID_mu_jpsi.push_back(muon::isLooseMuon(muon));
          highptID_mu_jpsi.push_back(muon::isHighPtMuon(muon,*pv));
          trackerhighptID_mu_jpsi.push_back(muon::isTrackerHighPtMuon(muon,*pv));
          mediumID_mu_jpsi.push_back(muon::isMediumMuon(muon));
          tightID_mu_jpsi.push_back(muon::isTightMuon(muon,*pv));
          isGlobal_mu_jpsi.push_back(muon.isGlobalMuon());
          isTracker_mu_jpsi.push_back(muon.isTrackerMuon());
          isPFMuon_mu_jpsi.push_back(muon.isPFMuon());
          isPFIsolationValid_mu_jpsi.push_back(muon.isPFIsolationValid());
          em_mu_jpsi.push_back(muon.calEnergy().em);
          emS9_mu_jpsi.push_back(muon.calEnergy().emS9);
          emS25_mu_jpsi.push_back(muon.calEnergy().emS25);
          emMax_mu_jpsi.push_back(muon.calEnergy().emMax);
          had_mu_jpsi.push_back(muon.calEnergy().had);
          hadS9_mu_jpsi.push_back(muon.calEnergy().hadS9);
          hadMax_mu_jpsi.push_back(muon.calEnergy().hadMax);
          ho_mu_jpsi.push_back(muon.calEnergy().ho);
          hoS9_mu_jpsi.push_back(muon.calEnergy().hoS9);
          numOfMatchedStations_mu_jpsi.push_back(muon.numberOfMatchedStations());
          numOfChambers_mu_jpsi.push_back(muon.numberOfChambers());
          segCompatibility_mu_jpsi.push_back(muon::segmentCompatibility(muon));
          caloCompatibility_mu_jpsi.push_back(muon.caloCompatibility());
          chi2LocalPosition_mu_jpsi.push_back(muon.combinedQuality().chi2LocalPosition);
          trkKink_mu_jpsi.push_back(muon.combinedQuality().trkKink);
          glbKink_mu_jpsi.push_back(muon.combinedQuality().glbKink);
          float pfrelIso04_val = -999.f;
          if (muon.isPFMuon() && muon.isPFIsolationValid()) {
            auto iso04 = muon.pfIsolationR04();
            pfrelIso04_val = ( iso04.sumChargedHadronPt
                             + std::max(0.f, iso04.sumNeutralHadronEt + iso04.sumPhotonEt - 0.5f*iso04.sumPUPt) )
                             / muon.pt();
          }
          pfrelIso04_mu_jpsi.push_back(pfrelIso04_val);
          const auto& glbTrk = muon.globalTrack();
          const auto& innTrk = muon.innerTrack();
          const auto& outTrk = muon.outerTrack();
          const auto& bestTrk = muon.muonBestTrack();
          const auto& tunePTrk = muon.tunePMuonBestTrack();
          nValidMuonHits_mu_jpsi.push_back(glbTrk.isNonnull() ? glbTrk->hitPattern().numberOfValidMuonHits() : -999);
          normChi2_global_mu_jpsi.push_back(glbTrk.isNonnull() ? (float)glbTrk->normalizedChi2() : -999.f);
          nValidPixelHits_mu_jpsi.push_back(innTrk.isNonnull() ? innTrk->hitPattern().numberOfValidPixelHits() : -999);
          trkLayers_mu_jpsi.push_back(innTrk.isNonnull() ? innTrk->hitPattern().trackerLayersWithMeasurement() : -999);
          pixelLayers_mu_jpsi.push_back(innTrk.isNonnull() ? innTrk->hitPattern().pixelLayersWithMeasurement() : -999);
          normChi2_inner_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->normalizedChi2() : -999.f);
          validFraction_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->validFraction() : -999.f);
          tunePtErrorOverPt_mu_jpsi.push_back(tunePTrk.isNonnull() && tunePTrk->pt() > 0.f ? (float)(tunePTrk->ptError()/tunePTrk->pt()) : -999.f);
          dxy_PV_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->dxy(pv->position()) : -999.f);
          dz_PV_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->dz(pv->position()) : -999.f);
          nMatches_mu_jpsi.push_back(muon.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
          innerTrk_pt_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->pt() : -999.f);
          innerTrk_eta_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->eta() : -999.f);
          innerTrk_phi_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->phi() : -999.f);
          outerTrk_pt_mu_jpsi.push_back(outTrk.isNonnull() ? (float)outTrk->pt() : -999.f);
          outerTrk_eta_mu_jpsi.push_back(outTrk.isNonnull() ? (float)outTrk->eta() : -999.f);
          outerTrk_phi_mu_jpsi.push_back(outTrk.isNonnull() ? (float)outTrk->phi() : -999.f);
          globalTrk_pt_mu_jpsi.push_back(glbTrk.isNonnull() ? (float)glbTrk->pt() : -999.f);
          globalTrk_eta_mu_jpsi.push_back(glbTrk.isNonnull() ? (float)glbTrk->eta() : -999.f);
          globalTrk_phi_mu_jpsi.push_back(glbTrk.isNonnull() ? (float)glbTrk->phi() : -999.f);
          bestTrk_pt_mu_jpsi.push_back(bestTrk.isNonnull() ? (float)bestTrk->pt() : -999.f);
          bestTrk_eta_mu_jpsi.push_back(bestTrk.isNonnull() ? (float)bestTrk->eta() : -999.f);
          bestTrk_phi_mu_jpsi.push_back(bestTrk.isNonnull() ? (float)bestTrk->phi() : -999.f);
          tunePTrk_pt_mu_jpsi.push_back(tunePTrk.isNonnull() ? (float)tunePTrk->pt() : -999.f);
          tunePTrk_eta_mu_jpsi.push_back(tunePTrk.isNonnull() ? (float)tunePTrk->eta() : -999.f);
          tunePTrk_phi_mu_jpsi.push_back(tunePTrk.isNonnull() ? (float)tunePTrk->phi() : -999.f);
          muon_charge_mu_jpsi.push_back(muon.charge());
          innerTrk_charge_mu_jpsi.push_back(innTrk.isNonnull() ? innTrk->charge() : -999);
          outerTrk_charge_mu_jpsi.push_back(outTrk.isNonnull() ? outTrk->charge() : -999);
          globalTrk_charge_mu_jpsi.push_back(glbTrk.isNonnull() ? glbTrk->charge() : -999);
          bestTrk_charge_mu_jpsi.push_back(bestTrk.isNonnull() ? bestTrk->charge() : -999);
          tunePTrk_charge_mu_jpsi.push_back(tunePTrk.isNonnull() ? tunePTrk->charge() : -999);
          innerTrk_ptError_mu_jpsi.push_back(innTrk.isNonnull() ? (float)innTrk->ptError() : -999.f);
          outerTrk_ptError_mu_jpsi.push_back(outTrk.isNonnull() ? (float)outTrk->ptError() : -999.f);
          globalTrk_ptError_mu_jpsi.push_back(glbTrk.isNonnull() ? (float)glbTrk->ptError() : -999.f);
          bestTrk_ptError_mu_jpsi.push_back(bestTrk.isNonnull() ? (float)bestTrk->ptError() : -999.f);
          tunePTrk_ptError_mu_jpsi.push_back(tunePTrk.isNonnull() ? (float)tunePTrk->ptError() : -999.f);
          dxy_bestTrk_PV_mu_jpsi.push_back(bestTrk.isNonnull() ? (float)bestTrk->dxy(pv->position()) : -999.f);
          dz_bestTrk_PV_mu_jpsi.push_back(bestTrk.isNonnull() ? (float)bestTrk->dz(pv->position()) : -999.f);

          pT_mu_sub_jpsi.push_back(muon_sub.pt());
          eta_mu_sub_jpsi.push_back(muon_sub.eta());
          phi_mu_sub_jpsi.push_back(muon_sub.phi());
          idx_mu_sub_jpsi.push_back(j);
          looseID_mu_sub_jpsi.push_back(muon::isLooseMuon(muon_sub));
          highptID_mu_sub_jpsi.push_back(muon::isHighPtMuon(muon_sub,*pv));
          trackerhighptID_mu_sub_jpsi.push_back(muon::isTrackerHighPtMuon(muon_sub,*pv));
          mediumID_mu_sub_jpsi.push_back(muon::isMediumMuon(muon_sub));
          tightID_mu_sub_jpsi.push_back(muon::isTightMuon(muon_sub,*pv));
          isGlobal_mu_sub_jpsi.push_back(muon_sub.isGlobalMuon());
          isTracker_mu_sub_jpsi.push_back(muon_sub.isTrackerMuon());
          isPFMuon_mu_sub_jpsi.push_back(muon_sub.isPFMuon());
          isPFIsolationValid_mu_sub_jpsi.push_back(muon_sub.isPFIsolationValid());
          em_mu_sub_jpsi.push_back(muon_sub.calEnergy().em);
          emS9_mu_sub_jpsi.push_back(muon_sub.calEnergy().emS9);
          emS25_mu_sub_jpsi.push_back(muon_sub.calEnergy().emS25);
          emMax_mu_sub_jpsi.push_back(muon_sub.calEnergy().emMax);
          had_mu_sub_jpsi.push_back(muon_sub.calEnergy().had);
          hadS9_mu_sub_jpsi.push_back(muon_sub.calEnergy().hadS9);
          hadMax_mu_sub_jpsi.push_back(muon_sub.calEnergy().hadMax);
          ho_mu_sub_jpsi.push_back(muon_sub.calEnergy().ho);
          hoS9_mu_sub_jpsi.push_back(muon_sub.calEnergy().hoS9);
          numOfMatchedStations_mu_sub_jpsi.push_back(muon_sub.numberOfMatchedStations());
          numOfChambers_mu_sub_jpsi.push_back(muon_sub.numberOfChambers());
          segCompatibility_mu_sub_jpsi.push_back(muon::segmentCompatibility(muon_sub));
          caloCompatibility_mu_sub_jpsi.push_back(muon_sub.caloCompatibility());
          chi2LocalPosition_mu_sub_jpsi.push_back(muon_sub.combinedQuality().chi2LocalPosition);
          trkKink_mu_sub_jpsi.push_back(muon_sub.combinedQuality().trkKink);
          glbKink_mu_sub_jpsi.push_back(muon_sub.combinedQuality().glbKink);
          float pfrelIso04_sub_val = -999.f;
          if (muon_sub.isPFMuon() && muon_sub.isPFIsolationValid()) {
            auto iso04s = muon_sub.pfIsolationR04();
            pfrelIso04_sub_val = ( iso04s.sumChargedHadronPt
                                 + std::max(0.f, iso04s.sumNeutralHadronEt + iso04s.sumPhotonEt - 0.5f*iso04s.sumPUPt) )
                                 / muon_sub.pt();
          }
          pfrelIso04_mu_sub_jpsi.push_back(pfrelIso04_sub_val);
          const auto& glbTrkS = muon_sub.globalTrack();
          const auto& innTrkS = muon_sub.innerTrack();
          const auto& outTrkS = muon_sub.outerTrack();
          const auto& bestTrkS = muon_sub.muonBestTrack();
          const auto& tunePTrkS = muon_sub.tunePMuonBestTrack();
          nValidMuonHits_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? glbTrkS->hitPattern().numberOfValidMuonHits() : -999);
          normChi2_global_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? (float)glbTrkS->normalizedChi2() : -999.f);
          nValidPixelHits_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? innTrkS->hitPattern().numberOfValidPixelHits() : -999);
          trkLayers_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? innTrkS->hitPattern().trackerLayersWithMeasurement() : -999);
          pixelLayers_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? innTrkS->hitPattern().pixelLayersWithMeasurement() : -999);
          normChi2_inner_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->normalizedChi2() : -999.f);
          validFraction_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->validFraction() : -999.f);
          tunePtErrorOverPt_mu_sub_jpsi.push_back(tunePTrkS.isNonnull() && tunePTrkS->pt() > 0.f ? (float)(tunePTrkS->ptError()/tunePTrkS->pt()) : -999.f);
          dxy_PV_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->dxy(pv->position()) : -999.f);
          dz_PV_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->dz(pv->position()) : -999.f);
          nMatches_mu_sub_jpsi.push_back(muon_sub.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
          innerTrk_pt_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->pt() : -999.f);
          innerTrk_eta_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->eta() : -999.f);
          innerTrk_phi_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->phi() : -999.f);
          outerTrk_pt_mu_sub_jpsi.push_back(outTrkS.isNonnull() ? (float)outTrkS->pt() : -999.f);
          outerTrk_eta_mu_sub_jpsi.push_back(outTrkS.isNonnull() ? (float)outTrkS->eta() : -999.f);
          outerTrk_phi_mu_sub_jpsi.push_back(outTrkS.isNonnull() ? (float)outTrkS->phi() : -999.f);
          globalTrk_pt_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? (float)glbTrkS->pt() : -999.f);
          globalTrk_eta_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? (float)glbTrkS->eta() : -999.f);
          globalTrk_phi_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? (float)glbTrkS->phi() : -999.f);
          bestTrk_pt_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? (float)bestTrkS->pt() : -999.f);
          bestTrk_eta_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? (float)bestTrkS->eta() : -999.f);
          bestTrk_phi_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? (float)bestTrkS->phi() : -999.f);
          tunePTrk_pt_mu_sub_jpsi.push_back(tunePTrkS.isNonnull() ? (float)tunePTrkS->pt() : -999.f);
          tunePTrk_eta_mu_sub_jpsi.push_back(tunePTrkS.isNonnull() ? (float)tunePTrkS->eta() : -999.f);
          tunePTrk_phi_mu_sub_jpsi.push_back(tunePTrkS.isNonnull() ? (float)tunePTrkS->phi() : -999.f);
          muon_charge_mu_sub_jpsi.push_back(muon_sub.charge());
          innerTrk_charge_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? innTrkS->charge() : -999);
          outerTrk_charge_mu_sub_jpsi.push_back(outTrkS.isNonnull() ? outTrkS->charge() : -999);
          globalTrk_charge_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? glbTrkS->charge() : -999);
          bestTrk_charge_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? bestTrkS->charge() : -999);
          tunePTrk_charge_mu_sub_jpsi.push_back(tunePTrkS.isNonnull() ? tunePTrkS->charge() : -999);
          innerTrk_ptError_mu_sub_jpsi.push_back(innTrkS.isNonnull() ? (float)innTrkS->ptError() : -999.f);
          outerTrk_ptError_mu_sub_jpsi.push_back(outTrkS.isNonnull() ? (float)outTrkS->ptError() : -999.f);
          globalTrk_ptError_mu_sub_jpsi.push_back(glbTrkS.isNonnull() ? (float)glbTrkS->ptError() : -999.f);
          bestTrk_ptError_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? (float)bestTrkS->ptError() : -999.f);
          tunePTrk_ptError_mu_sub_jpsi.push_back(tunePTrkS.isNonnull() ? (float)tunePTrkS->ptError() : -999.f);
          dxy_bestTrk_PV_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? (float)bestTrkS->dxy(pv->position()) : -999.f);
          dz_bestTrk_PV_mu_sub_jpsi.push_back(bestTrkS.isNonnull() ? (float)bestTrkS->dz(pv->position()) : -999.f);

          // 예시 출력
          //std::cout << "Invariant mass: " << invM << std::endl;
        }
      }
      //if (invM > 100) std::cout<<"wtf?"<<std::endl;

      // muon 기준 출력 등
      //std::cout << muon.pt() << " | " << muon.eta() << " | " << muon.phi() << " | " << close_muon << std::endl;
    }
    weight_mu_jpsi = aWeight;
    //if (aWeight < 0)std::cout<<aWeight<<std::endl;
    jpsi_->Fill();
    if (noIsoMuflag) return;

    if (aMET.corPt(pat::MET::Type1) > 20 || index_isomu != -1){
      for (size_t i = 0; i < muonHandle->size(); ++i) {
        const auto& muon = muonHandle->at(i);
        
        for (size_t j = i + 1; j < muonHandle->size(); ++j) {
          const auto& muon_sub = muonHandle->at(j);
          float dR = reco::deltaR(muon, muon_sub);
          if (dR < 0.15 && muon.pt() > 50) {
            const auto muon_vec     = math::PtEtaPhiMLorentzVector(muon.pt(), muon.eta(), muon.phi(), mumass_);
            const auto muon_sub_vec     = math::PtEtaPhiMLorentzVector(muon_sub.pt(), muon_sub.eta(), muon_sub.phi(), mumass_);
	    const auto sum = muon_vec + muon_sub_vec;
             
	    invM.push_back(sum.M());
            pT_mu.push_back(muon.pt());
            eta_mu.push_back(muon.eta());
            phi_mu.push_back(muon.phi());
            idx_mu.push_back(i);
	    looseID_mu.push_back(muon::isLooseMuon(muon));
	    highptID_mu.push_back(muon::isHighPtMuon(muon,*pv));
	    trackerhighptID_mu.push_back(muon::isTrackerHighPtMuon(muon,*pv));
            mediumID_mu.push_back(muon::isMediumMuon(muon));
            tightID_mu.push_back(muon::isTightMuon(muon,*pv));
            isGlobal_mu.push_back(muon.isGlobalMuon());
            isTracker_mu.push_back(muon.isTrackerMuon());
            isPFMuon_mu.push_back(muon.isPFMuon());
            isPFIsolationValid_mu.push_back(muon.isPFIsolationValid());
            em_mu.push_back(muon.calEnergy().em);
            emS9_mu.push_back(muon.calEnergy().emS9);
            emS25_mu.push_back(muon.calEnergy().emS25);
            emMax_mu.push_back(muon.calEnergy().emMax);
            had_mu.push_back(muon.calEnergy().had);
            hadS9_mu.push_back(muon.calEnergy().hadS9);
            hadMax_mu.push_back(muon.calEnergy().hadMax);
            ho_mu.push_back(muon.calEnergy().ho);
            hoS9_mu.push_back(muon.calEnergy().hoS9);
            numOfMatchedStations_mu.push_back(muon.numberOfMatchedStations());
            numOfChambers_mu.push_back(muon.numberOfChambers());
            segCompatibility_mu.push_back(muon::segmentCompatibility(muon));
            caloCompatibility_mu.push_back(muon.caloCompatibility());
            chi2LocalPosition_mu.push_back(muon.combinedQuality().chi2LocalPosition);
            trkKink_mu.push_back(muon.combinedQuality().trkKink);
            glbKink_mu.push_back(muon.combinedQuality().glbKink);
            float pfrelIso04_val_m = -999.f;
            if (muon.isPFMuon() && muon.isPFIsolationValid()) {
              auto iso04m = muon.pfIsolationR04();
              pfrelIso04_val_m = ( iso04m.sumChargedHadronPt
                                 + std::max(0.f, iso04m.sumNeutralHadronEt + iso04m.sumPhotonEt - 0.5f*iso04m.sumPUPt) )
                                 / muon.pt();
            }
            pfrelIso04_mu.push_back(pfrelIso04_val_m);
            const auto& glbTrkM = muon.globalTrack();
            const auto& innTrkM = muon.innerTrack();
            const auto& outTrkM = muon.outerTrack();
            const auto& bestTrkM = muon.muonBestTrack();
            const auto& tunePTrkM = muon.tunePMuonBestTrack();
            nValidMuonHits_mu.push_back(glbTrkM.isNonnull() ? glbTrkM->hitPattern().numberOfValidMuonHits() : -999);
            normChi2_global_mu.push_back(glbTrkM.isNonnull() ? (float)glbTrkM->normalizedChi2() : -999.f);
            nValidPixelHits_mu.push_back(innTrkM.isNonnull() ? innTrkM->hitPattern().numberOfValidPixelHits() : -999);
            trkLayers_mu.push_back(innTrkM.isNonnull() ? innTrkM->hitPattern().trackerLayersWithMeasurement() : -999);
            pixelLayers_mu.push_back(innTrkM.isNonnull() ? innTrkM->hitPattern().pixelLayersWithMeasurement() : -999);
            normChi2_inner_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->normalizedChi2() : -999.f);
            validFraction_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->validFraction() : -999.f);
            tunePtErrorOverPt_mu.push_back(tunePTrkM.isNonnull() && tunePTrkM->pt() > 0.f ? (float)(tunePTrkM->ptError()/tunePTrkM->pt()) : -999.f);
            dxy_PV_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->dxy(pv->position()) : -999.f);
            dz_PV_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->dz(pv->position()) : -999.f);
            nMatches_mu.push_back(muon.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
            innerTrk_pt_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->pt() : -999.f);
            innerTrk_eta_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->eta() : -999.f);
            innerTrk_phi_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->phi() : -999.f);
            outerTrk_pt_mu.push_back(outTrkM.isNonnull() ? (float)outTrkM->pt() : -999.f);
            outerTrk_eta_mu.push_back(outTrkM.isNonnull() ? (float)outTrkM->eta() : -999.f);
            outerTrk_phi_mu.push_back(outTrkM.isNonnull() ? (float)outTrkM->phi() : -999.f);
            globalTrk_pt_mu.push_back(glbTrkM.isNonnull() ? (float)glbTrkM->pt() : -999.f);
            globalTrk_eta_mu.push_back(glbTrkM.isNonnull() ? (float)glbTrkM->eta() : -999.f);
            globalTrk_phi_mu.push_back(glbTrkM.isNonnull() ? (float)glbTrkM->phi() : -999.f);
            bestTrk_pt_mu.push_back(bestTrkM.isNonnull() ? (float)bestTrkM->pt() : -999.f);
            bestTrk_eta_mu.push_back(bestTrkM.isNonnull() ? (float)bestTrkM->eta() : -999.f);
            bestTrk_phi_mu.push_back(bestTrkM.isNonnull() ? (float)bestTrkM->phi() : -999.f);
            tunePTrk_pt_mu.push_back(tunePTrkM.isNonnull() ? (float)tunePTrkM->pt() : -999.f);
            tunePTrk_eta_mu.push_back(tunePTrkM.isNonnull() ? (float)tunePTrkM->eta() : -999.f);
            tunePTrk_phi_mu.push_back(tunePTrkM.isNonnull() ? (float)tunePTrkM->phi() : -999.f);
            muon_charge_mu.push_back(muon.charge());
            innerTrk_charge_mu.push_back(innTrkM.isNonnull() ? innTrkM->charge() : -999);
            outerTrk_charge_mu.push_back(outTrkM.isNonnull() ? outTrkM->charge() : -999);
            globalTrk_charge_mu.push_back(glbTrkM.isNonnull() ? glbTrkM->charge() : -999);
            bestTrk_charge_mu.push_back(bestTrkM.isNonnull() ? bestTrkM->charge() : -999);
            tunePTrk_charge_mu.push_back(tunePTrkM.isNonnull() ? tunePTrkM->charge() : -999);
            innerTrk_ptError_mu.push_back(innTrkM.isNonnull() ? (float)innTrkM->ptError() : -999.f);
            outerTrk_ptError_mu.push_back(outTrkM.isNonnull() ? (float)outTrkM->ptError() : -999.f);
            globalTrk_ptError_mu.push_back(glbTrkM.isNonnull() ? (float)glbTrkM->ptError() : -999.f);
            bestTrk_ptError_mu.push_back(bestTrkM.isNonnull() ? (float)bestTrkM->ptError() : -999.f);
            tunePTrk_ptError_mu.push_back(tunePTrkM.isNonnull() ? (float)tunePTrkM->ptError() : -999.f);
            dxy_bestTrk_PV_mu.push_back(bestTrkM.isNonnull() ? (float)bestTrkM->dxy(pv->position()) : -999.f);
            dz_bestTrk_PV_mu.push_back(bestTrkM.isNonnull() ? (float)bestTrkM->dz(pv->position()) : -999.f);

	    pT_mu_sub.push_back(muon_sub.pt());
            eta_mu_sub.push_back(muon_sub.eta());
            phi_mu_sub.push_back(muon_sub.phi());
            idx_mu_sub.push_back(j);
	    looseID_mu_sub.push_back(muon::isLooseMuon(muon_sub));
	    highptID_mu_sub.push_back(muon::isHighPtMuon(muon_sub,*pv));
	    trackerhighptID_mu_sub.push_back(muon::isTrackerHighPtMuon(muon_sub,*pv));
            mediumID_mu_sub.push_back(muon::isMediumMuon(muon_sub));
            tightID_mu_sub.push_back(muon::isTightMuon(muon_sub,*pv));
            isGlobal_mu_sub.push_back(muon_sub.isGlobalMuon());
            isTracker_mu_sub.push_back(muon_sub.isTrackerMuon());
            isPFMuon_mu_sub.push_back(muon_sub.isPFMuon());
            isPFIsolationValid_mu_sub.push_back(muon_sub.isPFIsolationValid());
            em_mu_sub.push_back(muon_sub.calEnergy().em);
            emS9_mu_sub.push_back(muon_sub.calEnergy().emS9);
            emS25_mu_sub.push_back(muon_sub.calEnergy().emS25);
            emMax_mu_sub.push_back(muon_sub.calEnergy().emMax);
            had_mu_sub.push_back(muon_sub.calEnergy().had);
            hadS9_mu_sub.push_back(muon_sub.calEnergy().hadS9);
            hadMax_mu_sub.push_back(muon_sub.calEnergy().hadMax);
            ho_mu_sub.push_back(muon_sub.calEnergy().ho);
            hoS9_mu_sub.push_back(muon_sub.calEnergy().hoS9);
            numOfMatchedStations_mu_sub.push_back(muon_sub.numberOfMatchedStations());
            numOfChambers_mu_sub.push_back(muon_sub.numberOfChambers());
            segCompatibility_mu_sub.push_back(muon::segmentCompatibility(muon_sub));
            caloCompatibility_mu_sub.push_back(muon_sub.caloCompatibility());
            chi2LocalPosition_mu_sub.push_back(muon_sub.combinedQuality().chi2LocalPosition);
            trkKink_mu_sub.push_back(muon_sub.combinedQuality().trkKink);
            glbKink_mu_sub.push_back(muon_sub.combinedQuality().glbKink);
            float pfrelIso04_val_ms = -999.f;
            if (muon_sub.isPFMuon() && muon_sub.isPFIsolationValid()) {
              auto iso04ms = muon_sub.pfIsolationR04();
              pfrelIso04_val_ms = ( iso04ms.sumChargedHadronPt
                                  + std::max(0.f, iso04ms.sumNeutralHadronEt + iso04ms.sumPhotonEt - 0.5f*iso04ms.sumPUPt) )
                                  / muon_sub.pt();
            }
            pfrelIso04_mu_sub.push_back(pfrelIso04_val_ms);
            const auto& glbTrkMS = muon_sub.globalTrack();
            const auto& innTrkMS = muon_sub.innerTrack();
            const auto& outTrkMS = muon_sub.outerTrack();
            const auto& bestTrkMS = muon_sub.muonBestTrack();
            const auto& tunePTrkMS = muon_sub.tunePMuonBestTrack();
            nValidMuonHits_mu_sub.push_back(glbTrkMS.isNonnull() ? glbTrkMS->hitPattern().numberOfValidMuonHits() : -999);
            normChi2_global_mu_sub.push_back(glbTrkMS.isNonnull() ? (float)glbTrkMS->normalizedChi2() : -999.f);
            nValidPixelHits_mu_sub.push_back(innTrkMS.isNonnull() ? innTrkMS->hitPattern().numberOfValidPixelHits() : -999);
            trkLayers_mu_sub.push_back(innTrkMS.isNonnull() ? innTrkMS->hitPattern().trackerLayersWithMeasurement() : -999);
            pixelLayers_mu_sub.push_back(innTrkMS.isNonnull() ? innTrkMS->hitPattern().pixelLayersWithMeasurement() : -999);
            normChi2_inner_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->normalizedChi2() : -999.f);
            validFraction_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->validFraction() : -999.f);
            tunePtErrorOverPt_mu_sub.push_back(tunePTrkMS.isNonnull() && tunePTrkMS->pt() > 0.f ? (float)(tunePTrkMS->ptError()/tunePTrkMS->pt()) : -999.f);
            dxy_PV_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->dxy(pv->position()) : -999.f);
            dz_PV_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->dz(pv->position()) : -999.f);
            nMatches_mu_sub.push_back(muon_sub.numberOfMatches(reco::Muon::SegmentAndTrackArbitration));
            innerTrk_pt_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->pt() : -999.f);
            innerTrk_eta_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->eta() : -999.f);
            innerTrk_phi_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->phi() : -999.f);
            outerTrk_pt_mu_sub.push_back(outTrkMS.isNonnull() ? (float)outTrkMS->pt() : -999.f);
            outerTrk_eta_mu_sub.push_back(outTrkMS.isNonnull() ? (float)outTrkMS->eta() : -999.f);
            outerTrk_phi_mu_sub.push_back(outTrkMS.isNonnull() ? (float)outTrkMS->phi() : -999.f);
            globalTrk_pt_mu_sub.push_back(glbTrkMS.isNonnull() ? (float)glbTrkMS->pt() : -999.f);
            globalTrk_eta_mu_sub.push_back(glbTrkMS.isNonnull() ? (float)glbTrkMS->eta() : -999.f);
            globalTrk_phi_mu_sub.push_back(glbTrkMS.isNonnull() ? (float)glbTrkMS->phi() : -999.f);
            bestTrk_pt_mu_sub.push_back(bestTrkMS.isNonnull() ? (float)bestTrkMS->pt() : -999.f);
            bestTrk_eta_mu_sub.push_back(bestTrkMS.isNonnull() ? (float)bestTrkMS->eta() : -999.f);
            bestTrk_phi_mu_sub.push_back(bestTrkMS.isNonnull() ? (float)bestTrkMS->phi() : -999.f);
            tunePTrk_pt_mu_sub.push_back(tunePTrkMS.isNonnull() ? (float)tunePTrkMS->pt() : -999.f);
            tunePTrk_eta_mu_sub.push_back(tunePTrkMS.isNonnull() ? (float)tunePTrkMS->eta() : -999.f);
            tunePTrk_phi_mu_sub.push_back(tunePTrkMS.isNonnull() ? (float)tunePTrkMS->phi() : -999.f);
            muon_charge_mu_sub.push_back(muon_sub.charge());
            innerTrk_charge_mu_sub.push_back(innTrkMS.isNonnull() ? innTrkMS->charge() : -999);
            outerTrk_charge_mu_sub.push_back(outTrkMS.isNonnull() ? outTrkMS->charge() : -999);
            globalTrk_charge_mu_sub.push_back(glbTrkMS.isNonnull() ? glbTrkMS->charge() : -999);
            bestTrk_charge_mu_sub.push_back(bestTrkMS.isNonnull() ? bestTrkMS->charge() : -999);
            tunePTrk_charge_mu_sub.push_back(tunePTrkMS.isNonnull() ? tunePTrkMS->charge() : -999);
            innerTrk_ptError_mu_sub.push_back(innTrkMS.isNonnull() ? (float)innTrkMS->ptError() : -999.f);
            outerTrk_ptError_mu_sub.push_back(outTrkMS.isNonnull() ? (float)outTrkMS->ptError() : -999.f);
            globalTrk_ptError_mu_sub.push_back(glbTrkMS.isNonnull() ? (float)glbTrkMS->ptError() : -999.f);
            bestTrk_ptError_mu_sub.push_back(bestTrkMS.isNonnull() ? (float)bestTrkMS->ptError() : -999.f);
            tunePTrk_ptError_mu_sub.push_back(tunePTrkMS.isNonnull() ? (float)tunePTrkMS->ptError() : -999.f);
            dxy_bestTrk_PV_mu_sub.push_back(bestTrkMS.isNonnull() ? (float)bestTrkMS->dxy(pv->position()) : -999.f);
            dz_bestTrk_PV_mu_sub.push_back(bestTrkMS.isNonnull() ? (float)bestTrkMS->dz(pv->position()) : -999.f);

            // 예시 출력
            //std::cout << "Invariant mass: " << invM << std::endl;
          }
        }
        //if (invM > 100) std::cout<<"wtf?"<<std::endl;
        
        // muon 기준 출력 등
        //std::cout << muon.pt() << " | " << muon.eta() << " | " << muon.phi() << " | " << close_muon << std::endl;
      }
      weight_mu = aWeight; 
      //if (aWeight < 0)std::cout<<aWeight<<std::endl;
      muon_->Fill();
    }
  }
  
  return;
}

DEFINE_FWK_MODULE(MuonJpsi);
