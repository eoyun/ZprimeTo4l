#ifndef MergedLeptonHelper_H
#define MergedLeptonHelper_H 1

#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedDEtaInSeed.h"
#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedShowerShape.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "Geometry/CaloTopology/interface/CaloTopology.h"
//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"

#include "TTree.h"
#include "TString.h"
#include "TMath.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"

namespace MergedLeptonHelperFct {
  bool isNotMerged(const pat::ElectronRef& aEle,
                   const edm::Handle<edm::View<pat::Electron>>& eleHandle,
                   const reco::GsfTrackRef& addGsfTrk);
}

class MergedLeptonHelper {
public:
  typedef struct {
    float weight;
    float pt, eta, phi, en, et;
    int charge;
    float enSC, etSC, etaSC, phiSC, etaSeed, phiSeed, etaSCWidth, phiSCWidth;
    float full5x5_sigmaIetaIeta, full5x5_sigmaIphiIphi;
    float full5x5_E1x5, full5x5_E2x5, full5x5_E5x5, full5x5_hOverE, full5x5_r9;
    float dEtaIn, dPhiIn, dPhiSeed, dEtaEle, dPhiEle, dEtaSeed;
    float EseedOverP, EOverP;
    float ecalEn, ecalErr, trkErr, combErr, PFcombErr;
    float dr03TkSumPtHEEP, dr03EcalRecHitSumEt, dr03HcalDepth1TowerSumEt;
    float modTrkIso, modEcalIso;
    int lostHits, nValidHits, nValidPixelHits, GsfHits;
    float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz;
    float Gsfpt, Gsfeta, Gsfphi, GsfPtErr;
    int expectedMissingInnerHits;
    float convDist, convDcot, convRadius;
    int passConversionVeto, nbrem;
    float fbrem, fbremSC;
    float full5x5_e2x5Left, full5x5_e2x5Right, full5x5_e2x5Top, full5x5_e2x5Bottom;
    float full5x5_eLeft, full5x5_eRight, full5x5_eTop, full5x5_eBottom;
    float full5x5_eMax, full5x5_e2nd;
    float clus2ndMoment_sMaj, clus2ndMoment_sMin, clus2ndMoment_alpha;
    float dEtaSeedMiniAOD, dPhiInMiniAOD, sigIeIeMiniAOD;
    float union5x5dEtaIn, union5x5dPhiIn;
    float union5x5Energy, union5x5covIeIe, union5x5covIeIp, union5x5covIpIp;
    float union5x5covMaj, union5x5covMin;
    float alphaCalo;
    float GenPt, GenE;
    int u5x5numGood, u5x5numPoorReco, u5x5numOutOfTime, u5x5numFaultyHardware;
    int u5x5numNoisy, u5x5numPoorCalib, u5x5numSaturated, u5x5numLeadingEdgeRecovered;
    int u5x5NeighboursRecovered, u5x5numTowerRecovered, u5x5numDead, u5x5numKilled;
    int u5x5numTPSaturated, u5x5numL1SpikeFlag, u5x5numWeird, u5x5numDiWeird;
    int u5x5numHasSwitchToGain6, u5x5numHasSwitchToGain1, u5x5numUnknown;
    int u5x5numTPSaturatedAndTowerRecovered;
  } ElectronStruct;

  typedef struct {
    float weight;
    float pt, eta, phi;
    int lostHits, nValidHits, nValidPixelHits;
    float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz, ptErr;
    float deltaEtaSeedClusterAtVtx, deltaPhiSuperClusterAtVtx;
    float dPerpIn, normalizedDParaIn, alphaTrack;
  } AddTrkStruct;

  MergedLeptonHelper(edm::ConsumesCollector iC);
  virtual ~MergedLeptonHelper()=default;

public:
  void initElectronTree(const std::string& name,
                        const std::string& prefix,
                        const std::string& postfix);
  void initAddTrkTree(const std::string& name,
                      const std::string& prefix,
                      const std::string& postfix);

  void fillElectrons(const pat::ElectronRef& el,
                     const float& trkIso,
                     const float& ecalIso,
                     const ModifiedDEtaInSeed::variables& variablesDEtaIn,
                     const ModifiedShowerShape::variables& variables,
                     const EcalRecHitCollection* ecalRecHits,
                     const edm::EventSetup& iSetup,
                     const std::string& prefix,
		     edm::ConsumesCollector iC,
                     const float genPt = -1., // signal MC only
                     const float genE = -1.);

  void fillAddTracks(const pat::ElectronRef& el,
                     const reco::TrackBase* addTrk,
                     const ModifiedDEtaInSeed::variables& variables,
                     const EcalRecHitCollection* ecalRecHits,
                     const edm::EventSetup& iSetup,
                     const std::string& prefix);

  // non of these are owned by the helper
  void SetFileService(edm::Service<TFileService>* fs) { pFS_ = fs; }
  void SetPV(const reco::VertexRef& pv) { pPV_ = pv; }
  void SetBS(const reco::BeamSpot* bs) { pBS_ = bs; }
  void SetMCweight(double mcweight) { mcweight_ = mcweight; }
  void SetPositionCalcLog(PositionCalc& calc) { posCalcLog_ = calc; }

private:
  std::map<std::string,TTree*> tree_;
  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,ElectronStruct> elvalues_;
  std::map<std::string,AddTrkStruct> trkvalues_;
  
  edm::ESGetToken<CaloGeometry, CaloGeometryRecord> geometryToken_;
  edm::ESGetToken<CaloTopology, CaloTopologyRecord> topologyToken_;

  edm::Service<TFileService>* pFS_;
  reco::VertexRef pPV_;
  const reco::BeamSpot* pBS_;

  // Ecal position calculation algorithm
  PositionCalc posCalcLog_;

  double mcweight_;

  TString elstr_;
  TString addtrkstr_;
};

#endif
