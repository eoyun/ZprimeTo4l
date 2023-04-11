#ifndef MergedLeptonHelper_H
#define MergedLeptonHelper_H 1

#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"

#include "TTree.h"
#include "TString.h"
#include "TMath.h"

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
    float Gsfpt, Gsfeta, Gsfphi;
    int expectedMissingInnerHits;
    float convDist, convDcot, convRadius;
    int passConversionVeto, nbrem;
    float fbrem, fbremSC;
    float full5x5_e2x5Left, full5x5_e2x5Right, full5x5_e2x5Top, full5x5_e2x5Bottom, full5x5_eLeft, full5x5_eRight, full5x5_eTop, full5x5_eBottom;
    float full5x5_eMax, full5x5_e2nd;
    float etaE1st, phiE1st, etaE2nd, phiE2nd, etaSeedLinear, phiSeedLinear;
    float clus2ndMoment_sMaj, clus2ndMoment_sMin, clus2ndMoment_alpha;
    float full5x5_Em2m2, full5x5_Em2m1, full5x5_Em2p0, full5x5_Em2p1, full5x5_Em2p2;
    float full5x5_Em1m2, full5x5_Em1m1, full5x5_Em1p0, full5x5_Em1p1, full5x5_Em1p2;
    float full5x5_Ep0m2, full5x5_Ep0m1, full5x5_Ep0p0, full5x5_Ep0p1, full5x5_Ep0p2;
    float full5x5_Ep1m2, full5x5_Ep1m1, full5x5_Ep1p0, full5x5_Ep1p1, full5x5_Ep1p2;
    float full5x5_Ep2m2, full5x5_Ep2m1, full5x5_Ep2p0, full5x5_Ep2p1, full5x5_Ep2p2;
  } ElectronStruct;

  typedef struct {
    float weight;
    float Gsfpt, Gsfeta, Gsfphi;
    int lostHits, nValidHits, nValidPixelHits;
    float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz, ptErr;
    float deltaEtaSuperClusterAtVtx, deltaPhiSuperClusterAtVtx;
    float deltaEtaSeedClusterAtCalo, deltaPhiSeedClusterAtCalo;
    float deltaEtaSeedClusterAtVtx;
    float union3x3LogEta, union3x3LogPhi;
    float union3x3LinearEta, union3x3LinearPhi;
    float union3x3Energy;
  } AddGsfStruct;

  typedef struct {
    float weight;
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
    float weight;
    float PFphi, PFdPhi, PFpt, PFoGen, PFoMu, MuEta, PFSumEt, TPphi, TPdPhi, TPpt, TPoGen, TPoMu, TPSumEt, dSumEt, sumEtRatio, sumEtRatioTP;
    int nLepton;
    unsigned int bitmask;
  } METstruct;

  MergedLeptonHelper();
  virtual ~MergedLeptonHelper() {}

public:
  void initElectronTree(const std::string& name,
                        const std::string& prefix,
                        const std::string& postfix);
  void initAddGsfTree(const std::string& name,
                      const std::string& prefix,
                      const std::string& postfix);
  void initMuonTree(const std::string& name,
                    const std::string& prefix,
                    const std::string& postfix);
  void initMETTree(const std::string& name,
                   const std::string& prefix,
                   const std::string& postfix);

  void fillElectrons(const pat::ElectronRef& el,
                     const float& trkIso,
                     const float& ecalIso,
                     const reco::GsfTrackRef& addGsfTrk,
                     const edm::EventSetup& iSetup,
                     const EcalRecHitCollection* ecalRecHits,
                     const std::string& prefix);

  void fillGsfTracks(const pat::ElectronRef& el,
                     const reco::GsfTrackRef& addGsfTrk,
                     const edm::EventSetup& iSetup,
                     const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                     const EcalRecHitCollection* ecalRecHits,
                     const std::string& prefix);

  void fillMuons(const edm::Ptr<reco::Muon>&,
                 const edm::Ptr<reco::GenParticle>&,
                 const edm::Ptr<reco::GenParticle>&,
                 std::string);

  void fillMETs(const edm::Ptr<reco::Muon>&,
                const edm::Ptr<reco::GenParticle>&,
                const pat::MET&,
                const std::vector<edm::Ptr<reco::Muon>>&,
                const unsigned int&,
                const std::string&);
  void fillMETs(const edm::Ptr<reco::Candidate>&,
                const pat::MET&,
                const std::vector<edm::Ptr<reco::Candidate>>&,
                const unsigned int&,
                const std::string&);

  // non of these are owned by the helper
  void SetFileService(edm::Service<TFileService>* fs) { pFS_ = fs; }
  void SetPV(edm::Ptr<reco::Vertex> pv) { pPV_ = pv; }
  void SetMCweight(double mcweight) { mcweight_ = mcweight; }
  void SetPositionCalcLog(PositionCalc& calc) { posCalcLog_ = calc; }
  void SetPositionCalcLinear(PositionCalc& calc) { posCalcLinear_ = calc; }

private:
  std::map<std::string,TTree*> tree_;
  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,ElectronStruct> elvalues_;
  std::map<std::string,AddGsfStruct> gsfvalues_;
  std::map<std::string,MuonStruct> muonvalues_;
  std::map<std::string,METstruct> metvalues_;

  edm::Service<TFileService>* pFS_;
  edm::Ptr<reco::Vertex> pPV_;

  // Ecal position calculation algorithm
  PositionCalc posCalcLog_;
  PositionCalc posCalcLinear_;

  double mcweight_;

  TString elstr_;
  TString addgsfstr_;
  TString mustr_;
  TString metstr_;
};

#endif
