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
    float enSC, etSC, etaSC, phiSC, etaSCWidth, phiSCWidth;
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
    float convVtxFitProb, convVtxChi2, convDist, convDcot, convRadius;
    int passConversionVeto, nbrem;
    float fbrem, fbremSC;
    float full5x5_e2x5Left, full5x5_e2x5Right, full5x5_e2x5Top, full5x5_e2x5Bottom, full5x5_eLeft, full5x5_eRight, full5x5_eTop, full5x5_eBottom;
    float full5x5_eMax, full5x5_e2nd;
  } ElectronStruct;

  typedef struct {
    float weight;
    float Gsfpt, Gsfeta, Gsfphi;
    int lostHits, nValidHits, nValidPixelHits;
    float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz, ptErr;
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

  void fillElectrons(const edm::Ptr<reco::GsfElectron>& el,
                     const float& trkIso,
                     const float& ecalIso,
                     const reco::GsfTrackRef& addGsfTrk,
                     const edm::Handle<reco::ConversionCollection>& conversions,
                     const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                     const std::string& prefix);

  void fillGsfTracks(const reco::GsfTrackRef& addGsfTrk,
                     const reco::GsfTrackRef& orgGsfTrk,
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

private:
  std::map<std::string,TTree*> tree_;
  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,ElectronStruct> elvalues_;
  std::map<std::string,AddGsfStruct> gsfvalues_;
  std::map<std::string,MuonStruct> muonvalues_;
  std::map<std::string,METstruct> metvalues_;

  edm::Service<TFileService>* pFS_;
  edm::Ptr<reco::Vertex> pPV_;

  double mcweight_;

  TString elstr_;
  TString addgsfstr_;
  TString mustr_;
  TString metstr_;
};

#endif
