#ifndef Minituplizer_H
#define Minituplizer_H

#include <TROOT.h>
#include <TMath.h>
#include <vector>
#include <string>

class NtupleTrigger {
public:
  bool isFired;
  std::string name;

  NtupleTrigger(){};
  virtual ~NtupleTrigger(){};
};

class NtupleMuon {
public:
  bool isGlobalMuon;
  bool isStandAloneMuon;
  bool isTrackerMuon;
  bool isPFMuon;
  double px;
  double py;
  double pz;
  double pt;
  double phi;
  double eta;
  int charge;
  int nChambers;
  int stationMask;
  int nMatchedStations;
  double normalizedChi2;
  int nValidHits;
  int nValidTrackerHits;
  int nValidPixelHits;
  int nTrackerLayers;
  int nValidMuonHits;
  double qoverp;
  double theta;
  double lambda;
  double dxy;
  double d0;
  double dsz;
  double dz;
  double dxyBS;
  double dzBS;
  double dszBS;
  double dxyVTX;
  double dzVTX;
  double dszVTX;
  double vx;
  double vy;
  double vz;

  // muon tracks
  double muonBestTrack_px;
  double muonBestTrack_py;
  double muonBestTrack_pz;
  double muonBestTrack_pt;
  double muonBestTrack_ptError;
  int muonBestTrack_nValidPixelHits;
  int muonBestTrack_nTrackerLayers;
  int muonBestTrack_nValidMuonHits;
  double muonBestTrack_dxyVTX;
  double muonBestTrack_dzVTX;
  double muonBestTrack_dszVTX;
  double tunePMuonBestTrack_px;
  double tunePMuonBestTrack_py;
  double tunePMuonBestTrack_pz;
  double tunePMuonBestTrack_pt;
  double tunePMuonBestTrack_ptError;
  int tunePMuonBestTrack_nValidPixelHits;
  int tunePMuonBestTrack_nTrackerLayers;
  int tunePMuonBestTrack_nValidMuonHits;
  double tunePMuonBestTrack_dxyVTX;
  double tunePMuonBestTrack_dzVTX;
  double tunePMuonBestTrack_dszVTX;
  double innerTrack_px;
  double innerTrack_py;
  double innerTrack_pz;
  double innerTrack_pt;
  double innerTrack_ptError;
  int innerTrack_nValidPixelHits;
  int innerTrack_nTrackerLayers;
  double innerTrack_dxyVTX;
  double innerTrack_dzVTX;
  double innerTrack_dszVTX;
  double outerTrack_px;
  double outerTrack_py;
  double outerTrack_pz;
  double outerTrack_pt;
  double outerTrack_ptError;
  int outerTrack_nValidMuonHits;

  // muon isolation
  double miniChIso;
  double miniNhIso;
  double miniPhIso;
  double miniPUIso;
  double isolationR03_sumpt;
  double isolationR03_hadEt;
  double isolationR03_emEt;
  double isolationR05_sumpt;
  double isolationR05_hadEt;
  double isolationR05_emEt;
  double PfChargedHadronIsoR05;
  double PfNeutralHadronIsoR05;
  double PfGammaIsoR05;
  double PfChargedHadronIsoR04;
  double PfNeutralHadronIsoR04;
  double PfGammaIsoR04;
  double PfPUPtR04;
  double PfChargedHadronIsoR03;
  double PfNeutralHadronIsoR03;
  double PfGammaIsoR03;
  double PfPUPtR03;

  NtupleMuon(){};
  virtual ~NtupleMuon(){};
};

class NtupleElectron {
public:
  double pt;
  double eta;
  double rap;
  double phi;
  double en;
  int charge;
  double enSC;
  double preEnSC;
  double rawEnSC;
  double etSC;
  double etaSC;
  double phiSC;
  double etaSCWidth;
  double phiSCWidth;
  double full5x5_sigmaIetaIeta;
  double full5x5_sigmaIphiIphi;
  double full5x5_E1x5;
  double full5x5_E2x5;
  double full5x5_E5x5;
  double full5x5_hOverE;
  double full5x5_hOverEBC;
  double full5x5_r9;
  double dEtaIn;
  double dPhiIn;
  double dEtaSeed;
  double dPhiSeed;
  double dEtaEle;
  double dPhiEle;
  double EOverP;
  double EseedOverP;
  double EseedOverPout;
  double EeleOverPout;
  double pIn;
  double pOut;
  double ptIn;
  double ptOut;
  double ecalEn;
  double ecalErr;
  double trkErr;
  double combErr;
  double PFcombErr;
  double ooEmooP;
  double isoChargedHadrons;
  double isoNeutralHadrons;
  double isoPhotons;
  double isoChargedFromPU;
  double miniChIso;
  double miniNhIso;
  double miniPhIso;
  double miniPUIso;
  double dr03EcalRecHitSumEt;
  double dr03HcalDepth1TowerSumEt;
  bool passHEEP;
  bool passHEEPTrkIso;
  bool passHEEPEMHad1Iso;
  bool passHEEPEt;
  bool passHEEPdEtaInSeed;
  bool passHEEPdPhiIn;
  bool passHEEPE2x5OverE5x5;
  bool passHEEPHoverE;
  bool passHEEPdxy;
  bool passHEEPMissingHits;
  bool passHEEPEcalDriven;
  // bool passN1HEEPTrkIso;
  double HEEPTrkIso;
  double HEEPEMHad1Iso;
  int HEEPnrSatCrysValue;
  double HEEPTrkIsoMod;
  bool HEEPaddGsfTrkSel;
  double HEEPEcalRecHitIsoMod;

  int HEEPnoSelectedGsfTrk;
  double HEEPaddGsfTrk_Gsfpt;
  double HEEPaddGsfTrk_Gsfeta;
  double HEEPaddGsfTrk_Gsfphi;
  double HEEPaddGsfTrk_GsfptErr;
  int HEEPaddGsfTrk_lostHits;
  int HEEPaddGsfTrk_nValidHits;
  int HEEPaddGsfTrk_nValidPixelHits;
  double HEEPaddGsfTrk_chi2;
  double HEEPaddGsfTrk_d0;
  double HEEPaddGsfTrk_d0Err;
  double HEEPaddGsfTrk_dxy;
  double HEEPaddGsfTrk_dxyErr;
  double HEEPaddGsfTrk_dz;
  double HEEPaddGsfTrk_vz;
  double HEEPaddGsfTrk_dzErr;
  double HEEPaddGsfTrk_dsz;
  double HEEPaddGsfTrk_dszErr;

  double HEEPaddVtx_isValid;
  double HEEPaddVtx_dx;
  double HEEPaddVtx_dy;
  double HEEPaddVtx_dz;
  double HEEPaddVtx_chi2;
  double HEEPaddVtx_xErr;
  double HEEPaddVtx_yErr;
  double HEEPaddVtx_zErr;
  double HEEPaddVtx_pt;
  double HEEPaddVtx_rapidity;
  double HEEPaddVtx_phi;
  double HEEPaddVtx_M;
  int lostHits;
  int validHits;
  int nValidHits;
  int nValidPixelHits;
  double chi2;
  int GsfHits;
  int GsfCharge;
  double d0;
  double d0Err;
  double dxy;
  double dxyErr;
  double dz;
  double vz;
  double dzErr;
  double dsz;
  double dszErr;
  double Gsfpt;
  double Gsfeta;
  double Gsfphi;
  double GsfptErr;
  double Gsfpx;
  double Gsfpy;
  double Gsfpz;
  double KFchi2;
  int KFvalidHits;
  int KFpixHits;
  int KFstripHits;
  double KFpt;
  double KFeta;
  double KFphi;
  int KFcharge;
  double KFptErr;
  double KFpx;
  double KFpy;
  double KFpz;
  double KFd0;
  double KFdxy;
  double KFdz;
  double KFdsz;
  double expectedMissingInnerHits;
  bool passConversionVeto;
  double convVtxFitProb;
  double convVtxChi2;
  double convDist;
  double convDcot;
  double convRadius;
  double fbrem;
  double fbremSC;
  int nbrem;
  bool isEB;
  bool isEE;
  bool isEBEEGap;
  bool isEBEtaGap;
  bool isEBPhiGap;
  bool isEEDeeGap;
  bool isEERingGap;

  bool ecalDrivenSeed;
  bool trackerDrivenSeed;

  double mvaValue;
  int mvaCategory;

  NtupleElectron(){};
  virtual ~NtupleElectron(){};
};

class NtuplePhoton {
public:
  bool hasPixelSeed;
  double pt;
  double eta;
  double phi;
  double etaSC;
  double phiSC;
  double HoverE;
  double Full5x5_SigmaIEtaIEta;
  double ChIso;
  double NhIso;
  double PhIso;

  NtuplePhoton(){};
  virtual ~NtuplePhoton(){};
};

class NtupleJet {
public:
  float pt;
  float eta;
  float phi;
  float et;
  float energy;
  int charge;
  float chHadEn;
  float chHadFrac;
  float neHadEn;
  float neHadFrac;
  float phEn;
  float phFrac;
  float elEn;
  float elFrac;
  float muEn;
  float muFrac;
  int chHadMulti;
  int neHadMulti;
  int phMulti;
  int elMulti;
  int muMulti;
  float chEmEn;
  float chEmFrac;
  float chMuEn;
  float chMuFrac;
  float neEmEn;
  float neEmFrac;
  int chMulti;
  int neMulti;

  NtupleJet(){};
  virtual ~NtupleJet(){};
};

class NtupleMET {
public:
  double pt;
  double phi;
  double px;
  double py;
  double sumEt;
  double Type1_pt;
  double Type1_phi;
  double Type1_px;
  double Type1_py;
  double Type1_sumEt;

  NtupleMET(){};
  virtual ~NtupleMET(){};
};

class NtupleGenParticle {
public:
  double px;
  double py;
  double pz;
  double pt;
  double phi;
  double eta;
  double energy;
  double mass;
  int mother;
  double motherPt;
  int charge;
  int status;
  int id;
  bool fromHardProcessFinalState;
  bool fromHardProcessDecayed;
  bool fromHardProcessBeforeFSR;
  bool isHardProcess;
  bool isLastCopy;
  bool isLastCopyBeforeFSR;
  bool isPromptDecayed;
  bool isPromptFinalState;
  bool isDirectHardProcessTauDecayProductFinalState;
  bool isDirectPromptTauDecayProductFinalState;

  NtupleGenParticle(){};
  virtual ~NtupleGenParticle(){};
};

class NtupleGenJet {
public:
  float pt;
  float eta;
  float phi;
  float et;
  float energy;
  int charge;
  int nConst;
  float emEnergy;
  float hadEnergy;
  float invisibleEnergy;
  float auxiliaryEnergy;

  NtupleGenJet(){};
  virtual ~NtupleGenJet(){};
};

class NtupleGsfTrack {
public:
  int lostHits;
  int validHits;
  int nValidHits;
  int nValidPixelHits;
  double chi2;
  int GsfHits;
  int GsfCharge;
  double d0;
  double d0Err;
  double dxy;
  double dxyErr;
  double dz;
  double vz;
  double dzErr;
  double dsz;
  double dszErr;
  double Gsfpt;
  double Gsfeta;
  double Gsfphi;
  double GsfptErr;
  double Gsfpx;
  double Gsfpy;
  double Gsfpz;

  NtupleGsfTrack(){};
  virtual ~NtupleGsfTrack(){};
};

class NtuplePackedCand {
public:
  double p;
  double en;
  double et;
  double mass;
  double mt;
  double pt;
  double eta;
  double phi;
  double vz;
  double dxy;
  double dz;
  bool hasTrackDetails;
  bool trackHighPurity;
  int pdgId;
  double vertexChi2;
  double vertexNdof;
  double vertexNormalizedChi2;
  double trkPt;
  double trkEta;
  double trkPhi;
  double trkVz;
  int nValidHits;
  int nValidPixelHits;
  double trkPtError;
  double trkChi2;
  double trkNdof;
  double trkNormalizedChi2;
  int trkCharge;
  double trkDxy;
  double trkD0;
  double trkDz;
  double qoverp;
  bool isHighPurity;
  bool isJetCoreRegionalAlgo;
  int qualMask;
  int algoMask;

  NtuplePackedCand(){};
  virtual ~NtuplePackedCand(){};
};

class NtupleEvent {
public:
  int run;
  unsigned long long event;
  int lumi;
  int nVertices;
  int nPU;
  int nPUin;
  int nMuons;
  int nElectrons;
  int nPhotons;
  int nJets;
  int nMETs;
  int nGsfTrks;
  int nTrks;
  double rho;
  double weight;
  int nGenParticles;
  int nGenJets;

  std::vector<NtupleTrigger> triggers;
  std::vector<NtupleMuon> muons;
  std::vector<NtupleElectron> electrons;
  std::vector<NtuplePhoton> photons;
  std::vector<NtupleJet> jets;
  std::vector<NtupleMET> METs;
  std::vector<NtupleGenParticle> genparticles;
  std::vector<NtupleGenJet> genjets;
  std::vector<NtupleGsfTrack> GsfTracks;
  std::vector<NtuplePackedCand> lostTracks;
  std::vector<NtuplePackedCand> packedPFCands;

  NtupleEvent(){};
  virtual ~NtupleEvent(){};
};

#endif
