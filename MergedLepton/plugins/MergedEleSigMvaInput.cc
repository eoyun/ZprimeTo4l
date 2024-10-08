#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

// produce TTree for merged electron training with H->AA->4e events

class MergedEleSigMvaInput : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleSigMvaInput(const edm::ParameterSet&);
  virtual ~MergedEleSigMvaInput() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  void fillByGsfTrack(const edm::Event& iEvent,
                      const edm::EventSetup& iSetup,
                      const std::vector<pat::ElectronRef>& recoPair,
                      const std::vector<reco::GenParticleRef>& promptPair);

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dEtaInSeed2ndToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPhiInSC2ndToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaTrackToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaCaloToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> normDParaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5covIeIeToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5covIeIpToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5covIpIpToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dEtaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dPhiInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5EnergyToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EBrecHitToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EErecHitToken_;
  // Ecal position calculation algorithm
  PositionCalc posCalcLog_;

  const double ptThres_;
  const double ptThres2nd_;
  const double drThres_;
  const double drThres2_;
  const double drThres2nd_;
  const double drThres2nd2_;
  const bool allowNonPrompt_;

  std::map<std::string,TH1*> histo1d_;

  MergedLeptonHelper aHelper_;
};

MergedEleSigMvaInput::MergedEleSigMvaInput(const edm::ParameterSet& iConfig)
: srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
  srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
  pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
  trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
  ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
  dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
  dEtaInSeed2ndToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dEtaInSeed2nd"))),
  dPhiInSC2ndToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPhiInSC2nd"))),
  alphaTrackToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaTrack"))),
  alphaCaloToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaCalo"))),
  normDParaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("normalizedDParaIn"))),
  union5x5covIeIeToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5covIeIe"))),
  union5x5covIeIpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5covIeIp"))),
  union5x5covIpIpToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5covIpIp"))),
  union5x5dEtaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dEtaIn"))),
  union5x5dPhiInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dPhiIn"))),
  union5x5EnergyToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5Energy"))),
  addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
  addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
  beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
  EBrecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBrecHits"))),
  EErecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EErecHits"))),
  posCalcLog_(PositionCalc(iConfig.getParameter<edm::ParameterSet>("posCalcLog"))),
  ptThres_(iConfig.getParameter<double>("ptThres")),
  ptThres2nd_(iConfig.getParameter<double>("ptThres2nd")),
  drThres_(iConfig.getParameter<double>("drThres")),
  drThres2_(drThres_*drThres_),
  drThres2nd_(iConfig.getParameter<double>("drThres2nd")),
  drThres2nd2_(drThres2nd_*drThres2nd_),
  allowNonPrompt_(iConfig.getParameter<bool>("allowNonPrompt")),
  aHelper_(consumesCollector()) {
  usesResource("TFileService");
}

void MergedEleSigMvaInput::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  aHelper_.SetFileService(&fs);
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  aHelper_.initElectronTree("ElectronStruct","mergedEl1","el");
  aHelper_.initElectronTree("ElectronStruct","mergedEl2","el");
  aHelper_.initAddTrkTree("AddTrkStruct","mergedEl1Trk","addTrk");

  aHelper_.SetPositionCalcLog(posCalcLog_);
}

void MergedEleSigMvaInput::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  if (eleHandle->empty())
    return;

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken_, genInfo);
  double mcweight = genInfo->weight();

  double aWeight = mcweight/std::abs(mcweight);
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  aHelper_.SetMCweight(mcweight);

  reco::VertexRef primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = pvHandle->refAt(0).castTo<reco::VertexRef>();
    aHelper_.SetPV(primaryVertex);
  } else
    return;

  aHelper_.SetBS(beamSpotHandle.product());

  std::vector<reco::GenParticleRef> promptLeptons;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    const auto& genptc = genptcHandle->refAt(idx);

    if ( ( std::abs(genptc->pdgId())==11 ) &&
         genptc->fromHardProcessFinalState() &&
         ( std::abs(genptc->eta()) < 2.5 ) )
      promptLeptons.push_back(genptc.castTo<reco::GenParticleRef>());

    if (allowNonPrompt_) {
      if ( ( std::abs(genptc->pdgId())==11 ) &&
           genptc->statusFlags().isDirectHadronDecayProduct() &&
           ( std::abs(genptc->eta()) < 2.5 ) &&
           ( std::abs(genptc->mother()->pdgId())==443 ) )
        promptLeptons.push_back(genptc.castTo<reco::GenParticleRef>());
    }
  }

  if (promptLeptons.size() < 2)
    return;

  // sort by mother ptc, and then by et (motherA Et1, motherA Et2, motherB Et1, motherB Et2)
  auto sisterLambda = [](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
    return (a->mother() == b->mother()) ? (a->et() > b->et()) : (a->mother() < b->mother());
  };

  std::sort(promptLeptons.begin(),promptLeptons.end(),sisterLambda);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  std::vector<pat::ElectronRef> heeps1, heeps2;

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->refAt(idx);

    const float eta1stGSF = -(aEle->deltaEtaSeedClusterTrackAtVtx() - aEle->superCluster()->seed()->eta());
    const float u5x5Eta = (*union5x5dEtaInHandle)[aEle] + eta1stGSF;
    const float u5x5Et = (*union5x5EnergyHandle)[aEle]/std::cosh(u5x5Eta);

    if ( u5x5Et < ptThres_ )
      continue;

    bool matched = false;
    size_t igen = 0;

    // for given reco electron, check whether GEN-lv prompt electron exists within dR thres
    for (; igen < promptLeptons.size(); ++igen) {
      const auto& genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aEle->gsfTrack()->eta(),aEle->gsfTrack()->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2_) {
        matched = true;
        break;
      }
    }

    if ( !matched )
      continue;

    auto castEle = aEle.castTo<pat::ElectronRef>();

    (igen < 2) ? heeps1.push_back(castEle) : heeps2.push_back(castEle);
  }

  // for eemm & eeee final states
  if ( heeps1.empty() )
    return;

  if ( promptLeptons.at(0)->pt() + promptLeptons.at(1)->pt() < ptThres_ )
    return;

  if ( !allowNonPrompt_) {
    if ( promptLeptons.at(0)->pt() < ptThres_ || promptLeptons.at(1)->pt() < ptThres2nd_ )
      return;

    if ( heeps1.size() > 1 && promptLeptons.at(1)->pt() < ptThres_ )
      return;
  }

  // for eeee final state
  if ( promptLeptons.size() > 2 ) {
    if (promptLeptons.size() != 4)
      return;

    if ( heeps2.empty() )
      return;

    if ( promptLeptons.at(2)->pt() + promptLeptons.at(3)->pt() < ptThres_ )
      return;

    if ( !allowNonPrompt_) {
      if ( promptLeptons.at(2)->pt() < ptThres_ || promptLeptons.at(3)->pt() < ptThres2nd_ )
        return;

      if ( heeps2.size() > 1 && promptLeptons.at(3)->pt() < ptThres_ )
        return;
    }
  }

  auto sortByEt = [] (const pat::ElectronRef& a, const pat::ElectronRef& b) { return a->et() > b->et(); };

  if (!heeps1.empty()) {
    std::sort(heeps1.begin(),heeps1.end(),sortByEt);
    fillByGsfTrack(iEvent,iSetup,heeps1,std::vector<reco::GenParticleRef>({promptLeptons.at(0),promptLeptons.at(1)}));
  }

  if (!heeps2.empty()) {
    std::sort(heeps2.begin(),heeps2.end(),sortByEt);
    fillByGsfTrack(iEvent,iSetup,heeps2,std::vector<reco::GenParticleRef>({promptLeptons.at(2),promptLeptons.at(3)}));
  }
}

void MergedEleSigMvaInput::fillByGsfTrack(const edm::Event& iEvent,
                                          const edm::EventSetup& iSetup,
                                          const std::vector<pat::ElectronRef>& recoPair,
                                          const std::vector<reco::GenParticleRef>& promptPair) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> ecalIsoMapHandle;
  iEvent.getByToken(ecalIsoToken_, ecalIsoMapHandle);

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

  edm::Handle<edm::ValueMap<float>> union5x5covIpIpHandle;
  iEvent.getByToken(union5x5covIpIpToken_, union5x5covIpIpHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dPhiInHandle;
  iEvent.getByToken(union5x5dPhiInToken_, union5x5dPhiInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<EcalRecHitCollection> EBrecHitHandle;
  iEvent.getByToken(EBrecHitToken_, EBrecHitHandle);

  edm::Handle<EcalRecHitCollection> EErecHitHandle;
  iEvent.getByToken(EErecHitToken_, EErecHitHandle);

  if (recoPair.size() > 2)
    return;

  for (const auto aEle : recoPair) {
    auto addGsfTrk = (*addGsfTrkHandle)[aEle];
    auto addPackedCand = (*addPackedCandHandle)[aEle];
    auto orgGsfTrk = aEle->gsfTrack();

    const auto dEtaVariables = ModifiedDEtaInSeed::variables((*dPerpInHandle)[aEle],
                                                             (*dEtaInSeed2ndHandle)[aEle],
                                                             (*dPhiInSC2ndHandle)[aEle],
                                                             (*alphaTrackHandle)[aEle],
                                                             (*normDParaInHandle)[aEle]);
    const auto ssVariables = ModifiedShowerShape::variables((*union5x5covIeIeHandle)[aEle],
                                                             (*union5x5covIeIpHandle)[aEle],
                                                             (*union5x5covIpIpHandle)[aEle],
                                                             (*alphaCaloHandle)[aEle],
                                                             (*union5x5dEtaInHandle)[aEle],
                                                             (*union5x5dPhiInHandle)[aEle],
                                                             (*union5x5EnergyHandle)[aEle]);

    const EcalRecHitCollection* ecalRecHits = aEle->isEB() ? &(*EBrecHitHandle) : &(*EErecHitHandle);
    const float genPt = promptPair.at(0)->pt() + promptPair.at(1)->pt();
    const float genE = promptPair.at(0)->p() + promptPair.at(1)->p();

    if (allowNonPrompt_) {
      int32_t bitmap = aEle->userInt("modifiedHeepElectronID");
      int32_t mask = 0x00000381; // = 0011 1000 0001 - 1st for min Et, 7th for trk iso, 8th for EM+HadD1 iso, 9th for dxy
      int32_t pass = bitmap | mask;
      bool passMaskedId = pass==0x00000FFF; // HEEP ID has 12 cuts

      if ( !passMaskedId )
        continue;
    } else {
      if ( !aEle->electronID("modifiedHeepElectronID") )
        continue;
    }

    if ( addGsfTrk==orgGsfTrk && addPackedCand.isNull() ) { // ME w/o add GSF
      if (recoPair.size()==1)
        aHelper_.fillElectrons(aEle,
                               (*trkIsoMapHandle)[aEle],
                               (*ecalIsoMapHandle)[aEle],
                               dEtaVariables,
                               ssVariables,
                               ecalRecHits,
                               iSetup,
                               "mergedEl2",
			       consumesCollector(),
                               genPt,
                               genE);
    } else { // ME w/ GSF
      const bool isPackedCand = addGsfTrk==orgGsfTrk && addPackedCand.isNonnull();
      const reco::TrackBase* addTrk = isPackedCand ? addPackedCand->bestTrack() : addGsfTrk.get();

      bool matched1st = false;
      bool matched2nd = false;

      for (const auto& genptc : promptPair) {
        if ( reco::deltaR2(genptc->eta(),genptc->phi(),orgGsfTrk->eta(),orgGsfTrk->phi()) < drThres2nd2_ )
          matched1st = true;
        if ( reco::deltaR2(genptc->eta(),genptc->phi(),addTrk->eta(),addTrk->phi()) < drThres2nd2_ )
          matched2nd = true;
      }

      if ( !matched1st || !matched2nd )
        continue;

      aHelper_.fillElectrons(aEle,
                             (*trkIsoMapHandle)[aEle],
                             (*ecalIsoMapHandle)[aEle],
                             dEtaVariables,
                             ssVariables,
                             ecalRecHits,
                             iSetup,
                             "mergedEl1",
			     consumesCollector(),
                             genPt,
                             genE);
      aHelper_.fillAddTracks(aEle,
                             addTrk,
                             dEtaVariables,
                             ecalRecHits,
                             iSetup,
                             "mergedEl1Trk");
    }
  }

  return;
}

DEFINE_FWK_MODULE(MergedEleSigMvaInput);
