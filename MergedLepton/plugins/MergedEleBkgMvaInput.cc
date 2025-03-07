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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

// produce TTree for merged electron training with bkg MC (QCD, WJets, DY, TT)

class MergedEleBkgMvaInput : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleBkgMvaInput(const edm::ParameterSet&);
  virtual ~MergedEleBkgMvaInput() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<reco::GsfTrack>> GsfTrkToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> candToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTrackToken_;
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
  const edm::EDGetTokenT<EcalRecHitCollection> EBrecHitToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EErecHitToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  // Ecal position calculation algorithm
  PositionCalc posCalcLog_;

  const double ptThres_;
  const double ptThres2nd_;
  const double drThres_;
  const double drThres2_;
  const bool select0J_;
  const bool selectHT_;
  const double maxHT_;

  MergedLeptonHelper aHelper_;

  std::map<std::string,TH1*> histo1d_;
};

MergedEleBkgMvaInput::MergedEleBkgMvaInput(const edm::ParameterSet& iConfig)
: srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
  srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
  pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
  GsfTrkToken_(consumes<edm::View<reco::GsfTrack>>(iConfig.getParameter<edm::InputTag>("srcGsfTrack"))),
  candToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("srcPackedCand"))),
  lostTrackToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("srcLostTracks"))),
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
  EBrecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBrecHits"))),
  EErecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EErecHits"))),
  generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
  lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
  posCalcLog_(PositionCalc(iConfig.getParameter<edm::ParameterSet>("posCalcLog"))),
  ptThres_(iConfig.getParameter<double>("ptThres")),
  ptThres2nd_(iConfig.getParameter<double>("ptThres2nd")),
  drThres_(iConfig.getParameter<double>("drThres")),
  drThres2_(drThres_*drThres_),
  select0J_(iConfig.getParameter<bool>("select0J")),
  selectHT_(iConfig.getParameter<bool>("selectHT")),
  maxHT_(iConfig.getParameter<double>("maxHT")),
  aHelper_() {
  usesResource("TFileService");
}

void MergedEleBkgMvaInput::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  aHelper_.SetFileService(&fs);
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["lheNj"] = fs->make<TH1D>("lheNj","lheNj",5,0.,5.);
  histo1d_["lheNj_cut"] = fs->make<TH1D>("lheNj_cut","lheNj",5,0.,5.);
  histo1d_["lheHT"] = fs->make<TH1D>("lheHT","lheHT",400,0.,4000.);
  histo1d_["lheHT_cut"] = fs->make<TH1D>("lheHT_cut","lheHT",400,0.,4000.);

  aHelper_.initElectronTree("ElectronStruct","fake","el");
  aHelper_.initElectronTree("ElectronStruct","bkg","el");
  aHelper_.initAddTrkTree("AddTrkStruct","fakeTrk","addTrk");

  aHelper_.SetPositionCalcLog(posCalcLog_);
}

void MergedEleBkgMvaInput::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

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

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  edm::Handle<EcalRecHitCollection> EBrecHitHandle;
  iEvent.getByToken(EBrecHitToken_, EBrecHitHandle);

  edm::Handle<EcalRecHitCollection> EErecHitHandle;
  iEvent.getByToken(EErecHitToken_, EErecHitHandle);

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken_, genInfo);
  double mcweight = genInfo->weight();

  double aWeight = mcweight/std::abs(mcweight);
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);
  aHelper_.SetMCweight(mcweight);

  if (select0J_ || selectHT_) {
    edm::Handle<LHEEventProduct> lheEventHandle;
    iEvent.getByToken(lheToken_, lheEventHandle);

    const auto& hepeup = lheEventHandle->hepeup();
    const auto& pup = hepeup.PUP;
    unsigned int lheNj = 0;
    double lheHT = 0.;

    for (unsigned int i = 0, n = pup.size(); i < n; ++i) {
      int status = hepeup.ISTUP[i];
      int idabs = std::abs(hepeup.IDUP[i]);

      if ((status == 1) && ((idabs == 21) || (idabs > 0 && idabs < 7))) { //# gluons and quarks
        lheNj++;

        double pt = std::hypot(pup[i][0], pup[i][1]);  // first entry is px, second py
        lheHT += pt;
      }
    }

    histo1d_["lheNj"]->Fill( static_cast<float>(lheNj)+0.5 );
    histo1d_["lheHT"]->Fill( lheHT );

    if (select0J_) {
      if (lheNj > 0)
        return;
    }

    histo1d_["lheNj_cut"]->Fill( static_cast<float>(lheNj)+0.5 );

    if (selectHT_) {
      if (lheHT > maxHT_)
        return;
    }

    histo1d_["lheHT_cut"]->Fill( lheHT );
  }

  std::vector<reco::GenParticleRef> promptLeptons;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    const auto& genptc = genptcHandle->refAt(idx);

    if ( ( std::abs(genptc->pdgId())==11 ) &&
         genptc->fromHardProcessFinalState() &&
         ( std::abs(genptc->eta()) < 2.5 ) )
      promptLeptons.push_back(genptc.castTo<reco::GenParticleRef>());
  }

  reco::VertexRef primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = pvHandle->refAt(0).castTo<reco::VertexRef>();
    aHelper_.SetPV(primaryVertex);
  } else
    return;

  aHelper_.SetBS(beamSpotHandle.product());

  std::vector<pat::ElectronRef> heepElectrons;

  for (unsigned idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( !aEle->electronID("modifiedHeepElectronID") )
      continue;

    auto castEle = aEle.castTo<pat::ElectronRef>();

    heepElectrons.push_back(castEle);
  }

  if ( heepElectrons.size() < 2 )
    return;

  for (unsigned int idx = 0; idx < heepElectrons.size(); ++idx) {
    const auto& aEle = heepElectrons.at(idx);

    if ( aEle->et() < ptThres_ )
      continue;

    float genPt = -1.;
    float genE = -1.;

    for (unsigned igen = 0; igen < promptLeptons.size(); ++igen) {
      const auto& genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aEle->gsfTrack()->eta(),aEle->gsfTrack()->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2_) {
        genPt = genptc->pt();
        genE = genptc->p();

        break;
      }
    }

    edm::Handle<edm::View<reco::GsfTrack>> GsfTrkHandle;
    iEvent.getByToken(GsfTrkToken_, GsfTrkHandle);

    edm::Handle<edm::View<pat::PackedCandidate>> candHandle;
    iEvent.getByToken(candToken_, candHandle);

    edm::Handle<edm::View<pat::PackedCandidate>> lostTrackHandle;
    iEvent.getByToken(lostTrackToken_, lostTrackHandle);

    int nGSFtrk = 0;
    int nKFtrk = 0;

    for (unsigned idx = 0; idx < GsfTrkHandle->size(); idx++) {
      const auto& aTrk = GsfTrkHandle->refAt(idx);

      if (aTrk->pt() < ptThres2nd_)
        continue;

      if (aTrk.castTo<reco::GsfTrackRef>()==aEle->gsfTrack())
        continue;

      if (reco::deltaR2(aEle->eta(),aEle->phi(),aTrk->eta(),aTrk->phi()) < drThres2_)
        nGSFtrk++;
    }

    for (auto& handle : {candHandle,lostTrackHandle}) {
      for (unsigned idx = 0; idx < handle->size(); idx++) {
        const auto& aTrk = handle->refAt(idx);

        if (aTrk->pt() < ptThres2nd_)
          continue;

        if (std::abs(aTrk->pdgId())==11)
          continue;

        if ( aEle->closestCtfTrackRef().isNonnull() &&
             reco::deltaR2( aTrk->eta(), aTrk->phi(),
                            aEle->closestCtfTrackRef()->eta(), aEle->closestCtfTrackRef()->phi() ) < 0.001*0.001 )
          continue;

        if (reco::deltaR2(aEle->eta(),aEle->phi(),aTrk->eta(),aTrk->phi()) < drThres2_)
          nKFtrk++;
      }
    }

    const EcalRecHitCollection* ecalRecHits = aEle->isEB() ? &(*EBrecHitHandle) : &(*EErecHitHandle);

    // requirement for add GSF track
    const auto& addGsfTrk = (*addGsfTrkHandle)[aEle];
    const auto& addPackedCand = (*addPackedCandHandle)[aEle];
    const auto& orgGsfTrk = aEle->gsfTrack();

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

    std::string treename = "";

    if ( addGsfTrk==orgGsfTrk && addPackedCand.isNull() ) {
      treename = "bkg";
    } else {
      // find whether add GSF track has a corresponding electron
      treename = "fake";

      if ( MergedLeptonHelperFct::isNotMerged(aEle,eleHandle,addGsfTrk) )
        continue;

      const bool isPackedCand = addGsfTrk==orgGsfTrk && addPackedCand.isNonnull();
      const reco::Track* addTrk = isPackedCand ? addPackedCand->bestTrack() : addGsfTrk.get();

      aHelper_.fillAddTracks(aEle,
                             addTrk,
                             dEtaVariables,
                             ecalRecHits,
                             iSetup,
                             "fakeTrk",
                             isPackedCand);
    }

    aHelper_.fillElectrons(aEle,
                           (*trkIsoMapHandle)[aEle],
                           (*ecalIsoMapHandle)[aEle],
                           dEtaVariables,
                           ssVariables,
                           ecalRecHits,
                           iSetup,
                           treename,
                           genPt,genE,
                           nGSFtrk,nKFtrk);
  }
}

DEFINE_FWK_MODULE(MergedEleBkgMvaInput);
