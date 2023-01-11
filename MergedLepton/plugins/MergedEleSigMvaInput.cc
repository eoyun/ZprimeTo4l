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
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"
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

  void fillByGsfTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<edm::Ptr<pat::Electron>>& eles);

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EBrecHitToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EErecHitToken_;

  const double ptThres_;
  const double drThres_;

  std::map<std::string,TH1*> histo1d_;

  MergedLeptonHelper aHelper_;
};

MergedEleSigMvaInput::MergedEleSigMvaInput(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
EBrecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBrecHits"))),
EErecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EErecHits"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")),
aHelper_() {
  usesResource("TFileService");
}

void MergedEleSigMvaInput::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  aHelper_.SetFileService(&fs);
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  aHelper_.initElectronTree("ElectronStruct","heep1","el");
  aHelper_.initElectronTree("ElectronStruct","heep2","el");
  aHelper_.initElectronTree("ElectronStruct","mergedEl1","el");
  aHelper_.initElectronTree("ElectronStruct","mergedEl2","el");

  aHelper_.initAddGsfTree("AddGsfStruct","heep1Gsf","addGsf");
  aHelper_.initAddGsfTree("AddGsfStruct","heep2Gsf","addGsf");
  aHelper_.initAddGsfTree("AddGsfStruct","mergedEl1Gsf","addGsf");
}

void MergedEleSigMvaInput::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken_, genInfo);
  double mcweight = genInfo->weight();

  double aWeight = mcweight/std::abs(mcweight);
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  aHelper_.SetMCweight(mcweight);

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = pvHandle->ptrAt(0);
    aHelper_.SetPV(primaryVertex);
  } else
    return;

  std::vector<reco::GenParticleRef> promptLeptons;
  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    const auto& genptc = genptcHandle->refAt(idx);

    if ( ( std::abs(genptc->pdgId())==11 ) && genptc->isPromptFinalState() && genptc->pt() > ptThres_ )
      promptLeptons.push_back(genptc.castTo<reco::GenParticleRef>());
  }

  if (promptLeptons.size()!=4)
    return;

  // sort by mother ptc, and then by et (motherA Et1, motherA Et2, motherB Et1, motherB Et2)
  auto sisterLambda = [](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
    return (a->mother() == b->mother()) ? (a->et() > b->et()) : (a->mother() < b->mother());
  };

  std::sort(promptLeptons.begin(),promptLeptons.end(),sisterLambda);

  std::vector<edm::Ptr<pat::Electron>> heeps1, heeps2;

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->ptrAt(idx);

    if ( aEle->et() < ptThres_ )
      continue;

    bool matched = false;
    size_t igen = 0;

    // for given reco electron, check whether GEN-lv prompt electron exists within dR thres
    for (; igen < promptLeptons.size(); ++igen) {
      const auto& genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aEle->eta(),aEle->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if ( !matched )
      continue;

    if ( !aEle->electronID("modifiedHeepElectronID") )
      return;

    (igen < 2) ? heeps1.push_back(aEle) : heeps2.push_back(aEle);
  }

  // assume 4e only
  if ( heeps1.size() == 0 || heeps2.size() == 0 )
    return;

  fillByGsfTrack(iEvent,iSetup,heeps1);
  fillByGsfTrack(iEvent,iSetup,heeps2);
}

void MergedEleSigMvaInput::fillByGsfTrack(const edm::Event& iEvent, const edm::EventSetup& iSetup, std::vector<edm::Ptr<pat::Electron>>& eles) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> ecalIsoMapHandle;
  iEvent.getByToken(ecalIsoToken_, ecalIsoMapHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<reco::ConversionCollection> conversionsHandle;
  iEvent.getByToken(conversionsToken_, conversionsHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  edm::Handle<EcalRecHitCollection> EBrecHitHandle;
  iEvent.getByToken(EBrecHitToken_, EBrecHitHandle);

  edm::Handle<EcalRecHitCollection> EErecHitHandle;
  iEvent.getByToken(EErecHitToken_, EErecHitHandle);

  if ( eles.size()==1 ) {
    auto addGsfTrk = (*addGsfTrkHandle)[eles.front()];
    auto orgGsfTrk = eles.front()->gsfTrack();

    if ( addGsfTrk==orgGsfTrk ) { // ME w/o add GSF
      // veto any nearby electron over Et thres
      for (unsigned jdx = 0; jdx < eleHandle->size(); jdx++) {
        const auto& secEle = eleHandle->ptrAt(jdx);

        if (secEle==eles.front())
          continue;

        if (secEle->et() < ptThres_)
          continue;

        double dr2 = reco::deltaR2(eles.front()->eta(),eles.front()->phi(),secEle->eta(),secEle->phi());

        if (dr2 < drThres_*drThres_)
          return;
      }

      aHelper_.fillElectrons(eles.front(),
                             (*trkIsoMapHandle)[eles.front()],
                             (*ecalIsoMapHandle)[eles.front()],
                             addGsfTrk,
                             conversionsHandle,
                             beamSpotHandle,
                             iSetup,
                             EBrecHitHandle,
                             EErecHitHandle,
                             "mergedEl2");
    } else { // ME w/ GSF
      // find whether add GSF track has a corresponding electron
      bool notMerged = false;

      for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
        const auto& secEle = eleHandle->refAt(jdx);
        const auto& secGsfTrk = secEle->gsfTrack();

        if ( addGsfTrk==secGsfTrk ) {
          notMerged = true;
          break;
        }
      }

      if ( notMerged )
        return;

      aHelper_.fillElectrons(eles.front(),
                             (*trkIsoMapHandle)[eles.front()],
                             (*ecalIsoMapHandle)[eles.front()],
                             addGsfTrk,
                             conversionsHandle,
                             beamSpotHandle,
                             iSetup,
                             EBrecHitHandle,
                             EErecHitHandle,
                             "mergedEl1");
      aHelper_.fillGsfTracks((*addGsfTrkHandle)[eles.front()],eles.front()->gsfTrack(),"mergedEl1Gsf");
    }

    return;
  }

  // default case: eles.size() > 1
  std::sort(eles.begin(),eles.end(),[](const edm::Ptr<pat::Electron>& a, const edm::Ptr<pat::Electron>& b) { return a->et() > b->et(); });

  aHelper_.fillElectrons(eles.at(0),
                         (*trkIsoMapHandle)[eles.at(0)],
                         (*ecalIsoMapHandle)[eles.at(0)],
                         (*addGsfTrkHandle)[eles.at(0)],
                         conversionsHandle,
                         beamSpotHandle,
                         iSetup,
                         EBrecHitHandle,
                         EErecHitHandle,
                         "heep1");
  aHelper_.fillElectrons(eles.at(1),
                         (*trkIsoMapHandle)[eles.at(1)],
                         (*ecalIsoMapHandle)[eles.at(1)],
                         (*addGsfTrkHandle)[eles.at(1)],
                         conversionsHandle,
                         beamSpotHandle,
                         iSetup,
                         EBrecHitHandle,
                         EErecHitHandle,
                         "heep2");
  aHelper_.fillGsfTracks((*addGsfTrkHandle)[eles.at(0)],eles.at(0)->gsfTrack(),"heep1Gsf");
  aHelper_.fillGsfTracks((*addGsfTrkHandle)[eles.at(1)],eles.at(1)->gsfTrack(),"heep2Gsf");

  return;
}

DEFINE_FWK_MODULE(MergedEleSigMvaInput);
