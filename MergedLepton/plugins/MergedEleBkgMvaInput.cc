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

// produce TTree for merged electron training with bkg MC (QCD, WJets, DY, TT)

class MergedEleBkgMvaInput : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleBkgMvaInput(const edm::ParameterSet&);
  virtual ~MergedEleBkgMvaInput() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> nrSatCrysMapToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;

  const double ptThres_;
  const double drThres_;

  MergedLeptonHelper aHelper_;

  std::map<std::string,TH1*> histo1d_;
};

MergedEleBkgMvaInput::MergedEleBkgMvaInput(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
nrSatCrysMapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("nrSatCrysMap"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")),
aHelper_() {
  usesResource("TFileService");
}

void MergedEleBkgMvaInput::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  aHelper_.SetFileService(&fs);
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  aHelper_.initElectronTree("ElectronStruct","fake","el");
  aHelper_.initElectronTree("ElectronStruct","bkg","el");
  aHelper_.initAddGsfTree("AddGsfStruct","fakeGsf","addGsf");
}

void MergedEleBkgMvaInput::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> ecalIsoMapHandle;
  iEvent.getByToken(ecalIsoToken_, ecalIsoMapHandle);

  edm::Handle<edm::ValueMap<int>> nrSatCrysHandle;
  iEvent.getByToken(nrSatCrysMapToken_, nrSatCrysHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_,rhoHandle);

  edm::Handle<reco::ConversionCollection> conversionsHandle;
  iEvent.getByToken(conversionsToken_, conversionsHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

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

  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->ptrAt(idx);
    MergedLeptonIDs::cutflowElectron cutflow = MergedLeptonIDs::cutflowElectron::baseline;

    if ( aEle->et() < ptThres_ )
      continue;

    bool passModifiedHEEP =
      MergedLeptonIDs::isModifiedHEEP(*aEle,
                                      *primaryVertex,
                                      (*trkIsoMapHandle)[aEle],
                                      (*ecalIsoMapHandle)[aEle],
                                      (*nrSatCrysHandle)[aEle],
                                      *rhoHandle,
                                      cutflow);

    MergedLeptonIDs::cutflowElectron cutflow_HEEP = MergedLeptonIDs::cutflowElectron::baseline;
    bool passHEEP =
      MergedLeptonIDs::hasPassedHEEP(*aEle,
                                     *primaryVertex,
                                     (*nrSatCrysHandle)[aEle],
                                     *rhoHandle,
                                     cutflow_HEEP);

    if ( !passModifiedHEEP )
      continue;

    // requirement for add GSF track
    const auto& addGsfTrk = (*addGsfTrkHandle)[aEle];
    const auto& orgGsfTrk = aEle->gsfTrack();

    std::string treename = "";

    if ( MergedLeptonIDs::isSameGsfTrack(addGsfTrk,orgGsfTrk) ) {
      treename = "bkg";

      if ( !passHEEP )
        continue;

    } else {
      // find whether add GSF track has a corresponding electron
      treename = "fake";
      bool notMerged = false;

      for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
        const auto& secEle = eleHandle->refAt(jdx);

        double dr2 = reco::deltaR2(aEle->eta(),aEle->phi(),secEle->eta(),secEle->phi());

        if (dr2 > drThres2)
          continue;

        const auto& secGsfTrk = secEle->gsfTrack();

        if ( MergedLeptonIDs::isSameGsfTrack(addGsfTrk,secGsfTrk) ) {
          notMerged = true;
          break;
        }
      }

      if ( notMerged )
        continue;

      aHelper_.fillGsfTracks((*addGsfTrkHandle)[aEle],
                             aEle->gsfTrack(),
                             "fakeGsf");
    }

    aHelper_.fillElectrons(aEle,
                           (*trkIsoMapHandle)[aEle],
                           (*ecalIsoMapHandle)[aEle],
                           (*nrSatCrysHandle)[aEle],
                           conversionsHandle,
                           beamSpotHandle,
                           treename);
  }
}

DEFINE_FWK_MODULE(MergedEleBkgMvaInput);
