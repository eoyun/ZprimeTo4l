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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

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

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EBrecHitToken_;
  const edm::EDGetTokenT<EcalRecHitCollection> EErecHitToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  // Ecal position calculation algorithm
  PositionCalc posCalcLog_;
  PositionCalc posCalcLinear_;

  const double ptThres_;
  const bool select0J_;
  const bool selectHT_;
  const double maxHT_;

  MergedLeptonHelper aHelper_;

  std::map<std::string,TH1*> histo1d_;
};

MergedEleBkgMvaInput::MergedEleBkgMvaInput(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
EBrecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EBrecHits"))),
EErecHitToken_(consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("EErecHits"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
posCalcLog_(PositionCalc(iConfig.getParameter<edm::ParameterSet>("posCalcLog"))),
posCalcLinear_(PositionCalc(iConfig.getParameter<edm::ParameterSet>("posCalcLinear"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
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
  aHelper_.initAddGsfTree("AddGsfStruct","fakeGsf","addGsf");

  aHelper_.SetPositionCalcLog(posCalcLog_);
  aHelper_.SetPositionCalcLinear(posCalcLinear_);
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

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = pvHandle->ptrAt(0);
    aHelper_.SetPV(primaryVertex);
  } else
    return;

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

    const EcalRecHitCollection* ecalRecHits = aEle->isEB() ? &(*EBrecHitHandle) : &(*EErecHitHandle);

    // requirement for add GSF track
    const auto& addGsfTrk = (*addGsfTrkHandle)[aEle];
    const auto& orgGsfTrk = aEle->gsfTrack();

    std::string treename = "";

    if ( addGsfTrk==orgGsfTrk ) {
      treename = "bkg";
    } else {
      // find whether add GSF track has a corresponding electron
      treename = "fake";
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
        continue;

      aHelper_.fillGsfTracks(aEle,
                             (*addGsfTrkHandle)[aEle],
                             iSetup,
                             beamSpotHandle,
                             ecalRecHits,
                             "fakeGsf");
    }

    aHelper_.fillElectrons(aEle,
                           (*trkIsoMapHandle)[aEle],
                           (*ecalIsoMapHandle)[aEle],
                           addGsfTrk,
                           iSetup,
                           ecalRecHits,
                           treename);
  }
}

DEFINE_FWK_MODULE(MergedEleBkgMvaInput);
