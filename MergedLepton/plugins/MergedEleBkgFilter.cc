#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "CommonTools/Egamma/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

// produce TTree for merged electron training with bkg MC (QCD, WJets, DY, TT)

class MergedEleBkgFilter : public edm::one::EDFilter<edm::one::SharedResources> {
public:
  explicit MergedEleBkgFilter(const edm::ParameterSet&);
  virtual ~MergedEleBkgFilter() {}

private:
  bool filter(edm::Event&, const edm::EventSetup&) override;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;

  const double ptThres_;
  const bool select0J_;
  const bool selectHT_;
  const double maxHT_;
};

MergedEleBkgFilter::MergedEleBkgFilter(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
select0J_(iConfig.getParameter<bool>("select0J")),
selectHT_(iConfig.getParameter<bool>("selectHT")),
maxHT_(iConfig.getParameter<double>("maxHT"))
{}

bool MergedEleBkgFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
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

    if (select0J_) {
      if (lheNj > 0)
        return false;
    }

    if (selectHT_) {
      if (lheHT > maxHT_)
        return false;
    }
  }

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->ptrAt(0);
  else
    return false;

  std::vector<edm::Ptr<pat::Electron>> heepElectrons;

  for (unsigned idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->ptrAt(idx);

    if ( !aEle->electronID("modifiedHeepElectronID") )
      continue;

    heepElectrons.push_back(aEle);
  }

  if ( heepElectrons.size() < 2 )
    return false;

  for (unsigned int idx = 0; idx < heepElectrons.size(); ++idx) {
    const auto& aEle = heepElectrons.at(idx);

    if ( aEle->et() < ptThres_ )
      continue;
    // requirement for add GSF track
    const auto& addGsfTrk = (*addGsfTrkHandle)[aEle];
    const auto& orgGsfTrk = aEle->gsfTrack();

    if ( addGsfTrk==orgGsfTrk ) {
      continue;
    } else {
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
        continue;

      return true;
    }
  }

  return false;
}

DEFINE_FWK_MODULE(MergedEleBkgFilter);
