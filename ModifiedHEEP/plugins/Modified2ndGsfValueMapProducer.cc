#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"

#include <memory>
#include <vector>

class Modified2ndGsfValueMapProducer : public edm::stream::EDProducer<> {
public:
  explicit Modified2ndGsfValueMapProducer(const edm::ParameterSet&);
  ~Modified2ndGsfValueMapProducer() override {}

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  template<typename T>
  void writeValueMap( edm::Event& iEvent, const edm::Handle<edm::View<pat::Electron>>& handle,
                      const std::vector<T>& values, const std::string& label);

  edm::EDGetTokenT<edm::View<pat::Electron>> eleToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  edm::EDGetTokenT<std::vector<reco::Vertex>> pvToken_;
  edm::ESHandle<TransientTrackBuilder> trackBuilder;
};

Modified2ndGsfValueMapProducer::Modified2ndGsfValueMapProducer(const edm::ParameterSet& iConfig)
: eleToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("eleSrc"))),
  addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
  pvToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("pvSrc"))) {

  produces<edm::ValueMap<int>>("addGsfCharge");
  produces<edm::ValueMap<int>>("addGsfLostHits");
  produces<edm::ValueMap<int>>("addGsfIsHighPurityTrack");
  produces<edm::ValueMap<float>>("addGsfDxy");
  produces<edm::ValueMap<float>>("addGsfDxyErr");
  produces<edm::ValueMap<float>>("addGsfDz");
  produces<edm::ValueMap<float>>("addGsfDzErr");
  produces<edm::ValueMap<float>>("addGsfPt");
  produces<edm::ValueMap<float>>("addGsfEta");
  produces<edm::ValueMap<float>>("addGsfPhi");
  produces<edm::ValueMap<float>>("addGsfIp3d");
  produces<edm::ValueMap<float>>("addGsfSip3d");
}

template<typename T>
void Modified2ndGsfValueMapProducer::writeValueMap( edm::Event& iEvent, const edm::Handle<edm::View<pat::Electron>>& handle,
                                                    const std::vector<T>& values, const std::string& label) {
  auto valMap = std::make_unique<edm::ValueMap<T>>();
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valMap),label);
}

void Modified2ndGsfValueMapProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(eleToken_, eleHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<std::vector<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);

  reco::Vertex primaryVertex;
  bool pvValid = false;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = pvHandle->at(0);
    pvValid = true;
  }

  unsigned int nEle = eleHandle->size();

  std::vector<int> charge, lostHits, isHighPurityTrack;
  std::vector<float> dxy, dxyErr, dz, dzErr, pt, eta, phi, ip3d, sip3d;
  charge.reserve(nEle);
  lostHits.reserve(nEle);
  isHighPurityTrack.reserve(nEle);
  dxy.reserve(nEle);
  dxyErr.reserve(nEle);
  dz.reserve(nEle);
  dzErr.reserve(nEle);
  pt.reserve(nEle);
  eta.reserve(nEle);
  phi.reserve(nEle);
  ip3d.reserve(nEle);
  sip3d.reserve(nEle);

  for (size_t iEle = 0; iEle < nEle; ++iEle) {
    auto trk = (*addGsfTrkMap)[eleHandle->ptrAt(iEle)];

    charge.emplace_back( trk->charge() );
    lostHits.emplace_back( trk->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) );
    isHighPurityTrack.emplace_back( trk->quality(reco::TrackBase::TrackQuality::highPurity) );
    dxy.emplace_back( trk->dxy(primaryVertex.position()) );
    dxyErr.emplace_back( trk->dxyError(primaryVertex.position(), primaryVertex.covariance()) );
    dz.emplace_back( trk->dz(primaryVertex.position()) );
    dzErr.emplace_back( std::abs(std::hypot(trk->dzError(), primaryVertex.zError())) );
    pt.emplace_back( trk->pt() );
    eta.emplace_back( trk->eta() );
    phi.emplace_back( trk->phi() );

    reco::TransientTrack tt = trackBuilder->build(trk);
    std::pair<bool, Measurement1D> result = IPTools::signedImpactParameter3D(tt, GlobalVector(trk->px(), trk->py(), trk->pz()), primaryVertex);
    double d0_corr = result.second.value();
    double d0_err = pvValid ? result.second.error() : -1.0;

    ip3d.emplace_back( std::abs(d0_corr) );
    sip3d.emplace_back( std::abs(d0_corr)/d0_err );
  }

  writeValueMap(iEvent,eleHandle,charge,"addGsfCharge");
  writeValueMap(iEvent,eleHandle,lostHits,"addGsfLostHits");
  writeValueMap(iEvent,eleHandle,isHighPurityTrack,"addGsfIsHighPurityTrack");
  writeValueMap(iEvent,eleHandle,dxy,"addGsfDxy");
  writeValueMap(iEvent,eleHandle,dxyErr,"addGsfDxyErr");
  writeValueMap(iEvent,eleHandle,dz,"addGsfDz");
  writeValueMap(iEvent,eleHandle,dzErr,"addGsfDzErr");
  writeValueMap(iEvent,eleHandle,pt,"addGsfPt");
  writeValueMap(iEvent,eleHandle,eta,"addGsfEta");
  writeValueMap(iEvent,eleHandle,phi,"addGsfPhi");
  writeValueMap(iEvent,eleHandle,ip3d,"addGsfIp3d");
  writeValueMap(iEvent,eleHandle,sip3d,"addGsfSip3d");
}

void Modified2ndGsfValueMapProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.add<edm::InputTag>("eleSrc",edm::InputTag("slimmedElectrons"));
  desc.add<edm::InputTag>("addGsfTrkMap",edm::InputTag("ModifiedHEEPIDVarValueMaps","eleAddGsfTrk"));
  desc.add<edm::InputTag>("pvSrc",edm::InputTag("offlineSlimmedPrimaryVertices"));

  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Modified2ndGsfValueMapProducer);
