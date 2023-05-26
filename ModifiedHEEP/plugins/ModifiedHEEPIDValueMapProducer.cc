#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Handle.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedEleTkIsolFromCands.h"
#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedDEtaInSeed.h"

#include <memory>
#include <vector>

//Heavily inspired from ElectronIDValueMapProducer

class ModifiedHEEPIDValueMapProducer : public edm::stream::EDProducer<> {
private:
  //helper classes to handle AOD vs MiniAOD
  template<typename T>
  struct DualToken {
    edm::EDGetTokenT<T> aod;
    edm::EDGetTokenT<T> miniAOD;
  };

  class DataFormat {
  public:
    enum Format{AUTO=0,AOD=1,MINIAOD=2};
  private:
    int data_;
  public:
    DataFormat(int val) : data_(val) {}
    bool tryAOD() const { return data_==AUTO || data_==AOD; }
    bool tryMiniAOD() const { return data_==AUTO || data_==MINIAOD; }
    int operator()() const { return data_; }
  };

public:
  explicit ModifiedHEEPIDValueMapProducer(const edm::ParameterSet&);
  ~ModifiedHEEPIDValueMapProducer() override;

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  template<typename T>
  static void writeValueMap(edm::Event &iEvent,
                            const edm::Handle<edm::View<pat::Electron>>& handle,
                            const std::vector<T> & values,
                            const std::string& label);

  static float calTrkIso(const pat::Electron& ele,
                         const std::vector<edm::Handle<edm::View<pat::PackedCandidate>>>& handles,
                         const reco::TrackBase& gsfTrk,
                         const std::vector<ModifiedEleTkIsolFromCands::PIDVeto>& pidVetos,
                         const ModifiedEleTkIsolFromCands& trkIsoCalc);

  template <typename T> void setToken(edm::EDGetTokenT<T>& token, edm::InputTag tag) { token=consumes<T>(tag); }
  template <typename T> void setToken(edm::EDGetTokenT<T>& token, const edm::ParameterSet& iPara, const std::string& tag) { token=consumes<T>(iPara.getParameter<edm::InputTag>(tag)); }
  template <typename T> void setToken(std::vector<edm::EDGetTokenT<T>>& tokens, const edm::ParameterSet& iPara, const std::string& tagName) {
    auto tags = iPara.getParameter<std::vector<edm::InputTag> >(tagName);
    for(auto& tag : tags) {
      edm::EDGetTokenT<T> token;
      setToken(token,tag);
      tokens.push_back(token);
    }
  }

  template <typename T> void setToken(DualToken<T>& token, const edm::ParameterSet& iPara, const std::string& tagAOD, const std::string& tagMiniAOD, DataFormat format) {
    if(format.tryAOD())
      token.aod=consumes<T>(iPara.getParameter<edm::InputTag>(tagAOD));
    if(format.tryMiniAOD())
      token.miniAOD=consumes<T>(iPara.getParameter<edm::InputTag>(tagMiniAOD));
  }

  template <typename T> void setToken(std::vector<DualToken<T>>& tokens, const edm::ParameterSet& iPara, const std::string& tagAOD, const std::string& tagMiniAOD,DataFormat format) {
    auto tagsAOD = iPara.getParameter<std::vector<edm::InputTag>>(tagAOD);
    auto tagsMiniAOD = iPara.getParameter<std::vector<edm::InputTag>>(tagMiniAOD);
    size_t maxSize = std::max(tagsAOD.size(),tagsMiniAOD.size());
    tokens.clear();
    tokens.resize(maxSize);

    if (format.tryAOD()) {
      for(size_t tagNr=0; tagNr < tagsAOD.size(); tagNr++)
      	setToken(tokens[tagNr].aod,tagsAOD[tagNr]);
    }

    if (format.tryMiniAOD()) {
      for(size_t tagNr=0; tagNr < tagsMiniAOD.size(); tagNr++)
      	setToken(tokens[tagNr].miniAOD,tagsMiniAOD[tagNr]);
    }
  }

  template<typename T>
  static edm::Handle<T> getHandle(const edm::Event& iEvent,
                                  const edm::EDGetTokenT<T>& token) {
    edm::Handle<T> handle;
    iEvent.getByToken(token,handle);

    return handle;
  }

  template<typename T>
  static edm::Handle<T> getHandle(const edm::Event& iEvent,
                                  const DualToken<T>& token) {
    edm::Handle<T> handle;
    if (!token.aod.isUninitialized())
      iEvent.getByToken(token.aod,handle);
    if (!handle.isValid() && !token.miniAOD.isUninitialized())
      iEvent.getByToken(token.miniAOD,handle);

    return handle;
  }

  template<typename T>
  static std::vector<edm::Handle<T>> getHandles(const edm::Event& iEvent, const std::vector<DualToken<T>>& tokens) {
    std::vector<edm::Handle<T>> handles(tokens.size());
    if (tokens.empty())
      return handles;

    if (!tokens[0].aod.isUninitialized())
      iEvent.getByToken(tokens[0].aod,handles[0]);

    bool isAOD = handles[0].isValid();

    if ( !isAOD && !tokens[0].miniAOD.isUninitialized() )
      iEvent.getByToken(tokens[0].miniAOD,handles[0]);

    for (size_t tokenNr=1; tokenNr < tokens.size(); tokenNr++) {
      auto token = isAOD ? tokens[tokenNr].aod : tokens[tokenNr].miniAOD;

      if (!token.isUninitialized())
        iEvent.getByToken(token,handles[tokenNr]);
    }

    return handles;
  }

  template<typename T>
  static bool isEventAOD(const edm::Event& iEvent, const DualToken<T>& token) {
    edm::Handle<T> handle;

    if (!token.aod.isUninitialized())
      iEvent.getByToken(token.aod,handle);

    if (handle.isValid())
      return true;
    else
      return false;
  }

  DualToken<edm::View<pat::Electron>> eleToken_;
  std::vector<DualToken<edm::View<pat::PackedCandidate>>> candTokens_;
  DualToken<edm::View<reco::GsfTrack>> gsfTrkToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  DualToken<EcalRecHitCollection> EBrecHitToken_;
  DualToken<EcalRecHitCollection> EErecHitToken_;

  ModifiedEleTkIsolFromCands trkIsoCalc_;
  ModifiedEleTkIsolFromCands trkIso04Calc_;
  ModifiedDEtaInSeed dEtaInSeedCalc_;
  bool makeTrkIso04_;
  DataFormat dataFormat_;
  std::vector<ModifiedEleTkIsolFromCands::PIDVeto> candVetosAOD_;
  std::vector<ModifiedEleTkIsolFromCands::PIDVeto> candVetosMiniAOD_;
  const std::vector<edm::InputTag> candsTag_;

  const std::string eleTrkPtIsoLabel_ = "eleTrkPtIso";
  const std::string eleTrkPtIso04Label_ = "eleTrkPtIso04";
  const std::string eleAddGsfTrkLabel_ = "eleAddGsfTrk";
  const std::string eleAddPackedCandLabel_ = "eleAddPackedCand";
  const std::string modifiedDEtaInSeedLabel_ = "modifiedDEtaInSeed";
  const std::string alphaTrackLabel_ = "alphaTrack";
  const std::string alphaCaloLabel_ = "alphaCalo";
  const std::string union5x5ratioDRlabel_ = "union5x5ratioDR";
};

ModifiedHEEPIDValueMapProducer::ModifiedHEEPIDValueMapProducer(const edm::ParameterSet& iConfig) :
trkIsoCalc_(iConfig.getParameter<edm::ParameterSet>("trkIsoConfig")),
trkIso04Calc_(iConfig.getParameter<edm::ParameterSet>("trkIso04Config")),
dEtaInSeedCalc_(PositionCalc(iConfig.getParameter<edm::ParameterSet>("posCalcLog"))),
makeTrkIso04_(iConfig.getParameter<bool>("makeTrkIso04")),
dataFormat_(iConfig.getParameter<int>("dataFormat")),
candsTag_(iConfig.getParameter<std::vector<edm::InputTag>>("candsMiniAOD")) {
  setToken(eleToken_,iConfig,"elesAOD","elesMiniAOD",dataFormat_);
  setToken(candTokens_,iConfig,"candsAOD","candsMiniAOD",dataFormat_);
  setToken(gsfTrkToken_,iConfig,"gsfTrksAOD","gsfTrksMiniAOD",dataFormat_);
  setToken(beamSpotToken_,iConfig,"beamSpot");
  setToken(EBrecHitToken_,iConfig,"EBrecHitsAOD","EBrecHitsMiniAOD",dataFormat_);
  setToken(EErecHitToken_,iConfig,"EErecHitsAOD","EErecHitsMiniAOD",dataFormat_);

  auto fillVetos = [] (const auto& in, auto& out) {
    std::transform(in.begin(),in.end(),std::back_inserter(out),ModifiedEleTkIsolFromCands::pidVetoFromStr);
  };

  fillVetos(iConfig.getParameter<std::vector<std::string>>("candVetosAOD"),candVetosAOD_);
  if ( candVetosAOD_.size()!=iConfig.getParameter<std::vector<edm::InputTag>>("candsAOD").size() ) {
    throw cms::Exception("ConfigError") << " Error candVetosAOD should be the same size as candsAOD " << std::endl;
  }

  fillVetos(iConfig.getParameter<std::vector<std::string>>("candVetosMiniAOD"),candVetosMiniAOD_);
  if ( candVetosMiniAOD_.size()!=iConfig.getParameter<std::vector<edm::InputTag>>("candsMiniAOD").size() ) {
    throw cms::Exception("ConfigError") << " Error candVetosMiniAOD should be the same size as candsMiniAOD " << std::endl;
  }

  produces<edm::ValueMap<float> >(eleTrkPtIsoLabel_);
  produces<edm::ValueMap<float> >(modifiedDEtaInSeedLabel_);
  produces<edm::ValueMap<float> >(alphaTrackLabel_);
  produces<edm::ValueMap<float> >(alphaCaloLabel_);
  produces<edm::ValueMap<float> >(union5x5ratioDRlabel_);

  if (makeTrkIso04_)
    produces<edm::ValueMap<float> >(eleTrkPtIso04Label_);

  produces<edm::ValueMap<reco::GsfTrackRef>>(eleAddGsfTrkLabel_);
  produces<edm::ValueMap<pat::PackedCandidateRef>>(eleAddPackedCandLabel_);
}

ModifiedHEEPIDValueMapProducer::~ModifiedHEEPIDValueMapProducer() {}

void ModifiedHEEPIDValueMapProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  auto eleHandle = getHandle(iEvent,eleToken_);
  auto beamSpotHandle = getHandle(iEvent,beamSpotToken_);
  auto candHandles = getHandles(iEvent,candTokens_);
  auto gsfTrkHandle = getHandle(iEvent,gsfTrkToken_);
  auto EBrecHitHandle = getHandle(iEvent,EBrecHitToken_);
  auto EErecHitHandle = getHandle(iEvent,EErecHitToken_);

  bool isAOD = isEventAOD(iEvent,eleToken_);
  const auto& candVetos = isAOD ? candVetosAOD_ : candVetosMiniAOD_;

  std::vector<float> eleTrkPtIso;
  std::vector<float> eleTrkPtIso04;
  std::vector<float> modifiedDEtaInSeed;
  std::vector<float> alphaTrack;
  std::vector<float> alphaCalo;
  std::vector<float> union5x5ratioDR;
  std::vector<reco::GsfTrackRef> eleAddGsfTrk;
  std::vector<pat::PackedCandidateRef> eleAddPackedCand;

  for (const auto& ele : *eleHandle) {
    const auto& additionalGsfTrk = trkIsoCalc_.additionalGsfTrkSelector(ele,gsfTrkHandle);
    auto additionalCand = trkIsoCalc_.additionalPackedCandSelector(ele,candHandles,candVetos);
    auto addTrk = reco::TrackBase(*additionalGsfTrk);

    if ( additionalCand.isNonnull() && ele.gsfTrack()==additionalGsfTrk )
      addTrk = *(additionalCand->bestTrack());

    auto modifiedDEtaInSeedVal = ModifiedDEtaInSeed::variables();

    if ( ele.gsfTrack()!=additionalGsfTrk || additionalCand.isNonnull() )
      modifiedDEtaInSeedVal = dEtaInSeedCalc_.value(ele,
                                                    ele.isEB() ? &(*EBrecHitHandle) : &(*EErecHitHandle),
                                                    addTrk,
                                                    *beamSpotHandle,
                                                    iSetup);

    eleTrkPtIso.push_back(calTrkIso(ele,candHandles,addTrk,candVetos,trkIsoCalc_));
    modifiedDEtaInSeed.push_back(modifiedDEtaInSeedVal.modifiedDEtaInSeed);
    alphaTrack.push_back(modifiedDEtaInSeedVal.alphaTrack);
    alphaCalo.push_back(modifiedDEtaInSeedVal.alphaCalo);
    union5x5ratioDR.push_back(modifiedDEtaInSeedVal.union5x5ratioDR);

    if (makeTrkIso04_)
      eleTrkPtIso04.push_back(calTrkIso(ele,candHandles,addTrk,candVetos,trkIso04Calc_));

    eleAddGsfTrk.push_back(additionalGsfTrk);
    eleAddPackedCand.push_back(additionalCand);
  }

  writeValueMap(iEvent,eleHandle,eleTrkPtIso,eleTrkPtIsoLabel_);
  writeValueMap(iEvent,eleHandle,modifiedDEtaInSeed,modifiedDEtaInSeedLabel_);
  writeValueMap(iEvent,eleHandle,alphaTrack,alphaTrackLabel_);
  writeValueMap(iEvent,eleHandle,alphaCalo,alphaCaloLabel_);
  writeValueMap(iEvent,eleHandle,union5x5ratioDR,union5x5ratioDRlabel_);

  if (makeTrkIso04_)
    writeValueMap(iEvent,eleHandle,eleTrkPtIso04,eleTrkPtIso04Label_);

  writeValueMap(iEvent,eleHandle,eleAddGsfTrk,eleAddGsfTrkLabel_);
  writeValueMap(iEvent,eleHandle,eleAddPackedCand,eleAddPackedCandLabel_);
}

float ModifiedHEEPIDValueMapProducer::calTrkIso(const pat::Electron& ele,
                                                const std::vector<edm::Handle<edm::View<pat::PackedCandidate>>>& handles,
                                                const reco::TrackBase& addTrk,
                                                const std::vector<ModifiedEleTkIsolFromCands::PIDVeto>& pidVetos,
                                                const ModifiedEleTkIsolFromCands& trkIsoCalc) {
  if (ele.gsfTrack().isNull())
    return std::numeric_limits<float>::max();
  else {
    float trkIso=0.;

    for (size_t handleNr=0; handleNr < handles.size(); handleNr++) {
      const auto& handle = handles[handleNr];

      if (handle.isValid()) {
        if (handleNr < pidVetos.size())
          trkIso += trkIsoCalc.calIsolPt(*ele.gsfTrack(),handle,addTrk,pidVetos[handleNr]);
        else
          throw cms::Exception("LogicError") << " somehow the pidVetos and handles do not much, given this is checked at construction time, something has gone wrong in the code handle nr " << handleNr << " size of vetos " << pidVetos.size();
      }
    }

    return trkIso;
  }
}

template<typename T>
void ModifiedHEEPIDValueMapProducer::writeValueMap(edm::Event& iEvent,
                                                   const edm::Handle<edm::View<pat::Electron>>& handle,
                                                   const std::vector<T>& values,
                                                   const std::string& label) {
  std::unique_ptr<edm::ValueMap<T>> valMap(new edm::ValueMap<T>());
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valMap),label);
}

DEFINE_FWK_MODULE(ModifiedHEEPIDValueMapProducer);
