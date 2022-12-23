#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"

// a hack of RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDEtaInSeedCut.cc

class GsfEleModifiedDEtaInSeedCut : public CutApplicatorWithEventContentBase {
public:
  GsfEleModifiedDEtaInSeedCut(const edm::ParameterSet& c);

  result_type operator() (const reco::GsfElectronPtr&) const final;

  void setConsumes(edm::ConsumesCollector&) final;
  void getEventContent(const edm::EventBase&) final;

  double value(const reco::CandidatePtr& cand) const final;

  CandidateType candidateType() const final { return ELECTRON; }

private:
  const double _dEtaInSeedCutValueEB, _dEtaInSeedCutValueEE, _barrelCutOff;
  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addTrkHandle_;
};

DEFINE_EDM_PLUGIN(CutApplicatorFactory,
                  GsfEleModifiedDEtaInSeedCut,
                  "GsfEleModifiedDEtaInSeedCut");

GsfEleModifiedDEtaInSeedCut::GsfEleModifiedDEtaInSeedCut(const edm::ParameterSet& c) :
  CutApplicatorWithEventContentBase(c),
  _dEtaInSeedCutValueEB(c.getParameter<double>("dEtaInSeedCutValueEB")),
  _dEtaInSeedCutValueEE(c.getParameter<double>("dEtaInSeedCutValueEE")),
  _barrelCutOff(c.getParameter<double>("barrelCutOff"))
{
  contentTags_.emplace("addTrkRef",c.getParameter<edm::InputTag>("addTrkRef"));
}

void GsfEleModifiedDEtaInSeedCut::setConsumes(edm::ConsumesCollector& cc) {
  contentTokens_.emplace("addTrkRef",cc.consumes<edm::ValueMap<reco::GsfTrackRef>>(contentTags_["addTrkRef"]));
}

void GsfEleModifiedDEtaInSeedCut::getEventContent(const edm::EventBase& ev) {
  ev.getByLabel(contentTags_["addTrkRef"],addTrkHandle_);
}

//a little temporary 72X fix
float dEtaInSeed(const reco::GsfElectronPtr& ele) {
  return ele->superCluster().isNonnull() && ele->superCluster()->seed().isNonnull() ?
    ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta() : std::numeric_limits<float>::max();
}

CutApplicatorBase::result_type GsfEleModifiedDEtaInSeedCut::operator() (const reco::GsfElectronPtr& cand) const {
  bool has2ndGsf = (*addTrkHandle_)[cand]!=cand->gsfTrack();

  const float dEtaInSeedCutValue =
    ( std::abs(cand->superCluster()->eta()) < _barrelCutOff ? _dEtaInSeedCutValueEB : _dEtaInSeedCutValueEE );

  return (has2ndGsf) ? true : std::abs(dEtaInSeed(cand)) < dEtaInSeedCutValue;
}

double GsfEleModifiedDEtaInSeedCut::value(const reco::CandidatePtr& cand) const {
  reco::GsfElectronPtr ele(cand);
  return std::abs(dEtaInSeed(ele));
}
