#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include <cmath>

// a hack of RecoEgamma/ElectronIdentification/plugins/cuts/GsfEleDEtaInSeedCut.cc

class GsfEleModifiedDEtaInSeedCutRun3 : public CutApplicatorWithEventContentBase {
public:
  GsfEleModifiedDEtaInSeedCutRun3(const edm::ParameterSet& c);

  result_type operator() (const reco::GsfElectronPtr&) const final;

  void setConsumes(edm::ConsumesCollector&) final;
  void getEventContent(const edm::EventBase&) final;

  double value(const reco::CandidatePtr& cand) const final;

  CandidateType candidateType() const final { return ELECTRON; }

private:
  const double _dEtaInSeedCutValueEB, _dEtaInSeedCutValueEE, _modifiedDEtaInSeedCutValueEB, _barrelCutOff;
  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle_;
  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle_;
  edm::Handle<edm::ValueMap<float>> modifiedDEtaInSeedHandle_;
};

DEFINE_EDM_PLUGIN(CutApplicatorFactory,
                  GsfEleModifiedDEtaInSeedCutRun3,
                  "GsfEleModifiedDEtaInSeedCutRun3");

GsfEleModifiedDEtaInSeedCutRun3::GsfEleModifiedDEtaInSeedCutRun3(const edm::ParameterSet& c) :
  CutApplicatorWithEventContentBase(c),
  _dEtaInSeedCutValueEB(c.getParameter<double>("dEtaInSeedCutValueEB")),
  _dEtaInSeedCutValueEE(c.getParameter<double>("dEtaInSeedCutValueEE")),
  _modifiedDEtaInSeedCutValueEB(c.getParameter<double>("modifiedDEtaInSeedCutValueEB")),
  _barrelCutOff(c.getParameter<double>("barrelCutOff"))
{
  contentTags_.emplace("addGsfTrk",c.getParameter<edm::InputTag>("addGsfTrk"));
  contentTags_.emplace("addPackedCand",c.getParameter<edm::InputTag>("addPackedCand"));
  contentTags_.emplace("modifiedDEtaInSeed",c.getParameter<edm::InputTag>("modifiedDEtaInSeed"));
}

void GsfEleModifiedDEtaInSeedCutRun3::setConsumes(edm::ConsumesCollector& cc) {
  contentTokens_.emplace("addGsfTrk",cc.consumes<edm::ValueMap<reco::GsfTrackRef>>(contentTags_["addGsfTrk"]));
  contentTokens_.emplace("addPackedCand",cc.consumes<edm::ValueMap<pat::PackedCandidateRef>>(contentTags_["addPackedCand"]));
  contentTokens_.emplace("modifiedDEtaInSeed",cc.consumes<edm::ValueMap<float>>(contentTags_["modifiedDEtaInSeed"]));
}

void GsfEleModifiedDEtaInSeedCutRun3::getEventContent(const edm::EventBase& ev) {
  ev.getByLabel(contentTags_["addGsfTrk"],addGsfTrkHandle_);
  ev.getByLabel(contentTags_["addPackedCand"],addPackedCandHandle_);
  ev.getByLabel(contentTags_["modifiedDEtaInSeed"],modifiedDEtaInSeedHandle_);
}

CutApplicatorBase::result_type GsfEleModifiedDEtaInSeedCutRun3::operator() (const reco::GsfElectronPtr& cand) const {
  const bool has2ndTrk = (*addGsfTrkHandle_)[cand]!=cand->gsfTrack() || (*addPackedCandHandle_)[cand].isNonnull();
  const bool isEB = std::abs(cand->superCluster()->eta()) < _barrelCutOff;

  // default cut w/o modification
  double dEtaInSeedCutValue = ( isEB ? _dEtaInSeedCutValueEB : _dEtaInSeedCutValueEE );

  // modified cut w/ the 2nd track
  if (has2ndTrk && isEB)
    dEtaInSeedCutValue = _modifiedDEtaInSeedCutValueEB;

  return std::abs(value(cand)) < dEtaInSeedCutValue;
}

double GsfEleModifiedDEtaInSeedCutRun3::value(const reco::CandidatePtr& cand) const {
  reco::GsfElectronPtr ele(cand);

  if ( ele->superCluster().isNull() || ele->superCluster()->seed().isNull() )
    return std::numeric_limits<double>::max(); // shouldn't happen

  double val = ele->deltaEtaSuperClusterTrackAtVtx() - ele->superCluster()->eta() + ele->superCluster()->seed()->eta();

  const bool has2ndTrk = (*addGsfTrkHandle_)[ele]!=ele->gsfTrack() || (*addPackedCandHandle_)[ele].isNonnull();
  const bool isEB = std::abs(ele->superCluster()->eta()) < _barrelCutOff;

  if (has2ndTrk && isEB) {
    float modifiedDEtaInSeed = (*modifiedDEtaInSeedHandle_)[ele];

    if ( modifiedDEtaInSeed >= std::numeric_limits<float>::max() )
      return modifiedDEtaInSeed; // extrapolation failed, should be rejected

    val = (std::abs(val) < std::abs(modifiedDEtaInSeed)) ? val : modifiedDEtaInSeed;
  }

  if (has2ndTrk && !isEB) {
    const float dPerpIn = (*modifiedDEtaInSeedHandle_)[ele];
    // Both EE variables use dEtaInSeedCutValueEE: the smaller magnitude
    // implements OR. Failed extrapolation must not veto a passing dEtaInSeed.
    if (std::isfinite(dPerpIn) && std::abs(dPerpIn) < std::abs(val))
      val = dPerpIn;
  }

  return val;
}
