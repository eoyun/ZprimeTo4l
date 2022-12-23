#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/ElectronIdentification/interface/EBEECutValues.h"

class GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut : public CutApplicatorWithEventContentBase {
public:
  GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut(const edm::ParameterSet& c);

  result_type operator() (const reco::GsfElectronPtr&) const final;

  void setConsumes(edm::ConsumesCollector&) final;
  void getEventContent(const edm::EventBase&) final;

  double value(const reco::CandidatePtr& cand) const final;

  CandidateType candidateType() const final { return ELECTRON; }

private:
  EBEECutValues minE1x5OverE5x5Cut_;
  EBEECutValues minE2x5OverE5x5Cut_;
  EBEECutValuesInt maxNrSatCrysIn5x5Cut_;
  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addTrkHandle_;
};

DEFINE_EDM_PLUGIN(CutApplicatorFactory,
                  GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut,
                  "GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut");

GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut::GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut(const edm::ParameterSet& params) :
  CutApplicatorWithEventContentBase(params),
  minE1x5OverE5x5Cut_(params,"minE1x5OverE5x5"),
  minE2x5OverE5x5Cut_(params,"minE2x5OverE5x5"),
  maxNrSatCrysIn5x5Cut_(params,"maxNrSatCrysIn5x5")
{
  contentTags_.emplace("addTrkRef",params.getParameter<edm::InputTag>("addTrkRef"));
}

void GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut::setConsumes(edm::ConsumesCollector& cc) {
  contentTokens_.emplace("addTrkRef",cc.consumes<edm::ValueMap<reco::GsfTrackRef>>(contentTags_["addTrkRef"]));
}

void GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut::getEventContent(const edm::EventBase& ev) {
  ev.getByLabel(contentTags_["addTrkRef"],addTrkHandle_);
}

CutApplicatorBase::result_type GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut::operator() (const reco::GsfElectronPtr& cand) const {
  bool has2ndGsf = (*addTrkHandle_)[cand]!=cand->gsfTrack();

  if (has2ndGsf)
    return true;

  if ( cand->nSaturatedXtals()>maxNrSatCrysIn5x5Cut_(cand) )
    return true;

  const double e5x5 = cand->full5x5_e5x5();
  const double e2x5OverE5x5 = e5x5!=0 ? cand->full5x5_e2x5Max()/e5x5 : 0;
  const double e1x5OverE5x5 = e5x5!=0 ? cand->full5x5_e1x5()/e5x5 : 0;

  return e1x5OverE5x5 > minE1x5OverE5x5Cut_(cand) || e2x5OverE5x5 > minE2x5OverE5x5Cut_(cand);
}

double GsfEleModifiedFull5x5E2x5OverE5x5WithSatCut::
value(const reco::CandidatePtr& cand) const {
  reco::GsfElectronPtr ele(cand);
  //btw we broke somebodies nice model of assuming every cut is 1D....
  //what this is returning is fairly meaningless...
  return ele->full5x5_e1x5() ? ele->full5x5_e2x5Max()/ele->full5x5_e1x5() : 0.;
}
