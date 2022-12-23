#include "PhysicsTools/SelectorUtils/interface/CutApplicatorWithEventContentBase.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "RecoEgamma/ElectronIdentification/interface/EBEECutValues.h"

class GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut : public CutApplicatorWithEventContentBase {
public:
  GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut(const edm::ParameterSet& c);

  result_type operator() (const reco::GsfElectronPtr&) const final;

  void setConsumes(edm::ConsumesCollector&) final;
  void getEventContent(const edm::EventBase&) final;

  double value(const reco::CandidatePtr& cand) const final;

  CandidateType candidateType() const final { return ELECTRON; }

private:
  EBEECutValues maxSigmaIEtaIEtaCut_;
  EBEECutValuesInt maxNrSatCrysIn5x5Cut_;
  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addTrkHandle_;
};

DEFINE_EDM_PLUGIN(CutApplicatorFactory,
                  GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut,
                  "GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut");

GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut::GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut(const edm::ParameterSet& params) :
  CutApplicatorWithEventContentBase(params),
  maxSigmaIEtaIEtaCut_(params,"maxSigmaIEtaIEta"),
  maxNrSatCrysIn5x5Cut_(params,"maxNrSatCrysIn5x5")
{
  contentTags_.emplace("addTrkRef",params.getParameter<edm::InputTag>("addTrkRef"));
}

void GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut::setConsumes(edm::ConsumesCollector& cc) {
  contentTokens_.emplace("addTrkRef",cc.consumes<edm::ValueMap<reco::GsfTrackRef>>(contentTags_["addTrkRef"]));
}

void GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut::getEventContent(const edm::EventBase& ev) {
  ev.getByLabel(contentTags_["addTrkRef"],addTrkHandle_);
}

CutApplicatorBase::result_type GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut::operator() (const reco::GsfElectronPtr& cand) const {
  bool has2ndGsf = (*addTrkHandle_)[cand]!=cand->gsfTrack();

  if (has2ndGsf)
    return true;

  if (cand->nSaturatedXtals() > maxNrSatCrysIn5x5Cut_(cand))
    return true;
  else
    return cand->full5x5_sigmaIetaIeta() < maxSigmaIEtaIEtaCut_(cand);
}

double GsfEleModifiedFull5x5SigmaIEtaIEtaWithSatCut::value(const reco::CandidatePtr& cand) const {
  reco::GsfElectronPtr ele(cand);
  return ele->full5x5_sigmaIetaIeta();
}
