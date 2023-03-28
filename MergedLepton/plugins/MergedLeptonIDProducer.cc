#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/MET.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedMvaEstimator.h"

class MergedLeptonIDProducer : public edm::stream::EDProducer<> {
public:
  explicit MergedLeptonIDProducer(const edm::ParameterSet&);
  ~MergedLeptonIDProducer() override {}

private:
  void produce(edm::Event&, const edm::EventSetup&) override;

  template<typename T>
  void writeValueMap( edm::Event& iEvent,
                      const edm::Handle<edm::View<pat::Electron>>& handle,
                      const std::vector<T>& values,
                      const std::string& label );

  const edm::EDGetTokenT<edm::View<pat::Electron>> eleToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;

  const edm::FileInPath xgbPathDR1Et2EB_;
  const edm::FileInPath xgbPathDR2Et1EB_;
  const edm::FileInPath xgbPathDR2Et2EB_;
  const edm::FileInPath xgbPathDR2Et1EE_;
  const edm::FileInPath xgbPathDR2Et2EE_;
  const edm::FileInPath xgbPathBkgEt2EB_;
  const edm::FileInPath meanstdPathDR1Et2EB_;
  const edm::FileInPath meanstdPathDR2Et1EB_;
  const edm::FileInPath meanstdPathDR2Et2EB_;
  const edm::FileInPath meanstdPathDR2Et1EE_;
  const edm::FileInPath meanstdPathDR2Et2EE_;
  const edm::FileInPath meanstdPathBkgEt2EB_;

  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR1Et2EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et1EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et2EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et1EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et2EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorBkgEt2EB_;

  const double etThresEB_;
  const double etThresEE_;
  const double minEt_;
};

MergedLeptonIDProducer::MergedLeptonIDProducer(const edm::ParameterSet& iConfig)
: eleToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
  addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
  xgbPathDR1Et2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR1Et2EB")),
  xgbPathDR2Et1EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et1EB")),
  xgbPathDR2Et2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et2EB")),
  xgbPathDR2Et1EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et1EE")),
  xgbPathDR2Et2EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et2EE")),
  xgbPathBkgEt2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathBkgEt2EB")),
  meanstdPathDR1Et2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR1Et2EB")),
  meanstdPathDR2Et1EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et1EB")),
  meanstdPathDR2Et2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et2EB")),
  meanstdPathDR2Et1EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et1EE")),
  meanstdPathDR2Et2EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et2EE")),
  meanstdPathBkgEt2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathBkgEt2EB")),
  xgbEstimatorDR1Et2EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR1Et2EB_,meanstdPathDR1Et2EB_)),
  xgbEstimatorDR2Et1EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et1EB_,meanstdPathDR2Et1EB_)),
  xgbEstimatorDR2Et2EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et2EB_,meanstdPathDR2Et2EB_)),
  xgbEstimatorDR2Et1EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et1EE_,meanstdPathDR2Et1EE_)),
  xgbEstimatorDR2Et2EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et2EE_,meanstdPathDR2Et2EE_)),
  xgbEstimatorBkgEt2EB_(std::make_unique<MergedMvaEstimator>(xgbPathBkgEt2EB_,meanstdPathBkgEt2EB_)),
  etThresEB_(iConfig.getParameter<double>("etThresEB")),
  etThresEE_(iConfig.getParameter<double>("etThresEE")),
  minEt_(iConfig.getParameter<double>("minEt")) {
  produces<edm::ValueMap<float>>("mvaMergedElectronValues");
  produces<edm::ValueMap<int>>("mvaMergedElectronCategories");
}

template<typename T>
void MergedLeptonIDProducer::writeValueMap( edm::Event& iEvent,
                                            const edm::Handle<edm::View<pat::Electron>>& handle,
                                            const std::vector<T>& values,
                                            const std::string& label ) {
  auto valMap = std::make_unique<edm::ValueMap<T>>();
  typename edm::ValueMap<T>::Filler filler(*valMap);
  filler.insert(handle, values.begin(), values.end());
  filler.fill();
  iEvent.put(std::move(valMap),label);
}

void MergedLeptonIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(eleToken_, eleHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  unsigned int nEle = eleHandle->size();

  std::vector<float> mvas_mergedElectron;
  mvas_mergedElectron.reserve(nEle);
  std::vector<int> GSFtype_mergedElectron;
  GSFtype_mergedElectron.reserve(nEle);

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->refAt(idx);
    const auto& orgGsfTrk = aEle->gsfTrack();
    const auto& addGsfTrk = (*addGsfTrkMap)[aEle];

    // default category is bkgEt2EB (no 2nd GSF)
    MergedLeptonIDs::GSFtype aGSFtype_mergedElectron = MergedLeptonIDs::GSFtype::bkgEt2EB;
    double amvascore_mergedElectron = -1.;

    if ( addGsfTrk==orgGsfTrk ) {
      // no additional GSF track
      const auto castEle = aEle.castTo<pat::ElectronRef>();

      // apply lower Et cut for the bkg enriched CR
      // need to apply an additional offline Et cut afterwards
      if ( std::abs(aEle->superCluster()->eta()) < 1.5 && aEle->et() > minEt_ )
        amvascore_mergedElectron = xgbEstimatorBkgEt2EB_->computeMva(castEle);
    } else {
      // has additional GSF track
      // assume modified HEEP apriori
      // default category is DR1Et2EB
      aGSFtype_mergedElectron = MergedLeptonIDs::GSFtype::DR1Et2EB;

      // additional GSF track must not be an electron
      bool notMerged = false;

      for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
        const auto& secEle = eleHandle->refAt(jdx);
        const auto& secGsfTrk = secEle->gsfTrack();

        if ( addGsfTrk==secGsfTrk ) {
          notMerged = true;
          break;
        }
      }

      if ( !notMerged ) {
        const auto castEle = aEle.castTo<pat::ElectronRef>();
        aGSFtype_mergedElectron = MergedLeptonIDs::checkElectronGSFtype(castEle,orgGsfTrk,addGsfTrk,etThresEB_,etThresEE_,minEt_);

        switch ( aGSFtype_mergedElectron ) {
          case MergedLeptonIDs::GSFtype::DR1Et2EB:
            amvascore_mergedElectron = xgbEstimatorDR1Et2EB_->computeMva(castEle);
            break;
          case MergedLeptonIDs::GSFtype::DR2Et1EB:
            amvascore_mergedElectron = xgbEstimatorDR2Et1EB_->computeMva(castEle);
            break;
          case MergedLeptonIDs::GSFtype::DR2Et2EB:
            amvascore_mergedElectron = xgbEstimatorDR2Et2EB_->computeMva(castEle);
            break;
          case MergedLeptonIDs::GSFtype::DR2Et1EE:
            amvascore_mergedElectron = xgbEstimatorDR2Et1EE_->computeMva(castEle);
            break;
          case MergedLeptonIDs::GSFtype::DR2Et2EE:
            amvascore_mergedElectron = xgbEstimatorDR2Et2EE_->computeMva(castEle);
            break;
          case MergedLeptonIDs::GSFtype::extendedCR:
            // need to apply an additional offline Et cut afterwards
            aGSFtype_mergedElectron = MergedLeptonIDs::GSFtype::DR1Et2EB;
            amvascore_mergedElectron = xgbEstimatorDR1Et2EB_->computeMva(castEle);
            break;
          default:
            aGSFtype_mergedElectron = MergedLeptonIDs::GSFtype::DR1Et2EB;
            amvascore_mergedElectron = -1.;
            break;
        } // switch
      } // !notMerged
    } // end isSameGsfTrack

    mvas_mergedElectron.emplace_back(amvascore_mergedElectron);
    GSFtype_mergedElectron.emplace_back( static_cast<int>(aGSFtype_mergedElectron) );
  } // electron loop

  writeValueMap(iEvent,eleHandle,mvas_mergedElectron,"mvaMergedElectronValues");
  writeValueMap(iEvent,eleHandle,GSFtype_mergedElectron,"mvaMergedElectronCategories");

  return;
}

DEFINE_FWK_MODULE(MergedLeptonIDProducer);
