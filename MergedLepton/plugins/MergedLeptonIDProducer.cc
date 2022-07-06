#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
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
  void produceElectrons(edm::Event&, const reco::Vertex&, const double& rho);

  template<typename T>
  void writeValueMap( edm::Event& iEvent,
                      const edm::Handle<edm::View<pat::Electron>>& handle,
                      const std::vector<T>& values,
                      const std::string& label );

  enum openingAngle {
    nullAngle = -1,
    DR1EB,
    DR1EE,
    DR2EB,
    DR2EE,
    DR3EB,
    DR3EE
  };

  openingAngle checkElectronOpeningAngle( edm::Ptr<pat::Electron>& aEle,
                                          reco::GsfTrackRef& orgGsfTrk,
                                          reco::GsfTrackRef& addGsfTrk );

  const edm::EDGetTokenT<edm::View<pat::Electron>> eleToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> nrSatCrysMapToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<double> rhoToken_;

  const edm::FileInPath xgbPathDR1EB_;
  const edm::FileInPath xgbPathDR2EB_;
  const edm::FileInPath xgbPathDR3EB_;
  const edm::FileInPath xgbPathDR1EE_;
  const edm::FileInPath xgbPathDR2EE_;
  const edm::FileInPath xgbPathDR3EE_;
  const edm::FileInPath meanstdPathDR1EB_;
  const edm::FileInPath meanstdPathDR2EB_;
  const edm::FileInPath meanstdPathDR3EB_;
  const edm::FileInPath meanstdPathDR1EE_;
  const edm::FileInPath meanstdPathDR2EE_;
  const edm::FileInPath meanstdPathDR3EE_;

  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR1EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR3EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR1EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR3EE_;

  const double cutDR1EB_;
  const double cutDR2EB_;
  const double cutDR3EB_;
  const double cutDR1EE_;
  const double cutDR2EE_;
  const double cutDR3EE_;

  const edm::EDPutTokenT<int> nresolvedElectronToken_;
  const edm::EDPutTokenT<int> nmergedElectronToken_;
};

MergedLeptonIDProducer::MergedLeptonIDProducer(const edm::ParameterSet& iConfig)
: eleToken_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
  addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
  trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
  ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
  nrSatCrysMapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("nrSatCrysMap"))),
  pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
  beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
  rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
  xgbPathDR1EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR1EB")),
  xgbPathDR2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2EB")),
  xgbPathDR3EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR3EB")),
  xgbPathDR1EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR1EE")),
  xgbPathDR2EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2EE")),
  xgbPathDR3EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR3EE")),
  meanstdPathDR1EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR1EB")),
  meanstdPathDR2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2EB")),
  meanstdPathDR3EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR3EB")),
  meanstdPathDR1EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR1EE")),
  meanstdPathDR2EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2EE")),
  meanstdPathDR3EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR3EE")),
  xgbEstimatorDR1EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR1EB_,meanstdPathDR1EB_)),
  xgbEstimatorDR2EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR2EB_,meanstdPathDR2EB_)),
  xgbEstimatorDR3EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR3EB_,meanstdPathDR3EB_)),
  xgbEstimatorDR1EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR1EE_,meanstdPathDR1EE_)),
  xgbEstimatorDR2EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR2EE_,meanstdPathDR2EE_)),
  xgbEstimatorDR3EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR3EE_,meanstdPathDR3EE_)),
  cutDR1EB_(iConfig.getParameter<double>("cutDR1EB")),
  cutDR2EB_(iConfig.getParameter<double>("cutDR2EB")),
  cutDR3EB_(iConfig.getParameter<double>("cutDR3EB")),
  cutDR1EE_(iConfig.getParameter<double>("cutDR1EE")),
  cutDR2EE_(iConfig.getParameter<double>("cutDR2EE")),
  cutDR3EE_(iConfig.getParameter<double>("cutDR3EE")),
  nresolvedElectronToken_(produces<int>("nresolvedElectron")),
  nmergedElectronToken_(produces<int>("nmergedElectron")) {

  produces<edm::ValueMap<int>>("cutflow_modifiedHEEP");
  produces<edm::ValueMap<int>>("cutflow_HEEP");
  produces<edm::ValueMap<int>>("cutflow_mergedElectron");
  produces<edm::ValueMap<float>>("mva_mergedElectron");
  produces<edm::ValueMap<int>>("openingAngle_mergedElectron");
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

MergedLeptonIDProducer::openingAngle MergedLeptonIDProducer::checkElectronOpeningAngle( edm::Ptr<pat::Electron>& aEle,
                                                                                        reco::GsfTrackRef& orgGsfTrk,
                                                                                        reco::GsfTrackRef& addGsfTrk ) {
  double dr2 = reco::deltaR2(orgGsfTrk->eta(),orgGsfTrk->phi(),addGsfTrk->eta(),addGsfTrk->phi());
  auto square = [](double input) { return input*input; };

  if ( aEle->isEB() ) {
    if ( dr2 < square(0.0174) )
      return openingAngle::DR1EB;
    else if ( dr2 > square(0.0174) && dr2 < square(2.*0.0174) )
      return openingAngle::DR2EB;
    else
      return openingAngle::DR3EB;
  } else {
    if ( dr2 < square( 0.00864*std::abs(std::sinh(aEle->eta())) ) )
      return openingAngle::DR1EE;
    else if ( dr2 > square( 0.00864*std::abs(std::sinh(aEle->eta())) ) && dr2 < square( 2.*0.00864*std::abs(std::sinh(aEle->eta())) ) )
      return openingAngle::DR2EE;
    else
      return openingAngle::DR3EE;
  }

  return openingAngle::nullAngle;
}

void MergedLeptonIDProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_,rhoHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    primaryVertex = reco::Vertex(beamSpotHandle->position(),reco::Vertex::Error());

  produceElectrons(iEvent,primaryVertex,*rhoHandle);

  return;
}

void MergedLeptonIDProducer::produceElectrons(edm::Event& iEvent, const reco::Vertex& primaryVertex, const double& rho) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(eleToken_, eleHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkMap;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkMap);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> ecalIsoMapHandle;
  iEvent.getByToken(ecalIsoToken_, ecalIsoMapHandle);

  edm::Handle<edm::ValueMap<int>> nrSatCrysHandle;
  iEvent.getByToken(nrSatCrysMapToken_, nrSatCrysHandle);

  unsigned int nEle = eleHandle->size();

  std::vector<int> cutflows_modifiedHEEP;
  cutflows_modifiedHEEP.reserve(nEle);
  std::vector<int> cutflows_HEEP;
  cutflows_HEEP.reserve(nEle);
  std::vector<int> cutflows_mergedElectron;
  cutflows_mergedElectron.reserve(nEle);
  std::vector<float> mvas_mergedElectron;
  mvas_mergedElectron.reserve(nEle);
  std::vector<int> openingAngles_mergedElectron;
  openingAngles_mergedElectron.reserve(nEle);

  int nresolvedElectron = 0;
  int nmergedElectron = 0;

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    auto aEle = eleHandle->ptrAt(idx);
    auto orgGsfTrk = aEle->gsfTrack();
    auto addGsfTrk = (*addGsfTrkMap)[aEle];

    MergedLeptonIDs::cutflowElectron acutflow_modifiedHEEP = MergedLeptonIDs::cutflowElectron::baseline;

    MergedLeptonIDs::isModifiedHEEP(*aEle,
                                    primaryVertex,
                                    (*trkIsoMapHandle)[aEle],
                                    (*ecalIsoMapHandle)[aEle],
                                    (*nrSatCrysHandle)[aEle],
                                    rho,
                                    acutflow_modifiedHEEP);

    MergedLeptonIDs::cutflowElectron acutflow_HEEP = MergedLeptonIDs::cutflowElectron::baseline;
    MergedLeptonIDs::hasPassedHEEP(*aEle,
                                   primaryVertex,
                                   (*nrSatCrysHandle)[aEle],
                                   rho,
                                   acutflow_HEEP);

    cutflows_modifiedHEEP.emplace_back( static_cast<int>(acutflow_modifiedHEEP) );
    cutflows_HEEP.emplace_back( static_cast<int>(acutflow_HEEP) );

    MergedLeptonIDs::cutflowElectron acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::baseline;
    openingAngle aopenangle_mergedElectron = openingAngle::nullAngle;
    double amvascore_mergedElectron = -1.;

    if ( MergedLeptonIDs::isSameGsfTrack(addGsfTrk,orgGsfTrk) ) {
      // no additional GSF track
      if ( acutflow_HEEP!=MergedLeptonIDs::cutflowElectron::showerShape )
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::failedHEEP;
      else {
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::no2ndGsf;
        // TODO
      } // end hasPassedHEEP
    } else {
      // has additional GSF track
      if ( acutflow_modifiedHEEP!=MergedLeptonIDs::cutflowElectron::dxy )
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::failedHEEP;
      else {
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::has2ndGsf;

        // additional GSF track must not be an electron
        bool notMerged = false;

        for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
          auto secEle = eleHandle->ptrAt(jdx);
          const auto secGsfTrk = secEle->gsfTrack();

          if ( MergedLeptonIDs::isSameGsfTrack(addGsfTrk,secGsfTrk) ) {
            notMerged = true;
            break;
          }
        }

        if ( notMerged )
          acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::hasMatchedElectron;
        else {
          bool passedMVA = false;
          // apply BDT to further reject bkg contributions (e.g. brem)
          aopenangle_mergedElectron = checkElectronOpeningAngle(aEle,orgGsfTrk,addGsfTrk);

          switch ( aopenangle_mergedElectron ) {
            case openingAngle::DR1EB:
              amvascore_mergedElectron = xgbEstimatorDR1EB_->computeMva(aEle,primaryVertex);
              passedMVA = amvascore_mergedElectron > cutDR1EB_;
              break;
            case openingAngle::DR2EB:
              amvascore_mergedElectron = xgbEstimatorDR2EB_->computeMva(aEle,primaryVertex);
              passedMVA = amvascore_mergedElectron > cutDR2EB_;
              break;
            case openingAngle::DR3EB:
              amvascore_mergedElectron = xgbEstimatorDR3EB_->computeMva(aEle,primaryVertex);
              passedMVA = amvascore_mergedElectron > cutDR3EB_;
              break;
            case openingAngle::DR1EE:
              amvascore_mergedElectron = xgbEstimatorDR1EE_->computeMva(aEle,primaryVertex);
              passedMVA = amvascore_mergedElectron > cutDR1EE_;
              break;
            case openingAngle::DR2EE:
              amvascore_mergedElectron = xgbEstimatorDR2EE_->computeMva(aEle,primaryVertex);
              passedMVA = amvascore_mergedElectron > cutDR2EE_;
              break;
            case openingAngle::DR3EE:
              amvascore_mergedElectron = xgbEstimatorDR3EE_->computeMva(aEle,primaryVertex);
              passedMVA = amvascore_mergedElectron > cutDR3EE_;
              break;
            default:
              throw cms::Exception("LogicError") << "opening angle between the original and additional GSF track does not fit into any of MergedLeptonIDProducer::openingAngle categories" << std::endl;
              break;
          } // switch

          if (passedMVA)
            acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::passedMVA1;
        } // end hasMatchedElectron
      } // end isModifiedHEEP
    } // end isSameGsfTrack

    cutflows_mergedElectron.emplace_back( static_cast<int>(acutflow_mergedElectron) );
    mvas_mergedElectron.emplace_back(amvascore_mergedElectron);
    openingAngles_mergedElectron.emplace_back( static_cast<int>(aopenangle_mergedElectron) );

    MergedLeptonIDs::typeElectron atype = MergedLeptonIDs::checkTypeElectron(acutflow_modifiedHEEP,
                                                                             acutflow_HEEP,
                                                                             acutflow_mergedElectron);

    if ( atype==MergedLeptonIDs::typeElectron::resolved )
      nresolvedElectron++;
    if ( atype==MergedLeptonIDs::typeElectron::merged )
      nmergedElectron++;
  } // electron loop

  writeValueMap(iEvent,eleHandle,cutflows_modifiedHEEP,"cutflow_modifiedHEEP");
  writeValueMap(iEvent,eleHandle,cutflows_HEEP,"cutflow_HEEP");
  writeValueMap(iEvent,eleHandle,cutflows_mergedElectron,"cutflow_mergedElectron");
  writeValueMap(iEvent,eleHandle,mvas_mergedElectron,"mva_mergedElectron");
  writeValueMap(iEvent,eleHandle,openingAngles_mergedElectron,"openingAngle_mergedElectron");

  iEvent.emplace(nresolvedElectronToken_,nresolvedElectron);
  iEvent.emplace(nmergedElectronToken_,nmergedElectron);
}

DEFINE_FWK_MODULE(MergedLeptonIDProducer);
