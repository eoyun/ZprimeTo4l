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

  const edm::EDGetTokenT<edm::View<pat::Electron>> eleToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  const edm::EDGetTokenT<edm::ValueMap<int>> nrSatCrysMapToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<double> rhoToken_;

  const edm::FileInPath xgbPathDR1Et2EB_;
  const edm::FileInPath xgbPathDR2Et1EB_;
  const edm::FileInPath xgbPathDR2Et2EB_;
  const edm::FileInPath xgbPathDR1Et2EE_;
  const edm::FileInPath xgbPathDR2Et1EE_;
  const edm::FileInPath xgbPathDR2Et2EE_;
  const edm::FileInPath xgbPathBkgEt2EB_;
  const edm::FileInPath meanstdPathDR1Et2EB_;
  const edm::FileInPath meanstdPathDR2Et1EB_;
  const edm::FileInPath meanstdPathDR2Et2EB_;
  const edm::FileInPath meanstdPathDR1Et2EE_;
  const edm::FileInPath meanstdPathDR2Et1EE_;
  const edm::FileInPath meanstdPathDR2Et2EE_;
  const edm::FileInPath meanstdPathBkgEt2EB_;

  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR1Et2EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et1EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et2EB_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR1Et2EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et1EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorDR2Et2EE_;
  const std::unique_ptr<MergedMvaEstimator> xgbEstimatorBkgEt2EB_;

  const double cutDR1Et2EB_;
  const double cutDR2Et1EB_;
  const double cutDR2Et2EB_;
  const double cutDR1Et2EE_;
  const double cutDR2Et1EE_;
  const double cutDR2Et2EE_;
  const double cutBkgEt2EB_;
  const double etThresEB_;
  const double etThresEE_;
  const double minEt_;

  bool vetoConv_;
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
  xgbPathDR1Et2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR1Et2EB")),
  xgbPathDR2Et1EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et1EB")),
  xgbPathDR2Et2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et2EB")),
  xgbPathDR1Et2EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR1Et2EE")),
  xgbPathDR2Et1EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et1EE")),
  xgbPathDR2Et2EE_(iConfig.getParameter<edm::FileInPath>("xgbPathDR2Et2EE")),
  xgbPathBkgEt2EB_(iConfig.getParameter<edm::FileInPath>("xgbPathBkgEt2EB")),
  meanstdPathDR1Et2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR1Et2EB")),
  meanstdPathDR2Et1EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et1EB")),
  meanstdPathDR2Et2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et2EB")),
  meanstdPathDR1Et2EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR1Et2EE")),
  meanstdPathDR2Et1EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et1EE")),
  meanstdPathDR2Et2EE_(iConfig.getParameter<edm::FileInPath>("meanstdPathDR2Et2EE")),
  meanstdPathBkgEt2EB_(iConfig.getParameter<edm::FileInPath>("meanstdPathBkgEt2EB")),
  xgbEstimatorDR1Et2EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR1Et2EB_,meanstdPathDR1Et2EB_)),
  xgbEstimatorDR2Et1EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et1EB_,meanstdPathDR2Et1EB_)),
  xgbEstimatorDR2Et2EB_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et2EB_,meanstdPathDR2Et2EB_)),
  xgbEstimatorDR1Et2EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR1Et2EE_,meanstdPathDR1Et2EE_)),
  xgbEstimatorDR2Et1EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et1EE_,meanstdPathDR2Et1EE_)),
  xgbEstimatorDR2Et2EE_(std::make_unique<MergedMvaEstimator>(xgbPathDR2Et2EE_,meanstdPathDR2Et2EE_)),
  xgbEstimatorBkgEt2EB_(std::make_unique<MergedMvaEstimator>(xgbPathBkgEt2EB_,meanstdPathBkgEt2EB_)),
  cutDR1Et2EB_(iConfig.getParameter<double>("cutDR1Et2EB")),
  cutDR2Et1EB_(iConfig.getParameter<double>("cutDR2Et1EB")),
  cutDR2Et2EB_(iConfig.getParameter<double>("cutDR2Et2EB")),
  cutDR1Et2EE_(iConfig.getParameter<double>("cutDR1Et2EE")),
  cutDR2Et1EE_(iConfig.getParameter<double>("cutDR2Et1EE")),
  cutDR2Et2EE_(iConfig.getParameter<double>("cutDR2Et2EE")),
  cutBkgEt2EB_(iConfig.getParameter<double>("cutBkgEt2EB")),
  etThresEB_(iConfig.getParameter<double>("etThresEB")),
  etThresEE_(iConfig.getParameter<double>("etThresEE")),
  minEt_(iConfig.getParameter<double>("minEt")),
  vetoConv_(iConfig.getParameter<bool>("vetoConv")) {
  produces<edm::ValueMap<int>>("cutflowModifiedHEEP");
  produces<edm::ValueMap<int>>("cutflowHEEP");
  produces<edm::ValueMap<int>>("statusMergedElectron");
  produces<edm::ValueMap<float>>("mvaMergedElectron");
  produces<edm::ValueMap<int>>("GSFtypeMergedElectron");
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
  std::vector<int> status_mergedElectron;
  status_mergedElectron.reserve(nEle);
  std::vector<float> mvas_mergedElectron;
  mvas_mergedElectron.reserve(nEle);
  std::vector<int> GSFtype_mergedElectron;
  GSFtype_mergedElectron.reserve(nEle);

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    auto aEle = eleHandle->ptrAt(idx);
    auto orgGsfTrk = aEle->gsfTrack();
    auto addGsfTrk = (*addGsfTrkMap)[aEle];
    double et = MergedLeptonIDs::Et(aEle);

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
    MergedLeptonIDs::GSFtype aGSFtype_mergedElectron = MergedLeptonIDs::GSFtype::nulltype;
    double amvascore_mergedElectron = -1.;

    if ( MergedLeptonIDs::isSameGsfTrack(addGsfTrk,orgGsfTrk) ) {
      // no additional GSF track
      if ( acutflow_HEEP!=MergedLeptonIDs::cutflowElectron::showerShape )
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::failedHEEP;
      else if ( vetoConv_ && !aEle->passConversionVeto() )
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::failConvVeto;
      else {
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::no2ndGsf;

        bool passedMVA = false;

        if ( std::abs(aEle->eta()) < 1.5 && et > etThresEB_ ) {
          amvascore_mergedElectron = xgbEstimatorBkgEt2EB_->computeMva(aEle);
          passedMVA = amvascore_mergedElectron > cutBkgEt2EB_;
        } else
          acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::outsideDef;

        if (passedMVA)
          acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::passedMVA2;
      } // end hasPassedHEEP
    } else {
      // has additional GSF track
      if ( acutflow_modifiedHEEP!=MergedLeptonIDs::cutflowElectron::dxy )
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::failedHEEP;
      else if ( vetoConv_ && !aEle->passConversionVeto() )
        acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::failConvVeto;
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
          aGSFtype_mergedElectron = MergedLeptonIDs::checkElectronGSFtype(aEle,orgGsfTrk,addGsfTrk,etThresEB_,etThresEE_,minEt_);

          switch ( aGSFtype_mergedElectron ) {
            case MergedLeptonIDs::GSFtype::DR1Et2EB:
              amvascore_mergedElectron = xgbEstimatorDR1Et2EB_->computeMva(aEle);
              passedMVA = amvascore_mergedElectron > cutDR1Et2EB_;
              break;
            case MergedLeptonIDs::GSFtype::DR2Et1EB:
              amvascore_mergedElectron = xgbEstimatorDR2Et1EB_->computeMva(aEle);
              passedMVA = amvascore_mergedElectron > cutDR2Et1EB_;
              break;
            case MergedLeptonIDs::GSFtype::DR2Et2EB:
              amvascore_mergedElectron = xgbEstimatorDR2Et2EB_->computeMva(aEle);
              passedMVA = amvascore_mergedElectron > cutDR2Et2EB_;
              break;
            case MergedLeptonIDs::GSFtype::DR1Et2EE:
              amvascore_mergedElectron = xgbEstimatorDR1Et2EE_->computeMva(aEle);
              passedMVA = amvascore_mergedElectron > cutDR1Et2EE_;
              break;
            case MergedLeptonIDs::GSFtype::DR2Et1EE:
              amvascore_mergedElectron = xgbEstimatorDR2Et1EE_->computeMva(aEle);
              passedMVA = amvascore_mergedElectron > cutDR2Et1EE_;
              break;
            case MergedLeptonIDs::GSFtype::DR2Et2EE:
              amvascore_mergedElectron = xgbEstimatorDR2Et2EE_->computeMva(aEle);
              passedMVA = amvascore_mergedElectron > cutDR2Et2EE_;
              break;
            default:
              acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::outsideDef;
              break;
          } // switch

          if (passedMVA)
            acutflow_mergedElectron = MergedLeptonIDs::cutflowElectron::passedMVA1;
        } // end hasMatchedElectron
      } // end isModifiedHEEP
    } // end isSameGsfTrack

    status_mergedElectron.emplace_back( static_cast<int>(acutflow_mergedElectron) );
    mvas_mergedElectron.emplace_back(amvascore_mergedElectron);
    GSFtype_mergedElectron.emplace_back( static_cast<int>(aGSFtype_mergedElectron) );
  } // electron loop

  writeValueMap(iEvent,eleHandle,cutflows_modifiedHEEP,"cutflowModifiedHEEP");
  writeValueMap(iEvent,eleHandle,cutflows_HEEP,"cutflowHEEP");
  writeValueMap(iEvent,eleHandle,status_mergedElectron,"statusMergedElectron");
  writeValueMap(iEvent,eleHandle,mvas_mergedElectron,"mvaMergedElectron");
  writeValueMap(iEvent,eleHandle,GSFtype_mergedElectron,"GSFtypeMergedElectron");
}

DEFINE_FWK_MODULE(MergedLeptonIDProducer);
