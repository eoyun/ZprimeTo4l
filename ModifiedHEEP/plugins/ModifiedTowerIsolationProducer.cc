//*****************************************************************************
// File:      ModifiedTowerIsolationProducer.h
// ----------------------------------------------------------------------------
// OrigAuth:  Matthias Mozer
// Institute: IIHE-VUB
//=============================================================================
//*****************************************************************************

// -*- C++ -*-
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedTowerIsolation.h"

// Framework
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
#include "DataFormats/RecoCandidate/interface/RecoCandidate.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include "DataFormats/CaloTowers/interface/CaloTowerCollection.h"

#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedEleTkIsolFromCands.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/TrackReco/interface/TrackBase.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

//
// class declaration
//

class ModifiedTowerIsolationProducer : public edm::stream::EDProducer<> {
   public:
      explicit ModifiedTowerIsolationProducer(const edm::ParameterSet&);
      ~ModifiedTowerIsolationProducer() override;


      void produce(edm::Event&, const edm::EventSetup&) override;
   private:
     template <typename T> void setToken(edm::EDGetTokenT<T>& token,edm::InputTag tag){token=consumes<T>(tag);}
     template <typename T> void setToken(edm::EDGetTokenT<T>& token,const edm::ParameterSet& iPara,const std::string& tag){token=consumes<T>(iPara.getParameter<edm::InputTag>(tag));}
     template <typename T> void setToken(std::vector<edm::EDGetTokenT<T> >& tokens,const edm::ParameterSet& iPara,const std::string& tagName){
       auto tags =iPara.getParameter<std::vector<edm::InputTag> >(tagName);
       for(auto& tag : tags) {
         edm::EDGetTokenT<T> token;
         setToken(token,tag);
         tokens.push_back(token);
       }
     }
      // ----------member data ---------------------------

  edm::EDGetTokenT<edm::View<reco::GsfElectron>> emObjectProducer_;
  edm::EDGetTokenT<CaloTowerCollection> towerProducer_;
  edm::EDGetTokenT<reco::GsfTrackCollection> gsfTrkToken_;
  // edm::InputTag emObjecProducer_;
  // edm::InputTag towerProducer_;

  double egHcalIsoPtMin_;
  double egHcalIsoConeSizeOut_;
  double egHcalIsoConeSizeIn_;
  signed int egHcalDepth_;

  ModifiedEleTkIsolFromCands trkIsoCalc_;

  edm::ParameterSet conf_;

};

ModifiedTowerIsolationProducer::ModifiedTowerIsolationProducer(const edm::ParameterSet& config) :
trkIsoCalc_(config.getParameter<edm::ParameterSet>("trkIsoConfig")),
conf_(config)
{
 // use configuration file to setup input/output collection names
  setToken(emObjectProducer_,config,"emObjectProducer");
  setToken(towerProducer_,config,"towerProducer");
  setToken(gsfTrkToken_,config,"gsfTrks");
  // emObjecProducer_ = conf_.getParameter<edm::InputTag>("emObjectProducer");
  // towerProducer_ = conf_.getParameter<edm::InputTag>("towerProducer");

  egHcalIsoPtMin_               = conf_.getParameter<double>("etMin");
  egHcalIsoConeSizeIn_          = conf_.getParameter<double>("intRadius");
  egHcalIsoConeSizeOut_         = conf_.getParameter<double>("extRadius");
  egHcalDepth_                  = conf_.getParameter<int>("Depth");


  //register your products
  produces < edm::ValueMap<double> >("HcalDepth1TowerSumEt");
}


ModifiedTowerIsolationProducer::~ModifiedTowerIsolationProducer(){}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
ModifiedTowerIsolationProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get the  filtered objects
  edm::Handle< edm::View<reco::GsfElectron> > emObjectHandle;
  iEvent.getByToken(emObjectProducer_,emObjectHandle);
  // iEvent.getByLabel(emObjecProducer_, emObjectHandle);

  // Get the barrel hcal hits
  edm::Handle<CaloTowerCollection> towerHandle;
  iEvent.getByToken(towerProducer_, towerHandle);
  // iEvent.getByLabel(towerProducer_, towerHandle);

  edm::Handle<reco::GsfTrackCollection> gsfTrkHandle;
  iEvent.getByToken(gsfTrkToken_, gsfTrkHandle);

  const CaloTowerCollection* towers = towerHandle.product();

  auto isoMap = std::make_unique<edm::ValueMap<double>>();
  edm::ValueMap<double>::Filler filler(*isoMap);
  std::vector<double> retV(emObjectHandle->size(),0);

  ModifiedTowerIsolation myHadIsolation(egHcalIsoConeSizeOut_,
			      egHcalIsoConeSizeIn_,
			      egHcalIsoPtMin_,
			      egHcalDepth_,
			      towers) ;


  for( size_t i = 0 ; i < emObjectHandle->size(); ++i) {
    bool addGsfTrkSel = false;
    auto additionalGsfTrk = trkIsoCalc_.additionalGsfTrkSelector(*emObjectHandle->ptrAt(i),gsfTrkHandle, addGsfTrkSel);

    double isoValue = myHadIsolation.getTowerEtSum(&(emObjectHandle->at(i)),*(additionalGsfTrk.get()));
    retV[i]=isoValue;
  }

  filler.insert(emObjectHandle,retV.begin(),retV.end());
  filler.fill();
  iEvent.put(std::move(isoMap),"HcalDepth1TowerSumEt");

}

//define this as a plug-in
DEFINE_FWK_MODULE(ModifiedTowerIsolationProducer);
