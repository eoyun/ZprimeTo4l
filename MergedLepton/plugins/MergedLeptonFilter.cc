#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

class MergedLeptonFilter : public edm::EDFilter {
public:
  explicit MergedLeptonFilter(const edm::ParameterSet&);
  ~MergedLeptonFilter() {}

private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool hasMergedMuon(edm::Event&, const edm::EventSetup&);

  // edm::EDGetTokenT<edm::View<reco::GsfElectron>> srcEle_;
  edm::EDGetTokenT<edm::View<reco::Muon>> srcMuon_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;

  const double ptThres_;
  const double drThres_;

  int nPrompt_ = 0;
  int nHighpt_ = 0;
  int nLoose_ = 0;
  int nMerged_ = 0;
};

MergedLeptonFilter::MergedLeptonFilter(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")) {}

bool MergedLeptonFilter::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  return hasMergedMuon(iEvent,iSetup);
}

bool MergedLeptonFilter::hasMergedMuon(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return false;

  std::vector<edm::Ptr<reco::GenParticle>> promptMuons;
  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    auto genptc = genptcHandle->ptrAt(idx);

    if ( std::abs(genptc->pdgId())==13 && genptc->isPromptFinalState() && genptc->pt() > ptThres_ )
      promptMuons.push_back(genptc);
  }

  nPrompt_ += promptMuons.size();

  if (promptMuons.size()==0)
    return false;

  std::vector<edm::Ptr<reco::Muon>> highptMuons;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    auto aMuon = muonHandle->ptrAt(idx);

    if ( !(aMuon->pt() > ptThres_) || !muon::isHighPtMuon(*aMuon,primaryVertex) )
      continue;

    bool matched = false;

    for (const auto genptc : promptMuons) {
      double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if (matched)
      highptMuons.push_back(aMuon);
  }

  nHighpt_ += highptMuons.size();

  if (highptMuons.size()==0)
    return false;

  bool hasMerged = false;

  for (const auto highptmu : highptMuons) {
    int nMatched = 0;

    for (const auto genptc : promptMuons) {
      double dr2 = reco::deltaR2(highptmu->eta(),highptmu->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2)
        nMatched++;
    }

    bool hasOtherMuon = false;

    for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
      auto aMuon = muonHandle->ptrAt(idx);

      if ( aMuon->pt()==highptmu->pt() && aMuon->eta()==highptmu->eta() && aMuon->phi()==highptmu->phi() )
        continue;

      if ( aMuon->pt() > ptThres_ && muon::isLooseMuon(*aMuon) ) {
        double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),highptmu->eta(),highptmu->phi());

        if (dr2 < drThres2) {
          hasOtherMuon = true;
          nLoose_ ++;
          break;
        }
      }
    }

    if ( nMatched >=2 && !hasOtherMuon ) {
      hasMerged = true;
      nMerged_ ++;
    }
  }

  return hasMerged;
}

void MergedLeptonFilter::endJob() {
  std::cout << "nPrompt = " << nPrompt_ << " nHighpt = " << nHighpt_ << " nLoose = " << nLoose_ << " nMerged = " << nMerged_ << std::endl;
}

DEFINE_FWK_MODULE(MergedLeptonFilter);
