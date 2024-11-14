#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"

#include "ZprimeTo4l/Analysis/interface/MuonCorrectionHelper.h"

#include "TTree.h"
#include "TMath.h"

class DimuonControlAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit DimuonControlAnalyzer(const edm::ParameterSet&);
  virtual ~DimuonControlAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  math::PtEtaPhiMLorentzVector lvecFromTuneP(const pat::MuonRef& aMu);

  const bool isMC_;

  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> genptcToken_;
  const edm::EDGetTokenT<double> prefweight_token;
  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<edm::TriggerResults> METfilterToken_;

  const edm::EDGetTokenT<edm::View<pat::Muon>> muonToken_;
  const edm::EDGetTokenT<edm::View<pat::Photon>> photonToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> candToken_;
  const edm::EDGetTokenT<edm::View<pat::PackedCandidate>> lostTrackToken_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;

  const double ptThres_;
  const std::vector<std::string> METfilterList_;
  const std::vector<std::string> trigList_;
  const edm::FileInPath rochesterPath_;
  const edm::FileInPath purwgtPath_;

  MuonCorrectionHelper mucorrHelper_;

  const edm::ESGetToken<TransientTrackBuilder, TransientTrackRecord> ttbToken_;


  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::map<std::string,TH1*> histo1d_;

  TTree* tree_ = nullptr;
  float invM_ = -1.;
  float invMfsr_ = -1.;
  float pt_ = -1.;
  float ptll_ = -1.;
  float eta_ = std::numeric_limits<float>::max();
  float dr_ = -1.;
  float wgt_ = 0.;
  int charge_ = 0;
  int passIso_ = -1;
  int passHighPt_ = -1;
  int passTrkHighPt_ = -1;

  TTree* recoTree_ = nullptr;
  float recoInvM_ = -1.;
  float recoPt_ = -1.;
  float recoPtll_ = -1.;
  float recoEta_ = std::numeric_limits<float>::max();
  float recoDr_ = -1.;
  int recoCharge_ = 0;
  int nTrackerMuon_ = -1;

  const double mumass_ = 0.1056583745;
};

DimuonControlAnalyzer::DimuonControlAnalyzer(const edm::ParameterSet& iConfig) :
isMC_(iConfig.getParameter<bool>("isMC")),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
genptcToken_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genptc"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
METfilterToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("METfilters"))),
muonToken_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
photonToken_(consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("srcPhoton"))),
candToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("srcPackedCand"))),
lostTrackToken_(consumes<edm::View<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("srcLostTracks"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
METfilterList_(iConfig.getParameter<std::vector<std::string>>("METfilterList")),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
rochesterPath_(iConfig.getParameter<edm::FileInPath>("rochesterPath")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
mucorrHelper_(rochesterPath_),
ttbToken_(esConsumes(edm::ESInputTag("","TransientTrackBuilder"))) {
  usesResource("TFileService");
}

math::PtEtaPhiMLorentzVector DimuonControlAnalyzer::lvecFromTuneP(const pat::MuonRef& aMu) {
  return math::PtEtaPhiMLorentzVector(aMu->tunePMuonBestTrack()->pt(),
                                      aMu->tunePMuonBestTrack()->eta(),
                                      aMu->tunePMuonBestTrack()->phi(),
                                      mumass_);
}

void DimuonControlAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["pileup_nPU"] = fs->make<TH1D>("pileup_nPU","pileup_nPU",99,0.,99.);
  histo1d_["pileup_nPUtrue"] = fs->make<TH1D>("pileup_nPUtrue","pileup_nPU",99,0.,99.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);

  histo1d_["cutflow_2M"] = fs->make<TH1D>("cutflow_2M","cutflow",10,0.,10.);

  tree_ = fs->make<TTree>("dimuTree","dimuTree");
  tree_->Branch("invM",&invM_,"invM/F");
  tree_->Branch("invM_FSR",&invMfsr_,"invM_FSR/F");
  tree_->Branch("pt",&pt_,"pt/F");
  tree_->Branch("ptll",&ptll_,"ptll/F");
  tree_->Branch("eta",&eta_,"eta/F");
  tree_->Branch("dr",&dr_,"dr/F");
  tree_->Branch("wgt",&wgt_,"wgt/F");
  tree_->Branch("chargeProduct",&charge_,"chargeProduct/I");
  tree_->Branch("passIso",&passIso_,"passIso/I");
  tree_->Branch("passHighPt",&passHighPt_,"passHighPt/I");
  tree_->Branch("passTrkHighPt",&passTrkHighPt_,"passTrkHighPt/I");

  recoTree_ = fs->make<TTree>("recoTree","recoTree");
  recoTree_->Branch("invM",&recoInvM_,"invM/F");
  recoTree_->Branch("pt",&recoPt_,"pt/F");
  recoTree_->Branch("ptll",&recoPtll_,"ptll/F");
  recoTree_->Branch("eta",&recoEta_,"eta/F");
  recoTree_->Branch("dr",&recoDr_,"dr/F");
  recoTree_->Branch("wgt",&wgt_,"wgt/F");
  recoTree_->Branch("chargeProduct",&recoCharge_,"chargeProduct/I");
  recoTree_->Branch("nTrackerMuon",&nTrackerMuon_,"nTrackerMuon/I");
}

void DimuonControlAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(muonToken_, muonHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  double aWeight = 1.;

  if (isMC_) {
    edm::Handle<double> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    double prefiringweight = *theprefweight;

    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    aWeight = prefiringweight*mcweight/std::abs(mcweight);

    edm::Handle<edm::View<PileupSummaryInfo>> pusummary;
    iEvent.getByToken(pileupToken_, pusummary);

    for (unsigned int idx = 0; idx < pusummary->size(); ++idx) {
      const auto& apu = pusummary->refAt(idx);

      int bx = apu->getBunchCrossing();

      if (bx==0) { // in-time PU only
        aWeight *= purwgt_->GetBinContent( purwgt_->FindBin(apu->getTrueNumInteractions()) );
        histo1d_["pileup_nPU"]->Fill(static_cast<double>(apu->getPU_NumInteractions()),aWeight);
        histo1d_["pileup_nPUtrue"]->Fill(apu->getTrueNumInteractions(),aWeight);

        break;
      }
    }
  }

  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_, trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  const edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    const std::string& trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
        if (trigResultHandle.product()->accept(iTrig))
          isFired = true;
      }
    } // wanted triggers
  } // fired triggers

  std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjs;

  for (unsigned iTrig=0; iTrig < trigObjHandle->size(); iTrig++) {
    const auto& trigObj = trigObjHandle->refAt(iTrig);
    auto trigObjInst = trigObjHandle->at(iTrig); // workaround for copy
    trigObjInst.unpackPathNames(trigList);
    const auto& pathNames = trigObjInst.pathNames();

    for (const auto name : pathNames) {
      for (unsigned int jTrig = 0; jTrig < trigList_.size(); jTrig++) {
        if ( name.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos &&
             trigObjInst.hasPathName(name,true,true) ) {
          trigObjs.push_back(trigObj);
        }
      } // wanted triggers
    } // fired triggers
  } // trigger objs

  histo1d_["cutflow_2M"]->Fill( 0.5, aWeight );

  if (!isFired)
    return;

  histo1d_["cutflow_2M"]->Fill( 1.5, aWeight );

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  edm::Handle<edm::TriggerResults> METfilterHandle;
  iEvent.getByToken(METfilterToken_,METfilterHandle);
  edm::TriggerNames METfilters = iEvent.triggerNames(*METfilterHandle);

  unsigned int nPassedFilters = 0;

  for (unsigned int iTrig = 0; iTrig < METfilterHandle.product()->size(); iTrig++) {
    const std::string trigname = METfilters.triggerName(iTrig);

    if (METfilterHandle.product()->accept(iTrig)) {
      for (const auto& filterName : METfilterList_) {
        if (trigname.find(filterName) != std::string::npos)
          nPassedFilters++;
      }
    }
  }

  if (nPassedFilters!=METfilterList_.size())
    return;

  std::vector<pat::MuonRef> highPtMuons;
  std::vector<pat::MuonRef> highPtTrackerMuons;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    const auto& aMuon = muonHandle->refAt(idx);

    if ( aMuon->tunePMuonBestTrack()->pt() < ptThres_ || std::abs(aMuon->eta()) > 2.4 )
      continue;

    if ( muon::isHighPtMuon(*aMuon,primaryVertex) )
      highPtMuons.push_back(aMuon.castTo<pat::MuonRef>());
    else if ( muon::isTrackerHighPtMuon(*aMuon,primaryVertex) )
      highPtTrackerMuons.push_back(aMuon.castTo<pat::MuonRef>());
  }

  std::vector<pat::MuonRef> tagMuons;

  for (unsigned idx = 0; idx < trigObjs.size(); idx++) {
    const auto& trigObj = trigObjs.at(idx);

    for (const auto& aMu : highPtMuons) {
      if ( aMu->tunePMuonBestTrack()->pt() > 52.
           && reco::deltaR2(trigObj->eta(),trigObj->phi(),
                            aMu->tunePMuonBestTrack()->eta(),aMu->tunePMuonBestTrack()->phi()) < 0.01 )
        tagMuons.push_back(aMu);
    }

    for (const auto& aMu : highPtTrackerMuons) {
      if ( aMu->tunePMuonBestTrack()->pt() > 52.
           && reco::deltaR2(trigObj->eta(),trigObj->phi(),
                            aMu->tunePMuonBestTrack()->eta(),aMu->tunePMuonBestTrack()->phi()) < 0.01 )
        tagMuons.push_back(aMu);
    }
  }

  if (tagMuons.empty())
    return;

  histo1d_["cutflow_2M"]->Fill( 2.5, aWeight );

  std::vector<pat::MuonRef> isolatedHighPtMuons;
  std::vector<pat::MuonRef> isolatedHighPtTrackerMuons;

  auto sortByTuneP = [](const pat::MuonRef& a, const pat::MuonRef& b) {
    return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
  };

  for (unsigned idx = 0; idx < highPtMuons.size(); idx++) {
    const auto& firstMuon = highPtMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 )
      isolatedHighPtMuons.push_back(firstMuon);
  }

  for (unsigned idx = 0; idx < highPtTrackerMuons.size(); idx++) {
    const auto& firstMuon = highPtTrackerMuons.at(idx);
    double firstIso = firstMuon->trackIso();

    std::vector<pat::MuonRef> isoCands;

    for (unsigned jdx = 0; jdx < highPtTrackerMuons.size(); jdx++) {
      const auto& secondMuon = highPtTrackerMuons.at(jdx);

      if (firstMuon==secondMuon)
        continue;

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    for (unsigned jdx = 0; jdx < highPtMuons.size(); jdx++) {
      const auto& secondMuon = highPtMuons.at(jdx);

      if ( mucorrHelper_.checkIso(firstMuon,secondMuon->innerTrack(),*beamSpotHandle) )
        isoCands.push_back(secondMuon);
    }

    std::sort(isoCands.begin(),isoCands.end(),sortByTuneP);

    if (!isoCands.empty())
      firstIso -= isoCands.front()->innerTrack()->pt();

    if ( firstIso/firstMuon->tunePMuonBestTrack()->pt() < 0.1 )
      isolatedHighPtTrackerMuons.push_back(firstMuon);
  }

  std::sort(isolatedHighPtMuons.begin(),isolatedHighPtMuons.end(),sortByTuneP);
  std::sort(isolatedHighPtTrackerMuons.begin(),isolatedHighPtTrackerMuons.end(),sortByTuneP);

  edm::Handle<edm::View<pat::Photon>> photonHandle;
  iEvent.getByToken(photonToken_, photonHandle);

  // isolation loop
  for (const auto& aMu : tagMuons) {
    for (const auto& bMu : highPtMuons) {
      if (aMu==bMu || aMu->tunePMuonBestTrack()->charge()*bMu->tunePMuonBestTrack()->charge() > 0)
        continue;

      const auto lvec1 = lvecFromTuneP(aMu);
      const auto lvec2 = lvecFromTuneP(bMu);
      auto lvecFSR1 = lvec1;
      auto lvecFSR2 = lvec2;

      for (unsigned idx = 0; idx < photonHandle->size(); idx++) {
        const auto& aPho = photonHandle->refAt(idx);

        if ( aPho->pt() < 10. )
          continue;

        if ( !aPho->photonID("cutBasedPhotonID-Fall17-94X-V2-loose") )
          continue;

        if ( reco::deltaR2(lvec1.eta(),lvec1.phi(),aPho->eta(),aPho->phi()) < 0.01 )
          lvecFSR1 += aPho->p4();

        if ( reco::deltaR2(lvec2.eta(),lvec2.phi(),aPho->eta(),aPho->phi()) < 0.01 )
          lvecFSR2 += aPho->p4();
      }

      double invM = (lvec1 + lvec2).M();
      double invMfsr = (lvecFSR1 + lvecFSR2).M();

      if ( (2.5 < invM && invM < 4.5) || (60. < invM && invM < 200.) ) {
        if ( mucorrHelper_.checkIso(bMu,aMu->innerTrack(),*beamSpotHandle) ) {
          invM_ = invM;
          invMfsr_ = invMfsr;
          pt_ = bMu->tunePMuonBestTrack()->pt();
          ptll_ = (lvec1 + lvec2).pt();
          eta_ = bMu->tunePMuonBestTrack()->eta();
          dr_ = reco::deltaR(lvec1.eta(),lvec1.phi(),lvec2.eta(),lvec2.phi());
          charge_ = aMu->tunePMuonBestTrack()->charge()*bMu->tunePMuonBestTrack()->charge();
          wgt_ = aWeight;
          passIso_ = static_cast<int>(std::find(isolatedHighPtMuons.begin(),
                                                isolatedHighPtMuons.end(),
                                                bMu) != isolatedHighPtMuons.end());
          passHighPt_ = 1;
          passTrkHighPt_ = static_cast<int>(muon::isTrackerHighPtMuon(*bMu,primaryVertex));
          tree_->Fill();
        }
      }
    }

    for (const auto& bMu : highPtTrackerMuons) {
      if (aMu==bMu || aMu->tunePMuonBestTrack()->charge()*bMu->tunePMuonBestTrack()->charge() > 0)
        continue;

      const auto lvec1 = lvecFromTuneP(aMu);
      const auto lvec2 = lvecFromTuneP(bMu);
      double invM = (lvec1 + lvec2).M();

      if ( (2.5 < invM && invM < 4.5) || (60. < invM && invM < 200.) ) {
        if ( mucorrHelper_.checkIso(bMu,aMu->innerTrack(),*beamSpotHandle) ) {
          invM_ = invM;
          pt_ = bMu->tunePMuonBestTrack()->pt();
          ptll_ = (lvec1 + lvec2).pt();
          eta_ = bMu->tunePMuonBestTrack()->eta();
          dr_ = reco::deltaR(lvec1.eta(),lvec1.phi(),lvec2.eta(),lvec2.phi());
          charge_ = aMu->tunePMuonBestTrack()->charge()*bMu->tunePMuonBestTrack()->charge();
          wgt_ = aWeight;
          passIso_ = static_cast<int>(std::find(isolatedHighPtTrackerMuons.begin(),
                                                isolatedHighPtTrackerMuons.end(),
                                                bMu) != isolatedHighPtTrackerMuons.end());
          passHighPt_ = 0;
          passTrkHighPt_ = 1;
          tree_->Fill();
        }
      }
    }
  }

  edm::Handle<edm::View<pat::PackedCandidate>> candHandle;
  iEvent.getByToken(candToken_, candHandle);

  edm::Handle<edm::View<pat::PackedCandidate>> lostTrackHandle;
  iEvent.getByToken(lostTrackToken_, lostTrackHandle);

  auto lvec = [this] (const pat::MuonRef& aMu) {
    return math::PtEtaPhiMLorentzVector(aMu->innerTrack()->pt(),
                                        aMu->innerTrack()->eta(),
                                        aMu->innerTrack()->phi(),
                                        mumass_);
  };

  std::vector<pat::PackedCandidateRef> candTracks;

  for (unsigned icand = 0; icand < candHandle->size(); icand++) {
    const auto& acand = candHandle->refAt(icand);

    if ( !acand->hasTrackDetails() )
      continue;

    const reco::Track* atrack = acand->bestTrack();

    if ( atrack->pt() < 20. )
      continue;

    if ( std::abs(atrack->eta()) > 2.4 )
      continue;

    candTracks.push_back(acand.castTo<pat::PackedCandidateRef>());
  }

  for (unsigned icand = 0; icand < lostTrackHandle->size(); icand++) {
    const auto& acand = lostTrackHandle->refAt(icand);

    if ( !acand->hasTrackDetails() )
      continue;

    const reco::Track* atrack = acand->bestTrack();

    if ( atrack->pt() < 20. )
      continue;

    if ( std::abs(atrack->eta()) > 2.4 )
      continue;

    candTracks.push_back(acand.castTo<pat::PackedCandidateRef>());
  }

  // reco loop
  for (const auto& aMu : tagMuons) {
    bool isTagGlobal = false;

    if (!aMu->isTrackerMuon())
      continue;

    for (const auto& bMu : highPtMuons) {
      if (aMu==bMu) {
        isTagGlobal = true;
        break;
      }
    }

    if (!isTagGlobal)
      continue;

    std::vector<pat::MuonRef> trackerMuons;

    for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
      const auto& bMu = muonHandle->refAt(idx).castTo<pat::MuonRef>();

      if (aMu==bMu)
        continue;

      if (!bMu->isTrackerMuon())
        continue;

      if ( bMu->innerTrack()->pt() < 20. || std::abs(bMu->innerTrack()->eta()) > 2.4 )
        continue;

      if ( reco::deltaR2(aMu->innerTrack()->eta(),
                         aMu->innerTrack()->phi(),
                         bMu->innerTrack()->eta(),
                         bMu->innerTrack()->phi()) > 0.01 )
        continue;

      double mass = ( lvec(aMu)+lvec(bMu) ).M();

      if ( 2.5 < mass && mass < 4.0 )
        trackerMuons.push_back(bMu);
    }

    //edm::ESHandle<TransientTrackBuilder> TTbuilder;
    auto fitter = KalmanVertexFitter();
      const TransientTrackBuilder& TTbuilder = iSetup.getData(ttbToken_);

    //iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTbuilder);
    auto firstTrack = TTbuilder.build(aMu->innerTrack());

    std::vector<std::pair<double,pat::PackedCandidateRef>> candFinal;

    for (const auto& acand : candTracks) {
      const reco::Track* atrack = acand->bestTrack();

      if ( reco::deltaR2(aMu->eta(),aMu->phi(),atrack->eta(),atrack->phi()) > 0.01 )
        continue;

      if ( std::isnan(atrack->dzError()) || std::isinf(atrack->dzError()) )
        continue;

      if ( std::isnan(atrack->dxyError()) || std::isinf(atrack->dxyError()) )
        continue;

      if ( std::isnan(atrack->d0Error()) || std::isinf(atrack->d0Error()) )
        continue;

      std::vector<reco::TransientTrack> trackPair;
      trackPair.push_back(firstTrack);
      trackPair.push_back(TTbuilder.build(atrack));

      auto aVtx = fitter.vertex(trackPair);
      double prob = -1.;

      if (aVtx.isValid()) {
        prob = TMath::Prob(aVtx.totalChiSquared(),static_cast<int>(std::rint(aVtx.degreesOfFreedom())));

        if (prob > 0.01) {
          const auto lvecTag = lvec(aMu);
          const auto lvecCand = math::PtEtaPhiMLorentzVector(atrack->pt(),
                                                             atrack->eta(),
                                                             atrack->phi(),
                                                             mumass_);

          const auto p4 = lvecTag + lvecCand;
          double mass = p4.M();

          if (2.5 < mass && mass < 4.0)
            candFinal.push_back(std::make_pair(prob,acand));
        }
      }
    }

    auto sortByProb = [] (const std::pair<double,pat::PackedCandidateRef>& a, const std::pair<double,pat::PackedCandidateRef>&b) {
      return a.first > b.first;
    };

    std::sort(candFinal.begin(),candFinal.end(),sortByProb);

    if (!candFinal.empty()) {
      const auto& acand = candFinal.front();
      const auto lvecTag = lvec(aMu);
      const reco::Track* atrack = acand.second->bestTrack();

      const auto lvecCand = math::PtEtaPhiMLorentzVector(atrack->pt(),
                                                         atrack->eta(),
                                                         atrack->phi(),
                                                         mumass_);

      const auto p4 = lvecTag + lvecCand;

      recoInvM_ = p4.M();
      recoPt_ = atrack->pt();
      recoPtll_ = p4.Pt();
      recoEta_ = atrack->eta();
      recoDr_ = reco::deltaR(lvecTag.eta(),lvecTag.phi(),lvecCand.eta(),lvecCand.phi());
      recoCharge_ = aMu->innerTrack()->charge()*atrack->charge();
      nTrackerMuon_ = trackerMuons.size();
      recoTree_->Fill();
    }
  }

  histo1d_["cutflow_2M"]->Fill( 3.5, aWeight );
  histo1d_["nPV"]->Fill( static_cast<float>(pvHandle->size())+0.5, aWeight );

  return;
}

DEFINE_FWK_MODULE(DimuonControlAnalyzer);
