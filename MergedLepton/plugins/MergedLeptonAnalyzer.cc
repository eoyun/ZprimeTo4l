#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"

#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonIDs.h"
#include "ZprimeTo4l/MergedLepton/interface/MergedLeptonHelper.h"

#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

class MergedLeptonAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonAnalyzer(const edm::ParameterSet&);
  virtual ~MergedLeptonAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  static bool sortByTuneP(const edm::Ptr<reco::Muon> a, const edm::Ptr<reco::Muon> b);

  int categorize(std::vector<edm::Ptr<reco::Muon>>&, const std::vector<edm::Ptr<reco::GenParticle>>&);

  void fillByCategory(std::vector<edm::Ptr<reco::Muon>>& muons,
                      std::vector<edm::Ptr<reco::Muon>>& candidates,
                      const std::vector<edm::Ptr<reco::GenParticle>>& promptLeptons,
                      std::vector<edm::Ptr<reco::GenParticle>>::const_iterator muonItr,
                      const edm::Handle<edm::View<reco::Muon>>& muonHandle,
                      const reco::Vertex& vtx);

  void fillByGsfTrack(std::vector<edm::Ptr<reco::GsfElectron>>& eles,
                      std::vector<edm::Ptr<reco::GenParticle>>::const_iterator eleItr,
                      const edm::Handle<edm::ValueMap<float>>& trkIsoMapHandle,
                      const edm::Handle<edm::ValueMap<float>>& ecalIsoMapHandle,
                      const edm::Handle<edm::ValueMap<int>>& nrSatCrysHandle,
                      const edm::Handle<std::vector<reco::Conversion>>& conversions,
                      const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                      const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkHandle);

  edm::EDGetTokenT<edm::View<reco::Muon>> srcMuon_;
  edm::EDGetTokenT<edm::View<reco::GsfElectron>> srcEle_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<edm::View<pat::MET>> pfMETToken_;
  edm::EDGetTokenT<edm::View<pat::MET>> puppiMETToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;
  edm::EDGetTokenT<edm::ValueMap<float>> ecalIsoToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> nrSatCrysMapToken_;
  edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  edm::EDGetTokenT<double> rhoToken_;
  edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  edm::EDGetTokenT<double> prefweight_token;
  edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  edm::ESHandle<TransientTrackBuilder> TTBuilder_;
  edm::ESHandle<MagneticField> magField_;
  edm::ParameterSet vtxFitterPset_;

  const double ptThres_;
  const double drThres_;

  std::map<std::string,TH1*> histo1d_;

  MergedLeptonHelper aHelper_;
};

MergedLeptonAnalyzer::MergedLeptonAnalyzer(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcEle_(consumes<edm::View<reco::GsfElectron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pfMETToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcPFMET"))),
puppiMETToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcPuppiMET"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
ecalIsoToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("ecalIsoMap"))),
nrSatCrysMapToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("nrSatCrysMap"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
rhoToken_(consumes<double>(iConfig.getParameter<edm::InputTag>("rho"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
vtxFitterPset_(iConfig.getParameter<edm::ParameterSet>("KFParameters")),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")),
aHelper_(vtxFitterPset_) {
  usesResource("TFileService");
}

void MergedLeptonAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  aHelper_.SetFileService(&fs);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",10,0.,10.);
  histo1d_["cutflow_HEEP"] = fs->make<TH1D>("cutflow_HEEP","cutflow_HEEP",15,0.,15.);
  histo1d_["bitmask_METfilter_merged"] = fs->make<TH1D>("bitmask_METfilter_merged","bitmask_METfilter_merged",15,0,15);
  histo1d_["bitmask_METfilter_resolved"] = fs->make<TH1D>("bitmask_METfilter_resolved","bitmask_METfilter_resolved",15,0,15);
  histo1d_["pt_H"] = fs->make<TH1D>("pt_H","p_{T}(H)",200,0.,200.);
  histo1d_["isEcalDriven"] = fs->make<TH1D>("isEcalDriven","isEcalDriven",2,0.,2.);
  histo1d_["isTrackerDriven"] = fs->make<TH1D>("isTrackerDriven","isTrackerDriven",2,0.,2.);
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  aHelper_.initMuonTree("MuonStruct","solo","muon");
  aHelper_.initMuonTree("MuonStruct","tag","muon");
  aHelper_.initMuonTree("MuonStruct","probe","muon");
  aHelper_.initMuonTree("MuonStruct","highPt1","muon");
  aHelper_.initMuonTree("MuonStruct","highPt2","muon");

  aHelper_.initMETTree("METstruct","pf","MET");
  aHelper_.initMETTree("METstruct","puppi","MET");
  aHelper_.initMETTree("METstruct","stdPf","MET");
  aHelper_.initMETTree("METstruct","stdPuppi","MET");

  aHelper_.initElectronTree("ElectronStruct","heep1","el");
  aHelper_.initElectronTree("ElectronStruct","heep2","el");
  aHelper_.initElectronTree("ElectronStruct","mergedEl1","el");
  aHelper_.initElectronTree("ElectronStruct","mergedEl2","el");

  aHelper_.initAddGsfTree("AddGsfStruct","heep1Gsf","addGsf");
  aHelper_.initAddGsfTree("AddGsfStruct","heep2Gsf","addGsf");
  aHelper_.initAddGsfTree("AddGsfStruct","mergedEl1Gsf","addGsf");
}

void MergedLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<reco::GsfElectron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::View<pat::MET>> pfMEThandle;
  iEvent.getByToken(pfMETToken_, pfMEThandle);

  edm::Handle<edm::View<pat::MET>> puppiMEThandle;
  iEvent.getByToken(puppiMETToken_, puppiMEThandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> ecalIsoMapHandle;
  iEvent.getByToken(ecalIsoToken_, ecalIsoMapHandle);

  edm::Handle<edm::ValueMap<int>> nrSatCrysHandle;
  iEvent.getByToken(nrSatCrysMapToken_, nrSatCrysHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<double> rhoHandle;
  iEvent.getByToken(rhoToken_,rhoHandle);

  edm::Handle<reco::ConversionCollection> conversionsHandle;
  iEvent.getByToken(conversionsToken_, conversionsHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  edm::Handle<double> theprefweight;
  iEvent.getByToken(prefweight_token, theprefweight);
  double prefiringweight = *theprefweight;

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken_, genInfo);
  double mcweight = genInfo->weight();
  
  double aWeight = prefiringweight*mcweight/std::abs(mcweight);
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  edm::Handle<edm::TriggerResults> resultHandle;
  iEvent.getByToken(triggerToken_,resultHandle);

  aHelper_.SetMCweight(mcweight);
  aHelper_.SetPrefiringWeight(prefiringweight);

  iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",TTBuilder_);
  iSetup.get<IdealMagneticFieldRecord>().get(magField_);
  aHelper_.SetTTBuilder(&TTBuilder_);

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty()) {
    primaryVertex = pvHandle->ptrAt(0);
    aHelper_.SetPV(primaryVertex);
  } else
    return;

  std::vector<edm::Ptr<reco::GenParticle>> promptLeptons;
  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    auto genptc = genptcHandle->ptrAt(idx);

    if ( genptc->isHardProcess() && genptc->pdgId()==35 )
      histo1d_["pt_H"]->Fill(genptc->pt(),aWeight);

    if ( ( std::abs(genptc->pdgId())==11 || std::abs(genptc->pdgId())==13 ) && genptc->isPromptFinalState() && genptc->pt() > ptThres_ )
      promptLeptons.push_back(genptc);
  }

  if (promptLeptons.size()!=4)
    return;

  edm::TriggerNames triglist = iEvent.triggerNames(*resultHandle);

  unsigned int bitmask_METfilter = 0;

  for (unsigned int iTrig = 0; iTrig < resultHandle.product()->size(); iTrig++) {
    const std::string trigname = triglist.triggerName(iTrig);

    if (resultHandle.product()->accept(iTrig)) {
      if (trigname.find("Flag_goodVertices") != std::string::npos)
        bitmask_METfilter |= 1;
      if (trigname.find("Flag_globalSuperTightHalo2016Filter") != std::string::npos)
        bitmask_METfilter |= (1<<1);
      if (trigname.find("Flag_HBHENoiseFilter") != std::string::npos)
        bitmask_METfilter |= (1<<2);
      if (trigname.find("Flag_HBHENoiseIsoFilter") != std::string::npos)
        bitmask_METfilter |= (1<<3);
      if (trigname.find("EcalDeadCellTriggerPrimitiveFilter") != std::string::npos)
        bitmask_METfilter |= (1<<4);
      if (trigname.find("Flag_BadPFMuonFilter") != std::string::npos)
        bitmask_METfilter |= (1<<5);
      if (trigname.find("Flag_BadPFMuonDzFilter") != std::string::npos)
        bitmask_METfilter |= (1<<6);
      if (trigname.find("Flag_eeBadScFilter") != std::string::npos)
        bitmask_METfilter |= (1<<7);
    }
  }

  auto sisterLambda = [](const edm::Ptr<reco::GenParticle> a, const edm::Ptr<reco::GenParticle> b) {
    return (a->mother() == b->mother()) ? (a->pt() > b->pt()) : (a->mother() < b->mother());
  };

  std::sort(promptLeptons.begin(),promptLeptons.end(),sisterLambda);

  std::vector<edm::Ptr<reco::Muon>> highptMuons1, highptMuons2;
  std::vector<edm::Ptr<reco::GsfElectron>> heeps1, heeps2;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    auto aMuon = muonHandle->ptrAt(idx);

    if ( !(aMuon->pt() > ptThres_) || !muon::isHighPtMuon(*aMuon,*primaryVertex) )
      continue;

    bool matched = false;
    size_t igen = 0;

    for (; igen < promptLeptons.size(); ++igen) {
      auto genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if (matched)
      (igen < 2) ? highptMuons1.push_back(aMuon) : highptMuons2.push_back(aMuon);
  }

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    auto aEle = eleHandle->ptrAt(idx);
    int cutflow = 0;

    if ( !(aEle->pt() > ptThres_) )
      continue;

    bool matched = false;
    size_t igen = 0;

    for (; igen < promptLeptons.size(); ++igen) {
      auto genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aEle->eta(),aEle->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if ( !matched )
      continue;

    histo1d_["isEcalDriven"]->Fill(static_cast<float>(aEle->ecalDrivenSeed())+0.5,aWeight);
    histo1d_["isTrackerDriven"]->Fill(static_cast<float>(aEle->trackerDrivenSeed())+0.5,aWeight);

    bool passModifiedHEEP =
      MergedLeptonIDs::isModifiedHEEP(*aEle,
                                      *primaryVertex,
                                      (*trkIsoMapHandle)[aEle],
                                      (*ecalIsoMapHandle)[aEle],
                                      (*nrSatCrysHandle)[aEle],
                                      *rhoHandle,
                                      cutflow);
    histo1d_["cutflow"]->Fill(static_cast<float>(cutflow)+0.5,aWeight);

    int cutflow_HEEP = 0;
    MergedLeptonIDs::hasPassedHEEP(*aEle,
                                   *primaryVertex,
                                   (*nrSatCrysHandle)[aEle],
                                   *rhoHandle,
                                   cutflow_HEEP);
    histo1d_["cutflow_HEEP"]->Fill(static_cast<float>(cutflow_HEEP)+0.5,aWeight);

    if ( !passModifiedHEEP )
      return; // WARNING veto events

    if (matched)
      (igen < 2) ? heeps1.push_back(aEle) : heeps2.push_back(aEle);
  }

  bool noEle = true;
  bool noMu = true;

  // assume 4e or 4mu only
  if ( heeps1.size() > 0 && heeps2.size() > 0 )
    noEle = false;

  if ( highptMuons1.size() > 0 && highptMuons2.size() > 0 )
    noMu = false;

  if ( noEle && noMu )
    return;

  std::vector<edm::Ptr<reco::Muon>> candidates1, candidates2;
  std::vector<edm::Ptr<reco::Muon>> finalMuons;

  auto pushMuon = [](std::vector<edm::Ptr<reco::Muon>>& fnal,
                     const std::vector<edm::Ptr<reco::Muon>>& highpt,
                     const std::vector<edm::Ptr<reco::Muon>>& cands) {
    int ncand = fnal.size();
    fnal.push_back(highpt.front());

    if ( highpt.size() > 1 )
      fnal.push_back(highpt.at(1));
    else {
      if ( cands.size() > 0 )
        fnal.push_back(cands.front());
    }

    return static_cast<int>(fnal.size()) - ncand;
  };

  int nmuon = 0, nCand1 = 0, nCand2 = 0;

  if ( !noMu ) {
    fillByCategory(highptMuons1,candidates1,promptLeptons,promptLeptons.begin(),muonHandle,*primaryVertex);
    fillByCategory(highptMuons2,candidates2,promptLeptons,std::next(promptLeptons.begin(),2),muonHandle,*primaryVertex);
    nCand1 = pushMuon(finalMuons,highptMuons1, candidates1);
    nCand2 = pushMuon(finalMuons,highptMuons2, candidates2);
    nmuon = nCand1 + nCand2;
  }

  if ( !noEle ) {
    fillByGsfTrack(heeps1,
                   promptLeptons.begin(),
                   trkIsoMapHandle,
                   ecalIsoMapHandle,
                   nrSatCrysHandle,
                   conversionsHandle,
                   beamSpotHandle,
                   addGsfTrkHandle);
    fillByGsfTrack(heeps2,
                   std::next(promptLeptons.begin(),2),
                   trkIsoMapHandle,
                   ecalIsoMapHandle,
                   nrSatCrysHandle,
                   conversionsHandle,
                   beamSpotHandle,
                   addGsfTrkHandle);
  }

  if ( nmuon < 4 && !noMu ) {
    edm::Ptr<reco::Muon> mergedMuon;

    if ( nmuon < 3 )
      mergedMuon = ( highptMuons1.front()->tunePMuonBestTrack()->pt() < highptMuons2.front()->tunePMuonBestTrack()->pt() )
                   ? highptMuons1.front() : highptMuons2.front();
    else
      mergedMuon = ( nCand1 < 2 ) ? highptMuons1.front() : highptMuons2.front();

    float tunepPt = mergedMuon->tunePMuonBestTrack()->pt();
    auto matched = ( mergedMuon==highptMuons1.front() ) ? promptLeptons.begin() : std::next(promptLeptons.begin(),2);
    auto probe = std::next(matched,1);

    if ( std::abs(tunepPt-(*matched)->pt()) > std::abs(tunepPt-(*probe)->pt()) ) {
      probe = matched;
      matched = std::next(matched,1);
    }

    for (unsigned int idx = 0; idx < 8; idx++) {
      if ( bitmask_METfilter & (1<<idx) )
        histo1d_["bitmask_METfilter_merged"]->Fill(static_cast<float>(idx)+0.5,aWeight);
    }

    aHelper_.fillMETs(mergedMuon, *probe, pfMEThandle->at(0), finalMuons, bitmask_METfilter, "pf");
    aHelper_.fillMETs(mergedMuon, *probe, puppiMEThandle->at(0), finalMuons, bitmask_METfilter, "puppi");
  } else if ( nmuon==4 || noMu ) {
    std::vector<edm::Ptr<reco::Candidate>> finalCands;
    auto fillCandidate = [&finalCands](const edm::Ptr<reco::Candidate>& cand) { finalCands.push_back(cand); };
    std::for_each(finalMuons.begin(),finalMuons.end(),fillCandidate);
    std::for_each(heeps1.begin(),heeps1.end(),fillCandidate);
    std::for_each(heeps2.begin(),heeps2.end(),fillCandidate);

    if ( finalCands.empty() )
      return;

    auto findNN = [&finalCands](const pat::MET& met) {
      edm::Ptr<reco::Candidate> candidate;
      double dphicand = DBL_MAX;

      for (auto aLepton : finalCands) {
        double dphi = reco::deltaPhi(aLepton->phi(),met.phi());

        if (dphi < dphicand)
          candidate = aLepton;
      }

      return candidate;
    };

    for (unsigned int idx = 0; idx < 8; idx++) {
      if ( bitmask_METfilter & (1<<idx) )
        histo1d_["bitmask_METfilter_resolved"]->Fill(static_cast<float>(idx)+0.5,aWeight);
    }

    aHelper_.fillMETs(findNN(pfMEThandle->at(0)),pfMEThandle->at(0), finalCands, bitmask_METfilter, "stdPf");
    aHelper_.fillMETs(findNN(puppiMEThandle->at(0)),puppiMEThandle->at(0), finalCands, bitmask_METfilter, "stdPuppi");
  }
}

int MergedLeptonAnalyzer::categorize(std::vector<edm::Ptr<reco::Muon>>& muons, const std::vector<edm::Ptr<reco::GenParticle>>& promptLeptons) {
  if ( muons.size() > 1 )
    std::sort(muons.begin(),muons.end(),sortByTuneP);
  else if ( muons.size() == 1 ) {
    int nMatched = 0;
    const auto highptmu = muons.front();

    for (const auto genptc : promptLeptons) {
      double dr2 = reco::deltaR2(highptmu->eta(),highptmu->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres_*drThres_)
        nMatched++;
    }

    if ( nMatched < 2 )
      return -1; // not merged
  }

  return static_cast<int>(muons.size());
}

void MergedLeptonAnalyzer::fillByCategory(std::vector<edm::Ptr<reco::Muon>>& muons,
  std::vector<edm::Ptr<reco::Muon>>& candidates,
  const std::vector<edm::Ptr<reco::GenParticle>>& promptLeptons,
  std::vector<edm::Ptr<reco::GenParticle>>::const_iterator muonItr,
  const edm::Handle<edm::View<reco::Muon>>& muonHandle,
  const reco::Vertex& vtx) {
  switch (categorize(muons,promptLeptons)) {
    case -1:
      return;
    case 0:
      return;
    case 1: {
      auto aMuon = muons.front();
      float tunepPt = aMuon->tunePMuonBestTrack()->pt();

      auto matched = muonItr;
      auto probe = std::next(muonItr,1);

      if ( std::abs(tunepPt-(*matched)->pt()) > std::abs(tunepPt-(*probe)->pt()) ) {
        probe = matched;
        matched = std::next(matched,1);
      }

      for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
        auto otherMuon = muonHandle->ptrAt(idx);

        if ( aMuon==otherMuon )
          continue;

        if ( otherMuon->pt() > ptThres_ ) {
          double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),otherMuon->eta(),otherMuon->phi());

          if ( dr2 < drThres_*drThres_ )
            candidates.push_back(otherMuon);
        }
      }

      if (candidates.size()==0) {
        aHelper_.fillMuons(aMuon,*matched,*probe,"solo");
        return;
      }

      std::sort(candidates.begin(),candidates.end(),sortByTuneP);
      aHelper_.fillMuons(aMuon,*matched,*probe,"tag");
      aHelper_.fillMuons(candidates.front(),*probe,*matched,"probe");

      return;
    } default:
      aHelper_.fillMuons(muons.front(),*muonItr,*std::next(muonItr,1),"highPt1");
      aHelper_.fillMuons(muons.at(1),*std::next(muonItr,1),*muonItr,"highPt2");
      break;
  }

  return;
}

void MergedLeptonAnalyzer::fillByGsfTrack(std::vector<edm::Ptr<reco::GsfElectron>>& eles,
                                          std::vector<edm::Ptr<reco::GenParticle>>::const_iterator eleItr,
                                          const edm::Handle<edm::ValueMap<float>>& trkIsoMapHandle,
                                          const edm::Handle<edm::ValueMap<float>>& ecalIsoMapHandle,
                                          const edm::Handle<edm::ValueMap<int>>& nrSatCrysHandle,
                                          const edm::Handle<std::vector<reco::Conversion>>& conversions,
                                          const edm::Handle<reco::BeamSpot>& beamSpotHandle,
                                          const edm::Handle<edm::ValueMap<reco::GsfTrackRef>>& addGsfTrkHandle) {
  if ( eles.size()==1 ) {
    auto addGsfTrk = (*addGsfTrkHandle)[eles.front()];
    auto orgGsfTrk = eles.front()->gsfTrack();

    if ( addGsfTrk->pt()==orgGsfTrk->pt() && addGsfTrk->eta()==orgGsfTrk->eta() && addGsfTrk->phi()==orgGsfTrk->phi() ) {
      aHelper_.fillElectrons(eles.front(),
                             (*trkIsoMapHandle)[eles.front()],
                             (*ecalIsoMapHandle)[eles.front()],
                             (*nrSatCrysHandle)[eles.front()],
                             conversions,
                             beamSpotHandle,
                             (*eleItr)->pt()+(*std::next(eleItr,1))->pt(),
                             "mergedEl2");
    } else {
      aHelper_.fillElectrons(eles.front(),
                             (*trkIsoMapHandle)[eles.front()],
                             (*ecalIsoMapHandle)[eles.front()],
                             (*nrSatCrysHandle)[eles.front()],
                             conversions,
                             beamSpotHandle,
                             (*eleItr)->pt()+(*std::next(eleItr,1))->pt(),
                             "mergedEl1");
      aHelper_.fillGsfTracks((*addGsfTrkHandle)[eles.front()],eles.front()->gsfTrack(),"mergedEl1Gsf");
    }

    return;
  }

  // default case: eles.size() > 1
  std::sort(eles.begin(),eles.end(),[](const edm::Ptr<reco::GsfElectron>& a, const edm::Ptr<reco::GsfElectron>& b) { return a->pt() > b->pt(); });
  aHelper_.fillElectrons(eles.at(0),
                         (*trkIsoMapHandle)[eles.at(0)],
                         (*ecalIsoMapHandle)[eles.at(0)],
                         (*nrSatCrysHandle)[eles.at(0)],
                         conversions,
                         beamSpotHandle,
                         (*eleItr)->pt(),
                         "heep1");
  aHelper_.fillElectrons(eles.at(1),
                         (*trkIsoMapHandle)[eles.at(1)],
                         (*ecalIsoMapHandle)[eles.at(1)],
                         (*nrSatCrysHandle)[eles.at(1)],
                         conversions,
                         beamSpotHandle,
                         (*std::next(eleItr,1))->pt(),
                         "heep2");
  aHelper_.fillGsfTracks((*addGsfTrkHandle)[eles.at(0)],eles.at(0)->gsfTrack(),"heep1Gsf");
  aHelper_.fillGsfTracks((*addGsfTrkHandle)[eles.at(1)],eles.at(1)->gsfTrack(),"heep2Gsf");

  return;
}

bool MergedLeptonAnalyzer::sortByTuneP(const edm::Ptr<reco::Muon> a, const edm::Ptr<reco::Muon> b) {
  return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
}

DEFINE_FWK_MODULE(MergedLeptonAnalyzer);
