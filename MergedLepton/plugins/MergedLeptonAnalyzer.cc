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
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/PatCandidates/interface/MET.h"

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
  void initTObj(edm::Service<TFileService>&, std::string);
  int categorize(std::vector<edm::Ptr<reco::Muon>>&, const std::vector<edm::Ptr<reco::GenParticle>>&);
  bool isHighPtTrackerMuon(const reco::Muon& muon, const reco::Vertex& vtx);
  void fillByCategory(std::vector<edm::Ptr<reco::Muon>>& muons,
    std::vector<edm::Ptr<reco::Muon>>& candidates,
    const std::vector<edm::Ptr<reco::GenParticle>>& promptMuons,
    std::vector<edm::Ptr<reco::GenParticle>>::const_iterator muonItr,
    const edm::Handle<edm::View<reco::Muon>>& muonHandle,
    const reco::Vertex& vtx);
  void fillMuons(const edm::Ptr<reco::Muon>&, const edm::Ptr<reco::GenParticle>&, const reco::Vertex&, std::string);
  void fillMets(const edm::Ptr<reco::Muon>&, const edm::Ptr<reco::GenParticle>&, const pat::MET&, const std::vector<edm::Ptr<reco::Muon>>&, std::string);

  edm::EDGetTokenT<edm::View<reco::Muon>> srcMuon_;
  edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  edm::EDGetTokenT<edm::View<pat::MET>> pfMETToken_;
  edm::EDGetTokenT<edm::View<pat::MET>> puppiMETToken_;

  const double ptThres_;
  const double drThres_;

  typedef struct {
    int numberOfValidTrackerHits, numberOfValidPixelHits, numberOfValidStripHits,
    trackerLayersWithMeasurement, pixelLayersWithMeasurement, stripLayersWithMeasurement,
    trackerLayersWithoutMeasurement, pixelLayersWithoutMeasurement, stripLayersWithoutMeasurement;
    float trackerVoM, pixelVoM, stripVoM, pfPt, tunepPt, genPt, PFoGen, TPoGen;
  } MuonStruct;

  typedef struct {
    float PFphi, PFdPhi, PFpt, PFoGen, PFSumEt, TPphi, TPdPhi, TPpt, TPoGen, TPSumEt, dSumEt;
    int nMuon;
  } METstruct;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TTree*> tree_;
  std::map<std::string,MuonStruct> values_;
  std::map<std::string,METstruct> metvalues_;

  TString mustr_ = TString("numberOfValidTrackerHits/I:numberOfValidPixelHits:numberOfValidStripHits:")
  + "trackerLayersWithMeasurement:pixelLayersWithMeasurement:stripLayersWithMeasurement:"
  + "trackerLayersWithoutMeasurement:pixelLayersWithoutMeasurement:stripLayersWithoutMeasurement:"
  + "trackerVoM/F:pixelVoM:stripVoM:pfPt:tunepPt:genPt:PFoGen:TPoGen";

  TString metstr_ = "PFphi/F:PFdPhi:PFpt:PFoGen:PFSumEt:TPphi:TPdPhi:TPpt:TPoGen:TPSumEt:dSumEt:nMuon/I";
};

MergedLeptonAnalyzer::MergedLeptonAnalyzer(const edm::ParameterSet& iConfig) :
srcMuon_(consumes<edm::View<reco::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pfMETToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcPFMET"))),
puppiMETToken_(consumes<edm::View<pat::MET>>(iConfig.getParameter<edm::InputTag>("srcPuppiMET"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")) {
  usesResource("TFileService");
}

void MergedLeptonAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",10,0.,10.);

  initTObj(fs,"solo");
  initTObj(fs,"tag");
  initTObj(fs,"probe");

  METstruct pfmetStruct, puppimetStruct;
  metvalues_["pfMET"] = pfmetStruct;
  metvalues_["puppiMET"] = puppimetStruct;

  tree_["pfMETTree"] = fs->make<TTree>("pfMETTree","pfMETTree");
  tree_["pfMETTree"]->Branch("METstruct",&(metvalues_["pfMET"]),metstr_);
  tree_["puppiMETTree"] = fs->make<TTree>("puppiMETTree","puppiMETTree");
  tree_["puppiMETTree"]->Branch("METstruct",&(metvalues_["puppiMET"]),metstr_);
}

void MergedLeptonAnalyzer::initTObj(edm::Service<TFileService>& fs, std::string prefix) {
  histo1d_[prefix+"_muType"] = fs->make<TH1D>(TString(prefix)+"_muType","muon type",3,0.,3.);
  histo1d_[prefix+"_algo"] = fs->make<TH1D>(TString(prefix)+"_algo","muon algo",46,0.,46.);
  histo1d_[prefix+"_orgAlgo"] = fs->make<TH1D>(TString(prefix)+"_orgAlgo","muon original algo",46,0.,46.);
  histo1d_[prefix+"_algoMask"] = fs->make<TH1D>(TString(prefix)+"_algoMask","muon algo mask",46,0.,46.);
  histo1d_[prefix+"_quality"] = fs->make<TH1D>(TString(prefix)+"_quality","muon quality",8,0.,8.);
  histo1d_[prefix+"_passIDs"] = fs->make<TH1D>(TString(prefix)+"_passIDs","muon passed IDs",5,0.,5.);

  MuonStruct muStruct;
  values_[prefix+"_muon"] = muStruct;
  tree_[prefix+"_muonTree"] = fs->make<TTree>(TString(prefix)+"_muonTree","muonTree");
  tree_[prefix+"_muonTree"]->Branch("MuonStruct",&(values_[prefix+"_muon"]),mustr_);
}

void MergedLeptonAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::View<pat::MET>> pfMEThandle;
  iEvent.getByToken(pfMETToken_, pfMEThandle);

  edm::Handle<edm::View<pat::MET>> puppiMEThandle;
  iEvent.getByToken(puppiMETToken_, puppiMEThandle);

  reco::Vertex primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->at(0);
  else
    return;

  std::vector<edm::Ptr<reco::GenParticle>> promptMuons;
  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    auto genptc = genptcHandle->ptrAt(idx);

    if ( std::abs(genptc->pdgId())==13 && genptc->isPromptFinalState() && genptc->pt() > ptThres_ )
      promptMuons.push_back(genptc);
  }

  if (promptMuons.size()!=4)
    return;

  auto sisterLambda = [](const edm::Ptr<reco::GenParticle> a, const edm::Ptr<reco::GenParticle> b) {
    return (a->mother() == b->mother()) ? (a->pt() > b->pt()) : (a->mother() < b->mother());
  };

  std::sort(promptMuons.begin(),promptMuons.end(),sisterLambda);

  std::vector<edm::Ptr<reco::Muon>> highptMuons1, highptMuons2;

  for (unsigned int idx = 0; idx < muonHandle->size(); ++idx) {
    auto aMuon = muonHandle->ptrAt(idx);

    if ( !(aMuon->pt() > ptThres_) || !muon::isHighPtMuon(*aMuon,primaryVertex) )
      continue;

    bool matched = false;
    size_t igen = 0;

    for (; igen < promptMuons.size(); ++igen) {
      auto genptc = promptMuons.at(igen);
      double dr2 = reco::deltaR2(aMuon->eta(),aMuon->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if (matched)
      (igen < 2) ? highptMuons1.push_back(aMuon) : highptMuons2.push_back(aMuon);
  }

  if ( highptMuons1.size() < 1 || highptMuons2.size() < 1 )
    return;

  std::vector<edm::Ptr<reco::Muon>> candidates1, candidates2;

  fillByCategory(highptMuons1,candidates1,promptMuons,promptMuons.begin(),muonHandle,primaryVertex);
  fillByCategory(highptMuons2,candidates2,promptMuons,std::next(promptMuons.begin(),2),muonHandle,primaryVertex);

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

  int nCand1 = pushMuon(finalMuons,highptMuons1, candidates1);
  int nCand2 = pushMuon(finalMuons,highptMuons2, candidates2);
  int nmuon = nCand1 + nCand2;

  if ( nmuon < 4 ) {
    edm::Ptr<reco::Muon> mergedMuon;

    if ( nmuon < 3 )
      mergedMuon = ( highptMuons1.front()->tunePMuonBestTrack()->pt() < highptMuons2.front()->tunePMuonBestTrack()->pt() )
                   ? highptMuons1.front() : highptMuons2.front();
    else
      mergedMuon = ( nCand1 < 2 ) ? highptMuons1.front() : highptMuons2.front();

    float tunepPt = mergedMuon->tunePMuonBestTrack()->pt();
    auto matched = ( mergedMuon==highptMuons1.front() ) ? promptMuons.begin() : std::next(promptMuons.begin(),2);
    auto probe = std::next(matched,1);

    if ( std::abs(tunepPt-(*matched)->pt()) > std::abs(tunepPt-(*probe)->pt()) ) {
      probe = matched;
      matched = std::next(matched,1);
    }

    fillMets(mergedMuon, *probe, pfMEThandle->at(0), finalMuons, "pfMET");
    fillMets(mergedMuon, *probe, puppiMEThandle->at(0), finalMuons, "puppiMET");
  }
}

int MergedLeptonAnalyzer::categorize(std::vector<edm::Ptr<reco::Muon>>& muons, const std::vector<edm::Ptr<reco::GenParticle>>& promptMuons) {
  if ( muons.size() > 1 )
    std::sort(muons.begin(),muons.end(),sortByTuneP);
  else if ( muons.size() == 1 ) {
    int nMatched = 0;
    const auto highptmu = muons.front();

    for (const auto genptc : promptMuons) {
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
  const std::vector<edm::Ptr<reco::GenParticle>>& promptMuons,
  std::vector<edm::Ptr<reco::GenParticle>>::const_iterator muonItr,
  const edm::Handle<edm::View<reco::Muon>>& muonHandle,
  const reco::Vertex& vtx) {
  switch (categorize(muons,promptMuons)) {
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
        fillMuons(aMuon,*matched,vtx,"solo");
        return;
      }

      std::sort(candidates.begin(),candidates.end(),sortByTuneP);
      fillMuons(aMuon,*matched,vtx,"tag");
      fillMuons(candidates.front(),*probe,vtx,"probe");

      return;
    } default:
      fillMuons(muons.front(),*muonItr,vtx,"tag");
      fillMuons(muons.at(1),*std::next(muonItr,1),vtx,"probe");
      break;
  }

  return;
}

bool MergedLeptonAnalyzer::sortByTuneP(const edm::Ptr<reco::Muon> a, const edm::Ptr<reco::Muon> b) {
  return a->tunePMuonBestTrack()->pt() > b->tunePMuonBestTrack()->pt();
}

void MergedLeptonAnalyzer::fillMuons(const edm::Ptr<reco::Muon>& aMuon,
                                     const edm::Ptr<reco::GenParticle>& matched,
                                     const reco::Vertex& vtx,
                                     std::string prefix) {
  if (aMuon->isGlobalMuon())
    histo1d_[prefix+"_muType"]->Fill(0.5);
  if (aMuon->isTrackerMuon())
    histo1d_[prefix+"_muType"]->Fill(1.5);
  if (aMuon->isStandAloneMuon())
    histo1d_[prefix+"_muType"]->Fill(2.5);

  if (muon::isHighPtMuon(*aMuon,vtx))
    histo1d_[prefix+"_passIDs"]->Fill(0.5);
  if (isHighPtTrackerMuon(*aMuon,vtx))
    histo1d_[prefix+"_passIDs"]->Fill(1.5);
  if (muon::isLooseMuon(*aMuon))
    histo1d_[prefix+"_passIDs"]->Fill(2.5);
  if (muon::isMediumMuon(*aMuon))
    histo1d_[prefix+"_passIDs"]->Fill(3.5);
  if (muon::isTightMuon(*aMuon,vtx))
    histo1d_[prefix+"_passIDs"]->Fill(4.5);

  if (!aMuon->isTrackerMuon())
    return;

  auto innerTrack = aMuon->innerTrack();
  values_[prefix+"_muon"].numberOfValidTrackerHits = innerTrack->hitPattern().numberOfValidTrackerHits();
  values_[prefix+"_muon"].numberOfValidPixelHits = innerTrack->hitPattern().numberOfValidPixelHits();
  values_[prefix+"_muon"].numberOfValidStripHits = innerTrack->hitPattern().numberOfValidStripHits();
  values_[prefix+"_muon"].trackerLayersWithMeasurement = innerTrack->hitPattern().trackerLayersWithMeasurement();
  values_[prefix+"_muon"].pixelLayersWithMeasurement = innerTrack->hitPattern().pixelLayersWithMeasurement();
  values_[prefix+"_muon"].stripLayersWithMeasurement = innerTrack->hitPattern().stripLayersWithMeasurement();
  values_[prefix+"_muon"].trackerLayersWithoutMeasurement = innerTrack->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  values_[prefix+"_muon"].pixelLayersWithoutMeasurement = innerTrack->hitPattern().pixelLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  values_[prefix+"_muon"].stripLayersWithoutMeasurement = innerTrack->hitPattern().stripLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  values_[prefix+"_muon"].trackerVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidTrackerHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().trackerLayersWithMeasurement());
  values_[prefix+"_muon"].pixelVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidPixelHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().pixelLayersWithMeasurement());
  values_[prefix+"_muon"].stripVoM = static_cast<float>(innerTrack->hitPattern().numberOfValidStripHits())/static_cast<float>(aMuon->innerTrack()->hitPattern().stripLayersWithMeasurement());
  values_[prefix+"_muon"].pfPt = aMuon->pt();
  values_[prefix+"_muon"].tunepPt = aMuon->tunePMuonBestTrack()->pt();
  values_[prefix+"_muon"].genPt = matched->pt();
  values_[prefix+"_muon"].PFoGen = aMuon->pt()/matched->pt();
  values_[prefix+"_muon"].TPoGen = aMuon->tunePMuonBestTrack()->pt()/matched->pt();
  tree_[prefix+"_muonTree"]->Fill();

  histo1d_[prefix+"_algo"]->Fill(static_cast<float>(innerTrack->algo())+0.5);
  histo1d_[prefix+"_orgAlgo"]->Fill(static_cast<float>(innerTrack->originalAlgo())+0.5);

  for (int en = reco::TrackBase::undefAlgorithm; en != reco::TrackBase::algoSize; ++en) {
    if ( (innerTrack->algoMask())[en] )
      histo1d_[prefix+"_algoMask"]->Fill(static_cast<float>(en)+0.5);
  }

  for (int en = reco::TrackBase::loose; en != reco::TrackBase::qualitySize; ++en) {
    if ( innerTrack->quality( static_cast<reco::TrackBase::TrackQuality>(en) ) )
      histo1d_[prefix+"_quality"]->Fill(static_cast<float>(en)+0.5);
  }
}

void MergedLeptonAnalyzer::fillMets(const edm::Ptr<reco::Muon>& aMuon,
                                    const edm::Ptr<reco::GenParticle>& matched,
                                    const pat::MET& met,
                                    const std::vector<edm::Ptr<reco::Muon>>& finalMuons,
                                    std::string prefix) {
  metvalues_[prefix].PFphi = met.phi();
  metvalues_[prefix].PFdPhi = reco::deltaPhi(met.phi(),matched->phi());
  metvalues_[prefix].PFpt = met.pt();
  metvalues_[prefix].PFoGen = met.pt()/matched->pt();
  metvalues_[prefix].PFSumEt = met.sumEt();

  pat::MET::Vector2 tpvec;
  tpvec.px = met.px() + aMuon->px() - aMuon->tunePMuonBestTrack()->px();
  tpvec.py = met.py() + aMuon->py() - aMuon->tunePMuonBestTrack()->py();
  metvalues_[prefix].TPphi = tpvec.phi();
  metvalues_[prefix].TPdPhi = reco::deltaPhi(tpvec.phi(),matched->phi());
  metvalues_[prefix].TPpt = tpvec.pt();
  metvalues_[prefix].TPoGen = tpvec.pt()/matched->pt();

  float sumpt = 0.;
  float sumTPpt = 0.;

  for (auto mu : finalMuons) {
    sumpt += mu->pt();
    sumTPpt += mu->tunePMuonBestTrack()->pt();
  }

  metvalues_[prefix].TPSumEt = met.sumEt() + sumTPpt - sumpt;
  metvalues_[prefix].dSumEt = met.sumEt() - sumpt;
  metvalues_[prefix].nMuon = static_cast<int>(finalMuons.size());

  tree_[prefix+"Tree"]->Fill();
}

bool MergedLeptonAnalyzer::isHighPtTrackerMuon(const reco::Muon& muon, const reco::Vertex& vtx) {
  if (!muon.isTrackerMuon())
    return false;

  bool muMatchedSt = muon.numberOfMatchedStations() > 1;
  if (!muMatchedSt) {
    if (muon.isTrackerMuon() && muon.numberOfMatchedStations() == 1) {
      if (muon.expectedNnumberOfMatchedStations() < 2 || !(muon.stationMask() == 1 || muon.stationMask() == 16) ||
          muon.numberOfMatchedRPCLayers() > 2)
        muMatchedSt = true;
    }
  }

  bool muID = muMatchedSt;

  bool hits = muon.innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5 &&
              muon.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;

  bool momQuality = muon.tunePMuonBestTrack()->ptError() / muon.tunePMuonBestTrack()->pt() < 0.3;

  bool ip =
      std::abs(muon.innerTrack()->dxy(vtx.position())) < 0.2 && std::abs(muon.innerTrack()->dz(vtx.position())) < 0.5;

  return muID && hits && momQuality && ip;
}

DEFINE_FWK_MODULE(MergedLeptonAnalyzer);
