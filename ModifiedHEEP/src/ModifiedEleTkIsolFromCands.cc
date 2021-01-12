#include "ZprimeTo4l/ModifiedHEEP/interface/ModifiedEleTkIsolFromCands.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

ModifiedEleTkIsolFromCands::TrkCuts::TrkCuts(const edm::ParameterSet& para)
{
  minPt = para.getParameter<double>("minPt");
  auto sq = [](double val){return val*val;};
  minDR2 = sq(para.getParameter<double>("minDR"));
  maxDR2 = sq(para.getParameter<double>("maxDR"));
  minDEta = para.getParameter<double>("minDEta");
  dEta2nd = para.getParameter<double>("dEta2nd");
  dPhi2nd = para.getParameter<double>("dPhi2nd");
  maxDZ = para.getParameter<double>("maxDZ");
  minHits = para.getParameter<int>("minHits");
  minPixelHits = para.getParameter<int>("minPixelHits");
  maxDPtPt = para.getParameter<double>("maxDPtPt");
  addGsfminPt = para.getParameter<double>("addGsfminPt");
  addGsfDR2 = sq(para.getParameter<double>("maxDR"));

  auto qualNames = para.getParameter<std::vector<std::string> >("allowedQualities");
  auto algoNames = para.getParameter<std::vector<std::string> >("algosToReject");

  for(auto& qualName : qualNames){
    allowedQualities.push_back(reco::TrackBase::qualityByName(qualName));
  }
  for(auto& algoName : algoNames){
    algosToReject.push_back(reco::TrackBase::algoByName(algoName));
  }
  std::sort(algosToReject.begin(),algosToReject.end());

}

edm::ParameterSetDescription ModifiedEleTkIsolFromCands::TrkCuts::pSetDescript()
{
  edm::ParameterSetDescription desc;
  desc.add<double>("minPt",1.0);
  desc.add<double>("maxDR",0.3);
  desc.add<double>("minDR",0.000);
  desc.add<double>("minDEta",0.005);
  desc.add<double>("dEta2nd",0.005);
  desc.add<double>("dPhi2nd",0.05);
  desc.add<double>("maxDZ",0.1);
  desc.add<double>("maxDPtPt",-1);
  desc.add<double>("addGsfminPt",10.0);
  desc.add<int>("minHits",0);
  desc.add<int>("minPixelHits",0);
  desc.add<std::vector<std::string> >("allowedQualities",{"highPurity"});
  desc.add<std::vector<std::string> >("algosToReject");
  return desc;
}

ModifiedEleTkIsolFromCands::ModifiedEleTkIsolFromCands(const edm::ParameterSet& para):
  barrelCuts_(para.getParameter<edm::ParameterSet>("barrelCuts")),
  endcapCuts_(para.getParameter<edm::ParameterSet>("endcapCuts"))
{


}

edm::ParameterSetDescription ModifiedEleTkIsolFromCands::pSetDescript()
{
  edm::ParameterSetDescription desc;
  desc.add("barrelCuts",TrkCuts::pSetDescript());
  desc.add("endcapCuts",TrkCuts::pSetDescript());
  return desc;
}

std::pair<int,double>
ModifiedEleTkIsolFromCands::calIsol(const reco::TrackBase& eleTrk,
			    const pat::PackedCandidateCollection& cands,
          const reco::TrackBase& addTrk,
          int& nrMatchedTrk,
          float& rtMatchedTrk,
			    const PIDVeto pidVeto)const
{
  return calIsol(eleTrk.eta(),eleTrk.phi(),eleTrk.vz(),cands,addTrk,nrMatchedTrk,rtMatchedTrk,pidVeto);
}

std::pair<int,double>
ModifiedEleTkIsolFromCands::calIsol(const double eleEta,const double elePhi,
			    const double eleVZ,
			    const pat::PackedCandidateCollection& cands,
          const reco::TrackBase& addTrk,
          int& nrMatchedTrk,
          float& rtMatchedTrk,
			    const PIDVeto pidVeto)const
{
  double ptSum=0.;
  int nrTrks=0;

  const TrkCuts& cuts = std::abs(eleEta)<1.5 ? barrelCuts_ : endcapCuts_;

  for(auto& cand  : cands) {
    if(cand.hasTrackDetails() && cand.charge()!=0 && passPIDVeto(cand.pdgId(),pidVeto)){ // 94X
    // if(cand.charge()!=0 && passPIDVeto(cand.pdgId(),pidVeto)){ // 80X
      const reco::Track& trk = cand.pseudoTrack();
      if(passTrkSel(trk,trk.pt(),cuts,eleEta,elePhi,eleVZ)){
        if(std::abs(addTrk.eta()-trk.eta()) < cuts.dEta2nd && std::abs( reco::deltaPhi(addTrk.phi(),trk.phi()) ) < cuts.dPhi2nd) {
          nrMatchedTrk++;
          rtMatchedTrk = addTrk.pt()/trk.pt();
          continue;
        }

      	ptSum+=trk.pt();
      	nrTrks++;
      }
    }
  }

  return {nrTrks,ptSum};
}



std::pair<int,double>
ModifiedEleTkIsolFromCands::calIsol(const reco::TrackBase& eleTrk,
			    const reco::TrackCollection& tracks, const reco::TrackBase& addTrk, int& nrMatchedTrk, float& rtMatchedTrk)const
{
  return calIsol(eleTrk.eta(),eleTrk.phi(),eleTrk.vz(),tracks,addTrk,nrMatchedTrk,rtMatchedTrk);
}

std::pair<int,double>
ModifiedEleTkIsolFromCands::calIsol(const double eleEta,const double elePhi,
			    const double eleVZ,
			    const reco::TrackCollection& tracks,
          const reco::TrackBase& addTrk,
          int& nrMatchedTrk, float& rtMatchedTrk)const
{
  double ptSum=0.;
  int nrTrks=0;

  const TrkCuts& cuts = std::abs(eleEta)<1.5 ? barrelCuts_ : endcapCuts_;

  for(auto& trk : tracks){
    if(passTrkSel(trk,trk.pt(),cuts,eleEta,elePhi,eleVZ)){
      if(std::abs(addTrk.eta()-trk.eta()) < cuts.dEta2nd && std::abs( reco::deltaPhi(addTrk.phi(),trk.phi()) ) < cuts.dPhi2nd) {
        nrMatchedTrk++;
        rtMatchedTrk = addTrk.pt()/trk.pt();
        continue;
      }

      ptSum+=trk.pt();
      nrTrks++;
    }
  }
  return {nrTrks,ptSum};
}

bool ModifiedEleTkIsolFromCands::passPIDVeto(const int pdgId,const ModifiedEleTkIsolFromCands::PIDVeto veto)
{
  int pidAbs = std::abs(pdgId);
  switch (veto){
  case PIDVeto::NONE:
    return true;
  case PIDVeto::ELES:
    if(pidAbs==11) return false;
    else return true;
  case PIDVeto::NONELES:
    if(pidAbs==11) return true;
    else return false;
  }
  throw cms::Exception("CodeError") <<
    "invalid PIDVeto "<<static_cast<int>(veto)<<", "<<
    "this is likely due to some static casting of invalid ints somewhere";
}

ModifiedEleTkIsolFromCands::PIDVeto ModifiedEleTkIsolFromCands::pidVetoFromStr(const std::string& vetoStr)
{
  if(vetoStr=="NONE") return PIDVeto::NONE;
  else if(vetoStr=="ELES") return PIDVeto::ELES;
  else if(vetoStr=="NONELES") return PIDVeto::NONELES;
  else{
    throw cms::Exception("CodeError") <<"unrecognised string "<<vetoStr<<", either a typo or this function needs to be updated";
  }
}

bool ModifiedEleTkIsolFromCands::passTrkSel(const reco::TrackBase& trk,
				    const double trkPt,const TrkCuts& cuts,
				    const double eleEta,const double elePhi,
				    const double eleVZ)
{
  const float dR2 = reco::deltaR2(eleEta,elePhi,trk.eta(),trk.phi());
  const float dEta = trk.eta()-eleEta;
  const float dZ = eleVZ - trk.vz();

  return dR2>=cuts.minDR2 && dR2<=cuts.maxDR2 &&
    std::abs(dEta)>=cuts.minDEta &&
    std::abs(dZ)<cuts.maxDZ &&
    trk.hitPattern().numberOfValidHits() >= cuts.minHits &&
    trk.hitPattern().numberOfValidPixelHits() >=cuts.minPixelHits &&
    (trk.ptError()/trkPt < cuts.maxDPtPt || cuts.maxDPtPt<0) &&
    passQual(trk,cuts.allowedQualities) &&
    passAlgo(trk,cuts.algosToReject) &&
    trkPt > cuts.minPt;
}

bool ModifiedEleTkIsolFromCands::additionalGsfTrkSel(const reco::TrackBase& trk,
				    const double trkPt,const TrkCuts& cuts,
				    const double eleEta,const double elePhi,
				    const double eleVZ)
{
  const float dR2 = reco::deltaR2(eleEta,elePhi,trk.eta(),trk.phi());
  // const float dEta = trk.eta()-eleEta;
  const float dZ = eleVZ - trk.vz();

  return trk.eta() != eleEta &&
    dR2 <= cuts.addGsfDR2 &&
  //  std::abs(dEta)>=cuts.minDEta &&
    std::abs(dZ)<cuts.maxDZ &&
    trk.hitPattern().numberOfValidHits() >= cuts.minHits &&
    trk.hitPattern().numberOfValidPixelHits() >= cuts.minPixelHits &&
    //(trk.ptError()/trkPt < cuts.maxDPtPt || cuts.maxDPtPt<0) &&
    passQual(trk,cuts.allowedQualities) &&
    passAlgo(trk,cuts.algosToReject) &&
    trkPt > cuts.addGsfminPt;
}

const reco::GsfTrackRef ModifiedEleTkIsolFromCands::additionalGsfTrkSelector(const reco::GsfElectron& ele, edm::Handle<reco::GsfTrackCollection> gsfTrks, bool& addGsfTrkSel) {
  int noSelectedGsfTrk = 0;
  return additionalGsfTrkSelector(ele, gsfTrks, addGsfTrkSel, noSelectedGsfTrk);
}

const reco::GsfTrackRef ModifiedEleTkIsolFromCands::additionalGsfTrkSelector(const reco::GsfElectron& ele, edm::Handle<reco::GsfTrackCollection> gsfTrks, bool& addGsfTrkSel, int& noSelectedGsfTrk) {
  std::vector<reco::GsfTrackRef> additionalGsfTrks;
  const reco::GsfTrackRef eleTrk = ele.gsfTrack();
  for (size_t gsfTrkNr = 0; gsfTrkNr < gsfTrks->size(); gsfTrkNr++) {
    auto gsfTrk = gsfTrks->at(gsfTrkNr);
    if( gsfTrk.pt()==eleTrk->pt() && gsfTrk.eta()==eleTrk->eta() && gsfTrk.phi()==eleTrk->phi() ) continue;
    const ModifiedEleTkIsolFromCands::TrkCuts& cuts = std::abs(gsfTrk.eta())<1.5 ? barrelCuts_ : endcapCuts_;
    if(additionalGsfTrkSel(gsfTrk,gsfTrk.pt(),cuts,eleTrk->eta(),eleTrk->phi(),eleTrk->vz())) {
      additionalGsfTrks.push_back(edm::Ref<reco::GsfTrackCollection>(gsfTrks,gsfTrkNr));
      addGsfTrkSel = true;
    }
  }
  std::sort(additionalGsfTrks.begin(), additionalGsfTrks.end(), [](reco::GsfTrackRef a, reco::GsfTrackRef b) {return a->pt() > b->pt();} );
  if(additionalGsfTrks.empty()) {
    addGsfTrkSel = false;
    noSelectedGsfTrk = 0;
    return eleTrk;
  }

  noSelectedGsfTrk = additionalGsfTrks.size();

  return additionalGsfTrks.front();
}

bool ModifiedEleTkIsolFromCands::
passQual(const reco::TrackBase& trk,
	 const std::vector<reco::TrackBase::TrackQuality>& quals)
{
  if(quals.empty()) return true;

  for(auto qual : quals) {
    if(trk.quality(qual)) return true;
  }

  return false;
}

bool ModifiedEleTkIsolFromCands::
passAlgo(const reco::TrackBase& trk,
	 const std::vector<reco::TrackBase::TrackAlgorithm>& algosToRej)
{
  return algosToRej.empty() || !std::binary_search(algosToRej.begin(),algosToRej.end(),trk.algo());
}
