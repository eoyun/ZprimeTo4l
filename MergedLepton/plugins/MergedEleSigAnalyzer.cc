#include <memory>
#include <iostream>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "TH1D.h"
#include "TTree.h"
#include "TString.h"

// produce TTree for merged electron training with H->AA->4e events

class MergedEleSigAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedEleSigAnalyzer(const edm::ParameterSet&);
  virtual ~MergedEleSigAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override {}

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<reco::ConversionCollection> conversionsToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;

  const double ptThres_;
  const double ptThres2nd_;
  const double drThres_;

  std::map<std::string,TH1*> histo1d_;
};

MergedEleSigAnalyzer::MergedEleSigAnalyzer(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
conversionsToken_(consumes<reco::ConversionCollection>(iConfig.getParameter<edm::InputTag>("conversions"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
ptThres_(iConfig.getParameter<double>("ptThres")),
ptThres2nd_(iConfig.getParameter<double>("ptThres2nd")),
drThres_(iConfig.getParameter<double>("drThres")) {
  usesResource("TFileService");
}

void MergedEleSigAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;
  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);

  auto addHists = [&] (std::string postfix) {
    histo1d_["id_dEtaInSeed_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_merged_noGsf_"+postfix).c_str(),"Merged w/o GSF #Delta#eta_{in}(seed)",200,-0.025,0.025);
    histo1d_["id_dEtaInSeed_merged_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_merged_"+postfix).c_str(),"Merged w/ GSF #Delta#eta_{in}(seed)",200,-0.025,0.025);
    histo1d_["id_dEtaInSeed_resolved_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_resolved_"+postfix).c_str(),"Resolved #Delta#eta_{in}(seed)",200,-0.025,0.025);

    histo1d_["id_dR_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_dR_merged_noGsf_"+postfix).c_str(),"Merged w/o GSF GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_merged_"+postfix] = fs->make<TH1D>(("id_dR_merged_"+postfix).c_str(),"Merged w/ GSF GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_resolved_"+postfix] = fs->make<TH1D>(("id_dR_resolved_"+postfix).c_str(),"Resolved GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_dR_has2ndTrk_"+postfix).c_str(),"w/ 2nd trk GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_all_"+postfix] = fs->make<TH1D>(("id_dR_all_"+postfix).c_str(),"GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_sum_"+postfix] = fs->make<TH1D>(("id_dR_sum_"+postfix).c_str(),"GEN-lv #Delta R",1000,0.,1.);

    histo1d_["id_dR_GSF_merged_"+postfix] = fs->make<TH1D>(("id_dR_GSF_merged_"+postfix).c_str(),"Merged w/ GSF RECO-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_GSF_resolved_"+postfix] = fs->make<TH1D>(("id_dR_GSF_resolved_"+postfix).c_str(),"Resolved RECO-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_GSF_sum_"+postfix] = fs->make<TH1D>(("id_dR_GSF_sum_"+postfix).c_str(),"RECO-lv #Delta R",1000,0.,1.);

    histo1d_["id_showershape_"+postfix] = fs->make<TH1D>(("id_showershape_"+postfix).c_str(),"GEN-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_trackIsoHeep_"+postfix] = fs->make<TH1D>(("id_trackIsoHeep_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_trackIsoModified_"+postfix] = fs->make<TH1D>(("id_trackIsoModified_"+postfix).c_str(),"GEN-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_caloIsoHeep_"+postfix] = fs->make<TH1D>(("id_caloIsoHeep_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_caloIsoModified_"+postfix] = fs->make<TH1D>(("id_caloIsoModified_"+postfix).c_str(),"GEN-lv #Delta R passing modified calo iso",1000,0.,1.);

    histo1d_["id_conversionVeto_"+postfix] = fs->make<TH1D>(("id_conversionVeto_"+postfix).c_str(),"GEN-lv #Delta R passing conversion veto",1000,0.,1.);
    histo1d_["id_passHEEP_"+postfix] = fs->make<TH1D>(("id_passHEEP_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_passModifiedHEEP_"+postfix] = fs->make<TH1D>(("id_passModifiedHEEP_"+postfix).c_str(),"GEN-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_showershape_merged_"+postfix] = fs->make<TH1D>(("id_showershape_merged_"+postfix).c_str(),"GEN-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_trackIsoHeep_merged_"+postfix] = fs->make<TH1D>(("id_trackIsoHeep_merged_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_trackIsoModified_merged_"+postfix] = fs->make<TH1D>(("id_trackIsoModified_merged_"+postfix).c_str(),"GEN-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_caloIsoHeep_merged_"+postfix] = fs->make<TH1D>(("id_caloIsoHeep_merged_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_caloIsoModified_merged_"+postfix] = fs->make<TH1D>(("id_caloIsoModified_merged_"+postfix).c_str(),"GEN-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_passHEEP_merged_"+postfix] = fs->make<TH1D>(("id_passHEEP_merged_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_passModifiedHEEP_merged_"+postfix] = fs->make<TH1D>(("id_passModifiedHEEP_merged_"+postfix).c_str(),"GEN-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_showershape_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_showershape_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_trackIsoHeep_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_trackIsoHeep_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_trackIsoModified_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_trackIsoModified_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_caloIsoHeep_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_caloIsoHeep_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_caloIsoModified_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_caloIsoModified_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_passHEEP_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_passHEEP_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_passModifiedHEEP_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_passModifiedHEEP_merged_noGsf_"+postfix).c_str(),"GEN-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_showershape_resolved_"+postfix] = fs->make<TH1D>(("id_showershape_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_trackIsoHeep_resolved_"+postfix] = fs->make<TH1D>(("id_trackIsoHeep_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_trackIsoModified_resolved_"+postfix] = fs->make<TH1D>(("id_trackIsoModified_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_caloIsoHeep_resolved_"+postfix] = fs->make<TH1D>(("id_caloIsoHeep_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_caloIsoModified_resolved_"+postfix] = fs->make<TH1D>(("id_caloIsoModified_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_passHEEP_resolved_"+postfix] = fs->make<TH1D>(("id_passHEEP_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_passModifiedHEEP_resolved_"+postfix] = fs->make<TH1D>(("id_passModifiedHEEP_resolved_"+postfix).c_str(),"GEN-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_showershape_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_showershape_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_trackIsoHeep_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_trackIsoHeep_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_trackIsoModified_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_trackIsoModified_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_caloIsoHeep_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_caloIsoHeep_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_caloIsoModified_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_caloIsoModified_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_passHEEP_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_passHEEP_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_passModifiedHEEP_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_passModifiedHEEP_has2ndTrk_"+postfix).c_str(),"GEN-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_GSF_showershape_"+postfix] = fs->make<TH1D>(("id_GSF_showershape_"+postfix).c_str(),"RECO-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_GSF_trackIsoHeep_"+postfix] = fs->make<TH1D>(("id_GSF_trackIsoHeep_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_GSF_trackIsoModified_"+postfix] = fs->make<TH1D>(("id_GSF_trackIsoModified_"+postfix).c_str(),"RECO-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_GSF_caloIsoHeep_"+postfix] = fs->make<TH1D>(("id_GSF_caloIsoHeep_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_GSF_caloIsoModified_"+postfix] = fs->make<TH1D>(("id_GSF_caloIsoModified_"+postfix).c_str(),"RECO-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_GSF_passHEEP_"+postfix] = fs->make<TH1D>(("id_GSF_passHEEP_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_GSF_passModifiedHEEP_"+postfix] = fs->make<TH1D>(("id_GSF_passModifiedHEEP_"+postfix).c_str(),"RECO-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_GSF_showershape_merged_"+postfix] = fs->make<TH1D>(("id_GSF_showershape_merged_"+postfix).c_str(),"RECO-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_GSF_trackIsoHeep_merged_"+postfix] = fs->make<TH1D>(("id_GSF_trackIsoHeep_merged_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_GSF_trackIsoModified_merged_"+postfix] = fs->make<TH1D>(("id_GSF_trackIsoModified_merged_"+postfix).c_str(),"RECO-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_GSF_caloIsoHeep_merged_"+postfix] = fs->make<TH1D>(("id_GSF_caloIsoHeep_merged_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_GSF_caloIsoModified_merged_"+postfix] = fs->make<TH1D>(("id_GSF_caloIsoModified_merged_"+postfix).c_str(),"RECO-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_GSF_passHEEP_merged_"+postfix] = fs->make<TH1D>(("id_GSF_passHEEP_merged_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_GSF_passModifiedHEEP_merged_"+postfix] = fs->make<TH1D>(("id_GSF_passModifiedHEEP_merged_"+postfix).c_str(),"RECO-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["id_GSF_showershape_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_showershape_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing shower shape cuts",1000,0.,1.);
    histo1d_["id_GSF_trackIsoHeep_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_trackIsoHeep_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP track iso",1000,0.,1.);
    histo1d_["id_GSF_trackIsoModified_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_trackIsoModified_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing modified track iso",1000,0.,1.);
    histo1d_["id_GSF_caloIsoHeep_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_caloIsoHeep_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP calo iso",1000,0.,1.);
    histo1d_["id_GSF_caloIsoModified_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_caloIsoModified_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing modified calo iso",1000,0.,1.);
    histo1d_["id_GSF_passHEEP_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_passHEEP_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing HEEP",1000,0.,1.);
    histo1d_["id_GSF_passModifiedHEEP_resolved_"+postfix] = fs->make<TH1D>(("id_GSF_passModifiedHEEP_resolved_"+postfix).c_str(),"RECO-lv #Delta R passing modified HEEP",1000,0.,1.);

    histo1d_["GSF_pt_1st_"+postfix] = fs->make<TH1D>(("GSF_pt_1st_"+postfix).c_str(),"1st GSF track Pt of merged electron",500,0.,1000.);
    histo1d_["GSF_pt_2nd_"+postfix] = fs->make<TH1D>(("GSF_pt_2nd_"+postfix).c_str(),"2nd GSF track Pt of merged electron",500,0.,1000.);
    histo1d_["GSF_dPtOverPt_"+postfix] = fs->make<TH1D>(("GSF_dPtOverPt_"+postfix).c_str(),"GSF track dPt/Pt",500,0.,5.);
  };

  histo1d_["GEN_pt_1st"] = fs->make<TH1D>("GEN_pt_1st","GEN-lv Pt of leading lepton",500,0.,1000.);
  histo1d_["GEN_pt_2nd"] = fs->make<TH1D>("GEN_pt_2nd","GEN-lv Pt of subleading lepton",500,0.,1000.);

  addHists("all");
  addHists("EB-EB");
  addHists("EB-EE");
  addHists("EE-EE");
}

void MergedEleSigAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
  iEvent.getByToken(srcGenPtc_, genptcHandle);

  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken_, genInfo);
  double mcweight = genInfo->weight();

  double aWeight = mcweight/std::abs(mcweight);
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  edm::Ptr<reco::Vertex> primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->ptrAt(0);
  else
    return;

  std::vector<reco::GenParticleRef> promptLeptons;
  double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    const auto& genptc = genptcHandle->refAt(idx);

    if ( ( std::abs(genptc->pdgId())==11 ) &&
         genptc->fromHardProcessFinalState() &&
         ( std::abs(genptc->eta()) < 2.5 ) )
      promptLeptons.push_back(genptc.castTo<reco::GenParticleRef>());
  }

  if (promptLeptons.size()!=4)
    return;

  // sort by mother ptc, and then by et (motherA Et1, motherA Et2, motherB Et1, motherB Et2)
  auto sisterLambda = [](const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
    return (a->mother() == b->mother()) ? (a->et() > b->et()) : (a->mother() < b->mother());
  };

  std::sort(promptLeptons.begin(),promptLeptons.end(),sisterLambda);

  std::vector<pat::ElectronRef> pair1, pair2;

  for (unsigned int idx = 0; idx < eleHandle->size(); ++idx) {
    const auto& aEle = eleHandle->refAt(idx);

    if ( aEle->et() < ptThres2nd_ )
      continue;

    bool matched = false;
    size_t igen = 0;

    // for given reco electron, check whether GEN-lv prompt electron exists within dR thres
    for (; igen < promptLeptons.size(); ++igen) {
      const auto& genptc = promptLeptons.at(igen);
      double dr2 = reco::deltaR2(aEle->eta(),aEle->phi(),genptc->eta(),genptc->phi());

      if (dr2 < drThres2) {
        matched = true;
        break;
      }
    }

    if ( !matched )
      continue;

    auto castEle = aEle.castTo<pat::ElectronRef>();

    (igen < 2) ? pair1.push_back(castEle) : pair2.push_back(castEle);
  }

  // assume 4e only
  if ( pair1.size() == 0 || pair2.size() == 0 )
    return;

  const double drPrompt1 = std::sqrt( reco::deltaR2(promptLeptons.at(0)->eta(),promptLeptons.at(0)->phi(),
                                                    promptLeptons.at(1)->eta(),promptLeptons.at(1)->phi()) );
  const double drPrompt2 = std::sqrt( reco::deltaR2(promptLeptons.at(2)->eta(),promptLeptons.at(2)->phi(),
                                                    promptLeptons.at(3)->eta(),promptLeptons.at(3)->phi()) );

  histo1d_["GEN_pt_1st"]->Fill( promptLeptons.at(0)->pt(),aWeight );
  histo1d_["GEN_pt_1st"]->Fill( promptLeptons.at(2)->pt(),aWeight );
  histo1d_["GEN_pt_2nd"]->Fill( promptLeptons.at(1)->pt(),aWeight );
  histo1d_["GEN_pt_2nd"]->Fill( promptLeptons.at(3)->pt(),aWeight );

  if ( promptLeptons.at(0)->pt() < ptThres_ || promptLeptons.at(2)->pt() < ptThres_ )
    return;

  if ( promptLeptons.at(1)->pt() < ptThres2nd_ || promptLeptons.at(3)->pt() < ptThres2nd_ )
    return;

  if ( pair1.size() > 1 && promptLeptons.at(1)->pt() < ptThres_ )
    return;

  if ( pair2.size() > 1 && promptLeptons.at(3)->pt() < ptThres_ )
    return;

  auto passVIDbit = [] (const int32_t bitmap, unsigned cutIdx) -> bool {
    int32_t mask = 0x00000001 << cutIdx;
    int32_t result = mask & bitmap;

    return result!=0x00000000;
  };

  auto fillElectrons = [&] (std::vector<pat::ElectronRef>& apair, const double drPrompt, std::string postfix) {
    histo1d_["id_dR_all_"+postfix]->Fill(drPrompt,aWeight);

    if ( apair.size()==1 ) {
      auto addGsfTrk = (*addGsfTrkHandle)[apair.front()];
      auto orgGsfTrk = apair.front()->gsfTrack();

      if ( postfix=="EB-EB" ) {
        if ( std::abs(apair.at(0)->superCluster()->eta()) > 1.5 || std::abs(addGsfTrk->eta()) > 1.5 )
          return;
      }

      if ( postfix=="EB-EE" ) {
        if ( ( std::abs(apair.at(0)->superCluster()->eta()) > 1.5 && std::abs(addGsfTrk->eta()) > 1.5 ) ||
             ( std::abs(apair.at(0)->superCluster()->eta()) < 1.5 && std::abs(addGsfTrk->eta()) < 1.5 ) )
          return;
      }

      if ( postfix=="EE-EE" ) {
        if ( std::abs(apair.at(0)->superCluster()->eta()) < 1.5 || std::abs(addGsfTrk->eta()) < 1.5 )
          return;
      }

      if ( addGsfTrk==orgGsfTrk ) { // ME w/o add GSF
        // veto any nearby electron over Et thres
        for (unsigned jdx = 0; jdx < eleHandle->size(); jdx++) {
          const auto& secEleRef = eleHandle->refAt(jdx);
          auto secEle = secEleRef.castTo<pat::ElectronRef>();

          if (secEle==apair.front())
            continue;

          if (secEle->et() < ptThres2nd_)
            continue;

          double dr2 = reco::deltaR2(apair.front()->eta(),apair.front()->phi(),secEle->eta(),secEle->phi());

          if (dr2 < drThres_*drThres_)
            return;
        }

        // fill merged2 histo here
        histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dEtaInSeed_merged_noGsf_"+postfix]->Fill(apair.front()->deltaEtaSeedClusterTrackAtVtx(),aWeight);

        if ( passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),4) &&
             passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),5) ) {
          histo1d_["id_showershape_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_showershape_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),7) ) {
          histo1d_["id_trackIsoHeep_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_trackIsoHeep_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),8) ) {
          histo1d_["id_caloIsoHeep_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_caloIsoHeep_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("modifiedHeepElectronID"),7) ) {
          histo1d_["id_trackIsoModified_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_trackIsoModified_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("modifiedHeepElectronID"),8) ) {
          histo1d_["id_caloIsoModified_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_caloIsoModified_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }

        if ( apair.front()->passConversionVeto() )
          histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

        if ( apair.front()->electronID("heepElectronID-HEEPV70") ) {
          histo1d_["id_passHEEP_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_passHEEP_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }

        if ( apair.front()->electronID("modifiedHeepElectronID") ) {
          histo1d_["id_passModifiedHEEP_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_passModifiedHEEP_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        }
      } else { // ME w/ GSF
        // find whether add GSF track has a corresponding electron
        bool notMerged = false;

        for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
          const auto& secEle = eleHandle->refAt(jdx);
          const auto& secGsfTrk = secEle->gsfTrack();

          if ( addGsfTrk==secGsfTrk ) {
            notMerged = true;
            break;
          }
        }

        if ( notMerged )
          return;

        const double drReco = std::sqrt(reco::deltaR2(addGsfTrk->eta(),addGsfTrk->phi(),orgGsfTrk->eta(),orgGsfTrk->phi()));
        // fill merged1 histo here
        histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_merged_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_GSF_merged_"+postfix]->Fill(drReco,aWeight);
        histo1d_["id_dR_GSF_sum_"+postfix]->Fill(drReco,aWeight);
        histo1d_["id_dEtaInSeed_merged_"+postfix]->Fill(apair.front()->deltaEtaSeedClusterTrackAtVtx(),aWeight);
        histo1d_["GSF_pt_1st_"+postfix]->Fill(apair.front()->gsfTrack()->pt(),aWeight);
        histo1d_["GSF_pt_2nd_"+postfix]->Fill(addGsfTrk->pt(),aWeight);
        histo1d_["GSF_dPtOverPt_"+postfix]->Fill(addGsfTrk->ptError()/addGsfTrk->pt(),aWeight);

        if ( passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),4) &&
             passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),5) ) {
          histo1d_["id_showershape_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_showershape_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_showershape_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_showershape_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_showershape_merged_"+postfix]->Fill(drReco,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),7) ) {
          histo1d_["id_trackIsoHeep_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_trackIsoHeep_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_trackIsoHeep_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_trackIsoHeep_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_trackIsoHeep_merged_"+postfix]->Fill(drReco,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("heepElectronID-HEEPV70"),8) ) {
          histo1d_["id_caloIsoHeep_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_caloIsoHeep_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_caloIsoHeep_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_caloIsoHeep_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_caloIsoHeep_merged_"+postfix]->Fill(drReco,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("modifiedHeepElectronID"),7) ) {
          histo1d_["id_trackIsoModified_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_trackIsoModified_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_trackIsoModified_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_trackIsoModified_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_trackIsoModified_merged_"+postfix]->Fill(drReco,aWeight);
        }

        if ( passVIDbit(apair.front()->userInt("modifiedHeepElectronID"),8) ) {
          histo1d_["id_caloIsoModified_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_caloIsoModified_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_caloIsoModified_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_caloIsoModified_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_caloIsoModified_merged_"+postfix]->Fill(drReco,aWeight);
        }

        if ( apair.front()->passConversionVeto() )
          histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

        if ( apair.front()->electronID("heepElectronID-HEEPV70") ) {
          histo1d_["id_passHEEP_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_passHEEP_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_passHEEP_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_passHEEP_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_passHEEP_merged_"+postfix]->Fill(drReco,aWeight);
        }

        if ( apair.front()->electronID("modifiedHeepElectronID") ) {
          histo1d_["id_passModifiedHEEP_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_passModifiedHEEP_merged_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_passModifiedHEEP_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
          histo1d_["id_GSF_passModifiedHEEP_"+postfix]->Fill(drReco,aWeight);
          histo1d_["id_GSF_passModifiedHEEP_merged_"+postfix]->Fill(drReco,aWeight);
        }
      }

      return;
    } // apair.size()==1

    if ( postfix=="EB-EB" ) {
      if ( std::abs(apair.at(0)->superCluster()->eta()) > 1.5 || std::abs(apair.at(1)->superCluster()->eta()) > 1.5 )
        return;
    }

    if ( postfix=="EB-EE" ) {
      if ( ( std::abs(apair.at(0)->superCluster()->eta()) > 1.5 && std::abs(apair.at(1)->superCluster()->eta()) > 1.5 ) ||
           ( std::abs(apair.at(0)->superCluster()->eta()) < 1.5 && std::abs(apair.at(1)->superCluster()->eta()) < 1.5 ) )
        return;
    }

    if ( postfix=="EE-EE" ) {
      if ( std::abs(apair.at(0)->superCluster()->eta()) < 1.5 || std::abs(apair.at(1)->superCluster()->eta()) < 1.5 )
        return;
    }

    // default case: apair.size() > 1
    std::sort(apair.begin(),apair.end(),[](const pat::ElectronRef& a, const pat::ElectronRef& b) { return a->et() > b->et(); });

    const double drReco = std::sqrt(reco::deltaR2(apair.at(0)->gsfTrack()->eta(),apair.at(0)->gsfTrack()->phi(),
                                                  apair.at(1)->gsfTrack()->eta(),apair.at(1)->gsfTrack()->phi()));
    // fill heep 1&2 histo here
    histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_resolved_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_GSF_resolved_"+postfix]->Fill(drReco,aWeight);
    histo1d_["id_dR_GSF_sum_"+postfix]->Fill(drReco,aWeight);
    histo1d_["id_dEtaInSeed_resolved_"+postfix]->Fill(apair.at(0)->deltaEtaSeedClusterTrackAtVtx(),0.5*aWeight);
    histo1d_["id_dEtaInSeed_resolved_"+postfix]->Fill(apair.at(1)->deltaEtaSeedClusterTrackAtVtx(),0.5*aWeight);
    histo1d_["GSF_pt_1st_"+postfix]->Fill(apair.at(0)->gsfTrack()->pt(),aWeight);
    histo1d_["GSF_pt_2nd_"+postfix]->Fill(apair.at(1)->gsfTrack()->pt(),aWeight);
    histo1d_["GSF_dPtOverPt_"+postfix]->Fill((*addGsfTrkHandle)[apair.at(0)]->ptError()/(*addGsfTrkHandle)[apair.at(0)]->pt(),aWeight);
    histo1d_["GSF_dPtOverPt_"+postfix]->Fill((*addGsfTrkHandle)[apair.at(1)]->ptError()/(*addGsfTrkHandle)[apair.at(1)]->pt(),aWeight);

    if ( passVIDbit(apair.at(0)->userInt("heepElectronID-HEEPV70"),4) &&
         passVIDbit(apair.at(0)->userInt("heepElectronID-HEEPV70"),5) &&
         passVIDbit(apair.at(1)->userInt("heepElectronID-HEEPV70"),4) &&
         passVIDbit(apair.at(1)->userInt("heepElectronID-HEEPV70"),5) ) {
      histo1d_["id_showershape_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_showershape_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_showershape_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_showershape_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_showershape_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    if ( passVIDbit(apair.at(0)->userInt("heepElectronID-HEEPV70"),7) &&
         passVIDbit(apair.at(1)->userInt("heepElectronID-HEEPV70"),7) ) {
      histo1d_["id_trackIsoHeep_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_trackIsoHeep_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_trackIsoHeep_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_trackIsoHeep_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_trackIsoHeep_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    if ( passVIDbit(apair.at(0)->userInt("heepElectronID-HEEPV70"),8) &&
         passVIDbit(apair.at(1)->userInt("heepElectronID-HEEPV70"),8) ) {
      histo1d_["id_caloIsoHeep_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_caloIsoHeep_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_caloIsoHeep_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_caloIsoHeep_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_caloIsoHeep_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    if ( passVIDbit(apair.at(0)->userInt("modifiedHeepElectronID"),7) &&
         passVIDbit(apair.at(1)->userInt("modifiedHeepElectronID"),7) ) {
      histo1d_["id_trackIsoModified_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_trackIsoModified_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_trackIsoModified_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_trackIsoModified_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_trackIsoModified_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    if ( passVIDbit(apair.at(0)->userInt("modifiedHeepElectronID"),8) &&
         passVIDbit(apair.at(1)->userInt("modifiedHeepElectronID"),8) ) {
      histo1d_["id_caloIsoModified_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_caloIsoModified_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_caloIsoModified_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_caloIsoModified_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_caloIsoModified_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    if ( apair.at(0)->passConversionVeto() &&
         apair.at(1)->passConversionVeto() )
      histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

    if ( apair.at(0)->electronID("heepElectronID-HEEPV70") &&
         apair.at(1)->electronID("heepElectronID-HEEPV70") ) {
      histo1d_["id_passHEEP_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_passHEEP_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_passHEEP_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_passHEEP_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_passHEEP_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    if ( apair.at(0)->electronID("modifiedHeepElectronID") &&
         apair.at(1)->electronID("modifiedHeepElectronID") ) {
      histo1d_["id_passModifiedHEEP_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_passModifiedHEEP_resolved_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_passModifiedHEEP_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
      histo1d_["id_GSF_passModifiedHEEP_"+postfix]->Fill(drReco,aWeight);
      histo1d_["id_GSF_passModifiedHEEP_resolved_"+postfix]->Fill(drReco,aWeight);
    }

    return;
  };

  fillElectrons(pair1,drPrompt1,"all");
  fillElectrons(pair2,drPrompt2,"all");
  fillElectrons(pair1,drPrompt1,"EB-EB");
  fillElectrons(pair2,drPrompt2,"EB-EB");
  fillElectrons(pair1,drPrompt1,"EB-EE");
  fillElectrons(pair2,drPrompt2,"EB-EE");
  fillElectrons(pair1,drPrompt1,"EE-EE");
  fillElectrons(pair2,drPrompt2,"EE-EE");
}

DEFINE_FWK_MODULE(MergedEleSigAnalyzer);
