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
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
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
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dEtaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5EnergyToken_;

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
addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
union5x5dEtaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dEtaIn"))),
union5x5EnergyToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5Energy"))),
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
    histo1d_["id_dEtaInSeed_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_merged_noGsf_"+postfix).c_str(),";#Delta#eta_{in}^{seed};",200,-0.025,0.025);
    histo1d_["id_dEtaInSeed_merged_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_merged_"+postfix).c_str(),";#Delta#eta_{in}^{seed};",200,-0.025,0.025);
    histo1d_["id_dEtaInSeed_resolved_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_resolved_"+postfix).c_str(),";#Delta#eta_{in}^{seed};",200,-0.025,0.025);
    histo1d_["id_dEtaInSeed_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_dEtaInSeed_has2ndTrk_"+postfix).c_str(),";#Delta#eta_{in}^{seed};",200,-0.025,0.025);

    histo1d_["id_dPerpIn_merged_"+postfix] = fs->make<TH1D>(("id_dPerpIn_merged_"+postfix).c_str(),";#Delta u_{in}(union 5x5);",200,-0.025,0.025);
    histo1d_["id_dPerpIn_resolved_"+postfix] = fs->make<TH1D>(("id_dPerpIn_resolved_"+postfix).c_str(),";#Delta u_{in}(union 5x5);",200,-0.025,0.025);
    histo1d_["id_dPerpIn_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_dPerpIn_has2ndTrk_"+postfix).c_str(),";#Delta u_{in}(union 5x5);",200,-0.025,0.025);
    histo1d_["id_minDEtaDPerpIn_merged_"+postfix] = fs->make<TH1D>(("id_minDEtaDPerpIn_merged_"+postfix).c_str(),";min(#Delta#eta_{in}(seed),#Delta u_{in}(union 3x3));",200,-0.025,0.025);
    histo1d_["id_minDEtaDPerpIn_resolved_"+postfix] = fs->make<TH1D>(("id_minDEtaDPerpIn_resolved_"+postfix).c_str(),";min(#Delta#eta_{in}(seed),#Delta u_{in}(union 3x3));",200,-0.025,0.025);
    histo1d_["id_minDEtaDPerpIn_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_minDEtaDPerpIn_has2ndTrk_"+postfix).c_str(),";min(#Delta#eta_{in}(seed),#Delta u_{in}(union 3x3));",200,-0.025,0.025);

    histo1d_["id_E1x5oE5x5_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_E1x5oE5x5_merged_noGsf_"+postfix).c_str(),"Merged w/o GSF E1x5/E5x5",200,0.,1.);
    histo1d_["id_E1x5oE5x5_merged_"+postfix] = fs->make<TH1D>(("id_E1x5oE5x5_merged_"+postfix).c_str(),"Merged w/ GSF E1x5/E5x5",200,0.,1.);
    histo1d_["id_E1x5oE5x5_resolved_"+postfix] = fs->make<TH1D>(("id_E1x5oE5x5_resolved_"+postfix).c_str(),"Resolved E1x5/E5x5",200,0.,1.);

    histo1d_["id_sigIeIe_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_sigIeIe_merged_noGsf_"+postfix).c_str(),"Merged w/o GSF #sigma_{i#eta i#eta}",200,0.,0.05);
    histo1d_["id_sigIeIe_merged_"+postfix] = fs->make<TH1D>(("id_sigIeIe_merged_"+postfix).c_str(),"Merged w/ GSF #sigma_{i#eta i#eta}",200,0.,0.05);
    histo1d_["id_sigIeIe_resolved_"+postfix] = fs->make<TH1D>(("id_sigIeIe_resolved_"+postfix).c_str(),"Resolved #sigma_{i#eta i#eta}",200,0.,0.05);

    histo1d_["id_etaSC_merged_"+postfix] = fs->make<TH1D>(("id_etaSC_merged_"+postfix).c_str(),"Merged w/ GSF #eta_{SC}",200,-3.,3.);
    histo1d_["id_etaSC_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_etaSC_merged_noGsf_"+postfix).c_str(),"Merged w/o GSF #eta_{SC}",200,-3.,3.);

    histo1d_["id_dR_merged_noGsf_"+postfix] = fs->make<TH1D>(("id_dR_merged_noGsf_"+postfix).c_str(),"Merged w/o GSF GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_merged_"+postfix] = fs->make<TH1D>(("id_dR_merged_"+postfix).c_str(),"Merged w/ GSF GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_resolved_"+postfix] = fs->make<TH1D>(("id_dR_resolved_"+postfix).c_str(),"Resolved GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_has2ndTrk_"+postfix] = fs->make<TH1D>(("id_dR_has2ndTrk_"+postfix).c_str(),"w/ 2nd trk GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_all_"+postfix] = fs->make<TH1D>(("id_dR_all_"+postfix).c_str(),"GEN-lv #Delta R",1000,0.,1.);
    histo1d_["id_dR_sum_"+postfix] = fs->make<TH1D>(("id_dR_sum_"+postfix).c_str(),"GEN-lv #Delta R",1000,0.,1.);

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
  };

  histo1d_["GEN_pt_1st"] = fs->make<TH1D>("GEN_pt_1st","GEN-lv Pt of leading lepton",500,0.,1000.);
  histo1d_["GEN_pt_2nd"] = fs->make<TH1D>("GEN_pt_2nd","GEN-lv Pt of subleading lepton",500,0.,1000.);

  addHists("all");
  addHists("EB-EB");
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

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<edm::ValueMap<float>> dPerpInHandle;
  iEvent.getByToken(dPerpInToken_, dPerpInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  edm::Handle<GenEventInfoProduct> genInfo;
  iEvent.getByToken(generatorToken_, genInfo);
  double mcweight = genInfo->weight();

  double aWeight = mcweight/std::abs(mcweight);
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);

  if ( !pvHandle.isValid() || pvHandle->empty() )
    return;

  std::vector<reco::GenParticleRef> promptLeptons;
  const double drThres2 = drThres_*drThres_;

  for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
    const auto& genptc = genptcHandle->refAt(idx);

    if ( ( std::abs(genptc->pdgId())==11 ) &&
         genptc->fromHardProcessFinalState() &&
         ( std::abs(genptc->eta()) < 2.5 ) )
      promptLeptons.push_back(genptc.castTo<reco::GenParticleRef>());
  }

  if (promptLeptons.size() < 2)
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
      double dr2 = reco::deltaR2(aEle->gsfTrack()->eta(),aEle->gsfTrack()->phi(),genptc->eta(),genptc->phi());

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

  histo1d_["GEN_pt_1st"]->Fill( promptLeptons.at(0)->pt(),aWeight );
  histo1d_["GEN_pt_2nd"]->Fill( promptLeptons.at(1)->pt(),aWeight );

  // for eemm & eeee final states
  if ( pair1.empty() )
    return;

  if ( promptLeptons.at(0)->pt() < ptThres_ || promptLeptons.at(1)->pt() < ptThres2nd_ )
    return;

  if ( pair1.size() > 1 && promptLeptons.at(1)->pt() < ptThres_ )
    return;

  // for eeee final state
  if ( promptLeptons.size() > 2 ) {
    if (promptLeptons.size() != 4)
      return;

    histo1d_["GEN_pt_1st"]->Fill( promptLeptons.at(2)->pt(), aWeight );
    histo1d_["GEN_pt_2nd"]->Fill( promptLeptons.at(3)->pt(), aWeight );

    if ( pair2.empty() )
      return;

    if ( promptLeptons.at(2)->pt() < ptThres_ || promptLeptons.at(3)->pt() < ptThres2nd_ )
      return;

    if ( pair2.size() > 1 && promptLeptons.at(3)->pt() < ptThres_ )
      return;
  }

  auto passVIDbit = [] (const int32_t bitmap, unsigned cutIdx) -> bool {
    int32_t mask = 0x00000001 << cutIdx;
    int32_t result = mask & bitmap;

    return result!=0x00000000;
  };

  // auto passMergedElectronID = [&] (const pat::ElectronRef& aEle) {
  //   const float eta1stGSF = -(aEle->deltaEtaSeedClusterTrackAtVtx() - aEle->superCluster()->seed()->eta());
  //   const float u5x5Eta = (*union5x5dEtaInHandle)[aEle] + eta1stGSF;
  //   const float u5x5Et = (*union5x5EnergyHandle)[aEle]/std::cosh(u5x5Eta);
  //
  //   if ( aEle->electronID("mvaMergedElectron") ) {
  //     if ( u5x5Et > 50. )
  //       return true;
  //   }
  //
  //   return false;
  // };

  auto fillEfficiency = [&passVIDbit,&aWeight,this] (const std::vector<pat::ElectronRef>& apair,
                                            const double drPrompt,
                                            const std::string& postfix,
                                            std::initializer_list<std::string> histList,
                                            const std::string& idName,
                                            std::initializer_list<int> bitmaps = {}) -> void {
    for (const auto& aEle : apair) {
      bool pass = true;

      if (bitmaps.size()==0) {
        pass = static_cast<bool>(aEle->electronID(idName));
      } else {
        for (int abit : bitmaps) {
          if (!passVIDbit(aEle->userInt(idName),abit)) {
            pass = false;
            break;
          }
        }
      }

      if (pass) {
        for (const auto& name : histList) {
          histo1d_[name+"_"+postfix]->Fill(drPrompt,aWeight);
        }
      }
    }
  };

  auto fillElectrons = [&] (std::vector<pat::ElectronRef>& apair, const double drPrompt,
                            const std::vector<reco::GenParticleRef>& promptPair, const std::string& postfix) {
    if ( apair.size()!=1 && apair.size()!=2 )
      return;

    histo1d_["id_dR_all_"+postfix]->Fill(drPrompt,aWeight);

    if ( apair.size()==1 ) {
      auto addGsfTrk = (*addGsfTrkHandle)[apair.front()];
      auto addPackedCand = (*addPackedCandHandle)[apair.front()];
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

      if ( addGsfTrk==orgGsfTrk && addPackedCand.isNull() ) { // ME w/o add GSF
        // fill merged2 histo here
        histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_merged_noGsf_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dEtaInSeed_merged_noGsf_"+postfix]->Fill(apair.front()->deltaEtaSeedClusterTrackAtVtx(),aWeight);
        histo1d_["id_etaSC_merged_noGsf_"+postfix]->Fill(apair.front()->superCluster()->eta(),aWeight);

        const double e5x5 = apair.front()->full5x5_e5x5();
        histo1d_["id_E1x5oE5x5_merged_noGsf_"+postfix]->Fill(e5x5!=0 ? apair.front()->full5x5_e1x5()/e5x5 : 0, aWeight);
        histo1d_["id_sigIeIe_merged_noGsf_"+postfix]->Fill(apair.front()->full5x5_sigmaIetaIeta(), aWeight);

        if ( apair.front()->passConversionVeto() )
          histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

        fillEfficiency(apair,drPrompt,postfix,
                       {"id_showershape","id_showershape_merged_noGsf"},
                       "heepElectronID-HEEPV70",{4,5});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_trackIsoHeep","id_trackIsoHeep_merged_noGsf"},
                       "heepElectronID-HEEPV70",{7});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_caloIsoHeep","id_caloIsoHeep_merged_noGsf"},
                       "heepElectronID-HEEPV70",{8});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_passHEEP","id_passHEEP_merged_noGsf"},
                       "heepElectronID-HEEPV70");
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_trackIsoModified","id_trackIsoModified_merged_noGsf"},
                       "modifiedHeepElectronID",{7});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_caloIsoModified","id_caloIsoModified_merged_noGsf"},
                       "modifiedHeepElectronID",{8});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_passModifiedHEEP","id_passModifiedHEEP_merged_noGsf"},
                       "modifiedHeepElectronID");
      } else { // ME w/ GSF
        // find whether add GSF track has a corresponding electron
        bool notMerged = false;

        if (addPackedCand.isNull()) {
          for (unsigned int jdx = 0; jdx < eleHandle->size(); ++jdx) {
            const auto& secEle = eleHandle->refAt(jdx);
            const auto& secGsfTrk = secEle->gsfTrack();

            if (apair.front()==secEle.castTo<pat::ElectronRef>())
              continue;

            if ( addGsfTrk==secGsfTrk ) {
              notMerged = true;
              break;
            }
          }
        }

        if ( notMerged )
          return;

        // first look for GSF and then packedPFcand
        const reco::Track* addTrk = addGsfTrk.get();

        if ( addGsfTrk==orgGsfTrk && addPackedCand.isNonnull() )
          addTrk = addPackedCand->bestTrack();

        double drcand10 = reco::deltaR2(addTrk->eta(),addTrk->phi(),promptPair.at(0)->eta(),promptPair.at(0)->phi());
        double drcand11 = reco::deltaR2(addTrk->eta(),addTrk->phi(),promptPair.at(1)->eta(),promptPair.at(1)->phi());

        if ( drcand10 < drThres2 || drcand11 < drThres2 ) {
          double val = (*dPerpInHandle)[apair.front()];
          double dEtaInSeed = apair.front()->deltaEtaSeedClusterTrackAtVtx();
          val = (std::abs(val) < std::abs(dEtaInSeed)) ? val : dEtaInSeed;
          histo1d_["id_dPerpIn_merged_"+postfix]->Fill( (*dPerpInHandle)[apair.front()] , aWeight );
          histo1d_["id_dPerpIn_has2ndTrk_"+postfix]->Fill( (*dPerpInHandle)[apair.front()] , aWeight );
          histo1d_["id_minDEtaDPerpIn_merged_"+postfix]->Fill( val , aWeight );
          histo1d_["id_minDEtaDPerpIn_has2ndTrk_"+postfix]->Fill( val , aWeight );
        }

        // fill merged1 histo here
        histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_merged_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dR_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
        histo1d_["id_dEtaInSeed_merged_"+postfix]->Fill(apair.front()->deltaEtaSeedClusterTrackAtVtx(),aWeight);
        histo1d_["id_dEtaInSeed_has2ndTrk_"+postfix]->Fill(apair.front()->deltaEtaSeedClusterTrackAtVtx(),aWeight);
        histo1d_["id_etaSC_merged_"+postfix]->Fill(apair.front()->superCluster()->eta(),aWeight);

        const double e5x5 = apair.front()->full5x5_e5x5();
        histo1d_["id_E1x5oE5x5_merged_"+postfix]->Fill(e5x5!=0 ? apair.front()->full5x5_e1x5()/e5x5 : 0, aWeight);
        histo1d_["id_sigIeIe_merged_"+postfix]->Fill(apair.front()->full5x5_sigmaIetaIeta(), aWeight);

        if ( apair.front()->passConversionVeto() )
          histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

        fillEfficiency(apair,drPrompt,postfix,
                       {"id_showershape","id_showershape_merged","id_showershape_has2ndTrk"},
                       "heepElectronID-HEEPV70",{4,5});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_trackIsoHeep","id_trackIsoHeep_merged","id_trackIsoHeep_has2ndTrk"},
                       "heepElectronID-HEEPV70",{7});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_caloIsoHeep","id_caloIsoHeep_merged","id_caloIsoHeep_has2ndTrk"},
                       "heepElectronID-HEEPV70",{8});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_passHEEP","id_passHEEP_merged","id_passHEEP_has2ndTrk"},
                       "heepElectronID-HEEPV70");
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_trackIsoModified","id_trackIsoModified_merged","id_trackIsoModified_has2ndTrk"},
                       "modifiedHeepElectronID",{7});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_caloIsoModified","id_caloIsoModified_merged","id_caloIsoModified_has2ndTrk"},
                       "modifiedHeepElectronID",{8});
        fillEfficiency(apair,drPrompt,postfix,
                       {"id_passModifiedHEEP","id_passModifiedHEEP_merged","id_passModifiedHEEP_has2ndTrk"},
                       "modifiedHeepElectronID");
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

    // fill heep 1&2 histo here
    histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_sum_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_resolved_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_resolved_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dR_has2ndTrk_"+postfix]->Fill(drPrompt,aWeight);
    histo1d_["id_dEtaInSeed_resolved_"+postfix]->Fill(apair.at(0)->deltaEtaSeedClusterTrackAtVtx(),aWeight);
    histo1d_["id_dEtaInSeed_resolved_"+postfix]->Fill(apair.at(1)->deltaEtaSeedClusterTrackAtVtx(),aWeight);
    histo1d_["id_dEtaInSeed_has2ndTrk_"+postfix]->Fill(apair.at(0)->deltaEtaSeedClusterTrackAtVtx(),aWeight);
    histo1d_["id_dEtaInSeed_has2ndTrk_"+postfix]->Fill(apair.at(1)->deltaEtaSeedClusterTrackAtVtx(),aWeight);

    const double e5x5_1 = apair.at(0)->full5x5_e5x5();
    histo1d_["id_E1x5oE5x5_resolved_"+postfix]->Fill(e5x5_1!=0 ? apair.at(0)->full5x5_e1x5()/e5x5_1 : 0, aWeight);
    histo1d_["id_sigIeIe_resolved_"+postfix]->Fill(apair.at(0)->full5x5_sigmaIetaIeta(), aWeight);

    const double e5x5_2 = apair.at(1)->full5x5_e5x5();
    histo1d_["id_E1x5oE5x5_resolved_"+postfix]->Fill(e5x5_2!=0 ? apair.at(1)->full5x5_e1x5()/e5x5_2 : 0, aWeight);
    histo1d_["id_sigIeIe_resolved_"+postfix]->Fill(apair.at(1)->full5x5_sigmaIetaIeta(), aWeight);

    double val1 = (*dPerpInHandle)[apair.at(0)];
    const double dEtaInSeed1 = apair.at(0)->deltaEtaSeedClusterTrackAtVtx();
    val1 = (std::abs(val1) < std::abs(dEtaInSeed1)) ? val1 : dEtaInSeed1;
    histo1d_["id_dPerpIn_resolved_"+postfix]->Fill( (*dPerpInHandle)[apair.at(0)] , aWeight );
    histo1d_["id_dPerpIn_has2ndTrk_"+postfix]->Fill( (*dPerpInHandle)[apair.at(0)] , aWeight );
    histo1d_["id_minDEtaDPerpIn_resolved_"+postfix]->Fill( val1 , aWeight );
    histo1d_["id_minDEtaDPerpIn_has2ndTrk_"+postfix]->Fill( val1 , aWeight );

    double val2 = (*dPerpInHandle)[apair.at(1)];
    const double dEtaInSeed2 = apair.at(1)->deltaEtaSeedClusterTrackAtVtx();
    val2 = (std::abs(val2) < std::abs(dEtaInSeed2)) ? val2 : dEtaInSeed2;
    histo1d_["id_dPerpIn_resolved_"+postfix]->Fill( (*dPerpInHandle)[apair.at(1)] , aWeight );
    histo1d_["id_dPerpIn_has2ndTrk_"+postfix]->Fill( (*dPerpInHandle)[apair.at(1)] , aWeight );
    histo1d_["id_minDEtaDPerpIn_resolved_"+postfix]->Fill( val2 , aWeight );
    histo1d_["id_minDEtaDPerpIn_has2ndTrk_"+postfix]->Fill( val2 , aWeight );

    if ( apair.at(0)->passConversionVeto() )
      histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

    if ( apair.at(1)->passConversionVeto() )
      histo1d_["id_conversionVeto_"+postfix]->Fill(drPrompt,aWeight);

    fillEfficiency(apair,drPrompt,postfix,
                   {"id_showershape","id_showershape_resolved","id_showershape_has2ndTrk"},
                   "heepElectronID-HEEPV70",{4,5});
    fillEfficiency(apair,drPrompt,postfix,
                   {"id_trackIsoHeep","id_trackIsoHeep_resolved","id_trackIsoHeep_has2ndTrk"},
                   "heepElectronID-HEEPV70",{7});
    fillEfficiency(apair,drPrompt,postfix,
                   {"id_caloIsoHeep","id_caloIsoHeep_resolved","id_caloIsoHeep_has2ndTrk"},
                   "heepElectronID-HEEPV70",{8});
    fillEfficiency(apair,drPrompt,postfix,
                   {"id_passHEEP","id_passHEEP_resolved","id_passHEEP_has2ndTrk"},
                   "heepElectronID-HEEPV70");
    fillEfficiency(apair,drPrompt,postfix,
                   {"id_trackIsoModified","id_trackIsoModified_resolved","id_trackIsoModified_has2ndTrk"},
                   "modifiedHeepElectronID",{7});
    fillEfficiency(apair,drPrompt,postfix,
                   {"id_caloIsoModified","id_caloIsoModified_resolved","id_caloIsoModified_has2ndTrk"},
                   "modifiedHeepElectronID",{8});
    fillEfficiency(apair,drPrompt,postfix,
                   {"id_passModifiedHEEP","id_passModifiedHEEP_resolved","id_passModifiedHEEP_has2ndTrk"},
                   "modifiedHeepElectronID");

    return;
  };

  std::vector<reco::GenParticleRef> promptPair1 = {promptLeptons.at(0),promptLeptons.at(1)};
  const double drPrompt1 = std::sqrt( reco::deltaR2(promptLeptons.at(0)->eta(),promptLeptons.at(0)->phi(),
                                                    promptLeptons.at(1)->eta(),promptLeptons.at(1)->phi()) );
  fillElectrons(pair1,drPrompt1,promptPair1,"all");
  fillElectrons(pair1,drPrompt1,promptPair1,"EB-EB");
  fillElectrons(pair1,drPrompt1,promptPair1,"EE-EE");

  if ( promptLeptons.size()==4 ) {
    std::vector<reco::GenParticleRef> promptPair2 = {promptLeptons.at(2),promptLeptons.at(3)};
    const double drPrompt2 = std::sqrt( reco::deltaR2(promptLeptons.at(2)->eta(),promptLeptons.at(2)->phi(),
                                                      promptLeptons.at(3)->eta(),promptLeptons.at(3)->phi()) );
    fillElectrons(pair2,drPrompt2,promptPair2,"all");
    fillElectrons(pair2,drPrompt2,promptPair2,"EB-EB");
    fillElectrons(pair2,drPrompt2,promptPair2,"EE-EE");
  }
}

DEFINE_FWK_MODULE(MergedEleSigAnalyzer);
