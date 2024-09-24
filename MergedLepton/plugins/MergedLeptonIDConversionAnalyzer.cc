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
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"

//#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "CommonTools/Egamma/interface/ConversionTools.h"

#include "TH1D.h"
#include "TH2F.h"
#include "TTree.h"
#include "TString.h"
#include "TMath.h"

class MergedLeptonIDConversionAnalyzer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit MergedLeptonIDConversionAnalyzer(const edm::ParameterSet&);
  virtual ~MergedLeptonIDConversionAnalyzer() {}

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  const edm::EDGetTokenT<edm::View<pat::Electron>> srcEle_;
  const edm::EDGetTokenT<edm::View<pat::Muon>> srcMuon_;
  const edm::EDGetTokenT<edm::View<pat::Photon>> srcPhoton_;
  const edm::EDGetTokenT<edm::View<reco::Conversion>> conversionsToken_;
  const edm::EDGetTokenT<edm::View<reco::Conversion>> singleLegConversionsToken_;
  const edm::EDGetTokenT<edm::View<reco::GenParticle>> srcGenPtc_;
  const edm::EDGetTokenT<edm::View<reco::Vertex>> pvToken_;
  const edm::EDGetTokenT<edm::View<PileupSummaryInfo>> pileupToken_;
  const edm::EDGetTokenT<reco::BeamSpot> beamspotToken_;
  const edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
  const edm::EDGetTokenT<LHEEventProduct> lheToken_;
  const edm::EDGetTokenT<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkToken_;
  const edm::EDGetTokenT<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> trkIsoMapToken_;

  const edm::EDGetTokenT<edm::ValueMap<float>> alphaTrackToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> alphaCaloToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> dPerpInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> normDParaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dEtaInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5dPhiInToken_;
  const edm::EDGetTokenT<edm::ValueMap<float>> union5x5EnergyToken_;

  const edm::EDGetTokenT<edm::TriggerResults> triggerToken_;
  const edm::EDGetTokenT<edm::View<pat::TriggerObjectStandAlone>> triggerobjectsToken_;

  const edm::EDGetTokenT<double> prefweight_token;

  const std::vector<std::string> trigList_;

  const edm::FileInPath purwgtPath_;

  const bool isMC_;
  const double ptThresTrig_;
  const double ptThres_;
  const double drThres_;

  // PDG mass & error
  const double elmass_ = 0.0005109989461;
  const double elmassErr_ = 0.0000000000031;
  const double mumass_ = 0.1056583745;
  const double mumassErr_ = 0.0000000024;

  std::unique_ptr<TFile> purwgtFile_;
  TH1D* purwgt_;

  std::map<std::string,TH1*> histo1d_;
  std::map<std::string,TH2*> histo2d_;

  TTree* tree_ = nullptr;
  float invM_ = -1.;
  float u5x5Et_ = -1.;
  float wgt_ = 0.;
  int passME_ = -1;
};

MergedLeptonIDConversionAnalyzer::MergedLeptonIDConversionAnalyzer(const edm::ParameterSet& iConfig) :
srcEle_(consumes<edm::View<pat::Electron>>(iConfig.getParameter<edm::InputTag>("srcEle"))),
srcMuon_(consumes<edm::View<pat::Muon>>(iConfig.getParameter<edm::InputTag>("srcMuon"))),
srcPhoton_(consumes<edm::View<pat::Photon>>(iConfig.getParameter<edm::InputTag>("srcPhoton"))),
conversionsToken_(consumes<edm::View<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
singleLegConversionsToken_(consumes<edm::View<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("singleLegConversions"))),
srcGenPtc_(consumes<edm::View<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("srcGenPtc"))),
pvToken_(consumes<edm::View<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("srcPv"))),
pileupToken_(consumes<edm::View<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupSummary"))),
beamspotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamSpot"))),
generatorToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("generator"))),
lheToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvent"))),
addGsfTrkToken_(consumes<edm::ValueMap<reco::GsfTrackRef>>(iConfig.getParameter<edm::InputTag>("addGsfTrkMap"))),
addPackedCandToken_(consumes<edm::ValueMap<pat::PackedCandidateRef>>(iConfig.getParameter<edm::InputTag>("addPackedCandMap"))),
trkIsoMapToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("trkIsoMap"))),
alphaTrackToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaTrack"))),
alphaCaloToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("alphaCalo"))),
dPerpInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("dPerpIn"))),
normDParaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("normalizedDParaIn"))),
union5x5dEtaInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dEtaIn"))),
union5x5dPhiInToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5dPhiIn"))),
union5x5EnergyToken_(consumes<edm::ValueMap<float>>(iConfig.getParameter<edm::InputTag>("union5x5Energy"))),
triggerToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
triggerobjectsToken_(consumes<edm::View<pat::TriggerObjectStandAlone>>(iConfig.getParameter<edm::InputTag>("triggerObjects"))),
prefweight_token(consumes<double>(edm::InputTag("prefiringweight:nonPrefiringProb"))),
trigList_(iConfig.getParameter<std::vector<std::string>>("trigList")),
purwgtPath_(iConfig.getParameter<edm::FileInPath>("PUrwgt")),
isMC_(iConfig.getParameter<bool>("isMC")),
ptThresTrig_(iConfig.getParameter<double>("ptThresTrig")),
ptThres_(iConfig.getParameter<double>("ptThres")),
drThres_(iConfig.getParameter<double>("drThres")) {
  usesResource("TFileService");
}

void MergedLeptonIDConversionAnalyzer::beginJob() {
  TH1::SetDefaultSumw2();
  edm::Service<TFileService> fs;

  purwgtFile_ = std::make_unique<TFile>(purwgtPath_.fullPath().c_str(),"READ");
  purwgt_ = static_cast<TH1D*>(purwgtFile_->Get("PUrwgt"));

  histo1d_["totWeightedSum"] = fs->make<TH1D>("totWeightedSum","totWeightedSum",1,0.,1.);
  histo1d_["cutflow"] = fs->make<TH1D>("cutflow","cutflow",30,0.,30.);
  histo1d_["cutflow_tnpPair"] = fs->make<TH1D>("cutflow_tnpPair","cutflow_tnpPair",30,0.,30.);
  histo1d_["mva_HasTrkEB"] = fs->make<TH1D>("mva_HasTrkEB","MVA score",200,-1.,1.);
  histo1d_["mva_NoneEt2EB"] = fs->make<TH1D>("mva_NoneEt2EB","MVA score",200,-1.,1.);
  histo1d_["nPV"] = fs->make<TH1D>("nPV","nPV",99,0.,99.);
  histo1d_["PUsummary"] = fs->make<TH1D>("PUsummary","PUsummary",99,0.,99.);

  histo1d_["LHE_leadMuPt"] = fs->make<TH1D>("LHE_leadMuPt","LHE_leadMuPt",200,0.,200.);
  histo1d_["LHE_subleadMuPt"] = fs->make<TH1D>("LHE_subleadMuPt","LHE_subleadMuPt",200,0.,200.);
  histo1d_["GEN_leadMuPt"] = fs->make<TH1D>("GEN_leadMuPt","GEN_leadMuPt",200,0.,200.);
  histo1d_["GEN_subleadMuPt"] = fs->make<TH1D>("GEN_subleadMuPt","GEN_subleadMuPt",200,0.,200.);
  histo1d_["GEN_gammaPt"] = fs->make<TH1D>("GEN_gammaPt","GEN_gammaPt",200,0.,200.);
  histo1d_["GEN_dR"] = fs->make<TH1D>("GEN_dR","GEN_dR",300,0.,1.5);
  histo1d_["GEN_invM"] = fs->make<TH1D>("GEN_invM","GEN_invM",250,0.,250.);

  tree_ = fs->make<TTree>("dielTree","dielTree");
  tree_->Branch("invM",&invM_,"invM/F");
  tree_->Branch("u5x5Et",&u5x5Et_,"u5x5Et/F");
  tree_->Branch("wgt",&wgt_,"wgt/F");
  tree_->Branch("passME",&passME_,"passME/I");

  auto createProbeHisto = [&,this] (const std::string& prefix) {
    histo1d_[prefix+"_probeEle_MMG_invM"] = fs->make<TH1F>((prefix+"_probeEle_MMG_invM").c_str(),";Mass [GeV];",500,0.,500.);
    histo1d_[prefix+"_probeEle_MMG_rap"] = fs->make<TH1F>((prefix+"_probeEle_MMG_rap").c_str(),";rapidity;",200,-2.5,2.5);
    histo1d_[prefix+"_probeEle_MM_invM"] = fs->make<TH1F>((prefix+"_probeEle_MM_invM").c_str(),";Mass [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_mu1_pt"] = fs->make<TH1F>((prefix+"_probeEle_mu1_pt").c_str(),";p_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_mu1_eta"] = fs->make<TH1F>((prefix+"_probeEle_mu1_eta").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"_probeEle_mu2_pt"] = fs->make<TH1F>((prefix+"_probeEle_mu2_pt").c_str(),";p_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_mu2_eta"] = fs->make<TH1F>((prefix+"_probeEle_mu2_eta").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"_probeEle_probe_phoEt"] = fs->make<TH1F>((prefix+"_probeEle_probe_phoEt").c_str(),";E_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_probe_eleEt"] = fs->make<TH1F>((prefix+"_probeEle_probe_eleEt").c_str(),";E_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_probe_eleU5x5Et"] = fs->make<TH1F>((prefix+"_probeEle_probe_eleU5x5Et").c_str(),";E_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_probe_eta"] = fs->make<TH1F>((prefix+"_probeEle_probe_eta").c_str(),";#eta;",200,-2.5,2.5);
    histo1d_[prefix+"_probeEle_mu1_dR"] = fs->make<TH1F>((prefix+"_probeEle_mu1_dR").c_str(),";#Delta R(#mu,#gamma);",320,0.0,3.2);
    histo1d_[prefix+"_probeEle_mu2_dR"] = fs->make<TH1F>((prefix+"_probeEle_mu2_dR").c_str(),";#Delta R(#mu,#gamma);",320,0.0,3.2);

    histo2d_[prefix+"_probeEle_MMG_MM_invM"] = fs->make<TH2F>((prefix+"_probeEle_MMG_MM_invM").c_str(),"M_{#mu#mu#gamma};M_{#mu#mu};",100,0.,200.,100,0.,200.);

    histo1d_[prefix+"_probeEle_diel_dEta"] = fs->make<TH1D>((prefix+"_probeEle_diel_dEta").c_str(),";#Delta#eta;",400,-0.05,0.05);
    histo1d_[prefix+"_probeEle_diel_dPhi"] = fs->make<TH1D>((prefix+"_probeEle_diel_dPhi").c_str(),";#Delta#phi;",400,-0.1,0.1);
    histo1d_[prefix+"_probeEle_diel_dR"] = fs->make<TH1D>((prefix+"_probeEle_diel_dR").c_str(),";#Delta R;",300,0.,0.3);
    histo1d_[prefix+"_probeEle_diel_alphaCalo"] = fs->make<TH1D>((prefix+"_probeEle_diel_alphaCalo").c_str(),";#alpha_{calo};",360,-1.8,1.8);
    histo1d_[prefix+"_probeEle_diel_alphaTrk"] = fs->make<TH1D>((prefix+"_probeEle_diel_alphaTrk").c_str(),";#alpha_{trk};",400,-1.,1.);
    histo1d_[prefix+"_probeEle_diel_normDParaIn"] = fs->make<TH1D>((prefix+"_probeEle_diel_normDParaIn").c_str(),";#Delta v_{in}/#Delta R;",200,-0.5,1.5);
    histo1d_[prefix+"_probeEle_diel_HoE"] = fs->make<TH1D>((prefix+"_probeEle_diel_HoE").c_str(),";H/E;",150,0.,0.15);
    histo1d_[prefix+"_probeEle_diel_sigIeIe"] = fs->make<TH1D>((prefix+"_probeEle_diel_sigIeIe").c_str(),";#sigma_{i#eta i#eta};",300,0.,0.03);
    histo1d_[prefix+"_probeEle_diel_dEtaIn"] = fs->make<TH1D>((prefix+"_probeEle_diel_dEtaIn").c_str(),";#Delta#eta_{in};",400,-0.05,0.05);
    histo1d_[prefix+"_probeEle_diel_dEtaInSeed"] = fs->make<TH1D>((prefix+"_probeEle_diel_dEtaInSeed").c_str(),";#Delta#eta_{in seed};",400,-0.05,0.05);
    histo1d_[prefix+"_probeEle_diel_dPhiIn"] = fs->make<TH1D>((prefix+"_probeEle_diel_dPhiIn").c_str(),";#Delta#phi_{in};",400,-0.1,0.1);
    histo1d_[prefix+"_probeEle_diel_R9"] = fs->make<TH1D>((prefix+"_probeEle_diel_R9").c_str(),";R9;",200,0.,1.);
    histo1d_[prefix+"_probeEle_diel_EoverP"] = fs->make<TH1D>((prefix+"_probeEle_diel_EoverP").c_str(),";E/p;",250,0.,5.);
    histo1d_[prefix+"_probeEle_diel_trackIso"] = fs->make<TH1D>((prefix+"_probeEle_diel_trackIso").c_str(),";#Sigma p_{T iso};",300,-10.,50.);
    histo1d_[prefix+"_probeEle_diel_dPerpIn"] = fs->make<TH1D>((prefix+"_probeEle_diel_dPerpIn").c_str(),";#Delta u_{in};",400,-0.05,0.05);
    histo1d_[prefix+"_probeEle_diel_u5x5dEtaIn"] = fs->make<TH1D>((prefix+"_probeEle_diel_u5x5dEtaIn").c_str(),";#Delta #eta_{in}(u5x5);",500,-0.1,0.1);
    histo1d_[prefix+"_probeEle_diel_u5x5dPhiIn"] = fs->make<TH1D>((prefix+"_probeEle_diel_u5x5dPhiIn").c_str(),";#Delta #phi_{in}(u5x5);",500,-0.1,0.1);

    histo1d_[prefix+"_probeEle_diel_invM_ee"] = fs->make<TH1F>((prefix+"_probeEle_diel_invM_ee").c_str(),";Mass [GeV];",500,0.,5.);
    histo1d_[prefix+"_probeEle_diel_pt"] = fs->make<TH1F>((prefix+"_probeEle_diel_pt").c_str(),";p_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_diel_etSC"] = fs->make<TH1D>((prefix+"_probeEle_diel_etSC").c_str(),";E_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_1st_pt"] = fs->make<TH1F>((prefix+"_probeEle_1st_pt").c_str(),";p_{T} [GeV];",200,0.,200.);
    histo1d_[prefix+"_probeEle_2nd_pt"] = fs->make<TH1F>((prefix+"_probeEle_2nd_pt").c_str(),";p_{T} [GeV];",200,0.,200.);

    histo1d_[prefix+"_probeEle_mva_HasTrkEB"] = fs->make<TH1D>((prefix+"_probeEle_mva_HasTrkEB").c_str(),"MVA score",200,-1.,1.);
    histo1d_[prefix+"_probeEle_mva_NoneEt2EB"] = fs->make<TH1D>((prefix+"_probeEle_mva_NoneEt2EB").c_str(),"MVA score",200,-1.,1.);

    histo1d_[prefix+"_probeEle_conv_invM"] = fs->make<TH1F>((prefix+"_probeEle_conv_invM").c_str(),";Mass [GeV];",500,0.,5.);
    histo1d_[prefix+"_probeEle_conv_dPhiVtx"] = fs->make<TH1D>((prefix+"_probeEle_conv_dPhiVtx").c_str(),";#Delta#phi_{vtx};",400,-0.1,0.1);
    histo1d_[prefix+"_probeEle_conv_lxy"] = fs->make<TH1D>((prefix+"_probeEle_conv_lxy").c_str(),";lxy;",200,0.,100.);
    histo1d_[prefix+"_probeEle_conv_nHitsBeforeVtx_1st"] = fs->make<TH1D>((prefix+"_probeEle_conv_nHitsBeforeVtx_1st").c_str(),";nHits;",20,-1.,19.);
    histo1d_[prefix+"_probeEle_conv_nHitsBeforeVtx_2nd"] = fs->make<TH1D>((prefix+"_probeEle_conv_nHitsBeforeVtx_2nd").c_str(),";nHits;",20,-1.,19.);
    histo1d_[prefix+"_probeEle_conv_nSharedHits"] = fs->make<TH1D>((prefix+"_probeEle_conv_nSharedHits").c_str(),";nHits;",20,0.,20.);
  };

  createProbeHisto("MB_FSR");
  createProbeHisto("PassMergedEle");
  createProbeHisto("FailMergedEle");
  createProbeHisto("MB_ISR");

  createProbeHisto("singleLeg_MB_FSR");
  createProbeHisto("singleLeg_PassMergedEle");
  createProbeHisto("singleLeg_FailMergedEle");
  createProbeHisto("singleLeg_MB_ISR");

  createProbeHisto("singleLeg_Et10to15_MB_FSR");
  createProbeHisto("singleLeg_Et10to15_PassMergedEle");
  createProbeHisto("singleLeg_Et10to15_FailMergedEle");
  createProbeHisto("singleLeg_Et10to15_MB_ISR");

  createProbeHisto("singleLeg_Et15to20_MB_FSR");
  createProbeHisto("singleLeg_Et15to20_PassMergedEle");
  createProbeHisto("singleLeg_Et15to20_FailMergedEle");
  createProbeHisto("singleLeg_Et15to20_MB_ISR");

  createProbeHisto("singleLeg_Et20to30_MB_FSR");
  createProbeHisto("singleLeg_Et20to30_PassMergedEle");
  createProbeHisto("singleLeg_Et20to30_FailMergedEle");
  createProbeHisto("singleLeg_Et20to30_MB_ISR");

  createProbeHisto("singleLeg_Et30toInf_MB_FSR");
  createProbeHisto("singleLeg_Et30toInf_PassMergedEle");
  createProbeHisto("singleLeg_Et30toInf_FailMergedEle");
  createProbeHisto("singleLeg_Et30toInf_MB_ISR");
}

void MergedLeptonIDConversionAnalyzer::endJob() {
  purwgtFile_->Close();
}

void MergedLeptonIDConversionAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  edm::Handle<edm::View<reco::Vertex>> pvHandle;
  iEvent.getByToken(pvToken_, pvHandle);

  double aWeight = 1.;
  histo1d_["totWeightedSum"]->Fill(0.5,aWeight);
  histo1d_["cutflow"]->Fill(0.5,aWeight);

  if (isMC_) {
    edm::Handle<GenEventInfoProduct> genInfo;
    iEvent.getByToken(generatorToken_, genInfo);
    double mcweight = genInfo->weight();

    edm::Handle<LHEEventProduct> lheEventHandle;
    iEvent.getByToken(lheToken_, lheEventHandle);

    const auto& hepeup = lheEventHandle->hepeup();
    const auto& pup = hepeup.PUP;
    int lep = -1, lepBar = -1;

    for (unsigned int i = 0, n = pup.size(); i < n; ++i) {
      int idabs = std::abs(hepeup.IDUP[i]);

      if (idabs == 13)
        (hepeup.IDUP[i] > 0 ? lep : lepBar) = i;
    }

    std::pair<int, int> v(0, 0);
    if (lep != -1 && lepBar != -1)
      v = std::make_pair(lep, lepBar);
    else
      return;

    double lheMuPt1 = -1., lheMuPt2 = -1.;

    if (v.first != -1 && v.second != -1) {
      lheMuPt1 = std::hypot( pup[v.first][0], pup[v.first][1] );
      lheMuPt2 = std::hypot( pup[v.second][0], pup[v.second][1] );
    }

    histo1d_["LHE_leadMuPt"]->Fill( std::max(lheMuPt1,lheMuPt2) , aWeight);
    histo1d_["LHE_subleadMuPt"]->Fill( std::min(lheMuPt1,lheMuPt2) , aWeight);

    edm::Handle<edm::View<reco::GenParticle>> genptcHandle;
    iEvent.getByToken(srcGenPtc_, genptcHandle);

    std::vector<reco::GenParticleRef> hardprocMuons;
    std::vector<reco::GenParticleRef> promptPhotons;

    for (unsigned int idx = 0; idx < genptcHandle->size(); ++idx) {
      const auto& genptc = genptcHandle->refAt(idx);

      if ( std::abs(genptc->pdgId())==13 && genptc->fromHardProcessFinalState() )
        hardprocMuons.push_back( genptc.castTo<reco::GenParticleRef>() );

      if ( std::abs(genptc->pdgId())==22 && genptc->isPromptFinalState() &&
           genptc->pt() > 10. && std::abs(genptc->eta()) < 2.6 )
        promptPhotons.push_back( genptc.castTo<reco::GenParticleRef>() );
    }

    auto sortByPt = [] (const reco::GenParticleRef& a, const reco::GenParticleRef& b) {
      return a->pt() > b->pt();
    };

    std::sort(hardprocMuons.begin(),hardprocMuons.end(),sortByPt);
    std::sort(promptPhotons.begin(),promptPhotons.end(),sortByPt);

    if ( !promptPhotons.empty() && hardprocMuons.size() > 1 ) {
      const auto lvecM1 = math::PtEtaPhiMLorentzVector(hardprocMuons.at(0)->pt(),
                                                       hardprocMuons.at(0)->eta(),
                                                       hardprocMuons.at(0)->phi(),
                                                       mumass_);
      const auto lvecM2 = math::PtEtaPhiMLorentzVector(hardprocMuons.at(1)->pt(),
                                                       hardprocMuons.at(1)->eta(),
                                                       hardprocMuons.at(1)->phi(),
                                                       mumass_);
      const auto lvecG  = math::PtEtaPhiMLorentzVector(promptPhotons.front()->pt(),
                                                       promptPhotons.front()->eta(),
                                                       promptPhotons.front()->phi(),
                                                       0.);

      const double massMMG = ( lvecM1 + lvecM2 + lvecG ).M();

      const float dR1 = reco::deltaR( hardprocMuons.at(0)->eta(),hardprocMuons.at(0)->phi(),
                                      promptPhotons.front()->eta(),promptPhotons.front()->phi() );
      const float dR2 = reco::deltaR( hardprocMuons.at(1)->eta(),hardprocMuons.at(1)->phi(),
                                      promptPhotons.front()->eta(),promptPhotons.front()->phi() );

      histo1d_["GEN_leadMuPt"]->Fill(hardprocMuons.at(0)->pt(), aWeight);
      histo1d_["GEN_subleadMuPt"]->Fill(hardprocMuons.at(1)->pt(), aWeight);
      histo1d_["GEN_gammaPt"]->Fill(promptPhotons.front()->pt(), aWeight);
      histo1d_["GEN_dR"]->Fill(std::min(dR1,dR2), aWeight);
      histo1d_["GEN_invM"]->Fill(massMMG, aWeight);
    }

    edm::Handle<double> theprefweight;
    iEvent.getByToken(prefweight_token, theprefweight);
    double prefiringweight = *theprefweight;

    aWeight = prefiringweight*mcweight/std::abs(mcweight);

    edm::Handle<edm::View<PileupSummaryInfo>> pusummary;
    iEvent.getByToken(pileupToken_, pusummary);

    for (unsigned int idx = 0; idx < pusummary->size(); ++idx) {
      const auto& apu = pusummary->refAt(idx);

      int bx = apu->getBunchCrossing();

      if (bx==0) { // in-time PU only
        auto npu = apu->getTrueNumInteractions();
        aWeight *= purwgt_->GetBinContent( purwgt_->FindBin(apu->getTrueNumInteractions()) );
        histo1d_["PUsummary"]->Fill( static_cast<float>(npu)+0.5, aWeight );

        break;
      }
    }
  }

  histo1d_["cutflow"]->Fill(1.5,aWeight);

  reco::VertexRef primaryVertex;

  if (pvHandle.isValid() && !pvHandle->empty())
    primaryVertex = pvHandle->refAt(0).castTo<reco::VertexRef>();
  else
    return;

  edm::Handle<edm::TriggerResults> trigResultHandle;
  iEvent.getByToken(triggerToken_,trigResultHandle);

  edm::Handle<edm::View<pat::TriggerObjectStandAlone>> trigObjHandle;
  iEvent.getByToken(triggerobjectsToken_, trigObjHandle);

  const unsigned int nTrig = trigResultHandle.product()->size();
  std::vector<std::pair<std::string, int>> indices;
  edm::TriggerNames trigList = iEvent.triggerNames(*trigResultHandle);

  bool isFired = false;

  for (unsigned int iTrig = 0; iTrig != nTrig; iTrig++) {
    std::string trigName = trigList.triggerName(iTrig);
    for (unsigned int jTrig = 0; jTrig != trigList_.size(); jTrig++) {
      if (trigName.find(trigList_.at(jTrig).substr(0, trigList_.at(jTrig).find("*"))) != std::string::npos) {
        if (trigResultHandle.product()->accept(iTrig))
          isFired = true;
      }
    } // wanted triggers
  } // fired triggers

  if (!isFired)
    return;

  histo1d_["cutflow"]->Fill(2.5,aWeight);

  std::vector<edm::RefToBase<pat::TriggerObjectStandAlone>> trigObjs;

  for (unsigned iTrig = 0; iTrig < trigObjHandle->size(); iTrig++) {
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

  edm::Handle<edm::View<pat::Muon>> muonHandle;
  iEvent.getByToken(srcMuon_, muonHandle);

  if (muonHandle->empty())
    return;

  std::vector<pat::MuonRef> mediumMuons;

  for (unsigned iMu = 0; iMu < muonHandle->size(); iMu++) {
    const auto& aMu = muonHandle->refAt(iMu);

    if (aMu->innerTrack().isNull())
      continue;

    if ( aMu->pt() < ptThres_ || std::abs(aMu->eta()) > 2.4 )
      continue;

    if (!muon::isMediumMuon(*aMu) )
      continue;

    mediumMuons.push_back(aMu.castTo<pat::MuonRef>());
  }

  bool matched = false;
  std::vector<pat::MuonRef> trigMuons;

  for (const auto& trigObj : trigObjs) {
    for (const auto& aMu : mediumMuons) {
      if ( aMu->pt() < ptThresTrig_ || aMu->trackIso()/aMu->innerTrack()->pt() > 0.1 )
        continue;

      if ( reco::deltaR2(trigObj->eta(),trigObj->phi(),aMu->eta(),aMu->phi()) < 0.01 ) {
        trigMuons.push_back(aMu);
        matched = true;
      }
    }
  }

  if (!matched)
    return;

  histo1d_["cutflow"]->Fill(3.5,aWeight);
  histo1d_["nPV"]->Fill( static_cast<float>(pvHandle->size())+0.5, aWeight );

  edm::Handle<edm::View<pat::Electron>> eleHandle;
  iEvent.getByToken(srcEle_, eleHandle);

  edm::Handle<edm::View<pat::Photon>> phoHandle;
  iEvent.getByToken(srcPhoton_, phoHandle);

  edm::Handle<edm::View<reco::Conversion>> conversionHandle;
  iEvent.getByToken(conversionsToken_,conversionHandle);

  edm::Handle<edm::View<reco::Conversion>> singleLegConversionHandle;
  iEvent.getByToken(singleLegConversionsToken_,singleLegConversionHandle);

  edm::Handle<reco::BeamSpot> beamSpotHandle;
  iEvent.getByToken(beamspotToken_, beamSpotHandle);

  edm::Handle<edm::ValueMap<reco::GsfTrackRef>> addGsfTrkHandle;
  iEvent.getByToken(addGsfTrkToken_, addGsfTrkHandle);

  edm::Handle<edm::ValueMap<pat::PackedCandidateRef>> addPackedCandHandle;
  iEvent.getByToken(addPackedCandToken_, addPackedCandHandle);

  edm::Handle<edm::ValueMap<float>> trkIsoMapHandle;
  iEvent.getByToken(trkIsoMapToken_, trkIsoMapHandle);

  edm::Handle<edm::ValueMap<float>> alphaTrackHandle;
  iEvent.getByToken(alphaTrackToken_, alphaTrackHandle);

  edm::Handle<edm::ValueMap<float>> alphaCaloHandle;
  iEvent.getByToken(alphaCaloToken_, alphaCaloHandle);

  edm::Handle<edm::ValueMap<float>> dPerpInHandle;
  iEvent.getByToken(dPerpInToken_, dPerpInHandle);

  edm::Handle<edm::ValueMap<float>> normDParaInHandle;
  iEvent.getByToken(normDParaInToken_, normDParaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dEtaInHandle;
  iEvent.getByToken(union5x5dEtaInToken_, union5x5dEtaInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5dPhiInHandle;
  iEvent.getByToken(union5x5dPhiInToken_, union5x5dPhiInHandle);

  edm::Handle<edm::ValueMap<float>> union5x5EnergyHandle;
  iEvent.getByToken(union5x5EnergyToken_, union5x5EnergyHandle);

  if ( eleHandle->empty() || phoHandle->empty() )
    return;

  for (unsigned iMu = 0; iMu < trigMuons.size(); iMu++) {
    const auto& trigMuon = trigMuons.at(iMu);

    histo1d_["cutflow_tnpPair"]->Fill(0.5,aWeight);

    for (const auto& aMu : mediumMuons) {
      // prevent double counting
      bool counted = false;

      for (unsigned idx = 0; idx <= iMu; idx++) {
        if (aMu==trigMuons.at(idx)) {
          counted = true;
          break;
        }
      }

      if (counted)
        continue;

      histo1d_["cutflow_tnpPair"]->Fill(1.5,aWeight);

      for (unsigned iPho = 0; iPho < phoHandle->size(); iPho++) {
        const auto& aPho = phoHandle->refAt(iPho);

        if ( std::abs(aPho->superCluster()->eta()) > 1.4442 )
          continue;

        histo1d_["cutflow_tnpPair"]->Fill(2.5,aWeight);

        for (unsigned iEle = 0; iEle < eleHandle->size(); iEle++) {
          const auto& aEle = eleHandle->refAt(iEle);

          if ( !aEle->ecalDriven() || std::abs(aEle->superCluster()->eta()) > 1.4442 )
            continue;

          const float eta1stGSF = -(aEle->deltaEtaSeedClusterTrackAtVtx() - aEle->superCluster()->seed()->eta());
          const float u5x5Eta = (*union5x5dEtaInHandle)[aEle] + eta1stGSF;
          const float u5x5Et = (*union5x5EnergyHandle)[aEle]/std::cosh(u5x5Eta);

          if ( u5x5Et < 10. )
            continue;

          histo1d_["cutflow_tnpPair"]->Fill(3.5,aWeight);

          if ( aPho->superCluster()!=aEle->superCluster() )
            continue;

          histo1d_["cutflow_tnpPair"]->Fill(4.5,aWeight);

          const auto lvecMM = trigMuon->p4() + aMu->p4();
          const auto lvecMMG = lvecMM + aPho->p4();
          const double massMMG = lvecMMG.M();
          const double massMM = lvecMM.M();
          const double massSum = massMMG + massMM;

          const float dR1 = reco::deltaR( trigMuon->eta(),trigMuon->phi(),aPho->eta(),aPho->phi() );
          const float dR2 = reco::deltaR( aMu->eta(),aMu->phi(),aPho->eta(),aPho->phi() );

          if ( std::min(dR1,dR2) < drThres_ )
            continue;

          bool fsrCR = ( massMMG > 70. && massMMG < 110. && dR2 < 0.8 && massSum < 180. );
          bool isrCR = ( massMM > 70. && massMM < 110. );

          histo1d_["cutflow_tnpPair"]->Fill(5.5,aWeight);

          bool passMergedElectronID = aEle->electronID("mvaMergedElectron");
          const float mvascore = aEle->userFloat("mvaMergedElectronValues");

          if ( aEle->userInt("mvaMergedElectronCategories")==0 )
            histo1d_["mva_HasTrkEB"]->Fill(mvascore, aWeight);
          else
            histo1d_["mva_NoneEt2EB"]->Fill(mvascore, aWeight);

          // int32_t bitmap = aEle->userInt("modifiedHeepElectronID");
          // int32_t mask = 0x00000780; // = 0111 1000 0000 - 7th for trk iso, 8th for EM+HadD1 iso, 9th for dxy, 10 for missing inner hits
          // int32_t pass = bitmap | mask;
          // bool passMaskedId = pass==0x00000FFF; // HEEP ID has 12 cuts

          const auto& orgGsfTrk = aEle->gsfTrack();
          const auto& addGsfTrk = (*addGsfTrkHandle)[aEle];
          const auto& addPackedCand = (*addPackedCandHandle)[aEle];

          const float alphaTrack = (*alphaTrackHandle)[aEle];
          const float alphaCalo = (*alphaCaloHandle)[aEle];
          float trackIso = (*trkIsoMapHandle)[aEle];

          auto fillProbe = [&,this] (const std::string& prefix, const reco::Track* addTrkPtr, const reco::ConversionRef& theConv, float iso) {
            const auto lvecE1 = math::PtEtaPhiMLorentzVector( orgGsfTrk->pt(), orgGsfTrk->eta(), orgGsfTrk->phi(), elmass_ );
            const auto lvecEE = addTrkPtr ? lvecE1 + math::PtEtaPhiMLorentzVector( addTrkPtr->pt(), addTrkPtr->eta(), addTrkPtr->phi(), elmass_ ) : math::PtEtaPhiMLorentzVector(0.,0.,0.,0.);

            histo1d_[prefix+"_probeEle_MMG_invM"]->Fill( massMMG, aWeight );
            histo1d_[prefix+"_probeEle_MMG_rap"]->Fill( lvecMMG.Rapidity(), aWeight );
            histo1d_[prefix+"_probeEle_MM_invM"]->Fill( massMM, aWeight );
            histo1d_[prefix+"_probeEle_mu1_pt"]->Fill( trigMuon->pt(), aWeight );
            histo1d_[prefix+"_probeEle_mu1_eta"]->Fill( trigMuon->eta(), aWeight );
            histo1d_[prefix+"_probeEle_mu2_pt"]->Fill( aMu->pt(), aWeight );
            histo1d_[prefix+"_probeEle_mu2_eta"]->Fill( aMu->eta(), aWeight );
            histo1d_[prefix+"_probeEle_probe_phoEt"]->Fill( aPho->et(), aWeight );
            histo1d_[prefix+"_probeEle_probe_eleEt"]->Fill( aEle->et(), aWeight );
            histo1d_[prefix+"_probeEle_probe_eleU5x5Et"]->Fill( u5x5Et, aWeight );
            histo1d_[prefix+"_probeEle_probe_eta"]->Fill( aPho->superCluster()->eta(), aWeight );
            histo1d_[prefix+"_probeEle_mu1_dR"]->Fill( dR1, aWeight );
            histo1d_[prefix+"_probeEle_mu2_dR"]->Fill( dR2, aWeight );

            histo2d_[prefix+"_probeEle_MMG_MM_invM"]->Fill( massMMG, massMM, aWeight );

            histo1d_[prefix+"_probeEle_diel_dEta"]->Fill( addTrkPtr ? addTrkPtr->eta() - orgGsfTrk->eta() : std::numeric_limits<float>::max(), aWeight );
            histo1d_[prefix+"_probeEle_diel_dPhi"]->Fill( addTrkPtr ? reco::deltaPhi( addTrkPtr->phi(), orgGsfTrk->phi() ) : std::numeric_limits<float>::max(), aWeight );
            histo1d_[prefix+"_probeEle_diel_dR"]->Fill( addTrkPtr ? reco::deltaR( orgGsfTrk->eta(), orgGsfTrk->phi(), addTrkPtr->eta(), addTrkPtr->phi() ) : -1., aWeight );
            histo1d_[prefix+"_probeEle_diel_alphaCalo"]->Fill(alphaCalo,aWeight);
            histo1d_[prefix+"_probeEle_diel_alphaTrk"]->Fill(alphaTrack,aWeight);
            histo1d_[prefix+"_probeEle_diel_normDParaIn"]->Fill((*normDParaInHandle)[aEle],aWeight);
            histo1d_[prefix+"_probeEle_diel_HoE"]->Fill( aEle->hcalOverEcal(), aWeight );
            histo1d_[prefix+"_probeEle_diel_sigIeIe"]->Fill( aEle->full5x5_sigmaIetaIeta(), aWeight );
            histo1d_[prefix+"_probeEle_diel_dEtaIn"]->Fill( aEle->deltaEtaSuperClusterTrackAtVtx(), aWeight );
            histo1d_[prefix+"_probeEle_diel_dEtaInSeed"]->Fill( aEle->deltaEtaSuperClusterTrackAtVtx()
                                                                - aEle->superCluster()->eta()
                                                                + aEle->superCluster()->seed()->eta(), aWeight );
            histo1d_[prefix+"_probeEle_diel_dPhiIn"]->Fill( aEle->deltaPhiSuperClusterTrackAtVtx(), aWeight );
            histo1d_[prefix+"_probeEle_diel_R9"]->Fill( aEle->full5x5_r9(), aWeight );
            histo1d_[prefix+"_probeEle_diel_EoverP"]->Fill( aEle->eSuperClusterOverP(), aWeight );
            histo1d_[prefix+"_probeEle_diel_trackIso"]->Fill(iso, aWeight);
            histo1d_[prefix+"_probeEle_diel_dPerpIn"]->Fill((*dPerpInHandle)[aEle], aWeight);
            histo1d_[prefix+"_probeEle_diel_u5x5dEtaIn"]->Fill((*union5x5dEtaInHandle)[aEle], aWeight);
            histo1d_[prefix+"_probeEle_diel_u5x5dPhiIn"]->Fill((*union5x5dPhiInHandle)[aEle], aWeight);

            histo1d_[prefix+"_probeEle_diel_invM_ee"]->Fill( addTrkPtr ? lvecEE.M() : -1., aWeight );
            histo1d_[prefix+"_probeEle_diel_pt"]->Fill( addTrkPtr ? lvecEE.pt() : -1., aWeight );
            histo1d_[prefix+"_probeEle_diel_etSC"]->Fill( aEle->superCluster()->energy()*
                                                          aEle->superCluster()->position().rho()/
                                                          aEle->superCluster()->position().r(), aWeight );
            histo1d_[prefix+"_probeEle_1st_pt"]->Fill( orgGsfTrk->pt(), aWeight );
            histo1d_[prefix+"_probeEle_2nd_pt"]->Fill( addTrkPtr ? addTrkPtr->pt() : -1., aWeight );

            if ( aEle->userInt("mvaMergedElectronCategories")==0 )
              histo1d_[prefix+"_probeEle_mva_HasTrkEB"]->Fill(mvascore, aWeight);
            else
              histo1d_[prefix+"_probeEle_mva_NoneEt2EB"]->Fill(mvascore, aWeight);

            histo1d_[prefix+"_probeEle_conv_invM"]->Fill( theConv->pairInvariantMass(), aWeight );
            histo1d_[prefix+"_probeEle_conv_dPhiVtx"]->Fill( theConv->dPhiTracksAtVtx(), aWeight );
            histo1d_[prefix+"_probeEle_conv_lxy"]->Fill( theConv->lxy( beamSpotHandle->position() ), aWeight );
            histo1d_[prefix+"_probeEle_conv_nHitsBeforeVtx_1st"]->Fill( !theConv->nHitsBeforeVtx().empty() ? static_cast<float>(theConv->nHitsBeforeVtx().front())+0.5 : -0.5, aWeight );
            histo1d_[prefix+"_probeEle_conv_nHitsBeforeVtx_2nd"]->Fill( theConv->nHitsBeforeVtx().size() > 1 ? static_cast<float>(theConv->nHitsBeforeVtx().at(1))+0.5 : -0.5, aWeight );
            histo1d_[prefix+"_probeEle_conv_nSharedHits"]->Fill( static_cast<float>(theConv->nSharedHits())+0.5, aWeight );
          };

          auto fillProbes = [&,this] (const std::string& prefix, const reco::Track* addTrkPtr, const reco::ConversionRef& theConv, float iso) {
            if (fsrCR) {
              fillProbe(prefix+"MB_FSR",addTrkPtr,theConv,iso);

              if (passMergedElectronID)
                fillProbe(prefix+"PassMergedEle",addTrkPtr,theConv,iso);
              else
                fillProbe(prefix+"FailMergedEle",addTrkPtr,theConv,iso);
            }

            if (isrCR)
              fillProbe(prefix+"MB_ISR",addTrkPtr,theConv,iso);
          };

          if (addGsfTrk==orgGsfTrk && addPackedCand.isNull()) {
            bool isSingleLeg = false;
            auto conv = reco::ConversionRef();

            histo1d_["cutflow_tnpPair"]->Fill(6.5,aWeight);

            for (unsigned iConv = 0; iConv < singleLegConversionHandle->size(); iConv++) {
              const auto aConv = singleLegConversionHandle->refAt(iConv);

              if ( ConversionTools::matchesConversion(*aEle,*aConv) ) {
                isSingleLeg = true;
                conv = aConv.castTo<reco::ConversionRef>();
                break;
              }
            }

            // emulate track iso criteria
            if ( reco::deltaR2(orgGsfTrk->eta(),orgGsfTrk->phi(),aMu->innerTrack()->eta(),aMu->innerTrack()->phi()) < 0.09 ) {
              if ( std::abs( aMu->innerTrack()->eta() - orgGsfTrk->eta() ) > 0.005 ) {
                if ( std::abs( aMu->innerTrack()->vz() - orgGsfTrk->vz() ) < 0.1 &&
                     aMu->innerTrack()->hitPattern().numberOfValidHits() >= 8 &&
                     aMu->innerTrack()->hitPattern().numberOfValidPixelHits() >=1 &&
                     aMu->innerTrack()->ptError()/aMu->innerTrack()->pt() < 0.1 )
                  trackIso -= aMu->innerTrack()->pt();
              }
            }

            if (isSingleLeg) {
              histo1d_["cutflow_tnpPair"]->Fill(7.5,aWeight);

              fillProbes("singleLeg_",nullptr,conv,trackIso);

              if (fsrCR) {
                invM_ = massMMG;
                u5x5Et_ = u5x5Et;
                wgt_ = aWeight;
                passME_ = static_cast<int>(passMergedElectronID);
                tree_->Fill();
              }

              if (u5x5Et > 10. && u5x5Et < 15.)
                fillProbes("singleLeg_Et10to15_",nullptr,conv,trackIso);
              else if (u5x5Et > 15. && u5x5Et < 20.)
                fillProbes("singleLeg_Et15to20_",nullptr,conv,trackIso);
              else if (u5x5Et > 20. && u5x5Et < 30.)
                fillProbes("singleLeg_Et20to30_",nullptr,conv,trackIso);
              else if (u5x5Et > 30.)
                fillProbes("singleLeg_Et30toInf_",nullptr,conv,trackIso);
            }
          } else {
            // first look for GSF and then packedPFcand
            const reco::Track* addTrk = addGsfTrk.get();

            if ( addGsfTrk==orgGsfTrk && addPackedCand.isNonnull() )
              addTrk = addPackedCand->bestTrack();

            if ( orgGsfTrk->charge()*addTrk->charge() > 0 )
              continue;

            histo1d_["cutflow_tnpPair"]->Fill(8.5,aWeight);

            bool isDoubleLeg = false;
            auto conv = reco::ConversionRef();

            for (unsigned iConv = 0; iConv < conversionHandle->size(); iConv++) {
              const auto aConv = conversionHandle->refAt(iConv);

              if ( ConversionTools::matchesConversion(*aEle,*aConv) &&
                   ConversionTools::isGoodConversion(*aConv,beamSpotHandle->position()) ) {
                isDoubleLeg = true;
                conv = aConv.castTo<reco::ConversionRef>();
                break;
              }
            }

            // emulate track iso criteria
            if ( reco::deltaR2(orgGsfTrk->eta(),orgGsfTrk->phi(),aMu->innerTrack()->eta(),aMu->innerTrack()->phi()) < 0.09 ) {
              if ( std::abs( aMu->innerTrack()->eta() - orgGsfTrk->eta() ) > 0.005 ) {
                if ( !( std::abs( aMu->innerTrack()->eta() - addTrk->eta() ) < 0.005 &&
                        reco::deltaPhi( aMu->innerTrack()->phi(), addTrk->phi() ) < 0.05 ) ) {
                  if ( std::abs( aMu->innerTrack()->vz() - orgGsfTrk->vz() ) < 0.1 &&
                       aMu->innerTrack()->hitPattern().numberOfValidHits() >= 8 &&
                       aMu->innerTrack()->hitPattern().numberOfValidPixelHits() >=1 &&
                       aMu->innerTrack()->ptError()/aMu->innerTrack()->pt() < 0.1 )
                    trackIso -= aMu->innerTrack()->pt();
                }
              }
            }

            if (isDoubleLeg /*&& conv->nSharedHits()==0*/) {
              histo1d_["cutflow_tnpPair"]->Fill(9.5,aWeight);

              fillProbes("",addTrk,conv,trackIso);
            }
          } // HasTrk
        } // eleHandle
      } // phoHandle
    } // mediumMuons
  } // trigMuons
}

DEFINE_FWK_MODULE(MergedLeptonIDConversionAnalyzer);
