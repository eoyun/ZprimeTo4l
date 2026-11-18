// Shared core for the HEEP vs modified-HEEP efficiency study (vs dR_gen, EB/EE).
// Single source of truth for the selection + dR_gen binning, used by the condor
// worker (fillHeepEffi), the plotter (plotHeepEffi) and the single-file demo (runHeepEffi).
//
// Gen A->ee pairing: greedy nearest opposite-charge dR.
// Gen pair select : both |eta_gen|<2.5, leading gen pT>20 GeV.
// Reco            : pT_reco>10, dR(reco,gen)<0.001 to a gen e of a selected pair,
//                   closest-gen deterministic match, each reco counted once.
// EB |etaSC|<1.4442 ; EE 1.566<|etaSC|<2.5 ; transition excluded.
// Three IDs drawn: standard HEEP V70 (electron_passHeepV70 / heepV70Bitmap),
//   modified HEEP (electron_passModHeep / modHeepBitmap),
//   modified HEEP Run3 (electron_passModifiedHEEPRun3ID / bitmapModifiedHEEPRun3ID).
// Requires schema-v6 ntuples (26.10.08 production) that store ModifiedHEEPRun3ID.
#ifndef HEEPEFFICORE_H
#define HEEPEFFICORE_H
#include <TTree.h>
#include <TH1D.h>
#include <TH2D.h>
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>

namespace heep {

// dR_gen binning — filler and plotter MUST share these.
// Variable width for a log-x view: 1e-4..1e-3 step 1e-4, 1e-3..1e-2 step 1e-3,
// 1e-2..1e-1 step 1e-2, 1e-1..0.5 step 0.1. (entries < 1e-4 -> underflow)
static const double XMIN = 1e-4, XMAX = 1.0;   // log-x axis view range
inline std::vector<double> binEdges(){
  std::vector<double> e;
  for(int i=1;i<=9;i++)  e.push_back(1e-4*i);  // 0.0001 .. 0.0009
  for(int i=1;i<=9;i++)  e.push_back(1e-3*i);  // 0.001  .. 0.009
  for(int i=1;i<=9;i++)  e.push_back(1e-2*i);  // 0.01   .. 0.09
  for(int i=1;i<=10;i++) e.push_back(0.1*i);   // 0.1,0.2,...,1.0
  return e;
}

static const double ISOMAX = 50.;   // isolation histogram range [0, ISOMAX) GeV

// reco-pT binning (GeV) for efficiency vs pT: fine across the HEEP turn-on, coarse at high pT.
static const double PTMIN = 10., PTMAX = 1000.;   // log-x axis view range
inline std::vector<double> ptEdges(){
  std::vector<double> e;
  for(int x=10;x<50;  x+=5)  e.push_back(x);    // 10,15,...,45
  for(int x=50;x<100; x+=10) e.push_back(x);    // 50,60,...,90
  for(int x=100;x<200;x+=25) e.push_back(x);    // 100,125,150,175
  double hi[]={200,250,300,400,500,700,1000};
  for(double v:hi) e.push_back(v);
  return e;
}

struct Cutflow { long nEvt=0,nPair=0,nPairSel=0,nRecoPt=0,nDenEB=0,nDenEE=0; double drmin=1e9,drmax=-1e9; };
struct Hists   { TH1D *denEB=0,*befEB=0,*aftEB=0,*denEE=0,*befEE=0,*aftEE=0;             // vs dR_gen (bef=V70, aft=modHEEP)
                 TH1D *denEBpt=0,*befEBpt=0,*aftEBpt=0,*denEEpt=0,*befEEpt=0,*aftEEpt=0; // vs reco pT
                 TH1D *r3EB=0,*r3EE=0,*r3EBpt=0,*r3EEpt=0;   // numerator: modified HEEP Run3 (vs dR, vs pT)
                 TH1D *cfEB=0,*cfEE=0,        // modified-HEEP sequential cutflow (13 stages), denominator electrons
                      *cfV70EB=0,*cfV70EE=0,  // standard HEEP V70 sequential cutflow
                      *cfR3EB=0,*cfR3EE=0;    // modified HEEP Run3 sequential cutflow
                 // same cutflows restricted to the dR_gen > 0.1 subset (large-dR drop investigation)
                 TH1D *cfHiEB=0,*cfHiEE=0,*cfHiV70EB=0,*cfHiV70EE=0,*cfHiR3EB=0,*cfHiR3EE=0;
                 // per-cut marginal pass numerators: [cut bit 0..11][id 0=v70 1=mod 2=run3][region 0=EB 1=EE]
                 TH1D *cut[12][3][2]={}, *cutpt[12][3][2]={};    // vs dR_gen and vs reco pT
                 // isolation value distributions (denominator electrons; -999/invalid skipped)
                 TH1D *isoTrk[2][2]={};    // [0=modified(eleTrkPtIso) 1=std dr03TkSumPt][region]
                 TH1D *isoEcal[3][2]={};   // [0=modified 1=modified Run3 2=std dr03EcalRecHitSumEt][region]
                 TH1D *isoHcal[2]={};      // std dr03HcalTowerSumEt [region]
                 TH2D *drTrk[2][2]={}; };  // dR_gen vs track iso [0=mod 1=std][region]

// modified-HEEP cutflow stage labels (bin1=denominator, bin i+2 = cumulative pass of cuts 0..i)
static const char* kCutLabels[13] = {"denominator","ET>20","SCEta","dEtaInSeed","dPhiIn",
  "#sigma_{i#etai#eta}","E2x5/E5x5","H/E","trkIso","EcalHcalIso","dxy","missHits","ecalDriven"};
// filename-safe keys for the 12 cuts (bit 0..11), matching kCutLabels[bit+1]
static const char* kCutKeys[12] = {"et","sceta","detainseed","dphiin","sigmaieta","e2x5",
  "hoe","trkiso","caloiso","dxy","misshits","ecaldriven"};

inline double dphi(double a,double b){ double d=a-b; while(d> M_PI)d-=2*M_PI; while(d<-M_PI)d+=2*M_PI; return d; }
inline double dR(double e1,double p1,double e2,double p2){ double de=e1-e2, dp=dphi(p1,p2); return std::sqrt(de*de+dp*dp); }

inline Hists makeHists(){
  std::vector<double> e = binEdges(), p = ptEdges();
  auto mkdr=[&](const std::string& n){ TH1D* h=new TH1D(n.c_str(),";#DeltaR_{gen};electrons",(int)e.size()-1,e.data()); h->Sumw2(); return h; };
  auto mkpt=[&](const std::string& n){ TH1D* h=new TH1D(n.c_str(),";p_{T} [GeV];electrons",(int)p.size()-1,p.data()); h->Sumw2(); return h; };
  Hists h;
  h.denEB=mkdr("denominator_EB"); h.befEB=mkdr("numerator_before_EB"); h.aftEB=mkdr("numerator_after_EB");
  h.denEE=mkdr("denominator_EE"); h.befEE=mkdr("numerator_before_EE"); h.aftEE=mkdr("numerator_after_EE");
  h.denEBpt=mkpt("denominator_EB_pt"); h.befEBpt=mkpt("numerator_before_EB_pt"); h.aftEBpt=mkpt("numerator_after_EB_pt");
  h.denEEpt=mkpt("denominator_EE_pt"); h.befEEpt=mkpt("numerator_before_EE_pt"); h.aftEEpt=mkpt("numerator_after_EE_pt");
  h.r3EB=mkdr("numerator_run3_EB");       h.r3EE=mkdr("numerator_run3_EE");
  h.r3EBpt=mkpt("numerator_run3_EB_pt");  h.r3EEpt=mkpt("numerator_run3_EE_pt");
  auto mkcf=[&](const std::string& n, const char* etLabel){ TH1D* h=new TH1D(n.c_str(),";;electrons",13,0,13); h->Sumw2();
    for(int b=0;b<13;b++) h->GetXaxis()->SetBinLabel(b+1,kCutLabels[b]);
    h->GetXaxis()->SetBinLabel(2,etLabel); return h; };
  h.cfEB=mkcf("cutflow_EB","ET>20");       h.cfEE=mkcf("cutflow_EE","ET>20");
  h.cfV70EB=mkcf("cutflow_v70_EB","ET>35"); h.cfV70EE=mkcf("cutflow_v70_EE","ET>35");
  h.cfR3EB=mkcf("cutflow_r3_EB","ET>20");   h.cfR3EE=mkcf("cutflow_r3_EE","ET>20");
  h.cfHiEB=mkcf("cutflow_hidr_EB","ET>20");       h.cfHiEE=mkcf("cutflow_hidr_EE","ET>20");
  h.cfHiV70EB=mkcf("cutflow_v70_hidr_EB","ET>35"); h.cfHiV70EE=mkcf("cutflow_v70_hidr_EE","ET>35");
  h.cfHiR3EB=mkcf("cutflow_r3_hidr_EB","ET>20");    h.cfHiR3EE=mkcf("cutflow_r3_hidr_EE","ET>20");
  const char* idn[3]={"v70","mod","run3"}; const char* regn[2]={"EB","EE"};
  for(int c=0;c<12;c++) for(int d=0;d<3;d++) for(int rr=0;rr<2;rr++){
    std::string base=std::string(kCutKeys[c])+"_"+idn[d]+"_"+regn[rr];
    h.cut[c][d][rr]=mkdr(base); h.cutpt[c][d][rr]=mkpt(base+"_pt"); }
  // isolation distributions: 0-50 GeV (last bin includes overflow)
  auto mkiso=[&](const std::string& n,const char* xt){ TH1D* x=new TH1D(n.c_str(),(std::string(";")+xt+";electrons").c_str(),100,0,ISOMAX); x->Sumw2(); return x; };
  const char* trkn[2]={"mod","std"}; const char* ecn[3]={"mod","run3","std"};
  for(int rr=0;rr<2;rr++){
    for(int v=0;v<2;v++){
      h.isoTrk[v][rr]=mkiso(std::string("iso_trk_")+trkn[v]+"_"+regn[rr],"track iso [GeV]");
      h.drTrk[v][rr]=new TH2D((std::string("dr_vs_trkiso_")+trkn[v]+"_"+regn[rr]).c_str(),";#DeltaR_{gen};track iso [GeV]",(int)e.size()-1,e.data(),50,0,ISOMAX);
      h.drTrk[v][rr]->Sumw2(); }
    for(int v=0;v<3;v++) h.isoEcal[v][rr]=mkiso(std::string("iso_ecal_")+ecn[v]+"_"+regn[rr],"ECAL iso [GeV]");
    h.isoHcal[rr]=mkiso(std::string("iso_hcal_std_")+regn[rr],"HCAL iso [GeV]");
  }
  return h;
}

inline void processTree(TTree* t, Hists& h, Cutflow& cf, bool require2ndTrk=false){
  std::vector<float> *g_pt=0,*g_eta=0,*g_phi=0; std::vector<int>* g_pid=0;
  std::vector<float> *e_pt=0,*e_eta=0,*e_phi=0,*e_etaSC=0,*e_perp=0;
  std::vector<int>   *e_heep=0,*e_modheep=0,*e_modbm=0,*e_heepbm=0,*e_r3=0,*e_r3bm=0;
  t->SetBranchAddress("gen_pt",&g_pt);   t->SetBranchAddress("gen_eta",&g_eta);
  t->SetBranchAddress("gen_phi",&g_phi); t->SetBranchAddress("gen_pdgId",&g_pid);
  t->SetBranchAddress("electron_ptRaw",&e_pt);   t->SetBranchAddress("electron_etaRaw",&e_eta);
  t->SetBranchAddress("electron_phiRaw",&e_phi); t->SetBranchAddress("electron_etaSC",&e_etaSC);
  t->SetBranchAddress("electron_passHeepV70",&e_heep);
  t->SetBranchAddress("electron_passModHeep",&e_modheep);
  t->SetBranchAddress("electron_modHeepBitmap",&e_modbm);
  t->SetBranchAddress("electron_heepV70Bitmap",&e_heepbm);
  t->SetBranchAddress("electron_passModifiedHEEPRun3ID",&e_r3);
  t->SetBranchAddress("electron_bitmapModifiedHEEPRun3ID",&e_r3bm);
  t->SetBranchAddress("electron_dPerpIn",&e_perp);   // valid => modified HEEP found a 2nd (additional) track
  // isolation values; Run3 calo iso only exists in fixed (schema>=7) ntuples
  std::vector<float> *i_trkMod=0,*i_trkStd=0,*i_ecMod=0,*i_ecR3=0,*i_ecStd=0,*i_hcStd=0;
  t->SetBranchAddress("electron_modTrkIso",&i_trkMod);
  t->SetBranchAddress("electron_dr03TkSumPt",&i_trkStd);
  t->SetBranchAddress("electron_modEcalHcalIso",&i_ecMod);
  t->SetBranchAddress("electron_dr03EcalRecHitSumEt",&i_ecStd);
  t->SetBranchAddress("electron_dr03HcalTowerSumEt",&i_hcStd);
  const bool hasEcR3 = t->GetBranch("electron_modEcalHcalIsoRun3")!=nullptr;
  if(hasEcR3) t->SetBranchAddress("electron_modEcalHcalIsoRun3",&i_ecR3);
  auto ok=[](float v){ return std::isfinite(v) && v!=-999.f && std::fabs(v)<1e30; };
  auto cl=[](double v){ return v<ISOMAX ? v : ISOMAX-1e-6; };   // fold overflow into last bin

  Long64_t N=t->GetEntries();
  for(Long64_t i=0;i<N;i++){
    t->GetEntry(i); cf.nEvt++;
    std::vector<int> pos, neg;                       // pdg -11 = e+, 11 = e-
    for(size_t k=0;k<g_pid->size();k++){ if((*g_pid)[k]==-11) pos.push_back((int)k); else if((*g_pid)[k]==11) neg.push_back((int)k); }
    std::vector<bool> up(pos.size(),false), un(neg.size(),false);
    std::vector<double> genPairDR(g_pid->size(), -1.0);   // genIdx -> its selected pair's dR_gen (-1 if none)
    while(true){
      double best=1e9; int bp=-1,bn=-1;
      for(size_t a=0;a<pos.size();a++){ if(up[a])continue;
        for(size_t b=0;b<neg.size();b++){ if(un[b])continue;
          double d=dR((*g_eta)[pos[a]],(*g_phi)[pos[a]],(*g_eta)[neg[b]],(*g_phi)[neg[b]]);
          if(d<best){best=d;bp=(int)a;bn=(int)b;} } }
      if(bp<0) break;
      up[bp]=true; un[bn]=true; cf.nPair++;
      int ia=pos[bp], ib=neg[bn];
      bool etaOK  = std::fabs((*g_eta)[ia])<2.5 && std::fabs((*g_eta)[ib])<2.5;
      bool leadOK = std::max((*g_pt)[ia],(*g_pt)[ib])>20.0;
      if(etaOK && leadOK){ cf.nPairSel++; genPairDR[ia]=best; genPairDR[ib]=best;
        cf.drmin=std::min(cf.drmin,best); cf.drmax=std::max(cf.drmax,best); }
    }
    for(size_t r=0; r<e_pt->size(); r++){
      if((*e_pt)[r] <= 10.0) continue; cf.nRecoPt++;
      int gbest=-1; double dbest=1e9;
      for(size_t k=0;k<genPairDR.size();k++){ if(genPairDR[k]<0) continue;
        double d=dR((*e_eta)[r],(*e_phi)[r],(*g_eta)[k],(*g_phi)[k]);
        if(d<dbest){ dbest=d; gbest=(int)k; } }
      if(gbest<0 || dbest>=0.001) continue;
      double aeta=std::fabs((*e_etaSC)[r]);
      bool isEB = aeta<1.4442;
      bool isEE = aeta>1.566 && aeta<2.5;
      if(!isEB && !isEE) continue;
      if(require2ndTrk){ float pv=(*e_perp)[r]; if(!(std::isfinite(pv)&&pv!=-999.f&&std::fabs(pv)<1e30)) continue; }
      double drg = genPairDR[gbest];
      double ptv = (*e_pt)[r];
      { int ri=isEB?0:1;   // isolation value distributions for denominator electrons
        float tm=(*i_trkMod)[r], ts=(*i_trkStd)[r];
        if(ok(tm)){ h.isoTrk[0][ri]->Fill(cl(tm)); h.drTrk[0][ri]->Fill(drg,cl(tm)); }
        if(ok(ts)){ h.isoTrk[1][ri]->Fill(cl(ts)); h.drTrk[1][ri]->Fill(drg,cl(ts)); }
        if(ok((*i_ecMod)[r])) h.isoEcal[0][ri]->Fill(cl((*i_ecMod)[r]));
        if(hasEcR3 && ok((*i_ecR3)[r])) h.isoEcal[1][ri]->Fill(cl((*i_ecR3)[r]));
        if(ok((*i_ecStd)[r])) h.isoEcal[2][ri]->Fill(cl((*i_ecStd)[r]));
        if(ok((*i_hcStd)[r])) h.isoHcal[ri]->Fill(cl((*i_hcStd)[r])); }
      if(isEB){ cf.nDenEB++;
        h.denEB->Fill(drg);   if((*e_heep)[r]) h.befEB->Fill(drg);   if((*e_modheep)[r]) h.aftEB->Fill(drg);   if((*e_r3)[r]) h.r3EB->Fill(drg);
        h.denEBpt->Fill(ptv); if((*e_heep)[r]) h.befEBpt->Fill(ptv); if((*e_modheep)[r]) h.aftEBpt->Fill(ptv); if((*e_r3)[r]) h.r3EBpt->Fill(ptv); }
      else    { cf.nDenEE++;
        h.denEE->Fill(drg);   if((*e_heep)[r]) h.befEE->Fill(drg);   if((*e_modheep)[r]) h.aftEE->Fill(drg);   if((*e_r3)[r]) h.r3EE->Fill(drg);
        h.denEEpt->Fill(ptv); if((*e_heep)[r]) h.befEEpt->Fill(ptv); if((*e_modheep)[r]) h.aftEEpt->Fill(ptv); if((*e_r3)[r]) h.r3EEpt->Fill(ptv); }
      // sequential cutflow (cumulative pass of cuts 0..ci), per region, for all three IDs
      int bmod=(*e_modbm)[r], bheep=(*e_heepbm)[r], br3=(*e_r3bm)[r];
      TH1D* cfh  = isEB ? h.cfEB    : h.cfEE;
      TH1D* cfv  = isEB ? h.cfV70EB : h.cfV70EE;
      TH1D* cfr  = isEB ? h.cfR3EB  : h.cfR3EE;
      cfh->Fill(0.5); cfv->Fill(0.5); cfr->Fill(0.5);   // bin1: denominator
      for(int ci=0;ci<12;ci++){ int msk=(1<<(ci+1))-1;
        if((bmod&msk)==msk)  cfh->Fill(ci+1.5);
        if((bheep&msk)==msk) cfv->Fill(ci+1.5);
        if((br3&msk)==msk)   cfr->Fill(ci+1.5); }
      if(drg>0.1){   // same cutflow, restricted to large dR_gen (>0.1)
        TH1D* ch=isEB?h.cfHiEB:h.cfHiEE, *cv=isEB?h.cfHiV70EB:h.cfHiV70EE, *cr=isEB?h.cfHiR3EB:h.cfHiR3EE;
        ch->Fill(0.5); cv->Fill(0.5); cr->Fill(0.5);
        for(int ci=0;ci<12;ci++){ int msk=(1<<(ci+1))-1;
          if((bmod&msk)==msk)  ch->Fill(ci+1.5);
          if((bheep&msk)==msk) cv->Fill(ci+1.5);
          if((br3&msk)==msk)   cr->Fill(ci+1.5); } }
      // per-cut marginal pass (bit c set) for all 12 cuts, for v70/mod/run3
      { int ri=isEB?0:1, bits[3]={bheep,bmod,br3};
        for(int c=0;c<12;c++) for(int d=0;d<3;d++) if((bits[d]>>c)&1){ h.cut[c][d][ri]->Fill(drg); h.cutpt[c][d][ri]->Fill(ptv); } }
    }
  }
}

inline void printCutflow(const Cutflow& cf, const Hists& h, const char* tag=""){
  printf("=== cutflow %s ===\n", tag);
  printf("events                 : %ld\n", cf.nEvt);
  printf("gen pairs formed       : %ld\n", cf.nPair);
  printf("gen pairs selected     : %ld  (dR_gen min=%.4f max=%.4f)\n", cf.nPairSel, cf.drmin, cf.drmax);
  printf("reco e pT>10           : %ld\n", cf.nRecoPt);
  printf("denominator EB / EE    : %ld / %ld\n", cf.nDenEB, cf.nDenEE);
  printf("numerator before EB/EE : %.0f / %.0f\n", h.befEB->Integral(), h.befEE->Integral());
  printf("numerator after  EB/EE : %.0f / %.0f\n", h.aftEB->Integral(), h.aftEE->Integral());
}

} // namespace heep
#endif
