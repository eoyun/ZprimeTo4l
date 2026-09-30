// HEEP vs modified-HEEP electron ID efficiency vs dR_gen (EB/EE), Run-3 signal MC.
// Reads the flat ntuple (schema v5, with gen_* branches) directly with plain ROOT.
// Gen A->ee pairing: greedy nearest opposite-charge dR (per user request, first pass).
//
// Selections (per AN2023_143 Sec 4.2 baseline + stated implementation assumptions):
//   Gen pair : both |eta_gen| < 2.5, leading gen pT (max of the two) > 20 GeV.
//   Reco     : pT_reco > 10, dR(reco,gen) < 0.001 to a gen e of a selected pair,
//              closest-gen deterministic match, each reco counted once.
//   EB: |etaSC|<1.4442 ; EE: 1.566<|etaSC|<2.5 ; transition excluded.
//   Denominator filled once per matched reco e, at its pair's dR_gen (no ID).
//   Numerator before = standard HEEP V70 (electron_passHeepV70).
//   Numerator after  = modified HEEP     (electron_passModHeep).
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TLegend.h>
#include <TStyle.h>
#include <vector>
#include <string>
#include <cmath>
#include <cstdio>

static double dphi(double a, double b){ double d=a-b; while(d> M_PI)d-=2*M_PI; while(d<-M_PI)d+=2*M_PI; return d; }
static double dR(double e1,double p1,double e2,double p2){ double de=e1-e2, dp=dphi(p1,p2); return std::sqrt(de*de+dp*dp); }

int main(int argc, char** argv){
  std::string in     = (argc>1)? argv[1] : "Events.root";
  std::string outdir = (argc>2)? argv[2] : "imgs/signal";
  // dR_gen binning (configurable)
  const int    NB = 20;
  const double LO = 0.0, HI = 0.10;

  TFile* f = TFile::Open(in.c_str());
  if(!f || f->IsZombie()){ printf("cannot open %s\n", in.c_str()); return 1; }
  TTree* t = (TTree*)f->Get("resolvedNtuplizer/Events");
  if(!t){ printf("no tree resolvedNtuplizer/Events\n"); return 1; }

  std::vector<float> *g_pt=0,*g_eta=0,*g_phi=0; std::vector<int>* g_pid=0;
  std::vector<float> *e_pt=0,*e_eta=0,*e_phi=0,*e_etaSC=0;
  std::vector<int>   *e_heep=0,*e_modheep=0;
  t->SetBranchAddress("gen_pt",&g_pt);   t->SetBranchAddress("gen_eta",&g_eta);
  t->SetBranchAddress("gen_phi",&g_phi); t->SetBranchAddress("gen_pdgId",&g_pid);
  t->SetBranchAddress("electron_ptRaw",&e_pt);   t->SetBranchAddress("electron_etaRaw",&e_eta);
  t->SetBranchAddress("electron_phiRaw",&e_phi); t->SetBranchAddress("electron_etaSC",&e_etaSC);
  t->SetBranchAddress("electron_passHeepV70",&e_heep);
  t->SetBranchAddress("electron_passModHeep",&e_modheep);

  auto mk=[&](const char* n){ return new TH1D(n,";#DeltaR_{gen};electrons",NB,LO,HI); };
  TH1D* denEB=mk("denominator_EB"); TH1D* befEB=mk("numerator_before_EB"); TH1D* aftEB=mk("numerator_after_EB");
  TH1D* denEE=mk("denominator_EE"); TH1D* befEE=mk("numerator_before_EE"); TH1D* aftEE=mk("numerator_after_EE");
  for(TH1D* h : {denEB,befEB,aftEB,denEE,befEE,aftEE}) h->Sumw2();

  long nEvt=0, nPair=0, nPairSel=0, nRecoPt=0, nDenEB=0, nDenEE=0;
  double drmin=1e9, drmax=-1e9;

  Long64_t N=t->GetEntries();
  for(Long64_t i=0;i<N;i++){
    t->GetEntry(i); nEvt++;
    // gen electron charge lists (pdg 11 = e-, -11 = e+)
    std::vector<int> pos, neg;
    for(size_t k=0;k<g_pid->size();k++){ if((*g_pid)[k]==-11) pos.push_back((int)k); else if((*g_pid)[k]==11) neg.push_back((int)k); }
    // greedy nearest opposite-charge pairing -> per gen-electron pair dR_gen (only if pair selected)
    std::vector<bool> up(pos.size(),false), un(neg.size(),false);
    // genIdx -> dR_gen of its selected pair (-1 = not in a selected pair)
    std::vector<double> genPairDR(g_pid->size(), -1.0);
    while(true){
      double best=1e9; int bp=-1,bn=-1;
      for(size_t a=0;a<pos.size();a++){ if(up[a])continue;
        for(size_t b=0;b<neg.size();b++){ if(un[b])continue;
          double d=dR((*g_eta)[pos[a]],(*g_phi)[pos[a]],(*g_eta)[neg[b]],(*g_phi)[neg[b]]);
          if(d<best){best=d;bp=(int)a;bn=(int)b;} } }
      if(bp<0) break;
      up[bp]=true; un[bn]=true; nPair++;
      int ia=pos[bp], ib=neg[bn];
      // gen pair selection
      bool etaOK = std::fabs((*g_eta)[ia])<2.5 && std::fabs((*g_eta)[ib])<2.5;
      bool leadOK = std::max((*g_pt)[ia],(*g_pt)[ib])>20.0;
      if(etaOK && leadOK){ nPairSel++; genPairDR[ia]=best; genPairDR[ib]=best;
        drmin=std::min(drmin,best); drmax=std::max(drmax,best); }
    }
    // reco electrons -> denominator + numerators
    for(size_t r=0; r<e_pt->size(); r++){
      if((*e_pt)[r] <= 10.0) continue; nRecoPt++;
      // closest gen electron among selected-pair gen electrons, dR<0.001
      int gbest=-1; double dbest=1e9;
      for(size_t k=0;k<genPairDR.size();k++){ if(genPairDR[k]<0) continue;
        double d=dR((*e_eta)[r],(*e_phi)[r],(*g_eta)[k],(*g_phi)[k]);
        if(d<dbest){ dbest=d; gbest=(int)k; } }
      if(gbest<0 || dbest>=0.001) continue;
      double aeta=std::fabs((*e_etaSC)[r]);
      bool isEB = aeta<1.4442;
      bool isEE = aeta>1.566 && aeta<2.5;
      if(!isEB && !isEE) continue;
      double drg = genPairDR[gbest];
      if(isEB){ nDenEB++; denEB->Fill(drg); if((*e_heep)[r]) befEB->Fill(drg); if((*e_modheep)[r]) aftEB->Fill(drg); }
      else    { nDenEE++; denEE->Fill(drg); if((*e_heep)[r]) befEE->Fill(drg); if((*e_modheep)[r]) aftEE->Fill(drg); }
    }
  }

  // subset sanity check
  bool subsetOK=true;
  for(int b=1;b<=NB;b++){
    if(befEB->GetBinContent(b)>denEB->GetBinContent(b)+1e-6) subsetOK=false;
    if(aftEB->GetBinContent(b)>denEB->GetBinContent(b)+1e-6) subsetOK=false;
    if(befEE->GetBinContent(b)>denEE->GetBinContent(b)+1e-6) subsetOK=false;
    if(aftEE->GetBinContent(b)>denEE->GetBinContent(b)+1e-6) subsetOK=false;
  }

  printf("\n=== cutflow ===\n");
  printf("events                 : %ld\n", nEvt);
  printf("gen pairs formed       : %ld\n", nPair);
  printf("gen pairs selected     : %ld  (dR_gen min=%.4f max=%.4f)\n", nPairSel, drmin, drmax);
  printf("reco e pT>10           : %ld\n", nRecoPt);
  printf("denominator EB / EE    : %ld / %ld\n", nDenEB, nDenEE);
  printf("numerator before EB/EE : %.0f / %.0f\n", befEB->Integral(), befEE->Integral());
  printf("numerator after  EB/EE : %.0f / %.0f\n", aftEB->Integral(), aftEE->Integral());
  printf("numerator subset of denominator : %s\n", subsetOK?"YES":"NO");

  // write histograms
  TFile* fo=TFile::Open("hists_heepEffi.root","RECREATE");
  for(TH1D* h : {denEB,befEB,aftEB,denEE,befEE,aftEE}) h->Write();
  fo->Close();

  // draw efficiency (before vs after) for EB and EE
  gStyle->SetOptStat(0);
  std::string mkdircmd = "mkdir -p "+outdir; if(system(mkdircmd.c_str())){}
  TCanvas c("c","",1200,500); c.Divide(2,1);
  auto drawPad=[&](int pad,TH1D* den,TH1D* bef,TH1D* aft,const char* title){
    c.cd(pad);
    TEfficiency* eB=new TEfficiency(*bef,*den);
    TEfficiency* eA=new TEfficiency(*aft,*den);
    eB->SetTitle(Form("%s;#DeltaR_{gen};ID efficiency",title));
    eB->SetLineColor(kBlack); eB->SetMarkerColor(kBlack); eB->SetMarkerStyle(20);
    eA->SetLineColor(kRed);   eA->SetMarkerColor(kRed);   eA->SetMarkerStyle(24);
    eB->Draw("AP"); gPad->Update();
    if(eB->GetPaintedGraph()){ eB->GetPaintedGraph()->GetYaxis()->SetRangeUser(0.,1.05);
                               eB->GetPaintedGraph()->GetXaxis()->SetLimits(LO,HI); }
    eA->Draw("P SAME");
    TLegend* l=new TLegend(0.55,0.2,0.88,0.38);
    l->AddEntry(eB,"before (std HEEP V70)","lp");
    l->AddEntry(eA,"after (modified HEEP)","lp");
    l->Draw();
  };
  drawPad(1,denEB,befEB,aftEB,"EB");
  drawPad(2,denEE,befEE,aftEE,"EE");
  std::string png=outdir+"/heepEffi_vs_dRgen.png";
  c.SaveAs(png.c_str());
  printf("\nplot saved: %s\n", png.c_str());
  printf("hists saved: hists_heepEffi.root\n");
  return 0;
}
