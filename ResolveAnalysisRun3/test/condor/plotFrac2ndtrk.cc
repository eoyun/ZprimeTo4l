// Plotter: fraction of (matched reco) electrons that HAVE a 2nd (additional) track,
// vs dR_gen or reco pT, EB and EE on one canvas. 1 canvas = 1 plot.
//   numerator   = denominator electrons with a 2nd track  (out_addtrk/ denominator_<region>)
//   denominator = all matched reco electrons              (out/ denominator_<region>)
// Both are the same selection except the 2nd-track requirement, so the ratio is the
// 2nd-track-presence rate. For each: sample (15), H mass (3; summed), all (1).
// Usage: plotFrac2ndtrk <outdir> <dr|pt> <dendir> <numdir> <sample1> <sample2> ...
#include "HeepEffiCore.h"
#include <TFile.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <cstdio>

static TH1D* sumHist(const std::vector<std::string>& files, const std::string& hn){
  TH1D* out=0;
  for(const auto& fn : files){
    TFile* tf=TFile::Open(fn.c_str());
    if(!tf||tf->IsZombie()){ if(tf)tf->Close(); continue; }
    TH1D* h=(TH1D*)tf->Get(hn.c_str());
    if(h){ if(!out){ out=(TH1D*)h->Clone((hn+"_"+std::to_string((size_t)out)).c_str()); out->SetDirectory(0); } else out->Add(h); }
    tf->Close();
  }
  return out;
}
static std::string hmass(const std::string& s){
  auto p=s.find("_H"); if(p==std::string::npos) return "Hunknown";
  size_t i=p+2; std::string d="H"; while(i<s.size()&&std::isdigit((unsigned char)s[i])){ d+=s[i]; ++i; } return d;
}
struct Mode { std::string suf, xaxis; double xmin,xmax; };

static void draw(const Mode& m, const std::string& outdir, const std::string& name, const std::string& title,
                 const std::vector<std::string>& denf, const std::vector<std::string>& numf){
  gStyle->SetOptStat(0);
  TCanvas c("c","",900,700); c.SetLogx(); c.SetGridx(); c.SetGridy();
  TLegend* leg=new TLegend(0.15,0.74,0.30,0.87); leg->SetBorderSize(0); leg->SetFillStyle(0);
  const char* regs[2]={"EB","EE"}; int col[2]={kBlack,kRed+1}, mk[2]={20,24};
  bool first=true;
  for(int r=0;r<2;r++){
    TH1D* den=sumHist(denf,"denominator_"+std::string(regs[r])+m.suf);
    TH1D* num=sumHist(numf,"denominator_"+std::string(regs[r])+m.suf);
    if(!den||!num||den->Integral()<=0) continue;
    TEfficiency* e=new TEfficiency(*num,*den);
    e->SetLineColor(col[r]); e->SetMarkerColor(col[r]); e->SetMarkerStyle(mk[r]); e->SetMarkerSize(0.9);
    e->SetTitle(Form("%s  (fraction with 2nd track);%s;fraction with 2nd track", title.c_str(), m.xaxis.c_str()));
    e->Draw(first?"AP":"P SAME");
    if(first){ gPad->Update(); auto g=e->GetPaintedGraph(); if(g){ g->GetXaxis()->SetLimits(m.xmin,m.xmax); g->GetYaxis()->SetRangeUser(0.,1.05);} gPad->Update(); first=false; }
    leg->AddEntry(e, regs[r], "lp");
  }
  leg->Draw();
  std::string path=outdir+"/"+name+"_frac2ndtrk"+m.suf+".png";
  c.SaveAs(path.c_str());
  printf("plot saved: %s\n", path.c_str());
}

int main(int argc, char** argv){
  if(argc<6){ printf("usage: %s <outdir> <dr|pt> <dendir> <numdir> <sample...>\n", argv[0]); return 2; }
  std::string outdir=argv[1], modestr=argv[2], dendir=argv[3], numdir=argv[4];
  Mode m;
  if(modestr=="dr")      m={"",    "#DeltaR_{gen}", heep::XMIN, heep::XMAX};
  else if(modestr=="pt") m={"_pt", "p_{T} [GeV]",   heep::PTMIN, heep::PTMAX};
  else { printf("mode must be dr or pt\n"); return 2; }
  std::vector<std::string> samples;
  for(int i=5;i<argc;i++) samples.push_back(argv[i]);
  if(system(("mkdir -p "+outdir).c_str())){}
  auto den=[&](const std::string& s){ return dendir+"/hists_"+s+".root"; };
  auto num=[&](const std::string& s){ return numdir+"/hists_"+s+".root"; };

  for(auto& s : samples) draw(m,outdir,s,s,{den(s)},{num(s)});           // per sample
  std::map<std::string,std::vector<std::string>> byH;
  for(auto& s : samples) byH[hmass(s)].push_back(s);
  for(auto& kv : byH){ std::vector<std::string> d,n; for(auto& s:kv.second){d.push_back(den(s));n.push_back(num(s));} draw(m,outdir,kv.first,kv.first,d,n); }
  std::vector<std::string> da,na; for(auto& s:samples){da.push_back(den(s));na.push_back(num(s));}
  draw(m,outdir,"all_signal","all signal",da,na);
  printf("done [%s]: (%d samples + %d H-mass + 1 all) frac-2ndtrk\n", modestr.c_str(), (int)samples.size(), (int)byH.size());
  return 0;
}
