// Plotter: nominal HEEP (V70) vs modified HEEP ID efficiency, EB and EE as SEPARATE plots.
// Variable on x-axis selected by mode: "dr" (dR_gen) or "pt" (reco pT). 1 canvas = 1 plot.
// For each of: each sample (15), each H mass (3; num & den summed), all combined (1)
//   => (15 + 3 + 1) x {EB, EE} = 38 images per mode.
// Usage: plotHeepEffi <outdir> <dr|pt> <label1>=<file1.root> ...
#include "HeepEffiCore.h"
#include <TFile.h>
#include <TH1D.h>
#include <TEfficiency.h>
#include <TGraphAsymmErrors.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TROOT.h>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <cstdio>

struct Mode { std::string suf, xaxis, tag; double xmin, xmax; };

static TH1D* sumHist(const std::vector<std::string>& files, const std::string& hn){
  TH1D* out=0;
  for(const auto& fn : files){
    TFile* tf=TFile::Open(fn.c_str());
    if(!tf||tf->IsZombie()){ if(tf)tf->Close(); continue; }
    TH1D* h=(TH1D*)tf->Get(hn.c_str());
    if(h){ if(!out){ out=(TH1D*)h->Clone((hn+"_sum").c_str()); out->SetDirectory(0); } else out->Add(h); }
    tf->Close();
  }
  return out;
}

static std::string hmass(const std::string& s){
  auto p=s.find("_H"); if(p==std::string::npos) return "Hunknown";
  size_t i=p+2; std::string d="H";
  while(i<s.size() && std::isdigit((unsigned char)s[i])){ d+=s[i]; ++i; }
  return d;
}

// One region (EB or EE): nominal HEEP (before) vs modified HEEP (after) on one canvas.
static void drawRegion(const Mode& m, const std::string& outdir, const std::string& name,
                       const std::string& title, const std::vector<std::string>& files,
                       const std::string& region){
  TH1D* hden=sumHist(files,"denominator_"+region+m.suf);
  TH1D* hbef=sumHist(files,"numerator_before_"+region+m.suf);
  TH1D* haft=sumHist(files,"numerator_after_"+region+m.suf);
  TH1D* hr3 =sumHist(files,"numerator_run3_"+region+m.suf);
  if(!hden||!hbef||!haft||!hr3||hden->Integral()<=0){ printf("no data for %s %s (%s)\n", name.c_str(), region.c_str(), m.tag.c_str()); return; }

  TEfficiency* eN =new TEfficiency(*hbef,*hden);  // nominal HEEP V70
  TEfficiency* eM =new TEfficiency(*haft,*hden);  // modified HEEP
  TEfficiency* eR3=new TEfficiency(*hr3, *hden);  // modified HEEP Run3
  eN ->SetLineColor(kBlack);  eN ->SetMarkerColor(kBlack);  eN ->SetMarkerStyle(20); eN ->SetMarkerSize(0.9);
  eM ->SetLineColor(kBlue+1); eM ->SetMarkerColor(kBlue+1); eM ->SetMarkerStyle(24); eM ->SetMarkerSize(0.9);
  eR3->SetLineColor(kRed+1);  eR3->SetMarkerColor(kRed+1);  eR3->SetMarkerStyle(21); eR3->SetMarkerSize(0.9);

  gStyle->SetOptStat(0);
  TCanvas c("c","",900,700); c.SetLogx(); c.SetGridx(); c.SetGridy();
  eN->SetTitle(Form("%s  %s;%s;ID efficiency", title.c_str(), region.c_str(), m.xaxis.c_str()));
  eN->Draw("AP"); gPad->Update();
  if(auto g=eN->GetPaintedGraph()){ g->GetXaxis()->SetLimits(m.xmin,m.xmax); g->GetYaxis()->SetRangeUser(0.,1.05); }
  gPad->Update();
  eM->Draw("P SAME"); eR3->Draw("P SAME");
  TLegend* leg=new TLegend(0.16,0.16,0.46,0.33); leg->SetBorderSize(0); leg->SetFillStyle(0);
  leg->AddEntry(eN, "nominal HEEP (V70)","lp");
  leg->AddEntry(eM, "modified HEEP","lp");
  leg->AddEntry(eR3,"modified HEEP Run3","lp");
  leg->Draw();
  std::string path=outdir+"/"+name+"_"+region+m.suf+".png";
  c.SaveAs(path.c_str());
  printf("plot saved: %s\n", path.c_str());
}

static void drawBoth(const Mode& m, const std::string& outdir, const std::string& name,
                     const std::string& title, const std::vector<std::string>& files){
  drawRegion(m,outdir,name,title,files,"EB");
  drawRegion(m,outdir,name,title,files,"EE");
}

int main(int argc, char** argv){
  if(argc<4){ printf("usage: %s <outdir> <dr|pt> <label=file.root> ...\n", argv[0]); return 2; }
  std::string outdir=argv[1], modestr=argv[2];
  Mode m;
  if(modestr=="dr")      m={"",    "#DeltaR_{gen}", "dr", heep::XMIN, heep::XMAX};
  else if(modestr=="pt") m={"_pt", "p_{T} [GeV]",   "pt", heep::PTMIN, heep::PTMAX};
  else { printf("mode must be dr or pt (got %s)\n", modestr.c_str()); return 2; }

  std::vector<std::string> labels, files;
  for(int i=3;i<argc;i++){ std::string s=argv[i]; auto p=s.find('=');
    if(p==std::string::npos){ printf("skip malformed arg %s\n", s.c_str()); continue; }
    labels.push_back(s.substr(0,p)); files.push_back(s.substr(p+1)); }
  if(system(("mkdir -p "+outdir).c_str())){}

  for(size_t i=0;i<labels.size();i++) drawBoth(m, outdir, labels[i], labels[i], {files[i]});   // per sample
  std::map<std::string,std::vector<std::string>> byH;                                           // per H mass
  for(size_t i=0;i<labels.size();i++) byH[hmass(labels[i])].push_back(files[i]);
  for(auto& kv : byH) drawBoth(m, outdir, kv.first, kv.first, kv.second);
  drawBoth(m, outdir, "all_signal", "all signal", files);                                        // all

  printf("done [%s]: (%d samples + %d H-mass + 1 all) x {EB,EE}\n", m.tag.c_str(), (int)labels.size(), (int)byH.size());
  return 0;
}
