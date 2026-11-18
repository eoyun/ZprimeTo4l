// Plotter: modified-HEEP sequential cutflow, EB vs EE on one canvas (survival fraction
// normalized to the denominator). 1 canvas = 1 plot, for each of:
//   - each sample (15), each H mass (3; counts summed), all combined (1) => 19 images.
// Usage: plotCutflow <outdir> <label1>=<file1.root> ...
#include "HeepEffiCore.h"
#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <cmath>
#include <cstdio>

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
  size_t i=p+2; std::string d="H"; while(i<s.size()&&std::isdigit((unsigned char)s[i])){ d+=s[i]; ++i; } return d;
}
// normalize counts to the denominator (bin1) -> survival fraction with binomial errors
static TH1D* toFrac(TH1D* counts){
  if(!counts) return 0;
  TH1D* f=(TH1D*)counts->Clone((std::string(counts->GetName())+"_frac").c_str()); f->SetDirectory(0);
  double nd=counts->GetBinContent(1);
  for(int b=1;b<=f->GetNbinsX();b++){
    double n=counts->GetBinContent(b);
    double eff = nd>0? n/nd : 0.;
    double err = nd>0? std::sqrt(std::max(eff*(1-eff),0.)/nd) : 0.;
    f->SetBinContent(b,eff); f->SetBinError(b,err);
  }
  return f;
}

struct Mode { std::string hbase, idtitle, fsuf; };

static void drawCutflow(const Mode& m, const std::string& outdir, const std::string& name,
                        const std::string& title, const std::vector<std::string>& files){
  TH1D* cEB=sumHist(files,m.hbase+"_EB");
  TH1D* cEE=sumHist(files,m.hbase+"_EE");
  if((!cEB||cEB->GetBinContent(1)<=0) && (!cEE||cEE->GetBinContent(1)<=0)){ printf("no cutflow for %s (%s)\n", name.c_str(), m.idtitle.c_str()); return; }
  TH1D* fEB=toFrac(cEB); TH1D* fEE=toFrac(cEE);

  gStyle->SetOptStat(0);
  TCanvas c("c","",1000,700); c.SetGridy(); c.SetBottomMargin(0.18);
  std::string ttl=Form("%s  %s cutflow;;survival fraction", title.c_str(), m.idtitle.c_str());
  bool first=true;
  if(fEB){ fEB->SetLineColor(kBlack); fEB->SetMarkerColor(kBlack); fEB->SetMarkerStyle(20);
    fEB->SetTitle(ttl.c_str());
    fEB->GetYaxis()->SetRangeUser(0.,1.05); fEB->GetXaxis()->LabelsOption("v");
    fEB->Draw("P"); first=false; }
  if(fEE){ fEE->SetLineColor(kRed+1); fEE->SetMarkerColor(kRed+1); fEE->SetMarkerStyle(24);
    if(first){ fEE->SetTitle(ttl.c_str());
      fEE->GetYaxis()->SetRangeUser(0.,1.05); fEE->GetXaxis()->LabelsOption("v"); fEE->Draw("P"); first=false; }
    else fEE->Draw("P SAME"); }

  TLegend* leg=new TLegend(0.70,0.78,0.88,0.90); leg->SetBorderSize(0); leg->SetFillStyle(0);
  if(fEB) leg->AddEntry(fEB,Form("EB (N=%.0f)",cEB?cEB->GetBinContent(1):0),"lp");
  if(fEE) leg->AddEntry(fEE,Form("EE (N=%.0f)",cEE?cEE->GetBinContent(1):0),"lp");
  leg->Draw();
  std::string path=outdir+"/"+name+m.fsuf+".png";
  c.SaveAs(path.c_str());
  printf("plot saved: %s\n", path.c_str());
}

int main(int argc, char** argv){
  if(argc<4){ printf("usage: %s <outdir> <mod|v70|modrun3> <label=file.root> ...\n", argv[0]); return 2; }
  std::string outdir=argv[1], modestr=argv[2];
  Mode m;
  if(modestr=="mod")            m={"cutflow",          "modified-HEEP",      "_cutflow"};
  else if(modestr=="v70")       m={"cutflow_v70",      "HEEP V70",           "_cutflow_v70"};
  else if(modestr=="modrun3")   m={"cutflow_r3",       "modified-HEEP Run3", "_cutflow_run3"};
  else if(modestr=="mod_hidr")    m={"cutflow_hidr",     "modified-HEEP (#DeltaR_{gen}>0.1)",      "_cutflow_hidr"};
  else if(modestr=="v70_hidr")    m={"cutflow_v70_hidr", "HEEP V70 (#DeltaR_{gen}>0.1)",           "_cutflow_v70_hidr"};
  else if(modestr=="modrun3_hidr")m={"cutflow_r3_hidr",  "modified-HEEP Run3 (#DeltaR_{gen}>0.1)", "_cutflow_run3_hidr"};
  else { printf("mode must be mod, v70, modrun3, mod_hidr, v70_hidr or modrun3_hidr (got %s)\n", modestr.c_str()); return 2; }

  std::vector<std::string> labels, files;
  for(int i=3;i<argc;i++){ std::string s=argv[i]; auto p=s.find('=');
    if(p==std::string::npos){ printf("skip malformed arg %s\n", s.c_str()); continue; }
    labels.push_back(s.substr(0,p)); files.push_back(s.substr(p+1)); }
  if(system(("mkdir -p "+outdir).c_str())){}

  for(size_t i=0;i<labels.size();i++) drawCutflow(m, outdir, labels[i], labels[i], {files[i]});     // per sample
  std::map<std::string,std::vector<std::string>> byH;
  for(size_t i=0;i<labels.size();i++) byH[hmass(labels[i])].push_back(files[i]);
  for(auto& kv : byH) drawCutflow(m, outdir, kv.first, kv.first, kv.second);                         // per H mass
  drawCutflow(m, outdir, "all_signal", "all signal", files);                                         // all

  printf("done [%s]: %d samples + %d H-mass + 1 all = %d cutflow plots\n",
         m.idtitle.c_str(),(int)labels.size(),(int)byH.size(),(int)labels.size()+(int)byH.size()+1);
  return 0;
}
