// Plotter: per-cut marginal PASS efficiency vs dR_gen or reco pT, EB/EE separate,
// 3 IDs on one canvas (HEEP V70, modified HEEP, modified HEEP Run3). 1 canvas = 1 plot.
//   numerator = denominator electrons whose ID bitmap has that cut's bit set;
//   denominator = matched reco electrons.
// cut = one of the 12 keys: et sceta detainseed dphiin sigmaieta e2x5 hoe trkiso caloiso
//       dxy misshits ecaldriven   (bit 0..11)
// For each: each sample (15), each H mass (3; num & den summed), all combined (1).
// Usage: plotCut <outdir> <cutkey> <dr|pt> <label1>=<file1.root> ...
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
    if(h){ if(!out){ out=(TH1D*)h->Clone((hn+"_sum").c_str()); out->SetDirectory(0); } else out->Add(h); }
    tf->Close();
  }
  return out;
}
static std::string hmass(const std::string& s){
  auto p=s.find("_H"); if(p==std::string::npos) return "Hunknown";
  size_t i=p+2; std::string d="H"; while(i<s.size()&&std::isdigit((unsigned char)s[i])){ d+=s[i]; ++i; } return d;
}

struct Mode { std::string suf, xaxis, tag; double xmin,xmax; };
static const char* IDN[3]   = {"v70","mod","run3"};
static const char* IDLAB[3] = {"HEEP V70","modified HEEP","modified HEEP Run3"};
static int IDCOL[3], IDMK[3] = {20,24,21};

static void drawRegion(const Mode& m, const std::string& cut, const std::string& cutlab,
                       const std::string& outdir, const std::string& name, const std::string& title,
                       const std::vector<std::string>& files, const std::string& region){
  TH1D* hden=sumHist(files,"denominator_"+region+m.suf);
  if(!hden||hden->Integral()<=0){ printf("no denominator for %s %s\n", name.c_str(), region.c_str()); return; }
  std::vector<TEfficiency*> effs;
  for(int d=0;d<3;d++){
    TH1D* hn=sumHist(files,cut+"_"+IDN[d]+"_"+region+m.suf);
    if(!hn) continue;
    TEfficiency* e=new TEfficiency(*hn,*hden);
    e->SetLineColor(IDCOL[d]); e->SetMarkerColor(IDCOL[d]); e->SetMarkerStyle(IDMK[d]); e->SetMarkerSize(0.9);
    effs.push_back(e);
  }
  if(effs.empty()) return;
  gStyle->SetOptStat(0);
  TCanvas c("c","",900,700); c.SetLogx(); c.SetGridx(); c.SetGridy();
  for(size_t i=0;i<effs.size();i++){
    effs[i]->SetTitle(Form("%s  %s  (%s pass);%s;pass fraction", title.c_str(), region.c_str(), cutlab.c_str(), m.xaxis.c_str()));
    effs[i]->Draw(i==0?"AP":"P SAME");
    if(i==0){ gPad->Update(); auto g=effs[i]->GetPaintedGraph();
      if(g){ g->GetXaxis()->SetLimits(m.xmin,m.xmax); g->GetYaxis()->SetRangeUser(0.,1.05); } gPad->Update(); }
  }
  TLegend* leg=new TLegend(0.16,0.16,0.48,0.33); leg->SetBorderSize(0); leg->SetFillStyle(0);
  for(size_t i=0;i<effs.size();i++) leg->AddEntry(effs[i], IDLAB[i], "lp");
  leg->Draw();
  std::string path=outdir+"/"+name+"_"+cut+"_"+region+m.suf+".png";
  c.SaveAs(path.c_str());
  printf("plot saved: %s\n", path.c_str());
}

int main(int argc, char** argv){
  if(argc<5){ printf("usage: %s <outdir> <cutkey> <dr|pt> <label=file.root> ...\n", argv[0]); return 2; }
  IDCOL[0]=kBlack; IDCOL[1]=kBlue+1; IDCOL[2]=kRed+1;
  std::string outdir=argv[1], cut=argv[2], modestr=argv[3];
  // resolve cut key -> display label (bit index)
  int bit=-1; for(int c=0;c<12;c++) if(cut==heep::kCutKeys[c]) bit=c;
  if(bit<0){ printf("unknown cut '%s'. valid:", cut.c_str()); for(int c=0;c<12;c++) printf(" %s",heep::kCutKeys[c]); printf("\n"); return 2; }
  std::string cutlab=heep::kCutLabels[bit+1];
  Mode m;
  if(modestr=="dr")      m={"",    "#DeltaR_{gen}", "dr", heep::XMIN, heep::XMAX};
  else if(modestr=="pt") m={"_pt", "p_{T} [GeV]",   "pt", heep::PTMIN, heep::PTMAX};
  else { printf("mode must be dr or pt\n"); return 2; }

  std::vector<std::string> labels, files;
  for(int i=4;i<argc;i++){ std::string s=argv[i]; auto p=s.find('=');
    if(p==std::string::npos){ printf("skip malformed arg %s\n", s.c_str()); continue; }
    labels.push_back(s.substr(0,p)); files.push_back(s.substr(p+1)); }
  if(system(("mkdir -p "+outdir).c_str())){}

  auto both=[&](const std::string& nm,const std::string& ti,const std::vector<std::string>& fs){
    drawRegion(m,cut,cutlab,outdir,nm,ti,fs,"EB"); drawRegion(m,cut,cutlab,outdir,nm,ti,fs,"EE"); };
  for(size_t i=0;i<labels.size();i++) both(labels[i],labels[i],{files[i]});
  std::map<std::string,std::vector<std::string>> byH;
  for(size_t i=0;i<labels.size();i++) byH[hmass(labels[i])].push_back(files[i]);
  for(auto& kv : byH) both(kv.first,kv.first,kv.second);
  both("all_signal","all signal",files);
  printf("done [%s %s]: (%d samples + %d H-mass + 1 all) x {EB,EE}\n", cut.c_str(), m.tag.c_str(), (int)labels.size(), (int)byH.size());
  return 0;
}
