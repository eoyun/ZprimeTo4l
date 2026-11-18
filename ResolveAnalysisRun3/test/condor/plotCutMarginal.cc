// Plotter: single-cut (marginal) PASS efficiency for each of the 12 VID cuts, restricted to a
// dR_gen region (default dR_gen > 0.1). x-axis = the 12 cuts; 3 IDs on one canvas
// (HEEP V70, modified HEEP, modified HEEP Run3); EB and EE separate. 1 canvas = 1 plot.
//   eff(cut c) = (denominator electrons with bit c set, in dR region) / (denominator, in dR region)
// Reuses the existing per-cut histograms (<cutkey>_<id>_<region>) and denominator_<region>;
// NO re-fill needed. For each: sample (15), H mass (3; summed), all (1).
// Usage: plotCutMarginal <outdir> <drmin> <label1>=<file1.root> ...
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
// integral over bins with low edge >= drmin (i.e. the dR > drmin region), including overflow
static double integ(TH1D* h,double drmin){ if(!h) return 0; int b0=h->FindBin(drmin+1e-6); return h->Integral(b0,h->GetNbinsX()+1); }

static const char* IDN[3]   = {"v70","mod","run3"};
static const char* IDLAB[3] = {"HEEP V70","modified HEEP","modified HEEP Run3"};
static int IDCOL[3], IDMK[3] = {20,24,21};

static void drawRegion(double drmin, const std::string& outdir, const std::string& name,
                       const std::string& title, const std::vector<std::string>& files,
                       const std::string& region){
  TH1D* den=sumHist(files,"denominator_"+region);
  double dInt=integ(den,drmin);
  if(dInt<=0){ printf("no denominator (dR>%.2f) for %s %s\n", drmin, name.c_str(), region.c_str()); return; }
  std::vector<TH1D*> hs;
  for(int d=0;d<3;d++){
    TH1D* he=new TH1D((name+"_"+IDN[d]+"_"+region).c_str(),"",12,0,12);
    for(int c=0;c<12;c++){
      TH1D* hc=sumHist(files,std::string(heep::kCutKeys[c])+"_"+IDN[d]+"_"+region);
      double eff=integ(hc,drmin)/dInt;
      he->SetBinContent(c+1,eff); he->SetBinError(c+1,std::sqrt(std::max(eff*(1-eff),0.)/dInt));
      he->GetXaxis()->SetBinLabel(c+1,heep::kCutLabels[c+1]);
    }
    he->SetLineColor(IDCOL[d]); he->SetMarkerColor(IDCOL[d]); he->SetMarkerStyle(IDMK[d]); he->SetMarkerSize(1.0);
    hs.push_back(he);
  }
  gStyle->SetOptStat(0);
  TCanvas c("c","",1000,700); c.SetGridy(); c.SetBottomMargin(0.2);
  TLegend* leg=new TLegend(0.68,0.15,0.90,0.32); leg->SetBorderSize(0); leg->SetFillStyle(0);
  for(size_t i=0;i<hs.size();i++){
    hs[i]->SetTitle(Form("%s  %s  (#DeltaR_{gen}>%.2g, single-cut eff);;pass fraction", title.c_str(), region.c_str(), drmin));
    hs[i]->GetYaxis()->SetRangeUser(0.,1.05); hs[i]->GetXaxis()->LabelsOption("v");
    hs[i]->Draw(i==0?"P":"P SAME");
    leg->AddEntry(hs[i], IDLAB[i], "lp");
  }
  leg->Draw();
  std::string path=outdir+"/"+name+"_cutmarginal_"+region+".png";
  c.SaveAs(path.c_str());
  printf("plot saved: %s\n", path.c_str());
}

int main(int argc, char** argv){
  if(argc<4){ printf("usage: %s <outdir> <drmin> <label=file.root> ...\n", argv[0]); return 2; }
  IDCOL[0]=kBlack; IDCOL[1]=kBlue+1; IDCOL[2]=kRed+1;
  std::string outdir=argv[1]; double drmin=atof(argv[2]);
  std::vector<std::string> labels, files;
  for(int i=3;i<argc;i++){ std::string s=argv[i]; auto p=s.find('=');
    if(p==std::string::npos){ printf("skip malformed arg %s\n", s.c_str()); continue; }
    labels.push_back(s.substr(0,p)); files.push_back(s.substr(p+1)); }
  if(system(("mkdir -p "+outdir).c_str())){}
  auto both=[&](const std::string& nm,const std::string& ti,const std::vector<std::string>& fs){
    drawRegion(drmin,outdir,nm,ti,fs,"EB"); drawRegion(drmin,outdir,nm,ti,fs,"EE"); };
  for(size_t i=0;i<labels.size();i++) both(labels[i],labels[i],{files[i]});
  std::map<std::string,std::vector<std::string>> byH;
  for(size_t i=0;i<labels.size();i++) byH[hmass(labels[i])].push_back(files[i]);
  for(auto& kv : byH) both(kv.first,kv.first,kv.second);
  both("all_signal","all signal",files);
  printf("done [dR>%.2g]: (%d samples + %d H-mass + 1 all) x {EB,EE}\n", drmin, (int)labels.size(), (int)byH.size());
  return 0;
}
