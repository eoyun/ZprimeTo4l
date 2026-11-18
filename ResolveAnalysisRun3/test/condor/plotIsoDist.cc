// Plotter: isolation value distributions for denominator electrons, EB/EE separate.
//   mode trk : 1D track iso (modified eleTrkPtIso vs std dr03TkSumPt), unit-area, log-y,
//              line at the 5 GeV cut;  2D dR_gen vs track iso (one map per version).
//   mode calo: 1D ECAL iso (modified / modified Run3 / std dr03EcalRecHitSumEt) and
//              1D HCAL iso (std dr03HcalTowerSumEt). Empty (unavailable) curves are skipped.
// Last bin includes overflow (values >= 50 GeV). 1 canvas = 1 plot.
// For each: sample (15), H mass (3; summed), all (1).
// Usage: plotIsoDist <outdir> <trk|calo> <label1>=<file1.root> ...
#include "HeepEffiCore.h"
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <string>
#include <vector>
#include <map>
#include <cctype>
#include <cstdio>

static TH1* sumHist(const std::vector<std::string>& files, const std::string& hn){
  TH1* out=0;
  for(const auto& fn : files){
    TFile* tf=TFile::Open(fn.c_str());
    if(!tf||tf->IsZombie()){ if(tf)tf->Close(); continue; }
    TH1* h=(TH1*)tf->Get(hn.c_str());
    if(h){ if(!out){ out=(TH1*)h->Clone((hn+"_sum").c_str()); out->SetDirectory(0); } else out->Add(h); }
    tf->Close();
  }
  return out;
}
static std::string hmass(const std::string& s){
  auto p=s.find("_H"); if(p==std::string::npos) return "Hunknown";
  size_t i=p+2; std::string d="H"; while(i<s.size()&&std::isdigit((unsigned char)s[i])){ d+=s[i]; ++i; } return d;
}

struct Curve { std::string hname, label; int col, style; };

// overlay of unit-area 1D distributions; returns false if nothing drawn
static bool draw1D(const std::string& path, const std::string& title, const std::vector<Curve>& cs,
                   const std::vector<std::string>& files, double cutline){
  gStyle->SetOptStat(0);
  TCanvas c("c","",900,700); c.SetLogy(); c.SetGridx(); c.SetGridy();
  TLegend* leg=new TLegend(0.52,0.72,0.89,0.89); leg->SetBorderSize(0); leg->SetFillStyle(0);
  std::vector<TH1*> drawn; double ymax=0;
  for(const auto& cv : cs){
    TH1* h=sumHist(files,cv.hname);
    if(!h || h->GetEntries()<=0) continue;
    double n=h->GetEntries();
    h->Scale(1./h->Integral(0,h->GetNbinsX()+1));
    h->SetLineColor(cv.col); h->SetLineStyle(cv.style); h->SetLineWidth(2);
    ymax=std::max(ymax,h->GetMaximum());
    leg->AddEntry(h,Form("%s (N=%.0f)",cv.label.c_str(),n),"l");
    drawn.push_back(h);
  }
  if(drawn.empty()) return false;
  for(size_t i=0;i<drawn.size();i++){
    drawn[i]->SetTitle(title.c_str());
    drawn[i]->GetYaxis()->SetTitle("fraction / bin");
    drawn[i]->SetMinimum(1e-5); drawn[i]->SetMaximum(ymax*5);
    drawn[i]->Draw(i==0?"HIST":"HIST SAME");
  }
  if(cutline>0){ TLine* l=new TLine(cutline,1e-5,cutline,ymax*5); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
  leg->Draw();
  c.SaveAs(path.c_str()); printf("plot saved: %s\n", path.c_str());
  return true;
}

static void draw2D(const std::string& path, const std::string& title, const std::string& hname,
                   const std::vector<std::string>& files){
  TH1* h=sumHist(files,hname);
  if(!h || h->GetEntries()<=0) return;
  // normalize each dR_gen column to unit area -> track-iso distribution at fixed dR
  // (removes the variable-bin-width effect on colors)
  for(int ix=1; ix<=h->GetNbinsX(); ix++){
    double s=0; for(int iy=0; iy<=h->GetNbinsY()+1; iy++) s+=h->GetBinContent(ix,iy);
    if(s<=0) continue;
    for(int iy=0; iy<=h->GetNbinsY()+1; iy++) h->SetBinContent(ix,iy,h->GetBinContent(ix,iy)/s);
  }
  h->SetMinimum(1e-4); h->SetMaximum(1.);
  h->GetZaxis()->SetTitle("fraction within #DeltaR_{gen} bin");
  gStyle->SetOptStat(0);
  TCanvas c("c","",900,700); c.SetLogx(); c.SetLogz(); c.SetRightMargin(0.16);
  h->SetTitle(title.c_str());
  h->Draw("COLZ");
  TLine* l=new TLine(heep::XMIN,5.,heep::XMAX,5.); l->SetLineStyle(2); l->SetLineColor(kRed); l->Draw();
  c.SaveAs(path.c_str()); printf("plot saved: %s\n", path.c_str());
}

int main(int argc, char** argv){
  if(argc<4){ printf("usage: %s <outdir> <trk|calo> <label=file.root> ...\n", argv[0]); return 2; }
  std::string outdir=argv[1], mode=argv[2];
  if(mode!="trk" && mode!="calo"){ printf("mode must be trk or calo\n"); return 2; }
  std::vector<std::string> labels, files;
  for(int i=3;i<argc;i++){ std::string s=argv[i]; auto p=s.find('=');
    if(p==std::string::npos){ printf("skip malformed arg %s\n", s.c_str()); continue; }
    labels.push_back(s.substr(0,p)); files.push_back(s.substr(p+1)); }
  if(system(("mkdir -p "+outdir).c_str())){}

  auto one=[&](const std::string& nm,const std::string& ti,const std::vector<std::string>& fs){
    for(std::string r : {"EB","EE"}){
      if(mode=="trk"){
        draw1D(outdir+"/"+nm+"_trkiso_"+r+".png", ti+"  "+r+";track iso [GeV];",
               {{"iso_trk_std_"+r,"std dr03TkSumPt",kBlack,1},{"iso_trk_mod_"+r,"modified (eleTrkPtIso)",kRed+1,1}}, fs, 5.);
        draw2D(outdir+"/"+nm+"_dr_vs_trkiso_std_"+r+".png", ti+"  "+r+"  std dr03TkSumPt;#DeltaR_{gen};track iso [GeV]", "dr_vs_trkiso_std_"+r, fs);
        draw2D(outdir+"/"+nm+"_dr_vs_trkiso_mod_"+r+".png", ti+"  "+r+"  modified track iso;#DeltaR_{gen};track iso [GeV]", "dr_vs_trkiso_mod_"+r, fs);
      } else {
        draw1D(outdir+"/"+nm+"_ecaliso_"+r+".png", ti+"  "+r+";ECAL iso [GeV];",
               {{"iso_ecal_std_"+r,"std dr03EcalRecHitSumEt",kBlack,1},{"iso_ecal_mod_"+r,"modified",kBlue+1,1},
                {"iso_ecal_run3_"+r,"modified Run3",kRed+1,2}}, fs, -1);
        draw1D(outdir+"/"+nm+"_hcaliso_"+r+".png", ti+"  "+r+";HCAL iso [GeV];",
               {{"iso_hcal_std_"+r,"std dr03HcalTowerSumEt",kBlack,1}}, fs, -1);
      }
    }
  };
  for(size_t i=0;i<labels.size();i++) one(labels[i],labels[i],{files[i]});
  std::map<std::string,std::vector<std::string>> byH;
  for(size_t i=0;i<labels.size();i++) byH[hmass(labels[i])].push_back(files[i]);
  for(auto& kv : byH) one(kv.first,kv.first,kv.second);
  one("all_signal","all signal",files);
  printf("done [%s]: %d samples + %d H-mass + 1 all\n", mode.c_str(), (int)labels.size(), (int)byH.size());
  return 0;
}
