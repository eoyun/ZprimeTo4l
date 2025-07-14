#include "TROOT.h"
#include <cstdarg>
#include <cmath>
#include <algorithm>

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

static double retrieveLumi(const std::string& anlyzrEra) {
  if (anlyzrEra=="20UL16APV")
    return 19.5;
  else if (anlyzrEra=="20UL16" || anlyzrEra=="")
    return 16.8;
  else if (anlyzrEra=="20UL17")
    return 41.48;
  else if (anlyzrEra=="20UL18")
    return 59.83;

  return 0.;
}

static void fillHist(TTree* atree, TH1D* ahist, const TString& branch) {
  float val = 0.;
  float wgt = 0.;
  float e1pt, e1eta, e1phi, e2pt, e2eta, e2phi;

  if (branch!="e1e2Pt")
    atree->SetBranchAddress(branch.Data(),&val);

  atree->SetBranchAddress("e1Pt",&e1pt);
  atree->SetBranchAddress("e1Eta",&e1eta);
  atree->SetBranchAddress("e1Phi",&e1phi);
  atree->SetBranchAddress("e2Pt",&e2pt);
  atree->SetBranchAddress("e2Eta",&e2eta);
  atree->SetBranchAddress("e2Phi",&e2phi);
  atree->SetBranchAddress("wgt",&wgt);

  for (int idx=0; idx < atree->GetEntries(); ++idx) {
    atree->GetEntry(idx);

    if (branch=="e1e2Pt") {
      auto lvec1 = ROOT::Math::PtEtaPhiMVector(e1pt,e1eta,e1phi,0.);
      auto lvec2 = ROOT::Math::PtEtaPhiMVector(e2pt,e2eta,e2phi,0.);
      val = (lvec1+lvec2).Pt();
    }

    ahist->Fill(val,wgt);
  }

  return;
}

void drawStackFakeEle(TString era) {
  TH1::SetDefaultSumw2();

  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Internal";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
  } else if (era=="run2") {
    lumi_sqrtS = "";
    lumi_13TeV = "137.6 fb^{-1}";
  } else {
    std::cout << "check era..." << std::endl;
  }

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

  if( iPos==0 )
    relPosX = 0.12;

  int W = 800;
  int H = 800;

  int H_ref = 800;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  auto canvas = std::make_unique<TCanvas>("canvas","canvas",800,800);

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  TPad* p1 = new TPad("p1","",0,0.3,1,0.95);
  p1->SetFillColor(0);
  p1->SetFrameBorderSize(0);
  p1->SetBorderMode(0);
  p1->SetFrameFillStyle(0);
  p1->SetFrameBorderMode(0);
  p1->SetTickx(0);
  p1->SetTicky(0);
  p1->SetBottomMargin(0.00);
  p1->SetLeftMargin( L/W );
  p1->SetRightMargin( R/W );
  p1->Draw();

  TPad* p2 = new TPad("p2","",0,0,1,0.3);
  p2->SetTopMargin(0.0);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin( L/W );
  p2->SetRightMargin( R/W );
  p2->SetTickx(0);
  p2->SetTicky(0);
  p2->Draw();

  auto SaveAs = [&] (TCanvas* cv, const std::string& name, TPad* pad = nullptr) {
    cv->Update();

    // writing the lumi information and the CMS "logo"
    CMS_lumi( cv, iPeriod, iPos );

    if (pad) {
      pad->RedrawAxis();
      pad->GetFrame()->Draw();
    } else {
      cv->Update();
      cv->RedrawAxis();
      cv->GetFrame()->Draw();
    }

    cv->SaveAs(name.c_str());
  };

  class sample {
  public:
    TFile* file_;
    double xsec_;
    double lumi_;
    int color_;

    sample(TFile* afile, const double axsec, const double alumi, const int acolor) {
      file_ = afile;
      xsec_ = axsec;
      lumi_ = alumi;
      color_ = acolor;
    }
  };

  class samples {
  public:
    TString era_;
    double lumi_;
    std::vector<sample> samples_;
    std::vector<TH1D*> hists_;

    samples(TString aera, const double alumi) {
      era_ = aera;
      lumi_ = alumi;

      samples_ = {
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ZZ.root"),13.81,lumi_,TColor::GetColor("#7a21dd")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WZFXFX.root"),5.213,lumi_,TColor::GetColor("#9c9ca1")), // x0.65
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WW.root"),11.09	,lumi_,TColor::GetColor("#964a8b")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_tW_antitop.root"),79.3*0.5,lumi_,kGreen+2), // XSDB 32.51
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_tW_top.root"),79.3*0.5,lumi_,kGreen+2), // XSDB 32.45
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_t-ch_antitop.root"),80.0,lumi_,kGreen-3), // XSDB 71.75
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_t-ch_top.root"),134.2,lumi_,kGreen-3), // XSDB 119.7
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_s-ch.root"),6.839,lumi_,kGreen-6), // XSDB 3.549 // NNLO 6.839
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_TTTo2L2Nu.root"),833.9*0.1062,lumi_,TColor::GetColor("#e42536")), // 88.29
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_TTsemi.root"),833.9*0.4394,lumi_,kRed-7),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_TThad.root"),833.9*0.4544,lumi_,kRed-9),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_DY_2J.root"),353.6,lumi_,TColor::GetColor("#f89c20")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_DY_1J.root"),983.5,lumi_,TColor::GetColor("#f89c20")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_DY_0J.root"),5090.0,lumi_,TColor::GetColor("#f89c20")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WJets_2J.root"),3276.0,lumi_,TColor::GetColor("#5790fc")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WJets_1J.root"),8832.0,lumi_,TColor::GetColor("#5790fc")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WJets_0J.root"),52780.0	,lumi_,TColor::GetColor("#5790fc")),
      };
    }

    void retrieveHists(const TString& anlyzr, const TString& aname, const TString& branch, double xmin, double xmax, int rebin=1) {
      for (auto ahist : hists_)
        delete ahist;

      hists_.clear();

      TString postfix = era_;

      for (const auto& el : samples_) {
        const double sumwgt = ((TH1D*)el.file_->Get("evtCounter/h_sumW"))->GetBinContent(1);
        TTree* atree = (TTree*)el.file_->Get(anlyzr+postfix+"/"+aname);
        TH1D* ahist = new TH1D( aname + branch + TString(el.file_->GetName()), aname+branch, 50, xmin, xmax );

        fillHist(atree,ahist,branch);

        hists_.push_back(ahist);
        hists_.back()->SetFillColor(el.color_);
        hists_.back()->SetLineColor(el.color_);
        hists_.back()->Scale(lumi_*1000.*el.xsec_/sumwgt);
        hists_.back()->Rebin(rebin);

        delete atree;
      }
    }

    void addHists(const samples& other) {
      for (unsigned idx=0; idx<hists_.size(); idx++)
        hists_.at(idx)->Add(other.hists_.at(idx));
    }

    std::unique_ptr<THStack> extractStack() {
      TString aformat;
      aformat.Form("%s;%s;%s",hists_.back()->GetTitle(),hists_.back()->GetXaxis()->GetTitle(),hists_.back()->GetYaxis()->GetTitle());
      auto astack = std::make_unique<THStack>("stack",aformat);

      for (auto* ahist : hists_)
        astack->Add(ahist);

      return std::move(astack);
    }

    std::unique_ptr<TH1D> extractNominal() {
      auto ahist = std::unique_ptr<TH1D>((TH1D*)hists_.front()->Clone());

      for (unsigned idx=1; idx<hists_.size(); idx++)
        ahist->Add(hists_.at(idx));

      return std::move(ahist);
    }
  };

  TFile* datafile1 = new TFile("FakeEleAnalyzer_20UL16APV_data.root","READ");
  TFile* datafile2 = new TFile("FakeEleAnalyzer_20UL16_data.root","READ");
  TFile* datafile3 = new TFile("FakeEleAnalyzer_20UL17_data.root","READ");
  TFile* datafile4 = new TFile("FakeEleAnalyzer_20UL18_data.root","READ");

  auto retrieveDataHist = [] (TFile* afile, const TString& anlyzr, const TString& aname, const TString& branch, double xmin, double xmax, int rebin=1) {
    TTree* atree = (TTree*)afile->Get(anlyzr+"/"+aname);
    TH1D* ahist = new TH1D( aname + branch + TString(afile->GetName()), aname+branch, 50, xmin, xmax );

    fillHist(atree,ahist,branch);

    ahist->SetLineColor(kBlack);
    ahist->Rebin(rebin);

    delete atree;

    return ahist;
  };

  auto sample1 = samples("20UL16APV",19.5);
  auto sample2 = samples("20UL16",16.8);
  auto sample3 = samples("20UL17",41.48);
  auto sample4 = samples("20UL18",59.83);

  struct histo {
    TString name_;
    const double xmin_;
    const double xmax_;

    histo(TString name, double xmin, double xmax)
    : name_(name), xmin_(xmin), xmax_(xmax) {}
  };

  std::vector<TString> treeNames = {"denomTree1m","numerTree1m","denomTree2m","numerTree2m"};
  std::vector<histo> varNames = {histo("e1scEta",-3.,3.),histo("e1e2InvM",0.,20.),histo("e1e2Pt",0.,200.),histo("e1e2DR",0.,0.4),};

  for (auto atree : treeNames) {
    for (auto ahist : varNames) {

      auto* datahist = retrieveDataHist(datafile1,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_);
      datahist->Add( retrieveDataHist(datafile2,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_) );
      datahist->Add( retrieveDataHist(datafile3,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_) );
      datahist->Add( retrieveDataHist(datafile4,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_) );
      sample1.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_);
      sample2.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_);
      sample3.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_);
      sample4.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_);
      sample1.addHists(sample2);
      sample1.addHists(sample3);
      sample1.addHists(sample4);

      auto astack = sample1.extractStack();

      p1->cd();

      auto nominal = sample1.extractNominal();
      double theMax = nominal->GetMaximum();

      datahist->SetLineColor(kBlack);
      datahist->SetLineWidth(2);

      datahist->SetMaximum( 1.5*datahist->GetMaximum() );
      datahist->SetMinimum( 0.01 );

      datahist->Draw("E1");

      astack->Draw("hist&same");
      datahist->Draw("E1&same");

      auto legend = std::make_unique<TLegend>(0.6,0.55,0.95,0.9);
      legend->SetBorderSize(0);
      legend->SetNColumns(2);
      legend->AddEntry(datahist,"Data");
      legend->AddEntry(sample1.hists_.at(14),"WJets");
      legend->AddEntry(sample1.hists_.at(11),"DY");
      legend->AddEntry(sample1.hists_.at(10),"TThad");
      legend->AddEntry(sample1.hists_.at(9),"TTsemi");
      legend->AddEntry(sample1.hists_.at(8),"TTlep");
      legend->AddEntry(sample1.hists_.at(7),"ST_s-ch");
      legend->AddEntry(sample1.hists_.at(5),"ST_t-ch");
      legend->AddEntry(sample1.hists_.at(3),"ST_tW");
      legend->AddEntry(sample1.hists_.at(2),"WW");
      legend->AddEntry(sample1.hists_.at(1),"WZ");
      legend->AddEntry(sample1.hists_.at(0),"ZZ");
      legend->Draw();

      p2->cd();

      TH1* ratio = (TH1*)datahist->Clone("ratio");
      ratio->SetStats(0);
      ratio->SetTitle("");
      ratio->Divide(nominal.get());
      ratio->GetYaxis()->SetTitle("Data/MC");
      ratio->GetYaxis()->SetTitleSize(0.1);
      ratio->GetYaxis()->SetTitleOffset(0.4);
      ratio->GetXaxis()->SetLabelSize(0.1);
      ratio->GetYaxis()->SetLabelSize(0.1);
      ratio->GetXaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetRangeUser(0.01,1.99);
      ratio->GetYaxis()->SetNdivisions(505);
      ratio->GetXaxis()->SetTitleSize(0.12);
      ratio->GetXaxis()->SetTitleOffset(1.0);
      ratio->SetLineColor(kBlack);

      if (ahist.name_.Contains("e1e2Pt"))
        ratio->GetXaxis()->SetTitle("p_{T}(ee) [GeV]");

      ratio->Draw("E1");

      canvas->Update();

      CMS_lumi( canvas.get(), iPeriod, iPos );

      p1->RedrawAxis();
      p1->GetFrame()->Draw();

      canvas->SaveAs(atree+"_"+ahist.name_+".png");

      if (ahist.name_.Contains("e1e2Pt"))
        canvas->SaveAs(atree+"_"+ahist.name_+".pdf");
    }
  }

  return;
}
