#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

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

void runREMuFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static constexpr double WZxsec_ = 5.213; // 0.65*62.78;
  static constexpr double ZZxsec_ = 13.81;
  static double valLumi = retrieveLumi(era.Data());
  static TString postfix = era;
  TString fname = era;

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
    postfix = "";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
  } else if (era=="run2") {
    lumi_sqrtS = "";
    lumi_13TeV = "137.6 fb^{-1}";
    postfix = "";
    fname = "20UL16";
  } else {
    std::cout << "check era..." << std::endl;
  }

  static TString anlyzrMC = "resolvedEMuCRanalyzer"+postfix;
  static TString anlyzrData = "resolvedEMuCRanalyzerData";

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

  TFile* datafile = new TFile("MuAnalyzer_"+fname+"_data.root","READ");
  TFile* WZfile = new TFile("MuAnalyzer_"+fname+"_WZFXFX.root","READ");
  TFile* ZZfile = new TFile("MuAnalyzer_"+fname+"_ZZ.root","READ");

  TFile* datafile1 = new TFile("MuAnalyzer_20UL16APV_data.root","READ");
  TFile* WZfile1 = new TFile("MuAnalyzer_20UL16APV_WZFXFX.root","READ");
  TFile* ZZfile1 = new TFile("MuAnalyzer_20UL16APV_ZZ.root","READ");

  TFile* datafile2 = new TFile("MuAnalyzer_20UL17_data.root","READ");
  TFile* WZfile2 = new TFile("MuAnalyzer_20UL17_WZFXFX.root","READ");
  TFile* ZZfile2 = new TFile("MuAnalyzer_20UL17_ZZ.root","READ");

  TFile* datafile3 = new TFile("MuAnalyzer_20UL18_data.root","READ");
  TFile* WZfile3 = new TFile("MuAnalyzer_20UL18_WZFXFX.root","READ");
  TFile* ZZfile3 = new TFile("MuAnalyzer_20UL18_ZZ.root","READ");

  //TFile* H250A1file = new TFile("REMuCR_20UL16_H250A1.root","READ");
  // TFile* H750A1file = new TFile("REMuCR_"+era+"_H750A1.root","READ");
  //TFile* H2000A1file = new TFile("REMuCR_20UL16_H2000A1.root","READ");
  //TFile* H250A10file = new TFile("REMuCR_20UL16_H250A10.root","READ");
  //TFile* H750A10file = new TFile("REMuCR_20UL16_H750A10.root","READ");
  //TFile* H2000A10file = new TFile("REMuCR_20UL16_H2000A10.root","READ");

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  /*auto sigsamples = std::vector<SigSample>{
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A750.root","READ"),"H2000A750")
  };

  auto sigsamples1 = std::vector<SigSample>{
    SigSample(new TFile("MuAnalyzer_20UL16APV_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A750.root","READ"),"H2000A750")
  };

  auto sigsamples2 = std::vector<SigSample>{
    SigSample(new TFile("MuAnalyzer_20UL17_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("MuAnalyzer_20UL17_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("MuAnalyzer_20UL17_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("MuAnalyzer_20UL17_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("MuAnalyzer_20UL17_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("MuAnalyzer_20UL17_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("MuAnalyzer_20UL17_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("MuAnalyzer_20UL17_H2000A750.root","READ"),"H2000A750")
  };

  auto sigsamples3 = std::vector<SigSample>{
    SigSample(new TFile("MuAnalyzer_20UL18_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("MuAnalyzer_20UL18_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("MuAnalyzer_20UL18_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("MuAnalyzer_20UL18_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("MuAnalyzer_20UL18_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("MuAnalyzer_20UL18_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("MuAnalyzer_20UL18_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("MuAnalyzer_20UL18_H2000A750.root","READ"),"H2000A750")
  };*/

  auto sigsamples = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_"+fname+"_H750A100.root","READ"),"H750A100")};
  auto sigsamples1 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL16APV_H750A100.root","READ"),"H750A100")};
  auto sigsamples2 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL17_H750A100.root","READ"),"H750A100")};
  auto sigsamples3 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL18_H750A100.root","READ"),"H750A100")};

  auto* canvas_2 = new TCanvas("canvas_2","canvas_2",50,50,W,H);
  canvas_2->SetFillColor(0);
  canvas_2->SetBorderMode(0);
  canvas_2->SetFrameFillStyle(0);
  canvas_2->SetFrameBorderMode(0);
  canvas_2->SetLeftMargin( L/W );
  canvas_2->SetRightMargin( R/W );
  canvas_2->SetTopMargin( T/H );
  canvas_2->SetBottomMargin( B/H );
  canvas_2->SetTickx(0);
  canvas_2->SetTicky(0);

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
  p2->SetTickx(0);
  p2->SetTicky(0);
  p2->SetTopMargin(0.0);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin( L/W );
  p2->SetRightMargin( R/W );
  p2->Draw();

  auto* canvas_1 = new TCanvas("canvas_1","canvas_1",50,50,W,H);
  canvas_1->SetFillColor(0);
  canvas_1->SetBorderMode(0);
  canvas_1->SetFrameFillStyle(0);
  canvas_1->SetFrameBorderMode(0);
  canvas_1->SetLeftMargin( L/W );
  canvas_1->SetRightMargin( R/W );
  canvas_1->SetTopMargin( T/H );
  canvas_1->SetBottomMargin( B/H );
  canvas_1->SetTickx(0);
  canvas_1->SetTicky(0);

  auto SaveAs = [&] (TCanvas* canvas, const std::string& name, TPad* pad = nullptr) {
    canvas->Update();

    // writing the lumi information and the CMS "logo"
    CMS_lumi( canvas, iPeriod, iPos );

    if (pad) {
      pad->RedrawAxis();
      pad->GetFrame()->Draw();
    } else {
      canvas->Update();
      canvas->RedrawAxis();
      canvas->GetFrame()->Draw();
    }

    canvas->SaveAs(name.c_str());
  };

  class HistLoaderBase {
  protected:
    TFile* datafile_;
    TFile* WZfile_;
    TFile* ZZfile_;

  public:
    HistLoaderBase(TFile* adatafile, TFile* aWZfile, TFile* aZZfile) {
      datafile_ = adatafile;
      WZfile_ = aWZfile;
      ZZfile_ = aZZfile;
    };

    ~HistLoaderBase() {};

    TH1D* variateDn(TH1D* nominal, TH1D* up) const {
      TH1D* dn = (TH1D*)nominal->Clone();

      for (int ibin=0; ibin<nominal->GetNbinsX()+2; ibin++) {
        double val = nominal->GetBinContent(ibin);
        double var = up->GetBinContent(ibin);
        double valErr = nominal->GetBinError(ibin);
        double varErr = up->GetBinError(ibin);
        dn->SetBinContent( ibin, std::max( val - (var - val), 0.) );
        dn->SetBinError( ibin, std::hypot(2*valErr,varErr) );
      }

      return dn;
    }

    class SystVariation {
    public:
      SystVariation() : up_(nullptr), dn_(nullptr) {}
      SystVariation(TH1D* up, TH1D* dn) : up_(up), dn_(dn) {}
      ~SystVariation()=default;

      TH1D* up_;
      TH1D* dn_;
    };

  protected:
    std::map<std::string,SystVariation> syst_;
  };

  // 3P1F CR

  class HistLoader3P1F : public HistLoaderBase {
  public:
    HistLoader3P1F(TFile* adatafile, TFile* aWZfile, TFile* aZZfile) : HistLoaderBase(adatafile,aWZfile,aZZfile) {}

    ~HistLoader3P1F() {
      if (dataHist_)
        delete dataHist_, WZHist_, ZZHist_, FFHist_;

      if (FFHist_ffUp_)
        delete WZHist_idUp_, ZZHist_idUp_, WZHist_idDn_, ZZHist_idDn_, FFHist_ffUp_, FFHist_ffDn_, FFHist_ffUpM_, FFHist_ffDnM_;
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* WZHist_ = nullptr;
    TH1D* ZZHist_ = nullptr;
    TH1D* FFHist_ = nullptr;

    TH1D* WZHist_idUp_ = nullptr;
    TH1D* ZZHist_idUp_ = nullptr;
    TH1D* WZHist_idDn_ = nullptr;
    TH1D* ZZHist_idDn_ = nullptr;
    TH1D* FFHist_ffUp_ = nullptr;
    TH1D* FFHist_ffDn_ = nullptr;
    TH1D* FFHist_ffUpM_ = nullptr;
    TH1D* FFHist_ffDnM_ = nullptr;

  public:
    void add(const HistLoader3P1F& other) {
      this->dataHist_->Add(other.dataHist_);
      this->WZHist_->Add(other.WZHist_);
      this->ZZHist_->Add(other.ZZHist_);
      this->FFHist_->Add(other.FFHist_);

      if (FFHist_ffUp_) {
        this->WZHist_idUp_->Add(other.WZHist_idUp_);
        this->ZZHist_idUp_->Add(other.ZZHist_idUp_);
        this->WZHist_idDn_->Add(other.WZHist_idDn_);
        this->ZZHist_idDn_->Add(other.ZZHist_idDn_);
        this->FFHist_ffUp_->Add(other.FFHist_ffUp_);
        this->FFHist_ffDn_->Add(other.FFHist_ffDn_);
        this->FFHist_ffUpM_->Add(other.FFHist_ffUpM_);
        this->FFHist_ffDnM_->Add(other.FFHist_ffDnM_);
      }
    }

    void load(const std::string& nameNum, const std::string& name, std::string anlyzrEra="") {
      if (dataHist_)
        delete dataHist_, WZHist_, ZZHist_, FFHist_;

      if (FFHist_ffUp_) {
        delete WZHist_idUp_, ZZHist_idUp_, WZHist_idDn_, ZZHist_idDn_, FFHist_ffUp_, FFHist_ffDn_, FFHist_ffUpM_, FFHist_ffDnM_;

        FFHist_ffUp_ = nullptr;
      }

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum).c_str() )->Clone();
      WZHist_ = (TH1D*)WZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/")+nameNum).c_str() )->Clone();
      ZZHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/")+nameNum).c_str() )->Clone();
      FFHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name).c_str() )->Clone();

      const double lumi = retrieveLumi(anlyzrEra);

      if (TString(nameNum).Contains("llll_invM")) {
        WZHist_idUp_ = (TH1D*)WZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/")+nameNum+"_idUp").c_str() )->Clone();
        ZZHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/")+nameNum+"_idUp").c_str() )->Clone();
        WZHist_idDn_ = (TH1D*)WZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/")+nameNum+"_idDn").c_str() )->Clone();
        ZZHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/")+nameNum+"_idDn").c_str() )->Clone();
        FFHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffUpE").c_str() )->Clone();
        FFHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffDnE").c_str() )->Clone();
        FFHist_ffUpM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffUpM").c_str() )->Clone();
        FFHist_ffDnM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffDnM").c_str() )->Clone();

        WZHist_idUp_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        WZHist_idDn_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZHist_idUp_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZHist_idDn_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      }

      WZHist_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      WZHist_->SetFillColor(TColor::GetColor("#7a21dd"));
      ZZHist_->SetFillColor(TColor::GetColor("#5790fc"));
      WZHist_->SetLineWidth(0);
      ZZHist_->SetLineWidth(0);
      FFHist_->SetFillColor(TColor::GetColor("#9c9ca1"));
      FFHist_->SetLineWidth(0);
    };

    void compareAndRatio(TPad* padUp, TPad* padDn, int rebin=1) {
      if ( FFHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FFHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (FFHist_ffUp_) {
          FFHist_ffUp_->Rebin( FFHist_ffUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          FFHist_ffDn_->Rebin( FFHist_ffDn_->GetNbinsX()/dataHist_->GetNbinsX() );
          FFHist_ffUpM_->Rebin( FFHist_ffUpM_->GetNbinsX()/dataHist_->GetNbinsX() );
          FFHist_ffDnM_->Rebin( FFHist_ffDnM_->GetNbinsX()/dataHist_->GetNbinsX() );
        }
      }

      WZHist_->Rebin(rebin);
      ZZHist_->Rebin(rebin);
      FFHist_->Rebin(rebin);
      dataHist_->Rebin(rebin);

      if (FFHist_ffUp_) {
        WZHist_idUp_->Rebin(rebin);
        WZHist_idDn_->Rebin(rebin);
        ZZHist_idUp_->Rebin(rebin);
        ZZHist_idDn_->Rebin(rebin);
        FFHist_ffUp_->Rebin(rebin);
        FFHist_ffDn_->Rebin(rebin);
        FFHist_ffUpM_->Rebin(rebin);
        FFHist_ffDnM_->Rebin(rebin);
      }

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(ZZHist_);
      stack->Add(WZHist_);
      stack->Add(FFHist_);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));
      tmpHist->Add((TH1*)stack->GetHists()->At(2));

      padUp->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));
      dataHist_->SetMinimum(0.001);
      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        dataHist_->GetXaxis()->SetRangeUser(0.,500.);
      }

      TLegend* legend = new TLegend(0.75,0.7,0.98,0.92);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(WZHist_,"WZ");
      legend->AddEntry(ZZHist_,"ZZ");
      legend->AddEntry(FFHist_,"Nonprompt");
      legend->Draw();

      TH1D* ratio = (TH1D*)dataHist_->Clone();
      ratio->SetStats(0);
      ratio->SetTitle("");
      ratio->Divide(tmpHist);
      ratio->GetYaxis()->SetTitle("Obs/Exp");
      ratio->GetYaxis()->SetTitleSize(0.1);
      ratio->GetYaxis()->SetTitleOffset(0.4);
      ratio->GetXaxis()->SetLabelSize(0.1);
      ratio->GetYaxis()->SetLabelSize(0.08);
      ratio->GetXaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetLabelOffset(0.005);
      ratio->GetYaxis()->SetRangeUser(0.2,1.8);
      ratio->GetXaxis()->SetTitleSize(0.12);
      ratio->GetXaxis()->SetTitleOffset(0.75);
      ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
      ratio->SetLineColor(kBlack);

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        ratio->GetXaxis()->SetRangeUser(0.,500.);
      }

      padDn->cd();
      ratio->Draw("E1");

      if (FFHist_ffUp_) {
        TH1D* idUp = (TH1D*)WZHist_idUp_->Clone();
        idUp->Add(ZZHist_idUp_);
        idUp->Add(FFHist_);
        TH1D* idDn = (TH1D*)WZHist_idDn_->Clone();
        idDn->Add(ZZHist_idDn_);
        idDn->Add(FFHist_);

        FFHist_ffUp_->Add(WZHist_);
        FFHist_ffUp_->Add(ZZHist_);
        FFHist_ffDn_->Add(WZHist_);
        FFHist_ffDn_->Add(ZZHist_);
        FFHist_ffUpM_->Add(WZHist_);
        FFHist_ffUpM_->Add(ZZHist_);
        FFHist_ffDnM_->Add(WZHist_);
        FFHist_ffDnM_->Add(ZZHist_);

        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRdn, errRup;

        for (unsigned idx = 1; idx <= tmpHist->GetNbinsX(); idx++) {
          x0.push_back(tmpHist->GetBinCenter(idx));
          y0.push_back(tmpHist->GetBinContent(idx));
          errx.push_back(tmpHist->GetBinWidth(idx)/2.);

          double valIdUp = idUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valIdDn = tmpHist->GetBinContent(idx) - idDn->GetBinContent(idx);
          double valFFup = FFHist_ffUp_->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdn = tmpHist->GetBinContent(idx) - FFHist_ffDn_->GetBinContent(idx);
          double valFFupM = FFHist_ffUpM_->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdnM = tmpHist->GetBinContent(idx) - FFHist_ffDnM_->GetBinContent(idx);

          erryUp.push_back( std::sqrt(valIdUp*valIdUp + valFFup*valFFup + valFFupM*valFFupM) );
          erryDn.push_back( std::sqrt(valIdDn*valIdDn + valFFdn*valFFdn + valFFdnM*valFFdnM) );

          r0.push_back(1.);

          double rIdUp = valIdUp/tmpHist->GetBinContent(idx);
          double rIdDn = valIdDn/tmpHist->GetBinContent(idx);
          double rFFup = valFFup/tmpHist->GetBinContent(idx);
          double rFFdn = valFFdn/tmpHist->GetBinContent(idx);
          double rFFupM = valFFupM/tmpHist->GetBinContent(idx);
          double rFFdnM = valFFdnM/tmpHist->GetBinContent(idx);

          errRup.push_back( tmpHist->GetBinContent(idx) > 0. ? std::hypot(std::hypot(rIdUp,rFFup),rFFupM) : 0. );
          errRdn.push_back( tmpHist->GetBinContent(idx) > 0. ? std::hypot(std::hypot(rIdDn,rFFdn),rFFdnM) : 0. );
        }

        auto gr = new TGraphAsymmErrors(tmpHist->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kGray+2);
        gr->SetLineColor(kGray+2);
        gr->SetFillStyle(3004);
        padUp->cd();
        gr->Draw("2");

        auto rgr = new TGraphAsymmErrors(tmpHist->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
        rgr->SetFillColor(kGray+2);
        rgr->SetLineColor(kGray+2);
        rgr->SetFillStyle(3004);
        padDn->cd();
        rgr->Draw("2");
      }

      delete tmpHist;
    };
  };

  canvas_2->cd();
  auto aloader3P1F = HistLoader3P1F(datafile,WZfile,ZZfile);
  aloader3P1F.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF",postfix.Data());

  auto aloader3P1F1 = HistLoader3P1F(datafile1,WZfile1,ZZfile1);
  aloader3P1F1.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF","20UL16APV");

  auto aloader3P1F2 = HistLoader3P1F(datafile2,WZfile2,ZZfile2);
  aloader3P1F2.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF","20UL17");

  auto aloader3P1F3 = HistLoader3P1F(datafile3,WZfile3,ZZfile3);
  aloader3P1F3.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF","20UL18");

  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);

  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_llll_invM.pdf",p1);

/*  aloader3P1F.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll1ll2_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll1_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll1_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll2_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll2_dr.png",p1);*/

  // 4P0F CR

  class HistLoader4P0F : public HistLoaderBase {
  public:
    HistLoader4P0F(TFile* adatafile, TFile* aWZfile, TFile* aZZfile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile,aWZfile,aZZfile),
      sigFiles_(sigFiles) {}

    ~HistLoader4P0F() {
      if (dataHist_)
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FF2P2FHist_;
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* ZZ3P1FHist_ = nullptr;
    TH1D* ZZ4P0FHist_ = nullptr;
    TH1D* FF3P1FHist_ = nullptr;
    TH1D* FF2P2FHist_ = nullptr;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<std::map<std::string,SystVariation>> sigSyst_;

    TFile* datacard_ = nullptr;
    TDirectory* dir_ = nullptr;

  public:
    void preparecard(TString name, TString dirname) {
      datacard_ = new TFile(name,"RECREATE");
      dir_ = datacard_->mkdir(dirname);
      dir_->cd();
    }

    void close() {
      datacard_->Close();
      datacard_ = nullptr;
      dir_ = nullptr;
    }

    void add(const HistLoader4P0F& other) {
      this->dataHist_->Add(other.dataHist_);
      this->ZZ3P1FHist_->Add(other.ZZ3P1FHist_);
      this->ZZ4P0FHist_->Add(other.ZZ4P0FHist_);
      this->FF3P1FHist_->Add(other.FF3P1FHist_);
      this->FF2P2FHist_->Add(other.FF2P2FHist_);

      if (!syst_.empty()) {
        for (const auto& element : this->syst_) {
          this->syst_.at(element.first).up_->Add(other.syst_.at(element.first).up_);
          this->syst_.at(element.first).dn_->Add(other.syst_.at(element.first).dn_);
        }
      }

      if (!sigHist_.empty()) {
        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          this->sigHist_.at(idx)->Add(other.sigHist_.at(idx));

          if (!syst_.empty()) {
            for (const auto& element : this->sigSyst_.at(idx)) {
              this->sigSyst_.at(idx).at(element.first).up_->Add(other.sigSyst_.at(idx).at(element.first).up_);
              this->sigSyst_.at(idx).at(element.first).dn_->Add(other.sigSyst_.at(idx).at(element.first).dn_);
            }
          }
        }
      }
    }

    void load(const std::string& name, std::string anlyzrEra="") {
      if (dataHist_) {
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FF2P2FHist_;
      }

      if (!syst_.empty())
        syst_.clear();

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/4P0F_CR_")+name).c_str() )->Clone();
      ZZ3P1FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      ZZ4P0FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name).c_str() )->Clone();
      FF3P1FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      FF2P2FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2").c_str() )->Clone();

      const double lumi = retrieveLumi(anlyzrEra);
      sigHist_.clear();
      sigSyst_.clear();
      double sigLumi = 0.01;

      if ( TString(name).Contains("llll_invM") ) {
        TH1D* ZZ3P1FHist_idUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_idUp").c_str() )->Clone();
        TH1D* ZZ3P1FHist_idDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_idDn").c_str() )->Clone();
        syst_["ZZ3P1FHist_idEl"] = SystVariation(ZZ3P1FHist_idUp,ZZ3P1FHist_idDn);
        TH1D* ZZ4P0FHist_idUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_idDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_idEl"] = SystVariation(ZZ4P0FHist_idUp,ZZ4P0FHist_idDn);
        TH1D* ZZ3P1FHist_ffUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffUpE").c_str() )->Clone();
        TH1D* ZZ3P1FHist_ffDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffDnE").c_str() )->Clone();
        syst_["ZZ3P1FHist_ffEl"] = SystVariation(ZZ3P1FHist_ffUp,ZZ3P1FHist_ffDn);
        TH1D* FF3P1FHist_ffUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffUpE").c_str() )->Clone();
        TH1D* FF3P1FHist_ffDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffDnE").c_str() )->Clone();
        syst_["FF3P1FHist_ffEl"] = SystVariation(FF3P1FHist_ffUp,FF3P1FHist_ffDn);
        TH1D* FF2P2FHist_ffUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffUpE").c_str() )->Clone();
        TH1D* FF2P2FHist_ffDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffDnE").c_str() )->Clone();
        syst_["FF2P2FHist_ffEl"] = SystVariation(FF2P2FHist_ffUp,FF2P2FHist_ffDn);
        TH1D* ZZ3P1FHist_ffUpM = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffUpM").c_str() )->Clone();
        TH1D* ZZ3P1FHist_ffDnM = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffDnM").c_str() )->Clone();
        syst_["ZZ3P1FHist_ffMu"] = SystVariation(ZZ3P1FHist_ffUpM,ZZ3P1FHist_ffDnM);
        TH1D* FF3P1FHist_ffUpM = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffUpM").c_str() )->Clone();
        TH1D* FF3P1FHist_ffDnM = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffDnM").c_str() )->Clone();
        syst_["FF3P1FHist_ffMu"] = SystVariation(FF3P1FHist_ffUpM,FF3P1FHist_ffDnM);
        TH1D* FF2P2FHist_ffUpM = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffUpM").c_str() )->Clone();
        TH1D* FF2P2FHist_ffDnM = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffDnM").c_str() )->Clone();
        syst_["FF2P2FHist_ffMu"] = SystVariation(FF2P2FHist_ffUpM,FF2P2FHist_ffDnM);

        TH1D* ZZ3P1FHist_muScaleUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_altMuScale_xFF").c_str() )->Clone();
        TH1D* ZZ3P1FHist_muScaleDn = variateDn(ZZ3P1FHist_,ZZ3P1FHist_muScaleUp);
        syst_["ZZ3P1FHist_muScale"] = SystVariation(ZZ3P1FHist_muScaleUp,ZZ3P1FHist_muScaleDn);
        TH1D* FF3P1FHist_muScaleUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_altMuScale_xFF").c_str() )->Clone();
        TH1D* FF3P1FHist_muScaleDn = variateDn(FF3P1FHist_,FF3P1FHist_muScaleUp);
        syst_["FF3P1FHist_muScale"] = SystVariation(FF3P1FHist_muScaleUp,FF3P1FHist_muScaleDn);
        TH1D* FF2P2FHist_muScaleUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_altMuScale_xFF2").c_str() )->Clone();
        TH1D* FF2P2FHist_muScaleDn = variateDn(FF2P2FHist_,FF2P2FHist_muScaleUp);
        syst_["FF2P2FHist_muScale"] = SystVariation(FF2P2FHist_muScaleUp,FF2P2FHist_muScaleDn);
        TH1D* ZZ4P0FHist_muScaleUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuScale").c_str() )->Clone();
        TH1D* ZZ4P0FHist_muScaleDn = variateDn(ZZ4P0FHist_,ZZ4P0FHist_muScaleUp);
        syst_["ZZ4P0FHist_muScale"] = SystVariation(ZZ4P0FHist_muScaleUp,ZZ4P0FHist_muScaleDn);
        TH1D* ZZ4P0FHist_muSmearUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuSmear").c_str() )->Clone();
        TH1D* ZZ4P0FHist_muSmearDn = variateDn(ZZ4P0FHist_,ZZ4P0FHist_muSmearUp);
        syst_["ZZ4P0FHist_muSmear"] = SystVariation(ZZ4P0FHist_muSmearUp,ZZ4P0FHist_muSmearDn);
        TH1D* ZZ4P0FHist_muIdUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIdUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_muIdDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIdDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_muId"] = SystVariation(ZZ4P0FHist_muIdUp,ZZ4P0FHist_muIdDn);
        TH1D* ZZ4P0FHist_muIsoUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIsoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_muIsoDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIsoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_muIso"] = SystVariation(ZZ4P0FHist_muIsoUp,ZZ4P0FHist_muIsoDn);
        TH1D* ZZ4P0FHist_trigUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_trigDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_trig"] = SystVariation(ZZ4P0FHist_trigUp,ZZ4P0FHist_trigDn);
        TH1D* ZZ4P0FHist_muRecoUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muRecoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_muRecoDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muRecoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_muReco"] = SystVariation(ZZ4P0FHist_muRecoUp,ZZ4P0FHist_muRecoDn);
        TH1D* ZZ4P0FHist_elScaleUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elScaleUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_elScaleDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elScaleDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_elScale"] = SystVariation(ZZ4P0FHist_elScaleUp,ZZ4P0FHist_elScaleDn);
        TH1D* ZZ4P0FHist_elSigmaUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elSigmaUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_elSigmaDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elSigmaDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_elSigma"] = SystVariation(ZZ4P0FHist_elSigmaUp,ZZ4P0FHist_elSigmaDn);
        TH1D* ZZ4P0FHist_elRecoUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elRecoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_elRecoDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elRecoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_elReco"] = SystVariation(ZZ4P0FHist_elRecoUp_,ZZ4P0FHist_elRecoDn_);
        TH1D* ZZ4P0FHist_boostIsoUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_boostIsoDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_boostIso"] = SystVariation(ZZ4P0FHist_boostIsoUp,ZZ4P0FHist_boostIsoDn);
        TH1D* ZZ4P0FHist_PUrwgtUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_PUrwgtDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_PUrwgt"] = SystVariation(ZZ4P0FHist_PUrwgtUp,ZZ4P0FHist_PUrwgtDn);
        TH1D* ZZ4P0FHist_prefireUp = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_prefireDn = (TH1D*)ZZfile_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_prefire"] = SystVariation(ZZ4P0FHist_prefireUp,ZZ4P0FHist_prefireDn);

        for (const auto& element : syst_) {
          if (TString(element.first).Contains("ZZ")) {
            element.second.up_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
            element.second.dn_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
          }
        }

        for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
          sigHist_.push_back( (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name).c_str() )->Clone() );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          std::map<std::string,SystVariation> init;

          TH1D* sigHist_idUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
          TH1D* sigHist_idDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
          init["sigHist_idEl"] = SystVariation(sigHist_idUp,sigHist_idDn);
          TH1D* sigHist_muScaleUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuScale").c_str() )->Clone();
          TH1D* sigHist_muScaleDn = variateDn(sigHist_.back(),sigHist_muScaleUp);
          init["sigHist_muScale"] = SystVariation(sigHist_muScaleUp,sigHist_muScaleDn);
          TH1D* sigHist_muSmearUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuSmear").c_str() )->Clone();
          TH1D* sigHist_muSmearDn = variateDn(sigHist_.back(),sigHist_muSmearUp);
          init["sigHist_muSmear"] = SystVariation(sigHist_muSmearUp,sigHist_muSmearDn);
          TH1D* sigHist_muIdUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIdUp").c_str() )->Clone();
          TH1D* sigHist_muIdDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIdDn").c_str() )->Clone();
          init["sigHist_muId"] = SystVariation(sigHist_muIdUp,sigHist_muIdDn);
          TH1D* sigHist_muIsoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIsoUp").c_str() )->Clone();
          TH1D* sigHist_muIsoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muIsoDn").c_str() )->Clone();
          init["sigHist_muIso"] = SystVariation(sigHist_muIsoUp,sigHist_muIsoDn);
          TH1D* sigHist_trigUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigUp").c_str() )->Clone();
          TH1D* sigHist_trigDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigDn").c_str() )->Clone();
          init["sigHist_trig"] = SystVariation(sigHist_trigUp,sigHist_trigDn);
          TH1D* sigHist_muRecoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muRecoUp").c_str() )->Clone();
          TH1D* sigHist_muRecoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muRecoDn").c_str() )->Clone();
          init["sigHist_muReco"] = SystVariation(sigHist_muRecoUp,sigHist_muRecoDn);
          TH1D* sigHist_elScaleUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elScaleUp").c_str() )->Clone();
          TH1D* sigHist_elScaleDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elScaleDn").c_str() )->Clone();
          init["sigHist_elScale"] = SystVariation(sigHist_elScaleUp,sigHist_elScaleDn);
          TH1D* sigHist_elSigmaUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elSigmaUp").c_str() )->Clone();
          TH1D* sigHist_elSigmaDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elSigmaDn").c_str() )->Clone();
          init["sigHist_elSigma"] = SystVariation(sigHist_elSigmaUp,sigHist_elSigmaDn);
          TH1D* sigHist_elRecoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elRecoUp").c_str() )->Clone();
          TH1D* sigHist_elRecoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_elRecoDn").c_str() )->Clone();
          init["sigHist_elReco"] = SystVariation(sigHist_elRecoUp,sigHist_elRecoDn);
          TH1D* sigHist_boostIsoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoUp").c_str() )->Clone();
          TH1D* sigHist_boostIsoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoDn").c_str() )->Clone();
          init["sigHist_boostIso"] = SystVariation(sigHist_boostIsoUp,sigHist_boostIsoDn);
          TH1D* sigHist_PUrwgtUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtUp").c_str() )->Clone();
          TH1D* sigHist_PUrwgtDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtDn").c_str() )->Clone();
          init["sigHist_PUrwgt"] = SystVariation(sigHist_PUrwgtUp,sigHist_PUrwgtDn);
          TH1D* sigHist_prefireUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireUp").c_str() )->Clone();
          TH1D* sigHist_prefireDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireDn").c_str() )->Clone();
          init["sigHist_prefire"] = SystVariation(sigHist_prefireUp,sigHist_prefireDn);

          for (const auto& element : init) {
            element.second.up_->Scale( lumi*1000.*sigLumi / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );
            element.second.dn_->Scale( lumi*1000.*sigLumi / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );
          }

          sigHist_.at(idx)->Scale( lumi*1000.*sigLumi / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          sigSyst_.push_back(init);
        }
      }

      ZZ3P1FHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZ4P0FHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZ4P0FHist_->SetFillColor(TColor::GetColor("#5790fc"));
      ZZ4P0FHist_->SetLineWidth(0);
      FF3P1FHist_->SetFillColor(TColor::GetColor("#9c9ca1"));
      FF3P1FHist_->SetLineWidth(0);
      FF2P2FHist_->SetFillColor(TColor::GetColor("#9c9ca1"));
      FF2P2FHist_->SetLineWidth(0);
    };

    void compare(TPad* padUp, int rebin=1, TPad* padDn=nullptr) {
      if ( FF2P2FHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FF2P2FHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FF2P2FHist_->Rebin( FF2P2FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FF3P1FHist_->Rebin( FF3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        ZZ3P1FHist_->Rebin( ZZ3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( element.second.up_->GetNbinsX()/dataHist_->GetNbinsX() );
            element.second.dn_->Rebin( element.second.dn_->GetNbinsX()/dataHist_->GetNbinsX() );
          }
        }
      }

      if (!syst_.empty()) {
        for (const auto& element : syst_) {
          element.second.up_->Rebin( rebin );
          element.second.dn_->Rebin( rebin );
        }
      }

      ZZ3P1FHist_->Rebin(rebin);
      ZZ4P0FHist_->Rebin(rebin);
      FF3P1FHist_->Rebin(rebin);
      FF2P2FHist_->Rebin(rebin);
      dataHist_->Rebin(rebin);

      auto subtract = [] (TH1D* data3p1f, TH1D* zz3p1f) -> TH1D* {
        TH1D* subtracted = (TH1D*)data3p1f->Clone();

        for (unsigned ibin = 0; ibin < subtracted->GetNbinsX()+2; ibin++) {
          double val = data3p1f->GetBinContent(ibin) - zz3p1f->GetBinContent(ibin);
          double err = std::hypot( data3p1f->GetBinError(ibin), zz3p1f->GetBinError(ibin) );
          subtracted->SetBinContent(ibin,val);
          subtracted->SetBinError(ibin,err);
        }

        return subtracted;
      };

      auto truncateNegativeBin = [] (TH1D* ahist) {
        for (unsigned ibin = 0; ibin < ahist->GetNbinsX()+2; ibin++) {
          double val = ahist->GetBinContent(ibin);

          if (val < 0.) {
            ahist->SetBinContent(ibin,0.);
            ahist->SetBinError(ibin,0.);
          }
        }
      };

      TH1D* FF3P1Fsubtracted = subtract(FF3P1FHist_,ZZ3P1FHist_);
      TH1D* FFsubtractedFinal = subtract(FF3P1Fsubtracted,FF2P2FHist_);

      for (unsigned ibin = 0; ibin < FFsubtractedFinal->GetNbinsX()+2; ibin++) {
        double val = FFsubtractedFinal->GetBinContent(ibin);
        double err = FFsubtractedFinal->GetBinError(ibin);

        if ( val < 0.) {
          ZZ4P0FHist_->SetBinContent( ibin, std::max(ZZ4P0FHist_->GetBinContent(ibin) + val, 0.) );
          ZZ4P0FHist_->SetBinError( ibin, std::hypot( err, ZZ4P0FHist_->GetBinError(ibin) ) );
          FFsubtractedFinal->SetBinContent(ibin,0.);
          FFsubtractedFinal->SetBinError(ibin,0.);

          for (const auto& element : syst_) {
            if (element.first.find("ZZ4P0F")!=std::string::npos) {
              element.second.up_->SetBinContent( ibin, std::max(element.second.up_->GetBinContent(ibin) + val, 0.) );
              element.second.dn_->SetBinContent( ibin, std::max(element.second.dn_->GetBinContent(ibin) + val, 0.) );
              element.second.up_->SetBinError( ibin, std::hypot( err, element.second.up_->GetBinError(ibin) ) ); 
              element.second.dn_->SetBinError( ibin, std::hypot( err, element.second.dn_->GetBinError(ibin) ) ); 
            }
          }
        }
      }

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(FFsubtractedFinal);
      stack->Add(ZZ4P0FHist_);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));

      padUp->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        dataHist_->GetXaxis()->SetRangeUser(0.,1000.);
      }

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMinimum(0.001);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.65,0.65,0.95,0.9);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(FFsubtractedFinal,"Nonprompt");
      legend->AddEntry(ZZ4P0FHist_,"ZZ");

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);
          sigHist_.at(idx)->Draw("hist&same");

          for (const auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_->Rebin(rebin);
            sigSyst_.at(idx).at(element.first).dn_->Rebin(rebin);
          }
        }

        //legend->AddEntry(sigHist_.at(3),"H250/750/2000A1");
        legend->AddEntry(sigHist_.at(0),"X750Y100");
      }

      legend->Draw();

      if (padDn) {
        TH1D* ratio = (TH1D*)dataHist_->Clone();
        ratio->SetStats(0);
        ratio->SetTitle("");
        ratio->Divide(tmpHist);
        ratio->GetYaxis()->SetTitle("Obs/Exp");
        ratio->GetYaxis()->SetTitleSize(0.1);
        ratio->GetYaxis()->SetTitleOffset(0.4);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetXaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetLabelOffset(0.005);
        ratio->GetYaxis()->SetRangeUser(0.2,1.8);
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetTitleOffset(0.75);
        ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
        ratio->SetLineColor(kBlack);

        padDn->cd();
        ratio->Draw("E1");
        padUp->cd();
      }

      if (!syst_.empty()) {
        TH1D* FF3P1Fsubtracted_idUp = subtract(subtract(FF3P1FHist_,syst_.at("ZZ3P1FHist_idEl").up_),FF2P2FHist_);
        TH1D* FF3P1Fsubtracted_idDn = subtract(subtract(FF3P1FHist_,syst_.at("ZZ3P1FHist_idEl").dn_),FF2P2FHist_);
        truncateNegativeBin(FF3P1Fsubtracted_idUp);
        truncateNegativeBin(FF3P1Fsubtracted_idDn);
        TH1D* idUp = (TH1D*)syst_.at("ZZ4P0FHist_idEl").up_->Clone();
        TH1D* idDn = (TH1D*)syst_.at("ZZ4P0FHist_idEl").dn_->Clone();
        idUp->Add(FF3P1Fsubtracted_idUp);
        idDn->Add(FF3P1Fsubtracted_idDn);

        TH1D* FF3P1Fsubtracted_ffUp = subtract(subtract(syst_.at("FF3P1FHist_ffEl").up_,syst_.at("ZZ3P1FHist_ffEl").up_),syst_.at("FF2P2FHist_ffEl").up_);
        TH1D* FF3P1Fsubtracted_ffDn = subtract(subtract(syst_.at("FF3P1FHist_ffEl").dn_,syst_.at("ZZ3P1FHist_ffEl").dn_),syst_.at("FF2P2FHist_ffEl").dn_);
        truncateNegativeBin(FF3P1Fsubtracted_ffUp);
        truncateNegativeBin(FF3P1Fsubtracted_ffDn);
        TH1D* ffUpAdded = (TH1D*)FF3P1Fsubtracted_ffUp->Clone();
        TH1D* ffDnAdded = (TH1D*)FF3P1Fsubtracted_ffDn->Clone();
        ffUpAdded->Add(ZZ4P0FHist_);
        ffDnAdded->Add(ZZ4P0FHist_);

        TH1D* FF3P1Fsubtracted_ffUpM = subtract(subtract(syst_.at("FF3P1FHist_ffMu").up_,syst_.at("ZZ3P1FHist_ffMu").up_),syst_.at("FF2P2FHist_ffMu").up_);
        TH1D* FF3P1Fsubtracted_ffDnM = subtract(subtract(syst_.at("FF3P1FHist_ffMu").dn_,syst_.at("ZZ3P1FHist_ffMu").dn_),syst_.at("FF2P2FHist_ffMu").dn_);
        truncateNegativeBin(FF3P1Fsubtracted_ffUpM);
        truncateNegativeBin(FF3P1Fsubtracted_ffDnM);
        TH1D* ffUpMAdded = (TH1D*)FF3P1Fsubtracted_ffUpM->Clone();
        TH1D* ffDnMAdded = (TH1D*)FF3P1Fsubtracted_ffDnM->Clone();
        ffUpMAdded->Add(ZZ4P0FHist_);
        ffDnMAdded->Add(ZZ4P0FHist_);

        TH1D* FF3P1Fsubtracted_muScaleUp = subtract(subtract(syst_.at("FF3P1FHist_muScale").up_,syst_.at("ZZ3P1FHist_muScale").up_),syst_.at("FF2P2FHist_muScale").up_);
        TH1D* FF3P1Fsubtracted_muScaleDn = subtract(subtract(syst_.at("FF3P1FHist_muScale").dn_,syst_.at("ZZ3P1FHist_muScale").dn_),syst_.at("FF2P2FHist_muScale").dn_);
        truncateNegativeBin(FF3P1Fsubtracted_muScaleUp);
        truncateNegativeBin(FF3P1Fsubtracted_muScaleDn);
        TH1D* muScaleUpAdded = (TH1D*)FF3P1Fsubtracted_muScaleUp->Clone();
        TH1D* muScaleDnAdded = (TH1D*)FF3P1Fsubtracted_muScaleDn->Clone();
        muScaleUpAdded->Add(syst_.at("ZZ4P0FHist_muScale").up_);
        muScaleDnAdded->Add(syst_.at("ZZ4P0FHist_muScale").dn_);

        TH1D* muSmearUp = (TH1D*)syst_.at("ZZ4P0FHist_muSmear").up_->Clone();
        TH1D* muSmearDn = (TH1D*)syst_.at("ZZ4P0FHist_muSmear").dn_->Clone();
        muSmearUp->Add(FFsubtractedFinal);
        muSmearDn->Add(FFsubtractedFinal);

        TH1D* muIdUp = (TH1D*)syst_.at("ZZ4P0FHist_muId").up_->Clone();
        TH1D* muIdDn = (TH1D*)syst_.at("ZZ4P0FHist_muId").dn_->Clone();
        muIdUp->Add(FFsubtractedFinal);
        muIdDn->Add(FFsubtractedFinal);

        TH1D* muIsoUp = (TH1D*)syst_.at("ZZ4P0FHist_muIso").up_->Clone();
        TH1D* muIsoDn = (TH1D*)syst_.at("ZZ4P0FHist_muIso").dn_->Clone();
        muIsoUp->Add(FFsubtractedFinal);
        muIsoDn->Add(FFsubtractedFinal);

        TH1D* muTrigUp = (TH1D*)syst_.at("ZZ4P0FHist_trig").up_->Clone();
        TH1D* muTrigDn = (TH1D*)syst_.at("ZZ4P0FHist_trig").dn_->Clone();
        muTrigUp->Add(FFsubtractedFinal);
        muTrigDn->Add(FFsubtractedFinal);

        TH1D* muRecoUp = (TH1D*)syst_.at("ZZ4P0FHist_muReco").up_->Clone();
        TH1D* muRecoDn = (TH1D*)syst_.at("ZZ4P0FHist_muReco").dn_->Clone();
        muRecoUp->Add(FFsubtractedFinal);
        muRecoDn->Add(FFsubtractedFinal);

        TH1D* elScaleUp = (TH1D*)syst_.at("ZZ4P0FHist_elScale").up_->Clone();
        TH1D* elScaleDn = (TH1D*)syst_.at("ZZ4P0FHist_elScale").dn_->Clone();
        elScaleUp->Add(FFsubtractedFinal);
        elScaleDn->Add(FFsubtractedFinal);

        TH1D* elSigmaUp = (TH1D*)syst_.at("ZZ4P0FHist_elSigma").up_->Clone();
        TH1D* elSigmaDn = (TH1D*)syst_.at("ZZ4P0FHist_elSigma").dn_->Clone();
        elSigmaUp->Add(FFsubtractedFinal);
        elSigmaDn->Add(FFsubtractedFinal);

        TH1D* boostIsoUp = (TH1D*)syst_.at("ZZ4P0FHist_boostIso").up_->Clone();
        TH1D* boostIsoDn = (TH1D*)syst_.at("ZZ4P0FHist_boostIso").dn_->Clone();
        boostIsoUp->Add(FFsubtractedFinal);
        boostIsoDn->Add(FFsubtractedFinal);

        TH1D* elRecoUpAdded = (TH1D*)syst_.at("ZZ4P0FHist_elReco").up_->Clone();
        TH1D* elRecoDnAdded = (TH1D*)syst_.at("ZZ4P0FHist_elReco").dn_->Clone();
        elRecoUpAdded->Add(FFsubtractedFinal);
        elRecoDnAdded->Add(FFsubtractedFinal);

        TH1D* PUrwgtUp = (TH1D*)syst_.at("ZZ4P0FHist_PUrwgt").up_->Clone();
        TH1D* PUrwgtDn = (TH1D*)syst_.at("ZZ4P0FHist_PUrwgt").dn_->Clone();
        PUrwgtUp->Add(FFsubtractedFinal);
        PUrwgtDn->Add(FFsubtractedFinal);

        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRup, errRdn;

        auto sq = [] (const double val) { return val*val; };

        for (unsigned idx = 1; idx <= tmpHist->GetNbinsX(); idx++) {
          x0.push_back(tmpHist->GetBinCenter(idx));
          y0.push_back(tmpHist->GetBinContent(idx));
          errx.push_back(tmpHist->GetBinWidth(idx)/2.);

          double valIdUp = idUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valIdDn = tmpHist->GetBinContent(idx) - idDn->GetBinContent(idx);
          double valFFup = ffUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdn = tmpHist->GetBinContent(idx) - ffDnAdded->GetBinContent(idx);
          double valFFupM = ffUpMAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdnM = tmpHist->GetBinContent(idx) - ffDnMAdded->GetBinContent(idx);
          double valMuScaleUp = muScaleUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valMuScaleDn = tmpHist->GetBinContent(idx) - muScaleDnAdded->GetBinContent(idx);
          double valMuSmearUp = muSmearUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valMuSmearDn = tmpHist->GetBinContent(idx) - muSmearDn->GetBinContent(idx);
          double valMuIdUp = muIdUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valMuIdDn = tmpHist->GetBinContent(idx) - muIdDn->GetBinContent(idx);
          double valMuIsoUp = muIsoUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valMuIsoDn = tmpHist->GetBinContent(idx) - muIsoDn->GetBinContent(idx);
          double valTrigUp = muTrigUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valTrigDn = tmpHist->GetBinContent(idx) - muTrigDn->GetBinContent(idx);
          double valMuRecoUp = muRecoUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valMuRecoDn = tmpHist->GetBinContent(idx) - muRecoDn->GetBinContent(idx);
          double valElScaleUp = elScaleUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valElScaleDn = tmpHist->GetBinContent(idx) - elScaleDn->GetBinContent(idx);
          double valElSigmaUp = elSigmaUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valElSigmaDn = tmpHist->GetBinContent(idx) - elSigmaDn->GetBinContent(idx);
          double valElRecoUp = elRecoUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valElRecoDn = tmpHist->GetBinContent(idx) - elRecoDnAdded->GetBinContent(idx);
          double valboostIsoUp = boostIsoUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valboostIsoDn = tmpHist->GetBinContent(idx) - boostIsoDn->GetBinContent(idx);
          double valPUrwgtUp = PUrwgtUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valPUrwgtDn = tmpHist->GetBinContent(idx) - PUrwgtDn->GetBinContent(idx);
          double ZZnorm = 0.1*ZZ4P0FHist_->GetBinContent(idx);

          erryUp.push_back( std::sqrt( sq(valIdUp) + sq(valFFup) + sq(valFFupM) + sq(valMuScaleUp)
                                       + sq(valMuSmearUp) + sq(valMuIdUp) + sq(valMuIsoUp) + sq(valTrigUp)
                                       + sq(valMuRecoUp) + sq(valElScaleUp) + sq(valElSigmaUp) + sq(valElRecoUp)
                                       + sq(valboostIsoUp) + sq(valPUrwgtUp) + sq(ZZnorm) ) );
          erryDn.push_back( std::sqrt( sq(valIdDn) + sq(valFFdn) + sq(valFFdnM) + sq(valMuScaleDn)
                                       + sq(valMuSmearDn) + sq(valMuIdDn) + sq(valMuIsoDn) + sq(valTrigDn)
                                       + sq(valMuRecoDn) + sq(valElScaleDn) + sq(valElSigmaDn) + sq(valElRecoDn)
                                       + sq(valboostIsoDn) + sq(valPUrwgtDn) + sq(ZZnorm) ) );

          r0.push_back(1.);

          double rIdUp = valIdUp/tmpHist->GetBinContent(idx);
          double rIdDn = valIdDn/tmpHist->GetBinContent(idx);
          double rFFup = valFFup/tmpHist->GetBinContent(idx);
          double rFFdn = valFFdn/tmpHist->GetBinContent(idx);
          double rFFupM = valFFupM/tmpHist->GetBinContent(idx);
          double rFFdnM = valFFdnM/tmpHist->GetBinContent(idx);
          double rMuScaleUp = valMuScaleUp/tmpHist->GetBinContent(idx);
          double rMuScaleDn = valMuScaleDn/tmpHist->GetBinContent(idx);
          double rMuSmearUp = valMuSmearUp/tmpHist->GetBinContent(idx);
          double rMuSmearDn = valMuSmearDn/tmpHist->GetBinContent(idx);
          double rMuIdUp = valMuIdUp/tmpHist->GetBinContent(idx);
          double rMuIdDn = valMuIdDn/tmpHist->GetBinContent(idx);
          double rMuIsoUp = valMuIsoUp/tmpHist->GetBinContent(idx);
          double rMuIsoDn = valMuIsoDn/tmpHist->GetBinContent(idx);
          double rTrigUp = valTrigUp/tmpHist->GetBinContent(idx);
          double rTrigDn = valTrigDn/tmpHist->GetBinContent(idx);
          double rMuRecoUp = valMuRecoUp/tmpHist->GetBinContent(idx);
          double rMuRecoDn = valMuRecoDn/tmpHist->GetBinContent(idx);
          double rElScaleUp = valElScaleUp/tmpHist->GetBinContent(idx);
          double rElScaleDn = valElScaleDn/tmpHist->GetBinContent(idx);
          double rElSigmaUp = valElSigmaUp/tmpHist->GetBinContent(idx);
          double rElSigmaDn = valElSigmaDn/tmpHist->GetBinContent(idx);
          double rElRecoUp = valElRecoUp/tmpHist->GetBinContent(idx);
          double rElRecoDn = valElRecoDn/tmpHist->GetBinContent(idx);
          double rBoostIsoUp = valboostIsoUp/tmpHist->GetBinContent(idx);
          double rBoostIsoDn = valboostIsoDn/tmpHist->GetBinContent(idx);
          double rPUrwgtUp = valPUrwgtUp/tmpHist->GetBinContent(idx);
          double rPUrwgtDn = valPUrwgtDn/tmpHist->GetBinContent(idx);
          double rZZnorm = ZZnorm/tmpHist->GetBinContent(idx);

          double rUp = std::sqrt( sq(rIdUp) + sq(rFFup) + sq(rFFupM) + sq(rMuScaleUp)
                                     + sq(rMuSmearUp) + sq(rMuIdUp) + sq(rMuIsoUp) + sq(rTrigUp)
                                     + sq(rMuRecoUp) + sq(rElScaleUp) + sq(rElSigmaUp) + sq(rElRecoUp)
                                     + sq(rBoostIsoUp) + sq(rPUrwgtUp) + sq(rZZnorm) );
          double rDn = std::sqrt( sq(rIdDn) + sq(rFFdn) + sq(rFFdnM) + sq(rMuScaleDn)
                                     + sq(rMuSmearDn) + sq(rMuIdDn) + sq(rMuIsoDn) + sq(rTrigDn)
                                     + sq(rMuRecoDn) + sq(rElScaleDn) + sq(rElSigmaDn) + sq(rElRecoDn)
                                     + sq(rBoostIsoDn) + sq(rPUrwgtDn) + sq(rZZnorm) );

          errRup.push_back( tmpHist->GetBinContent(idx) > 0. ? rUp : 0. );
          errRdn.push_back( tmpHist->GetBinContent(idx) > 0. ? rDn : 0. );
        }

        auto gr = new TGraphAsymmErrors(tmpHist->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kGray+3);
        gr->SetLineColor(kGray+3);
        gr->SetFillStyle(3004);
        gr->Draw("2");

        if (padDn) {
          auto rgr = new TGraphAsymmErrors(tmpHist->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
          rgr->SetFillColor(kGray+2);
          rgr->SetLineColor(kGray+2);
          rgr->SetFillStyle(3004);
          padDn->cd();
          rgr->Draw("2");
        }

        if (dir_) {
          dir_->WriteTObject(dataHist_,"data_obs");
          dir_->WriteTObject(FFsubtractedFinal,"Nonprompt");
          dir_->WriteTObject(ZZ4P0FHist_,"ZZ");

          dir_->WriteTObject(FF3P1Fsubtracted_ffUp,"Nonprompt_resolvedEleFakeFactorUp");
          dir_->WriteTObject(FF3P1Fsubtracted_ffDn,"Nonprompt_resolvedEleFakeFactorDown");
          dir_->WriteTObject(FF3P1Fsubtracted_ffUpM,"Nonprompt_resolvedMuFakeFactorUp");
          dir_->WriteTObject(FF3P1Fsubtracted_ffDnM,"Nonprompt_resolvedMuFakeFactorDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_idEl").up_,"ZZ_modHeepIdUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_idEl").dn_,"ZZ_modHeepIdDown");
          dir_->WriteTObject(FF3P1Fsubtracted_idUp,"Nonprompt_modHeepIdUp");
          dir_->WriteTObject(FF3P1Fsubtracted_idDn,"Nonprompt_modHeepIdDown");
          dir_->WriteTObject(FF3P1Fsubtracted_muScaleUp,"Nonprompt_muMomentumScaleUp");
          dir_->WriteTObject(FF3P1Fsubtracted_muScaleDn,"Nonprompt_muMomentumScaleDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muScale").up_,"ZZ_muMomentumScaleUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muScale").dn_,"ZZ_muMomentumScaleDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muSmear").up_,"ZZ_muMomentumSmearUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muSmear").dn_,"ZZ_muMomentumSmearDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muId").up_,"ZZ_highPtIdUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muId").dn_,"ZZ_highPtIdDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muIso").up_,"ZZ_muLooseIsoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muIso").dn_,"ZZ_muLooseIsoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_trig").up_,"ZZ_muTrigUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_trig").dn_,"ZZ_muTrigDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muReco").up_,"ZZ_muRecoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_muReco").dn_,"ZZ_muRecoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_elScale").up_,"ZZ_elEnergyScaleUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_elScale").dn_,"ZZ_elEnergyScaleDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_elSigma").up_,"ZZ_elEnergySigmaUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_elSigma").dn_,"ZZ_elEnergySigmaDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_boostIso").up_,"ZZ_muBoostIsoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_boostIso").dn_,"ZZ_muBoostIsoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_elReco").up_,"ZZ_elRecoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_elReco").dn_,"ZZ_elRecoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_PUrwgt").up_,"ZZ_PUrwgtUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_PUrwgt").dn_,"ZZ_PUrwgtDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_prefire").up_,"ZZ_prefireUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_prefire").dn_,"ZZ_prefireDown");

          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_idEl").up_,sigFiles_.at(idx).name_+"_modHeepIdUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_idEl").dn_,sigFiles_.at(idx).name_+"_modHeepIdDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muScale").up_,sigFiles_.at(idx).name_+"_muMomentumScaleUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muScale").dn_,sigFiles_.at(idx).name_+"_muMomentumScaleDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muSmear").up_,sigFiles_.at(idx).name_+"_muMomentumSmearUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muSmear").dn_,sigFiles_.at(idx).name_+"_muMomentumSmearDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muId").up_,sigFiles_.at(idx).name_+"_highPtIdUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muId").dn_,sigFiles_.at(idx).name_+"_highPtIdDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muIso").up_,sigFiles_.at(idx).name_+"_muLooseIsoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muIso").dn_,sigFiles_.at(idx).name_+"_muLooseIsoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_trig").up_,sigFiles_.at(idx).name_+"_muTrigUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_trig").dn_,sigFiles_.at(idx).name_+"_muTrigDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muReco").up_,sigFiles_.at(idx).name_+"_muRecoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_muReco").dn_,sigFiles_.at(idx).name_+"_muRecoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_elScale").up_,sigFiles_.at(idx).name_+"_elEnergyScaleUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_elScale").dn_,sigFiles_.at(idx).name_+"_elEnergyScaleDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_elSigma").up_,sigFiles_.at(idx).name_+"_elEnergySigmaUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_elSigma").dn_,sigFiles_.at(idx).name_+"_elEnergySigmaDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_boostIso").up_,sigFiles_.at(idx).name_+"_muBoostIsoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_boostIso").dn_,sigFiles_.at(idx).name_+"_muBoostIsoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_elReco").up_,sigFiles_.at(idx).name_+"_elRecoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_elReco").dn_,sigFiles_.at(idx).name_+"_elRecoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_PUrwgt").up_,sigFiles_.at(idx).name_+"_PUrwgtUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_PUrwgt").dn_,sigFiles_.at(idx).name_+"_PUrwgtDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_prefire").up_,sigFiles_.at(idx).name_+"_prefireUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_prefire").dn_,sigFiles_.at(idx).name_+"_prefireDown");
          }
        }
      }
    };
  };

  auto aloader4P0F = HistLoader4P0F(datafile,WZfile,ZZfile,sigsamples);
  aloader4P0F.load("llll_invM",postfix.Data());

  auto aloader4P0F1 = HistLoader4P0F(datafile1,WZfile1,ZZfile1,sigsamples1);
  aloader4P0F1.load("llll_invM","20UL16APV");

  auto aloader4P0F2 = HistLoader4P0F(datafile2,WZfile2,ZZfile2,sigsamples2);
  aloader4P0F2.load("llll_invM","20UL17");

  auto aloader4P0F3 = HistLoader4P0F(datafile3,WZfile3,ZZfile3,sigsamples3);
  aloader4P0F3.load("llll_invM","20UL18");

  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);

  //aloader4P0F.preparecard("REMuFF_"+era+"_datacard.root","resolvedEMu");
  aloader4P0F.compare(p1,2,p2); // 2
  SaveAs(canvas_2,"REMuFF_4P0F_CR_llll_invM_zoomed.pdf",p1);
  //aloader4P0F.close();

  p1->SetLogy(0);
  canvas_1->SetLogy(0);

/*  aloader4P0F.load("ll1ll2_dr",postfix.Data());
  aloader4P0F1.load("ll1ll2_dr","20UL16APV");
  aloader4P0F2.load("ll1ll2_dr","20UL17");
  aloader4P0F3.load("ll1ll2_dr","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll1ll2_dr.png");

  aloader4P0F.load("ll1_invM",postfix.Data());
  aloader4P0F1.load("ll1_invM","20UL16APV");
  aloader4P0F2.load("ll1_invM","20UL17");
  aloader4P0F3.load("ll1_invM","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll1_invM.png");

  aloader4P0F.load("ll1_dr",postfix.Data());
  aloader4P0F1.load("ll1_dr","20UL16APV");
  aloader4P0F2.load("ll1_dr","20UL17");
  aloader4P0F3.load("ll1_dr","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll1_dr.png");

  aloader4P0F.load("ll2_invM",postfix.Data());
  aloader4P0F1.load("ll2_invM","20UL16APV");
  aloader4P0F2.load("ll2_invM","20UL17");
  aloader4P0F3.load("ll2_invM","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll2_invM.png");

  aloader4P0F.load("ll2_dr",postfix.Data());
  aloader4P0F1.load("ll2_dr","20UL16APV");
  aloader4P0F2.load("ll2_dr","20UL17");
  aloader4P0F3.load("ll2_dr","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll2_dr.png");*/

  return;
}
