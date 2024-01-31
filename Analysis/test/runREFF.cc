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

void runREFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static constexpr double WZxsec_ = 5.213; // 0.65*62.78;
  static constexpr double ZZxsec_ = 13.81;
  static double valLumi = retrieveLumi(era.Data());
  static TString postfix = era;

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
  } else {
    std::cout << "check era..." << std::endl;
  }

  static TString anlyzrMC = "resolvedEleCRanalyzer"+postfix;
  static TString anlyzrData = "resolvedEleCRanalyzerData";

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

  TString firstEra = era=="run2" ? "20UL16" : era;

  TFile* datafile = new TFile("EleAnalyzer_"+firstEra+"_data.root","READ");
  TFile* WZfile = new TFile("EleAnalyzer_"+firstEra+"_WZFXFX.root","READ");
  TFile* ZZfile = new TFile("EleAnalyzer_"+firstEra+"_ZZ.root","READ");

  TFile *datafile1, *WZfile1, *ZZfile1;
  TFile *datafile2, *WZfile2, *ZZfile2;
  TFile *datafile3, *WZfile3, *ZZfile3;

  if (era.Contains("run2")) {
    datafile1 = new TFile("EleAnalyzer_20UL16APV_data.root","READ");
    WZfile1 = new TFile("EleAnalyzer_20UL16APV_WZ.root","READ");
    ZZfile1 = new TFile("EleAnalyzer_20UL16APV_ZZ.root","READ");

    datafile2 = new TFile("EleAnalyzer_20UL17_data.root","READ");
    WZfile2 = new TFile("EleAnalyzer_20UL17_WZ.root","READ");
    ZZfile2 = new TFile("EleAnalyzer_20UL17_ZZ.root","READ");

    datafile3 = new TFile("EleAnalyzer_20UL18_data.root","READ");
    WZfile3 = new TFile("EleAnalyzer_20UL18_WZ.root","READ");
    ZZfile3 = new TFile("EleAnalyzer_20UL18_ZZ.root","READ");
  }

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  std::vector<SigSample> sigsamples = {
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A750.root","READ"),"H2000A750")
  };

  std::vector<SigSample> sigsamples1, sigsamples2, sigsamples3;

  if (era.Contains("run2")) {
    sigsamples1 = {
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples2 = {
      SigSample(new TFile("EleAnalyzer_20UL17_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples3 = {
      SigSample(new TFile("EleAnalyzer_20UL18_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A750.root","READ"),"H2000A750")
    };
  }

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
    canvas->cd();

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
        delete WZHist_idUp_, ZZHist_idUp_, WZHist_idDn_, ZZHist_idDn_, FFHist_ffUp_, FFHist_ffDn_;
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
      }
    }

    void load(const std::string& nameNum, const std::string& name, std::string anlyzrEra="") {
      if (dataHist_)
        delete dataHist_, WZHist_, ZZHist_, FFHist_;

      if (FFHist_ffUp_) {
        delete WZHist_idUp_, ZZHist_idUp_, WZHist_idDn_, ZZHist_idDn_, FFHist_ffUp_, FFHist_ffDn_;

        FFHist_ffUp_ = nullptr;
      }

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum).c_str() )->Clone();
      WZHist_ = (TH1D*)WZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/")+nameNum).c_str() )->Clone();
      ZZHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/")+nameNum).c_str() )->Clone();
      FFHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name).c_str() )->Clone();

      const double lumi = retrieveLumi(anlyzrEra);

      if (TString(nameNum).Contains("llll_invM")) {
        WZHist_idUp_ = (TH1D*)WZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/")+nameNum+"_idUp").c_str() )->Clone();
        ZZHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/")+nameNum+"_idUp").c_str() )->Clone();
        WZHist_idDn_ = (TH1D*)WZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/")+nameNum+"_idDn").c_str() )->Clone();
        ZZHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/")+nameNum+"_idDn").c_str() )->Clone();
        FFHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffUp").c_str() )->Clone();
        FFHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffDn").c_str() )->Clone();

        WZHist_idUp_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        WZHist_idDn_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZHist_idUp_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZHist_idDn_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      }

      WZHist_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) ); // 13.81
      WZHist_->SetFillColor(kViolet+1);
      ZZHist_->SetFillColor(kBlue-4);
      WZHist_->SetLineWidth(0);
      ZZHist_->SetLineWidth(0);
      FFHist_->SetFillColor(33);
      FFHist_->SetLineWidth(0);
    };

    void compareAndRatio(TPad* padUp, TPad* padDn, int rebin=1) {
      if ( FFHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FFHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (FFHist_ffUp_) {
          FFHist_ffUp_->Rebin( FFHist_ffUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          FFHist_ffDn_->Rebin( FFHist_ffDn_->GetNbinsX()/dataHist_->GetNbinsX() );
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
        dataHist_->SetMaximum(1.*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));
        dataHist_->GetXaxis()->SetRangeUser(0.,500.);
      }

      // TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
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

          erryUp.push_back(std::hypot(valIdUp,valFFup));
          erryDn.push_back(std::hypot(valIdDn,valFFdn));

          r0.push_back(1.);

          double rIdUp = valIdUp/tmpHist->GetBinContent(idx);
          double rIdDn = valIdDn/tmpHist->GetBinContent(idx);
          double rFFup = valFFup/tmpHist->GetBinContent(idx);
          double rFFdn = valFFdn/tmpHist->GetBinContent(idx);

          errRup.push_back(std::hypot(rIdUp,rFFup));
          errRdn.push_back(std::hypot(rIdDn,rFFdn));
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
  auto aloader3P1F1 = HistLoader3P1F(datafile1,WZfile1,ZZfile1);
  auto aloader3P1F2 = HistLoader3P1F(datafile2,WZfile2,ZZfile2);
  auto aloader3P1F3 = HistLoader3P1F(datafile3,WZfile3,ZZfile3);

  aloader3P1F.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF",postfix.Data());

  if (era.Contains("run2")) {
    aloader3P1F1.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF","20UL16APV");
    aloader3P1F2.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF","20UL17");
    aloader3P1F3.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF","20UL18");
    aloader3P1F.add(aloader3P1F1);
    aloader3P1F.add(aloader3P1F2);
    aloader3P1F.add(aloader3P1F3);
  }

  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REFF_3P1F_CR_llll_invM.eps",p1);

  aloader3P1F.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF",postfix.Data());

  if (era.Contains("run2")) {
    aloader3P1F1.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL16APV");
    aloader3P1F2.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL17");
    aloader3P1F3.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL18");
    aloader3P1F.add(aloader3P1F1);
    aloader3P1F.add(aloader3P1F2);
    aloader3P1F.add(aloader3P1F3);
  }

  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll1ll2_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF",postfix.Data());

  if (era.Contains("run2")) {
    aloader3P1F1.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL16APV");
    aloader3P1F2.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL17");
    aloader3P1F3.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL18");
    aloader3P1F.add(aloader3P1F1);
    aloader3P1F.add(aloader3P1F2);
    aloader3P1F.add(aloader3P1F3);
  }

  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll1_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF",postfix.Data());

  if (era.Contains("run2")) {
    aloader3P1F1.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL16APV");
    aloader3P1F2.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL17");
    aloader3P1F3.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL18");
    aloader3P1F.add(aloader3P1F1);
    aloader3P1F.add(aloader3P1F2);
    aloader3P1F.add(aloader3P1F3);
  }

  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll1_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF",postfix.Data());

  if (era.Contains("run2")) {
    aloader3P1F1.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL16APV");
    aloader3P1F2.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL17");
    aloader3P1F3.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL18");
    aloader3P1F.add(aloader3P1F1);
    aloader3P1F.add(aloader3P1F2);
    aloader3P1F.add(aloader3P1F3);
  }

  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll2_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF",postfix.Data());

  if (era.Contains("run2")) {
    aloader3P1F1.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL16APV");
    aloader3P1F2.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL17");
    aloader3P1F3.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL18");
    aloader3P1F.add(aloader3P1F1);
    aloader3P1F.add(aloader3P1F2);
    aloader3P1F.add(aloader3P1F3);
  }

  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll2_dr.png",p1);

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
      ZZ3P1FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      ZZ4P0FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name).c_str() )->Clone();
      FF3P1FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      FF2P2FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2").c_str() )->Clone();

      const double lumi = retrieveLumi(anlyzrEra);
      sigHist_.clear();
      sigSyst_.clear();
      double sigLumi = 0.001;

      for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
        sigHist_.push_back( (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name).c_str() )->Clone() );
        sigHist_.back()->Scale( lumi*1000.*sigLumi / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) ); // h_sum4E
        sigHist_.back()->SetLineWidth(2);
        sigHist_.back()->SetLineColor(kRed);
      }

      ZZ3P1FHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) ); // 13.81
      ZZ4P0FHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZ4P0FHist_->SetFillColor(kBlue-4);
      ZZ4P0FHist_->SetLineWidth(0);
      FF3P1FHist_->SetFillColor(33);
      FF3P1FHist_->SetLineWidth(0);
      FF2P2FHist_->SetFillColor(33);
      FF2P2FHist_->SetLineWidth(0);

      if ( TString(name).Contains("llll_invM") ) {
        TH1D* ZZ3P1FHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_idUp").c_str() )->Clone();
        TH1D* ZZ3P1FHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_idDn").c_str() )->Clone();
        syst_["ZZ3P1FHist_id"] = SystVariation(ZZ3P1FHist_idUp_,ZZ3P1FHist_idDn_);
        TH1D* ZZ4P0FHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_id"] = SystVariation(ZZ4P0FHist_idUp_,ZZ4P0FHist_idDn_);
        TH1D* ZZ3P1FHist_ffUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffUp").c_str() )->Clone();
        TH1D* ZZ3P1FHist_ffDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffDn").c_str() )->Clone();
        syst_["ZZ3P1FHist_ff"] = SystVariation(ZZ3P1FHist_ffUp_,ZZ3P1FHist_ffDn_);
        TH1D* FF3P1FHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffUp").c_str() )->Clone();
        TH1D* FF3P1FHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffDn").c_str() )->Clone();
        syst_["FF3P1FHist_ff"] = SystVariation(FF3P1FHist_ffUp_,FF3P1FHist_ffDn_);
        TH1D* FF2P2FHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffUp").c_str() )->Clone();
        TH1D* FF2P2FHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffDn").c_str() )->Clone();
        syst_["FF2P2FHist_ff"] = SystVariation(FF2P2FHist_ffUp_,FF2P2FHist_ffDn_);
        TH1D* ZZ4P0FHist_scaleUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_scaleUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_scaleDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_scaleDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_scale"] = SystVariation(ZZ4P0FHist_scaleUp_,ZZ4P0FHist_scaleDn_);
        TH1D* ZZ4P0FHist_sigmaUp_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_sigmaUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_sigmaDn_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_sigmaDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_sigma"] = SystVariation(ZZ4P0FHist_sigmaUp_,ZZ4P0FHist_sigmaDn_);

        for (const auto& element : syst_) {
          if (TString(element.first).Contains("ZZ")) {
            element.second.up_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
            element.second.dn_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
          }
        }

        for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
          std::map<std::string,SystVariation> init;

          TH1D* sigHist_idUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
          TH1D* sigHist_idDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
          init["sigHist_id"] = SystVariation(sigHist_idUp,sigHist_idDn);
          TH1D* sigHist_scaleUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_scaleUp").c_str() )->Clone();
          TH1D* sigHist_scaleDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_scaleDn").c_str() )->Clone();
          init["sigHist_scale"] = SystVariation(sigHist_scaleUp,sigHist_scaleDn);
          TH1D* sigHist_sigmaUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_sigmaUp").c_str() )->Clone();
          TH1D* sigHist_sigmaDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedEleCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_sigmaDn").c_str() )->Clone();
          init["sigHist_sigma"] = SystVariation(sigHist_sigmaUp,sigHist_sigmaDn);

          for (const auto& element : init) {
            element.second.up_->Scale( lumi*1000.*sigLumi / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );
            element.second.dn_->Scale( lumi*1000.*sigLumi / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );
          }

          sigSyst_.push_back(init);
        }
      }
    };

    void compare(TPad* pad, int rebin=1) {
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

      TH1D* FF3P1Fsubtracted = subtract(FF3P1FHist_,ZZ3P1FHist_);
      TH1D* FFsubtractedFinal = subtract(FF3P1Fsubtracted,FF2P2FHist_);
      TH1D* FFsubtractedFinal_copy = (TH1D*)FFsubtractedFinal->Clone();

      for (unsigned ibin = 0; ibin < FFsubtractedFinal->GetNbinsX()+2; ibin++) {
        double val = FFsubtractedFinal->GetBinContent(ibin);
        double err = FFsubtractedFinal->GetBinError(ibin);

        if ( val < 0.) {
          ZZ4P0FHist_->SetBinContent( ibin, ZZ4P0FHist_->GetBinContent(ibin) - val );
          ZZ4P0FHist_->SetBinError( ibin, std::hypot( err, ZZ4P0FHist_->GetBinError(ibin) ) );
          FFsubtractedFinal->SetBinContent(ibin,0.);
          FFsubtractedFinal->SetBinError(ibin,0.);
        }
      }

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(FFsubtractedFinal);
      stack->Add(ZZ4P0FHist_);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));

      pad->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        dataHist_->GetXaxis()->SetRangeUser(0.,1000.);
      }

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(FFsubtractedFinal,"Nonprompt");
      legend->AddEntry(ZZ4P0FHist_,"ZZ");

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);

          for (const auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_->Rebin(rebin);
            sigSyst_.at(idx).at(element.first).dn_->Rebin(rebin);
          }

          sigHist_.at(idx)->Draw("hist&same");
        }

        //legend->AddEntry(sigHist_.at(3),"H250/750/2000A1");
        legend->AddEntry(sigHist_.at(0),"X750Y100");
      }

      legend->Draw();

      if (!syst_.empty()) {
        TH1D* FF3P1Fsubtracted_idUp = subtract(subtract(FF3P1FHist_,syst_.at("ZZ3P1FHist_id").up_),FF2P2FHist_);
        TH1D* FF3P1Fsubtracted_idDn = subtract(subtract(FF3P1FHist_,syst_.at("ZZ3P1FHist_id").dn_),FF2P2FHist_);
        TH1D* FF3P1Fsubtracted_ffUp = subtract(subtract(syst_.at("FF3P1FHist_ff").up_,syst_.at("ZZ3P1FHist_ff").up_),syst_.at("FF2P2FHist_ff").up_);
        TH1D* FF3P1Fsubtracted_ffDn = subtract(subtract(syst_.at("FF3P1FHist_ff").dn_,syst_.at("ZZ3P1FHist_ff").dn_),syst_.at("FF2P2FHist_ff").dn_);

        TH1D* idUp = (TH1D*)syst_.at("ZZ4P0FHist_id").up_->Clone();
        TH1D* idDn = (TH1D*)syst_.at("ZZ4P0FHist_id").dn_->Clone();
        idUp->Add(FF3P1Fsubtracted_idUp);
        idDn->Add(FF3P1Fsubtracted_idDn);

        TH1D* ffUpAdded = (TH1D*)FF3P1Fsubtracted_ffUp->Clone();
        TH1D* ffDnAdded = (TH1D*)FF3P1Fsubtracted_ffDn->Clone();
        ffUpAdded->Add(ZZ4P0FHist_);
        ffDnAdded->Add(ZZ4P0FHist_);     

        TH1D* scaleUpAdded = (TH1D*)FFsubtractedFinal_copy->Clone();
        TH1D* scaleDnAdded = (TH1D*)FFsubtractedFinal_copy->Clone();
        scaleUpAdded->Add(syst_.at("ZZ4P0FHist_scale").up_);
        scaleDnAdded->Add(syst_.at("ZZ4P0FHist_scale").dn_);
        TH1D* sigmaUpAdded = (TH1D*)FFsubtractedFinal_copy->Clone();
        TH1D* sigmaDnAdded = (TH1D*)FFsubtractedFinal_copy->Clone();
        sigmaUpAdded->Add(syst_.at("ZZ4P0FHist_sigma").up_);
        sigmaDnAdded->Add(syst_.at("ZZ4P0FHist_sigma").dn_);

        std::vector<double> x0, y0, errx, erryDn, erryUp;

        for (unsigned idx = 1; idx <= tmpHist->GetNbinsX(); idx++) {
          x0.push_back(tmpHist->GetBinCenter(idx));
          y0.push_back(tmpHist->GetBinContent(idx));
          errx.push_back(tmpHist->GetBinWidth(idx)/2.);

          double valIdUp = idUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valIdDn = tmpHist->GetBinContent(idx) - idDn->GetBinContent(idx);
          double valFFup = ffUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdn = tmpHist->GetBinContent(idx) - ffDnAdded->GetBinContent(idx);
          double valScaleUp = scaleUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valScaleDn = tmpHist->GetBinContent(idx) - scaleDnAdded->GetBinContent(idx);
          double valSigmaUp = sigmaUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valSigmaDn = tmpHist->GetBinContent(idx) - sigmaDnAdded->GetBinContent(idx);

          erryUp.push_back(std::sqrt(valIdUp*valIdUp + valFFup*valFFup + valScaleUp*valScaleUp + valSigmaUp*valSigmaUp));
          erryDn.push_back(std::sqrt(valIdDn*valIdDn + valFFdn*valFFdn + valScaleDn*valScaleDn + valSigmaDn*valSigmaDn));
        }

        auto gr = new TGraphAsymmErrors(tmpHist->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kGray+3);
        gr->SetLineColor(kGray+3);
        gr->SetFillStyle(3004);
        gr->Draw("2");

        if (dir_) {
          dir_->WriteTObject(dataHist_,"data_obs");
          dir_->WriteTObject(FF2P2FHist_,"2P2F");
          dir_->WriteTObject(FF3P1Fsubtracted,"3P1F");
          dir_->WriteTObject(ZZ4P0FHist_,"ZZ");

          dir_->WriteTObject(FF3P1Fsubtracted_ffUp,"3P1F_resolvedEleFakeFactorUp");
          dir_->WriteTObject(FF3P1Fsubtracted_ffDn,"3P1F_resolvedEleFakeFactorDown");
          dir_->WriteTObject(syst_.at("FF2P2FHist_ff").up_,"2P2F_resolvedEleFakeFactorUp");
          dir_->WriteTObject(syst_.at("FF2P2FHist_ff").dn_,"2P2F_resolvedEleFakeFactorDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_id").up_,"ZZ_modHeepIdUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_id").dn_,"ZZ_modHeepIdDown");
          dir_->WriteTObject(FF3P1Fsubtracted_idUp,"3P1F_modHeepIdUp");
          dir_->WriteTObject(FF3P1Fsubtracted_idDn,"3P1F_modHeepIdDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_scale").up_,"ZZ_elEnergyScaleUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_scale").dn_,"ZZ_elEnergyScaleDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_sigma").up_,"ZZ_elEnergySigmaUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_sigma").dn_,"ZZ_elEnergySigmaDown");

          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_id").up_,sigFiles_.at(idx).name_+"_modHeepIdUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_id").dn_,sigFiles_.at(idx).name_+"_modHeepIdDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_scale").up_,sigFiles_.at(idx).name_+"_elEnergyScaleUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_scale").dn_,sigFiles_.at(idx).name_+"_elEnergyScaleDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_sigma").up_,sigFiles_.at(idx).name_+"_elEnergySigmaUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_sigma").dn_,sigFiles_.at(idx).name_+"_elEnergySigmaDown");
          }
        }
      }
    };
  };

  auto aloader4P0F = HistLoader4P0F(datafile,WZfile,ZZfile,sigsamples);
  auto aloader4P0F1 = HistLoader4P0F(datafile1,WZfile1,ZZfile1,sigsamples1);
  auto aloader4P0F2 = HistLoader4P0F(datafile2,WZfile2,ZZfile2,sigsamples2);
  auto aloader4P0F3 = HistLoader4P0F(datafile3,WZfile3,ZZfile3,sigsamples3);
  aloader4P0F.load("llll_invM",postfix.Data());

  if (era.Contains("run2")) {
    aloader4P0F1.load("llll_invM","20UL16APV");
    aloader4P0F2.load("llll_invM","20UL17");
    aloader4P0F3.load("llll_invM","20UL18");
    aloader4P0F.add(aloader4P0F1);
    aloader4P0F.add(aloader4P0F2);
    aloader4P0F.add(aloader4P0F3);
  }

  aloader4P0F.preparecard("REFF_"+era+"_datacard.root","resolvedEle");
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REFF_4P0F_CR_llll_invM_zoomed.eps");
  aloader4P0F.close();

  p1->SetLogy(0);
  canvas_1->SetLogy(0);

  aloader4P0F.load("ll1ll2_dr",postfix.Data());

  if (era.Contains("run2")) {
    aloader4P0F1.load("ll1ll2_dr","20UL16APV");
    aloader4P0F2.load("ll1ll2_dr","20UL17");
    aloader4P0F3.load("ll1ll2_dr","20UL18");
    aloader4P0F.add(aloader4P0F1);
    aloader4P0F.add(aloader4P0F2);
    aloader4P0F.add(aloader4P0F3);
  }

  aloader4P0F.compare(canvas_1);
  SaveAs(canvas_1,"REFF_4P0F_CR_ll1ll2_dr.png");

  aloader4P0F.load("ll1_invM",postfix.Data());

  if (era.Contains("run2")) {
    aloader4P0F1.load("ll1_invM","20UL16APV");
    aloader4P0F2.load("ll1_invM","20UL17");
    aloader4P0F3.load("ll1_invM","20UL18");
    aloader4P0F.add(aloader4P0F1);
    aloader4P0F.add(aloader4P0F2);
    aloader4P0F.add(aloader4P0F3);
  }

  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REFF_4P0F_CR_ll1_invM.png");

  aloader4P0F.load("ll1_dr",postfix.Data());

  if (era.Contains("run2")) {
    aloader4P0F1.load("ll1_dr","20UL16APV");
    aloader4P0F2.load("ll1_dr","20UL17");
    aloader4P0F3.load("ll1_dr","20UL18");
    aloader4P0F.add(aloader4P0F1);
    aloader4P0F.add(aloader4P0F2);
    aloader4P0F.add(aloader4P0F3);
  }

  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REFF_4P0F_CR_ll1_dr.png");

  aloader4P0F.load("ll2_invM",postfix.Data());

  if (era.Contains("run2")) {
    aloader4P0F1.load("ll2_invM","20UL16APV");
    aloader4P0F2.load("ll2_invM","20UL17");
    aloader4P0F3.load("ll2_invM","20UL18");
    aloader4P0F.add(aloader4P0F1);
    aloader4P0F.add(aloader4P0F2);
    aloader4P0F.add(aloader4P0F3);
  }

  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REFF_4P0F_CR_ll2_invM.png");

  aloader4P0F.load("ll2_dr",postfix.Data());

  if (era.Contains("run2")) {
    aloader4P0F1.load("ll2_dr","20UL16APV");
    aloader4P0F2.load("ll2_dr","20UL17");
    aloader4P0F3.load("ll2_dr","20UL18");
    aloader4P0F.add(aloader4P0F1);
    aloader4P0F.add(aloader4P0F2);
    aloader4P0F.add(aloader4P0F3);
  }

  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REFF_4P0F_CR_ll2_dr.png");

  return;
}
