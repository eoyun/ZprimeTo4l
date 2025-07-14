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

void runRMFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Preliminary";  // default extra text is "Preliminary"

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

  static TString anlyzrMC = "resolvedMuCRanalyzer"+postfix;
  static TString anlyzrData = "resolvedMuCRanalyzerData";

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

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  auto sampleNames = std::vector<TString> {"H750A100"};

  /*auto sampleNames = std::vector<TString> {
    "H250A0p4",
    "H250A0p6",
    "H250A0p8",
    "H250A1",
    "H250A1p5",
    "H250A2",
    "H250A5",
    "H250A10",
    "H250A50",
    "H250A100",
    "H500A0p4",
    "H500A0p6",
    "H500A0p8",
    "H500A1",
    "H500A1p5",
    "H500A2",
    "H500A5",
    "H500A10",
    "H500A50",
    "H500A100",
    "H750A0p4",
    "H750A0p6",
    "H750A0p8",
    "H750A1",
    "H750A1p5",
    "H750A2",
    "H750A5",
    "H750A10",
    "H750A50",
    "H750A100",
    "H750A250",
    "H1000A0p4",
    "H1000A0p6",
    "H1000A0p8",
    "H1000A1",
    "H1000A1p5",
    "H1000A2",
    "H1000A5",
    "H1000A10",
    "H1000A50",
    "H1000A100",
    "H1000A250",
    "H1500A0p4",
    "H1500A0p6",
    "H1500A0p8",
    "H1500A1",
    "H1500A1p5",
    "H1500A2",
    "H1500A5",
    "H1500A10",
    "H1500A50",
    "H1500A100",
    "H1500A250",
    "H1500A500",
    "H2000A0p4",
    "H2000A0p6",
    "H2000A0p8",
    "H2000A1",
    "H2000A1p5",
    "H2000A2",
    "H2000A5",
    "H2000A10",
    "H2000A50",
    "H2000A100",
    "H2000A750",
    "H5000A1",
    "H5000A2",
    "H5000A5",
    "H5000A10",
    "H5000A50",
    "H5000A100",
    "H5000A750",
    "H5000A2000",
    "H250Z0p4",
    "H250Z0p6",
    "H250Z0p8",
    "H250Z1",
    "H250Z1p5",
    "H250Z2",
    "H250Z5",
    "H250Z10",
    "H250Z50",
    "H250Z100",
    "H500Z0p4",
    "H500Z0p6",
    "H500Z0p8",
    "H500Z1",
    "H500Z1p5",
    "H500Z2",
    "H500Z5",
    "H500Z10",
    "H500Z50",
    "H500Z100",
    "H500Z250",
    "H750Z0p4",
    "H750Z0p6",
    "H750Z0p8",
    "H750Z1",
    "H750Z1p5",
    "H750Z2",
    "H750Z5",
    "H750Z10",
    "H750Z50",
    "H750Z100",
    "H750Z250",
    "H750Z500",
    "H1000Z0p4",
    "H1000Z0p6",
    "H1000Z0p8",
    "H1000Z1",
    "H1000Z1p5",
    "H1000Z2",
    "H1000Z5",
    "H1000Z10",
    "H1000Z50",
    "H1000Z100",
    "H1000Z500",
    "H1000Z750",
    "H1500Z0p4",
    "H1500Z0p6",
    "H1500Z0p8",
    "H1500Z1",
    "H1500Z1p5",
    "H1500Z2",
    "H1500Z5",
    "H1500Z10",
    "H1500Z50",
    "H1500Z100",
    "H1500Z500",
    "H1500Z1000",
    "H2000Z0p4",
    "H2000Z0p6",
    "H2000Z0p8",
    "H2000Z1",
    "H2000Z1p5",
    "H2000Z2",
    "H2000Z5",
    "H2000Z10",
    "H2000Z50",
    "H2000Z100",
    "H2000Z500",
    "H2000Z1000",
    "H2000Z1500"
  };*/

  std::vector<SigSample> sigsamples, sigsamples1, sigsamples2, sigsamples3;

  for (const auto& name : sampleNames) {
    sigsamples.push_back(SigSample(new TFile("MuAnalyzer_"+fname+"_"+name+".root","READ"),name));
    sigsamples1.push_back(SigSample(new TFile("MuAnalyzer_20UL16APV_"+name+".root","READ"),name));
    sigsamples2.push_back(SigSample(new TFile("MuAnalyzer_20UL17_"+name+".root","READ"),name));
    sigsamples3.push_back(SigSample(new TFile("MuAnalyzer_20UL18_"+name+".root","READ"),name));
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
  p1->SetBottomMargin(0.0);
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
    }

    ~HistLoaderBase()=default;

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

      if (FFHist_ffUp_) {
        delete FFHist_ffUp_, FFHist_ffDn_;
        FFHist_ffUp_ = nullptr;
      }
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* WZHist_ = nullptr;
    TH1D* ZZHist_ = nullptr;
    TH1D* FFHist_ = nullptr;
    TH1D* FFHist_ffUp_ = nullptr;
    TH1D* FFHist_ffDn_ = nullptr;

  public:
    void add(const HistLoader3P1F& other) {
      this->dataHist_->Add(other.dataHist_);
      this->WZHist_->Add(other.WZHist_);
      this->ZZHist_->Add(other.ZZHist_);
      this->FFHist_->Add(other.FFHist_);

      if (FFHist_ffUp_) {
        this->FFHist_ffUp_->Add(other.FFHist_ffUp_);
        this->FFHist_ffDn_->Add(other.FFHist_ffDn_);
      }
    }

    void load(const std::string& nameNum, const std::string& name, std::string anlyzrEra="") {
      if (dataHist_)
        delete dataHist_, WZHist_, ZZHist_, FFHist_;

      if (FFHist_ffUp_) {
        delete FFHist_ffUp_, FFHist_ffDn_;
        FFHist_ffUp_ = nullptr;
      }

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum).c_str() )->Clone();
      WZHist_ = (TH1D*)WZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/")+nameNum).c_str() )->Clone();
      ZZHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/")+nameNum).c_str() )->Clone();
      FFHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name).c_str() )->Clone();

      const double lumi = retrieveLumi(anlyzrEra);
      WZHist_->Scale( WZxsec_*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZHist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      WZHist_->SetFillColor(TColor::GetColor("#7a21dd"));
      ZZHist_->SetFillColor(TColor::GetColor("#5790fc"));
      WZHist_->SetLineWidth(0);
      ZZHist_->SetLineWidth(0);
      FFHist_->SetFillColor(TColor::GetColor("#9c9ca1"));
      FFHist_->SetLineWidth(0);

      if (TString(nameNum).Contains("llll_invM")) {
        FFHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffUp").c_str() )->Clone();
        FFHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffDn").c_str() )->Clone();
      }
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
          double valFFup = FFHist_ffUp_->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdn = tmpHist->GetBinContent(idx) - FFHist_ffDn_->GetBinContent(idx);

          erryUp.push_back(valFFup);
          erryDn.push_back(valFFdn);

          r0.push_back(1.);

          double rFFup = valFFup/tmpHist->GetBinContent(idx);
          double rFFdn = valFFdn/tmpHist->GetBinContent(idx);

          errRup.push_back( tmpHist->GetBinContent(idx) > 0. ? rFFup : 0. );
          errRdn.push_back( tmpHist->GetBinContent(idx) > 0. ? rFFdn : 0. );
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
  /*auto aloader3P1F = HistLoader3P1F(datafile,WZfile,ZZfile);
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
  SaveAs(canvas_2,"RMFF_3P1F_CR_llll_invM.pdf",p1);*/

/*  aloader3P1F.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"RMFF_3P1F_CR_ll1ll2_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"RMFF_3P1F_CR_ll1_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"RMFF_3P1F_CR_ll1_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"RMFF_3P1F_CR_ll2_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF",postfix.Data());
  aloader3P1F1.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL16APV");
  aloader3P1F2.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL17");
  aloader3P1F3.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF","20UL18");
  aloader3P1F.add(aloader3P1F1);
  aloader3P1F.add(aloader3P1F2);
  aloader3P1F.add(aloader3P1F3);
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"RMFF_3P1F_CR_ll2_dr.png",p1);*/

  // 4P0F CR

  class HistLoader4P0F : public HistLoaderBase {
  public:
    HistLoader4P0F(TFile* adatafile, TFile* aWZfile, TFile* aZZfile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile,aWZfile,aZZfile),
      sigFiles_(sigFiles) {}

    ~HistLoader4P0F() {
      if (dataHist_)
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FF2P2FHist_, FFdr03Hist_, ZZdr03Hist_;
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* ZZ3P1FHist_ = nullptr;
    TH1D* ZZ4P0FHist_ = nullptr;
    TH1D* FF3P1FHist_ = nullptr;
    TH1D* FF2P2FHist_ = nullptr;
    TH1D* FFdr03Hist_ = nullptr;
    TH1D* ZZdr03Hist_ = nullptr;

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
      this->FFdr03Hist_->Add(other.FFdr03Hist_);
      this->ZZdr03Hist_->Add(other.ZZdr03Hist_);

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
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FF2P2FHist_, FFdr03Hist_, ZZdr03Hist_;
      }

      if (!syst_.empty())
        syst_.clear();

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/4P0F_CR_")+name).c_str() )->Clone();
      ZZ3P1FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      ZZ4P0FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name).c_str() )->Clone();
      FF3P1FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      FF2P2FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2").c_str() )->Clone();
      FFdr03Hist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CRdr03_")+name+"_xFF").c_str() )->Clone();
      ZZdr03Hist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/2P2F_CRdr03_")+name+"_xFF").c_str() )->Clone();

      const double lumi = retrieveLumi(anlyzrEra);
      sigHist_.clear();
      sigSyst_.clear();
      double sigLumi = 0.01; // 0.0001;

      if ( TString(name).Contains("llll_invM") ) {
        TH1D* ZZ3P1FHist_ffUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffUp").c_str() )->Clone();
        TH1D* ZZ3P1FHist_ffDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_xFF_ffDn").c_str() )->Clone();
        syst_["ZZ3P1FHist_ff"] = SystVariation(ZZ3P1FHist_ffUp,ZZ3P1FHist_ffDn);
        TH1D* FF3P1FHist_ffUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffUp").c_str() )->Clone();
        TH1D* FF3P1FHist_ffDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffDn").c_str() )->Clone();
        syst_["FF3P1FHist_ff"] = SystVariation(FF3P1FHist_ffUp,FF3P1FHist_ffDn);
        TH1D* FF2P2FHist_ffUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffUp").c_str() )->Clone();
        TH1D* FF2P2FHist_ffDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffDn").c_str() )->Clone();
        syst_["FF2P2FHist_ff"] = SystVariation(FF2P2FHist_ffUp,FF2P2FHist_ffDn);
        TH1D* ZZ3P1FHist_scaleUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/3P1F_CR_")+name+"_altMuScale_xFF").c_str() )->Clone();
        TH1D* ZZ3P1FHist_scaleDn = variateDn(ZZ3P1FHist_,ZZ3P1FHist_scaleUp);
        syst_["ZZ3P1FHist_scale"] = SystVariation(ZZ3P1FHist_scaleUp,ZZ3P1FHist_scaleDn);
        TH1D* FF3P1FHist_scaleUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_altMuScale_xFF").c_str() )->Clone();
        TH1D* FF3P1FHist_scaleDn = variateDn(FF3P1FHist_,FF3P1FHist_scaleUp);
        syst_["FF3P1FHist_scale"] = SystVariation(FF3P1FHist_scaleUp,FF3P1FHist_scaleDn);
        TH1D* FF2P2FHist_scaleUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_altMuScale_xFF2").c_str() )->Clone();
        TH1D* FF2P2FHist_scaleDn = variateDn(FF2P2FHist_,FF2P2FHist_scaleUp);
        syst_["FF2P2FHist_scale"] = SystVariation(FF2P2FHist_scaleUp,FF2P2FHist_scaleDn);
        TH1D* ZZ4P0FHist_scaleUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuScale").c_str() )->Clone();
        TH1D* ZZ4P0FHist_scaleDn = variateDn(ZZ4P0FHist_,ZZ4P0FHist_scaleUp);
        syst_["ZZ4P0FHist_scale"] = SystVariation(ZZ4P0FHist_scaleUp,ZZ4P0FHist_scaleDn);
        TH1D* ZZ4P0FHist_smearUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuSmear").c_str() )->Clone();
        TH1D* ZZ4P0FHist_smearDn = variateDn(ZZ4P0FHist_,ZZ4P0FHist_smearUp);
        syst_["ZZ4P0FHist_smear"] = SystVariation(ZZ4P0FHist_smearUp,ZZ4P0FHist_smearDn);
        TH1D* ZZ4P0FHist_idUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_idDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_id"] = SystVariation(ZZ4P0FHist_idUp,ZZ4P0FHist_idDn);
        TH1D* ZZ4P0FHist_isoUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_isoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_isoDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_isoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_iso"] = SystVariation(ZZ4P0FHist_isoUp,ZZ4P0FHist_isoDn);
        TH1D* ZZ4P0FHist_boostIsoUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_boostIsoDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_boostIso"] = SystVariation(ZZ4P0FHist_boostIsoUp,ZZ4P0FHist_boostIsoDn);
        TH1D* ZZ4P0FHist_trigUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_trigDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_trig"] = SystVariation(ZZ4P0FHist_trigUp,ZZ4P0FHist_trigDn);
        TH1D* ZZ4P0FHist_recoUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_recoUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_recoDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_recoDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_reco"] = SystVariation(ZZ4P0FHist_recoUp,ZZ4P0FHist_recoDn);
        TH1D* ZZ4P0FHist_PUrwgtUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_PUrwgtDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_PUrwgt"] = SystVariation(ZZ4P0FHist_PUrwgtUp,ZZ4P0FHist_PUrwgtDn);
        TH1D* ZZ4P0FHist_prefireUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireUp").c_str() )->Clone();
        TH1D* ZZ4P0FHist_prefireDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireDn").c_str() )->Clone();
        syst_["ZZ4P0FHist_prefire"] = SystVariation(ZZ4P0FHist_prefireUp,ZZ4P0FHist_prefireDn);

        TH1D* FFdr03Hist_ffUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CRdr03_")+name+"_xFF_ffUp").c_str() )->Clone();
        TH1D* FFdr03Hist_ffDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CRdr03_")+name+"_xFF_ffDn").c_str() )->Clone();
        syst_["FFdr03Hist_ff"] = SystVariation(FFdr03Hist_ffUp,FFdr03Hist_ffDn);
        TH1D* ZZdr03Hist_ffUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/2P2F_CRdr03_")+name+"_xFF_ffUp").c_str() )->Clone();
        TH1D* ZZdr03Hist_ffDn = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/2P2F_CRdr03_")+name+"_xFF_ffDn").c_str() )->Clone();
        syst_["ZZdr03Hist_ff"] = SystVariation(ZZdr03Hist_ffUp,ZZdr03Hist_ffDn);

        TH1D* FFdr03Hist_muScaleUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CRdr03_")+name+"_altMuScale_xFF").c_str() )->Clone();
        TH1D* FFdr03Hist_muScaleDn = variateDn(FFdr03Hist_,FFdr03Hist_muScaleUp);
        syst_["FFdr03Hist_muScale"] = SystVariation(FFdr03Hist_muScaleUp,FFdr03Hist_muScaleDn);
        TH1D* ZZdr03Hist_muScaleUp = (TH1D*)ZZfile_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/2P2F_CRdr03_")+name+"_altMuScale_xFF").c_str() )->Clone();
        TH1D* ZZdr03Hist_muScaleDn = variateDn(ZZdr03Hist_,ZZdr03Hist_muScaleUp);
        syst_["ZZdr03Hist_muScale"] = SystVariation(ZZdr03Hist_muScaleUp,ZZdr03Hist_muScaleDn);

        for (const auto& element : syst_) {
          if (TString(element.first).Contains("ZZ")) {
            element.second.up_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
            element.second.dn_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
          }
        }

        for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
          sigHist_.push_back( (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name).c_str() )->Clone() );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          std::map<std::string,SystVariation> init;

          TH1D* sigHist_scaleUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuScale").c_str() )->Clone();
          TH1D* sigHist_scaleDn = variateDn(sigHist_.back(),sigHist_scaleUp);
          init["sigHist_scale"] = SystVariation(sigHist_scaleUp,sigHist_scaleDn);
          TH1D* sigHist_smearUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_altMuSmear").c_str() )->Clone();
          TH1D* sigHist_smearDn = variateDn(sigHist_.back(),sigHist_smearUp);
          init["sigHist_smear"] = SystVariation(sigHist_smearUp,sigHist_smearDn);
          TH1D* sigHist_idUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
          TH1D* sigHist_idDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
          init["sigHist_id"] = SystVariation(sigHist_idUp,sigHist_idDn);
          TH1D* sigHist_isoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_isoUp").c_str() )->Clone();
          TH1D* sigHist_isoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_isoDn").c_str() )->Clone();
          init["sigHist_iso"] = SystVariation(sigHist_isoUp,sigHist_isoDn);
          TH1D* sigHist_boostIsoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoUp").c_str() )->Clone();
          TH1D* sigHist_boostIsoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_muBoostIsoDn").c_str() )->Clone();
          init["sigHist_boostIso"] = SystVariation(sigHist_boostIsoUp,sigHist_boostIsoDn);
          TH1D* sigHist_trigUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigUp").c_str() )->Clone();
          TH1D* sigHist_trigDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_trigDn").c_str() )->Clone();
          init["sigHist_trig"] = SystVariation(sigHist_trigUp,sigHist_trigDn);
          TH1D* sigHist_recoUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_recoUp").c_str() )->Clone();
          TH1D* sigHist_recoDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_recoDn").c_str() )->Clone();
          init["sigHist_reco"] = SystVariation(sigHist_recoUp,sigHist_recoDn);
          TH1D* sigHist_PUrwgtUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtUp").c_str() )->Clone();
          TH1D* sigHist_PUrwgtDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_PUrwgtDn").c_str() )->Clone();
          init["sigHist_PUrwgt"] = SystVariation(sigHist_PUrwgtUp,sigHist_PUrwgtDn);
          TH1D* sigHist_prefireUp = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireUp").c_str() )->Clone();
          TH1D* sigHist_prefireDn = (TH1D*)sigFiles_.at(idx).file_->Get( (std::string("resolvedMuCRanalyzer"+anlyzrEra+"/4P0F_CR_")+name+"_prefireDn").c_str() )->Clone();
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
      ZZdr03Hist_->Scale( ZZxsec_*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZ4P0FHist_->SetFillColor(TColor::GetColor("#5790fc"));
      ZZ4P0FHist_->SetLineWidth(0);
      FF3P1FHist_->SetFillColor(TColor::GetColor("#9c9ca1"));
      FF3P1FHist_->SetLineWidth(0);
      FF2P2FHist_->SetFillColor(TColor::GetColor("#9c9ca1"));
      FF2P2FHist_->SetLineWidth(0);
      FFdr03Hist_->SetFillColor(TColor::GetColor("#92dadd"));
      FFdr03Hist_->SetLineWidth(0);
    };

    void compare(TPad* padUp, int rebin=1, TPad* padDn=nullptr) {
      if ( FF2P2FHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FF2P2FHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FF2P2FHist_->Rebin( FF2P2FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FF3P1FHist_->Rebin( FF3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        ZZ3P1FHist_->Rebin( ZZ3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FFdr03Hist_->Rebin( FFdr03Hist_->GetNbinsX()/dataHist_->GetNbinsX() );
        ZZdr03Hist_->Rebin( ZZdr03Hist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( element.second.up_->GetNbinsX()/dataHist_->GetNbinsX() );
            element.second.dn_->Rebin( element.second.dn_->GetNbinsX()/dataHist_->GetNbinsX() );
          }
        }
      }

      if (rebin!=1) {
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
        FFdr03Hist_->Rebin(rebin);
        ZZdr03Hist_->Rebin(rebin);
        dataHist_->Rebin(rebin);
      }

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
      TH1D* FFsubtractedFinal = subtract(subtract(FF3P1Fsubtracted,FF2P2FHist_),FF2P2FHist_);
      truncateNegativeBin(FFsubtractedFinal);
      FFsubtractedFinal->Add(FF2P2FHist_);

      TH1D* FFdr03Final = subtract(FFdr03Hist_,ZZdr03Hist_);
      truncateNegativeBin(FFdr03Final);

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(FFsubtractedFinal);
      stack->Add(FFdr03Final);
      stack->Add(ZZ4P0FHist_);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));
      tmpHist->Add((TH1*)stack->GetHists()->At(2));

      padUp->cd();
      dataHist_->SetMaximum(200.*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));

      if ( TString(dataHist_->GetName()).Contains("llll_invM") && rebin!=1 ) {
        dataHist_->GetXaxis()->SetRangeUser(0.,1500.);
        dataHist_->GetXaxis()->SetTitle("M(4#mu) [GeV]");
      }

      padUp->SetLogy();
      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMinimum(0.2);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.5,0.52,0.95,0.93);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(FFsubtractedFinal,"Nonprompt (#DeltaR > 0.3)");
      legend->AddEntry(FFdr03Final,"Nonprompt (#DeltaR < 0.3)");
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

        // legend->AddEntry(sigHist_.at(3),"H250/750/2000A1");
        legend->AddEntry(sigHist_.at(0),"M_{X} = 750 GeV, M_{Y} = 100 GeV");
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
        TH1D* FF3P1Fsubtracted_ffUp = subtract(subtract(subtract(syst_.at("FF3P1FHist_ff").up_,syst_.at("ZZ3P1FHist_ff").up_),syst_.at("FF2P2FHist_ff").up_),syst_.at("FF2P2FHist_ff").up_);
        TH1D* FF3P1Fsubtracted_ffDn = subtract(subtract(subtract(syst_.at("FF3P1FHist_ff").dn_,syst_.at("ZZ3P1FHist_ff").dn_),syst_.at("FF2P2FHist_ff").dn_),syst_.at("FF2P2FHist_ff").dn_);
        truncateNegativeBin(FF3P1Fsubtracted_ffUp);
        truncateNegativeBin(FF3P1Fsubtracted_ffDn);
        FF3P1Fsubtracted_ffUp->Add(syst_.at("FF2P2FHist_ff").up_);
        FF3P1Fsubtracted_ffDn->Add(syst_.at("FF2P2FHist_ff").dn_);
        TH1D* ffUpAdded = (TH1D*)FF3P1Fsubtracted_ffUp->Clone();
        TH1D* ffDnAdded = (TH1D*)FF3P1Fsubtracted_ffDn->Clone();
        ffUpAdded->Add(ZZ4P0FHist_);
        ffDnAdded->Add(ZZ4P0FHist_);
        ffUpAdded->Add(FFdr03Final);
        ffDnAdded->Add(FFdr03Final);

        TH1D* FFdr03subtracted_ffUp = subtract(syst_.at("FFdr03Hist_ff").up_,syst_.at("ZZdr03Hist_ff").up_);
        TH1D* FFdr03subtracted_ffDn = subtract(syst_.at("FFdr03Hist_ff").dn_,syst_.at("ZZdr03Hist_ff").dn_);
        truncateNegativeBin(FFdr03subtracted_ffUp);
        truncateNegativeBin(FFdr03subtracted_ffDn);
        TH1D* ffdr03upAdded = (TH1D*)FFdr03subtracted_ffUp->Clone();
        TH1D* ffdr03dnAdded = (TH1D*)FFdr03subtracted_ffDn->Clone();
        ffdr03upAdded->Add(ZZ4P0FHist_);
        ffdr03dnAdded->Add(ZZ4P0FHist_);
        ffdr03upAdded->Add(FFsubtractedFinal);
        ffdr03dnAdded->Add(FFsubtractedFinal);

        TH1D* FF3P1Fsubtracted_scaleUp = subtract(subtract(subtract(syst_.at("FF3P1FHist_scale").up_,syst_.at("ZZ3P1FHist_scale").up_),syst_.at("FF2P2FHist_scale").up_),syst_.at("FF2P2FHist_scale").up_);
        TH1D* FF3P1Fsubtracted_scaleDn = subtract(subtract(subtract(syst_.at("FF3P1FHist_scale").dn_,syst_.at("ZZ3P1FHist_scale").dn_),syst_.at("FF2P2FHist_scale").dn_),syst_.at("FF2P2FHist_scale").dn_);
        truncateNegativeBin(FF3P1Fsubtracted_scaleUp);
        truncateNegativeBin(FF3P1Fsubtracted_scaleDn);
        TH1D* FFdr03Final_muScaleUp = subtract(syst_.at("FFdr03Hist_muScale").up_,syst_.at("ZZdr03Hist_muScale").up_);
        TH1D* FFdr03Final_muScaleDn = subtract(syst_.at("FFdr03Hist_muScale").dn_,syst_.at("ZZdr03Hist_muScale").dn_);
        truncateNegativeBin(FFdr03Final_muScaleUp);
        truncateNegativeBin(FFdr03Final_muScaleDn); 
        TH1D* scaleUpAdded = (TH1D*)FF3P1Fsubtracted_scaleUp->Clone();
        TH1D* scaleDnAdded = (TH1D*)FF3P1Fsubtracted_scaleDn->Clone();
        scaleUpAdded->Add(syst_.at("FF2P2FHist_scale").up_);
        scaleDnAdded->Add(syst_.at("FF2P2FHist_scale").dn_);
        scaleUpAdded->Add(syst_.at("ZZ4P0FHist_scale").up_);
        scaleDnAdded->Add(syst_.at("ZZ4P0FHist_scale").dn_);
        scaleUpAdded->Add(FFdr03Final_muScaleUp);
        scaleDnAdded->Add(FFdr03Final_muScaleDn);

        TH1D* smearUp = (TH1D*)syst_.at("ZZ4P0FHist_smear").up_->Clone();
        TH1D* smearDn = (TH1D*)syst_.at("ZZ4P0FHist_smear").dn_->Clone();
        smearUp->Add(FFsubtractedFinal);
        smearDn->Add(FFsubtractedFinal);
        smearUp->Add(FFdr03Final);
        smearDn->Add(FFdr03Final);

        TH1D* idUp = (TH1D*)syst_.at("ZZ4P0FHist_id").up_->Clone();
        TH1D* idDn = (TH1D*)syst_.at("ZZ4P0FHist_id").dn_->Clone();
        idUp->Add(FFsubtractedFinal);
        idDn->Add(FFsubtractedFinal);
        idUp->Add(FFdr03Final);
        idDn->Add(FFdr03Final);

        TH1D* isoUp = (TH1D*)syst_.at("ZZ4P0FHist_iso").up_->Clone();
        TH1D* isoDn = (TH1D*)syst_.at("ZZ4P0FHist_iso").dn_->Clone();
        isoUp->Add(FFsubtractedFinal);
        isoDn->Add(FFsubtractedFinal);
        isoUp->Add(FFdr03Final);
        isoDn->Add(FFdr03Final);

        TH1D* boostIsoUp = (TH1D*)syst_.at("ZZ4P0FHist_boostIso").up_->Clone();
        TH1D* boostIsoDn = (TH1D*)syst_.at("ZZ4P0FHist_boostIso").dn_->Clone();
        boostIsoUp->Add(FFsubtractedFinal);
        boostIsoDn->Add(FFsubtractedFinal);
        boostIsoUp->Add(FFdr03Final);
        boostIsoDn->Add(FFdr03Final);

        TH1D* trigUp = (TH1D*)syst_.at("ZZ4P0FHist_trig").up_->Clone();
        TH1D* trigDn = (TH1D*)syst_.at("ZZ4P0FHist_trig").dn_->Clone();
        trigUp->Add(FFsubtractedFinal);
        trigDn->Add(FFsubtractedFinal);
        trigUp->Add(FFdr03Final);
        trigDn->Add(FFdr03Final);

        TH1D* recoUp = (TH1D*)syst_.at("ZZ4P0FHist_reco").up_->Clone();
        TH1D* recoDn = (TH1D*)syst_.at("ZZ4P0FHist_reco").dn_->Clone();
        recoUp->Add(FFsubtractedFinal);
        recoDn->Add(FFsubtractedFinal);
        recoUp->Add(FFdr03Final);
        recoDn->Add(FFdr03Final);

        TH1D* PUrwgtUp = (TH1D*)syst_.at("ZZ4P0FHist_PUrwgt").up_->Clone();
        TH1D* PUrwgtDn = (TH1D*)syst_.at("ZZ4P0FHist_PUrwgt").dn_->Clone();
        PUrwgtUp->Add(FFsubtractedFinal);
        PUrwgtDn->Add(FFsubtractedFinal);
        PUrwgtUp->Add(FFdr03Final);
        PUrwgtDn->Add(FFdr03Final);

        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRup, errRdn;

        for (unsigned idx = 1; idx <= tmpHist->GetNbinsX(); idx++) {
          x0.push_back(tmpHist->GetBinCenter(idx));
          y0.push_back(tmpHist->GetBinContent(idx));
          errx.push_back(tmpHist->GetBinWidth(idx)/2.);

          double valFFup = ffUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdn = tmpHist->GetBinContent(idx) - ffDnAdded->GetBinContent(idx);
          double scaleUp = scaleUpAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double scaleDn = tmpHist->GetBinContent(idx) - scaleDnAdded->GetBinContent(idx);
          double valSmearUp = smearUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valSmearDn = tmpHist->GetBinContent(idx) - smearDn->GetBinContent(idx);
          double valIdUp = idUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valIdDn = tmpHist->GetBinContent(idx) - idDn->GetBinContent(idx);
          double valIsoUp = isoUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valIsoDn = tmpHist->GetBinContent(idx) - isoDn->GetBinContent(idx);
          double valboostIsoUp = boostIsoUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valboostIsoDn = tmpHist->GetBinContent(idx) - boostIsoDn->GetBinContent(idx);
          double valTrigUp = trigUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valTrigDn = tmpHist->GetBinContent(idx) - trigDn->GetBinContent(idx);
          double valRecoUp = recoUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valRecoDn = tmpHist->GetBinContent(idx) - recoDn->GetBinContent(idx);
          double valPUrwgtUp = PUrwgtUp->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valPUrwgtDn = tmpHist->GetBinContent(idx) - PUrwgtDn->GetBinContent(idx);
          double valFFdr03up = ffdr03upAdded->GetBinContent(idx) - tmpHist->GetBinContent(idx);
          double valFFdr03dn = tmpHist->GetBinContent(idx) - ffdr03dnAdded->GetBinContent(idx);
          double ZZnorm = 0.1*ZZ4P0FHist_->GetBinContent(idx);

          erryUp.push_back( std::sqrt(valFFup*valFFup + scaleUp*scaleUp + valSmearUp*valSmearUp
                                      + valIdUp*valIdUp + valIsoUp*valIsoUp + valTrigUp*valTrigUp + valRecoUp*valRecoUp + valboostIsoUp*valboostIsoUp + valPUrwgtUp*valPUrwgtUp + valFFdr03up*valFFdr03up + ZZnorm*ZZnorm) );
          erryDn.push_back( std::sqrt(valFFdn*valFFdn + scaleDn*scaleDn + valSmearDn*valSmearDn
                                      + valIdDn*valIdDn + valIsoDn*valIsoDn + valTrigDn*valTrigDn + valRecoDn*valRecoDn + valboostIsoDn*valboostIsoDn + valPUrwgtDn*valPUrwgtDn + valFFdr03dn*valFFdr03dn + ZZnorm*ZZnorm) );

          r0.push_back(1.);

          double rIdUp = valIdUp/tmpHist->GetBinContent(idx);
          double rIdDn = valIdDn/tmpHist->GetBinContent(idx);
          double rFFup = valFFup/tmpHist->GetBinContent(idx);
          double rFFdn = valFFdn/tmpHist->GetBinContent(idx);
          double rScaleUp = scaleUp/tmpHist->GetBinContent(idx);
          double rScaleDn = scaleDn/tmpHist->GetBinContent(idx);
          double rSmearUp = valSmearUp/tmpHist->GetBinContent(idx);
          double rSmearDn = valSmearDn/tmpHist->GetBinContent(idx);
          double rIsoUp = valIsoUp/tmpHist->GetBinContent(idx);
          double rIsoDn = valIsoDn/tmpHist->GetBinContent(idx);
          double rBoostIsoUp = valboostIsoUp/tmpHist->GetBinContent(idx);
          double rBoostIsoDn = valboostIsoDn/tmpHist->GetBinContent(idx);
          double rTrigUp = valTrigUp/tmpHist->GetBinContent(idx);
          double rTrigDn = valTrigDn/tmpHist->GetBinContent(idx);
          double rRecoUp = valRecoUp/tmpHist->GetBinContent(idx);
          double rRecoDn = valRecoDn/tmpHist->GetBinContent(idx);
          double rPUrwgtUp = valPUrwgtUp/tmpHist->GetBinContent(idx);
          double rPUrwgtDn = valPUrwgtDn/tmpHist->GetBinContent(idx);
          double rFFdr03up = valFFdr03up/tmpHist->GetBinContent(idx);
          double rFFdr03dn = valFFdr03dn/tmpHist->GetBinContent(idx);
          double rZZnorm = ZZnorm/tmpHist->GetBinContent(idx);

          double rUp = std::sqrt(rIdUp*rIdUp + rFFup*rFFup + rScaleUp*rScaleUp + rSmearUp*rSmearUp + rIsoUp*rIsoUp
                                      + rBoostIsoUp*rBoostIsoUp + rTrigUp*rTrigUp + rRecoUp*rRecoUp + rPUrwgtUp*rPUrwgtUp + rFFdr03up*rFFdr03up + rZZnorm*rZZnorm);
          double rDn = std::sqrt(rIdDn*rIdDn + rFFdn*rFFdn + rScaleDn*rScaleDn + rSmearDn*rSmearDn + rIsoDn*rIsoDn
                                      + rBoostIsoDn*rBoostIsoDn + rTrigDn*rTrigDn + rRecoDn*rRecoDn + rPUrwgtDn*rPUrwgtDn + rFFdr03dn*rFFdr03dn + rZZnorm*rZZnorm);

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
          dir_->WriteTObject(FFdr03Final,"NonpromptDR03");
          dir_->WriteTObject(ZZ4P0FHist_,"ZZ");

          dir_->WriteTObject(FF3P1Fsubtracted_ffUp,"Nonprompt_resolvedMuFakeFactorUp");
          dir_->WriteTObject(FF3P1Fsubtracted_ffDn,"Nonprompt_resolvedMuFakeFactorDown");
          dir_->WriteTObject(FFdr03subtracted_ffUp,"NonpromptDR03_resolvedMuDR03FakeFactorUp");
          dir_->WriteTObject(FFdr03subtracted_ffDn,"NonpromptDR03_resolvedMuDR03FakeFactorDown");
          dir_->WriteTObject(FF3P1Fsubtracted_scaleUp,"Nonprompt_muMomentumScaleUp");
          dir_->WriteTObject(FF3P1Fsubtracted_scaleDn,"Nonprompt_muMomentumScaleDown");
          dir_->WriteTObject(FFdr03Final_muScaleUp,"NonpromptDR03_muMomentumScaleUp");
          dir_->WriteTObject(FFdr03Final_muScaleDn,"NonpromptDR03_muMomentumScaleDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_scale").up_,"ZZ_muMomentumScaleUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_scale").dn_,"ZZ_muMomentumScaleDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_smear").up_,"ZZ_muMomentumSmearUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_smear").dn_,"ZZ_muMomentumSmearDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_id").up_,"ZZ_highPtIdUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_id").dn_,"ZZ_highPtIdDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_iso").up_,"ZZ_muLooseIsoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_iso").dn_,"ZZ_muLooseIsoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_boostIso").up_,"ZZ_muBoostIsoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_boostIso").dn_,"ZZ_muBoostIsoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_trig").up_,"ZZ_muTrigUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_trig").dn_,"ZZ_muTrigDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_reco").up_,"ZZ_muRecoUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_reco").dn_,"ZZ_muRecoDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_PUrwgt").up_,"ZZ_PUrwgtUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_PUrwgt").dn_,"ZZ_PUrwgtDown");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_prefire").up_,"ZZ_prefireUp");
          dir_->WriteTObject(syst_.at("ZZ4P0FHist_prefire").dn_,"ZZ_prefireDown");

          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_scale").up_,sigFiles_.at(idx).name_+"_muMomentumScaleUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_scale").dn_,sigFiles_.at(idx).name_+"_muMomentumScaleDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_smear").up_,sigFiles_.at(idx).name_+"_muMomentumSmearUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_smear").dn_,sigFiles_.at(idx).name_+"_muMomentumSmearDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_id").up_,sigFiles_.at(idx).name_+"_highPtIdUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_id").dn_,sigFiles_.at(idx).name_+"_highPtIdDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_iso").up_,sigFiles_.at(idx).name_+"_muLooseIsoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_iso").dn_,sigFiles_.at(idx).name_+"_muLooseIsoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_boostIso").up_,sigFiles_.at(idx).name_+"_muBoostIsoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_boostIso").dn_,sigFiles_.at(idx).name_+"_muBoostIsoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_trig").up_,sigFiles_.at(idx).name_+"_muTrigUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_trig").dn_,sigFiles_.at(idx).name_+"_muTrigDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_reco").up_,sigFiles_.at(idx).name_+"_muRecoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigHist_reco").dn_,sigFiles_.at(idx).name_+"_muRecoDown");
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

  //aloader4P0F.preparecard("RMFF_"+era+"_datacard.root","resolvedMu");
  aloader4P0F.compare(p1,4,p2); // 4
  SaveAs(canvas_2,"RMFF_4P0F_CR_llll_invM_zoomed_log.png",p1);
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
  SaveAs(canvas_1,"RMFF_4P0F_CR_ll1ll2_dr.png");

  aloader4P0F.load("ll1_invM",postfix.Data());
  aloader4P0F1.load("ll1_invM","20UL16APV");
  aloader4P0F2.load("ll1_invM","20UL17");
  aloader4P0F3.load("ll1_invM","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"RMFF_4P0F_CR_ll1_invM.png");

  aloader4P0F.load("ll1_dr",postfix.Data());
  aloader4P0F1.load("ll1_dr","20UL16APV");
  aloader4P0F2.load("ll1_dr","20UL17");
  aloader4P0F3.load("ll1_dr","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"RMFF_4P0F_CR_ll1_dr.png");

  aloader4P0F.load("ll2_invM",postfix.Data());
  aloader4P0F1.load("ll2_invM","20UL16APV");
  aloader4P0F2.load("ll2_invM","20UL17");
  aloader4P0F3.load("ll2_invM","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"RMFF_4P0F_CR_ll2_invM.png");

  aloader4P0F.load("ll2_dr",postfix.Data());
  aloader4P0F1.load("ll2_dr","20UL16APV");
  aloader4P0F2.load("ll2_dr","20UL17");
  aloader4P0F3.load("ll2_dr","20UL18");
  aloader4P0F.add(aloader4P0F1);
  aloader4P0F.add(aloader4P0F2);
  aloader4P0F.add(aloader4P0F3);
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"RMFF_4P0F_CR_ll2_dr.png");*/

  return;
}
