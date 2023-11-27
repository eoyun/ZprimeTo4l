#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runREMuFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static double valLumi = 0.;
  static TString postfix = era;

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
    valLumi = 19.5;
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
    valLumi = 16.8;
    postfix = "";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
    valLumi = 41.48;
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
    valLumi = 59.83;
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

  TFile* datafile = new TFile("REMuCR_"+era+"_data.root","READ");
  TFile* WZfile = new TFile("REMuCR_"+era+"_WZ.root","READ");
  TFile* ZZfile = new TFile("REMuCR_"+era+"_ZZ.root","READ");

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

  auto H750A1sample = SigSample(new TFile("REMuCR_"+era+"_H750A1.root","READ"),"H750A1");

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
    void load(const std::string& nameNum, const std::string& name) {
      if (dataHist_)
        delete dataHist_, WZHist_, ZZHist_, FFHist_;

      if (FFHist_ffUp_) {
        delete WZHist_idUp_, ZZHist_idUp_, WZHist_idDn_, ZZHist_idDn_, FFHist_ffUp_, FFHist_ffDn_, FFHist_ffUpM_, FFHist_ffDnM_;

        FFHist_ffUp_ = nullptr;
      }

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum).c_str() )->Clone();
      WZHist_ = (TH1D*)WZfile_->Get( (std::string(anlyzrMC+"/")+nameNum).c_str() )->Clone();
      ZZHist_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/")+nameNum).c_str() )->Clone();
      FFHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name).c_str() )->Clone();

      const double lumi = valLumi;

      if (TString(nameNum).Contains("llll_invM")) {
        WZHist_idUp_ = (TH1D*)WZfile_->Get( (std::string(anlyzrMC+"/")+nameNum+"_idUp").c_str() )->Clone();
        ZZHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/")+nameNum+"_idUp").c_str() )->Clone();
        WZHist_idDn_ = (TH1D*)WZfile_->Get( (std::string(anlyzrMC+"/")+nameNum+"_idDn").c_str() )->Clone();
        ZZHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/")+nameNum+"_idDn").c_str() )->Clone();
        FFHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffUpE").c_str() )->Clone();
        FFHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffDnE").c_str() )->Clone();
        FFHist_ffUpM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffUpM").c_str() )->Clone();
        FFHist_ffDnM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_ffDnM").c_str() )->Clone();

        WZHist_idUp_->Scale( 0.65*62.78*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        WZHist_idDn_->Scale( 0.65*62.78*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZHist_idUp_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZHist_idDn_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      }

      WZHist_->Scale( 0.65*62.78*lumi*1000./ ((TH1D*)WZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZHist_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
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

      TLegend* legend = new TLegend(0.15,0.7,0.4,0.93);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(WZHist_,"WZ");
      legend->AddEntry(ZZHist_,"ZZ");
      legend->AddEntry(FFHist_,"Data-driven");
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

          r0.push_back(ratio->GetBinContent(idx));

          double rIdUp = dataHist_->GetBinContent(idx)/idDn->GetBinContent(idx) - ratio->GetBinContent(idx);
          double rIdDn = ratio->GetBinContent(idx) - dataHist_->GetBinContent(idx)/idUp->GetBinContent(idx);
          double rFFup = dataHist_->GetBinContent(idx)/FFHist_ffDn_->GetBinContent(idx) - ratio->GetBinContent(idx);
          double rFFdn = ratio->GetBinContent(idx) - dataHist_->GetBinContent(idx)/FFHist_ffUp_->GetBinContent(idx);
          /*add rFFupM & rFFdnM*/

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
        //rgr->Draw("2");
      }

      delete tmpHist;
    };
  };

  canvas_2->cd();
  auto aloader3P1F = HistLoader3P1F(datafile,WZfile,ZZfile);

  aloader3P1F.load("3P1F_CR_llll_invM","2P2F_CR_llll_invM_xFF");
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_llll_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1ll2_dr","2P2F_CR_ll1ll2_dr_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll1ll2_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_invM","2P2F_CR_ll1_invM_xFF");
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll1_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_dr","2P2F_CR_ll1_dr_xFF");
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll1_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_invM","2P2F_CR_ll2_invM_xFF");
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll2_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_dr","2P2F_CR_ll2_dr_xFF");
  aloader3P1F.compareAndRatio(p1,p2,2);
  SaveAs(canvas_2,"REMuFF_3P1F_CR_ll2_dr.png",p1);

  // 4P0F CR

  class HistLoader4P0F : public HistLoaderBase {
  public:
    HistLoader4P0F(TFile* adatafile, TFile* aWZfile, TFile* aZZfile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile,aWZfile,aZZfile),
      sigFiles_(sigFiles) {}

    ~HistLoader4P0F() {
      if (dataHist_)
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FF2P2FHist_;

      if (FF3P1FHist_ffUp_) {
        delete ZZ3P1FHist_idUp_, ZZ4P0FHist_idUp_, ZZ3P1FHist_idDn_, ZZ4P0FHist_idDn_;
        delete FF3P1FHist_ffUp_, FF3P1FHist_ffDn_, FF2P2FHist_ffUp_, FF2P2FHist_ffDn_, ZZ3P1FHist_ffUp_, ZZ3P1FHist_ffDn_;
        delete FF3P1FHist_ffUpM_, FF3P1FHist_ffDnM_, FF2P2FHist_ffUpM_, FF2P2FHist_ffDnM_, ZZ3P1FHist_ffUpM_, ZZ3P1FHist_ffDnM_;

        FF3P1FHist_ffUp_ = nullptr;
      }
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* ZZ3P1FHist_ = nullptr;
    TH1D* ZZ4P0FHist_ = nullptr;
    TH1D* FF3P1FHist_ = nullptr;
    TH1D* FF2P2FHist_ = nullptr;

    TH1D* ZZ3P1FHist_idUp_ = nullptr;
    TH1D* ZZ4P0FHist_idUp_ = nullptr;
    TH1D* ZZ3P1FHist_idDn_ = nullptr;
    TH1D* ZZ4P0FHist_idDn_ = nullptr;
    TH1D* ZZ3P1FHist_ffUp_ = nullptr;
    TH1D* ZZ3P1FHist_ffDn_ = nullptr;
    TH1D* FF3P1FHist_ffUp_ = nullptr;
    TH1D* FF3P1FHist_ffDn_ = nullptr;
    TH1D* FF2P2FHist_ffUp_ = nullptr;
    TH1D* FF2P2FHist_ffDn_ = nullptr;
    TH1D* ZZ3P1FHist_ffUpM_ = nullptr;
    TH1D* ZZ3P1FHist_ffDnM_ = nullptr;
    TH1D* FF3P1FHist_ffUpM_ = nullptr;
    TH1D* FF3P1FHist_ffDnM_ = nullptr;
    TH1D* FF2P2FHist_ffUpM_ = nullptr;
    TH1D* FF2P2FHist_ffDnM_ = nullptr;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<TH1D*> sigHist_idUp_;
    std::vector<TH1D*> sigHist_idDn_;

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

    void load(const std::string& name) {
      if (dataHist_) {
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FF2P2FHist_;
      }

      if (FF3P1FHist_ffUp_) {
        delete ZZ3P1FHist_idUp_, ZZ4P0FHist_idUp_, ZZ3P1FHist_idDn_, ZZ4P0FHist_idDn_;
        delete FF3P1FHist_ffUp_, FF3P1FHist_ffDn_, FF2P2FHist_ffUp_, FF2P2FHist_ffDn_, ZZ3P1FHist_ffUp_, ZZ3P1FHist_ffDn_;
        delete FF3P1FHist_ffUpM_, FF3P1FHist_ffDnM_, FF2P2FHist_ffUpM_, FF2P2FHist_ffDnM_, ZZ3P1FHist_ffUpM_, ZZ3P1FHist_ffDnM_;

        FF3P1FHist_ffUp_ = nullptr;
      }

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/4P0F_CR_")+name).c_str() )->Clone();
      ZZ3P1FHist_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      ZZ4P0FHist_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/4P0F_CR_")+name).c_str() )->Clone();
      FF3P1FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      FF2P2FHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2").c_str() )->Clone();

      const double lumi = valLumi;
      sigHist_.clear();
      sigHist_idUp_.clear();
      sigHist_idDn_.clear();

      for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
        sigHist_.push_back( (TH1D*)sigFiles_.at(idx).file_->Get( (std::string(anlyzrMC+"/4P0F_CR_")+name).c_str() )->Clone() );
        sigHist_.back()->Scale( lumi*1000.*0.001 / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) ); // h_sum2E2M
        sigHist_.back()->SetLineWidth(2);
        int color = (idx > 2) ? kRed : kOrange;
        sigHist_.back()->SetLineColor(color);
      }

      ZZ3P1FHist_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZ4P0FHist_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
      ZZ4P0FHist_->SetFillColor(kBlue-4);
      ZZ4P0FHist_->SetLineWidth(0);
      FF3P1FHist_->SetFillColor(38);
      FF3P1FHist_->SetLineWidth(0);
      FF2P2FHist_->SetFillColor(33);
      FF2P2FHist_->SetLineWidth(0);

      if ( TString(name).Contains("llll_invM") ) {
        ZZ3P1FHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF_idUp").c_str() )->Clone();
        ZZ4P0FHist_idUp_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone();
        ZZ3P1FHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF_idDn").c_str() )->Clone();
        ZZ4P0FHist_idDn_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone();
        ZZ3P1FHist_ffUp_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF_ffUpE").c_str() )->Clone();
        ZZ3P1FHist_ffDn_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF_ffDnE").c_str() )->Clone();
        FF3P1FHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffUpE").c_str() )->Clone();
        FF2P2FHist_ffUp_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffUpE").c_str() )->Clone();
        FF3P1FHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffDnE").c_str() )->Clone();
        FF2P2FHist_ffDn_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffDnE").c_str() )->Clone();
        ZZ3P1FHist_ffUpM_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF_ffUpM").c_str() )->Clone();
        ZZ3P1FHist_ffDnM_ = (TH1D*)ZZfile_->Get( (std::string(anlyzrMC+"/3P1F_CR_")+name+"_xFF_ffDnM").c_str() )->Clone();
        FF3P1FHist_ffUpM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffUpM").c_str() )->Clone();
        FF2P2FHist_ffUpM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffUpM").c_str() )->Clone();
        FF3P1FHist_ffDnM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/3P1F_CR_")+name+"_xFF_ffDnM").c_str() )->Clone();
        FF2P2FHist_ffDnM_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/2P2F_CR_")+name+"_xFF2_ffDnM").c_str() )->Clone();

        ZZ3P1FHist_idUp_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ4P0FHist_idUp_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ3P1FHist_idDn_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ4P0FHist_idDn_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ3P1FHist_ffUp_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ3P1FHist_ffDn_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ3P1FHist_ffUpM_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );
        ZZ3P1FHist_ffDnM_->Scale( 13.81*lumi*1000./ ((TH1D*)ZZfile_->Get("evtCounter/h_sumW"))->GetBinContent(1) );

        for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
          sigHist_idUp_.push_back( (TH1D*)sigFiles_.at(idx).file_->Get( (std::string(anlyzrMC+"/4P0F_CR_")+name+"_idUp").c_str() )->Clone() );
          sigHist_idUp_.back()->Scale( lumi*1000.*0.001 / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );
          sigHist_idDn_.push_back( (TH1D*)sigFiles_.at(idx).file_->Get( (std::string(anlyzrMC+"/4P0F_CR_")+name+"_idDn").c_str() )->Clone() );
          sigHist_idDn_.back()->Scale( lumi*1000.*0.001 / ( (TH1D*)sigFiles_.at(idx).file_->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );
        }
      }
    };

    void compare(TPad* pad, int rebin=1) {
      if ( FF2P2FHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FF2P2FHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FF2P2FHist_->Rebin( FF2P2FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FF3P1FHist_->Rebin( FF3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        ZZ3P1FHist_->Rebin( ZZ3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (FF3P1FHist_ffUp_) {
          FF2P2FHist_ffUp_->Rebin( FF2P2FHist_ffUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF2P2FHist_ffDn_->Rebin( FF2P2FHist_ffDn_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF3P1FHist_ffUp_->Rebin( FF3P1FHist_ffUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF3P1FHist_ffDn_->Rebin( FF3P1FHist_ffDn_->GetNbinsX()/dataHist_->GetNbinsX() );
          ZZ3P1FHist_idUp_->Rebin( ZZ3P1FHist_idUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          ZZ3P1FHist_idDn_->Rebin( ZZ3P1FHist_idDn_->GetNbinsX()/dataHist_->GetNbinsX() );
          ZZ3P1FHist_ffUp_->Rebin( ZZ3P1FHist_ffUp_->GetNbinsX()/dataHist_->GetNbinsX() );
          ZZ3P1FHist_ffDn_->Rebin( ZZ3P1FHist_ffDn_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF2P2FHist_ffUpM_->Rebin( FF2P2FHist_ffUpM_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF2P2FHist_ffDnM_->Rebin( FF2P2FHist_ffDnM_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF3P1FHist_ffUpM_->Rebin( FF3P1FHist_ffUpM_->GetNbinsX()/dataHist_->GetNbinsX() );
          FF3P1FHist_ffDnM_->Rebin( FF3P1FHist_ffDnM_->GetNbinsX()/dataHist_->GetNbinsX() );
          ZZ3P1FHist_ffUpM_->Rebin( ZZ3P1FHist_ffUpM_->GetNbinsX()/dataHist_->GetNbinsX() );
          ZZ3P1FHist_ffDnM_->Rebin( ZZ3P1FHist_ffDnM_->GetNbinsX()/dataHist_->GetNbinsX() );
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

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(FF2P2FHist_);
      stack->Add(FF3P1Fsubtracted);
      stack->Add(ZZ4P0FHist_);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));
      tmpHist->Add((TH1*)stack->GetHists()->At(2));

      pad->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        dataHist_->GetXaxis()->SetRangeUser(100.,300.);
      }

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(FF3P1Fsubtracted,"Data-driven (3P1F)");
      legend->AddEntry(FF2P2FHist_,"Data-driven (2P2F)");
      legend->AddEntry(ZZ4P0FHist_,"ZZ");

      if ( TString(dataHist_->GetName()).Contains("llll_invM") ) {
        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);
          sigHist_idUp_.at(idx)->Rebin(rebin);
          sigHist_idDn_.at(idx)->Rebin(rebin);
          sigHist_.at(idx)->Draw("hist&same");
        }

        //legend->AddEntry(sigHist_.at(3),"H250/750/2000A1");
        legend->AddEntry(sigHist_.at(0),"H750A1");
      }

      legend->Draw();

      if (FF3P1FHist_ffUp_) {
        TH1D* FF3P1Fsubtracted_idUp = subtract(FF3P1FHist_,ZZ3P1FHist_idUp_);
        TH1D* FF3P1Fsubtracted_idDn = subtract(FF3P1FHist_,ZZ3P1FHist_idDn_);
        TH1D* FF3P1Fsubtracted_ffUp = subtract(FF3P1FHist_ffUp_,ZZ3P1FHist_ffUp_);
        TH1D* FF3P1Fsubtracted_ffDn = subtract(FF3P1FHist_ffDn_,ZZ3P1FHist_ffDn_);
        TH1D* FF3P1Fsubtracted_ffUpM = subtract(FF3P1FHist_ffUpM_,ZZ3P1FHist_ffUpM_);
        TH1D* FF3P1Fsubtracted_ffDnM = subtract(FF3P1FHist_ffDnM_,ZZ3P1FHist_ffDnM_);

        TH1D* idUp = (TH1D*)ZZ4P0FHist_idUp_->Clone();
        idUp->Add(FF2P2FHist_);
        idUp->Add(FF3P1Fsubtracted_idUp);
        TH1D* idDn = (TH1D*)ZZ4P0FHist_idDn_->Clone();
        idDn->Add(FF2P2FHist_);
        idDn->Add(FF3P1Fsubtracted_idDn);

        TH1D* ffUpAdded = (TH1D*)FF2P2FHist_ffUp_->Clone();
        TH1D* ffDnAdded = (TH1D*)FF2P2FHist_ffDn_->Clone();

        ffUpAdded->Add(FF3P1Fsubtracted_ffUp);
        ffUpAdded->Add(ZZ4P0FHist_);
        ffDnAdded->Add(FF3P1Fsubtracted_ffDn);
        ffDnAdded->Add(ZZ4P0FHist_);

        TH1D* ffUpMAdded = (TH1D*)FF2P2FHist_ffUpM_->Clone();
        TH1D* ffDnMAdded = (TH1D*)FF2P2FHist_ffDnM_->Clone();

        ffUpMAdded->Add(FF3P1Fsubtracted_ffUpM);
        ffUpMAdded->Add(ZZ4P0FHist_);
        ffDnMAdded->Add(FF3P1Fsubtracted_ffDnM);
        ffDnMAdded->Add(ZZ4P0FHist_);

        std::vector<double> x0, y0, errx, erryDn, erryUp;

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

          erryUp.push_back( std::sqrt(valIdUp*valIdUp + valFFup*valFFup + valFFupM*valFFupM) );
          erryDn.push_back( std::sqrt(valIdDn*valIdDn + valFFdn*valFFdn + valFFdnM*valFFdnM) );
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
          dir_->WriteTObject(FF2P2FHist_ffUp_,"2P2F_resolvedEleFakeFactorUp");
          dir_->WriteTObject(FF2P2FHist_ffDn_,"2P2F_resolvedEleFakeFactorDown");
          dir_->WriteTObject(FF3P1Fsubtracted_ffUpM,"3P1F_resolvedMuFakeFactorUp");
          dir_->WriteTObject(FF3P1Fsubtracted_ffDnM,"3P1F_resolvedMuFakeFactorDown");
          dir_->WriteTObject(FF2P2FHist_ffUpM_,"2P2F_resolvedMuFakeFactorUp");
          dir_->WriteTObject(FF2P2FHist_ffDnM_,"2P2F_resolvedMuFakeFactorDown");
          dir_->WriteTObject(ZZ4P0FHist_idUp_,"ZZ_modHeepIdUp");
          dir_->WriteTObject(ZZ4P0FHist_idDn_,"ZZ_modHeepIdDown");
          dir_->WriteTObject(FF3P1Fsubtracted_idUp,"3P1F_modHeepIdUp");
          dir_->WriteTObject(FF3P1Fsubtracted_idDn,"3P1F_modHeepIdDown");

          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
            dir_->WriteTObject(sigHist_idUp_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdUp");
            dir_->WriteTObject(sigHist_idDn_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdDown");
          }
        }
      }
    };
  };

  auto aloader4P0F = HistLoader4P0F(datafile,WZfile,ZZfile,std::vector<SigSample>({H750A1sample}));

  aloader4P0F.load("llll_invM");
  aloader4P0F.preparecard("REMuFF_"+era+"_datacard.root","resolvedEMu");
  aloader4P0F.compare(canvas_1,1);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_llll_invM_zoomed.png");
  aloader4P0F.close();

  p1->SetLogy(0);
  canvas_1->SetLogy(0);

  aloader4P0F.load("ll1ll2_dr");
  aloader4P0F.compare(canvas_1);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll1ll2_dr.png");

  aloader4P0F.load("ll1_invM");
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll1_invM.png");

  aloader4P0F.load("ll1_dr");
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll1_dr.png");

  aloader4P0F.load("ll2_invM");
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll2_invM.png");

  aloader4P0F.load("ll2_dr");
  aloader4P0F.compare(canvas_1,2);
  SaveAs(canvas_1,"REMuFF_4P0F_CR_ll2_dr.png");

  return;
}
