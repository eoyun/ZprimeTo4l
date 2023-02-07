#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runREFF() {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  lumi_13TeV = "16.8 fb^{-1}";

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

  if( iPos==0 )
    relPosX = 0.12;

  int W = 800;
  int H = 600;

  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  TFile* datafile = new TFile("ResolvedEleCR_20UL16_data.root","READ");
  TFile* WZfile = new TFile("ResolvedEleCR_20UL16_WZ.root","READ");
  TFile* ZZfile = new TFile("ResolvedEleCR_20UL16_ZZ.root","READ");

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
  p1->SetBottomMargin(0.005);
  p1->SetLeftMargin( L/W );
  p1->SetRightMargin( R/W );
  p1->Draw();

  TPad* p2 = new TPad("p2","",0,0,1,0.3);
  p2->SetTickx(0);
  p2->SetTicky(0);
  p2->SetTopMargin(0.03);
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
        delete dataHist_, WZHist_, ZZHist_, FFHist_PFPF_, FFHist_PPFF_;
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* WZHist_ = nullptr;
    TH1D* ZZHist_ = nullptr;
    TH1D* FFHist_PFPF_ = nullptr;
    TH1D* FFHist_PPFF_ = nullptr;

  public:
    void load(const std::string& nameNum, const std::string& namePFPF, const std::string& namePPFF) {
      if (dataHist_)
        delete dataHist_, WZHist_, ZZHist_, FFHist_PFPF_, FFHist_PPFF_;

      dataHist_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/")+nameNum).c_str() )->Clone();
      WZHist_ = (TH1D*)WZfile_->Get( (std::string("resolvedEleCRanalyzer/")+nameNum).c_str() )->Clone();
      ZZHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer/")+nameNum).c_str() )->Clone();
      FFHist_PFPF_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/")+namePFPF).c_str() )->Clone();
      FFHist_PPFF_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/")+namePPFF).c_str() )->Clone();

      const double lumi = 16.8;
      WZHist_->Scale( 5.213*lumi*1000./ ((TH1D*)WZfile_->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );
      ZZHist_->Scale( 1.325*lumi*1000./ ((TH1D*)ZZfile_->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );
      WZHist_->SetFillColor(kViolet+1);
      ZZHist_->SetFillColor(kBlue+1);
      WZHist_->SetLineWidth(0);
      ZZHist_->SetLineWidth(0);
      FFHist_PFPF_->SetFillColor(33);
      FFHist_PFPF_->SetLineWidth(0);
      FFHist_PPFF_->SetFillColor(36);
      FFHist_PPFF_->SetLineWidth(0);
    };

    void compareAndRatio(TPad* padUp, TPad* padDn) {
      if ( FFHist_PFPF_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FFHist_PFPF_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FFHist_PFPF_->Rebin( FFHist_PFPF_->GetNbinsX()/dataHist_->GetNbinsX() );
        FFHist_PPFF_->Rebin( FFHist_PPFF_->GetNbinsX()/dataHist_->GetNbinsX() );
      }

      WZHist_->Rebin(2);
      ZZHist_->Rebin(2);
      FFHist_PFPF_->Rebin(2);
      FFHist_PPFF_->Rebin(2);
      dataHist_->Rebin(2);

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(ZZHist_);
      stack->Add(WZHist_);
      stack->Add(FFHist_PFPF_);
      stack->Add(FFHist_PPFF_);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));
      tmpHist->Add((TH1*)stack->GetHists()->At(2));
      tmpHist->Add((TH1*)stack->GetHists()->At(3));

      padUp->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));
      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(WZHist_,"WZ");
      legend->AddEntry(ZZHist_,"ZZ");
      legend->AddEntry(FFHist_PFPF_,"Data-driven (PFPF)");
      legend->AddEntry(FFHist_PPFF_,"Data-driven (PPFF)");
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
      ratio->GetYaxis()->SetRangeUser(0.,2.);
      ratio->GetXaxis()->SetTitleSize(0.12);
      ratio->GetXaxis()->SetTitleOffset(0.75);
      ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
      ratio->SetLineColor(kBlack);

      padDn->cd();
      ratio->Draw("E1");

      delete tmpHist;
    };
  };

  canvas_2->cd();
  auto aloader3P1F = HistLoader3P1F(datafile,WZfile,ZZfile);

  aloader3P1F.load("3P1F_CR_llll_invM","2P2F_CR_PFPF_llll_invM_xFF","2P2F_CR_PPFF_llll_invM_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_llll_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1ll2_dr","2P2F_CR_PFPF_ll1ll2_dr_xFF","2P2F_CR_PPFF_ll1ll2_dr_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll1ll2_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_invM","2P2F_CR_PFPF_ll1_invM_xFF","2P2F_CR_PPFF_ll1_invM_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll1_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll1_dr","2P2F_CR_PFPF_ll1_dr_xFF","2P2F_CR_PPFF_ll1_dr_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll1_dr.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_invM","2P2F_CR_PFPF_ll2_invM_xFF","2P2F_CR_PPFF_ll2_invM_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll2_invM.png",p1);

  aloader3P1F.load("3P1F_CR_ll2_dr","2P2F_CR_PFPF_ll2_dr_xFF","2P2F_CR_PPFF_ll2_dr_xFF");
  aloader3P1F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_3P1F_CR_ll2_dr.png",p1);

  // 4P0F CR

  class HistLoader4P0F : public HistLoaderBase {
  public:
    HistLoader4P0F(TFile* adatafile, TFile* aWZfile, TFile* aZZfile) : HistLoaderBase(adatafile,aWZfile,aZZfile) {}

    ~HistLoader4P0F() {
      if (dataHist_)
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FFPFPFHist_, FFPPFFHist_;
    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* ZZ3P1FHist_ = nullptr;
    TH1D* ZZ4P0FHist_ = nullptr;
    TH1D* FF3P1FHist_ = nullptr;
    TH1D* FFPFPFHist_ = nullptr;
    TH1D* FFPPFFHist_ = nullptr;

  public:
    void load(const std::string& name) {
      if (dataHist_)
        delete dataHist_, ZZ3P1FHist_, ZZ4P0FHist_, FF3P1FHist_, FFPFPFHist_, FFPPFFHist_;

      dataHist_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/4P0F_CR_")+name).c_str() )->Clone();
      ZZ3P1FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      ZZ4P0FHist_ = (TH1D*)ZZfile_->Get( (std::string("resolvedEleCRanalyzer/4P0F_CR_")+name).c_str() )->Clone();
      FF3P1FHist_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/3P1F_CR_")+name+"_xFF").c_str() )->Clone();
      FFPFPFHist_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/2P2F_CR_PFPF_")+name+"_xFF2").c_str() )->Clone();
      FFPPFFHist_ = (TH1D*)datafile_->Get( (std::string("resolvedEleCRanalyzerData/2P2F_CR_PPFF_")+name+"_xFF2").c_str() )->Clone();

      const double lumi = 16.8;
      ZZ3P1FHist_->Scale( 1.325*lumi*1000./ ((TH1D*)ZZfile_->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );
      ZZ4P0FHist_->Scale( 1.325*lumi*1000./ ((TH1D*)ZZfile_->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );
      ZZ4P0FHist_->SetFillColor(kBlue+1);
      ZZ4P0FHist_->SetLineWidth(0);
      FF3P1FHist_->SetFillColor(38);
      FF3P1FHist_->SetLineWidth(0);
      FFPFPFHist_->SetFillColor(33);
      FFPFPFHist_->SetLineWidth(0);
      FFPPFFHist_->SetFillColor(36);
      FFPPFFHist_->SetLineWidth(0);
    };

    void compare(TPad* pad) {
      if ( FFPFPFHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FFPFPFHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FFPFPFHist_->Rebin( FFPFPFHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FFPPFFHist_->Rebin( FFPPFFHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FF3P1FHist_->Rebin( FF3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        ZZ3P1FHist_->Rebin( ZZ3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
      }

      TH1D* FF3P1Fsubtracted = (TH1D*)FF3P1FHist_->Clone();

      for (unsigned ibin = 0; ibin < FF3P1Fsubtracted->GetNbinsX()+2; ibin++) {
        double val = FF3P1FHist_->GetBinContent(ibin) - ZZ3P1FHist_->GetBinContent(ibin);
        auto square = [] (const double a) { return a*a; };
        double err = std::sqrt( square(FF3P1FHist_->GetBinError(ibin)) + square(ZZ3P1FHist_->GetBinError(ibin)) );
        FF3P1Fsubtracted->SetBinContent(ibin,val);
        FF3P1Fsubtracted->SetBinError(ibin,err);
      }

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(ZZ4P0FHist_);
      stack->Add(FFPFPFHist_);
      stack->Add(FFPPFFHist_);
      stack->Add(FF3P1Fsubtracted);

      pad->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));
      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(FF3P1Fsubtracted,"Data-driven (3P1F)");
      legend->AddEntry(FFPFPFHist_,"Data-driven (PFPF)");
      legend->AddEntry(FFPPFFHist_,"Data-driven (PPFF)");
      legend->AddEntry(ZZ4P0FHist_,"ZZ");
      legend->Draw();

      delete stack;
    };

    void compareAndRatio(TPad* padUp, TPad* padDn) {
      if ( FFPFPFHist_->GetNbinsX() % dataHist_->GetNbinsX() == 0 &&
           FFPFPFHist_->GetNbinsX() != dataHist_->GetNbinsX() ) {
        FFPFPFHist_->Rebin( FFPFPFHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FFPPFFHist_->Rebin( FFPPFFHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        FF3P1FHist_->Rebin( FF3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
        ZZ3P1FHist_->Rebin( ZZ3P1FHist_->GetNbinsX()/dataHist_->GetNbinsX() );
      }

      ZZ3P1FHist_->Rebin(2);
      ZZ4P0FHist_->Rebin(2);
      FF3P1FHist_->Rebin(2);
      FFPFPFHist_->Rebin(2);
      FFPPFFHist_->Rebin(2);
      dataHist_->Rebin(2);

      TH1D* FF3P1Fsubtracted = (TH1D*)FF3P1FHist_->Clone();

      for (unsigned ibin = 0; ibin < FF3P1Fsubtracted->GetNbinsX()+2; ibin++) {
        double val = FF3P1FHist_->GetBinContent(ibin) - ZZ3P1FHist_->GetBinContent(ibin);
        auto square = [] (const double a) { return a*a; };
        double err = std::sqrt( square(FF3P1FHist_->GetBinError(ibin)) + square(ZZ3P1FHist_->GetBinError(ibin)) );
        FF3P1Fsubtracted->SetBinContent(ibin,val);
        FF3P1Fsubtracted->SetBinError(ibin,err);
      }

      THStack* stack = new THStack( (std::string(dataHist_->GetName())+"_stack").c_str() , dataHist_->GetTitle() );
      stack->Add(ZZ4P0FHist_);
      stack->Add(FFPFPFHist_);
      stack->Add(FFPPFFHist_);
      stack->Add(FF3P1Fsubtracted);

      TH1* tmpHist = (TH1*)stack->GetHists()->At(0)->Clone();
      tmpHist->Add((TH1*)stack->GetHists()->At(1));
      tmpHist->Add((TH1*)stack->GetHists()->At(2));
      tmpHist->Add((TH1*)stack->GetHists()->At(3));

      padUp->cd();
      dataHist_->SetMaximum(1.5*std::max(dataHist_->GetMaximum(),stack->GetMaximum()));
      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->Draw("E1");
      stack->Draw("hist&same");
      dataHist_->Draw("E1&same");

      TLegend* legend = new TLegend(0.65,0.7,0.9,0.9);
      legend->SetBorderSize(0);
      legend->AddEntry(dataHist_,"Data");
      legend->AddEntry(FF3P1Fsubtracted,"Data-driven (3P1F)");
      legend->AddEntry(FFPFPFHist_,"Data-driven (PFPF)");
      legend->AddEntry(FFPPFFHist_,"Data-driven (PPFF)");
      legend->AddEntry(ZZ4P0FHist_,"ZZ");
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
      ratio->GetYaxis()->SetRangeUser(0.,2.);
      ratio->GetXaxis()->SetTitleSize(0.12);
      ratio->GetXaxis()->SetTitleOffset(0.75);
      ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
      ratio->SetLineColor(kBlack);

      padDn->cd();
      ratio->Draw("E1");

      delete tmpHist;
    };
  };

  auto aloader4P0F = HistLoader4P0F(datafile,WZfile,ZZfile);

  aloader4P0F.load("llll_invM");
  aloader4P0F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_4P0F_CR_llll_invM.png",p1);

  aloader4P0F.load("ll1ll2_dr");
  aloader4P0F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_4P0F_CR_ll1ll2_dr.png",p1);

  aloader4P0F.load("ll1_invM");
  aloader4P0F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_4P0F_CR_ll1_invM.png",p1);

  aloader4P0F.load("ll1_dr");
  aloader4P0F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_4P0F_CR_ll1_dr.png",p1);

  aloader4P0F.load("ll2_invM");
  aloader4P0F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_4P0F_CR_ll2_invM.png",p1);

  aloader4P0F.load("ll2_dr");
  aloader4P0F.compareAndRatio(p1,p2);
  SaveAs(canvas_2,"REFF_4P0F_CR_ll2_dr.png",p1);

  return;
}
