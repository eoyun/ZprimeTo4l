#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void combineFlavor() {
  TString era = "run2";
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "";  // default extra text is "Preliminary"

  lumi_sqrtS = "";
  lumi_13TeV = "138 fb^{-1}";

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 11;

  if( iPos==0 )
    relPosX = 0.12;

  int W = 800;
  int H = 800;

  int H_ref = 800;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref; // 0.1
  float R = 0.05*W_ref; // 0.02

  auto* canvas_2 = new TCanvas("canvas_2","canvas_2",800,800,W,H);
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

  auto SaveAs = [&] (TCanvas* canvas, const std::string& name, TPad* pad = nullptr) {
    canvas->Update();

    // writing the lumi information and the CMS "logo"
    CMS_lumi( canvas, iPeriod, iPos );

    if (pad) {
      pad->cd();
      pad->RedrawAxis();
      pad->GetFrame()->Draw();
    } else {
      canvas->Update();
      canvas->RedrawAxis();
      canvas->GetFrame()->Draw();
    }

    canvas->SaveAs(name.c_str());
  };

  class SignalRegion {
  public:
    TH1D* data_;
    TH1D* sig_;

    std::vector<TH1D*> nominals_;
    std::vector<std::map<TString,TH1D*>> up_;
    std::vector<std::map<TString,TH1D*>> dn_;

  protected:
    TFile* file_;
    TString strSR_;

    TString signame_;
    std::vector<int> colors_;
    std::vector<TString> bkgs_;
    std::vector<TString> systs_;
    int sigColor_ = kRed;
    int rebin_ = 1;

    void init() {
      data_ = (TH1D*)file_->Get(strSR_+"/data_obs");
      sig_ = (TH1D*)file_->Get(strSR_+"/"+signame_);
      sig_->SetLineColor(sigColor_);
      sig_->SetLineWidth(2);
      sig_->Scale(100.);

      data_->Rebin(rebin_);
      sig_->Rebin(rebin_);

      for (unsigned idx=0; idx<bkgs_.size(); idx++) {
        const auto& bkgname = bkgs_.at(idx);
        TH1D* nominal = (TH1D*)file_->Get(strSR_+"/"+bkgname);
        nominal->SetFillColor(colors_.at(idx));
        nominal->SetMarkerColor(colors_.at(idx));
        nominal->SetLineWidth(0);
        nominal->Rebin(rebin_);

        nominals_.push_back(nominal);

        std::map<TString,TH1D*> up;
        std::map<TString,TH1D*> dn;

        for (const auto& systname : systs_) {
          TH1D* h1 = (TH1D*)file_->Get(strSR_+"/"+bkgname+"_"+systname+"Up");
          TH1D* h2 = (TH1D*)file_->Get(strSR_+"/"+bkgname+"_"+systname+"Down");

          if (h1) {
            h1->Rebin(rebin_);
            h2->Rebin(rebin_);
            up[systname] = h1;
            dn[systname] = h2;
          } else {
            up[systname] = (TH1D*)nominal->Clone();
            dn[systname] = (TH1D*)nominal->Clone();
          }
        }

        up_.push_back(up);
        dn_.push_back(dn);
      }
    }

  public:
    SignalRegion(TFile* afile, TString astr)
    : file_(afile), strSR_(astr) {}

    ~SignalRegion()=default;

    void add(SignalRegion& other) {
      data_->Add(other.data_);
      sig_->Add(other.sig_);

      for (unsigned idx=0; idx<nominals_.size(); idx++) {
        nominals_.at(idx)->Add(other.nominals_.at(idx));

        for (const auto& systname : systs_) {
          up_.at(idx).at(systname)->Add(other.up_.at(idx).at(systname));
          dn_.at(idx).at(systname)->Add(other.dn_.at(idx).at(systname));
        }
      }
    }

    std::unique_ptr<THStack> returnStack() {
      auto astack = std::make_unique<THStack>("stack",data_->GetTitle());

      for (auto* x : nominals_)
        astack->Add(x);

      return std::move(astack);
    }

    std::unique_ptr<TH1D> returnNominal() {
      auto x = std::unique_ptr<TH1D>((TH1D*)nominals_.front()->Clone());

      for (unsigned idx=1; idx<nominals_.size(); idx++)
        x->Add(nominals_.at(idx));

      return std::move(x);
    }

    std::unique_ptr<TH1D> returnRatio() {
      auto ratio = std::unique_ptr<TH1D>((TH1D*)data_->Clone());
      auto norm = returnNominal();
      ratio->Divide(norm.get());

      ratio->SetStats(0);
      ratio->SetTitle("");
      ratio->GetYaxis()->SetTitle("Obs/Exp");
      ratio->GetYaxis()->SetTitleSize(0.14);
      ratio->GetYaxis()->SetTitleOffset(0.34);
      ratio->GetXaxis()->SetLabelSize(0.12);
      ratio->GetYaxis()->SetLabelSize(0.12);
      ratio->GetXaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetLabelOffset(0.005);
      ratio->GetYaxis()->SetRangeUser(0.2,1.8);
      ratio->GetXaxis()->SetTitleSize(0.15);
      ratio->GetXaxis()->SetTitleOffset(0.8);
      ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
      ratio->SetLineColor(kBlack);

      return std::move(ratio);
    }

    std::unique_ptr<TH1D> returnVarianceUp(TString systname) {
      auto x = std::unique_ptr<TH1D>((TH1D*)up_.front().at(systname)->Clone());

      for (unsigned idx=1; idx<up_.size(); idx++)
        x->Add((TH1D*)up_.at(idx).at(systname));

      return std::move(x);
    }

    std::unique_ptr<TH1D> returnVarianceDn(TString systname) {
      auto x = std::unique_ptr<TH1D>((TH1D*)dn_.front().at(systname)->Clone());

      for (unsigned idx=1; idx<dn_.size(); idx++)
        x->Add((TH1D*)dn_.at(idx).at(systname));

      return std::move(x);
    }

    std::unique_ptr<TGraphAsymmErrors> returnError() {
      auto nominal = returnNominal();
      std::vector<double> x0, y0, errx, erryDn, erryUp;

      std::vector<std::unique_ptr<TH1D>> up, dn;

      for (const auto& systname : systs_) {
        up.push_back(std::move(returnVarianceUp(systname)));
        dn.push_back(std::move(returnVarianceDn(systname)));
      }

      for (int idx=1; idx<=nominal->GetNbinsX(); idx++) {
        x0.push_back(nominal->GetBinCenter(idx));
        y0.push_back(nominal->GetBinContent(idx));
        errx.push_back(nominal->GetBinWidth(idx)/2.);

        double errUpSqr = 0.;
        double errDnSqr = 0.;

        for (unsigned isyst=0; isyst<up.size(); isyst++) {
          double valUp = up.at(isyst)->GetBinContent(idx) - nominal->GetBinContent(idx);
          double valDn = dn.at(isyst)->GetBinContent(idx) - nominal->GetBinContent(idx);

          errUpSqr += valUp*valUp;
          errDnSqr += valDn*valDn;
        }

        erryUp.push_back(std::sqrt(errUpSqr));
        erryDn.push_back(std::sqrt(errDnSqr));        
      }

      auto gr = std::make_unique<TGraphAsymmErrors>(nominal->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+3);
      gr->SetLineColor(kGray+3);
      gr->SetFillStyle(3004);

      return std::move(gr);
    }

    std::unique_ptr<TGraphAsymmErrors> returnRatioError() {
      auto nominal = returnNominal();
      std::vector<double> x0, y0, errx, erryDn, erryUp;

      std::vector<std::unique_ptr<TH1D>> up, dn;

      for (const auto& systname : systs_) {
        up.push_back(std::move(returnVarianceUp(systname)));
        dn.push_back(std::move(returnVarianceDn(systname)));
      }

      for (int idx=1; idx<=nominal->GetNbinsX(); idx++) {
        x0.push_back(nominal->GetBinCenter(idx));
        y0.push_back(1.);
        errx.push_back(nominal->GetBinWidth(idx)/2.);

        double errUpSqr = 0.;
        double errDnSqr = 0.;
        double valNom = nominal->GetBinContent(idx);

        for (unsigned isyst=0; isyst<up.size(); isyst++) {
          double valUp = up.at(isyst)->GetBinContent(idx) - nominal->GetBinContent(idx);
          double valDn = dn.at(isyst)->GetBinContent(idx) - nominal->GetBinContent(idx);

          errUpSqr += valUp*valUp;
          errDnSqr += valDn*valDn;
        }

        erryUp.push_back(valNom > 0. ? std::sqrt(errUpSqr)/valNom : 0.);
        erryDn.push_back(valNom > 0. ? std::sqrt(errDnSqr)/valNom : 0.);        
      }

      auto gr = std::make_unique<TGraphAsymmErrors>(nominal->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+3);
      gr->SetLineColor(kGray+3);
      gr->SetFillStyle(3004);

      return std::move(gr);
    }
  };

  class ResolvedSR : public SignalRegion {
  public:
    ResolvedSR(TFile* afile, TString astr)
    : SignalRegion(afile,astr) {

      signame_ = "H750A100";
      colors_ = {TColor::GetColor("#9c9ca1"),TColor::GetColor("#92dadd"),TColor::GetColor("#5790fc")};
      bkgs_ = {"Nonprompt","NonpromptDR03","ZZ"};
      systs_ = {"resolvedEleFakeFactor",
                "resolvedMuFakeFactor",
                "modHeepId",
                "muMomentumScale",
                "muMomentumSmear",
                "highPtId",
                "muLooseIso",
                "muTrig",
                "elTrig",
                "muReco",
                "elEnergyScale",
                "elEnergySigma",
                "muBoostIso",
                "elReco",
                "PUrwgt",
                "prefire",
                "ZZnorm"};
      rebin_ = 4;

      init();

      unsigned ZZidx = 1;
      TH1D* ZZup = (TH1D*)nominals_.at(ZZidx)->Clone();
      TH1D* ZZdn = (TH1D*)nominals_.at(ZZidx)->Clone();
      ZZup->Scale(1.1);
      ZZdn->Scale(0.9);
      (up_.at(ZZidx))["ZZnorm"] = ZZup;
      (dn_.at(ZZidx))["ZZnorm"] = ZZdn;

      data_->GetXaxis()->SetRangeUser(0.,1500.);
      data_->GetXaxis()->SetTitle("M(4l) [GeV]");
      data_->GetYaxis()->SetLabelSize(0.06);
      data_->SetLineWidth(2);
      data_->SetLineColor(kBlack);
      data_->SetMinimum(0.001);
    }
  };

  class MergedEleSR : public SignalRegion {
  public:
    MergedEleSR(TFile* afile, TString astr)
    : SignalRegion(afile,astr) {

      signame_ = "H750Z1";
      colors_ = {TColor::GetColor("#f89c20"),TColor::GetColor("#9c9ca1")};
      bkgs_ = {"OS","SS"};
      systs_ = {"mergedEleFakeFactorSS",
                "mergedEleFakeFactorOS",
                "modHeepId",
                "mergedEleEnCorr",
                "muMomentumScale",
                "elTrig",
                "elReco"};

      rebin_ = 10;

      init();

      unsigned OSidx = 0;
      TH1D* OSup = (TH1D*)nominals_.at(OSidx)->Clone();
      TH1D* OSdn = (TH1D*)nominals_.at(OSidx)->Clone();
      OSup->Scale(1.2);
      OSdn->Scale(0.8);
      (up_.at(OSidx))["OSnorm"] = OSup;
      (dn_.at(OSidx))["OSnorm"] = OSdn;

      data_->GetXaxis()->SetRangeUser(0.,1500.);
      data_->GetXaxis()->SetTitle("M(e_{ME}ll) [GeV]");
      data_->SetLineWidth(2);
      data_->SetLineColor(kBlack);
      data_->GetXaxis()->SetTitleSize(0.07);
      data_->GetYaxis()->SetLabelSize(0.06);
      data_->SetMinimum(0.001);
    }
  };

  class CleanedMuSR : public SignalRegion {
  public:
    CleanedMuSR(TFile* afile, TString astr)
    : SignalRegion(afile,astr) {

      signame_ = "H2000Z1";
      colors_ = {kGray};
      bkgs_ = {"MM"};
      systs_ = {"mergedMuFakeFactor",
                "JES"};

      rebin_ = 10;

      init();

      data_->GetXaxis()->SetRangeUser(0.,2500.);
      data_->GetXaxis()->SetTitle("M_{T}(ll#mup_{T}^{miss}) [GeV]");
      data_->SetLineWidth(2);
      data_->SetLineColor(kBlack);
      data_->GetYaxis()->SetLabelSize(0.06);
      data_->SetMinimum(0.001);
    }
  };

  TFile* eeeeFile = new TFile("REFF_run2_datacard.root","READ");
  TFile* eemmFile = new TFile("REMuFF_run2_datacard.root","READ");
  TFile* mmmmFile = new TFile("RMFF_run2_datacard.root","READ");

  auto resolvedSR = ResolvedSR(eemmFile,"resolvedEMu");
  auto eeeeSR = ResolvedSR(eeeeFile,"resolvedEle");
  auto mmmmSR = ResolvedSR(mmmmFile,"resolvedMu");

  resolvedSR.add(eeeeSR);
  resolvedSR.add(mmmmSR);

  auto resolvedStack = resolvedSR.returnStack();
  auto resolvedRatio = resolvedSR.returnRatio();
  auto resolvedError = resolvedSR.returnError();
  auto resolvedRatioError = resolvedSR.returnRatioError();

  p1->cd();
  p1->SetLogy();
  resolvedSR.data_->SetMaximum(700.*std::max(resolvedSR.data_->GetMaximum(),resolvedSR.returnNominal()->GetMaximum()));
  resolvedSR.data_->SetMinimum(0.2);
  resolvedSR.data_->GetYaxis()->SetTitle("Events");
  resolvedSR.data_->GetYaxis()->SetTitleSize(0.08);
  resolvedSR.data_->GetYaxis()->SetTitleOffset(0.7);
  resolvedSR.data_->Draw("E1");
  resolvedStack->Draw("hist&same");
  resolvedSR.data_->Draw("E1&same");
  resolvedSR.sig_->Draw("hist&same");
  resolvedError->Draw("2");

  TLegend* resolvedLegend = new TLegend(0.4,0.52,0.95,0.93);
  resolvedLegend->SetMargin(0.2);
  resolvedLegend->SetBorderSize(0);
  resolvedLegend->SetFillStyle(0);
  resolvedLegend->AddEntry(resolvedSR.data_,"Data");
  resolvedLegend->AddEntry(resolvedSR.nominals_.at(0),"Nonprompt (#Delta R > 0.3)");
  resolvedLegend->AddEntry(resolvedSR.nominals_.at(1),"Nonprompt (#Delta R < 0.3)");
  resolvedLegend->AddEntry(resolvedSR.nominals_.at(2),"ZZ");
  resolvedLegend->AddEntry(resolvedSR.sig_,"M_{X} = 750 GeV, M_{Y} = 100 GeV");
  resolvedLegend->SetTextSize(0.055);

  resolvedLegend->Draw();

  p2->cd();
  resolvedRatio->Draw("E1");
  resolvedRatioError->Draw("2");

  p1->cd();

  //SaveAs(canvas_2,"resolved.pdf",p1);

  //return;

  p1->SetLogy(0);

  TFile* eeeFile = new TFile("ME3E_run2_datacard.root","READ");
  TFile* emmFile = new TFile("MEMu2M_run2_datacard.root","READ");

  auto mergedEleSR = MergedEleSR(emmFile,"mergedEMu2M");
  auto eeeSR = MergedEleSR(eeeFile,"mergedEle3E");

  mergedEleSR.add(eeeSR);

  auto mergedStack = mergedEleSR.returnStack();
  auto mergedRatio = mergedEleSR.returnRatio();
  auto mergedError = mergedEleSR.returnError();
  auto mergedRatioError = mergedEleSR.returnRatioError();

  p1->cd();
  mergedEleSR.data_->SetMaximum(1.5*std::max(mergedEleSR.data_->GetMaximum(),mergedEleSR.returnNominal()->GetMaximum()));
  mergedEleSR.data_->SetMinimum(0.001);
  mergedEleSR.data_->GetYaxis()->SetTitle("Events");
  mergedEleSR.data_->GetYaxis()->SetTitleSize(0.08);
  mergedEleSR.data_->GetYaxis()->SetTitleOffset(0.7);
  mergedEleSR.data_->Draw("E1");
  mergedStack->Draw("hist&same");
  mergedEleSR.data_->Draw("E1&same");
  mergedEleSR.sig_->Draw("hist&same");
  mergedError->Draw("2");

  TLegend* mergedLegend = new TLegend(0.43,0.55,0.97,0.93);
  mergedLegend->SetBorderSize(0);
  mergedLegend->SetFillStyle(0);
  mergedLegend->SetMargin(0.2);
  mergedLegend->AddEntry(mergedEleSR.data_,"Data");
  mergedLegend->AddEntry(mergedEleSR.nominals_.at(0),"Prompt bkg");
  mergedLegend->AddEntry(mergedEleSR.nominals_.at(1),"Nonprompt bkg");
  mergedLegend->AddEntry(mergedEleSR.sig_,"M_{X} = 750 GeV, M_{Y} = 1 GeV");
  mergedLegend->SetTextSize(0.055);

  mergedLegend->Draw();

  p2->cd();

  mergedRatio->GetXaxis()->SetTitleSize(0.15);
  mergedRatio->GetXaxis()->SetTitleOffset(0.85);
  mergedRatio->Draw("E1");
  mergedRatioError->Draw("2");

  p1->cd();

  //SaveAs(canvas_2,"merged.pdf",p1);

  //return;

  //p1->SetRightMargin( 0.04 );
  //p2->SetRightMargin( 0.04 );

  TFile* ceeFile = new TFile("MM2E_run2_datacard.root","READ");
  TFile* cmmFile = new TFile("MMFF_run2_datacard.root","READ");

  auto cleanedMuSR = CleanedMuSR(cmmFile,"mergedMu3M");
  auto ceeSR = CleanedMuSR(ceeFile,"mergedMu2E");

  cleanedMuSR.add(ceeSR);

  auto cleanedStack = cleanedMuSR.returnStack();
  auto cleanedRatio = cleanedMuSR.returnRatio();
  auto cleanedError = cleanedMuSR.returnError();
  auto cleanedRatioError = cleanedMuSR.returnRatioError();

  p1->cd();
  cleanedMuSR.data_->SetMaximum(1.5*std::max(cleanedMuSR.data_->GetMaximum(),cleanedMuSR.returnNominal()->GetMaximum()));
  cleanedMuSR.data_->SetMinimum(0.001);
  cleanedMuSR.data_->GetYaxis()->SetTitle("Events");
  cleanedMuSR.data_->GetYaxis()->SetTitleSize(0.08);
  cleanedMuSR.data_->GetYaxis()->SetTitleOffset(0.7);
  cleanedMuSR.data_->Draw("E1");
  cleanedStack->Draw("hist&same");
  cleanedMuSR.data_->Draw("E1&same");
  cleanedMuSR.sig_->Draw("hist&same");
  cleanedError->Draw("2");

  TLegend* cleanedLegend = new TLegend(0.4,0.6,0.93,0.93);
  cleanedLegend->SetBorderSize(0);
  cleanedLegend->AddEntry(cleanedMuSR.data_,"Data");
  cleanedLegend->AddEntry(cleanedMuSR.nominals_.at(0),"Bkg");
  cleanedLegend->AddEntry(cleanedMuSR.sig_,"M_{X} = 2 TeV, M_{Y} = 1 GeV");
  cleanedLegend->SetTextSize(0.055);

  cleanedLegend->Draw();

  p2->cd();
  //cleanedRatio->GetYaxis()->SetTitle("");
  cleanedRatio->GetXaxis()->SetTitleOffset(0.85);
  cleanedRatio->Draw("E1");
  cleanedRatioError->Draw("2");

  p1->cd();

  SaveAs(canvas_2,"cleaned.pdf",p1);
}
