#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  double valLumi = 0.;
  TString postfix = era;

  if (era=="20UL16APV") {
    lumi_13TeV = "19.5 fb^{-1}";
    valLumi = 19.5;
  } else if (era=="20UL16") {
    lumi_13TeV = "16.8 fb^{-1}";
    valLumi = 16.8;
    postfix = "";
  } else if (era=="20UL17") {
    lumi_13TeV = "41.48 fb^{-1}";
    valLumi = 41.48;
  } else if (era=="20UL18") {
    lumi_13TeV = "59.83 fb^{-1}";
    valLumi = 59.83;
  } else {
    std::cout << "check params!" << std::endl;
    throw;
  }

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

  TFile* datafile = new TFile("MergedEleCR_"+era+"_data.root","READ");
  TFile* WZfile = new TFile("MergedEleCR_"+era+"_WZ.root","READ");
  TFile* ZZfile = new TFile("MergedEleCR_"+era+"_ZZ.root","READ");

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

  auto drawRatio = [] (const TH1D* numer, const TH1D* denom, TPad* pad) {
    TH1D* ratio = (TH1D*)numer->Clone();
    ratio->SetStats(0);
    ratio->SetTitle("");
    ratio->Divide(denom);
    ratio->GetYaxis()->SetTitle("Obs/Exp");
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetTitleOffset(0.75);
    ratio->SetLineColor(kBlack);

    pad->cd();
    ratio->Draw("E1");
  };

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

  auto compare = [&datafile] (TString numName, TString denomName, const int color, TPad* pad) -> std::pair<TH1D*,TH1D*> {
    TH1D* num = (TH1D*)datafile->Get(numName)->Clone();
    TH1D* denom = (TH1D*)datafile->Get(denomName)->Clone();

    if ( denom->GetNbinsX() % num->GetNbinsX()!=0 ) {
      // hardcode (rebin to GCD)
      denom->Rebin( 5 );
      num->Rebin( 2 );
    } else {
      denom->Rebin( denom->GetNbinsX()/num->GetNbinsX() );
    }

    denom->SetFillColor(color);
    pad->cd();
    num->SetMaximum( 1.2*std::max(num->GetMaximum(),denom->GetMaximum()) );
    num->SetLineWidth(2);
    num->SetLineColor(kBlack);
    num->Draw("E1");
    denom->SetLineWidth(0);
    denom->Draw("hist&same");
    num->Draw("E1&same");

    return std::make_pair(num,denom);
  };

  TString anlyzrMC = "mergedEleCRanalyzer"+postfix;
  TString anlyzrData = "mergedEleCRanalyzerData";

  TH1D* invM_SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_CRME_SSll_invM")->Clone();
  TH1D* invM_SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_SSll_invM_xFF")->Clone();
  TH1D* invM_SSmixed = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_SSll_invM")->Clone();
  TH1D* invM_SSanti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_SSll_invM_CR_xFF")->Clone();

  TH1D* invM_OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_CRME_OSll_invM")->Clone();
  TH1D* invM_OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_OSll_invM_xFF")->Clone();
  TH1D* invM_OSmixed = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_OSll_invM")->Clone();
  TH1D* invM_OSanti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_OSll_invM_CR_xFF")->Clone();

  // invM
  compare(anlyzrData+"/2E_CRME_SSll_invM",anlyzrData+"/2E_mixedME_SSll_invM_xFF",kCyan+1,canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedSS.png");
  auto apair = compare(anlyzrData+"/2E_mixedME_SSll_invM",anlyzrData+"/2E_antiME_SSll_invM_CR_xFF",kCyan+1,p1);
  drawRatio(apair.first,apair.second,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiSS.png",p1);
  compare(anlyzrData+"/2E_CRME_OSll_invM",anlyzrData+"/2E_mixedME_OSll_invM_xFF",kOrange,canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedOS.png");
  apair = compare(anlyzrData+"/2E_mixedME_OSll_invM",anlyzrData+"/2E_antiME_OSll_invM_CR_xFF",kOrange,p1);
  drawRatio(apair.first,apair.second,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiOS.png",p1);

  // Et
  compare(anlyzrData+"/2E_CRME_SSll_Et",anlyzrData+"/2E_mixedAntiME_SSll_Et_xFF",kCyan+1,canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedSS.png");
  apair = compare(anlyzrData+"/2E_mixedME_SSll_Et_noCorr",anlyzrData+"/2E_antiME_SSll_Et_xFF",kCyan+1,p1);
  drawRatio(apair.first,apair.second,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiSS.png",p1);
  compare(anlyzrData+"/2E_CRME_OSll_Et",anlyzrData+"/2E_mixedAntiME_OSll_Et_xFF",kOrange,canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedOS.png");
  apair = compare(anlyzrData+"/2E_mixedME_OSll_Et_noCorr",anlyzrData+"/2E_antiME_OSll_Et_xFF",kOrange,p1);
  drawRatio(apair.first,apair.second,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiOS.png",p1);

  // eta
  compare(anlyzrData+"/2E_CRME_SSll_eta",anlyzrData+"/2E_mixedAntiME_SSll_eta_xFF",kCyan+1,canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedSS.png");
  apair = compare(anlyzrData+"/2E_mixedME_SSll_eta",anlyzrData+"/2E_antiME_SSll_eta_xFF",kCyan+1,p1);
  drawRatio(apair.first,apair.second,p2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiSS.png",p1);
  compare(anlyzrData+"/2E_CRME_OSll_eta",anlyzrData+"/2E_mixedAntiME_OSll_eta_xFF",kOrange,canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedOS.png");
  apair = compare(anlyzrData+"/2E_mixedME_OSll_eta",anlyzrData+"/2E_antiME_OSll_eta_xFF",kOrange,p1);
  drawRatio(apair.first,apair.second,p2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiOS.png",p1);

  auto square = [](double x) { return x*x; };

  auto subtractHist = [&square](const TH1D* denom, const TH1D* denom_prompt) -> TH1D* {
    auto* cloned = (TH1D*)denom->Clone();

    for (unsigned idx = 0; idx < cloned->GetNbinsX()+2; idx++) {
      cloned->SetBinContent(idx,cloned->GetBinContent(idx)-denom_prompt->GetBinContent(idx));
      cloned->SetBinError(idx,std::sqrt(square(cloned->GetBinError(idx))+square(denom_prompt->GetBinError(idx))));
    }

    return cloned;
  };

  double intlumi = valLumi;
  double WZxsec = 5.213;
  double ZZxsec = 1.325;
  double WZsumwgt = ((TH1D*)WZfile->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1);
  double ZZsumwgt = ((TH1D*)ZZfile->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1);

  auto compare3E = [&](TString numName, TString denomName) {
    TH1D* denomSS = (TH1D*)datafile->Get(anlyzrData+"/"+denomName+"_xSSFF")->Clone();
    TH1D* numer = (TH1D*)datafile->Get(anlyzrData+"/"+numName)->Clone();

    TH1D* denomSS_WZ = (TH1D*)WZfile->Get(anlyzrMC+"/"+denomName+"_xSSFF")->Clone();
    denomSS_WZ->Scale( intlumi*1000.*WZxsec/WZsumwgt );
    TH1D* denomSS_ZZ = (TH1D*)ZZfile->Get(anlyzrMC+"/"+denomName+"_xSSFF")->Clone();
    denomSS_ZZ->Scale( intlumi*1000.*ZZxsec/ZZsumwgt );
    TH1D* denomOS_WZ = (TH1D*)WZfile->Get(anlyzrMC+"/"+denomName+"_xOSFF")->Clone();
    denomOS_WZ->Scale( intlumi*1000.*WZxsec/WZsumwgt );
    TH1D* denomOS_ZZ = (TH1D*)ZZfile->Get(anlyzrMC+"/"+denomName+"_xOSFF")->Clone();
    denomOS_ZZ->Scale( intlumi*1000.*ZZxsec/ZZsumwgt );

    TH1D* denomSSfinal = subtractHist( subtractHist(denomSS,denomSS_WZ), denomSS_ZZ);
    denomSSfinal->SetFillColor(kCyan+1);
    denomSSfinal->Rebin( denomSSfinal->GetNbinsX()/numer->GetNbinsX() );
    TH1D* denomOSfinal = (TH1D*)denomOS_WZ->Clone();
    denomOSfinal->Add(denomOS_ZZ);
    denomOSfinal->SetFillColor(kOrange);
    denomOSfinal->Rebin( denomOSfinal->GetNbinsX()/numer->GetNbinsX() );

    THStack* denomFinal = new THStack(numName,";GeV;");
    denomOSfinal->SetLineWidth(0);
    denomSSfinal->SetLineWidth(0);
    denomFinal->Add(denomOSfinal);
    denomFinal->Add(denomSSfinal);

    numer->SetLineWidth(2);
    numer->SetLineColor(kBlack);
    numer->SetMaximum(1.5*numer->GetMaximum());

    numer->Draw("E1");
    denomFinal->Draw("hist&same");
    numer->Draw("E1&same");

    if (!numName.Contains("Et")) {
      TLegend* legend_left = new TLegend(0.15,0.7,0.4,0.9);
      legend_left->SetBorderSize(0);
      legend_left->AddEntry(numer,"Data");
      legend_left->AddEntry(denomSSfinal,"Data-driven (X #rightarrow ME)");
      legend_left->AddEntry(denomOSfinal,"Data-driven (e #rightarrow ME)");
      legend_left->Draw();
    }
  };

  canvas_2->cd();

  compare3E("3E_CRME_lll_invM","3E_antiME_lll_invM_CR");
  SaveAs(canvas_2,"FF_3E_invM.png");

  compare3E("3E_CRME_all_Et","3E_antiME_Et");
  SaveAs(canvas_2,"FF_3E_Et.png");

  compare3E("3E_CRME_all_eta","3E_antiME_eta");
  SaveAs(canvas_2,"FF_3E_eta.png");

  return;
}
