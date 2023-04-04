#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void estimateMEFF(TString era) {
  setTDRStyle();
  // gStyle->SetLineWidth(2);

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"

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
  } else {
    std::cout << "check era..." << std::endl;
  }

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 11;

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

  // EB
  TFile* datafile = new TFile("MergedEleCR_"+era+"_data.root","READ");

  TH1D* SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_SSll_Et_noCorr")->Clone();
  TH1D* SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_SSll_Et_noCorr")->Clone();
  int nbinsSS = 11;
  double xbinsSS[12] = {0, 50, 60, 70, 80, 100, 120, 150, 200, 300, 500, 1000};
  double xcenSS[11] = {55,65,75,90,110,135,175,250,400,750,1000};
  TH1D* SSnum_rebin = (TH1D*)SSnum->Rebin(nbinsSS, "2E_mixedME_SSll_Et_noCorr_rebin", xbinsSS);
  TH1D* SSdenom_rebin = (TH1D*)SSdenom->Rebin(nbinsSS, "2E_antiME_SSll_Et_noCorr_rebin", xbinsSS);

  SSnum_rebin->Divide( SSdenom_rebin );

  TH1D* OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_OSll_Et_noCorr")->Clone();
  TH1D* OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_OSll_Et_noCorr")->Clone();
  int nbinsOS = 26;
  double xbinsOS[27] = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                        225, 250, 275, 300, 350, 400, 500, 600, 700, 1000};
  double xcenOS[26] = {55,65,75,85,95,105,115,125,135,145,155,165,175,185,195,212.5,237.5,262.5,287.5,325,375,450,550,650,850,1000};
  TH1D* OSnum_rebin = (TH1D*)OSnum->Rebin(nbinsOS, "2E_mixedME_OSll_Et_noCorr_rebin", xbinsOS);
  TH1D* OSdenom_rebin = (TH1D*)OSdenom->Rebin(nbinsOS, "2E_antiME_OSll_Et_noCorr_rebin", xbinsOS);

  OSnum_rebin->Divide( OSdenom_rebin );

  TF1* ssboth = new TF1("ssboth","[0]*x+[1]",50,1000);
  ssboth->SetLineColor(kBlue);
  ssboth->SetLineWidth(2);
  ssboth->SetLineStyle(2);
  TFitResultPtr fitSS = SSnum_rebin->Fit(ssboth,"RS");
  fitSS->SetName("fitSS");
  double ciSS[11];
  fitSS->GetConfidenceIntervals(11,1,0,xcenSS,ciSS,0.6827,false);
  double xbinwSS[11] = {5,5,5,10,10,15,25,50,100,250,0};
  double ybinSS[11];

  for (unsigned idx = 0; idx<11; idx++) {
    ybinSS[idx] = ssboth->Eval(xcenSS[idx]);
  }

  auto errSS = new TGraphErrors(11,xcenSS,ybinSS,xbinwSS,ciSS);
  errSS->SetFillColor(kBlue);
  errSS->SetFillStyle(3003);

  TF1* osboth = new TF1("osboth","[0]*x+[1]",50,1000);
  osboth->SetLineColor(kRed);
  osboth->SetLineWidth(2);
  osboth->SetLineStyle(2);
  TFitResultPtr fitOS = OSnum_rebin->Fit(osboth,"RS");
  fitOS->SetName("fitOS");
  double ciOS[26];
  fitOS->GetConfidenceIntervals(26,1,0,xcenOS,ciOS,0.6827,false);
  double xbinwOS[26] = {5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,12.5,12.5,12.5,12.5,25,25,50,50,50,150,0};
  double ybinOS[26];

  for (unsigned idx = 0; idx<26; idx++) {
    ybinOS[idx] = osboth->Eval(xcenOS[idx]);
  }

  auto errOS = new TGraphErrors(26,xcenOS,ybinOS,xbinwOS,ciOS);
  errOS->SetFillColor(kRed);
  errOS->SetFillStyle(3001);

  // save file
  TFile* outfile = new TFile("MEFF_"+era+".root","RECREATE");
  SSnum_rebin->Write();
  OSnum_rebin->Write();
  ssboth->Write();
  osboth->Write();
  fitSS->Write();
  fitOS->Write();
  outfile->Close();

  auto* canvas = new TCanvas("canvas","canvas",50,50,W,H);

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

  // EB
  auto legend = std::make_unique<TLegend>(0.82,0.78,0.95,0.9);
  legend->SetBorderSize(0);
  legend->AddEntry(SSnum_rebin,"SS");
  legend->AddEntry(OSnum_rebin,"OS");

  OSnum_rebin->GetYaxis()->SetRangeUser(0.,0.7);
  OSnum_rebin->SetLineWidth(2);
  OSnum_rebin->GetYaxis()->SetTitle("Fake factor");
  OSnum_rebin->GetXaxis()->SetTitle("E_{T} [GeV]");
  OSnum_rebin->SetStats(0);
  OSnum_rebin->SetLineColor(kRed);
  OSnum_rebin->Draw("E1");

  SSnum_rebin->SetLineColor(kBlue);
  SSnum_rebin->SetLineWidth(2);
  SSnum_rebin->Draw("same&E1");
  errSS->Draw("3");
  errOS->Draw("3");
  legend->Draw();

  TPaveText* textlow = new TPaveText(0.12,0.63,0.55,0.7,"NDC");
  textlow->SetBorderSize(0);
  textlow->SetFillStyle(3025);
  textlow->SetFillColor(0);
  // x < 100 ? [0]*x+[1] : (x < 200 ? 100*[0]+[1]+[2]*(x-100) : 100*[0]+[1]+100*[2]+[3]*(x-200))
  TString textsslow;
  textsslow.Form("(%.3g#pm%.3g) #times E_{T} + %.3f#pm%.3f", ssboth->GetParameter(0), ssboth->GetParError(0),ssboth->GetParameter(1), ssboth->GetParError(1));
  textlow->AddText(textsslow);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(12);
  //TString textsshigh;
  //textsshigh.Form("(%.3g#pm%.3g) #times E_{T} + %.3f#pm%.3f (< 100 GeV)", ssboth->GetParameter(2), ssboth->GetParError(2), ssboth->GetParameter(3), ssboth->GetParError(3));
  //textlow->AddText(textsshigh);
  //((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  //((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(32);

  TPaveText* texthigh = new TPaveText(0.55,0.63,0.95,0.7,"NDC");
  texthigh->SetBorderSize(0);
  texthigh->SetFillColor(0);
  texthigh->SetFillStyle(3025);
  TString textoslow;
  textoslow.Form("(%.3g#pm%.3g) #times E_{T} + %.3f#pm%.3f", osboth->GetParameter(0), osboth->GetParError(0),osboth->GetParameter(1), osboth->GetParError(1));
  texthigh->AddText(textoslow);
  ((TText*)texthigh->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)texthigh->GetListOfLines()->Last())->SetTextAlign(12);

  textlow->Draw();
  texthigh->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("FF_2E.png");

  return;
}
