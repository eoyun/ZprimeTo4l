#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void estimateMEFF() {
  setTDRStyle();
  // gStyle->SetLineWidth(2);

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "2016 (13 TeV)";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  lumi_13TeV = "16.8 fb^{-1}";

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

  TFile* datafile = new TFile("MergedEleCR_20UL16_data.root","READ");

  TH1D* SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_mixedME")->Clone();
  TH1D* SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_antiME")->Clone();
  int nbinsSS = 12;
  double xbinsSS[13] = {0, 50, 60, 70, 80, 90, 100, 125, 150, 200, 300, 500, 1000};
  TH1D* SSnum_rebin = (TH1D*)SSnum->Rebin(nbinsSS, "2E_Et_SSCR_mixedME_rebin", xbinsSS);
  TH1D* SSdenom_rebin = (TH1D*)SSdenom->Rebin(nbinsSS, "2E_Et_SSCR_antiME_rebin", xbinsSS);

  SSnum_rebin->Divide( SSdenom_rebin );

  TH1D* OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_mixedME")->Clone();
  TH1D* OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_antiME")->Clone();
  int nbinsOS = 25;
  double xbinsOS[26] = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                        225, 250, 275, 300, 350, 400, 500, 700, 1000};
  TH1D* OSnum_rebin = (TH1D*)OSnum->Rebin(nbinsOS, "2E_Et_OSCR_mixedME_rebin", xbinsOS);
  TH1D* OSdenom_rebin = (TH1D*)OSdenom->Rebin(nbinsOS, "2E_Et_OSCR_antiME_rebin", xbinsOS);

  OSnum_rebin->Divide( OSdenom_rebin );

  TF1* sslow = new TF1("sslow","[0]",50,200);
  TF1* sshigh = new TF1("sshigh","[0]",200,1000);
  sslow->SetLineColor(kBlue);
  sslow->SetLineWidth(2);
  sslow->SetLineStyle(2);
  sshigh->SetLineColor(kBlue);
  sshigh->SetLineWidth(2);
  sshigh->SetLineStyle(2);
  SSnum_rebin->Fit(sslow,"R");
  SSnum_rebin->Fit(sshigh,"R+");

  TF1* oshigh = new TF1("oshigh","[0]",500,1000);
  TF1* oslow = new TF1("oslow","[0]*(x-500.)+[1]",50,500);
  oslow->SetLineColor(kRed);
  oslow->SetLineWidth(2);
  oslow->SetLineStyle(2);
  oshigh->SetLineColor(kRed);
  oshigh->SetLineWidth(2);
  oshigh->SetLineStyle(2);
  OSnum_rebin->Fit(oshigh,"R");
  oslow->FixParameter(1,oshigh->GetParameter(0));
  OSnum_rebin->Fit(oslow,"R+");

  TFile* outfile = new TFile("MEFF_20UL16.root","RECREATE");
  SSnum_rebin->Write();
  OSnum_rebin->Write();
  sslow->Write();
  sshigh->Write();
  oslow->Write();
  oshigh->Write();
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

  // auto legend = std::make_unique<TLegend>(0.15,0.62,0.35,0.77);
  auto legend = std::make_unique<TLegend>(0.82,0.78,1.,0.9);
  legend->SetBorderSize(0);
  legend->AddEntry(SSnum_rebin,"SS");
  legend->AddEntry(OSnum_rebin,"OS");

  OSnum_rebin->GetYaxis()->SetRangeUser(0.,0.5);
  OSnum_rebin->SetLineWidth(2);
  OSnum_rebin->GetYaxis()->SetTitle("Fake factor");
  OSnum_rebin->GetXaxis()->SetTitle("E_{T} [GeV]");
  OSnum_rebin->SetStats(0);
  OSnum_rebin->SetLineColor(kRed);
  OSnum_rebin->Draw("E1");

  SSnum_rebin->SetLineColor(kBlue);
  SSnum_rebin->SetLineWidth(2);
  SSnum_rebin->Draw("same&E1");
  legend->Draw();

  TPaveText* textlow = new TPaveText(0.15,0.6,0.65,0.75,"NDC");
  textlow->SetBorderSize(0);
  textlow->SetFillColor(0);
  TString textsslow;
  textsslow.Form("%.3f #pm %.3f (< 200 GeV)", sslow->GetParameter(0), sslow->GetParError(0));
  textlow->AddText(textsslow);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(12);
  TString textoslow;
  textoslow.Form("(%.3g#pm%.3g) #times (E_{T}-500) + %.3f (< 500 GeV)", oslow->GetParameter(0), oslow->GetParError(0), oslow->GetParameter(1));
  textlow->AddText(textoslow);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* texthigh = new TPaveText(0.68,0.6,0.95,0.75,"NDC");
  texthigh->SetBorderSize(0);
  texthigh->SetFillColor(0);
  TString textsshigh;
  textsshigh.Form("%.3f #pm %.3f (> 200 GeV)", sshigh->GetParameter(0), sshigh->GetParError(0));
  texthigh->AddText(textsshigh);
  ((TText*)texthigh->GetListOfLines()->Last())->SetTextColor(kBlue);
  TString textoshigh;
  textoshigh.Form("%.3f #pm %.3f (> 500 GeV)", oshigh->GetParameter(0), oshigh->GetParError(0));
  texthigh->AddText(textoshigh);
  ((TText*)texthigh->GetListOfLines()->Last())->SetTextColor(kRed);

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
