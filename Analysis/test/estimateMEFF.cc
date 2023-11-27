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
  int H = 800;

  int H_ref = 800;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  // EB
  TFile* datafile = new TFile("EleAnalyzer_"+era+"_data.root","READ");

  auto estimateCenter = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    for (unsigned idx = 1; idx < vec.size()-1; idx++)
      out.push_back( (vec.at(idx) + vec.at(idx+1) ) / 2. );

    out.push_back(vec.back());

    return std::move(out);
  };

  auto estimateWidth = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    for (unsigned idx = 1; idx < vec.size()-1; idx++)
      out.push_back( ( vec.at(idx+1) - vec.at(idx) ) / 2. );

    out.push_back(0.);

    return std::move(out);
  };

  TH1D* SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_EB_mixedME")->Clone();
  TH1D* SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_EB_antiME")->Clone();
  std::vector<double> xbinsSS = {0, 50, 60, 70, 100, 150, 250, 500, 1000};
  const int nbinsSS = xbinsSS.size()-1;
  std::vector<double> xcenSS = estimateCenter(xbinsSS);
  TH1D* SSnum_rebin = (TH1D*)SSnum->Rebin(nbinsSS, "2E_Et_SSCR_EB_mixedME", &(xbinsSS[0]));
  TH1D* SSdenom_rebin = (TH1D*)SSdenom->Rebin(nbinsSS, "2E_Et_SSCR_EB_antiME", &(xbinsSS[0]));

  SSnum_rebin->Divide( SSdenom_rebin );

  TH1D* OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_EB_mixedME")->Clone();
  TH1D* OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_EB_antiME")->Clone();
  std::vector<double> xbinsOS = {0, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200,
                                 225, 250, 275, 300, 400, 500, 700, 1000};
  const int nbinsOS = xbinsOS.size()-1;
  std::vector<double> xcenOS = estimateCenter(xbinsOS);
  TH1D* OSnum_rebin = (TH1D*)OSnum->Rebin(nbinsOS, "2E_Et_OSCR_EB_mixedME", &(xbinsOS[0]));
  TH1D* OSdenom_rebin = (TH1D*)OSdenom->Rebin(nbinsOS, "2E_Et_OSCR_EB_antiME", &(xbinsOS[0]));

  OSnum_rebin->Divide( OSdenom_rebin );

  TF1* ssboth = new TF1("ssboth","[0]",50,1000);
  ssboth->SetLineColor(kBlue);
  ssboth->SetLineWidth(2);
  ssboth->SetLineStyle(2);
  TFitResultPtr fitSS = SSnum_rebin->Fit(ssboth,"RS");
  fitSS->SetName("fitSS");
  double ciSS[nbinsSS];
  fitSS->GetConfidenceIntervals(nbinsSS,1,0,&(xcenSS[0]),ciSS,0.95,false); // 0.6827
  std::vector<double> xbinwSS = estimateWidth(xbinsSS);
  double ybinSS[nbinsSS];

  for (unsigned idx = 0; idx < nbinsSS; idx++) {
    ybinSS[idx] = ssboth->Eval(xcenSS[idx]);
  }

  auto errSS = new TGraphErrors(nbinsSS,&(xcenSS[0]),ybinSS,&(xbinwSS[0]),ciSS);
  errSS->SetFillColor(kBlue);
  errSS->SetFillStyle(3003);

  TF1* osboth = new TF1("osboth","[0]+[1]*x+[2]/sqrt(x)",50,1000);
  osboth->SetLineColor(kRed);
  osboth->SetLineWidth(2);
  osboth->SetLineStyle(2);
  TFitResultPtr fitOS = OSnum_rebin->Fit(osboth,"RS");
  fitOS->SetName("fitOS");
  double ciOS[nbinsOS];
  fitOS->GetConfidenceIntervals(nbinsOS,1,0,&(xcenOS[0]),ciOS,0.95,false);
  std::vector<double> xbinwOS = estimateWidth(xbinsOS);
  double ybinOS[nbinsOS];

  for (unsigned idx = 0; idx < nbinsOS; idx++) {
    ybinOS[idx] = osboth->Eval(xcenOS[idx]);
  }

  auto errOS = new TGraphErrors(nbinsOS,&(xcenOS[0]),ybinOS,&(xbinwOS[0]),ciOS);
  errOS->SetFillColor(kRed);
  errOS->SetFillStyle(3001);

  TH1D* OSnumEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_eta_OSCR_EB_mixedME")->Clone();
  TH1D* OSdenomEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_eta_OSCR_EB_antiME")->Clone();

  OSnumEta->Divide( OSdenomEta );
  TF1* osEta = new TF1("osEta","[0]",-1.5,1.5);
  OSnumEta->Fit(osEta,"RS");

  TH2D* OSnum2d = (TH2D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_eta_OSCR_EB_mixedME")->Clone();
  TH2D* OSdenom2d = (TH2D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_eta_OSCR_EB_antiME")->Clone();

  OSdenom2d->RebinX(2);
  OSnum2d->RebinX(2);

  OSnum2d->Divide( OSdenom2d );

  TObjArray aSlices;
  TF1* os2d = new TF1("os2d","[0]+[1]*x+[2]/sqrt(x)",50,1000);
  os2d->SetParLimits(0,-0.5,0.5);
  os2d->SetParLimits(1,-0.001,0.001);
  OSnum2d->FitSlicesY(os2d,0,-1,0,"QNR",&aSlices);

  TF1* par2 = new TF1("par2","[0]",-1.5,1.5);
  ((TH1D*)aSlices.At(2))->Fit(par2,"RS");
  OSnumEta->Scale(par2->GetParameter(0)/osEta->GetParameter(0)); // par2->GetParameter(0)

  // save file
  TFile* outfile = new TFile("MEFF_"+era+".root","RECREATE");
  SSnum_rebin->Write();
  OSnum_rebin->Write();
  OSnumEta->Write();
  ssboth->Write();
  osboth->Write();
  fitSS->Write();
  fitOS->Write();
  aSlices.At(2)->Write();
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
  auto legend = std::make_unique<TLegend>(0.85,0.8,0.95,0.9);
  legend->SetBorderSize(0);
  legend->AddEntry(SSnum_rebin,"SS");
  legend->AddEntry(OSnum_rebin,"OS");

  OSnum_rebin->GetYaxis()->SetRangeUser(0.,1.0);
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

  TPaveText* textlow = new TPaveText(0.12,0.65,0.3,0.69,"NDC");
  textlow->SetBorderSize(0);
  textlow->SetFillStyle(3025);
  textlow->SetFillColor(0);
  TString textsslow;
  textsslow.Form(" %.3f#pm%.3f", ssboth->GetParameter(0), ssboth->GetParError(0));
  textlow->AddText(textsslow);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(12);
  //TString textsshigh;
  //textsshigh.Form("(%.3g#pm%.3g) #times E_{T} + %.3f#pm%.3f (< 100 GeV)", ssboth->GetParameter(2), ssboth->GetParError(2), ssboth->GetParameter(3), ssboth->GetParError(3));
  //textlow->AddText(textsshigh);
  //((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  //((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(32);

  TPaveText* texthigh = new TPaveText(0.3,0.645,0.95,0.685,"NDC");
  texthigh->SetBorderSize(0);
  texthigh->SetFillColor(0);
  texthigh->SetFillStyle(3025);
  TString textoslow;
  textoslow.Form("%.3f#pm%.3f + (%.3g#pm%.3g)#timesE_{T} + (%.3g#pm%.3g)/#surd E_{T}", osboth->GetParameter(0), osboth->GetParError(0),osboth->GetParameter(1), osboth->GetParError(1), osboth->GetParameter(2), osboth->GetParError(2));
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

  OSnumEta->GetYaxis()->SetRangeUser(0.,1.0);
  OSnumEta->SetLineWidth(2);
  OSnumEta->GetYaxis()->SetTitle("Fake factor");
  OSnumEta->GetXaxis()->SetTitle("#eta_{SC}");
  OSnumEta->SetLineColor(kRed);
  OSnumEta->Draw("E1");

  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("FF_eta.png");

  OSnum2d->Draw("colz");
  canvas->SaveAs("FF_2d.png");

  TF1* par0 = new TF1("par0","[0]",-1.5,1.5);
  ((TH1D*)aSlices.At(0))->Fit(par0,"RS");
  aSlices.At(0)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("FF_slice0.png");

  TF1* par1 = new TF1("par1","[0]",-1.5,1.5);
  ((TH1D*)aSlices.At(1))->Fit(par1,"RS");
  aSlices.At(1)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("FF_slice1.png");

  aSlices.At(2)->Draw("E1");
  OSnumEta->Draw("E1&same");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("FF_slice2.png");

  aSlices.At(3)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("FF_slice3.png");

  return;
}
