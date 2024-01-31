#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void estimateMEFF(TString era) {
  setTDRStyle();
  gStyle->SetOptFit(0);

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
  TFile* datafile = new TFile("MergedEleCR_"+era+"_data.root","READ");

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
  errSS->SetFillStyle(3004);

  TH1D* OSnumEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_eta_OSCR_EB_mixedME")->Clone();
  TH1D* OSdenomEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_eta_OSCR_EB_antiME")->Clone();

  OSnumEta->Divide( OSdenomEta );

  TH1D* SSnumEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_SSll_eta")->Clone();
  TH1D* SSdenomEta = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_SSll_eta")->Clone();
  SSdenomEta->Rebin(SSdenomEta->GetNbinsX()/SSnumEta->GetNbinsX());
  SSdenomEta->Rebin(5);
  SSnumEta->Rebin(5);

  SSnumEta->Divide( SSdenomEta );

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

  TF2* osboth = new TF2("osboth","[0]+[1]*y+[2]*sqrt(cosh(x)/y)",-1.4442,1.4442,50.,1000.);
  osboth->SetLineColor(kRed);
  osboth->SetLineWidth(2);
  osboth->SetLineStyle(2);
  TFitResultPtr fitOS = OSnum2d->Fit(osboth,"RS+");
  fitOS->SetName("fitOS");
  double ciOS[nbinsOS];
  std::vector<double> xycenOS;
  const double projectX = 0.75;
  TString formula_projected;
  formula_projected.Form("[0]+[1]*x+[2]*sqrt(cosh(%.2g)/x)",projectX);
  TF1* osboth_projected = new TF1("osboth_projected",formula_projected,50.,1000.);
  osboth_projected->SetLineColor(kRed);
  osboth_projected->SetLineWidth(2);
  osboth_projected->SetLineStyle(2);
  osboth_projected->SetParameters(osboth->GetParameters());

  for (const auto& xcen : xcenOS) {
    xycenOS.push_back(projectX);
    xycenOS.push_back(xcen);
  }

  fitOS->GetConfidenceIntervals(nbinsOS,2,1,&(xycenOS[0]),ciOS,0.95,false);
  std::vector<double> xbinwOS = estimateWidth(xbinsOS);
  double ybinOS[nbinsOS];

  for (unsigned idx = 0; idx < nbinsOS; idx++) {
    ybinOS[idx] = osboth->Eval(projectX,xcenOS[idx]);
  }

  auto errOS = new TGraphErrors(nbinsOS,&(xcenOS[0]),ybinOS,&(xbinwOS[0]),ciOS);
  errOS->SetFillColor(kRed);
  errOS->SetFillStyle(3005);

  // save file
  TFile* outfile = new TFile("MEFF_"+era+".root","RECREATE");
  SSnum_rebin->Write();
  OSnum2d->Write();
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
  auto legend = std::make_unique<TLegend>(0.8,0.75,0.95,0.9);
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
  osboth_projected->Draw("same");

  SSnum_rebin->SetLineColor(kBlue);
  SSnum_rebin->SetLineWidth(2);
  SSnum_rebin->Draw("same&E1");
  errSS->Draw("3");
  errOS->Draw("3");
  legend->Draw();

  TPaveText* textlow = new TPaveText(0.12,0.65,0.28,0.69,"NDC");
  textlow->SetBorderSize(0);
  textlow->SetFillStyle(3025);
  textlow->SetFillColor(0);
  TString textsslow;
  textsslow.Form(" %.3f#pm%.3f", ssboth->GetParameter(0), ssboth->GetParError(0));
  textlow->AddText(textsslow);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textlow->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* texthigh = new TPaveText(0.28,0.645,0.965,0.685,"NDC");
  texthigh->SetBorderSize(0);
  texthigh->SetFillColor(0);
  texthigh->SetFillStyle(3025);
  TString textoslow;
  textoslow.Form("%.3f#pm%.3f + (%.3g#pm%.3g)#timesp_{T} + (%.3g#pm%.3g)/#surdE", osboth->GetParameter(0), osboth->GetParError(0),osboth->GetParameter(1), osboth->GetParError(1), osboth->GetParameter(2), osboth->GetParError(2));
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
  canvas->SaveAs("MEFF_Et_"+era+".pdf");

  OSnumEta->GetYaxis()->SetRangeUser(0.,1.0);
  OSnumEta->SetLineWidth(2);
  OSnumEta->GetYaxis()->SetTitle("Fake factor");
  OSnumEta->GetXaxis()->SetTitle("#eta_{SC}");
  OSnumEta->SetLineColor(kRed);
  OSnumEta->Draw("E1");
  SSnumEta->SetLineWidth(2);
  SSnumEta->SetLineColor(kBlue);
  SSnumEta->Draw("E1&same");
  OSnumEta->Draw("E1&same");
  legend->Draw();

  double projectY = 65.;
  TString formula_eta;
  formula_eta.Form("[0]+[1]*%.3g+[2]*sqrt(cosh(x)/%.3g)",projectY,projectY);
  TF1* osboth_eta = new TF1("osboth_eta",formula_eta,-1.4442,1.4442);
  osboth_eta->SetLineColor(kRed);
  osboth_eta->SetLineWidth(2);
  osboth_eta->SetLineStyle(2);
  osboth_eta->SetParameters(osboth->GetParameters());
  osboth_eta->Draw("same");

  const int binstart = OSnumEta->FindFixBin(-1.5);
  const int binend = OSnumEta->FindFixBin(1.5);
  const unsigned nbinsEta = binend-binstart+1;
  std::vector<double> etaCenter, etaCenterXY, xbinwEta;
  double ciOSeta[nbinsEta];

  for (int bin=binstart; bin<=binend; bin++) {
    etaCenter.push_back(OSnumEta->GetBinCenter(bin));
    etaCenterXY.push_back(OSnumEta->GetBinCenter(bin));
    etaCenterXY.push_back(projectY);
    xbinwEta.push_back(0.);
  }

  fitOS->GetConfidenceIntervals(nbinsEta,2,1,&(etaCenterXY[0]),ciOSeta,0.95,false);
  double ybinOSeta[nbinsEta];

  for (unsigned idx = 0; idx < nbinsEta; idx++) {
    ybinOSeta[idx] = osboth->Eval(etaCenter.at(idx),projectY);
  }

  auto errOSeta = new TGraphErrors(nbinsEta,&(etaCenter[0]),ybinOSeta,&(xbinwEta[0]),ciOSeta);
  errOSeta->SetFillColor(kRed);
  errOSeta->SetFillStyle(3005);
  errOSeta->Draw("3");

  TF1* ssEta = new TF1("ssEta","[0]",-1.4442,1.4442);
  ssEta->SetLineColor(kBlue);
  ssEta->SetLineWidth(2);
  ssEta->SetLineStyle(2);
  ssEta->SetParameters(ssboth->GetParameters());
  ssEta->Draw("same");

  double etaSS[2] = {-1.4442,1.4442};
  double ybinSSeta[2] = {ssEta->GetParameter(0),ssEta->GetParameter(0)};
  double ciSSeta[2] = {ciSS[0],ciSS[0]};
  auto errSSeta = new TGraphErrors(2,etaSS,ybinSSeta,&(xbinwEta[0]),ciSSeta);
  errSSeta->SetFillColor(kBlue);
  errSSeta->SetFillStyle(3004);
  errSSeta->Draw("3");

  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_eta_"+era+".pdf");

  CMS_lumi( canvas, iPeriod, iPos );

  extraText = "Internal";

  OSnum2d->Draw("col");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->SaveAs("FF_2d.png");

  TF1* par0 = new TF1("par0","[0]",-1.5,1.5);
  ((TH1D*)aSlices.At(0))->Fit(par0,"RS");
  aSlices.At(0)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_slice0_"+era+".pdf");

  TF1* par1 = new TF1("par1","[0]",-1.5,1.5);
  ((TH1D*)aSlices.At(1))->Fit(par1,"RS");
  aSlices.At(1)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_slice1_"+era+".pdf");

  TF1* par2cosh = new TF1("par2cosh","[0]*sqrt(cosh(x))",-1.5,1.5);
  ((TH1D*)aSlices.At(2))->Fit(par2cosh,"RS");
  aSlices.At(2)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_slice2_"+era+".pdf");

  aSlices.At(3)->Draw("E1");
  CMS_lumi( canvas, iPeriod, iPos );
  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("MEFF_slice3_"+era+".pdf");

  return;
}
