#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void estimateREFF(TString era) {
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

  TFile* datafile = new TFile("EleAnalyzer_"+era+"_data.root","READ");

  auto estimateCenter = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    for (unsigned idx = 0; idx < vec.size()-1; idx++)
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

  TH1D* EB_3P0F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_Et_EB")->Clone();
  TH1D* EB_2P1F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_Et_EB")->Clone();
  std::vector<double> xbinsEB = {0,10,20,35,50,75,100,150,200,300,500,1000};
  const int nbinsEB = xbinsEB.size()-1;
  std::vector<double> xcenEB = estimateCenter(xbinsEB);
  TH1D* EB_3P0F_rebin = (TH1D*)EB_3P0F->Rebin(nbinsEB, "EB_3P0F_rebin", &(xbinsEB[0]));
  TH1D* EB_2P1F_rebin = (TH1D*)EB_2P1F->Rebin(nbinsEB, "EB_2P1F_rebin", &(xbinsEB[0]));

  TH1D* all_3P0F_rebin = (TH1D*)EB_3P0F_rebin->Clone();
  TH1D* all_2P1F_rebin = (TH1D*)EB_2P1F_rebin->Clone();

  EB_3P0F_rebin->Divide( EB_2P1F_rebin );

  TH1D* EE_3P0F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_Et_EE")->Clone();
  TH1D* EE_2P1F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_Et_EE")->Clone();
  std::vector<double> xbinsEE = {0,10,20,35,50,75,100,150,200,300,500,1000};
  const int nbinsEE = xbinsEE.size()-1;
  std::vector<double> xcenEE = estimateCenter(xbinsEE);
  TH1D* EE_3P0F_rebin = (TH1D*)EE_3P0F->Rebin(nbinsEE, "EE_3P0F_rebin", &(xbinsEE[0]));
  TH1D* EE_2P1F_rebin = (TH1D*)EE_2P1F->Rebin(nbinsEE, "EE_2P1F_rebin", &(xbinsEE[0]));

  all_3P0F_rebin->Add(EE_3P0F_rebin);
  all_2P1F_rebin->Add(EE_2P1F_rebin);
  all_3P0F_rebin->Divide(all_2P1F_rebin);

  EE_3P0F_rebin->Divide( EE_2P1F_rebin );

  TF1* fitEB = new TF1("REFF_EB","[0]",0,1000);
  TF1* fitEE = new TF1("REFF_EE","[0]",0,1000);
  fitEB->SetLineColor(kRed);
  fitEB->SetLineWidth(2);
  fitEB->SetLineStyle(2);
  fitEE->SetLineColor(kBlue);
  fitEE->SetLineWidth(2);
  fitEE->SetLineStyle(2);
  TFitResultPtr fitResultEB = EB_3P0F_rebin->Fit(fitEB,"RS");
  TFitResultPtr fitResultEE = EE_3P0F_rebin->Fit(fitEE,"RS");

  double ciEB[nbinsEB+1];
  fitResultEB->GetConfidenceIntervals(nbinsEB+1,1,0,&(xbinsEB[0]),ciEB,0.95,false); // 0.6827
  std::vector<double> xbinwEB = estimateWidth(xbinsEB);
  double ybinEB[nbinsEB+1];

  for (unsigned idx = 0; idx < nbinsEB+1; idx++) {
    ybinEB[idx] = fitEB->Eval(xcenEB[idx]);
  }

  auto errEB = new TGraphErrors(nbinsEB+1,&(xcenEB[0]),ybinEB,&(xbinwEB[0]),ciEB);
  errEB->SetFillColor(kRed);
  errEB->SetFillStyle(3004);

  double ciEE[nbinsEE+1];
  fitResultEE->GetConfidenceIntervals(nbinsEE+1,1,0,&(xbinsEE[0]),ciEE,0.95,false); // 0.6827
  std::vector<double> xbinwEE = estimateWidth(xbinsEE);
  double ybinEE[nbinsEE+1];

  for (unsigned idx = 0; idx < nbinsEE+1; idx++) {
    ybinEE[idx] = fitEE->Eval(xcenEE[idx]);
  }

  auto errEE = new TGraphErrors(nbinsEE+1,&(xcenEE[0]),ybinEE,&(xbinwEE[0]),ciEE);
  errEE->SetFillColor(kBlue);
  errEE->SetFillStyle(3005);

  TF1* fitAll = new TF1("REFF_all","[0]",0,1000);
  fitAll->SetLineColor(kRed);
  fitAll->SetLineWidth(2);
  fitAll->SetLineStyle(2);
  all_3P0F_rebin->Fit(fitAll,"R");

  // try dR

  TH1D* dr_3P0F_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_dr_EB")->Clone();
  TH1D* dr_2P1F_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_dr_EB")->Clone();
  std::vector<double> xbinsdrEB = {0.0,0.1,0.2,0.25,0.3,0.35,0.4,0.5,0.7,1.0,1.5,2.0,3.15,6.4};
  const int nbinsdrEB = xbinsdrEB.size()-1;
  std::vector<double> xcendrEB = estimateCenter(xbinsdrEB);
  TH1D* dr_3P0F_EB_rebin = (TH1D*)dr_3P0F_EB->Rebin(nbinsdrEB,"dr_3P0F_EB_rebin",&(xbinsdrEB[0]));
  TH1D* dr_2P1F_EB_rebin = (TH1D*)dr_2P1F_EB->Rebin(nbinsdrEB,"dr_2P1F_EB_rebin",&(xbinsdrEB[0]));

  TH1D* all_dr_3P0F_rebin = (TH1D*)dr_3P0F_EB_rebin->Clone("all_dr_3P0F_rebin");
  TH1D* all_dr_2P1F_rebin = (TH1D*)dr_2P1F_EB_rebin->Clone();

  dr_3P0F_EB_rebin->Divide( dr_2P1F_EB_rebin );

  TH1D* dr_3P0F_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_dr_EE")->Clone();
  TH1D* dr_2P1F_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_dr_EE")->Clone();
  TH1D* dr_3P0F_EE_rebin = (TH1D*)dr_3P0F_EE->Rebin(nbinsdrEB,"dr_3P0F_EE_rebin",&(xbinsdrEB[0]));
  TH1D* dr_2P1F_EE_rebin = (TH1D*)dr_2P1F_EE->Rebin(nbinsdrEB,"dr_2P1F_EE_rebin",&(xbinsdrEB[0]));

  all_dr_3P0F_rebin->Add(dr_3P0F_EE_rebin);
  all_dr_2P1F_rebin->Add(dr_2P1F_EE_rebin);
  all_dr_3P0F_rebin->Divide(all_dr_2P1F_rebin);

  dr_3P0F_EE_rebin->Divide( dr_2P1F_EE_rebin );

  TF1* fitDREB = new TF1("REFF_dr_EB","x > 0.4 ? [0]+[1]*x : [0]+0.4*[1]+[2]*(x-0.4)",0,6.4);
  TF1* fitDREE = new TF1("REFF_dr_EE","x > 0.4 ? [0]+[1]*x : [0]+0.4*[1]+[2]*(x-0.4)",0,6.4);
  fitDREB->SetLineColor(kRed);
  fitDREB->SetLineWidth(2);
  fitDREB->SetLineStyle(2);
  fitDREE->SetLineColor(kBlue);
  fitDREE->SetLineWidth(2);
  fitDREE->SetLineStyle(2);
  TFitResultPtr fitResultDREB = dr_3P0F_EB_rebin->Fit(fitDREB,"RS");
  TFitResultPtr fitResultDREE = dr_3P0F_EE_rebin->Fit(fitDREE,"RS");

  double cidrEB[nbinsdrEB+1];
  fitResultDREB->GetConfidenceIntervals(nbinsdrEB+1,1,0,&(xbinsdrEB[0]),cidrEB,0.95,false);
  std::vector<double> xbinwdrEB = estimateWidth(xbinsdrEB);
  double ybindrEB[nbinsdrEB+1];

  for (unsigned idx = 0; idx < nbinsdrEB+1; idx++) {
    ybindrEB[idx] = fitDREB->Eval(xcendrEB[idx]);
  }

  auto errDREB = new TGraphErrors(nbinsdrEB+1,&(xcendrEB[0]),ybindrEB,&(xbinwdrEB[0]),cidrEB);
  errDREB->SetFillColor(kRed);
  errDREB->SetFillStyle(3004);

  double cidrEE[nbinsdrEB+1];
  fitResultDREE->GetConfidenceIntervals(nbinsdrEB+1,1,0,&(xbinsdrEB[0]),cidrEE,0.95,false);
  double ybindrEE[nbinsdrEB+1];

  for (unsigned idx = 0; idx < nbinsdrEB+1; idx++) {
    ybindrEE[idx] = fitDREE->Eval(xcendrEB[idx]);
  }

  auto errDREE = new TGraphErrors(nbinsdrEB+1,&(xcendrEB[0]),ybindrEE,&(xbinwdrEB[0]),cidrEE);
  errDREE->SetFillColor(kBlue);
  errDREE->SetFillStyle(3005);

  TF1* fitDRall = new TF1("REFF_dr_all","x > 0.4 ? [0]+[1]*x : [0]+0.4*[1]+[2]*(x-0.4)",0,6.4);
  fitDRall->SetLineColor(kRed);
  fitDRall->SetLineWidth(2);
  fitDRall->SetLineStyle(2);
  TFitResultPtr fitResultDR = all_dr_3P0F_rebin->Fit(fitDRall,"RS");

  double cidr[nbinsdrEB+1];
  fitResultDR->GetConfidenceIntervals(nbinsdrEB+1,1,0,&(xbinsdrEB[0]),cidr,0.95,false);
  double ybindr[nbinsdrEB+1];

  for (unsigned idx = 0; idx < nbinsdrEB+1; idx++) {
    ybindr[idx] = fitDRall->Eval(xcendrEB[idx]);
  }

  auto errDR = new TGraphErrors(nbinsdrEB+1,&(xcendrEB[0]),ybindr,&(xbinwdrEB[0]),cidr);
  errDR->SetFillColor(kRed);
  errDR->SetFillStyle(3004);

  // save file
  TFile* outfile = new TFile("REFF_"+era+".root","RECREATE");
  EB_3P0F_rebin->Write();
  EE_3P0F_rebin->Write();
  all_3P0F_rebin->Write();
  dr_3P0F_EB_rebin->Write();
  dr_3P0F_EE_rebin->Write();
  all_dr_3P0F_rebin->Write();
  fitEB->Write();
  fitEE->Write();
  fitAll->Write();
  fitDREB->Write();
  fitDREE->Write();
  fitDRall->Write();
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

  EB_3P0F_rebin->GetYaxis()->SetRangeUser(0.,0.3);
  EB_3P0F_rebin->SetLineWidth(2);
  EB_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  EB_3P0F_rebin->GetXaxis()->SetTitle("E_{T} [GeV]");
  EB_3P0F_rebin->SetStats(0);
  EB_3P0F_rebin->SetLineColor(kRed);
  EB_3P0F_rebin->Draw("E1");

  EE_3P0F_rebin->SetLineColor(kBlue);
  EE_3P0F_rebin->SetLineWidth(2);
  EE_3P0F_rebin->Draw("same&E1");

  errEB->Draw("3");
  errEE->Draw("3");

  TPaveText* textEB = new TPaveText(0.15,0.6,0.5,0.75,"NDC");
  textEB->SetBorderSize(0);
  textEB->SetFillColor(0);
  TString textEBstr;
  textEBstr.Form("%.3g #pm %.3g (EB)", fitEB->GetParameter(0), fitEB->GetParError(0));
  textEB->AddText(textEBstr);
  ((TText*)textEB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textEB->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* textEE = new TPaveText(0.6,0.6,0.92,0.75,"NDC");
  textEE->SetBorderSize(0);
  textEE->SetFillColor(0);
  TString textEEstr;
  textEEstr.Form("%.3g #pm %.3g (EE)", fitEE->GetParameter(0), fitEE->GetParError(0));
  textEE->AddText(textEEstr);
  ((TText*)textEE->GetListOfLines()->Last())->SetTextColor(kBlue);

  textEB->Draw();
  textEE->Draw();

  TLegend* legend = new TLegend(0.8,0.75,0.95,0.9);
  legend->SetBorderSize(0);
  legend->AddEntry(EB_3P0F_rebin,"EB");
  legend->AddEntry(EE_3P0F_rebin,"EE");
  legend->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("REFF_"+era+".pdf");

  // all
  all_3P0F_rebin->GetYaxis()->SetRangeUser(0.,0.3);
  all_3P0F_rebin->SetLineWidth(2);
  all_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  all_3P0F_rebin->GetXaxis()->SetTitle("E_{T} [GeV]");
  all_3P0F_rebin->SetStats(0);
  all_3P0F_rebin->SetLineColor(kRed);
  all_3P0F_rebin->Draw("E1");

  TPaveText* textAll = new TPaveText(0.15,0.6,0.4,0.75,"NDC");
  textAll->SetBorderSize(0);
  textAll->SetFillColor(0);
  TString textAllStr;
  textAllStr.Form("%.3g #pm %.3g", fitAll->GetParameter(0), fitAll->GetParError(0));
  textAll->AddText(textAllStr);
  ((TText*)textAll->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textAll->GetListOfLines()->Last())->SetTextAlign(12);

  textAll->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("REFF_all_"+era+".png");

  // draw dr

  dr_3P0F_EB_rebin->Draw("E1");

  dr_3P0F_EB_rebin->GetYaxis()->SetRangeUser(0.,2.0);
  dr_3P0F_EB_rebin->SetLineWidth(2);
  dr_3P0F_EB_rebin->GetYaxis()->SetTitle("Fake factor");
  dr_3P0F_EB_rebin->GetXaxis()->SetTitle("#Delta R");
  dr_3P0F_EB_rebin->SetStats(0);
  dr_3P0F_EB_rebin->SetLineColor(kRed);
  dr_3P0F_EB_rebin->Draw("E1");

  dr_3P0F_EE_rebin->SetLineColor(kBlue);
  dr_3P0F_EE_rebin->SetLineWidth(2);
  dr_3P0F_EE_rebin->Draw("same&E1");

  errDREE->Draw("3");
  errDREB->Draw("3");

  TPaveText* textDREB = new TPaveText(0.15,0.6,0.5,0.75,"NDC");
  textDREB->SetBorderSize(0);
  textDREB->SetFillColor(0);
  TString textDREBstrBelow;
  textDREBstrBelow.Form("%.3g#DeltaR + %.3g (EB) (#Delta R < 0.4)", fitDREB->GetParameter(2), fitDREB->GetParameter(0)+0.4*fitDREB->GetParameter(1)-0.4*fitDREB->GetParameter(2));
  TString textDREBstrAbove;
  textDREBstrAbove.Form("%.3g#DeltaR + %.3g (EB) (#Delta R > 0.4)", fitDREB->GetParameter(1), fitDREB->GetParameter(0));
  textDREB->AddText(textDREBstrBelow);
  ((TText*)textDREB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDREB->GetListOfLines()->Last())->SetTextAlign(12);
  textDREB->AddText(textDREBstrAbove);
  ((TText*)textDREB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDREB->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* textDREE = new TPaveText(0.55,0.6,0.92,0.75,"NDC");
  textDREE->SetBorderSize(0);
  textDREE->SetFillColor(0);
  TString textDREEstrBelow;
  textDREEstrBelow.Form("%.3g#DeltaR + %.3g (EE) (#Delta R < 0.4)", fitDREE->GetParameter(2), fitDREE->GetParameter(0)+0.4*fitDREE->GetParameter(1)-0.4*fitDREE->GetParameter(2));
  TString textDREEstrAbove;
  textDREEstrAbove.Form("%.3g#DeltaR + %.3g (EE) (#Delta R > 0.4)", fitDREE->GetParameter(1), fitDREE->GetParameter(0));
  textDREE->AddText(textDREEstrBelow);
  ((TText*)textDREE->GetListOfLines()->Last())->SetTextColor(kBlue);
  textDREE->AddText(textDREEstrAbove);
  ((TText*)textDREE->GetListOfLines()->Last())->SetTextColor(kBlue);

  textDREB->Draw();
  textDREE->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();

  canvas->SaveAs("REFF_dr_"+era+".png");

  // all
  all_dr_3P0F_rebin->Draw("E1");

  all_dr_3P0F_rebin->GetYaxis()->SetRangeUser(0.,2.0);
  all_dr_3P0F_rebin->SetLineWidth(2);
  all_dr_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  all_dr_3P0F_rebin->GetXaxis()->SetTitle("#Delta R");
  all_dr_3P0F_rebin->SetStats(0);
  all_dr_3P0F_rebin->SetLineColor(kRed);
  all_dr_3P0F_rebin->Draw("E1");

  errDR->Draw("3");

  TPaveText* textDRall = new TPaveText(0.5,0.6,0.9,0.75,"NDC");
  textDRall->SetBorderSize(0);
  textDRall->SetFillColor(0);
  TString textDRAllstrBelow;
  textDRAllstrBelow.Form("%.3g #times #DeltaR + %.3g (#Delta R < 0.4)", fitDRall->GetParameter(2), fitDRall->GetParameter(0)+0.4*fitDRall->GetParameter(1)-0.4*fitDRall->GetParameter(2));
  TString textDRallAbove;
  textDRallAbove.Form("%.3g #times #DeltaR + %.3g (#Delta R > 0.4)", fitDRall->GetParameter(1), fitDRall->GetParameter(0));
  textDRall->AddText(textDRAllstrBelow);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextAlign(12);
  textDRall->AddText(textDRallAbove);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextAlign(12);

  textDRall->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();

  canvas->SaveAs("REFF_all_dr_"+era+".pdf");

  return;
}
