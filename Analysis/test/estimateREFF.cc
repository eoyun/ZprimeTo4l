#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void estimateREFF() {
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

  TFile* datafile = new TFile("ResolvedEleCR_20UL16_data.root","READ");
  TFile* WZfile = new TFile("ResolvedEleCR_20UL16_WZ.root","READ");
  TFile* ZZfile = new TFile("ResolvedEleCR_20UL16_ZZ.root","READ");

  TH1D* EB_3P0F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_Et_EB")->Clone();
  TH1D* EB_2P1F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_Et_EB")->Clone();
  int nbinsEB = 7;
  double xbinsEB[8] = {0,10,20,35,50,100,200,500};
  TH1D* EB_3P0F_rebin = (TH1D*)EB_3P0F->Rebin(nbinsEB, "EB_3P0F_rebin", xbinsEB);
  TH1D* EB_2P1F_rebin = (TH1D*)EB_2P1F->Rebin(nbinsEB, "EB_2P1F_rebin", xbinsEB);

  EB_3P0F_rebin->Divide( EB_2P1F_rebin );

  TH1D* EE_3P0F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_Et_EE")->Clone();
  TH1D* EE_2P1F = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_Et_EE")->Clone();
  int nbinsEE = 7;
  double xbinsEE[8] = {0,10,20,35,50,100,200,500};
  TH1D* EE_3P0F_rebin = (TH1D*)EE_3P0F->Rebin(nbinsEE, "EE_3P0F_rebin", xbinsEE);
  TH1D* EE_2P1F_rebin = (TH1D*)EE_2P1F->Rebin(nbinsEE, "EE_2P1F_rebin", xbinsEE);

  EE_3P0F_rebin->Divide( EE_2P1F_rebin );

  TF1* fitEB = new TF1("REFF_EB","[0]",0,500);
  TF1* fitEE = new TF1("REFF_EE","[0]",0,500);
  fitEB->SetLineColor(kRed);
  fitEB->SetLineWidth(2);
  fitEB->SetLineStyle(2);
  fitEE->SetLineColor(kBlue);
  fitEE->SetLineWidth(2);
  fitEE->SetLineStyle(2);
  EB_3P0F_rebin->Fit(fitEB,"R");
  EE_3P0F_rebin->Fit(fitEE,"R");

  // try dR

  TH1D* dr_3P0F_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_dr_EB")->Clone();
  TH1D* dr_2P1F_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_dr_EB")->Clone();
  int nbinsdrEB = 12;
  double xbinsdrEB[13] = {0.0,0.1,0.2,0.3,0.5,1.0,1.5,2.0,2.5,3.15,3.5,4.0,6.4};
  TH1D* dr_3P0F_EB_rebin = (TH1D*)dr_3P0F_EB->Rebin(nbinsdrEB,"dr_3P0F_EB_rebin",xbinsdrEB);
  TH1D* dr_2P1F_EB_rebin = (TH1D*)dr_2P1F_EB->Rebin(nbinsdrEB,"dr_2P1F_EB_rebin",xbinsdrEB);

  dr_3P0F_EB_rebin->Divide( dr_2P1F_EB_rebin );

  TH1D* dr_3P0F_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P0F_dr_EE")->Clone();
  TH1D* dr_2P1F_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P1F_dr_EE")->Clone();
  TH1D* dr_3P0F_EE_rebin = (TH1D*)dr_3P0F_EE->Rebin(nbinsdrEB,"dr_3P0F_EE_rebin",xbinsdrEB);
  TH1D* dr_2P1F_EE_rebin = (TH1D*)dr_2P1F_EE->Rebin(nbinsdrEB,"dr_2P1F_EE_rebin",xbinsdrEB);

  dr_3P0F_EE_rebin->Divide( dr_2P1F_EE_rebin );

  TF1* fitDREB_below = new TF1("REFF_dr_below_EB","[0]",0,0.3);
  TF1* fitDREE_below = new TF1("REFF_dr_below_EE","[0]",0,0.3);
  TF1* fitDREB_above = new TF1("REFF_dr_above_EB","[0]",0.3,6.4);
  TF1* fitDREE_above = new TF1("REFF_dr_above_EE","[0]",0.3,6.4);
  fitDREB_below->SetLineColor(kRed);
  fitDREB_below->SetLineWidth(2);
  fitDREB_below->SetLineStyle(2);
  fitDREB_above->SetLineColor(kRed);
  fitDREB_above->SetLineWidth(2);
  fitDREB_above->SetLineStyle(2);
  fitDREE_below->SetLineColor(kBlue);
  fitDREE_below->SetLineWidth(2);
  fitDREE_below->SetLineStyle(2);
  fitDREE_above->SetLineColor(kBlue);
  fitDREE_above->SetLineWidth(2);
  fitDREE_above->SetLineStyle(2);
  dr_3P0F_EB_rebin->Fit(fitDREB_below,"R");
  dr_3P0F_EB_rebin->Fit(fitDREB_above,"R+");
  dr_3P0F_EE_rebin->Fit(fitDREE_below,"R");
  dr_3P0F_EE_rebin->Fit(fitDREE_above,"R+");

  // PPFF

  TH1D* dr_3P1F_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P1F_CR_dr_EB")->Clone();
  TH1D* dr_PPFF_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P2F_CR_PPFF_dr_EB")->Clone();
  TH1D* dr_PFPFxFF_EB = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P2F_CR_PFPF_dr_EB_xFF")->Clone();
  TH1D* dr_3P1F_WZ_EB = (TH1D*)WZfile->Get("resolvedEleCRanalyzer/3P1F_CR_dr_EB")->Clone();
  TH1D* dr_3P1F_ZZ_EB = (TH1D*)ZZfile->Get("resolvedEleCRanalyzer/3P1F_CR_dr_EB")->Clone();
  int nbinsPPFF_drEB = 6;
  double xbinsPPFF_drEB[7] = {0.0,0.1,0.2,0.3,1.0,3.15,6.4};
  TH1D* dr_3P1F_EB_rebin = (TH1D*)dr_3P1F_EB->Rebin(nbinsPPFF_drEB,"dr_3P1F_EB_rebin",xbinsPPFF_drEB);
  TH1D* dr_PPFF_EB_rebin = (TH1D*)dr_PPFF_EB->Rebin(nbinsPPFF_drEB,"dr_PPFF_EB_rebin",xbinsPPFF_drEB);
  TH1D* dr_PFPFxFF_EB_rebin = (TH1D*)dr_PFPFxFF_EB->Rebin(nbinsPPFF_drEB,"dr_PFPFxFF_EB_rebin",xbinsPPFF_drEB);
  TH1D* dr_3P1F_WZ_EB_rebin = (TH1D*)dr_3P1F_WZ_EB->Rebin(nbinsPPFF_drEB,"dr_3P1F_WZ_EB_rebin",xbinsPPFF_drEB);
  TH1D* dr_3P1F_ZZ_EB_rebin = (TH1D*)dr_3P1F_ZZ_EB->Rebin(nbinsPPFF_drEB,"dr_3P1F_ZZ_EB_rebin",xbinsPPFF_drEB);
  const double lumi = 16.8;
  dr_3P1F_WZ_EB_rebin->Scale( 5.213*lumi*1000./ ((TH1D*)WZfile->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );
  dr_3P1F_ZZ_EB_rebin->Scale( 1.325*lumi*1000./ ((TH1D*)ZZfile->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );

  for (unsigned ibin = 0; ibin < dr_3P1F_EB_rebin->GetNbinsX()+2; ibin++) {
    if ( dr_3P1F_EB_rebin->GetBinContent(ibin)==0. )
      continue;

    auto square = [] (double x) { return x*x; };
    double val = std::max(dr_3P1F_EB_rebin->GetBinContent(ibin) - dr_PFPFxFF_EB_rebin->GetBinContent(ibin)
                          - dr_3P1F_WZ_EB_rebin->GetBinContent(ibin) - dr_3P1F_ZZ_EB_rebin->GetBinContent(ibin) ,0.);
    double err = std::sqrt( square(dr_3P1F_EB_rebin->GetBinError(ibin)) + square(dr_PFPFxFF_EB_rebin->GetBinError(ibin))
                            + square(dr_3P1F_WZ_EB_rebin->GetBinError(ibin)) + square(dr_3P1F_ZZ_EB_rebin->GetBinError(ibin)) );
    dr_3P1F_EB_rebin->SetBinContent(ibin,val);
    dr_3P1F_EB_rebin->SetBinError(ibin,err);
  }

  dr_3P1F_EB_rebin->Divide( dr_PPFF_EB_rebin );

  TH1D* dr_3P1F_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/3P1F_CR_dr_EE")->Clone();
  TH1D* dr_PPFF_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P2F_CR_PPFF_dr_EE")->Clone();
  TH1D* dr_PFPFxFF_EE = (TH1D*)datafile->Get("resolvedEleCRanalyzerData/2P2F_CR_PFPF_dr_EE_xFF")->Clone();
  TH1D* dr_3P1F_EE_rebin = (TH1D*)dr_3P1F_EE->Rebin(nbinsPPFF_drEB,"dr_3P1F_EE_rebin",xbinsPPFF_drEB);
  TH1D* dr_PPFF_EE_rebin = (TH1D*)dr_PPFF_EE->Rebin(nbinsPPFF_drEB,"dr_PPFF_EE_rebin",xbinsPPFF_drEB);
  TH1D* dr_PFPFxFF_EE_rebin = (TH1D*)dr_PFPFxFF_EE->Rebin(nbinsPPFF_drEB,"dr_PFPFxFF_EE_rebin",xbinsPPFF_drEB);
  TH1D* dr_3P1F_WZ_EE = (TH1D*)WZfile->Get("resolvedEleCRanalyzer/3P1F_CR_dr_EE")->Clone();
  TH1D* dr_3P1F_ZZ_EE = (TH1D*)ZZfile->Get("resolvedEleCRanalyzer/3P1F_CR_dr_EE")->Clone();
  TH1D* dr_3P1F_WZ_EE_rebin = (TH1D*)dr_3P1F_WZ_EE->Rebin(nbinsPPFF_drEB,"dr_3P1F_WZ_EE_rebin",xbinsPPFF_drEB);
  TH1D* dr_3P1F_ZZ_EE_rebin = (TH1D*)dr_3P1F_ZZ_EE->Rebin(nbinsPPFF_drEB,"dr_3P1F_ZZ_EE_rebin",xbinsPPFF_drEB);
  dr_3P1F_WZ_EE_rebin->Scale( 5.213*lumi*1000./ ((TH1D*)WZfile->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );
  dr_3P1F_ZZ_EE_rebin->Scale( 1.325*lumi*1000./ ((TH1D*)ZZfile->Get("resolvedEleCRanalyzer/totWeightedSum"))->GetBinContent(1) );

  for (unsigned ibin = 0; ibin < dr_3P1F_EE_rebin->GetNbinsX()+2; ibin++) {
    if ( dr_3P1F_EE_rebin->GetBinContent(ibin)==0. )
      continue;

    auto square = [] (double x) { return x*x; };
    double val = std::max(dr_3P1F_EE_rebin->GetBinContent(ibin) - dr_PFPFxFF_EE_rebin->GetBinContent(ibin)
                          - dr_3P1F_WZ_EE_rebin->GetBinContent(ibin) - dr_3P1F_ZZ_EE_rebin->GetBinContent(ibin), 0.);
    double err = std::sqrt( square(dr_3P1F_EE_rebin->GetBinError(ibin)) + square(dr_PFPFxFF_EE_rebin->GetBinError(ibin))
                            + square(dr_3P1F_WZ_EE_rebin->GetBinError(ibin)) + square(dr_3P1F_ZZ_EE_rebin->GetBinError(ibin)) );
    dr_3P1F_EE_rebin->SetBinContent(ibin,val);
    dr_3P1F_EE_rebin->SetBinError(ibin,err);
  }

  dr_3P1F_EE_rebin->Divide( dr_PPFF_EE_rebin );

  TF1* fitPPFF_DREB_below = new TF1("REFF_PPFF_dr_below_EB","[0]",0,0.3);
  TF1* fitPPFF_DREE_below = new TF1("REFF_PPFF_dr_below_EE","[0]",0,0.3);
  TF1* fitPPFF_DREB_above = new TF1("REFF_PPFF_dr_above_EB","[0]",0.3,6.4);
  TF1* fitPPFF_DREE_above = new TF1("REFF_PPFF_dr_above_EE","[0]",0.3,6.4);
  fitPPFF_DREB_below->SetLineColor(kRed);
  fitPPFF_DREB_below->SetLineWidth(2);
  fitPPFF_DREB_below->SetLineStyle(2);
  fitPPFF_DREB_above->SetLineColor(kRed);
  fitPPFF_DREB_above->SetLineWidth(2);
  fitPPFF_DREB_above->SetLineStyle(2);
  fitPPFF_DREE_below->SetLineColor(kBlue);
  fitPPFF_DREE_below->SetLineWidth(2);
  fitPPFF_DREE_below->SetLineStyle(2);
  fitPPFF_DREE_above->SetLineColor(kBlue);
  fitPPFF_DREE_above->SetLineWidth(2);
  fitPPFF_DREE_above->SetLineStyle(2);
  dr_3P1F_EB_rebin->Fit(fitPPFF_DREB_below,"R");
  dr_3P1F_EB_rebin->Fit(fitPPFF_DREB_above,"R+");
  dr_3P1F_EE_rebin->Fit(fitPPFF_DREE_below,"R");
  dr_3P1F_EE_rebin->Fit(fitPPFF_DREE_above,"R+");

  // save file
  TFile* outfile = new TFile("REFF_20UL16.root","RECREATE");
  EB_3P0F_rebin->Write();
  EE_3P0F_rebin->Write();
  dr_3P0F_EB_rebin->Write();
  dr_3P0F_EE_rebin->Write();
  fitEB->Write();
  fitEE->Write();
  fitDREB_below->Write();
  fitDREE_below->Write();
  fitDREB_above->Write();
  fitDREE_above->Write();
  dr_3P1F_EB_rebin->Write();
  dr_3P1F_EE_rebin->Write();
  fitPPFF_DREB_below->Write();
  fitPPFF_DREE_below->Write();
  fitPPFF_DREB_above->Write();
  fitPPFF_DREE_above->Write();
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

  EB_3P0F_rebin->GetYaxis()->SetRangeUser(0.,0.5);
  EB_3P0F_rebin->SetLineWidth(2);
  EB_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  EB_3P0F_rebin->GetXaxis()->SetTitle("E_{T} [GeV]");
  EB_3P0F_rebin->SetStats(0);
  EB_3P0F_rebin->SetLineColor(kRed);
  EB_3P0F_rebin->Draw("E1");

  EE_3P0F_rebin->SetLineColor(kBlue);
  EE_3P0F_rebin->SetLineWidth(2);
  EE_3P0F_rebin->Draw("same&E1");

  TPaveText* textEB = new TPaveText(0.15,0.6,0.4,0.75,"NDC");
  textEB->SetBorderSize(0);
  textEB->SetFillColor(0);
  TString textEBstr;
  textEBstr.Form("%.3g #pm %.3g (EB)", fitEB->GetParameter(0), fitEB->GetParError(0));
  textEB->AddText(textEBstr);
  ((TText*)textEB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textEB->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* textEE = new TPaveText(0.7,0.6,0.92,0.75,"NDC");
  textEE->SetBorderSize(0);
  textEE->SetFillColor(0);
  TString textEEstr;
  textEEstr.Form("%.3g #pm %.3g (EE)", fitEE->GetParameter(0), fitEE->GetParError(0));
  textEE->AddText(textEEstr);
  ((TText*)textEE->GetListOfLines()->Last())->SetTextColor(kBlue);

  textEB->Draw();
  textEE->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("REFF_20UL16.png");

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

  TPaveText* textDREB = new TPaveText(0.15,0.6,0.5,0.75,"NDC");
  textDREB->SetBorderSize(0);
  textDREB->SetFillColor(0);
  TString textDREBstrBelow;
  textDREBstrBelow.Form("%.3g #pm %.3g (EB) (#Delta R < 0.3)", fitDREB_below->GetParameter(0), fitDREB_below->GetParError(0));
  TString textDREBstrAbove;
  textDREBstrAbove.Form("%.3g #pm %.3g (EB) (#Delta R > 0.3)", fitDREB_above->GetParameter(0), fitDREB_above->GetParError(0));
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
  textDREEstrBelow.Form("%.3g #pm %.3g (EE) (#Delta R < 0.3)", fitDREE_below->GetParameter(0), fitDREE_below->GetParError(0));
  TString textDREEstrAbove;
  textDREEstrAbove.Form("%.3g #pm %.3g (EE) (#Delta R > 0.3)", fitDREE_above->GetParameter(0), fitDREE_above->GetParError(0));
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

  canvas->SaveAs("REFF_dr_20UL16.png");

  // draw PPFF

  dr_3P1F_EB_rebin->Draw("E1");

  dr_3P1F_EB_rebin->GetYaxis()->SetRangeUser(0.,1.0);
  dr_3P1F_EB_rebin->SetLineWidth(2);
  dr_3P1F_EB_rebin->GetYaxis()->SetTitle("Fake factor");
  dr_3P1F_EB_rebin->GetXaxis()->SetTitle("#Delta R");
  dr_3P1F_EB_rebin->SetStats(0);
  dr_3P1F_EB_rebin->SetLineColor(kRed);
  dr_3P1F_EB_rebin->Draw("E1");

  dr_3P1F_EE_rebin->SetLineColor(kBlue);
  dr_3P1F_EE_rebin->SetLineWidth(2);
  dr_3P1F_EE_rebin->Draw("same&E1");

  TPaveText* textPPFF_DREB = new TPaveText(0.15,0.6,0.5,0.75,"NDC");
  textPPFF_DREB->SetBorderSize(0);
  textPPFF_DREB->SetFillColor(0);
  TString textPPFF_DREBstrBelow;
  textPPFF_DREBstrBelow.Form("%.3g #pm %.3g (EB) (#Delta R < 0.3)", fitPPFF_DREB_below->GetParameter(0), fitPPFF_DREB_below->GetParError(0));
  TString textPPFF_DREBstrAbove;
  textPPFF_DREBstrAbove.Form("%.3g #pm %.3g (EB) (#Delta R > 0.3)", fitPPFF_DREB_above->GetParameter(0), fitPPFF_DREB_above->GetParError(0));
  textPPFF_DREB->AddText(textPPFF_DREBstrBelow);
  ((TText*)textPPFF_DREB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textPPFF_DREB->GetListOfLines()->Last())->SetTextAlign(12);
  textPPFF_DREB->AddText(textPPFF_DREBstrAbove);
  ((TText*)textPPFF_DREB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textPPFF_DREB->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* textPPFF_DREE = new TPaveText(0.55,0.6,0.92,0.75,"NDC");
  textPPFF_DREE->SetBorderSize(0);
  textPPFF_DREE->SetFillColor(0);
  TString textPPFF_DREEstrBelow;
  textPPFF_DREEstrBelow.Form("%.3g #pm %.3g (EE) (#Delta R < 0.3)", fitPPFF_DREE_below->GetParameter(0), fitPPFF_DREE_below->GetParError(0));
  TString textPPFF_DREEstrAbove;
  textPPFF_DREEstrAbove.Form("%.3g #pm %.3g (EE) (#Delta R > 0.3)", fitPPFF_DREE_above->GetParameter(0), fitPPFF_DREE_above->GetParError(0));
  textPPFF_DREE->AddText(textPPFF_DREEstrBelow);
  ((TText*)textPPFF_DREE->GetListOfLines()->Last())->SetTextColor(kBlue);
  textPPFF_DREE->AddText(textPPFF_DREEstrAbove);
  ((TText*)textPPFF_DREE->GetListOfLines()->Last())->SetTextColor(kBlue);

  textPPFF_DREB->Draw();
  textPPFF_DREE->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();

  canvas->SaveAs("REFF_PPFF_dr_20UL16.png");

  return;
}
