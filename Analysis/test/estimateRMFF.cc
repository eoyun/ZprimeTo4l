#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void estimateRMFF(TString era) {
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

  TFile* datafile = new TFile("RMCR_"+era+"_data.root","READ");

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

  TH1D* MB_3P0F = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/3P0F_pt_MB")->Clone();
  TH1D* MB_2P1F = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/2P1F_pt_MB")->Clone();
  int nbinsMB = 12;
  double xbinsMB[13] = {0,10,20,35,50,75,100,150,200,300,500,700,1000};
  TH1D* MB_3P0F_rebin = (TH1D*)MB_3P0F->Rebin(nbinsMB, "MB_3P0F_rebin", xbinsMB);
  TH1D* MB_2P1F_rebin = (TH1D*)MB_2P1F->Rebin(nbinsMB, "MB_2P1F_rebin", xbinsMB);

  TH1D* ME_3P0F = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/3P0F_pt_ME")->Clone();
  TH1D* ME_2P1F = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/2P1F_pt_ME")->Clone();
  int nbinsME = 12;
  double xbinsME[13] = {0,10,20,35,50,75,100,150,200,300,500,700,1000};
  TH1D* ME_3P0F_rebin = (TH1D*)ME_3P0F->Rebin(nbinsME, "ME_3P0F_rebin", xbinsME);
  TH1D* ME_2P1F_rebin = (TH1D*)ME_2P1F->Rebin(nbinsME, "ME_2P1F_rebin", xbinsME);

  TH1D* all_3P0F_rebin = (TH1D*)MB_3P0F_rebin->Clone();
  all_3P0F_rebin->Add(ME_3P0F_rebin);

  TH1D* all_2P1F_rebin = (TH1D*)MB_2P1F_rebin->Clone();
  all_2P1F_rebin->Add(ME_2P1F_rebin);

  MB_3P0F_rebin->Divide( MB_2P1F_rebin );
  ME_3P0F_rebin->Divide( ME_2P1F_rebin );

  all_3P0F_rebin->Divide( all_2P1F_rebin );

  TF1* fitMB = new TF1("RMFF_MB","[0]",0,1000);
  TF1* fitME = new TF1("RMFF_ME","[0]",0,1000);
  fitMB->SetLineColor(kRed);
  fitMB->SetLineWidth(2);
  fitMB->SetLineStyle(2);
  fitME->SetLineColor(kBlue);
  fitME->SetLineWidth(2);
  fitME->SetLineStyle(2);
  MB_3P0F_rebin->Fit(fitMB,"R");
  ME_3P0F_rebin->Fit(fitME,"R");

  TF1* fitAll = new TF1("RMFF_all","[0]",0,1000);
  fitAll->SetLineColor(kRed);
  fitAll->SetLineWidth(2);
  fitAll->SetLineStyle(2);
  // all_3P0F_rebin->Fit(fitAll,"R");

  // try dR

  TH1D* dr_3P0F_MB = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/3P0F_dr_MB")->Clone();
  TH1D* dr_2P1F_MB = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/2P1F_dr_MB")->Clone();
  std::vector<double> xbinsdrMB = {0.0,0.1,0.2,0.25,0.3,0.35,0.4,0.5,0.7,1.0,1.5,2.0,3.15,6.4};
  const int nbinsdrMB = xbinsdrMB.size()-1;
  std::vector<double> xcendrMB = estimateCenter(xbinsdrMB);
  TH1D* dr_3P0F_MB_rebin = (TH1D*)dr_3P0F_MB->Rebin(nbinsdrMB,"dr_3P0F_MB_rebin",&(xbinsdrMB[0]));
  TH1D* dr_2P1F_MB_rebin = (TH1D*)dr_2P1F_MB->Rebin(nbinsdrMB,"dr_2P1F_MB_rebin",&(xbinsdrMB[0]));

  TH1D* dr_3P0F_ME = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/3P0F_dr_ME")->Clone();
  TH1D* dr_2P1F_ME = (TH1D*)datafile->Get("resolvedMuCRanalyzerData/2P1F_dr_ME")->Clone();
  TH1D* dr_3P0F_ME_rebin = (TH1D*)dr_3P0F_ME->Rebin(nbinsdrMB,"dr_3P0F_ME_rebin",&(xbinsdrMB[0]));
  TH1D* dr_2P1F_ME_rebin = (TH1D*)dr_2P1F_ME->Rebin(nbinsdrMB,"dr_2P1F_ME_rebin",&(xbinsdrMB[0]));

  TH1D* all_dr_3P0F_rebin = (TH1D*)dr_3P0F_MB_rebin->Clone("all_dr_3P0F_rebin");
  all_dr_3P0F_rebin->Add(dr_3P0F_ME_rebin);

  TH1D* all_dr_2P1F_rebin = (TH1D*)dr_2P1F_MB_rebin->Clone();
  all_dr_2P1F_rebin->Add(dr_2P1F_ME_rebin);

  dr_3P0F_MB_rebin->Divide( dr_2P1F_MB_rebin );
  dr_3P0F_ME_rebin->Divide( dr_2P1F_ME_rebin );

  all_dr_3P0F_rebin->Divide(all_dr_2P1F_rebin);

  TF1* fitDRMB_above = new TF1("RMFF_dr_MB","x > 0.8 ? [0]*x+[1] : [0]*0.8+[1] + [2]*(x-0.8)",0.0,6.4);
  TF1* fitDRME_above = new TF1("RMFF_dr_ME","x > 0.8 ? [0]*x+[1] : [0]*0.8+[1] + [2]*(x-0.8)",0.0,6.4);
  TF1* fitDRall = new TF1("RMFF_dr_all","x > 0.8 ? [0]*x+[1] : [0]*0.8+[1] + [2]*(x-0.8)",0.0,6.4);
  fitDRMB_above->SetLineColor(kRed);
  fitDRMB_above->SetLineWidth(2);
  fitDRMB_above->SetLineStyle(2);
  fitDRME_above->SetLineColor(kBlue);
  fitDRME_above->SetLineWidth(2);
  fitDRME_above->SetLineStyle(2);

  fitDRall->SetLineColor(kRed);
  fitDRall->SetLineWidth(2);
  fitDRall->SetLineStyle(2);

  TFitResultPtr fitResultMB = dr_3P0F_MB_rebin->Fit(fitDRMB_above,"RS");
  fitResultMB->SetName("fitMB");
  double ciMB[nbinsdrMB+1];
  fitResultMB->GetConfidenceIntervals(nbinsdrMB+1,1,0,&(xcendrMB[0]),ciMB,0.95,false);

  std::vector<double> xbinwdrMB = estimateWidth(xbinsdrMB);
  double ybindrMB[nbinsdrMB+1];

  for (unsigned idx = 0; idx < nbinsdrMB+1; idx++) {
    ybindrMB[idx] = fitDRMB_above->Eval(xcendrMB[idx]);
  }

  auto errMB = new TGraphErrors(nbinsdrMB+1,&(xcendrMB[0]),&(ybindrMB[0]),&(xbinwdrMB[0]),ciMB);
  errMB->SetFillColor(kRed);
  errMB->SetFillStyle(3003);

  TFitResultPtr fitResultME = dr_3P0F_ME_rebin->Fit(fitDRME_above,"RS");
  fitResultME->SetName("fitME");
  double ciME[nbinsdrMB+1];
  fitResultME->GetConfidenceIntervals(nbinsdrMB+1,1,0,&(xcendrMB[0]),ciME,0.95,false);
  double ybindrME[nbinsdrMB+1];

  for (unsigned idx = 0; idx < nbinsdrMB+1; idx++) {
    ybindrME[idx] = fitDRME_above->Eval(xcendrMB[idx]);
  }

  auto errME = new TGraphErrors(nbinsdrMB+1,&(xcendrMB[0]),&(ybindrME[0]),&(xbinwdrMB[0]),ciME);
  errME->SetFillColor(kBlue);
  errME->SetFillStyle(3003);

  TFitResultPtr fitResultAll = all_dr_3P0F_rebin->Fit(fitDRall,"RS");
  fitResultAll->SetName("fitAll");
  double ciAll[nbinsdrMB+1];
  fitResultAll->GetConfidenceIntervals(nbinsdrMB+1,1,0,&(xcendrMB[0]),ciAll,0.95,false);
  double ybinAll[nbinsdrMB+1];

  for (unsigned idx = 0; idx < nbinsdrMB+1; idx++) {
    ybinAll[idx] = fitDRall->Eval(xcendrMB[idx]);
  }

  auto errAll = new TGraphErrors(nbinsdrMB+1,&(xcendrMB[0]),&(ybinAll[0]),&(xbinwdrMB[0]),ciAll);
  errAll->SetFillColor(kRed);
  errAll->SetFillStyle(3003);

  // save file
  TFile* outfile = new TFile("RMFF_"+era+".root","RECREATE");
  MB_3P0F_rebin->Write();
  ME_3P0F_rebin->Write();
  all_3P0F_rebin->Write();
  dr_3P0F_MB_rebin->Write();
  dr_3P0F_ME_rebin->Write();
  all_dr_3P0F_rebin->Write();
  //fitMB->Write();
  //fitME->Write();
  //fitAll->Write();
  fitDRMB_above->Write();
  fitDRME_above->Write();
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

  MB_3P0F_rebin->GetYaxis()->SetRangeUser(0.,0.1);
  MB_3P0F_rebin->SetLineWidth(2);
  MB_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  MB_3P0F_rebin->GetXaxis()->SetTitle("p_{T} [GeV]");
  MB_3P0F_rebin->SetStats(0);
  MB_3P0F_rebin->SetLineColor(kRed);
  MB_3P0F_rebin->Draw("E1");

  ME_3P0F_rebin->SetLineColor(kBlue);
  ME_3P0F_rebin->SetLineWidth(2);
  ME_3P0F_rebin->Draw("same&E1");

  TPaveText* textMB = new TPaveText(0.15,0.6,0.4,0.75,"NDC");
  textMB->SetBorderSize(0);
  textMB->SetFillColor(0);
  TString textMBstr;
  textMBstr.Form("%.3g #pm %.3g (MB)", fitMB->GetParameter(0), fitMB->GetParError(0));
  textMB->AddText(textMBstr);
  ((TText*)textMB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textMB->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* textME = new TPaveText(0.7,0.6,0.92,0.75,"NDC");
  textME->SetBorderSize(0);
  textME->SetFillColor(0);
  TString textMEstr;
  textMEstr.Form("%.3g #pm %.3g (ME)", fitME->GetParameter(0), fitME->GetParError(0));
  textME->AddText(textMEstr);
  ((TText*)textME->GetListOfLines()->Last())->SetTextColor(kBlue);

  TLegend* legend = new TLegend(0.8,0.75,0.95,0.9);
  legend->SetBorderSize(0);
  legend->AddEntry(MB_3P0F_rebin,"MB");
  legend->AddEntry(ME_3P0F_rebin,"ME");
  legend->Draw();

  // textMB->Draw();
  // textME->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("RMFF_"+era+".png");

  // all
  all_3P0F_rebin->GetYaxis()->SetRangeUser(0.,0.1);
  all_3P0F_rebin->SetLineWidth(2);
  all_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  all_3P0F_rebin->GetXaxis()->SetTitle("p_{T} [GeV]");
  all_3P0F_rebin->SetStats(0);
  all_3P0F_rebin->SetLineColor(kRed);
  all_3P0F_rebin->Draw("E1");

  TPaveText* textAll = new TPaveText(0.15,0.6,0.4,0.75,"NDC");
  textAll->SetBorderSize(0);
  textAll->SetFillColor(0);
  TString textAllStr;
  textAllStr.Form("%.3g #pm %.3g (MB)", fitAll->GetParameter(0), fitAll->GetParError(0));
  textAll->AddText(textAllStr);
  ((TText*)textAll->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textAll->GetListOfLines()->Last())->SetTextAlign(12);

  // textAll->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();
  canvas->SaveAs("RMFF_all_"+era+".png");

  // draw dr

  dr_3P0F_MB_rebin->Draw("E1");

  dr_3P0F_MB_rebin->GetYaxis()->SetRangeUser(0.,0.1);
  dr_3P0F_MB_rebin->SetLineWidth(2);
  dr_3P0F_MB_rebin->GetYaxis()->SetTitle("Fake factor");
  dr_3P0F_MB_rebin->GetXaxis()->SetTitle("#Delta R");
  dr_3P0F_MB_rebin->SetStats(0);
  dr_3P0F_MB_rebin->SetLineColor(kRed);
  dr_3P0F_MB_rebin->Draw("E1");
  errMB->Draw("3");

  dr_3P0F_ME_rebin->SetLineColor(kBlue);
  dr_3P0F_ME_rebin->SetLineWidth(2);
  dr_3P0F_ME_rebin->Draw("same&E1");
  errME->Draw("3");

  TPaveText* textDRMB = new TPaveText(0.15,0.6,0.5,0.75,"NDC");
  textDRMB->SetBorderSize(0);
  textDRMB->SetFillColor(0);
  TString textDRMBstrBelow;
  textDRMBstrBelow.Form("%.3g#DeltaR + %.3g (MB) (#DeltaR < 0.8)", fitDRMB_above->GetParameter(2), fitDRMB_above->GetParameter(0)*0.8+fitDRMB_above->GetParameter(1)-fitDRMB_above->GetParameter(2)*0.8);
  TString textDRMBstrAbove;
  textDRMBstrAbove.Form("%.3g#DeltaR + %.3g (MB) (#DeltaR > 0.5)", fitDRMB_above->GetParameter(0), fitDRMB_above->GetParameter(1));
  textDRMB->AddText(textDRMBstrBelow);
  ((TText*)textDRMB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDRMB->GetListOfLines()->Last())->SetTextAlign(12);
  textDRMB->AddText(textDRMBstrAbove);
  ((TText*)textDRMB->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDRMB->GetListOfLines()->Last())->SetTextAlign(12);

  TPaveText* textDRME = new TPaveText(0.55,0.6,0.92,0.75,"NDC");
  textDRME->SetBorderSize(0);
  textDRME->SetFillColor(0);
  TString textDRMEstrBelow;
  textDRMEstrBelow.Form("%.3g#DeltaR + %.3g (ME) (#DeltaR < 0.8)", fitDRME_above->GetParameter(2), fitDRME_above->GetParameter(0)*0.8+fitDRME_above->GetParameter(1)-fitDRME_above->GetParameter(2)*0.8);
  TString textDRMEstrAbove;
  textDRMEstrAbove.Form("%.3g#DeltaR + %.3g (ME) (#DeltaR > 0.5)", fitDRME_above->GetParameter(0), fitDRME_above->GetParameter(1));
  textDRME->AddText(textDRMEstrBelow);
  ((TText*)textDRME->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textDRME->GetListOfLines()->Last())->SetTextAlign(12);
  textDRME->AddText(textDRMEstrAbove);
  ((TText*)textDRME->GetListOfLines()->Last())->SetTextColor(kBlue);
  ((TText*)textDRME->GetListOfLines()->Last())->SetTextAlign(12);

  textDRMB->Draw();
  textDRME->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();

  canvas->SaveAs("RMFF_dr_"+era+".png");

  // all
  all_dr_3P0F_rebin->Draw("E1");

  all_dr_3P0F_rebin->GetYaxis()->SetRangeUser(0.,0.1);
  all_dr_3P0F_rebin->SetLineWidth(2);
  all_dr_3P0F_rebin->GetYaxis()->SetTitle("Fake factor");
  all_dr_3P0F_rebin->GetXaxis()->SetTitle("#Delta R");
  all_dr_3P0F_rebin->SetStats(0);
  all_dr_3P0F_rebin->SetLineColor(kRed);
  all_dr_3P0F_rebin->Draw("E1");
  errAll->Draw("3");

  TPaveText* textDRall = new TPaveText(0.15,0.6,0.5,0.75,"NDC");
  textDRall->SetBorderSize(0);
  textDRall->SetFillColor(0);
  TString textDRallStrBelow;
  textDRallStrBelow.Form("%.3g#DeltaR + %.3g (#DeltaR < 0.8)", fitDRall->GetParameter(2), fitDRall->GetParameter(0)*0.8+fitDRall->GetParameter(1)-fitDRall->GetParameter(2)*0.8);
  TString textDRallStrAbove;
  textDRallStrAbove.Form("%.3g#DeltaR + %.3g (#DeltaR > 0.5)", fitDRall->GetParameter(0), fitDRall->GetParameter(1));
  textDRall->AddText(textDRallStrBelow);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextAlign(12);
  textDRall->AddText(textDRallStrAbove);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)textDRall->GetListOfLines()->Last())->SetTextAlign(12);

  textDRall->Draw();

  canvas->Update();

  // writing the lumi information and the CMS "logo"
  CMS_lumi( canvas, iPeriod, iPos );

  canvas->Update();
  canvas->RedrawAxis();
  canvas->GetFrame()->Draw();

  canvas->SaveAs("RMFF_all_dr_"+era+".png");

  return;
}
