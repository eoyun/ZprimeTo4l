#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runSF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  double valLumi = 0.;
  TString postfixAn = era;

  if (era=="20UL16APV") {
    lumi_13TeV = "19.5 fb^{-1}";
    valLumi = 19.5;
  } else if (era=="20UL16") {
    lumi_13TeV = "16.8 fb^{-1}";
    valLumi = 16.8;
    postfixAn = "";
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

  class sample {
  public:
    TFile* file;
    double xsec;
    int color;

    sample(TFile* afile, const double axsec, const int acolor) {
      file = afile;
      xsec = axsec;
      color = acolor;
    }
  };

  TFile* datafile = new TFile("MergedEleCR_"+era+"_data.root","READ");

  const double normNLO = 6077.22/6424.0;
  auto sample_DY_PtZ0To50 = sample(new TFile("MergedEleCR_"+era+"_DY_FXFX_PtZ-0To50.root","READ"),6424.0*normNLO,kOrange);
  auto sample_DY_PtZ50To100 = sample(new TFile("MergedEleCR_"+era+"_DY_FXFX_PtZ-50To100.root","READ"),397.4*normNLO,kOrange-4);
  auto sample_DY_PtZ100To250 = sample(new TFile("MergedEleCR_"+era+"_DY_FXFX_PtZ-100To250.root","READ"),97.2*normNLO,kOrange+6);
  auto sample_DY_PtZ250To400 = sample(new TFile("MergedEleCR_"+era+"_DY_FXFX_PtZ-250To400.root","READ"),3.701*normNLO,kOrange-3);
  auto sample_DY_PtZ400To650 = sample(new TFile("MergedEleCR_"+era+"_DY_FXFX_PtZ-400To650.root","READ"),0.5086*normNLO,kOrange+7);
  auto sample_DY_PtZ650ToInf = sample(new TFile("MergedEleCR_"+era+"_DY_FXFX_PtZ-650ToInf.root","READ"),0.04728*normNLO,kOrange-5);

  std::vector<sample> samples_DY_NLO = {
    sample_DY_PtZ0To50,
    sample_DY_PtZ50To100,
    sample_DY_PtZ100To250,
    sample_DY_PtZ250To400,
    sample_DY_PtZ400To650,
    sample_DY_PtZ650ToInf
  };

  const double normLO = 6077.22/5379.0;
  auto sample_DY_HT0To70 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-0To70.root","READ"),5379.0*normLO,kOrange);
  auto sample_DY_HT70To100 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-70To100.root","READ"),140.0*normLO,kOrange-4);
  auto sample_DY_HT100To200 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-100To200.root","READ"),139.2*normLO,kOrange+6);
  auto sample_DY_HT200To400 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-200To400.root","READ"),38.4*normLO,kOrange-3);
  auto sample_DY_HT400To600 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-400To600.root","READ"),5.174*normLO,kOrange+7);
  auto sample_DY_HT600To800 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-600To800.root","READ"),1.258*normLO,kOrange-5);
  auto sample_DY_HT800To1200 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-800To1200.root","READ"),0.5598*normLO,kOrange+5);
  auto sample_DY_HT1200To2500 = sample(new TFile("MergedEleCR_"+era+"_DY_HT-1200To2500.root","READ"),0.1305*normLO,kOrange-6);
  auto sample_DY_HT2500ToInf = sample(new TFile("MergedEleCR_"+era+"_DY_HT-2500ToInf.root","READ"),0.002997*normLO,kOrange+4);

  std::vector<sample> samples_DY_LO = {
    sample_DY_HT0To70,
    sample_DY_HT70To100,
    sample_DY_HT100To200,
    sample_DY_HT200To400,
    sample_DY_HT400To600,
    sample_DY_HT600To800,
    sample_DY_HT800To1200,
    sample_DY_HT1200To2500,
    sample_DY_HT2500ToInf
  };

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

  TString anlyzrMC = "mergedEleCRanalyzer"+postfixAn;
  TString anlyzrData = "mergedEleCRanalyzerData";

  // THStack* stack_VpT = new THStack("stack_VpT",";GeV;");
  //
  // for (auto& asample : samples_DY_NLO) {
  //   TH1D* ahist = (TH1D*)asample.file->Get(anlyzrMC+"/lheVpT_cut");
  //   ahist->SetFillColor(asample.color);
  //   ahist->SetLineColor(asample.color);
  //   stack_VpT->Add(ahist);
  // }

  // canvas_1->cd();
  // stack_VpT->Draw("hist");
  // SaveAs(canvas_1,"SF_VpT.png");

  // THStack* stack_HT = new THStack("stack_HT",";GeV;");
  //
  // for (auto& asample : samples_DY_LO) {
  //   TH1D* ahist = (TH1D*)asample.file->Get(anlyzrMC+"/lheHT_cut");
  //   ahist->SetFillColor(asample.color);
  //   ahist->SetLineColor(asample.color);
  //   stack_HT->Add(ahist);
  // }

  // stack_HT->Draw("hist");
  // SaveAs(canvas_1,"SF_HT.png");

  auto drawRatio = [] (const TH1D* numer, const TH1D* denom, TPad* pad, TString postfix="") {
    TH1D* ratio = (TH1D*)numer->Clone();
    ratio->SetStats(0);
    ratio->SetTitle("");
    ratio->Divide(denom);
    ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);

    if (postfix=="_EB")
      ratio->GetYaxis()->SetRangeUser(0.7,1.3);

    if (postfix=="_EE")
      ratio->GetYaxis()->SetRangeUser(0.,3.);

    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetTitleOffset(0.75);
    ratio->SetLineColor(kBlack);

    pad->cd();
    ratio->Draw("E1");
  };

  auto drawMass = [&datafile,&anlyzrData,&anlyzrMC,&valLumi,&canvas_2,&p1,&p2,&SaveAs] (std::vector<sample> samples, TString name,TString postfix) {
    TH1D* h_data = (TH1D*)((TH1D*)datafile->Get(anlyzrData+"/"+name))->Clone();
    h_data->SetLineWidth(2);
    h_data->SetLineColor(kBlack);

    THStack* astack = new THStack("stack_"+name,";GeV;");

    for (auto asample = samples.rbegin(); asample!=samples.rend(); ++asample) {
      TH1D* ahist = (TH1D*)asample->file->Get(anlyzrMC+"/"+name);
      ahist->SetFillColor(asample->color);
      ahist->SetLineColor(asample->color);
      ahist->Scale(valLumi*1000.*asample->xsec/((TH1D*)asample->file->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1));
      astack->Add(ahist);
    }

    p1->cd();
    // p1->SetLogy();
    h_data->GetXaxis()->SetRangeUser(60,120);
    h_data->SetMaximum(1.5*h_data->GetMaximum());
    h_data->Draw("E1");
    astack->Draw("hist&sames");
    h_data->Draw("E1&sames");

    p2->cd();

    TList* stackHists = astack->GetHists();
    TH1D* tmpHist = (TH1D*)stackHists->At(0)->Clone();

    for (int idx = 1; idx < stackHists->GetSize(); ++idx)
      tmpHist->Add((TH1D*)stackHists->At(idx));

    TH1D* ratio = (TH1D*)h_data->Clone();
    ratio->Divide(tmpHist);
    ratio->GetXaxis()->SetRangeUser(60,120);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetYaxis()->SetTitle("Data/MC");
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetRangeUser(0.5,1.5);
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetTitleOffset(0.75);
    ratio->Draw("E1");

    SaveAs(canvas_2,("SF_"+name+"_"+postfix+".png").Data(),p1);
  };

  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_pass","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_fail","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_pass_EB","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_fail_EB","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_pass_EE","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_fail_EE","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_pass_Et-50to100","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_fail_Et-50to100","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_pass_Et-100to300","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_fail_Et-100to300","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_pass_Et-300to1000","VpT");
  drawMass(samples_DY_NLO,"2E_SFCR_OSll_invM_fail_Et-300to1000","VpT");

  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_pass","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_fail","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_pass_EB","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_fail_EB","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_pass_EE","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_fail_EE","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_pass_Et-50to100","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_fail_Et-50to100","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_pass_Et-100to300","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_fail_Et-100to300","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_pass_Et-300to1000","HT");
  drawMass(samples_DY_LO,"2E_SFCR_OSll_invM_fail_Et-300to1000","HT");

  const unsigned nbinsOS = 11;
  double xbinsOS[nbinsOS+1] = {0, 50, 60, 70, 80, 90, 100, 125, 150, 200,
                               300, 1000};

  const unsigned nbinsOS_EE = 4;
  double xbinsOS_EE[nbinsOS_EE+1] = {0, 50, 70, 100, 1000};

  auto addHists = [&] (TH1D* ahist, std::vector<sample> samples, TString name) {
    for (unsigned idx = 1; idx < samples.size(); idx++) {
      TH1D* cloned = (TH1D*)((TH1D*)samples.at(idx).file->Get(anlyzrMC+"/"+name))->Clone();
      cloned->Scale(valLumi*1000.*samples.at(idx).xsec/((TH1D*)samples.at(idx).file->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1));
      ahist->Add(cloned);
    }
  };

  auto drawSF = [&] (TString postfix) {
    TH1D* numerData = (TH1D*)((TH1D*)datafile->Get(anlyzrData+"/2E_SFCR_OSll_Et_numer"+postfix))->Clone();
    TH1D* denomData = (TH1D*)((TH1D*)datafile->Get(anlyzrData+"/2E_SFCR_OSll_Et_denom"+postfix))->Clone();
    TH1D* numerData_rebin = (TH1D*)numerData->Rebin(nbinsOS, "2E_SFCR_OSll_Et_numer_rebin"+postfix, xbinsOS);
    TH1D* denomData_rebin = (TH1D*)denomData->Rebin(nbinsOS, "2E_SFCR_OSll_Et_denom_rebin"+postfix, xbinsOS);

    if (postfix=="_EE") {
      numerData_rebin = (TH1D*)numerData->Rebin(nbinsOS_EE, "2E_SFCR_OSll_Et_numer_rebin"+postfix, xbinsOS_EE);
      denomData_rebin = (TH1D*)denomData->Rebin(nbinsOS_EE, "2E_SFCR_OSll_Et_denom_rebin"+postfix, xbinsOS_EE);
    }

    numerData_rebin->Divide(denomData_rebin);

    TH1D* numerMC_NLO = (TH1D*)((TH1D*)samples_DY_NLO.at(0).file->Get(anlyzrMC+"/2E_SFCR_OSll_Et_numer"+postfix))->Clone();
    numerMC_NLO->Scale(valLumi*1000.*samples_DY_NLO.at(0).xsec/((TH1D*)samples_DY_NLO.at(0).file->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1));
    TH1D* denomMC_NLO = (TH1D*)((TH1D*)samples_DY_NLO.at(0).file->Get(anlyzrMC+"/2E_SFCR_OSll_Et_denom"+postfix))->Clone();
    denomMC_NLO->Scale(valLumi*1000.*samples_DY_NLO.at(0).xsec/((TH1D*)samples_DY_NLO.at(0).file->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1));

    addHists(numerMC_NLO,samples_DY_NLO,"2E_SFCR_OSll_Et_numer"+postfix);
    addHists(denomMC_NLO,samples_DY_NLO,"2E_SFCR_OSll_Et_denom"+postfix);

    TH1D* numerMC_LO = (TH1D*)((TH1D*)samples_DY_LO.at(0).file->Get(anlyzrMC+"/2E_SFCR_OSll_Et_numer"+postfix))->Clone();
    numerMC_LO->Scale(valLumi*1000.*samples_DY_LO.at(0).xsec/((TH1D*)samples_DY_LO.at(0).file->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1));
    TH1D* denomMC_LO = (TH1D*)((TH1D*)samples_DY_LO.at(0).file->Get(anlyzrMC+"/2E_SFCR_OSll_Et_denom"+postfix))->Clone();
    denomMC_LO->Scale(valLumi*1000.*samples_DY_LO.at(0).xsec/((TH1D*)samples_DY_LO.at(0).file->Get(anlyzrMC+"/totWeightedSum"))->GetBinContent(1));

    addHists(numerMC_LO,samples_DY_LO,"2E_SFCR_OSll_Et_numer"+postfix);
    addHists(denomMC_LO,samples_DY_LO,"2E_SFCR_OSll_Et_denom"+postfix);

    TH1D* numerMC_NLO_rebin = (TH1D*)numerMC_NLO->Rebin(nbinsOS, "2E_SFCR_OSll_Et_numer_rebin_NLO"+postfix, xbinsOS);
    TH1D* denomMC_NLO_rebin = (TH1D*)denomMC_NLO->Rebin(nbinsOS, "2E_SFCR_OSll_Et_denom_rebin_NLO"+postfix, xbinsOS);

    if (postfix=="_EE") {
      numerMC_NLO_rebin = (TH1D*)numerMC_NLO->Rebin(nbinsOS_EE, "2E_SFCR_OSll_Et_numer_rebin_NLO"+postfix, xbinsOS_EE);
      denomMC_NLO_rebin = (TH1D*)denomMC_NLO->Rebin(nbinsOS_EE, "2E_SFCR_OSll_Et_denom_rebin_NLO"+postfix, xbinsOS_EE);
    }

    numerMC_NLO_rebin->Divide(denomMC_NLO_rebin);

    TH1D* numerMC_LO_rebin = (TH1D*)numerMC_LO->Rebin(nbinsOS, "2E_SFCR_OSll_Et_numer_rebin_LO"+postfix, xbinsOS);
    TH1D* denomMC_LO_rebin = (TH1D*)denomMC_LO->Rebin(nbinsOS, "2E_SFCR_OSll_Et_denom_rebin_LO"+postfix, xbinsOS);

    if (postfix=="_EE") {
      numerMC_LO_rebin = (TH1D*)numerMC_LO->Rebin(nbinsOS_EE, "2E_SFCR_OSll_Et_numer_rebin_LO"+postfix, xbinsOS_EE);
      denomMC_LO_rebin = (TH1D*)denomMC_LO->Rebin(nbinsOS_EE, "2E_SFCR_OSll_Et_denom_rebin_LO"+postfix, xbinsOS_EE);
    }

    numerMC_LO_rebin->Divide(denomMC_LO_rebin);

    p1->cd();
    p1->SetLogx();
    numerData_rebin->GetYaxis()->SetTitle("Efficiency");
    numerData_rebin->SetLineColor(kBlack);
    numerData_rebin->SetLineWidth(2);
    numerData_rebin->GetYaxis()->SetTitleSize(0.05);
    numerData_rebin->GetYaxis()->SetTitleOffset(0.7);
    numerData_rebin->Draw("E1");
    numerMC_NLO_rebin->SetLineColor(kGray);
    numerMC_NLO_rebin->SetLineWidth(2);
    numerMC_NLO_rebin->Draw("E1&same");
    numerMC_LO_rebin->SetLineColor(kRed-3);
    numerMC_LO_rebin->SetLineWidth(2);
    numerMC_LO_rebin->Draw("E1&same");

    auto legend = std::make_unique<TLegend>(0.15,0.75,0.35,0.9);
    legend->SetBorderSize(0);
    legend->AddEntry(numerData_rebin,"Data");
    legend->AddEntry(numerMC_NLO_rebin,"DY NLO");
    legend->AddEntry(numerMC_LO_rebin,"DY LO");
    legend->Draw();

    p2->cd();
    p2->SetLogx();
    drawRatio(numerData_rebin,numerMC_NLO_rebin,p2,postfix);
    TH1D* cloned_LO = (TH1D*)numerData_rebin->Clone();
    cloned_LO->Divide(numerMC_LO_rebin);
    cloned_LO->SetLineColor(kRed-3);
    cloned_LO->SetLineWidth(2);
    cloned_LO->Draw("E1&sames");
    SaveAs(canvas_2,("SF"+postfix+".png").Data(),p1);
  };

  drawSF("");
  drawSF("_EB");
  drawSF("_EE");

  return;
}
