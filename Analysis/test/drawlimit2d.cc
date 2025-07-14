#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/tdrstyle.C"
#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/CMS_lumi.C"

void drawlimit2d(TString isZA) {
  TString era="run2";
  setTDRStyle();
  // gStyle->SetLineWidth(2);

  writeExtraText = true;       // if extra text
  extraText  = "";  // default extra text is "Preliminary"

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
  } else if (era=="run2") {
    lumi_sqrtS = "";
    lumi_13TeV = "138 fb^{-1}";
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
  float R = 0.06*W_ref;

  auto* canvas1 = new TCanvas("canvas","canvas",800,800,W,H);
  canvas1->SetFillColor(0);
  canvas1->SetBorderMode(0);
  canvas1->SetFrameFillStyle(0);
  canvas1->SetFrameBorderMode(0);
  canvas1->SetLeftMargin( L/W );
  canvas1->SetRightMargin( 0.15 );
  canvas1->SetTopMargin( T/H );
  canvas1->SetBottomMargin( B/H );
  canvas1->SetTickx(0);
  canvas1->SetTicky(0);

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

  auto retrieve = [&isZA] (const TString& hmass, const TString& amass) -> double {
    TString filename = TString("./")+isZA+amass+"/"+hmass+"/higgsCombine.X"+hmass+isZA+amass+".HybridNew.mH"+hmass+".root";

    auto afile = std::make_unique<TFile>(filename,"READ");

    if (afile->IsZombie() || afile->TestBit(TFile::kRecovered))
      return -1.;

    TTree* atree = (TTree*)afile->Get("limit");
    double val;
    atree->SetBranchAddress("limit",&val);

    unsigned nentry = atree->GetEntries();

    if (nentry==0)
      return 9999.;

    atree->GetEntry(0);
    double valFinal = val;

    return 0.1*valFinal;
  };

  std::vector<int> massvec = {250,275,300,325,350,375,400,425,450,500,550,650,750,850,1000,1250,1500,1750,2000};
  std::vector<TString> avec = {"0p4","0p6","0p8","1","1p5","2","5","10","50","100","250","500","750"};

  std::vector<float> xbin = {250.,275.,300.,325.,350.,375.,400.,425.,450.,500.,550.,650.,750.,850.,1000.,1250.,1500.,1750.,2000.,2500.};
  std::vector<float> ybin = {0.4,0.6,0.8,1.,1.5,2.,5.,10.,50.,100.,250.,500.,750.,1250.};

  if (isZA=="Z") {
    avec = {"0p4","0p6","0p8","1","1p5","2","5","10","50","100","250","500","750","1000","1500"};
    ybin = {0.4,0.6,0.8,1.,1.5,2.,5.,10.,50.,100.,250.,500.,750.,1000.,1500.,2000.};
  }

  TH2D* customHist = new TH2D("hist","hist",xbin.size()-1,&(xbin[0]),ybin.size()-1,&(ybin[0]));

  for (const auto hmass : massvec) {
    for (int idx = 0; idx<avec.size(); idx++) {
      auto ymassStr = avec.at(idx);
      double ymass = ybin.at(idx);
      double val = retrieve(std::to_string(hmass),ymassStr);

      if (val < 0.)
        continue;

      int ibin = customHist->FindFixBin(static_cast<float>(hmass)+1.,ymass+0.1);
      customHist->SetBinContent(ibin, val ); 
    }
  }

  auto fillGap = [&customHist] (bool fillY) {
    for (int idx=2; idx<=customHist->GetNbinsX(); idx++) {
      for (int jdx=2; jdx<=customHist->GetNbinsY(); jdx++) {
        double con = customHist->GetBinContent(idx,jdx);

        double conp1 = fillY ? customHist->GetBinContent(idx,jdx+1) : customHist->GetBinContent(idx+1,jdx);
        double conm1 = fillY ? customHist->GetBinContent(idx,jdx-1) : customHist->GetBinContent(idx-1,jdx);


        if (con == 0. && conp1 > 0. && conm1 > 0.)
          customHist->SetBinContent( customHist->GetBin(idx,jdx), (conp1+conm1)/2. );
      }
    }
  };

  if (isZA=="A") {
    int ibinX2000Y750 = customHist->FindFixBin(2000.+1.,750.+0.1);
    int ibinX2000Y250 = customHist->FindFixBin(2000.+1.,250.+0.1);
    int ibinX2000Y500 = customHist->FindFixBin(2000.+1.,500.+0.1);
    int ibinX2000Y100 = customHist->FindFixBin(2000.+1.,100.+0.1);
    customHist->SetBinContent(ibinX2000Y250, customHist->GetBinContent(ibinX2000Y100) + ( customHist->GetBinContent(ibinX2000Y750) - customHist->GetBinContent(ibinX2000Y100) )/3. );
    customHist->SetBinContent(ibinX2000Y500, customHist->GetBinContent(ibinX2000Y100) + 2.*( customHist->GetBinContent(ibinX2000Y750) - customHist->GetBinContent(ibinX2000Y100) )/3. );

    int ibinX1750Y500 = customHist->FindFixBin(1750.+1.,500.+0.1);
    int ibinX1750Y250 = customHist->FindFixBin(1750.+1.,250.+0.1);
    int ibinX1500Y500 = customHist->FindFixBin(1500.+1.,500.+0.1);
    int ibinX1500Y250 = customHist->FindFixBin(1500.+1.,250.+0.1);
    customHist->SetBinContent(ibinX1750Y250, (customHist->GetBinContent(ibinX1500Y250) + customHist->GetBinContent(ibinX2000Y250))/2. );
    customHist->SetBinContent(ibinX1750Y500, (customHist->GetBinContent(ibinX1500Y500) + customHist->GetBinContent(ibinX2000Y500))/2. );

    for (int idx=0; idx<massvec.size()-1; idx++) {
      int hmass = massvec.at(idx);

      if (hmass < 850)
        continue;

      int ibin = customHist->FindFixBin(static_cast<float>(hmass)+1.,250.+0.1);

      if ( customHist->GetBinContent(ibin)!=0. )
        continue;

      int ibinp1 = customHist->FindFixBin(static_cast<float>(massvec.at(idx+1))+1.,250.+0.1);
      int ibinm1 = customHist->FindFixBin(static_cast<float>(massvec.at(idx-1))+1.,250.+0.1);

      double val = (customHist->GetBinContent(ibinp1) + customHist->GetBinContent(ibinm1))/2.;
      customHist->SetBinContent(ibin, val);
    }
  } else if (isZA=="Z") {
    int ibinX2000Y1500 = customHist->FindFixBin(2000.+1.,1500.+0.1);
    int ibinX2000Y500 = customHist->FindFixBin(2000.+1.,500.+0.1);
    int ibinX2000Y1000 = customHist->FindFixBin(2000.+1.,1000.+0.1);
    int ibinX2000Y750 = customHist->FindFixBin(2000.+1.,750.+0.1);
    customHist->SetBinContent(ibinX2000Y750, customHist->GetBinContent(ibinX2000Y500) + ( customHist->GetBinContent(ibinX2000Y1500) - customHist->GetBinContent(ibinX2000Y500) )/3. );
    customHist->SetBinContent(ibinX2000Y1000, customHist->GetBinContent(ibinX2000Y500) + 2.*( customHist->GetBinContent(ibinX2000Y1500) - customHist->GetBinContent(ibinX2000Y500) )/3. );

    int ibinX750Y250 = customHist->FindFixBin(750.+1.,250.+0.1);
    int ibinX500Y250 = customHist->FindFixBin(500.+1.,250.+0.1);
    int ibinX550Y250 = customHist->FindFixBin(550.+1.,250.+0.1);
    int ibinX650Y250 = customHist->FindFixBin(650.+1.,250.+0.1);
    customHist->SetBinContent(ibinX550Y250, customHist->GetBinContent(ibinX500Y250) + ( customHist->GetBinContent(ibinX750Y250) - customHist->GetBinContent(ibinX500Y250) )/3. );
    customHist->SetBinContent(ibinX650Y250, customHist->GetBinContent(ibinX500Y250) + 2.*( customHist->GetBinContent(ibinX750Y250) - customHist->GetBinContent(ibinX500Y250) )/3. );

    fillGap(true);
    fillGap(false);
  }

  canvas1->SetLogy();
  canvas1->SetLogz();

  TString Ychar = isZA=="A" ? "Y" : "Z";

  customHist->GetXaxis()->SetLabelSize(0.035);
  customHist->GetXaxis()->SetTitleSize(0.05);
  customHist->GetXaxis()->SetTitleOffset(1.05);
  customHist->GetXaxis()->SetTitle("M_{X} [GeV]");
  customHist->GetYaxis()->SetLabelSize(0.04);
  customHist->GetYaxis()->SetTitleSize(0.05);
  customHist->GetYaxis()->SetTitleOffset(1.2);
  customHist->GetYaxis()->SetTitle("M_{Y} [GeV]");
  customHist->GetZaxis()->SetLabelSize(0.035);
  customHist->GetZaxis()->SetTitleSize(0.05);
  customHist->GetZaxis()->SetTitle("#sigma(pp#rightarrowX) #times B(X#rightarrow"+Ychar+"Y#rightarrow4l) [fb]");
  customHist->GetZaxis()->SetTitleOffset(0.98);
  customHist->GetZaxis()->SetRangeUser(0.02,10.);
  customHist->Draw("colz");
  canvas1->Update();

  canvas1->cd();

  SaveAs(canvas1,std::string(TString("limit2D_"+isZA+".pdf").Data()));

  return;
}
