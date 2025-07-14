#include "TFile.h"
#include "TH1D.h"

#include <regex>

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void drawEff(TString filename, TString strZA) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "    Simulation Preliminary";  // default extra text is "Preliminary"

  lumi_sqrtS = "13 TeV";
  lumi_13TeV = "";

  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

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

  auto canvas = std::make_unique<TCanvas>("canvas","canvas",800,800,W,H);
  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( 0.15 );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

  TFile* afile = TFile::Open(filename,"READ");
  auto keyList = afile->GetListOfKeys();

  for (auto keyObj : *keyList) {
    auto key = (TKey*)keyObj;

    if (TString(key->GetClassName())==TString("TDirectoryFile")) {
      TDirectory* adir = (TDirectory*)key->ReadObj();
      TString dirname = adir->GetName();
      auto keyList2 = adir->GetListOfKeys();

      std::vector<float> mxs;
      std::vector<float> mys;
      std::vector<float> effs;

      std::vector<int> massvec = {250,275,300,325,350,375,400,425,450,500,550,650,750,850,1000,1250,1500,1750,2000}; 
      std::vector<float> xbin = {250.,275.,300.,325.,350.,375.,400.,425.,450.,500.,550.,650.,750.,850.,1000.,1250.,1500.,1750.,2000.,2500.};
      std::vector<float> ybin = {0.4,0.6,0.8,1.,1.5,2.,5.,10.,50.,100.,250.,500.,750.,1250.};

      TH2D* customHist = new TH2D("hist","hist",xbin.size()-1,&(xbin[0]),ybin.size()-1,&(ybin[0]));

      TH1D* customHist1d = new TH1D("hist1d","hist1d",xbin.size()-1,&(xbin[0]));

      for (auto keyObj2 : *keyList2) {
        auto key2 = (TKey*)keyObj2;

        if (TString(key2->GetClassName()).Contains("TH1")) {
          std::string keyname = std::string(key2->GetName());
          std::regex sigPattern("H([0-9]+)([AZ])([0-9]*p*[0-9]+).*");
          std::regex sigSystPattern("H[0-9]+[AZ][0-9]*p*[0-9]+_(.*)");
          std::smatch XYmass, systName;

          bool isSig = std::regex_match(keyname, XYmass, sigPattern);
          bool isSyst = std::regex_match(keyname, systName, sigSystPattern);
          std::string ymassStr = std::regex_replace( XYmass[3].str(), std::regex("p"), "." );

          double mx = 0, my = 0;
          bool isZA = false;

          if (isSig) {
            mx = std::stod(XYmass[1].str());
            my = std::stod(ymassStr);
            isZA = XYmass[2].str()=="Z";
          }

          if (isSig && !isSyst && ( (isZA&&strZA=="Z") || (!isZA&&strZA=="A") )) {
            TH1* ahist = (TH1*)key2->ReadObj();
            mxs.push_back(static_cast<float>(mx));
            mys.push_back(static_cast<float>(my));

            int ibin = customHist->FindFixBin(static_cast<float>(mx)+1.,static_cast<float>(my)+0.1);
            customHist->SetBinContent(ibin, std::max(ahist->Integral()/(137.6*0.1),0.0001) );
            effs.push_back( ahist->Integral()/(137.6*0.1) );
            std::cout << "MX = " << mx << " MY = " << my << " axE = " << ahist->Integral()/(137.6*0.1) << std::endl;

            if (my==0.8) {
              int ibin1d = customHist1d->FindFixBin(static_cast<float>(mx)+1.);
              customHist1d->SetBinContent(ibin1d, std::max(ahist->Integral()/(137.6*0.1),0.0001));
            }

            if (my==250.) {
              for (auto x : xbin) {
                if (x < 800.)
                  continue;

                int jbin = customHist->FindFixBin(static_cast<float>(x)+1.,static_cast<float>(my)+0.1);
                customHist->SetBinContent(jbin, ahist->Integral()/(137.6*0.1) ); 
              }
            }
          }
        }
      } // loop hists

      for (const auto& xpt : xbin) {
        for (const auto& ypt : ybin) {
          int bin = customHist->FindFixBin(xpt,ypt);

          if ( ( strZA=="A" && ypt < 100. ) || ( strZA=="Z" && ypt < 100. ) || ( ypt==100. && ( filename.Contains("ME2E") || filename.Contains("MMFF") || filename.Contains("MEMu1M") || filename.Contains("MM2E") ) ) ) {
            if (customHist->GetBinContent(bin)==0.0)
              customHist->SetBinContent(bin, 0.0001);
          }
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

      TString isZA = strZA;

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

      auto effHisto = std::make_unique<TGraph2D>(effs.size(),&(mxs[0]),&(mys[0]),&(effs[0]));
      effHisto->SetNpy(375);
      effHisto->SetNpx(35);

      canvas->SetLogy();

      customHist->GetXaxis()->SetLabelSize(0.035);
      customHist->GetXaxis()->SetTitleSize(0.035);
      customHist->GetXaxis()->SetTitleOffset(1.2);
      customHist->GetXaxis()->SetTitle("M_{X} [GeV]");
      customHist->GetYaxis()->SetLabelSize(0.035);
      customHist->GetYaxis()->SetTitleSize(0.035);
      customHist->GetYaxis()->SetTitleOffset(1.4);
      customHist->GetYaxis()->SetTitle("M_{Y} [GeV]");
      customHist->GetZaxis()->SetLabelSize(0.035);
      customHist->GetZaxis()->SetTitleSize(0.035);
      customHist->GetZaxis()->SetTitleOffset(1.2);
      customHist->Draw("colz");
      canvas->Update();

      //effHisto->GetXaxis()->SetLimits(0.,2000.);
      //effHisto->GetHistogram()->GetYaxis()->SetRangeUser(1.,750.);
      customHist->GetZaxis()->SetRangeUser(0.,0.3);

      canvas->Update();

      CMS_lumi( canvas.get(), iPeriod, iPos );

      canvas->Update();
      canvas->RedrawAxis();
      canvas->GetFrame()->Draw();
      canvas->SaveAs(dirname+"_"+strZA+".png");

      customHist1d->Draw("");

      canvas->Update();

      CMS_lumi( canvas.get(), iPeriod, iPos );

      canvas->Update();
      canvas->RedrawAxis();
      canvas->GetFrame()->Draw();

      // canvas->SaveAs(dirname+"_"+strZA+"_1d.png");
    }
  }
}
