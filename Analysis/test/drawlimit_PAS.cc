#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
// #include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/tdrstyle.C"
#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/CMS_lumi.C"

void drawlimit_PAS(TString isZA) {
  setTDRStyle();
  // gStyle->SetLineWidth(2);
  TString era = "run2";

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"

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
  int iPos = 0;

  if( iPos==0 )
    relPosX = 0.11;

  int W = 1000;
  int H = 1200;

  int H_ref = 1200;
  int W_ref = 1000;

  // references for T, B, L, R
  float T = 0.05*H_ref;
  float B = 0.08*H_ref;
  float L = 0.1*W_ref;
  float R = 0.04*W_ref;

  auto* canvas1 = new TCanvas("canvas","canvas",W,H,W,H);
  canvas1->SetFillColor(0);
  canvas1->SetBorderMode(0);
  canvas1->SetFrameFillStyle(0);
  canvas1->SetFrameBorderMode(0);
  canvas1->SetLeftMargin( L/W+0.01 );
  canvas1->SetRightMargin( R/W );
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
      //canvas->RedrawAxis();
      //canvas->GetFrame()->Draw();
    }

    canvas->SaveAs(name.c_str());
  };

  auto retrieve = [&isZA] (const TString& atype, const TString& hmass, const TString& amass, const TString& quant) -> double {
    TString filename = TString("./")+isZA+amass+"/"+hmass+"/higgsCombine.X"+hmass+isZA+amass+".HybridNew.mH"+hmass;

    if (quant!="")
      filename += ".quant"+quant;

    filename += ".root";

    //if (atype!="")
    //  filename = TString("higgsCombine.")+atype+"_X"+hmass+"Y"+amass+".HybridNew.mH"+hmass+".quant"+quant+".root";

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

  auto graphNominal = [&retrieve,&massvec] (const TString& atype, const TString& amass, const TString& obs) -> std::unique_ptr<TGraph> {
    std::vector<double> limits, massdouble;

    for (unsigned idx=0; idx<massvec.size(); idx++) {
      double val = retrieve(atype,std::to_string(massvec.at(idx)),amass,obs);

      if (val < 0.)
        continue;

      limits.push_back(val);
      massdouble.push_back(static_cast<double>(massvec.at(idx)));
    }

    auto gr = std::make_unique<TGraph>(limits.size(),&(massdouble[0]),&(limits[0]));
    return std::move(gr);
  };

  auto graph1sig = [&retrieve,&massvec] (const TString& amass) -> std::unique_ptr<TGraphAsymmErrors> {
    std::vector<double> x0, y0, errydn1sig, erryup1sig, errx;

    for (unsigned idx=0; idx<massvec.size(); idx++) {
      double nominal = retrieve("",std::to_string(massvec.at(idx)),amass,"0.500");
      double quant0p16 = retrieve("",std::to_string(massvec.at(idx)),amass,"0.160");
      double quant0p84 = retrieve("",std::to_string(massvec.at(idx)),amass,"0.840");

      if (nominal < 0. || quant0p16 < 0. || quant0p84 < 0.)
        continue;

      y0.push_back(nominal);
      errydn1sig.push_back(nominal-quant0p16);
      erryup1sig.push_back(quant0p84-nominal);
      errx.push_back(0.);
      x0.push_back(static_cast<double>(massvec.at(idx)));
    }

    y0.push_back(y0.back());
    errydn1sig.push_back(errydn1sig.back());
    erryup1sig.push_back(erryup1sig.back());
    errx.push_back(errx.back());
    x0.push_back(2500.);

    auto gr = std::make_unique<TGraphAsymmErrors>(x0.size(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(errydn1sig[0]),&(erryup1sig[0]));

    return std::move(gr);    
  };

  auto graph2sig = [&retrieve,&massvec] (const TString& amass) -> std::unique_ptr<TGraphAsymmErrors> {
    std::vector<double> x0, y0, errydn2sig, erryup2sig, errx;

    for (unsigned idx=0; idx<massvec.size(); idx++) {
      double nominal = retrieve("",std::to_string(massvec.at(idx)),amass,"0.500");
      double quant0p025 = retrieve("",std::to_string(massvec.at(idx)),amass,"0.025");
      double quant0p975 = retrieve("",std::to_string(massvec.at(idx)),amass,"0.975");

      if (nominal < 0. || quant0p025 < 0. || quant0p975 < 0.)
        continue;

      y0.push_back(nominal);
      errydn2sig.push_back(nominal-quant0p025);
      erryup2sig.push_back(quant0p975-nominal);
      errx.push_back(0.);
      x0.push_back(static_cast<double>(massvec.at(idx)));
    }

    y0.push_back(y0.back());
    errydn2sig.push_back(errydn2sig.back());
    erryup2sig.push_back(erryup2sig.back());
    errx.push_back(errx.back());
    x0.push_back(2500.);

    auto gr = std::make_unique<TGraphAsymmErrors>(x0.size(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(errydn2sig[0]),&(erryup2sig[0]));

    return std::move(gr);    
  };

  // Define the rows and columns for subpads
  int rows = 4;
  int cols = 3;

  // Create subpads and divide the canvas

  double leftmargin = L/W;
  double rightmargin = R/W;
  double topmargin = T/H;
  double bottommargin = B/H;

  std::vector<TString> amassvec = {"0p4","0p6","0p8","1","1p5","2","5","10","50","100"};
  int seq = 0;

  TGraph* grobs = nullptr;
  TGraph* grexp = nullptr;
  TGraphAsymmErrors* gr1sig = nullptr;
  TGraphAsymmErrors* gr2sig = nullptr;

  for (int i = 0; i < cols; ++i) {
    for (int j = 0; j < rows; ++j) {
      double xmin = leftmargin + static_cast<double>(i)*(1.-leftmargin-rightmargin)/cols;
      double xmax = leftmargin + static_cast<double>(i+1)*(1.-leftmargin-rightmargin)/cols;
      double ymin = bottommargin + static_cast<double>(rows-j-1)*(1.-topmargin-bottommargin)/rows;
      double ymax = bottommargin + static_cast<double>(rows-j)*(1.-topmargin-bottommargin)/rows;

      if (i==0)
        xmin -= 0.3*leftmargin;
      if (j==rows-1)
        ymin -= 0.3*bottommargin;

      if (j==0 && (i==cols-1 || i==cols-2))
        continue;

      // Create subpad
      TPad *pad = new TPad(Form("pad%d", i * cols + j + 1), Form("Subpad %d", i * cols + j + 1), xmin, ymin, xmax, ymax);

      pad->SetLeftMargin( 0. );
      pad->SetRightMargin( 0.0 );
      pad->SetTopMargin( 0.0 );
      pad->SetBottomMargin( 0.0 );
      pad->SetFillColor(0);
      pad->SetBorderMode(0);
      pad->SetFrameFillStyle(0);
      pad->SetFrameBorderMode(0);

      if (i==0)
        pad->SetLeftMargin(1.3*leftmargin);
      //if (i==cols-1)
      //  pad->SetRightMargin(rightmargin);

      if (j==rows-1)
        pad->SetBottomMargin(1.1*bottommargin);
      //if (j==0)
      //  pad->SetTopMargin(topmargin);

      pad->Draw();
      pad->cd();

      auto gr1 = graph1sig(amassvec.at(seq));
      auto gr2 = graph2sig(amassvec.at(seq));
      auto grNominal = graphNominal("",amassvec.at(seq),"");
      auto grNominalExp = graphNominal("",amassvec.at(seq),"0.500");

      pad->SetLogy();

      TString Ychar = isZA=="A" ? "Y" : "Z";

      gr2->SetFillColor(kOrange);
      gr2->SetLineColor(kOrange);
      gr2->SetMarkerColor(kOrange);
      gr2->SetMarkerSize(0);
      //gr2->GetXaxis()->SetTitle("M_{X} [GeV]");
      //gr2->GetXaxis()->SetTitleSize(0.05);
      //gr2->GetYaxis()->SetTitle("#sigma(pp#rightarrowX) #times B(X#rightarrow"+Ychar+"Y#rightarrow4l) [fb]");
      //gr2->GetYaxis()->SetTitleOffset(1.1);
      //gr2->GetXaxis()->SetTitleOffset(1.0);
      //gr2->GetYaxis()->SetTitleSize(0.04);
      gr2->GetYaxis()->SetLabelSize(0.09);
      gr2->GetXaxis()->SetLabelSize(0.1);

      gr2->SetMaximum(20.);
      gr2->SetMinimum(0.02);

      TAxis *axis = gr2->GetXaxis();
      axis->SetLimits(250.,1999.9);
      axis->SetNdivisions(504);

      gr1->SetFillColor(kGreen+1);
      gr1->SetLineColor(kGreen+1);
      gr1->SetMarkerColor(kGreen+1);
      gr1->SetMarkerSize(0);
      gr2->Draw("a4"); // a3
      gr1->Draw("same4"); // samea3l

      grNominal->SetLineColor(kBlack);
      grNominal->SetMarkerSize(0);
      grNominal->SetLineWidth(2);
      grNominal->Draw("sameC");
      grNominalExp->SetLineStyle(kDashed);
      grNominalExp->SetMarkerSize(0);
      grNominalExp->SetLineWidth(2);
      grNominalExp->SetLineColor(kBlack);
      grNominalExp->Draw("sameC");

      gr1sig = gr1.release();
      gr2sig = gr2.release();
      grobs = grNominal.release();
      grexp = grNominalExp.release();

      pad->RedrawAxis();

      TString amassDecimal = amassvec.at(seq);
      amassDecimal.ReplaceAll("p",".");

      auto atext = std::make_unique<TPaveText>(0.55,0.75,0.95,0.95,"NDC");
      atext->SetBorderSize(0);
      atext->SetFillColor(0);
      atext->SetFillStyle(0);
      TString astr = "M_{Y} = "+amassDecimal+" GeV";
      atext->AddText(astr);
      ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
      // ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(12);
      atext->Draw();
      atext.release();

      canvas1->cd();  // Go back to the main canvas

      seq++;
    }
  }

  canvas1->cd();

  auto legend = std::make_unique<TLegend>(0.54,0.8,0.99,0.9);
  legend->SetNColumns(2);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(grobs,"Observed");
  legend->AddEntry(grexp,"Median expected");
  //legend->AddEntry(grMerged.get(),"merged median exp.");
  //legend->AddEntry(grResolved.get(),"resolved median exp.");
  legend->AddEntry(gr1sig,"68% expected");
  legend->AddEntry(gr2sig,"95% expected");
  legend->Draw();

  auto btext = std::make_unique<TPaveText>(0.68,0.9,0.95,0.95,"NDC");
  btext->SetBorderSize(0);
  btext->SetFillColor(0);
  btext->SetFillStyle(0);
  TString bstr = "95% CL upper limits";
  btext->AddText(bstr);
  ((TText*)btext->GetListOfLines()->Last())->SetTextColor(kBlack);
  btext->Draw();

  auto ytext = std::make_unique<TPaveText>(0.0,0.4,0.4,0.92,"NDC");
  ytext->SetBorderSize(0);
  ytext->SetFillColor(0);
  ytext->SetFillStyle(0);
  TString Ychar = isZA=="A" ? "Y" : "Z";
  TString astr = "#sigma(pp#rightarrowX) #times B(X#rightarrow"+Ychar+"Y#rightarrow4l) [fb]";
  ytext->AddText(astr);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextFont(42);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAngle(90);
  ((TText*)ytext->GetListOfLines()->Last())->SetTextAlign(13);
  ytext->Draw();

  auto xtext = std::make_unique<TPaveText>(0.8,0.01,0.97,0.055,"NDC");
  xtext->SetBorderSize(0);
  xtext->SetFillColor(0);
  xtext->SetFillStyle(0);
  TString xstr = "M_{X} [GeV]";
  xtext->AddText(xstr);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextColor(kBlack);
  ((TText*)xtext->GetListOfLines()->Last())->SetTextFont(42);
  xtext->Draw();

  canvas1->Update();  // Refresh canvas

  SaveAs(canvas1,std::string((TString("limit_")+isZA+"_PAS_"+era+".pdf").Data()));

  return;







  /*auto gr1 = graph1sig();
  auto gr2 = graph2sig();
  //auto grMerged = graphNominal("merged");
  //auto grResolved = graphNominal("resolved");
  auto grNominal = graphNominal("","");
  auto grNominalExp = graphNominal("","0.500");

//  canvas1->SetLogx();
  canvas1->SetLogy();

  TString Ychar = isZA=="A" ? "Y" : "Z";

  gr2->SetFillColor(kOrange);
  gr2->SetLineColor(kOrange);
  gr2->SetMarkerColor(kOrange);
  gr2->SetMarkerSize(0);
  gr2->GetXaxis()->SetTitle("M_{X} [GeV]");
  gr2->GetXaxis()->SetTitleSize(0.05);
  gr2->GetYaxis()->SetTitle("#sigma(pp#rightarrowX) #times B(X#rightarrow"+Ychar+"Y#rightarrow4l) [fb]");
  gr2->GetYaxis()->SetTitleOffset(1.1);
  gr2->GetXaxis()->SetTitleOffset(1.0);
  gr2->GetYaxis()->SetTitleSize(0.04);
  gr2->GetYaxis()->SetLabelSize(0.04);
  gr2->GetXaxis()->SetLabelSize(0.04);

  gr2->SetMaximum(20.);
  gr2->SetMinimum(0.01);

  TAxis *axis = gr2->GetXaxis();
  axis->SetLimits(250.,2000.);

  gr1->SetFillColor(kGreen+1);
  gr1->SetLineColor(kGreen+1);
  gr1->SetMarkerColor(kGreen+1);
  gr1->SetMarkerSize(0);
  gr2->Draw("a4"); // a3
  gr1->Draw("same4"); // samea3l

  grNominal->SetLineColor(kBlack);
  grNominal->SetMarkerSize(0);
  grNominal->SetLineWidth(2);
  grNominal->Draw("sameC");
  grNominalExp->SetLineStyle(kDashed);
  grNominalExp->SetMarkerSize(0);
  grNominalExp->SetLineWidth(2);
  grNominalExp->SetLineColor(kBlack);
  grNominalExp->Draw("sameC");

  auto legend = std::make_unique<TLegend>(0.6,0.62,0.92,0.85);
  legend->SetBorderSize(0);
  legend->SetFillStyle(0);
  legend->AddEntry(grNominal.get(),"Median observed");
  legend->AddEntry(grNominalExp.get(),"Median expected");
  //legend->AddEntry(grMerged.get(),"merged median exp.");
  //legend->AddEntry(grResolved.get(),"resolved median exp.");
  legend->AddEntry(gr1.get(),"68% expected");
  legend->AddEntry(gr2.get(),"95% expected");
  legend->Draw();

  auto btext = std::make_unique<TPaveText>(0.6,0.85,0.92,0.9,"NDC");
  btext->SetBorderSize(0);
  btext->SetFillColor(0);
  btext->SetFillStyle(0);
  TString bstr = "95% CL upper limits";
  btext->AddText(bstr);
  ((TText*)btext->GetListOfLines()->Last())->SetTextColor(kBlack);
  btext->Draw();

  TString amassDecimal = amass;
  amassDecimal.ReplaceAll("p",".");

  auto atext = std::make_unique<TPaveText>(0.68,0.15,0.92,0.2,"NDC");
  atext->SetBorderSize(0);
  atext->SetFillColor(0);
  atext->SetFillStyle(0);
  TString astr = "M_{Y} = "+amassDecimal+" GeV";
  atext->AddText(astr);
  ((TText*)atext->GetListOfLines()->Last())->SetTextColor(kBlack);
  // ((TText*)atext->GetListOfLines()->Last())->SetTextAlign(12);
  atext->Draw();

  canvas1->cd();

  SaveAs(canvas1,std::string((TString("limit_")+isZA+"_PAS_"+era+".pdf").Data()));

  return;*/
}
