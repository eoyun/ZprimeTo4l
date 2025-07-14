#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

static double retrieveLumi(const std::string& anlyzrEra) {
  if (anlyzrEra=="20UL16APV")
    return 19.5;
  else if (anlyzrEra=="20UL16" || anlyzrEra=="")
    return 16.8;
  else if (anlyzrEra=="20UL17")
    return 41.48;
  else if (anlyzrEra=="20UL18")
    return 59.83;

  return 0.;
}

void estimateREFFfromTTree(TString era) {
  setTDRStyle();
  // gStyle->SetLineWidth(2);

  writeExtraText = true;       // if extra text
  extraText  = "Preliminary";  // default extra text is "Preliminary"

  static constexpr double WZxsec = 5.213; // 0.65*62.78; // 5.213
  static constexpr double ZZxsec = 13.81;
  static TString postfix = era;
  TString fname = era;

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
    postfix = "";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
  } else if (era=="run2") {
    lumi_sqrtS = "";
    lumi_13TeV = "137.6 fb^{-1}";
    postfix = "";
    fname = "20UL16";
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

  class RMFFera {
  public:
    RMFFera(const TString& aera)
    : era_(aera) {
      dataFile_ = new TFile("EleAnalyzer_"+aera+"_data.root","READ");
      WZFile_ = new TFile("EleAnalyzer_"+aera+"_WZFXFX.root","READ");
      ZZFile_ = new TFile("EleAnalyzer_"+aera+"_ZZ.root","READ");
    }

    ~RMFFera()=default;

    TH1D* rebinnedHisto(TH1D* ahist, std::vector<double>& binning) {
      const int nbin = binning.size()-1;
      TH1D* rebin = (TH1D*)ahist->Rebin(nbin,TString(ahist->GetName())+"_rebin",&(binning[0]));

      return rebin;
    }

    TH1D* subtractHist(const TH1D* ahist, const TH1D* bhist) {
      TH1D* result = (TH1D*)ahist->Clone();

      for (unsigned ibin = 0; ibin < ahist->GetNbinsX()+2; ibin++) {
        if (ahist->GetBinContent(ibin)==0.)
          continue;

        result->SetBinContent(ibin, std::max( ahist->GetBinContent(ibin) - bhist->GetBinContent(ibin), 0.) );
        result->SetBinError(ibin, result->GetBinContent(ibin)==0. ? 0. : std::hypot(ahist->GetBinError(ibin), bhist->GetBinError(ibin)) );
      }

      return result;
    }

    TH1D* fillHist(TFile* afile, const TString& treeName, std::vector<double>& binning, const double xsec=0.) {
      const TString aera = (era_=="20UL16") ? "" : era_;
      const TString anlyzrName = (xsec==0.) ? "resolvedEleCRanalyzerData" : "resolvedEleCRanalyzer"+aera;
      TTree* atree = (TTree*)afile->Get(anlyzrName+"/"+treeName);

      const int nbin = binning.size()-1;
      const TString title = treeName + era_ + TString(afile->GetName());
      TH1D* ahist = new TH1D(title,title,nbin,&(binning[0]));

      float pt = -1., dr = -1., wgt = 0., invM = -1., iso=-1.;
      int passMva = 0;
      float dEtaInSeed = 0., dPerpIn = 0.;
      atree->SetBranchAddress("pt",&pt);
      atree->SetBranchAddress("dr",&dr);
      atree->SetBranchAddress("wgt",&wgt);
      atree->SetBranchAddress("invM",&invM);

      if (treeName=="denomTree") {
        atree->SetBranchAddress("iso",&iso);
        atree->SetBranchAddress("passMva",&passMva);
        atree->SetBranchAddress("dEtaInSeed",&dEtaInSeed);
        atree->SetBranchAddress("dPerpIn",&dPerpIn);
      }

      for (int idx=0; idx < atree->GetEntries(); ++idx) {
        atree->GetEntry(idx);

        if ((fabs(dEtaInSeed) < 0.01 || fabs(dPerpIn) < 0.01))
          ahist->Fill(pt,wgt);
      }

      if (xsec > 0.) {
        const double sumwgt = ((TH1D*)afile->Get("evtCounter/h_sumW"))->GetBinContent(1);
        ahist->Scale(xsec*retrieveLumi(era_.Data())*1000./sumwgt);
      }

      return ahist;
    }

    void fill(std::vector<double>& binning) {
      dataNumerHist_ = fillHist(dataFile_,"numerTree",binning);
      WZNumerHist_ = fillHist(WZFile_,"numerTree",binning,WZxsec);
      ZZNumerHist_ = fillHist(ZZFile_,"numerTree",binning,ZZxsec);

      dataDenomHist_ = fillHist(dataFile_,"denomTree",binning);
      WZDenomHist_ = fillHist(WZFile_,"denomTree",binning,WZxsec);
      ZZDenomHist_ = fillHist(ZZFile_,"denomTree",binning,ZZxsec);
    }

    void add(const RMFFera& other) {
      dataNumerHist_->Add(other.dataNumerHist_);
      WZNumerHist_->Add(other.WZNumerHist_);
      ZZNumerHist_->Add(other.ZZNumerHist_);

      dataDenomHist_->Add(other.dataDenomHist_);
      WZDenomHist_->Add(other.WZDenomHist_);
      ZZDenomHist_->Add(other.ZZDenomHist_);
    }

    void load(std::vector<double>& xbins) {
      fill(xbins);
    }

    TH1D* run() {
      TH1D* numerHist = subtractHist( subtractHist(dataNumerHist_,WZNumerHist_), ZZNumerHist_ );
      TH1D* denomHist = subtractHist( subtractHist(dataDenomHist_,WZDenomHist_), ZZDenomHist_ );

      numerHist->Divide(denomHist);

      return numerHist;
    }

    const TString era_;

    TFile* dataFile_;
    TFile* WZFile_;
    TFile* ZZFile_;

    TH1D* dataNumerHist_;
    TH1D* WZNumerHist_;
    TH1D* ZZNumerHist_;

    TH1D* dataDenomHist_;
    TH1D* WZDenomHist_;
    TH1D* ZZDenomHist_;
  };

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

  std::vector<double> xbins = {0.,20.,25.,30.,40.,50.,70.,100.,150.,250.,500.};

  auto rmff1 = RMFFera(fname);
  auto rmff2 = RMFFera("20UL16APV");
  auto rmff3 = RMFFera("20UL17");
  auto rmff4 = RMFFera("20UL18");
  rmff1.load(xbins);
  rmff2.load(xbins);
  rmff3.load(xbins);
  rmff4.load(xbins);
  rmff1.add(rmff2);
  rmff1.add(rmff3);
  rmff1.add(rmff4);

  auto estimateCenter = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    out.push_back(vec.front());

    for (unsigned idx = 0; idx < vec.size()-1; idx++)
      out.push_back( (vec.at(idx) + vec.at(idx+1) ) / 2. );

    out.push_back(vec.back());

    return std::move(out);
  };

  auto estimateWidth = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;

    out.push_back(0.);

    for (unsigned idx = 1; idx < vec.size()-1; idx++)
      out.push_back( ( vec.at(idx+1) - vec.at(idx) ) / 2. );

    out.push_back(0.);

    return std::move(out);
  };

  TH1D* result = rmff1.run();
  result->SetMinimum(0.);
  result->GetXaxis()->SetTitle("p_{T}(e) [GeV]");
  result->Draw("E1");

  TF1* func = new TF1("REFFfunc","[0]+[1]*TMath::Exp(-[2]*x)",20.,500.);
  func->SetParameter(2,0.0001);
  func->SetParameter(0,0.01);
  func->SetLineColor(kRed);
  func->SetLineWidth(2);
  func->SetLineStyle(2);
  func->SetParLimits(0,0.,1.);
  TFitResultPtr fitResult = result->Fit(func,"RS");

  canvas_2->Update();
  TPaveStats* stats = (TPaveStats*)canvas_2->GetPrimitive("stats");
  stats->SetTextColor(kRed);
  stats->SetX1NDC(.65);
  stats->SetX2NDC(.95);
  stats->SetY1NDC(.7);
  stats->SetY2NDC(.9);

  const unsigned nbinSel = result->GetNbinsX();
  std::vector<double> xcen = estimateCenter(xbins);
  std::vector<double> xbinw = estimateWidth(xbins);
  double ci[nbinSel+2];
  fitResult->GetConfidenceIntervals(nbinSel+2,1,0,&(xcen[0]),ci,0.68,false);
  double ybin[nbinSel+2];

  for (unsigned idx = 0; idx < nbinSel+2; idx++) {
    ybin[idx] = func->Eval(xcen[idx]);
  }

  auto errGr = new TGraphErrors(nbinSel+2,&(xcen[0]),ybin,&(xbinw[0]),ci);
  errGr->SetFillColor(kRed);
  errGr->SetFillStyle(3004);
  errGr->Draw("3"); 

  SaveAs(canvas_2,"reff_run2.pdf");

  // return;

  // save file
  TH1D* outhist = (TH1D*)result->Clone("REFF_pt");
  const TString outname = "REFF_"+era+".root";
  TFile* outfile = new TFile(outname,"RECREATE");

  std::cout << "...writing " << outname << std::endl;
  outhist->Write();
  func->Write();
  outfile->Close();
}
