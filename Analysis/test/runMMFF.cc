#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runMMFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  static double valLumi = 0.;
  static TString postfix = era;

  if (era=="20UL16APV") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "19.5 fb^{-1}";
    valLumi = 19.5;
  } else if (era=="20UL16") {
    lumi_sqrtS = "2016 (13 TeV)";
    lumi_13TeV = "16.8 fb^{-1}";
    valLumi = 16.8;
    postfix = "";
  } else if (era=="20UL17") {
    lumi_sqrtS = "2017 (13 TeV)";
    lumi_13TeV = "41.48 fb^{-1}";
    valLumi = 41.48;
  } else if (era=="20UL18") {
    lumi_sqrtS = "2018 (13 TeV)";
    lumi_13TeV = "59.83 fb^{-1}";
    valLumi = 59.83;
  } else {
    std::cout << "check era..." << std::endl;
  }

  static TString anlyzrMC = "mergedMuCRanalyzer"+postfix;
  static TString anlyzrData = "mergedMuCRanalyzerData";

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

  TFile* datafile = new TFile("MMCR_"+era+"_data.root","READ");

  //TFile* H250A1file = new TFile("MMCR_"+era+"_H250A1.root","READ");
  //TFile* H750A1file = new TFile("MMCR_"+era+"_H750A1.root","READ");
  //TFile* H2000A1file = new TFile("MMCR_"+era+"_H2000A1.root","READ");

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  auto H750A1sample = SigSample(new TFile("MMCR_"+era+"_H750A1.root","READ"),"H750A1");
  std::vector<SigSample> sigsamples = {H750A1sample};

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

  auto drawRatio = [] (const TH1D* numer, const TH1D* denom, TPad* pad) {
    TH1D* ratio = (TH1D*)numer->Clone();
    ratio->SetStats(0);
    ratio->SetTitle("");
    ratio->Divide(denom);
    ratio->GetYaxis()->SetTitle("Obs/Exp");
    ratio->GetYaxis()->SetTitleSize(0.1);
    ratio->GetYaxis()->SetTitleOffset(0.4);
    ratio->GetXaxis()->SetLabelSize(0.1);
    ratio->GetYaxis()->SetLabelSize(0.1);
    ratio->GetXaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetLabelOffset(0.01);
    ratio->GetYaxis()->SetRangeUser(0.,2.);
    ratio->GetXaxis()->SetTitleSize(0.12);
    ratio->GetXaxis()->SetTitleOffset(0.75);
    ratio->SetLineColor(kBlack);

    pad->cd();
    ratio->Draw("E1");
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

  auto density = [] (TH1D* ahist) {
    for (unsigned ibin = 0; ibin < ahist->GetNbinsX(); ibin++) {
      ahist->SetBinContent(ibin, ahist->GetBinContent(ibin)/ahist->GetBinWidth(ibin));
      ahist->SetBinError(ibin, ahist->GetBinError(ibin)/ahist->GetBinWidth(ibin));
    }
  };

  TH1D* nm = (TH1D*)datafile->Get("mergedMuCRanalyzerData/3M_antiRpt_mt")->Clone();
  TH1D* dm = (TH1D*)datafile->Get("mergedMuCRanalyzerData/3M_antiDphi_antiRpt_mt")->Clone();

  nm->Rebin( nm->GetNbinsX() );
  dm->Rebin( dm->GetNbinsX() );

  nm->Divide(nm,dm,1,1,"B");
  const double ff = nm->GetBinContent(1);
  const double ffErr = dm->GetBinError(1);

  class HistLoaderBase {
  protected:
    TFile* datafile_ = nullptr;
    TFile* WZfile_ = nullptr;
    TFile* ZZfile_ = nullptr;

    TFile* datacard_ = nullptr;
    TDirectory* dir_ = nullptr;

  public:
    HistLoaderBase(TFile* adatafile, TFile* aWZfile, TFile* aZZfile) {
      datafile_ = adatafile;
      WZfile_ = aWZfile;
      ZZfile_ = aZZfile;
    }

    HistLoaderBase(TFile* adatafile) {
      datafile_ = adatafile;
    }

    ~HistLoaderBase()=default;

    void preparecard(TString name, TString dirname) {
      datacard_ = new TFile(name,"RECREATE");
      dir_ = datacard_->mkdir(dirname);
      dir_->cd();
    }

    void close() {
      datacard_->Close();
      datacard_ = nullptr;
      dir_ = nullptr;
    }
  };

  class HistLoaderMM : public HistLoaderBase {
  public:
    HistLoaderMM(TFile* adatafile, std::vector<SigSample> sigsamples)
    : HistLoaderBase(adatafile,nullptr,nullptr), sigFiles_(sigsamples) {}

    ~HistLoaderMM()=default;

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* FFHist_ = nullptr;
    TH1D* FFHist_up_ = nullptr;
    TH1D* FFHist_dn_ = nullptr;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;

  public:
    void load(TString numName, TString denomName) {
      if (dataHist_) {
        delete dataHist_, FFHist_;
        dataHist_ = nullptr;
      }

      if (FFHist_up_) {
        delete FFHist_up_, FFHist_dn_;
        FFHist_up_ = nullptr;
      }

      dataHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+numName)->Clone();
      FFHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+denomName)->Clone();
      FFHist_->SetFillColor(kGray);
      FFHist_->SetLineWidth(0);

      sigHist_.clear();

      if (numName.Contains("3M_mt")) {
        const unsigned isigDiv = 3;

        auto retrieveSigHist = [this,&numName] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( anlyzrMC+"/"+numName+systName )->Clone();
          ahist->Scale( valLumi*1000.*0.001 / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned isig = 0; isig < sigFiles_.size(); isig++) {
          sigHist_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"") );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed-3);

          if (isig > isigDiv-1)
            sigHist_.back()->SetLineColor(kBlue-3);
        }

        FFHist_up_ = (TH1D*)datafile_->Get(anlyzrData+"/"+denomName+"_up")->Clone();
        FFHist_dn_ = (TH1D*)datafile_->Get(anlyzrData+"/"+denomName+"_dn")->Clone();
      }
    }

    void compare(TPad* pad, int rebin=1, double scale=1.) {
      if ( dataHist_->GetNbinsX()!=FFHist_->GetNbinsX() ) {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (FFHist_up_) {
          FFHist_up_->Rebin( FFHist_up_->GetNbinsX()/dataHist_->GetNbinsX() );
          FFHist_dn_->Rebin( FFHist_dn_->GetNbinsX()/dataHist_->GetNbinsX() );
        }
      }

      if (rebin!=1) {
        dataHist_->Rebin(rebin);
        FFHist_->Rebin(rebin);

        if (FFHist_up_) {
          FFHist_up_->Rebin(rebin);
          FFHist_dn_->Rebin(rebin);
        }

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);
        }
      }

      if (scale!=1.) {
        FFHist_->Scale(scale);
        FFHist_up_->Scale(scale);
        FFHist_dn_->Scale(scale);
      }

      pad->cd();

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMaximum(1.5*dataHist_->GetMaximum());

      if (TString(dataHist_->GetName()).Contains("mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,2500.);
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      dataHist_->Draw("E1");
      FFHist_->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("3M_mt")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (TString(dataHist_->GetName()).Contains("3M_mt")) {
        TLegend* legend = new TLegend(0.63,0.77,0.93,0.9);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(FFHist_,"Data-driven bkg");
        legend->AddEntry(sigHist_.at(0),"M_{A} = 1 GeV (#sigma = 1 fb)");

        legend->Draw();
      }

      if (FFHist_up_) {
        std::vector<double> x0, y0, errx, erryDn, erryUp;

        for (unsigned idx = 1; idx <= FFHist_->GetNbinsX(); idx++) {
          x0.push_back(FFHist_->GetBinCenter(idx));
          y0.push_back(FFHist_->GetBinContent(idx));
          errx.push_back(FFHist_->GetBinWidth(idx)/2.);

          double valFFup = FFHist_up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
          double valFFdn = FFHist_->GetBinContent(idx) - FFHist_dn_->GetBinContent(idx);

          erryUp.push_back( valFFup );
          erryDn.push_back( valFFdn );
        }

        auto gr = new TGraphAsymmErrors(FFHist_->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kGray+2);
        gr->SetLineColor(kGray+2);
        gr->SetFillStyle(3004);
        gr->Draw("2");
      }

      if (dir_) {
        dir_->WriteTObject(dataHist_,"data_obs");
        dir_->WriteTObject(FFHist_,"MM");

        dir_->WriteTObject(FFHist_up_,"MM_mergedMuFakeFactorUp");
        dir_->WriteTObject(FFHist_dn_,"MM_mergedMuFakeFactorDown");

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
        }
      }
    }
  }; // class

  auto aloaderMM = HistLoaderMM(datafile,sigsamples);

  aloaderMM.load("3M_mt","3M_antiDphi_antiRpt_mt_xFF"); // 3M_antiRpt_mt_xFF
  aloaderMM.preparecard("MMFF_"+era+"_datacard.root","mergedMu3M");
  aloaderMM.compare(canvas_1,10,ff);
  SaveAs(canvas_1,"FF_3M_mt.png");
  aloaderMM.close();

  aloaderMM.load("3M_CRdphi_mt","3M_CRdphi_antiRpt_mt_xFF");
  aloaderMM.compare(canvas_1,10);
  SaveAs(canvas_1,"FF_3M_CRdphi_mt.png");

  aloaderMM.load("3M_antiDphi_mt","3M_antiDphi_antiRpt_mt_xFF");
  aloaderMM.compare(canvas_1,5);
  SaveAs(canvas_1,"FF_3M_antiDphi_mt.png");

  aloaderMM.load("3M_antiDphi_MET_pt","3M_antiDphi_antiRpt_MET_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_antiDphi_MET_pt.png");

  aloaderMM.load("3M_CRdphi_MM_pt","3M_CRdphi_antiRpt_MM_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_CRdphi_MET_pt.png");

  aloaderMM.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_antiDphi_MM_pt.png");

  aloaderMM.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_CRdphi_MM_pt.png");

  return;
}
