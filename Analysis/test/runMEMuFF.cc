#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runMEMuFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static double valLumi = 0.;
  static constexpr double WZxsec_ = 0.65*62.78;
  static constexpr double ZZxsec_ = 13.81;
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

  static TString anlyzrMC = "mergedEMuCRanalyzer"+postfix;
  static TString anlyzrData = "mergedEMuCRanalyzerData";

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

  TFile* datafile = new TFile("MuAnalyzer_"+era+"_data.root","READ");
  TFile* WZfile = new TFile("MuAnalyzer_"+era+"_WZ.root","READ");
  TFile* ZZfile = new TFile("MuAnalyzer_"+era+"_ZZ.root","READ");

  //TFile* H250A1file = new TFile("MEMuCR_"+era+"_H250A1.root","READ");
  //TFile* H750A1file = new TFile("MEMuCR_"+era+"_H750A1.root","READ");
  //TFile* H2000A1file = new TFile("MEMuCR_"+era+"_H2000A1.root","READ");
  //TFile* H250A10file = new TFile("MEMuCR_"+era+"_H250A10.root","READ");
  //TFile* H750A10file = new TFile("MEMuCR_"+era+"_H750A10.root","READ");
  //TFile* H2000A10file = new TFile("MEMuCR_"+era+"_H2000A10.root","READ");

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  auto H750A1sample = SigSample(new TFile("MuAnalyzer_"+era+"_H750A1.root","READ"),"H750A1");
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

  TH1D* nm = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_FF_CRME_MET_pt")->Clone();
  TH1D* dm = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_FF_CRME_antiDphi_MET_pt")->Clone();

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
    std::vector<TH1D*> sigHist_modHeepIdUp_;
    std::vector<TH1D*> sigHist_modHeepIdDn_;
    std::vector<TH1D*> sigHist_mergedEleIdUp_;
    std::vector<TH1D*> sigHist_mergedEleIdDn_;

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
      sigHist_modHeepIdUp_.clear();
      sigHist_modHeepIdDn_.clear();
      sigHist_mergedEleIdUp_.clear();
      sigHist_mergedEleIdDn_.clear();

      if (numName.Contains("CRME_mt") || numName.Contains("CRME_CRdphi_mt")) {
        FFHist_up_ = (TH1D*)datafile_->Get(anlyzrData+"/"+denomName+"_up")->Clone();
        FFHist_dn_ = (TH1D*)datafile_->Get(anlyzrData+"/"+denomName+"_dn")->Clone();
      }

      if (numName.Contains("CRME_mt")) {
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

          sigHist_modHeepIdUp_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_heepIdUp") );
          sigHist_modHeepIdDn_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_heepIdDn") );
          sigHist_mergedEleIdUp_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_mergedEleIdUp") );
          sigHist_mergedEleIdDn_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_mergedEleIdDn") );

          if (isig > isigDiv-1)
            sigHist_.back()->SetLineColor(kBlue-3);
        }
      }
    }

    void compare(TPad* pad, int rebin=1) {
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
          sigHist_modHeepIdUp_.at(idx)->Rebin(rebin);
          sigHist_modHeepIdDn_.at(idx)->Rebin(rebin);
          sigHist_mergedEleIdUp_.at(idx)->Rebin(rebin);
          sigHist_mergedEleIdDn_.at(idx)->Rebin(rebin);
        }
      }

      pad->cd();

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMaximum(1.5*dataHist_->GetMaximum());

      if (TString(dataHist_->GetName()).Contains("CRME_mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,2500.);
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      if (TString(dataHist_->GetName()).Contains("CRME_CRdphi_mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,1500.);
      }

      dataHist_->Draw("E1");
      FFHist_->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("CRME_mt")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (TString(dataHist_->GetName()).Contains("CRME_mt")) {
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
          dir_->WriteTObject(sigHist_modHeepIdUp_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdUp");
          dir_->WriteTObject(sigHist_modHeepIdDn_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdDown");
          dir_->WriteTObject(sigHist_mergedEleIdUp_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdUp");
          dir_->WriteTObject(sigHist_mergedEleIdDn_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdDown");
        }
      }
    }
  }; // class

  auto aloaderMM = HistLoaderMM(datafile,sigsamples);

  aloaderMM.load("1M_CRME_mt","1M_CRME_antiRpt_mt_xFF");
  aloaderMM.preparecard("MEMu1M_"+era+"_datacard.root","mergedEMu1M");
  aloaderMM.compare(canvas_1,5);
  SaveAs(canvas_1,"FF_1M_CRME_mt.png");
  aloaderMM.close();

  aloaderMM.load("1M_CRME_antiDphi_mt","1M_CRME_antiDphi_antiRpt_mt_xFF");
  aloaderMM.compare(canvas_1,5);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_mt.png");

  aloaderMM.load("1M_CRME_CRdphi_mt","1M_CRME_CRdphi_antiRpt_mt_xFF");
  aloaderMM.compare(canvas_1,5);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_mt.png");

  aloaderMM.load("1M_CRME_antiDphi_MET_pt","1M_CRME_antiDphi_antiRpt_MET_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_MET_pt.png");

  aloaderMM.load("1M_CRME_CRdphi_MET_pt","1M_CRME_CRdphi_antiRpt_MET_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_MET_pt.png");

  aloaderMM.load("1M_CRME_antiDphi_MM_pt","1M_CRME_antiDphi_antiRpt_MM_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_MM_pt.png");

  aloaderMM.load("1M_CRME_CRdphi_MM_pt","1M_CRME_CRdphi_antiRpt_MM_pt_xFF");
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_MM_pt.png");

  aloaderMM.load("1M_CRME_antiDphi_all_Et","1M_CRME_antiDphi_antiRpt_all_Et_xFF");
  aloaderMM.compare(canvas_1);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_all_Et.png");

  aloaderMM.load("1M_CRME_CRdphi_all_Et","1M_CRME_CRdphi_antiRpt_all_Et_xFF");
  aloaderMM.compare(canvas_1);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_all_Et.png");

  class HistLoader3E : public HistLoaderBase {
  public:
    HistLoader3E(TFile* adatafile, TFile* aWZfile, TFile* aZZfile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile,aWZfile,aZZfile),
      sigFiles_(sigFiles) {}

    ~HistLoader3E() {
      delete dataHist_;
    }

    class denomHists {
    public:
      denomHists() {}
      ~denomHists()=default;

    private:
      TH1D* SSFFHist_ = nullptr;
      TH1D* SSWZHist_ = nullptr;
      TH1D* SSZZHist_ = nullptr;
      TH1D* OSWZHist_ = nullptr;
      TH1D* OSZZHist_ = nullptr;

    public:
      void load(TString denomNameSS, TString denomNameOS, TFile* datafile, TFile* WZfile, TFile* ZZfile) {
        SSFFHist_ = (TH1D*)datafile->Get(anlyzrData+"/"+denomNameSS)->Clone();

        const double WZsumwgt = ((TH1D*)WZfile->Get("evtCounter/h_sumW"))->GetBinContent(1);
        const double ZZsumwgt = ((TH1D*)ZZfile->Get("evtCounter/h_sumW"))->GetBinContent(1);

        SSWZHist_ = (TH1D*)WZfile->Get(anlyzrMC+"/"+denomNameSS)->Clone();
        SSWZHist_->Scale( valLumi*1000.*WZxsec_/WZsumwgt );
        SSZZHist_ = (TH1D*)ZZfile->Get(anlyzrMC+"/"+denomNameSS)->Clone();
        SSZZHist_->Scale( valLumi*1000.*ZZxsec_/ZZsumwgt );
        OSWZHist_ = (TH1D*)WZfile->Get(anlyzrMC+"/"+denomNameOS)->Clone();
        OSWZHist_->Scale( valLumi*1000.*WZxsec_/WZsumwgt );
        OSZZHist_ = (TH1D*)ZZfile->Get(anlyzrMC+"/"+denomNameOS)->Clone();
        OSZZHist_->Scale( valLumi*1000.*ZZxsec_/ZZsumwgt );
      }

      void rebin(int r) {
        SSFFHist_->Rebin(r);
        SSWZHist_->Rebin(r);
        SSZZHist_->Rebin(r);
        OSWZHist_->Rebin(r);
        OSZZHist_->Rebin(r);
      }

      TH1D* SSFFHist() { return SSFFHist_; }

      std::unique_ptr<TH1D> subtractHist(const TH1D* denom, const TH1D* denom_prompt) {
        auto cloned = std::unique_ptr<TH1D>((TH1D*)denom->Clone());

        for (unsigned idx = 0; idx < cloned->GetNbinsX()+2; idx++) {
          cloned->SetBinContent(idx,std::max(cloned->GetBinContent(idx)-denom_prompt->GetBinContent(idx),0.));
          cloned->SetBinError(idx,std::hypot(cloned->GetBinError(idx),denom_prompt->GetBinError(idx)));
        }

        return std::move(cloned);
      };

      std::unique_ptr<TH1D> returnSS() {
        auto tempSubtract = subtractHist( SSFFHist_, SSWZHist_ );
        auto denomSSfinal = subtractHist( tempSubtract.get(), SSZZHist_ );
        denomSSfinal->SetFillColor(kCyan+1);

        return std::move(denomSSfinal);
      }

      std::unique_ptr<TH1D> returnOS() {
        auto denomOSfinal = std::unique_ptr<TH1D>((TH1D*)OSWZHist_->Clone());
        denomOSfinal->Add(OSZZHist_);
        denomOSfinal->SetFillColor(kOrange);

        return std::move(denomOSfinal);
      }

      std::unique_ptr<TH1D> returnAdded() {
        auto ss = returnSS();
        auto os = returnOS();
        ss->Add( os.get() );

        return std::move(ss);
      }
    }; // class denomHists

  private:
    TH1D* dataHist_ = nullptr;
    denomHists nominal_;
    denomHists SSFFup_;
    denomHists SSFFdn_;
    denomHists OSFFup_;
    denomHists OSFFdn_;
    denomHists heepIdUp_;
    denomHists heepIdDn_;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<TH1D*> sigHist_modHeepIdUp_;
    std::vector<TH1D*> sigHist_modHeepIdDn_;
    std::vector<TH1D*> sigHist_mergedEleIdUp_;
    std::vector<TH1D*> sigHist_mergedEleIdDn_;

  public:
    void load(TString numName, TString denomName) {
      if (dataHist_) {
        delete dataHist_;
      }

      dataHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+numName)->Clone();
      nominal_.load(denomName+"_xSSFF",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_);
      SSFFup_.load(denomName+"_xSSFF_up",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_);
      SSFFdn_.load(denomName+"_xSSFF_dn",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_);
      OSFFup_.load(denomName+"_xSSFF",denomName+"_xOSFF_up",datafile_,WZfile_,ZZfile_);
      OSFFdn_.load(denomName+"_xSSFF",denomName+"_xOSFF_dn",datafile_,WZfile_,ZZfile_);

      sigHist_.clear();
      sigHist_modHeepIdUp_.clear();
      sigHist_modHeepIdDn_.clear();
      sigHist_mergedEleIdUp_.clear();
      sigHist_mergedEleIdDn_.clear();

      if (numName.Contains("invM")) {
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

          sigHist_modHeepIdUp_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_heepIdUp") );
          sigHist_modHeepIdDn_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_heepIdDn") );
          sigHist_mergedEleIdUp_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_mergedEleIdUp") );
          sigHist_mergedEleIdDn_.push_back( retrieveSigHist(sigFiles_.at(isig).file_,"_mergedEleIdDn") );

          if (isig > isigDiv-1)
            sigHist_.back()->SetLineColor(kBlue-3);
        }

        heepIdUp_.load(denomName+"_xSSFF_heepIdUp",denomName+"_xOSFF_heepIdUp",datafile_,WZfile_,ZZfile_);
        heepIdDn_.load(denomName+"_xSSFF_heepIdDn",denomName+"_xOSFF_heepIdDn",datafile_,WZfile_,ZZfile_);
      }
    }

    void compare(TPad* pad, int rebin=1) {
      if (rebin!=1) {
        dataHist_->Rebin(rebin);

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);
          sigHist_modHeepIdUp_.at(idx)->Rebin(rebin);
          sigHist_modHeepIdDn_.at(idx)->Rebin(rebin);
          sigHist_mergedEleIdUp_.at(idx)->Rebin(rebin);
          sigHist_mergedEleIdDn_.at(idx)->Rebin(rebin);
        }
      }

      nominal_.rebin(nominal_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      SSFFup_.rebin(SSFFup_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      SSFFdn_.rebin(SSFFdn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      OSFFup_.rebin(OSFFup_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      OSFFdn_.rebin(OSFFdn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());

      auto denomSSfinal = nominal_.returnSS();
      auto denomOSfinal = nominal_.returnOS();

      auto denomfinal_nominal = nominal_.returnAdded();
      auto denomfinal_SSFFup = SSFFup_.returnAdded();
      auto denomfinal_SSFFdn = SSFFdn_.returnAdded();
      auto denomfinal_OSFFup = OSFFup_.returnAdded();
      auto denomfinal_OSFFdn = OSFFdn_.returnAdded();
      std::unique_ptr<TH1D> denomfinal_heepIdUp, denomfinal_heepIdDn;

      if ( TString(dataHist_->GetName()).Contains("invM") ) {
        heepIdUp_.rebin(heepIdUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        heepIdDn_.rebin(heepIdDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        denomfinal_heepIdUp = heepIdUp_.returnAdded();
        denomfinal_heepIdDn = heepIdDn_.returnAdded();
      }

      THStack* denomFinal = new THStack("final",";GeV;");
      denomOSfinal->SetLineWidth(0);
      denomSSfinal->SetLineWidth(0);
      denomFinal->Add(denomOSfinal.get());
      denomFinal->Add(denomSSfinal.get());

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMaximum(1.5*dataHist_->GetMaximum());

      if (TString(dataHist_->GetName()).Contains("invM")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,1000.);
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      dataHist_->Draw("E1");
      denomFinal->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("invM")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (!TString(dataHist_->GetName()).Contains("Et")) {
        TLegend* legend = new TLegend(0.55,0.7,0.95,0.9);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(denomSSfinal.get(),"Data-driven (X #rightarrow ME)");
        legend->AddEntry(denomOSfinal.get(),"Data-driven (e #rightarrow ME)");

        if (TString(dataHist_->GetName()).Contains("invM")) {
          legend->AddEntry(sigHist_.at(0),"M_{A} = 1 GeV (#sigma = 1 fb)");
          //legend_left->AddEntry(sigNums.at(isigDiv),"M_{A} = 10 GeV (#sigma = 10 fb)");
        }

        legend->Draw();
      }

      std::vector<double> x0, y0, errx, erryDn, erryUp;

      for (unsigned idx = 1; idx <= denomfinal_nominal->GetNbinsX(); idx++) {
        x0.push_back(denomfinal_nominal->GetBinCenter(idx));
        y0.push_back(denomfinal_nominal->GetBinContent(idx));
        errx.push_back(denomfinal_nominal->GetBinWidth(idx)/2.);

        double valSSFFup = denomfinal_SSFFup->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
        double valSSFFdn = denomfinal_nominal->GetBinContent(idx) - denomfinal_SSFFdn->GetBinContent(idx);
        double valOSFFup = denomfinal_OSFFup->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
        double valOSFFdn = denomfinal_nominal->GetBinContent(idx) - denomfinal_OSFFdn->GetBinContent(idx);

        double valHeepIdUp = 0., valHeepIdDn = 0.;

        if ( TString(dataHist_->GetName()).Contains("invM") ) {
          valHeepIdUp = denomfinal_heepIdUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valHeepIdDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_heepIdDn->GetBinContent(idx);
        }

        erryUp.push_back( std::sqrt( valSSFFup*valSSFFup + valOSFFup*valOSFFup + valHeepIdUp*valHeepIdUp ) );
        erryDn.push_back( std::sqrt( valSSFFdn*valSSFFdn + valOSFFdn*valOSFFdn + valHeepIdDn*valHeepIdDn ) );
      }

      auto gr = new TGraphAsymmErrors(denomfinal_nominal->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+2);
      gr->SetLineColor(kGray+2);
      gr->SetFillStyle(3004);
      gr->Draw("2");

      if (dir_) {
        dir_->WriteTObject(dataHist_,"data_obs");
        dir_->WriteTObject(nominal_.returnSS().release(),"SS");
        dir_->WriteTObject(nominal_.returnOS().release(),"OS");

        dir_->WriteTObject(SSFFup_.returnSS().release(),"SS_mergedEleFakeFactorSSUp");
        dir_->WriteTObject(SSFFdn_.returnSS().release(),"SS_mergedEleFakeFactorSSDown");
        dir_->WriteTObject(OSFFup_.returnOS().release(),"OS_mergedEleFakeFactorOSUp");
        dir_->WriteTObject(OSFFdn_.returnOS().release(),"OS_mergedEleFakeFactorOSDown");
        dir_->WriteTObject(heepIdUp_.returnSS().release(),"SS_modHeepIdUp");
        dir_->WriteTObject(heepIdUp_.returnOS().release(),"OS_modHeepIdUp");
        dir_->WriteTObject(heepIdDn_.returnSS().release(),"SS_modHeepIdDown");
        dir_->WriteTObject(heepIdDn_.returnOS().release(),"OS_modHeepIdDown");

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
          dir_->WriteTObject(sigHist_modHeepIdUp_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdUp");
          dir_->WriteTObject(sigHist_modHeepIdDn_.at(idx),sigFiles_.at(idx).name_+"_modHeepIdDown");
          dir_->WriteTObject(sigHist_mergedEleIdUp_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdUp");
          dir_->WriteTObject(sigHist_mergedEleIdDn_.at(idx),sigFiles_.at(idx).name_+"_mergedEleIdDown");
        }
      }

      denomSSfinal.release();
      denomOSfinal.release();
    }
  };

  canvas_1->cd();

  auto aloader3E = HistLoader3E(datafile,WZfile,ZZfile,sigsamples);

  //canvas_1->SetLogy();
  aloader3E.load("2M_CRME_lll_invM","2M_antiME_lll_invM_CR");
  aloader3E.preparecard("MEMu2M_"+era+"_datacard.root","mergedEMu2M");
  aloader3E.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_2M_invM.png");
  aloader3E.close();
  canvas_1->SetLogy(0);

  aloader3E.load("2M_Et_CR_EB_CRME","2M_Et_CR_EB_antiME");
  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2M_Et.png");

  aloader3E.load("2M_CRME_all_eta","2M_antiME_eta");
  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2M_eta.png");

  canvas_2->cd();

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

  TFile* afile = new TFile("MMFF_"+era+".root","RECREATE");

  std::vector<double> xbins = {0,150,200,250,300,350,400,500,600,800,1000,1200,1500};
  //std::vector<double> xbins = {0,100,150,200,250,275,300,350,400,450,500,600,700,800,900,1000,1200,1500,2000};
  const int nbins = xbins.size()-1;
  std::vector<double> xcen = estimateCenter(xbins);

  TH1D* numer = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_CRME_antiDphi_mt")->Clone();
  TH1D* denom = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_CRME_antiDphi_antiRpt_mt")->Clone();
  //TH1D* numer = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_CRME_antiDphi_MET_sumEt")->Clone();
  //TH1D* denom = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_CRME_antiDphi_antiRpt_MET_sumEt")->Clone();

  numer = (TH1D*)numer->Rebin( nbins,"1M_MMFF_numer_rebin",&(xbins[0]) );
  denom = (TH1D*)denom->Rebin( nbins,"denom_rebin",&(xbins[0]) );

  numer->Divide(denom);
  //numer->SetMaximum(0.7);
  numer->SetMinimum(0.0);
  numer->GetYaxis()->SetTitle("Fake factor");
  //numer->GetXaxis()->SetTitle("p_{T}");
  numer->SetLineColor(kBlack);
  numer->SetLineWidth(2);
  numer->SetStats(0);

  TF1* ffFunc = new TF1("MMFF","[0]*x+[1]",150,2000);
  ffFunc->SetLineColor(kRed);
  ffFunc->SetLineWidth(2);
  ffFunc->SetLineStyle(2);

  TFitResultPtr fitResult = numer->Fit(ffFunc,"RS");
  fitResult->SetName("fitResult");
  double ci[nbins];
  fitResult->GetConfidenceIntervals(nbins,1,0,&(xcen[0]),ci,0.95,false); // 0.6827
  std::vector<double> xbinw = estimateWidth(xbins);
  double ybin[nbins];

  for (unsigned idx = 0; idx < nbins; idx++) {
    ybin[idx] = ffFunc->Eval(xcen[idx]);
  }

  auto errGr = new TGraphErrors(nbins,&(xcen[0]),ybin,&(xbinw[0]),ci);
  errGr->SetFillColor(kRed);
  errGr->SetFillStyle(3003);

  numer->Draw("E1");
  errGr->Draw("3");

  TPaveText* texthigh = new TPaveText(0.65,0.82,0.95,0.9,"NDC");
  texthigh->SetBorderSize(0);
  texthigh->SetFillColor(0);
  texthigh->SetFillStyle(3025);
  TString textoslow;
  textoslow.Form("%.3g #times x + %.4f", ffFunc->GetParameter(0), ffFunc->GetParameter(1));
  texthigh->AddText(textoslow);
  ((TText*)texthigh->GetListOfLines()->Last())->SetTextColor(kRed);
  ((TText*)texthigh->GetListOfLines()->Last())->SetTextAlign(12);
  texthigh->Draw();

  SaveAs(canvas_2,"FF_1M_FF.png");

  numer->Write();
  ffFunc->Write();
  fitResult->Write();
  afile->Close();

  /*TH2D* srHist = (TH2D*)H2000A1file->Get("mergedEMuCRanalyzer20UL18/1M_dphi_ratioPt");
  canvas_1->cd();
  srHist->Draw("col");
  lumi_13TeV = "";
  extraText  = "   Simulation Work in progress";
  SaveAs(canvas_1,"srHist_H2000A1.png");*/

  return;
}
