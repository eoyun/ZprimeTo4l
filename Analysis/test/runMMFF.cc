#include "TFile.h"
#include "TH1D.h"
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
  } else {
    std::cout << "check era..." << std::endl;
  }

  static TString anlyzrMC = "mergedMuCRanalyzer"+postfix;
  static TString anlyzrData = "mergedMuCRanalyzerData";

  static TString anlyzrEMuMC = "mergedEMuCRanalyzer"+postfix;
  static TString anlyzrEMuData = "mergedEMuCRanalyzerData";

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

  TString firstEra = era=="run2" ? "20UL16" : era;

  TFile* datafile = new TFile("MuAnalyzer_"+firstEra+"_data.root","READ");

  TFile *datafile1;
  TFile *datafile2;
  TFile *datafile3;

  if (era.Contains("run2")) {
    datafile1 = new TFile("MuAnalyzer_20UL16APV_data.root","READ");
    datafile2 = new TFile("MuAnalyzer_20UL17_data.root","READ");
    datafile3 = new TFile("MuAnalyzer_20UL18_data.root","READ");
  }

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  std::vector<SigSample> sigsamples = {
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A750.root","READ"),"H2000A750")
  };

  std::vector<SigSample> sigsamples1, sigsamples2, sigsamples3;

  if (era.Contains("run2")) {
    sigsamples1 = {
      SigSample(new TFile("MuAnalyzer_20UL16APV_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples2 = {
      SigSample(new TFile("MuAnalyzer_20UL17_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("MuAnalyzer_20UL17_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("MuAnalyzer_20UL17_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("MuAnalyzer_20UL17_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("MuAnalyzer_20UL17_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("MuAnalyzer_20UL17_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("MuAnalyzer_20UL17_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("MuAnalyzer_20UL17_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples3 = {
      SigSample(new TFile("MuAnalyzer_20UL18_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("MuAnalyzer_20UL18_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("MuAnalyzer_20UL18_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("MuAnalyzer_20UL18_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("MuAnalyzer_20UL18_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("MuAnalyzer_20UL18_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("MuAnalyzer_20UL18_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("MuAnalyzer_20UL18_H2000A750.root","READ"),"H2000A750")
    };
  }

  sigsamples = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_"+firstEra+"_H2000A1.root","READ"),"H2000A1")};
  sigsamples1 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A1.root","READ"),"H2000A1")};
  sigsamples2 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL17_H2000A1.root","READ"),"H2000A1")};
  sigsamples3 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL18_H2000A1.root","READ"),"H2000A1")};

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
  p1->SetBottomMargin(0.00);
  p1->SetLeftMargin( L/W );
  p1->SetRightMargin( R/W );
  p1->Draw();

  TPad* p2 = new TPad("p2","",0,0,1,0.3);
  p2->SetTickx(0);
  p2->SetTicky(0);
  p2->SetTopMargin(0.0);
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

/*  TH1D* nm = (TH1D*)datafile->Get("mergedMuCRanalyzerData/3M_antiRpt_mt")->Clone();
  TH1D* dm = (TH1D*)datafile->Get("mergedMuCRanalyzerData/3M_antiDphi_antiRpt_mt")->Clone();

  nm->Rebin( nm->GetNbinsX() );
  dm->Rebin( dm->GetNbinsX() );

  nm->Divide(nm,dm,1,1,"B");
  const double ff = nm->GetBinContent(1);
  const double ffErr = dm->GetBinError(1);*/

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

    TH1D* variateDn(TH1D* nominal, TH1D* up) const {
      TH1D* dn = (TH1D*)nominal->Clone();

      for (int ibin=0; ibin<nominal->GetNbinsX()+2; ibin++) {
        double val = nominal->GetBinContent(ibin);
        double var = up->GetBinContent(ibin);
        double valErr = nominal->GetBinError(ibin);
        double varErr = up->GetBinError(ibin);
        dn->SetBinContent( ibin, std::max( val - (var - val), 0.) );
        dn->SetBinError( ibin, std::hypot(2*valErr,varErr) );
      }

      return dn;
    }

    class SystVariation {
    public:
      SystVariation() : up_(nullptr), dn_(nullptr) {}
      SystVariation(TH1D* up, TH1D* dn) : up_(up), dn_(dn) {}
      ~SystVariation()=default;

      TH1D* up_;
      TH1D* dn_;
    };

  protected:
    std::map<std::string,SystVariation> syst_;

    TH1D* dataHist_ = nullptr;
    TH1D* FFHist_ = nullptr;

    std::vector<TH1D*> sigHist_;
    std::vector<std::map<std::string,SystVariation>> sigSyst_;

  public:
    void add(const HistLoaderBase& other) {
      this->dataHist_->Add(other.dataHist_);
      this->FFHist_->Add(other.FFHist_);

      if (!syst_.empty()) {
        for (const auto& element : this->syst_) {
          this->syst_.at(element.first).up_->Add(other.syst_.at(element.first).up_);
          this->syst_.at(element.first).dn_->Add(other.syst_.at(element.first).dn_);
        }
      }

      if (!sigHist_.empty()) {
        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          this->sigHist_.at(idx)->Add(other.sigHist_.at(idx));

          if (!syst_.empty()) {
            for (const auto& element : this->sigSyst_.at(idx)) {
              this->sigSyst_.at(idx).at(element.first).up_->Add(other.sigSyst_.at(idx).at(element.first).up_);
              this->sigSyst_.at(idx).at(element.first).dn_->Add(other.sigSyst_.at(idx).at(element.first).dn_);
            }
          }
        }
      }
    }
  };

  class HistLoaderMM : public HistLoaderBase {
  public:
    HistLoaderMM(TFile* adatafile, std::vector<SigSample> sigsamples)
    : HistLoaderBase(adatafile,nullptr,nullptr), sigFiles_(sigsamples) {}

    ~HistLoaderMM()=default;

  private:
    std::vector<SigSample> sigFiles_;

  public:
    void load(TString numName, TString denomName, std::string anlyzrEra="") {
      if (dataHist_) {
        delete dataHist_, FFHist_;
        dataHist_ = nullptr;
      }

      if (!syst_.empty())
        syst_.clear();

      dataHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+numName)->Clone();
      FFHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+denomName)->Clone();
      FFHist_->SetFillColor(kGray);
      FFHist_->SetLineWidth(5);
      FFHist_->SetLineColor(kGray);

      sigHist_.clear();
      sigSyst_.clear();

      if ( numName.Contains("3M_mt") || numName.Contains("3M_CRdphi_mt") ) {
        TH1D* FFHistUp = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_up" )->Clone();
        TH1D* FFHistDn = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_dn" )->Clone();
        syst_["FFHist"] = SystVariation(FFHistUp,FFHistDn);

        TH1D* JESup = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_JESup" )->Clone();
        TH1D* JESdn = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_JESdn" )->Clone();
        syst_["JES"] = SystVariation(JESup,JESdn);
      }

      if (numName.Contains("3M_mt")) {
        const double lumi = retrieveLumi(anlyzrEra);
        const double sigLumi = 0.01;

        auto retrieveSigHist = [this,&numName,&sigLumi,&lumi,&anlyzrEra] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( TString("mergedMuCRanalyzer")+anlyzrEra.c_str()+"/"+numName+systName )->Clone();
          ahist->Scale( lumi*1000.*sigLumi / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned idx = 0; idx < sigFiles_.size(); idx++) {
          sigHist_.push_back( retrieveSigHist(sigFiles_.at(idx).file_,"") );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          std::map<std::string,SystVariation> init;

          TH1D* sigJESup = retrieveSigHist(sigFiles_.at(idx).file_,"_JESup");
          TH1D* sigJESdn = retrieveSigHist(sigFiles_.at(idx).file_,"_JESdn");
          init["sigJES"] = SystVariation(sigJESup,sigJESdn);
          TH1D* sigJERup = retrieveSigHist(sigFiles_.at(idx).file_,"_JERup");
          TH1D* sigJERdn = retrieveSigHist(sigFiles_.at(idx).file_,"_JERdn");
          init["sigJER"] = SystVariation(sigJERup,sigJERdn);
          TH1D* sigMuIdUp = retrieveSigHist(sigFiles_.at(idx).file_,"_idUp");
          TH1D* sigMuIdDn = retrieveSigHist(sigFiles_.at(idx).file_,"_idDn");
          init["sigMuId"] = SystVariation(sigMuIdUp,sigMuIdDn);
          TH1D* sigMuIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_isoUp");
          TH1D* sigMuIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_isoDn");
          init["sigMuIso"] = SystVariation(sigMuIsoUp,sigMuIsoDn);
          TH1D* sigMuBoostIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoUp");
          TH1D* sigMuBoostIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoDn");
          init["sigMuBoostIso"] = SystVariation(sigMuBoostIsoUp,sigMuBoostIsoDn);
          TH1D* sigTrigUp = retrieveSigHist(sigFiles_.at(idx).file_,"_trigUp");
          TH1D* sigTrigDn = retrieveSigHist(sigFiles_.at(idx).file_,"_trigDn");
          init["sigTrig"] = SystVariation(sigTrigUp,sigTrigDn);
          TH1D* sigMuRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_recoUp");
          TH1D* sigMuRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_recoDn");
          init["sigMuReco"] = SystVariation(sigMuRecoUp,sigMuRecoDn);
          TH1D* sigPUrwgtUp = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtUp");
          TH1D* sigPUrwgtDn = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtDn");
          init["sigPUrwgt"] = SystVariation(sigPUrwgtUp,sigPUrwgtDn);
          TH1D* sigPrefireUp = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireUp");
          TH1D* sigPrefireDn = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireDn");
          init["sigPrefire"] = SystVariation(sigPrefireUp,sigPrefireDn);

          sigSyst_.push_back(init);
        }
      }
    }

    void compare(TPad* padUp, int rebin=1, double scale=1., TPad* padDn=nullptr) {
      if ( dataHist_->GetNbinsX()!=FFHist_->GetNbinsX() ) {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( element.second.up_->GetNbinsX()/dataHist_->GetNbinsX() );
            element.second.dn_->Rebin( element.second.dn_->GetNbinsX()/dataHist_->GetNbinsX() );
          }
        }
      }

      if (rebin!=1) {
        dataHist_->Rebin(rebin);
        FFHist_->Rebin(rebin);

       if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( rebin );
            element.second.dn_->Rebin( rebin );
          }
        }

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);

          for (const auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_->Rebin(rebin);
            sigSyst_.at(idx).at(element.first).dn_->Rebin(rebin);
          }
        }
      } else if (true) {
        std::vector<double> binEdges = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,2500.};

        dataHist_ = (TH1D*)dataHist_->Rebin(binEdges.size()-1,TString(dataHist_->GetName())+"_rebin",&(binEdges[0]));
        FFHist_ = (TH1D*)FFHist_->Rebin(binEdges.size()-1,TString(FFHist_->GetName())+"_rebin",&(binEdges[0]));

        if (!syst_.empty()) {
          for (auto& element : syst_) {
            element.second.up_ = (TH1D*)element.second.up_->Rebin( binEdges.size()-1,TString(element.second.up_->GetName())+"_rebin",&(binEdges[0]) );
            element.second.dn_ = (TH1D*)element.second.dn_->Rebin( binEdges.size()-1,TString(element.second.dn_->GetName())+"_rebin",&(binEdges[0]) );
          }
        }

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx) = (TH1D*)sigHist_.at(idx)->Rebin( binEdges.size()-1,TString(sigHist_.at(idx)->GetName())+"_rebin",&(binEdges[0]) );

          for (auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_ = (TH1D*)sigSyst_.at(idx).at(element.first).up_->Rebin( binEdges.size()-1,TString(sigSyst_.at(idx).at(element.first).up_->GetName())+"_rebin",&(binEdges[0]) );
            sigSyst_.at(idx).at(element.first).dn_ = (TH1D*)sigSyst_.at(idx).at(element.first).dn_->Rebin( binEdges.size()-1,TString(sigSyst_.at(idx).at(element.first).dn_->GetName())+"_rebin",&(binEdges[0]) );
          }
        }
      }

      if (scale!=1.) {
        FFHist_->Scale(scale);

        for (const auto& element : syst_) {
          element.second.up_->Scale(scale);
          element.second.dn_->Scale(scale);
        }
      }

      padUp->cd();

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMaximum(2.*dataHist_->GetMaximum());

      if (TString(dataHist_->GetName()).Contains("mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,2500.);
        dataHist_->SetMinimum(0.001);
        dataHist_->SetMaximum(3.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      dataHist_->Draw("E1");
      FFHist_->Draw("hist&same");
      FFHist_->Draw("E0X0&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("3M_mt")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (TString(dataHist_->GetName()).Contains("3M_mt") || TString(dataHist_->GetName()).Contains("3M_CRdphi_mt")) {
        double ylow = TString(dataHist_->GetName()).Contains("3M_CRdphi_mt") ? 0.75 : 0.68;
        TLegend* legend = new TLegend(0.75,ylow,0.93,0.93);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(FFHist_,"Bkg");

        if (TString(dataHist_->GetName()).Contains("3M_mt"))
          legend->AddEntry(sigHist_.at(0),"X2000Y1");

        legend->Draw();
      }

      if (padDn) {
        TH1D* ratio = (TH1D*)dataHist_->Clone();
        ratio->SetStats(0);
        ratio->SetTitle("");
        ratio->Divide(FFHist_);
        ratio->GetYaxis()->SetTitle("Obs/Exp");
        ratio->GetYaxis()->SetTitleSize(0.1);
        ratio->GetYaxis()->SetTitleOffset(0.4);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetXaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetLabelOffset(0.005);
        ratio->GetYaxis()->SetRangeUser(0.2,1.8);
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetTitleOffset(0.75);
        ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
        ratio->SetLineColor(kBlack);

        padDn->cd();
        ratio->Draw("E1");
        padUp->cd();
      }

      if (!syst_.empty()) {
        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRup, errRdn;

        for (unsigned idx = 1; idx <= FFHist_->GetNbinsX(); idx++) {
          x0.push_back(FFHist_->GetBinCenter(idx));
          y0.push_back(FFHist_->GetBinContent(idx));
          errx.push_back(FFHist_->GetBinWidth(idx)/2.);

          double valFFup = syst_.at("FFHist").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
          double valFFdn = FFHist_->GetBinContent(idx) - syst_.at("FFHist").dn_->GetBinContent(idx);
          double valJESup = 0., valJESdn = 0.;

          if (syst_.count("JES")) {
            valJESup = syst_.at("JES").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
            valJESdn = FFHist_->GetBinContent(idx) - syst_.at("JES").dn_->GetBinContent(idx);
          }

          erryUp.push_back( std::sqrt(valFFup*valFFup + valJESup*valJESup) );
          erryDn.push_back( std::sqrt(valFFdn*valFFdn + valJESdn*valJESdn) );

          r0.push_back(1.);

          double rFFup = valFFup/FFHist_->GetBinContent(idx);
          double rFFdn = valFFdn/FFHist_->GetBinContent(idx);
          double rJESup = valJESup/FFHist_->GetBinContent(idx);
          double rJESdn = valJESdn/FFHist_->GetBinContent(idx);

          errRup.push_back( FFHist_->GetBinContent(idx) > 0. ? std::sqrt( rFFup*rFFup + rJESup*rJESup ) : 0. );
          errRdn.push_back( FFHist_->GetBinContent(idx) > 0. ? std::sqrt( rFFdn*rFFdn + rJESdn*rJESdn ) : 0. );
        }

        auto gr = new TGraphAsymmErrors(FFHist_->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kBlack);
        gr->SetLineColor(kBlack);
        gr->SetFillStyle(3004);
        gr->Draw("2");

        if (padDn) {
          auto rgr = new TGraphAsymmErrors(FFHist_->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
          rgr->SetFillColor(kGray+2);
          rgr->SetLineColor(kGray+2);
          rgr->SetFillStyle(3004);
          padDn->cd();
          rgr->Draw("2");
        }
      }

      if (dir_) {
        dir_->WriteTObject(dataHist_,"data_obs");
        dir_->WriteTObject(FFHist_,"MM");

        dir_->WriteTObject(syst_["FFHist"].up_,"MM_mergedMuFakeFactorUp");
        dir_->WriteTObject(syst_["FFHist"].dn_,"MM_mergedMuFakeFactorDown");
        dir_->WriteTObject(syst_["JES"].up_,"MM_JESUp");
        dir_->WriteTObject(syst_["JES"].dn_,"MM_JESDown");

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJES").up_,sigFiles_.at(idx).name_+"_JESUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJES").dn_,sigFiles_.at(idx).name_+"_JESDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJER").up_,sigFiles_.at(idx).name_+"_JERUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJER").dn_,sigFiles_.at(idx).name_+"_JERDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").up_,sigFiles_.at(idx).name_+"_highPtIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").dn_,sigFiles_.at(idx).name_+"_highPtIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").up_,sigFiles_.at(idx).name_+"_muLooseIsoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").dn_,sigFiles_.at(idx).name_+"_muLooseIsoDown");
          //dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").up_,sigFiles_.at(idx).name_+"_muBoostIsoUp");
          //dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").dn_,sigFiles_.at(idx).name_+"_muBoostIsoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").up_,sigFiles_.at(idx).name_+"_muTrigUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").dn_,sigFiles_.at(idx).name_+"_muTrigDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").up_,sigFiles_.at(idx).name_+"_muRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").dn_,sigFiles_.at(idx).name_+"_muRecoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").up_,sigFiles_.at(idx).name_+"_PUrwgtUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").dn_,sigFiles_.at(idx).name_+"_PUrwgtDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").up_,sigFiles_.at(idx).name_+"_prefireUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").dn_,sigFiles_.at(idx).name_+"_prefireDown");
        }
      }
    }
  }; // class

  auto aloaderMM = HistLoaderMM(datafile,sigsamples);
  auto aloaderMM1 = HistLoaderMM(datafile1,sigsamples1);
  auto aloaderMM2 = HistLoaderMM(datafile2,sigsamples2);
  auto aloaderMM3 = HistLoaderMM(datafile3,sigsamples3);

  aloaderMM.load("3M_mt","3M_antiRpt_mt_xFF",postfix.Data()); // 3M_antiRpt_mt_xFF 3M_antiDphi_antiRpt_mt_xFF

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_mt","3M_antiRpt_mt_xFF","20UL16APV");
    aloaderMM2.load("3M_mt","3M_antiRpt_mt_xFF","20UL17");
    aloaderMM3.load("3M_mt","3M_antiRpt_mt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  //aloaderMM.preparecard("MMFF_"+era+"_datacard.root","mergedMu3M");
  aloaderMM.compare(p1,1,1.,p2); // 10
  SaveAs(canvas_2,"MMFF_3M_mt.pdf",p1);
  //aloaderMM.close();

  /*aloaderMM.load("3M_CRdphi_mt","3M_CRdphi_antiRpt_mt_xFF");

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_CRdphi_mt","3M_CRdphi_antiRpt_mt_xFF","20UL16APV");
    aloaderMM2.load("3M_CRdphi_mt","3M_CRdphi_antiRpt_mt_xFF","20UL17");
    aloaderMM3.load("3M_CRdphi_mt","3M_CRdphi_antiRpt_mt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  aloaderMM.compare(p1,10,1.,p2); // ff
  SaveAs(canvas_2,"MMFF_3M_CRdphi_mt.pdf",p1);*/

/*  aloaderMM.load("3M_antiDphi_mt","3M_antiDphi_antiRpt_mt_xFF");

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_antiDphi_mt","3M_antiDphi_antiRpt_mt_xFF","20UL16APV");
    aloaderMM2.load("3M_antiDphi_mt","3M_antiDphi_antiRpt_mt_xFF","20UL17");
    aloaderMM3.load("3M_antiDphi_mt","3M_antiDphi_antiRpt_mt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  aloaderMM.compare(canvas_1,5);
  SaveAs(canvas_1,"FF_3M_antiDphi_mt.png");

  aloaderMM.load("3M_antiDphi_MET_pt","3M_antiDphi_antiRpt_MET_pt_xFF");

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_antiDphi_MET_pt","3M_antiDphi_antiRpt_MET_pt_xFF","20UL16APV");
    aloaderMM2.load("3M_antiDphi_MET_pt","3M_antiDphi_antiRpt_MET_pt_xFF","20UL17");
    aloaderMM3.load("3M_antiDphi_MET_pt","3M_antiDphi_antiRpt_MET_pt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_antiDphi_MET_pt.png");

  aloaderMM.load("3M_CRdphi_MM_pt","3M_CRdphi_antiRpt_MM_pt_xFF");

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_CRdphi_MM_pt","3M_CRdphi_antiRpt_MM_pt_xFF","20UL16APV");
    aloaderMM2.load("3M_CRdphi_MM_pt","3M_CRdphi_antiRpt_MM_pt_xFF","20UL17");
    aloaderMM3.load("3M_CRdphi_MM_pt","3M_CRdphi_antiRpt_MM_pt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_CRdphi_MET_pt.png");

  aloaderMM.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF");

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF","20UL16APV");
    aloaderMM2.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF","20UL17");
    aloaderMM3.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_antiDphi_MM_pt.png");

  aloaderMM.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF");

  if (era.Contains("run2")) {
    aloaderMM1.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF","20UL16APV");
    aloaderMM2.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF","20UL17");
    aloaderMM3.load("3M_antiDphi_MM_pt","3M_antiDphi_antiRpt_MM_pt_xFF","20UL18");
    aloaderMM.add(aloaderMM1);
    aloaderMM.add(aloaderMM2);
    aloaderMM.add(aloaderMM3);
  }

  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_3M_CRdphi_MM_pt.png");*/

  class HistLoader2E : public HistLoaderBase {
  public:
    HistLoader2E(TFile* adatafile, std::vector<SigSample> sigsamples)
    : HistLoaderBase(adatafile,nullptr,nullptr), sigFiles_(sigsamples) {}

    ~HistLoader2E()=default;

  private:
    std::vector<SigSample> sigFiles_;

  public:
    void load(TString numName, TString denomName, std::string anlyzrEra="") {
      if (dataHist_) {
        delete dataHist_, FFHist_;
        dataHist_ = nullptr;
      }

      if (!syst_.empty())
        syst_.clear();

      dataHist_ = (TH1D*)datafile_->Get(anlyzrEMuData+"/"+numName)->Clone();
      FFHist_ = (TH1D*)datafile_->Get(anlyzrEMuData+"/"+denomName)->Clone();
      FFHist_->SetFillColor(kGray);
      FFHist_->SetLineWidth(5);
      FFHist_->SetLineColor(kGray);

      sigHist_.clear();
      sigSyst_.clear();

      if ( numName.Contains("2E_mt") || numName.Contains("2E_CRdphi_mt") ) {
        TH1D* FFHistUp = (TH1D*)datafile_->Get( std::string(anlyzrEMuData+"/")+denomName+"_up" )->Clone();
        TH1D* FFHistDn = (TH1D*)datafile_->Get( std::string(anlyzrEMuData+"/")+denomName+"_dn" )->Clone();
        syst_["FFHist"] = SystVariation(FFHistUp,FFHistDn);

        TH1D* JESup = (TH1D*)datafile_->Get( std::string(anlyzrEMuData+"/"+denomName+"_JESup").c_str() )->Clone();
        TH1D* JESdn = (TH1D*)datafile_->Get( std::string(anlyzrEMuData+"/"+denomName+"_JESdn").c_str() )->Clone();
        syst_["JES"] = SystVariation(JESup,JESdn);
      }

      if (numName.Contains("2E_mt")) {
        const double lumi = retrieveLumi(anlyzrEra);
        const double sigLumi = 0.0001;

        auto retrieveSigHist = [this,&numName,&sigLumi,&lumi,&anlyzrEra] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( TString("mergedEMuCRanalyzer")+anlyzrEra.c_str()+"/"+numName+systName )->Clone();
          ahist->Scale( lumi*1000.*sigLumi / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

/*        for (unsigned idx = 0; idx < sigFiles_.size(); idx++) {
          sigHist_.push_back( retrieveSigHist(sigFiles_.at(idx).file_,"") );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          std::map<std::string,SystVariation> init;

          TH1D* sigJESup = retrieveSigHist(sigFiles_.at(idx).file_,"_JESup");
          TH1D* sigJESdn = retrieveSigHist(sigFiles_.at(idx).file_,"_JESdn");
          init["sigJES"] = SystVariation(sigJESup,sigJESdn);
          TH1D* sigJERup = retrieveSigHist(sigFiles_.at(idx).file_,"_JERup");
          TH1D* sigJERdn = retrieveSigHist(sigFiles_.at(idx).file_,"_JERdn");
          init["sigJER"] = SystVariation(sigJERup,sigJERdn);
          TH1D* sigMuIdUp = retrieveSigHist(sigFiles_.at(idx).file_,"_idUp");
          TH1D* sigMuIdDn = retrieveSigHist(sigFiles_.at(idx).file_,"_idDn");
          init["sigMuId"] = SystVariation(sigMuIdUp,sigMuIdDn);
          TH1D* sigMuIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_isoUp");
          TH1D* sigMuIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_isoDn");
          init["sigMuIso"] = SystVariation(sigMuIsoUp,sigMuIsoDn);
          TH1D* sigMuBoostIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoUp");
          TH1D* sigMuBoostIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoDn");
          init["sigMuBoostIso"] = SystVariation(sigMuIsoUp,sigMuIsoDn);
          TH1D* sigTrigUp = retrieveSigHist(sigFiles_.at(idx).file_,"_trigUp");
          TH1D* sigTrigDn = retrieveSigHist(sigFiles_.at(idx).file_,"_trigDn");
          init["sigTrig"] = SystVariation(sigTrigUp,sigTrigDn);
          TH1D* sigMuRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_recoUp");
          TH1D* sigMuRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_recoDn");
          init["sigMuReco"] = SystVariation(sigMuRecoUp,sigMuRecoDn);

          sigSyst_.push_back(init);
        }*/
      }
    }

    void compare(TPad* padUp, int rebin=1, double scale=1., TPad* padDn=nullptr) {
      if ( dataHist_->GetNbinsX()!=FFHist_->GetNbinsX() ) {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( element.second.up_->GetNbinsX()/dataHist_->GetNbinsX() );
            element.second.dn_->Rebin( element.second.dn_->GetNbinsX()/dataHist_->GetNbinsX() );
          }
        }
      }

      if (rebin!=1) {
        dataHist_->Rebin(rebin);
        FFHist_->Rebin(rebin);

       if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( rebin );
            element.second.dn_->Rebin( rebin );
          }
        }

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);

          for (const auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_->Rebin(rebin);
            sigSyst_.at(idx).at(element.first).dn_->Rebin(rebin);
          }
        }
      } else if (true) {
        std::vector<double> binEdges = {0.,50.,100.,150.,200.,250.,300.,350.,400.,450.,500.,600.,700.,800.,2500.};

        dataHist_ = (TH1D*)dataHist_->Rebin(binEdges.size()-1,TString(dataHist_->GetName())+"_rebin",&(binEdges[0]));
        FFHist_ = (TH1D*)FFHist_->Rebin(binEdges.size()-1,TString(FFHist_->GetName())+"_rebin",&(binEdges[0]));

        if (!syst_.empty()) {
          for (auto& element : syst_) {
            element.second.up_ = (TH1D*)element.second.up_->Rebin( binEdges.size()-1,TString(element.second.up_->GetName())+"_rebin",&(binEdges[0]) );
            element.second.dn_ = (TH1D*)element.second.dn_->Rebin( binEdges.size()-1,TString(element.second.dn_->GetName())+"_rebin",&(binEdges[0]) );
          }
        }

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx) = (TH1D*)sigHist_.at(idx)->Rebin( binEdges.size()-1,TString(sigHist_.at(idx)->GetName())+"_rebin",&(binEdges[0]) );

          for (auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_ = (TH1D*)sigSyst_.at(idx).at(element.first).up_->Rebin( binEdges.size()-1,TString(sigSyst_.at(idx).at(element.first).up_->GetName())+"_rebin",&(binEdges[0]) );
            sigSyst_.at(idx).at(element.first).dn_ = (TH1D*)sigSyst_.at(idx).at(element.first).dn_->Rebin( binEdges.size()-1,TString(sigSyst_.at(idx).at(element.first).dn_->GetName())+"_rebin",&(binEdges[0]) );
          }
        }
      }

      if (scale!=1.) {
        FFHist_->Scale(scale);

        for (const auto& element : syst_) {
          element.second.up_->Scale(scale);
          element.second.dn_->Scale(scale);
        }
      }

      padUp->cd();

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMaximum(3.*dataHist_->GetMaximum());

      if (TString(dataHist_->GetName()).Contains("2E_mt")) {
        dataHist_->SetMaximum(1.5*dataHist_->GetMaximum());
      }

      if (TString(dataHist_->GetName()).Contains("mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,2500.);
        dataHist_->SetMinimum(0.001);
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      dataHist_->Draw("E1");
      FFHist_->Draw("hist&same");
      FFHist_->Draw("E0X0&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("2E_mt")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (TString(dataHist_->GetName()).Contains("2E_mt") || TString(dataHist_->GetName()).Contains("2E_CRdphi_mt")) {
        double ylow = TString(dataHist_->GetName()).Contains("2E_CRdphi_mt") ? 0.75 : 0.75; // 0.68
        TLegend* legend = new TLegend(0.75,ylow,0.93,0.93);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(FFHist_,"Bkg");

//        if (TString(dataHist_->GetName()).Contains("2E_mt"))
//          legend->AddEntry(sigHist_.at(0),"X2000Y1");

        legend->Draw();
      }

      if (padDn) {
        TH1D* ratio = (TH1D*)dataHist_->Clone();
        ratio->SetStats(0);
        ratio->SetTitle("");
        ratio->Divide(FFHist_);
        ratio->GetYaxis()->SetTitle("Obs/Exp");
        ratio->GetYaxis()->SetTitleSize(0.1);
        ratio->GetYaxis()->SetTitleOffset(0.4);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetLabelSize(0.08);
        ratio->GetXaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetLabelOffset(0.005);
        ratio->GetYaxis()->SetRangeUser(0.2,1.8);
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetTitleOffset(0.75);
        ratio->GetYaxis()->SetNdivisions(5, 5, 0, kTRUE);
        ratio->SetLineColor(kBlack);

        padDn->cd();
        ratio->Draw("E1");
        padUp->cd();
      }

      if (!syst_.empty()) {
        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRup, errRdn;

        for (unsigned idx = 1; idx <= FFHist_->GetNbinsX(); idx++) {
          x0.push_back(FFHist_->GetBinCenter(idx));
          y0.push_back(FFHist_->GetBinContent(idx));
          errx.push_back(FFHist_->GetBinWidth(idx)/2.);

          double valFFup = syst_.at("FFHist").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
          double valFFdn = FFHist_->GetBinContent(idx) - syst_.at("FFHist").dn_->GetBinContent(idx);
          double valJESup = 0., valJESdn = 0.;

          if (syst_.count("JES")) {
            valJESup = syst_.at("JES").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
            valJESdn = FFHist_->GetBinContent(idx) - syst_.at("JES").dn_->GetBinContent(idx);
          }

          erryUp.push_back( std::sqrt(valFFup*valFFup + valJESup*valJESup) );
          erryDn.push_back( std::sqrt(valFFdn*valFFdn + valJESdn*valJESdn) );

          r0.push_back(1.);

          double rFFup = valFFup/FFHist_->GetBinContent(idx);
          double rFFdn = valFFdn/FFHist_->GetBinContent(idx);
          double rJESup = valJESup/FFHist_->GetBinContent(idx);
          double rJESdn = valJESdn/FFHist_->GetBinContent(idx);

          errRup.push_back( FFHist_->GetBinContent(idx) > 0. ? std::sqrt( rFFup*rFFup + rJESup*rJESup ) : 0. );
          errRdn.push_back( FFHist_->GetBinContent(idx) > 0. ? std::sqrt( rFFdn*rFFdn + rJESdn*rJESdn ) : 0. );
        }

        auto gr = new TGraphAsymmErrors(FFHist_->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kBlack);
        gr->SetLineColor(kBlack);
        gr->SetFillStyle(3004);
        gr->Draw("2");

        if (padDn) {
          auto rgr = new TGraphAsymmErrors(FFHist_->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
          rgr->SetFillColor(kGray+2);
          rgr->SetLineColor(kGray+2);
          rgr->SetFillStyle(3004);
          padDn->cd();
          rgr->Draw("2");
        }
      }

      if (dir_) {
        dir_->WriteTObject(dataHist_,"data_obs");
        dir_->WriteTObject(FFHist_,"MM");

        dir_->WriteTObject(syst_["FFHist"].up_,"MM_mergedMuFakeFactorUp");
        dir_->WriteTObject(syst_["FFHist"].dn_,"MM_mergedMuFakeFactorDown");
        dir_->WriteTObject(syst_["JES"].up_,"MM_JESUp");
        dir_->WriteTObject(syst_["JES"].dn_,"MM_JESDown");

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJES").up_,sigFiles_.at(idx).name_+"_JESUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJES").dn_,sigFiles_.at(idx).name_+"_JESDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJER").up_,sigFiles_.at(idx).name_+"_JERUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJER").dn_,sigFiles_.at(idx).name_+"_JERDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").up_,sigFiles_.at(idx).name_+"_highPtIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").dn_,sigFiles_.at(idx).name_+"_highPtIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").up_,sigFiles_.at(idx).name_+"_muLooseIsoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").dn_,sigFiles_.at(idx).name_+"_muLooseIsoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").up_,sigFiles_.at(idx).name_+"_muBoostIsoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").dn_,sigFiles_.at(idx).name_+"_muBoostIsoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").up_,sigFiles_.at(idx).name_+"_muTrigUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").dn_,sigFiles_.at(idx).name_+"_muTrigDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").up_,sigFiles_.at(idx).name_+"_muRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").dn_,sigFiles_.at(idx).name_+"_muRecoDown");
        }
      }
    }
  }; // class

  auto aloader2E = HistLoader2E(datafile,sigsamples);
  auto aloader2E1 = HistLoader2E(datafile1,sigsamples1);
  auto aloader2E2 = HistLoader2E(datafile2,sigsamples2);
  auto aloader2E3 = HistLoader2E(datafile3,sigsamples3);

  aloader2E.load("2E_mt","2E_antiRpt_mt_xFF"); // 3M_antiRpt_mt_xFF 3M_antiDphi_antiRpt_mt_xFF

  if (era.Contains("run2")) {
    aloader2E1.load("2E_mt","2E_antiRpt_mt_xFF","20UL16APV");
    aloader2E2.load("2E_mt","2E_antiRpt_mt_xFF","20UL17");
    aloader2E3.load("2E_mt","2E_antiRpt_mt_xFF","20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  //aloader2E.preparecard("MM2E_"+era+"_datacard.root","mergedMu2E");
  aloader2E.compare(p1,1,1.,p2); // ff
  SaveAs(canvas_2,"MMFF_2E_mt.pdf",p1);
  //aloader2E.close();

  return;

  aloader2E.load("2E_CRdphi_mt","2E_CRdphi_antiRpt_mt_xFF");

  if (era.Contains("run2")) {
    aloader2E1.load("2E_CRdphi_mt","2E_CRdphi_antiRpt_mt_xFF","20UL16APV");
    aloader2E2.load("2E_CRdphi_mt","2E_CRdphi_antiRpt_mt_xFF","20UL17");
    aloader2E3.load("2E_CRdphi_mt","2E_CRdphi_antiRpt_mt_xFF","20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,10,1.,p2); // ff
  SaveAs(canvas_2,"MMFF_2E_CRdphi_mt.pdf",p1);

  return;
}
