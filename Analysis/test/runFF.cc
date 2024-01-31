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

void runFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static constexpr double WZxsec_ = 5.213; // 0.65*62.78;
  static constexpr double ZZxsec_ = 13.81;
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

  static TString anlyzrMC = "mergedEleCRanalyzer"+postfix;
  static TString anlyzrData = "mergedEleCRanalyzerData";

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

  TFile* datafile = new TFile("EleAnalyzer_"+firstEra+"_data.root","READ");
  TFile* WZfile = new TFile("EleAnalyzer_"+firstEra+"_WZFXFX.root","READ");
  TFile* ZZfile = new TFile("EleAnalyzer_"+firstEra+"_ZZ.root","READ");

  TFile *datafile1, *WZfile1, *ZZfile1;
  TFile *datafile2, *WZfile2, *ZZfile2;
  TFile *datafile3, *WZfile3, *ZZfile3;

  if (era.Contains("run2")) {
    datafile1 = new TFile("EleAnalyzer_20UL16APV_data.root","READ");
    WZfile1 = new TFile("EleAnalyzer_20UL16APV_WZ.root","READ");
    ZZfile1 = new TFile("EleAnalyzer_20UL16APV_ZZ.root","READ");

    datafile2 = new TFile("EleAnalyzer_20UL17_data.root","READ");
    WZfile2 = new TFile("EleAnalyzer_20UL17_WZ.root","READ");
    ZZfile2 = new TFile("EleAnalyzer_20UL17_ZZ.root","READ");

    datafile3 = new TFile("EleAnalyzer_20UL18_data.root","READ");
    WZfile3 = new TFile("EleAnalyzer_20UL18_WZ.root","READ");
    ZZfile3 = new TFile("EleAnalyzer_20UL18_ZZ.root","READ");
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
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("EleAnalyzer_"+firstEra+"_H2000A750.root","READ"),"H2000A750")
  };

  std::vector<SigSample> sigsamples3E = sigsamples; // {SigSample(new TFile("EleAnalyzer_"+firstEra+"_H750A5.root","READ"),"H750A5")};

  std::vector<SigSample> sigsamples1, sigsamples2, sigsamples3;
  std::vector<SigSample> sigsamples3E1, sigsamples3E2, sigsamples3E3;

  if (era.Contains("run2")) {
    sigsamples1 = {
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("EleAnalyzer_20UL16APV_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples2 = {
      SigSample(new TFile("EleAnalyzer_20UL17_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("EleAnalyzer_20UL17_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("EleAnalyzer_20UL17_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("EleAnalyzer_20UL17_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples3 = {
      SigSample(new TFile("EleAnalyzer_20UL18_H250A1.root","READ"),"H250A1"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A2.root","READ"),"H250A2"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A5.root","READ"),"H250A5"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A10.root","READ"),"H250A10"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A50.root","READ"),"H250A50"),
      SigSample(new TFile("EleAnalyzer_20UL18_H250A100.root","READ"),"H250A100"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A1.root","READ"),"H750A1"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A2.root","READ"),"H750A2"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A5.root","READ"),"H750A5"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A10.root","READ"),"H750A10"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A50.root","READ"),"H750A50"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A100.root","READ"),"H750A100"),
      SigSample(new TFile("EleAnalyzer_20UL18_H750A250.root","READ"),"H750A250"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A1.root","READ"),"H2000A1"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A2.root","READ"),"H2000A2"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A5.root","READ"),"H2000A5"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A10.root","READ"),"H2000A10"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A50.root","READ"),"H2000A50"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A100.root","READ"),"H2000A100"),
      SigSample(new TFile("EleAnalyzer_20UL18_H2000A750.root","READ"),"H2000A750")
    };

    sigsamples3E1 = sigsamples1; // {SigSample(new TFile("EleAnalyzer_20UL16APV_H750A5.root","READ"),"H750A5")};
    sigsamples3E2 = sigsamples2; // {SigSample(new TFile("EleAnalyzer_20UL17_H750A5.root","READ"),"H750A5")};
    sigsamples3E3 = sigsamples3; // {SigSample(new TFile("EleAnalyzer_20UL18_H750A5.root","READ"),"H750A5")};
  }

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
  };

  class HistLoader2E : public HistLoaderBase {
  public:
    HistLoader2E(TFile* adatafile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile),
      sigFiles_(sigFiles) {}

    ~HistLoader2E() {
      delete dataHist_, FFHist_;

      if (dataHist2_) {
        delete dataHist2_, FFHist2_;
        dataHist2_ = nullptr;
      }

    }

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* FFHist_ = nullptr;

    TH1D* dataHist2_ = nullptr;
    TH1D* FFHist2_ = nullptr;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<std::map<std::string,SystVariation>> sigSyst_;

  public:
    void add(const HistLoader2E& other) {
      this->dataHist_->Add(other.dataHist_);
      this->FFHist_->Add(other.FFHist_);

      if (dataHist2_) {
        this->dataHist2_->Add(other.dataHist2_);
        this->FFHist2_->Add(other.FFHist2_);
      }

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

    void load(const std::string& nameNum, const std::string& name, int color, std::string anlyzrEra="", const std::string& nameNum2="", const std::string& name2="", int color2=0) {
      if (dataHist_)
        delete dataHist_, FFHist_;

      if (dataHist2_) {
        delete dataHist2_, FFHist2_;
        dataHist2_ = nullptr;
      }

      if (!syst_.empty())
        syst_.clear();

      dataHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum).c_str() )->Clone();
      FFHist_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name).c_str() )->Clone();

      FFHist_->SetLineWidth(0);
      FFHist_->SetFillColor(color);

      if ( TString(name).Contains("invM") || TString(name).Contains("eta") ) {
        TH1D* FFHistUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_up").c_str() )->Clone();
        TH1D* FFHistDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_dn").c_str() )->Clone();
        syst_["FFHist"] = SystVariation(FFHistUp,FFHistDn);
      }

      const double lumi = retrieveLumi(anlyzrEra);
      const double sigLumi = 0.001;
      sigHist_.clear();
      sigSyst_.clear();

      if (nameNum2!="") {
        dataHist2_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum2).c_str() )->Clone();
        FFHist2_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2).c_str() )->Clone();
        FFHist2_->SetLineWidth(0);
        FFHist2_->SetFillColor(color2);

        if ( TString(name2).Contains("invM") ) {
          TH1D* FFHist2Up = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_up").c_str() )->Clone();
          TH1D* FFHist2Dn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_dn").c_str() )->Clone();
          syst_["FFHist2"] = SystVariation(FFHist2Up,FFHist2Dn);
        }

        auto retrieveSigHist = [this,&nameNum,&nameNum2,&lumi,&sigLumi,&anlyzrEra] (TFile* afile, const std::string& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( (std::string("mergedEleCRanalyzer"+anlyzrEra+"/")+nameNum+systName).c_str() )->Clone();
          ahist->Add( (TH1D*)afile->Get( (std::string("mergedEleCRanalyzer"+anlyzrEra+"/")+nameNum2+systName).c_str() ) );
          ahist->Scale( lumi*1000.*sigLumi / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned idx=0; idx<sigFiles_.size(); idx++) {
          TH1D* asigHist = retrieveSigHist(sigFiles_.at(idx).file_,"");
          sigHist_.push_back( asigHist );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          if ( TString(name2).Contains("invM") ) {
            std::map<std::string,SystVariation> init;

            TH1D* sigModHeepUp = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdUp");
            TH1D* sigModHeepDn = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdDn");
            init["sigModHeep"] = SystVariation(sigModHeepUp,sigModHeepDn);
            TH1D* sigMergedEleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdUp");
            TH1D* sigMergedEleDn = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdDn");
            init["sigMergedEle"] = SystVariation(sigMergedEleUp,sigMergedEleDn);

            if ( TString(nameNum2).Contains("CRME") ) {
              TH1D* sigScaleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleScale");
              TH1D* sigScaleDn = variateDn(asigHist,sigScaleUp);
              init["sigScale"] = SystVariation(sigScaleUp,sigScaleDn);
              TH1D* sigSmearUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleSmear");
              TH1D* sigSmearDn = variateDn(asigHist,sigSmearUp);
              init["sigSmear"] = SystVariation(sigSmearUp,sigSmearDn);
            }

            sigSyst_.push_back(init);
          }
        }
      }
    }

    void compare(TPad* padUp, bool removeZero=false, TPad* padDn=nullptr, int rebin=1) {
      if ( FFHist_->GetNbinsX() % dataHist_->GetNbinsX()!=0 ) {
        // hardcode (rebin to GCD)
        FFHist_->Rebin( 5 );
        dataHist_->Rebin( 2 );

        if (dataHist2_) {
          FFHist2_->Rebin( 5 );
          dataHist2_->Rebin( 2 );
        }

        if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( 5 );
            element.second.dn_->Rebin( 5 );
          }
        }
      } else {
        FFHist_->Rebin( FFHist_->GetNbinsX()/dataHist_->GetNbinsX() );

        if (dataHist2_)
          FFHist2_->Rebin( FFHist2_->GetNbinsX()/dataHist2_->GetNbinsX() );

        if (!syst_.empty()) {
          for (const auto& element : syst_) {
            element.second.up_->Rebin( element.second.up_->GetNbinsX()/dataHist_->GetNbinsX() );
            element.second.dn_->Rebin( element.second.dn_->GetNbinsX()/dataHist_->GetNbinsX() );
          }
        }
      }

      dataHist_->Rebin(rebin);
      FFHist_->Rebin(rebin);

      if (dataHist2_) {
        dataHist2_->Rebin(rebin);
        FFHist2_->Rebin(rebin);
      }

      if (!syst_.empty()) {
        for (const auto& element : syst_) {
          element.second.up_->Rebin( rebin );
          element.second.dn_->Rebin( rebin );
        }
      }

      TH1D* dataHist = dataHist_;
      TH1D* FFadded = (TH1D*)FFHist_->Clone();
      TNamed* FFobj = FFHist_;

      if (dataHist2_) {
        dataHist->Add(dataHist2_);
        FFadded->Add(FFHist2_);
        auto* stack = new THStack("stack",";GeV;");
        stack->Add(FFHist_);
        stack->Add(FFHist2_);
        FFobj = stack;
      }

      if ( TString(dataHist_->GetName()).Contains("invM") && TString(dataHist_->GetName()).Contains("mixedME") )
        dataHist->GetXaxis()->SetRangeUser(0.,1000.);
      else if ( TString(dataHist_->GetName()).Contains("invM") && TString(dataHist_->GetName()).Contains("CRME") ) {
        if (dataHist2_) {
          dataHist->GetXaxis()->SetRangeUser(0.,1000.);
          dataHist->SetMinimum( 1.1 );
        } else
          dataHist->GetXaxis()->SetRangeUser(0.,250.);
      }

      padUp->cd();

      if ( padUp->GetLogy() )
        dataHist->SetMaximum( 10.*std::max(dataHist_->GetMaximum(),FFadded->GetMaximum()) );
      else
        dataHist->SetMaximum( 1.2*std::max(dataHist_->GetMaximum(),FFadded->GetMaximum()) );

      if (removeZero)
        dataHist->SetMinimum( 0.001 );

      dataHist->SetLineWidth(2);
      dataHist->SetLineColor(kBlack);
      dataHist->Draw("E1");
      FFobj->Draw("hist&same");
      dataHist->Draw("E1&same");

      if (dataHist2_) {
        TLegend* legend = new TLegend(0.55,0.65,0.95,0.93);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist,"Data");

        legend->AddEntry(FFHist_,"Nonprompt bkg");
        legend->AddEntry(FFHist2_,"Prompt bkg");

        if ( TString(dataHist_->GetName()).Contains("CRME") && TString(dataHist_->GetName()).Contains("invM") ) {
          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            sigHist_.at(idx)->Rebin(rebin);

            for (const auto& element : sigSyst_.at(idx)) {
              sigSyst_.at(idx).at(element.first).up_->Rebin(rebin);
              sigSyst_.at(idx).at(element.first).dn_->Rebin(rebin);
            }

            sigHist_.at(idx)->Draw("hist&same");
          }

          //legend->AddEntry(sigHist_.at(3),"H250/750/2000A1");
          legend->AddEntry(sigHist_.at(0),"X750Y1");
        }

        legend->Draw();
      }

      TH1D* ratio = nullptr;

      if (padDn) {
        ratio = (TH1D*)dataHist->Clone();
        ratio->SetStats(0);
        ratio->SetTitle("");
        ratio->Divide(FFadded);
        ratio->GetYaxis()->SetTitle("Obs/Exp");
        ratio->GetYaxis()->SetTitleSize(0.1);
        ratio->GetYaxis()->SetTitleOffset(0.4);
        ratio->GetXaxis()->SetLabelSize(0.1);
        ratio->GetYaxis()->SetLabelSize(0.1);
        ratio->GetXaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetLabelOffset(0.01);
        ratio->GetYaxis()->SetRangeUser(0.2,1.8);
        ratio->GetYaxis()->SetNdivisions(505);
        ratio->GetXaxis()->SetTitleSize(0.12);
        ratio->GetXaxis()->SetTitleOffset(0.75);
        ratio->SetLineColor(kBlack);

        padDn->cd();
        ratio->Draw("E1");
      }

      if (!syst_.empty()) {
        std::vector<double> x0, y0, errx, erryDn, erryUp;
        std::vector<double> r0, errRdn, errRup;

        TH1D* FFHistUpAdded = (TH1D*)syst_.at("FFHist").up_->Clone();
        TH1D* FFHistDnAdded = (TH1D*)syst_.at("FFHist").dn_->Clone();
        TH1D* FFHist2UpAdded = nullptr;
        TH1D* FFHist2DnAdded = nullptr;

        if (syst_.count("FFHist2")) {
          FFHistUpAdded->Add(FFHist2_);
          FFHistDnAdded->Add(FFHist2_);

          FFHist2UpAdded = (TH1D*)syst_.at("FFHist2").up_->Clone();
          FFHist2DnAdded = (TH1D*)syst_.at("FFHist2").dn_->Clone();
          FFHist2UpAdded->Add(FFHist_);
          FFHist2DnAdded->Add(FFHist_);
        }

        for (unsigned idx = 1; idx <= FFadded->GetNbinsX(); idx++) {
          x0.push_back(FFadded->GetBinCenter(idx));
          y0.push_back(FFadded->GetBinContent(idx));
          errx.push_back(FFadded->GetBinWidth(idx)/2.);

          double valFFUp = syst_.at("FFHist").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
          double valFFDn = FFHist_->GetBinContent(idx) - syst_.at("FFHist").dn_->GetBinContent(idx);
          double valFF2Up = 0.;
          double valFF2Dn = 0.;

          if (syst_.count("FFHist2")) {
            valFF2Up = syst_.at("FFHist2").up_->GetBinContent(idx) - FFHist2_->GetBinContent(idx);
            valFF2Dn = FFHist2_->GetBinContent(idx) - syst_.at("FFHist2").dn_->GetBinContent(idx);
          }

          erryUp.push_back( std::sqrt( valFFUp*valFFUp + valFF2Up*valFF2Up ) );
          erryDn.push_back( std::sqrt( valFFDn*valFFDn + valFF2Dn*valFF2Dn ) );

          if (ratio) {
            r0.push_back(1.);

            double rvalFFUp = FFHistUpAdded->GetBinContent(idx)/FFadded->GetBinContent(idx) - 1.;
            double rvalFFDn = 1. - FFHistDnAdded->GetBinContent(idx)/FFadded->GetBinContent(idx);
            double rvalFF2Up = 0.;
            double rvalFF2Dn = 0.;

            if (syst_.count("FFHist2")) {
              rvalFF2Up = FFHist2UpAdded->GetBinContent(idx)/FFadded->GetBinContent(idx) - 1.;
              rvalFF2Dn = 1. - FFHist2DnAdded->GetBinContent(idx)/FFadded->GetBinContent(idx);
            }

            errRup.push_back( std::sqrt( rvalFFUp*rvalFFUp + rvalFF2Up*rvalFF2Up ) );
            errRdn.push_back( std::sqrt( rvalFFDn*rvalFFDn + rvalFF2Dn*rvalFF2Dn ) );
          }
        }

        auto gr = new TGraphAsymmErrors(dataHist->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
        gr->SetFillColor(kGray+2);
        gr->SetLineColor(kGray+2);
        gr->SetFillStyle(3004);
        padUp->cd();
        gr->Draw("2");

        if (padDn) {
          auto rgr = new TGraphAsymmErrors(dataHist->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
          rgr->SetFillColor(kGray+2);
          rgr->SetLineColor(kGray+2);
          rgr->SetFillStyle(3004);
          padDn->cd();
          rgr->Draw("2");
        }

        if (dir_) {
          dir_->WriteTObject(dataHist_,"data_obs");
          dir_->WriteTObject(FFHist_,"SS");
          dir_->WriteTObject(FFHist2_,"OS");

          dir_->WriteTObject(syst_.at("FFHist").up_,"SS_mergedEleFakeFactorSSUp");
          dir_->WriteTObject(syst_.at("FFHist").dn_,"SS_mergedEleFakeFactorSSDown");
          dir_->WriteTObject(syst_.at("FFHist2").up_,"OS_mergedEleFakeFactorOSUp");
          dir_->WriteTObject(syst_.at("FFHist2").dn_,"OS_mergedEleFakeFactorOSDown");

          for (unsigned idx=0; idx<sigHist_.size(); idx++) {
            dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
            dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").up_,sigFiles_.at(idx).name_+"_modHeepIdUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").dn_,sigFiles_.at(idx).name_+"_modHeepIdDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigMergedEle").up_,sigFiles_.at(idx).name_+"_mergedEleIdUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigMergedEle").dn_,sigFiles_.at(idx).name_+"_mergedEleIdDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigScale").up_,sigFiles_.at(idx).name_+"_mergedEleEnScaleUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigScale").dn_,sigFiles_.at(idx).name_+"_mergedEleEnScaleDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigSmear").up_,sigFiles_.at(idx).name_+"_mergedEleEnSmearUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigSmear").dn_,sigFiles_.at(idx).name_+"_mergedEleEnSmearDown");
          }
        }
      }
    }
  };

  // invM
  auto aloader2E = HistLoader2E(datafile,sigsamples);
  aloader2E.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,postfix.Data());

  auto aloader2E1 = HistLoader2E(datafile1,sigsamples1);
  auto aloader2E2 = HistLoader2E(datafile2,sigsamples2);
  auto aloader2E3 = HistLoader2E(datafile3,sigsamples3);

  if (era.Contains("run2")) {
    aloader2E1.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL16APV");
    aloader2E2.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL17"); 
    aloader2E3.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,true,p2,4);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiSS.png",p1);

  aloader2E.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,"20UL16APV");
    aloader2E2.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,"20UL17");
    aloader2E3.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_mixedOS.png",p1);

  p1->SetLogy();
  aloader2E.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,"20UL16APV");
    aloader2E2.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,"20UL17");
    aloader2E3.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,false,p2,4);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiOS.png",p1);

  //p1->SetLogx();
  //p2->SetLogx();
  p1->SetLogy();
  aloader2E.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,postfix.Data(),"2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange);

  if (era.Contains("run2")) {
    aloader2E1.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL16APV","2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange);
    aloader2E2.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL17","2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange);
    aloader2E3.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL18","2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange);
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.preparecard("ME2E_"+era+"_datacard.root","mergedEle2E");
  aloader2E.compare(p1,false,p2,2);
  SaveAs(canvas_2,"MEFF_2E_invM_denom_CR.eps",p1);
  aloader2E.close();

  p1->SetLogy();
  aloader2E.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,postfix.Data(),"2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange);

  if (era.Contains("run2")) {
    aloader2E1.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL16APV","2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange);
    aloader2E2.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL17","2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange);
    aloader2E3.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL18","2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange);
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,false,p2,4);
  SaveAs(canvas_2,"MEFF_2E_invM_mixed_CR.eps",p1);

  //p1->SetLogx(0);
  //p2->SetLogx(0);

  aloader2E.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL16APV");
    aloader2E2.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL17");
    aloader2E3.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedSS.png");

  // Et
  canvas_1->SetLogy();
  aloader2E.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,"20UL16APV");
    aloader2E2.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,"20UL17");
    aloader2E3.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedSS.png");

  aloader2E.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,"20UL16APV");
    aloader2E2.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,"20UL17");
    aloader2E3.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiSS.png",p1);

  aloader2E.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,"20UL16APV");
    aloader2E2.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,"20UL17");
    aloader2E3.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedOS.png");

  aloader2E.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,"20UL16APV");
    aloader2E2.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,"20UL17");
    aloader2E3.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiOS.png",p1);
  canvas_1->SetLogy(0);

  // eta
  p1->SetLogy(0);
  aloader2E.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,"20UL16APV");
    aloader2E2.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,"20UL17");
    aloader2E3.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedSS.png");

  aloader2E.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,"20UL16APV");
    aloader2E2.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,"20UL17");
    aloader2E3.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,true,p2,2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiSS.png",p1);

  aloader2E.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,"20UL16APV");
    aloader2E2.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,"20UL17");
    aloader2E3.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedOS.png");

  aloader2E.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,postfix.Data());

  if (era.Contains("run2")) {
    aloader2E1.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,"20UL16APV");
    aloader2E2.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,"20UL17");
    aloader2E3.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,"20UL18");
    aloader2E.add(aloader2E1);
    aloader2E.add(aloader2E2);
    aloader2E.add(aloader2E3);
  }

  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiOS.png",p1);

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
      void add(const denomHists& other) {
        this->SSFFHist_->Add(other.SSFFHist_);
        this->SSWZHist_->Add(other.SSWZHist_);
        this->SSZZHist_->Add(other.SSZZHist_);
        this->OSWZHist_->Add(other.OSWZHist_);
        this->OSZZHist_->Add(other.OSZZHist_);
      }

      void load(TString denomNameSS, TString denomNameOS, TFile* datafile, TFile* WZfile, TFile* ZZfile, std::string anlyzrEra="") {
        SSFFHist_ = (TH1D*)datafile->Get(anlyzrData+"/"+denomNameSS)->Clone();

        const double WZsumwgt = ((TH1D*)WZfile->Get("evtCounter/h_sumW"))->GetBinContent(1);
        const double ZZsumwgt = ((TH1D*)ZZfile->Get("evtCounter/h_sumW"))->GetBinContent(1);
        const double lumi = retrieveLumi(anlyzrEra);

        SSWZHist_ = (TH1D*)WZfile->Get("mergedEleCRanalyzer"+anlyzrEra+"/"+denomNameSS)->Clone();
        SSWZHist_->Scale( lumi*1000.*WZxsec_/WZsumwgt );
        SSZZHist_ = (TH1D*)ZZfile->Get("mergedEleCRanalyzer"+anlyzrEra+"/"+denomNameSS)->Clone();
        SSZZHist_->Scale( lumi*1000.*ZZxsec_/ZZsumwgt );
        OSWZHist_ = (TH1D*)WZfile->Get("mergedEleCRanalyzer"+anlyzrEra+"/"+denomNameOS)->Clone();
        OSWZHist_->Scale( lumi*1000.*WZxsec_/WZsumwgt );
        OSZZHist_ = (TH1D*)ZZfile->Get("mergedEleCRanalyzer"+anlyzrEra+"/"+denomNameOS)->Clone();
        OSZZHist_->Scale( lumi*1000.*ZZxsec_/ZZsumwgt );
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
    std::vector<std::map<std::string,SystVariation>> sigSyst_;

  public:
    void add(const HistLoader3E& other) {
      this->dataHist_->Add(other.dataHist_);
      this->nominal_.add(other.nominal_);
      this->SSFFup_.add(other.SSFFup_);
      this->SSFFdn_.add(other.SSFFdn_);
      this->OSFFup_.add(other.OSFFup_);
      this->OSFFdn_.add(other.OSFFdn_);

      if (!sigSyst_.empty()) {
        this->heepIdUp_.add(other.heepIdUp_);
        this->heepIdDn_.add(other.heepIdDn_);
      }

      if (!sigHist_.empty()) {
        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          this->sigHist_.at(idx)->Add(other.sigHist_.at(idx));

          if (!sigSyst_.empty()) {
            for (const auto& element : this->sigSyst_.at(idx)) {
              this->sigSyst_.at(idx).at(element.first).up_->Add(other.sigSyst_.at(idx).at(element.first).up_);
              this->sigSyst_.at(idx).at(element.first).dn_->Add(other.sigSyst_.at(idx).at(element.first).dn_);
            }
          }
        }
      }
    }

    void load(TString numName, TString denomName, std::string anlyzrEra="") {
      if (dataHist_) {
        delete dataHist_;
      }

      dataHist_ = (TH1D*)datafile_->Get(anlyzrData+"/"+numName)->Clone();
      nominal_.load(denomName+"_xSSFF",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
      SSFFup_.load(denomName+"_xSSFF_up",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
      SSFFdn_.load(denomName+"_xSSFF_dn",denomName+"_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
      OSFFup_.load(denomName+"_xSSFF",denomName+"_xOSFF_up",datafile_,WZfile_,ZZfile_,anlyzrEra);
      OSFFdn_.load(denomName+"_xSSFF",denomName+"_xOSFF_dn",datafile_,WZfile_,ZZfile_,anlyzrEra);

      sigHist_.clear();
      sigSyst_.clear();

      if (numName.Contains("invM")) {
        const double lumi = retrieveLumi(anlyzrEra);
        const double sigLumi = 0.001;

        auto retrieveSigHist = [this,&numName,&anlyzrEra,&lumi,&sigLumi] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( "mergedEleCRanalyzer"+anlyzrEra+"/"+numName+systName )->Clone();
          ahist->Scale( lumi*1000.*sigLumi / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned idx = 0; idx < sigFiles_.size(); idx++) {
          TH1D* asigHist = retrieveSigHist(sigFiles_.at(idx).file_,"");
          sigHist_.push_back( asigHist );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          std::map<std::string,SystVariation> init;

          TH1D* sigModHeepUp = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdUp");
          TH1D* sigModHeepDn = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdDn");
          init["sigModHeep"] = SystVariation(sigModHeepUp,sigModHeepDn);
          TH1D* sigMergedEleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdUp");
          TH1D* sigMergedEleDn = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdDn");
          init["sigMergedEle"] = SystVariation(sigMergedEleUp,sigMergedEleDn);

          TH1D* sigScaleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleScale");
          TH1D* sigScaleDn = variateDn(asigHist,sigScaleUp);
          init["sigScale"] = SystVariation(sigScaleUp,sigScaleDn);
          TH1D* sigSmearUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleSmear");
          TH1D* sigSmearDn = variateDn(asigHist,sigSmearUp);
          init["sigSmear"] = SystVariation(sigSmearUp,sigSmearDn);

          TH1D* sigResolvedScaleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_scaleUp");
          TH1D* sigResolvedScaleDn = retrieveSigHist(sigFiles_.at(idx).file_,"_scaleDn");
          init["sigResolvedEnScale"] = SystVariation(sigResolvedScaleUp,sigResolvedScaleDn);

          TH1D* sigResolvedSigmaUp = retrieveSigHist(sigFiles_.at(idx).file_,"_sigmaUp");
          TH1D* sigResolvedSigmaDn = retrieveSigHist(sigFiles_.at(idx).file_,"_sigmaDn");
          init["sigResolvedEnSigma"] = SystVariation(sigResolvedSigmaUp,sigResolvedSigmaDn);

          sigSyst_.push_back(init);
        }

        heepIdUp_.load(denomName+"_xSSFF_heepIdUp",denomName+"_xOSFF_heepIdUp",datafile_,WZfile_,ZZfile_,anlyzrEra);
        heepIdDn_.load(denomName+"_xSSFF_heepIdDn",denomName+"_xOSFF_heepIdDn",datafile_,WZfile_,ZZfile_,anlyzrEra);
      }
    }

    void compare(TPad* pad, int rebin=1) {
      if (rebin!=1) {
        dataHist_->Rebin(rebin);
        nominal_.rebin(rebin);
        SSFFup_.rebin(rebin);
        SSFFdn_.rebin(rebin);
        OSFFup_.rebin(rebin);
        OSFFdn_.rebin(rebin);

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);

          for (const auto& element : sigSyst_.at(idx)) {
            element.second.up_->Rebin( rebin );
            element.second.dn_->Rebin( rebin );
          }
        }
      }

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

        if (pad->GetLogy())
          dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      } else if (TString(dataHist_->GetName()).Contains("Et")) {
        dataHist_->GetXaxis()->SetRangeUser(50.,300.);
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
        legend->AddEntry(denomSSfinal.get(),"Nonprompt bkg");
        legend->AddEntry(denomOSfinal.get(),"Prompt bkg");

        if (TString(dataHist_->GetName()).Contains("invM")) {
          legend->AddEntry(sigHist_.at(0),"X750A5");
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
          dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").up_,sigFiles_.at(idx).name_+"_modHeepIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").dn_,sigFiles_.at(idx).name_+"_modHeepIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMergedEle").up_,sigFiles_.at(idx).name_+"_mergedEleIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMergedEle").dn_,sigFiles_.at(idx).name_+"_mergedEleIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigScale").up_,sigFiles_.at(idx).name_+"_mergedEleEnScaleUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigScale").dn_,sigFiles_.at(idx).name_+"_mergedEleEnScaleDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigSmear").up_,sigFiles_.at(idx).name_+"_mergedEleEnSmearUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigSmear").dn_,sigFiles_.at(idx).name_+"_mergedEleEnSmearDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigResolvedEnScale").up_,sigFiles_.at(idx).name_+"_elEnergyScaleUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigResolvedEnScale").dn_,sigFiles_.at(idx).name_+"_elEnergyScaleDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigResolvedEnSigma").up_,sigFiles_.at(idx).name_+"_elEnergySigmaUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigResolvedEnSigma").dn_,sigFiles_.at(idx).name_+"_elEnergySigmaDown");
        }
      }

      denomSSfinal.release();
      denomOSfinal.release();
    }
  };

  canvas_1->cd();

  auto aloader3E = HistLoader3E(datafile,WZfile,ZZfile,sigsamples3E);
  auto aloader3E1 = HistLoader3E(datafile1,WZfile1,ZZfile1,sigsamples3E1);
  auto aloader3E2 = HistLoader3E(datafile2,WZfile2,ZZfile2,sigsamples3E2);
  auto aloader3E3 = HistLoader3E(datafile3,WZfile3,ZZfile3,sigsamples3E3);

  canvas_1->SetLogy();
  aloader3E.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR",postfix.Data());

  if (era.Contains("run2")) {
    aloader3E1.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR","20UL16APV");
    aloader3E2.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR","20UL17");
    aloader3E3.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR","20UL18");
    aloader3E.add(aloader3E1);
    aloader3E.add(aloader3E2);
    aloader3E.add(aloader3E3);
  }

  aloader3E.preparecard("ME3E_"+era+"_datacard.root","mergedEle3E");
  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"MEFF_3E_invM.eps");
  aloader3E.close();
  canvas_1->SetLogy(0);

  canvas_1->SetLogy();
  aloader3E.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME",postfix.Data());

  if (era.Contains("run2")) {
    aloader3E1.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME","20UL16APV");
    aloader3E2.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME","20UL17");
    aloader3E3.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME","20UL18");
    aloader3E.add(aloader3E1);
    aloader3E.add(aloader3E2);
    aloader3E.add(aloader3E3);
  }

  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"FF_3E_Et.png");
  canvas_1->SetLogy(0);

  aloader3E.load("3E_eta_CR_EB_CRME","3E_antiME_eta",postfix.Data());

  if (era.Contains("run2")) {
    aloader3E1.load("3E_eta_CR_EB_CRME","3E_antiME_eta","20UL16APV");
    aloader3E2.load("3E_eta_CR_EB_CRME","3E_antiME_eta","20UL17");
    aloader3E3.load("3E_eta_CR_EB_CRME","3E_antiME_eta","20UL18");
    aloader3E.add(aloader3E1);
    aloader3E.add(aloader3E2);
    aloader3E.add(aloader3E3);
  }

  aloader3E.compare(canvas_1,2);
  SaveAs(canvas_1,"FF_3E_eta.png");

  return;
}
