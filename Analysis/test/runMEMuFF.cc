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

void runMEMuFF(TString era) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Work in progress";  // default extra text is "Preliminary"

  static constexpr double WZxsec_ = 5.213; // 0.65*62.78;
  static constexpr double ZZxsec_ = 13.81;
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

  TFile* datafile = new TFile("MuAnalyzer_"+fname+"_data.root","READ");
  TFile* WZfile = new TFile("MuAnalyzer_"+fname+"_WZFXFX.root","READ");
  TFile* ZZfile = new TFile("MuAnalyzer_"+fname+"_ZZ.root","READ");

  TFile* datafile1 = new TFile("MuAnalyzer_20UL16APV_data.root","READ");
  TFile* WZfile1 = new TFile("MuAnalyzer_20UL16APV_WZFXFX.root","READ");
  TFile* ZZfile1 = new TFile("MuAnalyzer_20UL16APV_ZZ.root","READ");

  TFile* datafile2 = new TFile("MuAnalyzer_20UL17_data.root","READ");
  TFile* WZfile2 = new TFile("MuAnalyzer_20UL17_WZFXFX.root","READ");
  TFile* ZZfile2 = new TFile("MuAnalyzer_20UL17_ZZ.root","READ");

  TFile* datafile3 = new TFile("MuAnalyzer_20UL18_data.root","READ");
  TFile* WZfile3 = new TFile("MuAnalyzer_20UL18_WZFXFX.root","READ");
  TFile* ZZfile3 = new TFile("MuAnalyzer_20UL18_ZZ.root","READ");

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

  /*auto sigsamples = std::vector<SigSample>{
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A750.root","READ"),"H2000A750")
  };

  auto sigsamples1 = std::vector<SigSample>{
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

  auto sigsamples2 = std::vector<SigSample>{
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

  auto sigsamples3 = std::vector<SigSample>{
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
  };*/

  auto sigsamples = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_"+fname+"_H750A1.root","READ"),"H750A1")};
  auto sigsamples1 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL16APV_H750A1.root","READ"),"H750A1")};
  auto sigsamples2 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL17_H750A1.root","READ"),"H750A1")};
  auto sigsamples3 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL18_H750A1.root","READ"),"H750A1")};

  auto sigsamples4 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_"+fname+"_H2000A1.root","READ"),"H2000A1")};
  auto sigsamples5 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL16APV_H2000A1.root","READ"),"H2000A1")};
  auto sigsamples6 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL17_H2000A1.root","READ"),"H2000A1")};
  auto sigsamples7 = std::vector<SigSample>{SigSample(new TFile("MuAnalyzer_20UL18_H2000A1.root","READ"),"H2000A1")};

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



/*  TH1D* nm = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_FF_CRME_MET_pt")->Clone();
  TH1D* dm = (TH1D*)datafile->Get("mergedEMuCRanalyzerData/1M_FF_CRME_antiDphi_MET_pt")->Clone();

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
  };

  class HistLoaderMM : public HistLoaderBase {
  public:
    HistLoaderMM(TFile* adatafile, std::vector<SigSample> sigsamples)
    : HistLoaderBase(adatafile,nullptr,nullptr), sigFiles_(sigsamples) {}

    ~HistLoaderMM()=default;

  private:
    TH1D* dataHist_ = nullptr;
    TH1D* FFHist_ = nullptr;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<std::map<std::string,SystVariation>> sigSyst_;

  public:
    void add(const HistLoaderMM& other) {
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
      FFHist_->SetLineWidth(0);

      sigHist_.clear();
      sigSyst_.clear();

      if (numName.Contains("CRME_mt") || numName.Contains("CRME_CRdphi_mt")) {
        TH1D* FFHistUp = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_up" )->Clone();
        TH1D* FFHistDn = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_dn" )->Clone();
        syst_["FFHist"] = SystVariation(FFHistUp,FFHistDn);
        TH1D* JESup = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_JESup" )->Clone();
        TH1D* JESdn = (TH1D*)datafile_->Get( std::string(anlyzrData+"/")+denomName+"_JESdn" )->Clone();
        syst_["JES"] = SystVariation(JESup,JESdn);
      }

      if (numName.Contains("CRME_mt")) {
        const double lumi = retrieveLumi(anlyzrEra);
        const double sigLumi = 0.01;

        auto retrieveSigHist = [this,&numName,&anlyzrEra,&lumi,&sigLumi] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( "mergedEMuCRanalyzer"+anlyzrEra+"/"+numName+systName )->Clone();
          ahist->Scale( lumi*1000.*sigLumi / ( (TH1D*)afile->Get( std::string("evtCounter/h_sumW").c_str() ) )->GetBinContent(1) );

          return ahist;
        };

        for (unsigned idx = 0; idx < sigFiles_.size(); idx++) {
          sigHist_.push_back( retrieveSigHist(sigFiles_.at(idx).file_,"") );
          sigHist_.back()->SetLineWidth(2);
          sigHist_.back()->SetLineColor(kRed);

          std::map<std::string,SystVariation> init;

          TH1D* sigModHeepUp = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdUp");
          TH1D* sigModHeepDn = retrieveSigHist(sigFiles_.at(idx).file_,"_heepIdDn");
          init["sigModHeep"] = SystVariation(sigModHeepUp,sigModHeepDn);
          TH1D* sigMergedEleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdUp");
          TH1D* sigMergedEleDn = retrieveSigHist(sigFiles_.at(idx).file_,"_mergedEleIdDn");
          init["sigMergedEle"] = SystVariation(sigMergedEleUp,sigMergedEleDn);
          TH1D* sigJESup = retrieveSigHist(sigFiles_.at(idx).file_,"_JESup");
          TH1D* sigJESdn = retrieveSigHist(sigFiles_.at(idx).file_,"_JESdn");
          init["sigJES"] = SystVariation(sigJESup,sigJESdn);
          TH1D* sigJERup = retrieveSigHist(sigFiles_.at(idx).file_,"_JERup");
          TH1D* sigJERdn = retrieveSigHist(sigFiles_.at(idx).file_,"_JERdn");
          init["sigJER"] = SystVariation(sigJERup,sigJERdn);
          TH1D* sigMuIdUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muIdUp");
          TH1D* sigMuIdDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muIdDn");
          init["sigMuId"] = SystVariation(sigMuIdUp,sigMuIdDn);
          TH1D* sigMuIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muIsoUp");
          TH1D* sigMuIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muIsoDn");
          init["sigMuIso"] = SystVariation(sigMuIsoUp,sigMuIsoDn);
          TH1D* sigMuBoostIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoUp");
          TH1D* sigMuBoostIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoDn");
          init["sigMuBoostIso"] = SystVariation(sigMuBoostIsoUp,sigMuBoostIsoDn);
          TH1D* sigTrigUp = retrieveSigHist(sigFiles_.at(idx).file_,"_trigUp");
          TH1D* sigTrigDn = retrieveSigHist(sigFiles_.at(idx).file_,"_trigDn");
          init["sigTrig"] = SystVariation(sigTrigUp,sigTrigDn);
          TH1D* sigMuRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muRecoUp");
          TH1D* sigMuRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muRecoDn");
          init["sigMuReco"] = SystVariation(sigMuRecoUp,sigMuRecoDn);
          TH1D* sigElRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoUp");
          TH1D* sigElRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoDn");
          init["sigElReco"] = SystVariation(sigElRecoUp,sigElRecoDn);
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

    void compare(TPad* padUp, int rebin=1, TPad* padDn=nullptr) {
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
        std::vector<double> binEdges = {0.,25.,50.,75.,100.,125.,150.,175.,200.,225.,250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,550.,600.,700.,800.,1000.,2500.};

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

      padUp->cd();

      dataHist_->SetLineWidth(2);
      dataHist_->SetLineColor(kBlack);
      dataHist_->SetMinimum(0.001);
      dataHist_->SetMaximum( std::max(1.5*dataHist_->GetMaximum(),1.5*FFHist_->GetMaximum()) );

      if (TString(dataHist_->GetName()).Contains("CRME_mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,2500.);
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      }

      if (TString(dataHist_->GetName()).Contains("CRME_CRdphi_mt")) {
        dataHist_->GetXaxis()->SetRangeUser(0.,2500.);
      }

      dataHist_->Draw("E1");
      FFHist_->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("CRME_mt")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (TString(dataHist_->GetName()).Contains("CRME_mt") || TString(dataHist_->GetName()).Contains("CRME_CRdphi_mt")) {
        double ylow = TString(dataHist_->GetName()).Contains("CRME_CRdphi_mt") ? 0.75 : 0.68;
        TLegend* legend = new TLegend(0.75,ylow,0.93,0.93);
        legend->SetBorderSize(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(FFHist_,"Bkg");

        if (TString(dataHist_->GetName()).Contains("CRME_mt"))
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
        gr->SetFillColor(kGray+2);
        gr->SetLineColor(kGray+2);
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
          dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").up_,sigFiles_.at(idx).name_+"_modHeepIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").dn_,sigFiles_.at(idx).name_+"_modHeepIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMergedEle").up_,sigFiles_.at(idx).name_+"_mergedEleIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMergedEle").dn_,sigFiles_.at(idx).name_+"_mergedEleIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJES").up_,sigFiles_.at(idx).name_+"_JESUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJES").dn_,sigFiles_.at(idx).name_+"_JESDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJER").up_,sigFiles_.at(idx).name_+"_JERUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigJER").dn_,sigFiles_.at(idx).name_+"_JERDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").up_,sigFiles_.at(idx).name_+"_highPtIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").dn_,sigFiles_.at(idx).name_+"_highPtIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").up_,sigFiles_.at(idx).name_+"_muLooseIsoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").dn_,sigFiles_.at(idx).name_+"_muLooseIsoDown");
          // dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").up_,sigFiles_.at(idx).name_+"_muBoostIsoUp");
          // dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").dn_,sigFiles_.at(idx).name_+"_muBoostIsoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").up_,sigFiles_.at(idx).name_+"_muTrigUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").dn_,sigFiles_.at(idx).name_+"_muTrigDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").up_,sigFiles_.at(idx).name_+"_muRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").dn_,sigFiles_.at(idx).name_+"_muRecoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").up_,sigFiles_.at(idx).name_+"_elRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").dn_,sigFiles_.at(idx).name_+"_elRecoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").up_,sigFiles_.at(idx).name_+"_PUrwgtUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").dn_,sigFiles_.at(idx).name_+"_PUrwgtDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").up_,sigFiles_.at(idx).name_+"_prefireUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").dn_,sigFiles_.at(idx).name_+"_prefireDown");
        }
      }
    }
  }; // class

  auto aloaderMM = HistLoaderMM(datafile,sigsamples4);
  auto aloaderMM1 = HistLoaderMM(datafile1,sigsamples5);
  auto aloaderMM2 = HistLoaderMM(datafile2,sigsamples6);
  auto aloaderMM3 = HistLoaderMM(datafile3,sigsamples7);

  aloaderMM.load("1M_CRME_mt","1M_CRME_antiRpt_mt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_mt","1M_CRME_antiRpt_mt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_mt","1M_CRME_antiRpt_mt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_mt","1M_CRME_antiRpt_mt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  //aloaderMM.preparecard("MEMu1M_"+era+"_datacard.root","mergedEMu1M");
  aloaderMM.compare(p1,1,p2); // 5
  SaveAs(canvas_2,"MEMuFF_1M_CRME_mt.pdf",p1);
  //aloaderMM.close();

/*  aloaderMM.load("1M_CRME_antiDphi_mt","1M_CRME_antiDphi_antiRpt_mt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_antiDphi_mt","1M_CRME_antiDphi_antiRpt_mt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_antiDphi_mt","1M_CRME_antiDphi_antiRpt_mt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_antiDphi_mt","1M_CRME_antiDphi_antiRpt_mt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1,5);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_mt.png");*/

  /*aloaderMM.load("1M_CRME_CRdphi_mt","1M_CRME_CRdphi_antiRpt_mt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_CRdphi_mt","1M_CRME_CRdphi_antiRpt_mt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_CRdphi_mt","1M_CRME_CRdphi_antiRpt_mt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_CRdphi_mt","1M_CRME_CRdphi_antiRpt_mt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(p1,5,p2);
  SaveAs(canvas_2,"MEMuFF_1M_CRME_CRdphi_mt.pdf",p1);*/

/*  aloaderMM.load("1M_CRME_antiDphi_MET_pt","1M_CRME_antiDphi_antiRpt_MET_pt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_antiDphi_MET_pt","1M_CRME_antiDphi_antiRpt_MET_pt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_antiDphi_MET_pt","1M_CRME_antiDphi_antiRpt_MET_pt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_antiDphi_MET_pt","1M_CRME_antiDphi_antiRpt_MET_pt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_MET_pt.png");

  aloaderMM.load("1M_CRME_CRdphi_MET_pt","1M_CRME_CRdphi_antiRpt_MET_pt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_CRdphi_MET_pt","1M_CRME_CRdphi_antiRpt_MET_pt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_CRdphi_MET_pt","1M_CRME_CRdphi_antiRpt_MET_pt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_CRdphi_MET_pt","1M_CRME_CRdphi_antiRpt_MET_pt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_MET_pt.png");

  aloaderMM.load("1M_CRME_antiDphi_MM_pt","1M_CRME_antiDphi_antiRpt_MM_pt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_antiDphi_MM_pt","1M_CRME_antiDphi_antiRpt_MM_pt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_antiDphi_MM_pt","1M_CRME_antiDphi_antiRpt_MM_pt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_antiDphi_MM_pt","1M_CRME_antiDphi_antiRpt_MM_pt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_MM_pt.png");

  aloaderMM.load("1M_CRME_CRdphi_MM_pt","1M_CRME_CRdphi_antiRpt_MM_pt_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_CRdphi_MM_pt","1M_CRME_CRdphi_antiRpt_MM_pt_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_CRdphi_MM_pt","1M_CRME_CRdphi_antiRpt_MM_pt_xFF","20UL17");
  aloaderMM3.load("1M_CRME_CRdphi_MM_pt","1M_CRME_CRdphi_antiRpt_MM_pt_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1,4);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_MM_pt.png");

  aloaderMM.load("1M_CRME_antiDphi_all_Et","1M_CRME_antiDphi_antiRpt_all_Et_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_antiDphi_all_Et","1M_CRME_antiDphi_antiRpt_all_Et_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_antiDphi_all_Et","1M_CRME_antiDphi_antiRpt_all_Et_xFF","20UL17");
  aloaderMM3.load("1M_CRME_antiDphi_all_Et","1M_CRME_antiDphi_antiRpt_all_Et_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1);
  SaveAs(canvas_1,"FF_1M_CRME_antiDphi_all_Et.png");

  aloaderMM.load("1M_CRME_CRdphi_all_Et","1M_CRME_CRdphi_antiRpt_all_Et_xFF",postfix.Data());
  aloaderMM1.load("1M_CRME_CRdphi_all_Et","1M_CRME_CRdphi_antiRpt_all_Et_xFF","20UL16APV");
  aloaderMM2.load("1M_CRME_CRdphi_all_Et","1M_CRME_CRdphi_antiRpt_all_Et_xFF","20UL17");
  aloaderMM3.load("1M_CRME_CRdphi_all_Et","1M_CRME_CRdphi_antiRpt_all_Et_xFF","20UL18");
  aloaderMM.add(aloaderMM1);
  aloaderMM.add(aloaderMM2);
  aloaderMM.add(aloaderMM3);
  aloaderMM.compare(canvas_1);
  SaveAs(canvas_1,"FF_1M_CRME_CRdphi_all_Et.png");*/

  class HistLoader2M : public HistLoaderBase {
  public:
    HistLoader2M(TFile* adatafile, TFile* aWZfile, TFile* aZZfile, std::vector<SigSample> sigFiles)
    : HistLoaderBase(adatafile,aWZfile,aZZfile),
      sigFiles_(sigFiles) {}

    ~HistLoader2M() {
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
      void setSSFFHist(TH1D* h) { SSFFHist_ = h; }
      void setSSWZHist(TH1D* h) { SSWZHist_ = h; }
      void setSSZZHist(TH1D* h) { SSZZHist_ = h; }
      void setOSWZHist(TH1D* h) { OSWZHist_ = h; }
      void setOSZZHist(TH1D* h) { OSZZHist_ = h; }
      TH1D* SSFFHist() { return SSFFHist_; }
      TH1D* SSWZHist() { return SSWZHist_; }
      TH1D* SSZZHist() { return SSZZHist_; }
      TH1D* OSWZHist() { return OSWZHist_; }
      TH1D* OSZZHist() { return OSZZHist_; }

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

        SSWZHist_ = (TH1D*)WZfile->Get("mergedEMuCRanalyzer"+anlyzrEra+"/"+denomNameSS)->Clone();
        SSWZHist_->Scale( lumi*1000.*WZxsec_/WZsumwgt );
        SSZZHist_ = (TH1D*)ZZfile->Get("mergedEMuCRanalyzer"+anlyzrEra+"/"+denomNameSS)->Clone();
        SSZZHist_->Scale( lumi*1000.*ZZxsec_/ZZsumwgt );
        OSWZHist_ = (TH1D*)WZfile->Get("mergedEMuCRanalyzer"+anlyzrEra+"/"+denomNameOS)->Clone();
        OSWZHist_->Scale( lumi*1000.*WZxsec_/WZsumwgt );
        OSZZHist_ = (TH1D*)ZZfile->Get("mergedEMuCRanalyzer"+anlyzrEra+"/"+denomNameOS)->Clone();
        OSZZHist_->Scale( lumi*1000.*ZZxsec_/ZZsumwgt );
      }

      void rebin(int r) {
        SSFFHist_->Rebin(r);
        SSWZHist_->Rebin(r);
        SSZZHist_->Rebin(r);
        OSWZHist_->Rebin(r);
        OSZZHist_->Rebin(r);
      }

      void rebin(int r, TString name, double* bin) {
        SSFFHist_ = (TH1D*)SSFFHist_->Rebin(r,name+"_SSFF",bin);
        SSWZHist_ = (TH1D*)SSWZHist_->Rebin(r,name+"_SSWZ",bin);
        SSZZHist_ = (TH1D*)SSZZHist_->Rebin(r,name+"_SSZZ",bin);
        OSWZHist_ = (TH1D*)OSWZHist_->Rebin(r,name+"_OSWZ",bin);
        OSZZHist_ = (TH1D*)OSZZHist_->Rebin(r,name+"_OSZZ",bin);
      }

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
        denomSSfinal->SetFillColor(TColor::GetColor("#9c9ca1"));

        return std::move(denomSSfinal);
      }

      std::unique_ptr<TH1D> returnOS() {
        auto denomOSfinal = std::unique_ptr<TH1D>((TH1D*)OSWZHist_->Clone());
        denomOSfinal->Add(OSZZHist_);
        denomOSfinal->SetFillColor(TColor::GetColor("#f89c20"));

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
    denomHists altMuScaleUp_;
    denomHists altMuScaleDn_;
    denomHists preCorrUp_;
    denomHists preCorrDn_;
    denomHists elRecoUp_;
    denomHists elRecoDn_;

    std::vector<SigSample> sigFiles_;
    std::vector<TH1D*> sigHist_;
    std::vector<std::map<std::string,SystVariation>> sigSyst_;

  public:
    void add(const HistLoader2M& other) {
      this->dataHist_->Add(other.dataHist_);
      this->nominal_.add(other.nominal_);
      this->SSFFup_.add(other.SSFFup_);
      this->SSFFdn_.add(other.SSFFdn_);
      this->OSFFup_.add(other.OSFFup_);
      this->OSFFdn_.add(other.OSFFdn_);

      if (!sigSyst_.empty()) {
        this->heepIdUp_.add(other.heepIdUp_);
        this->heepIdDn_.add(other.heepIdDn_);
        this->altMuScaleUp_.add(other.altMuScaleUp_);
        this->altMuScaleDn_.add(other.altMuScaleDn_);
        this->preCorrUp_.add(other.preCorrUp_);
        this->preCorrDn_.add(other.preCorrDn_);
        this->elRecoUp_.add(other.elRecoUp_);
        this->elRecoDn_.add(other.elRecoDn_);
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
        const double sigLumi = 0.01;

        auto retrieveSigHist = [this,&numName,&anlyzrEra,&lumi,&sigLumi] (TFile* afile, const TString& systName) -> TH1D* {
          TH1D* ahist = (TH1D*)afile->Get( "mergedEMuCRanalyzer"+anlyzrEra+"/"+numName+systName )->Clone();
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
          TH1D* sigMuScaleUp = retrieveSigHist(sigFiles_.at(idx).file_,"_altMuScale");
          TH1D* sigMuScaleDn = variateDn(asigHist,sigMuScaleUp);
          init["sigMuScale"] = SystVariation(sigMuScaleUp,sigMuScaleDn);
          TH1D* sigMuSmearUp = retrieveSigHist(sigFiles_.at(idx).file_,"_altMuSmear");
          TH1D* sigMuSmearDn = variateDn(asigHist,sigMuSmearUp);
          init["sigMuSmear"] = SystVariation(sigMuSmearUp,sigMuSmearDn);
          TH1D* sigMuIdUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muIdUp");
          TH1D* sigMuIdDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muIdDn");
          init["sigMuId"] = SystVariation(sigMuIdUp,sigMuIdDn);
          TH1D* sigMuIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muIsoUp");
          TH1D* sigMuIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muIsoDn");
          init["sigMuIso"] = SystVariation(sigMuIsoUp,sigMuIsoDn);
          TH1D* sigMuBoostIsoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoUp");
          TH1D* sigMuBoostIsoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muBoostIsoDn");
          init["sigMuBoostIso"] = SystVariation(sigMuBoostIsoUp,sigMuBoostIsoDn);
          TH1D* sigTrigUp = retrieveSigHist(sigFiles_.at(idx).file_,"_trigUp");
          TH1D* sigTrigDn = retrieveSigHist(sigFiles_.at(idx).file_,"_trigDn");
          init["sigTrig"] = SystVariation(sigTrigUp,sigTrigDn);
          TH1D* sigMuRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_muRecoUp");
          TH1D* sigMuRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_muRecoDn");
          init["sigMuReco"] = SystVariation(sigMuRecoUp,sigMuRecoDn);
          TH1D* sigElRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoUp");
          TH1D* sigElRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoDn");
          init["sigElReco"] = SystVariation(sigElRecoUp,sigElRecoDn);
          TH1D* sigPUrwgtUp = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtUp");
          TH1D* sigPUrwgtDn = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtDn");
          init["sigPUrwgt"] = SystVariation(sigPUrwgtUp,sigPUrwgtDn);
          TH1D* sigPrefireUp = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireUp");
          TH1D* sigPrefireDn = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireDn");
          init["sigPrefire"] = SystVariation(sigPrefireUp,sigPrefireDn);

          sigSyst_.push_back(init);
        }

        heepIdUp_.load(denomName+"_xSSFF_heepIdUp",denomName+"_xOSFF_heepIdUp",datafile_,WZfile_,ZZfile_,anlyzrEra);
        heepIdDn_.load(denomName+"_xSSFF_heepIdDn",denomName+"_xOSFF_heepIdDn",datafile_,WZfile_,ZZfile_,anlyzrEra);
        elRecoUp_.load(denomName+"_xSSFF_elRecoUp",denomName+"_xOSFF_elRecoUp",datafile_,WZfile_,ZZfile_,anlyzrEra);
        elRecoDn_.load(denomName+"_xSSFF_elRecoDn",denomName+"_xOSFF_elRecoDn",datafile_,WZfile_,ZZfile_,anlyzrEra);
        altMuScaleUp_.load(denomName+"_altMuScale_xSSFF",denomName+"_altMuScale_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
        altMuScaleDn_.setSSFFHist( variateDn(nominal_.SSFFHist(),altMuScaleUp_.SSFFHist()) );
        altMuScaleDn_.setSSWZHist( variateDn(nominal_.SSWZHist(),altMuScaleUp_.SSWZHist()) );
        altMuScaleDn_.setSSZZHist( variateDn(nominal_.SSZZHist(),altMuScaleUp_.SSZZHist()) );
        altMuScaleDn_.setOSWZHist( variateDn(nominal_.OSWZHist(),altMuScaleUp_.OSWZHist()) );
        altMuScaleDn_.setOSZZHist( variateDn(nominal_.OSZZHist(),altMuScaleUp_.OSZZHist()) );
        preCorrUp_.load(denomName+"_xSSFF_preCorr",denomName+"_xOSFF_preCorr",datafile_,WZfile_,ZZfile_,anlyzrEra);
        preCorrDn_.setSSFFHist( variateDn(nominal_.SSFFHist(),preCorrUp_.SSFFHist()) );
        preCorrDn_.setSSWZHist( variateDn(nominal_.SSWZHist(),preCorrUp_.SSWZHist()) );
        preCorrDn_.setSSZZHist( variateDn(nominal_.SSZZHist(),preCorrUp_.SSZZHist()) );
        preCorrDn_.setOSWZHist( variateDn(nominal_.OSWZHist(),preCorrUp_.OSWZHist()) );
        preCorrDn_.setOSZZHist( variateDn(nominal_.OSZZHist(),preCorrUp_.OSZZHist()) );
      }
    }

    void compare(TPad* padUp, int rebin=1, TPad* padDn=nullptr) {
      if (rebin!=1) {
        dataHist_->Rebin(rebin);

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx)->Rebin(rebin);

          for (const auto& element : sigSyst_.at(idx)) {
            element.second.up_->Rebin( rebin );
            element.second.dn_->Rebin( rebin );
          }
        }
      }

      nominal_.rebin(nominal_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      SSFFup_.rebin(SSFFup_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      SSFFdn_.rebin(SSFFdn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      OSFFup_.rebin(OSFFup_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      OSFFdn_.rebin(OSFFdn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());

      if ( TString(dataHist_->GetName()).Contains("invM") ) {
        heepIdUp_.rebin(heepIdUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        heepIdDn_.rebin(heepIdDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        altMuScaleUp_.rebin(altMuScaleUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        altMuScaleDn_.rebin(altMuScaleDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        elRecoUp_.rebin(elRecoUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        elRecoDn_.rebin(elRecoDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        preCorrUp_.rebin(preCorrUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        preCorrDn_.rebin(preCorrDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
      }

      if (rebin==1 && true) {
        std::vector<double> binEdges = {0.,50.,60.,70.,80.,90.,100.,110.,120.,130.,140.,150.,160.,170.,180.,190.,200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,550.,600.,700.,800.,1000.,2500.};

        dataHist_ = (TH1D*)dataHist_->Rebin(binEdges.size()-1,TString(dataHist_->GetName())+"_rebin",&(binEdges[0]));

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          sigHist_.at(idx) = (TH1D*)sigHist_.at(idx)->Rebin( binEdges.size()-1,TString(sigHist_.at(idx)->GetName())+"_rebin",&(binEdges[0]) );

          for (const auto& element : sigSyst_.at(idx)) {
            sigSyst_.at(idx).at(element.first).up_ = (TH1D*)sigSyst_.at(idx).at(element.first).up_->Rebin( binEdges.size()-1,TString(sigSyst_.at(idx).at(element.first).up_->GetName())+"_rebin",&(binEdges[0]) );
            sigSyst_.at(idx).at(element.first).dn_ = (TH1D*)sigSyst_.at(idx).at(element.first).dn_->Rebin( binEdges.size()-1,TString(sigSyst_.at(idx).at(element.first).dn_->GetName())+"_rebin",&(binEdges[0]) );
          }
        }

        nominal_.rebin( binEdges.size()-1,TString("nominal_rebin"),&(binEdges[0]) );
        SSFFup_.rebin( binEdges.size()-1,TString("SSFFup_rebin"),&(binEdges[0]) );
        SSFFdn_.rebin( binEdges.size()-1,TString("SSFFdn_rebin"),&(binEdges[0]) );
        OSFFup_.rebin( binEdges.size()-1,TString("OSFFup_rebin"),&(binEdges[0]) );
        OSFFdn_.rebin( binEdges.size()-1,TString("OSFFdn_rebin"),&(binEdges[0]) );

        heepIdUp_.rebin( binEdges.size()-1,TString("heepIdUp_rebin"),&(binEdges[0]) );
        heepIdDn_.rebin( binEdges.size()-1,TString("heepIdDn_rebin"),&(binEdges[0]) );
        altMuScaleUp_.rebin( binEdges.size()-1,TString("altMuScaleUp_rebin"),&(binEdges[0]) );
        altMuScaleDn_.rebin( binEdges.size()-1,TString("altMuScaleDn_rebin"),&(binEdges[0]) );
        elRecoUp_.rebin( binEdges.size()-1,TString("elRecoUp_rebin"),&(binEdges[0]) );
        elRecoDn_.rebin( binEdges.size()-1,TString("elRecoDn_rebin"),&(binEdges[0]) );
        preCorrUp_.rebin( binEdges.size()-1,TString("preCorrUp_rebin"),&(binEdges[0]) );
        preCorrDn_.rebin( binEdges.size()-1,TString("preCorrDn_rebin"),&(binEdges[0]) );
      }

      auto denomSSfinal = nominal_.returnSS();
      auto denomOSfinal = nominal_.returnOS();

      auto denomfinal_nominal = nominal_.returnAdded();
      auto denomfinal_SSFFup = SSFFup_.returnAdded();
      auto denomfinal_SSFFdn = SSFFdn_.returnAdded();
      auto denomfinal_OSFFup = OSFFup_.returnAdded();
      auto denomfinal_OSFFdn = OSFFdn_.returnAdded();
      std::unique_ptr<TH1D> denomfinal_heepIdUp, denomfinal_heepIdDn;
      std::unique_ptr<TH1D> denomfinal_muScaleUp, denomfinal_muScaleDn;
      std::unique_ptr<TH1D> denomfinal_precorrUp, denomfinal_precorrDn;
      std::unique_ptr<TH1D> denomfinal_elRecoUp, denomfinal_elRecoDn;

      if ( TString(dataHist_->GetName()).Contains("invM") ) {
        denomfinal_heepIdUp = heepIdUp_.returnAdded();
        denomfinal_heepIdDn = heepIdDn_.returnAdded();
        denomfinal_muScaleUp = altMuScaleUp_.returnAdded();
        denomfinal_muScaleDn = altMuScaleDn_.returnAdded();
        denomfinal_elRecoUp = elRecoUp_.returnAdded();
        denomfinal_elRecoDn = elRecoDn_.returnAdded();
        denomfinal_precorrUp = preCorrUp_.returnAdded();
        denomfinal_precorrDn = preCorrDn_.returnAdded();
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
        dataHist_->GetXaxis()->SetTitle("M_{e#mu#mu} [GeV]");
        dataHist_->SetMaximum(3.*dataHist_->GetMaximum());
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      } else if (TString(dataHist_->GetName()).Contains("Et")) {
        dataHist_->GetXaxis()->SetTitle("E_{T}^{5#times5} [GeV]");
        dataHist_->GetXaxis()->SetRangeUser(0.,250.);
        dataHist_->SetMaximum(1.3*dataHist_->GetMaximum());
      } else if (TString(dataHist_->GetName()).Contains("eta")) {
        dataHist_->GetXaxis()->SetTitle("#eta^{5#times5}");
        dataHist_->SetMaximum(1.5*dataHist_->GetMaximum());
      }

      padUp->cd();
      dataHist_->SetMinimum(0.001);
      dataHist_->Draw("E1");
      denomFinal->Draw("hist&same");
      dataHist_->Draw("E1&same");

      if (TString(dataHist_->GetName()).Contains("invM")) {
        for (auto* sigNum : sigHist_)
          sigNum->Draw("hist&same");
      }

      if (true) {
        TLegend* legend = nullptr;

        if (TString(dataHist_->GetName()).Contains("invM"))
          legend = new TLegend(0.55,0.58,0.95,0.93);
        else
          legend = new TLegend(0.55,0.65,0.95,0.93);

        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->AddEntry(dataHist_,"Data");
        legend->AddEntry(denomSSfinal.get(),"Nonprompt bkg");
        legend->AddEntry(denomOSfinal.get(),"Prompt bkg");

        if (TString(dataHist_->GetName()).Contains("invM")) {
          legend->AddEntry(sigHist_.at(0),"X750Y1");
          //legend_left->AddEntry(sigNums.at(isigDiv),"M_{A} = 10 GeV (#sigma = 10 fb)");
        }

        legend->Draw();
      }

      if (padDn) {
        TH1D* ratio = (TH1D*)dataHist_->Clone();
        ratio->SetStats(0);
        ratio->SetTitle("");
        ratio->Divide(denomfinal_nominal.get());
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
      }

      std::vector<double> x0, y0, errx, erryDn, erryUp;
      std::vector<double> r0, errRdn, errRup;

      for (unsigned idx = 1; idx <= denomfinal_nominal->GetNbinsX(); idx++) {
        x0.push_back(denomfinal_nominal->GetBinCenter(idx));
        y0.push_back(denomfinal_nominal->GetBinContent(idx));
        errx.push_back(denomfinal_nominal->GetBinWidth(idx)/2.);

        double valSSFFup = denomfinal_SSFFup->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
        double valSSFFdn = denomfinal_nominal->GetBinContent(idx) - denomfinal_SSFFdn->GetBinContent(idx);
        double valOSFFup = denomfinal_OSFFup->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
        double valOSFFdn = denomfinal_nominal->GetBinContent(idx) - denomfinal_OSFFdn->GetBinContent(idx);
        double normSSFF = 0.;//0.2*nominal_.returnSS()->GetBinContent(idx);
        double normOSFF = 0.2*nominal_.returnOS()->GetBinContent(idx);

        double valHeepIdUp = 0., valHeepIdDn = 0., valMuScaleUp = 0., valMuScaleDn = 0.;

        double valPreCorrUp = 0., valPreCorrDn = 0.;
        double valElRecoUp = 0., valElRecoDn = 0.;

        if ( TString(dataHist_->GetName()).Contains("invM") ) {
          valHeepIdUp = denomfinal_heepIdUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valHeepIdDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_heepIdDn->GetBinContent(idx);
          valMuScaleUp = denomfinal_muScaleUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valMuScaleDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_muScaleDn->GetBinContent(idx);
          valPreCorrUp = denomfinal_precorrUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valPreCorrDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_precorrDn->GetBinContent(idx);
          valElRecoUp = denomfinal_elRecoUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valElRecoDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_elRecoDn->GetBinContent(idx);
        }

        erryUp.push_back( std::sqrt( valSSFFup*valSSFFup + valOSFFup*valOSFFup + valHeepIdUp*valHeepIdUp + valMuScaleUp*valMuScaleUp
                                   + valPreCorrUp*valPreCorrUp + valElRecoUp*valElRecoUp + normSSFF*normSSFF + normOSFF*normOSFF ) );
        erryDn.push_back( std::sqrt( valSSFFdn*valSSFFdn + valOSFFdn*valOSFFdn + valHeepIdDn*valHeepIdDn + valMuScaleDn*valMuScaleDn
                                   + valPreCorrDn*valPreCorrDn + valElRecoDn*valElRecoDn + normSSFF*normSSFF + normOSFF*normOSFF ) );

        if (padDn) {
          r0.push_back(1.);

          double rvalSSFFup = denomfinal_SSFFup->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
          double rvalSSFFdn = 1. - denomfinal_SSFFdn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
          double rvalOSFFup = denomfinal_OSFFup->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
          double rvalOSFFdn = 1. - denomfinal_OSFFdn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
          double rnormSSFF = normSSFF/denomfinal_nominal->GetBinContent(idx);
          double rnormOSFF = normOSFF/denomfinal_nominal->GetBinContent(idx);

          double rHeepIdUp = 0., rHeepIdDn = 0., rPreCorrUp = 0., rPreCorrDn = 0.;
          double rElRecoUp = 0., rElRecoDn = 0., rMuScaleUp = 0., rMuScaleDn = 0.;

          if ( TString(dataHist_->GetName()).Contains("invM") ) {
            rHeepIdUp = denomfinal_heepIdUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rHeepIdDn = 1. - denomfinal_heepIdDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
            rPreCorrUp = denomfinal_precorrUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rPreCorrDn = 1. - denomfinal_precorrDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
            rElRecoUp = denomfinal_elRecoUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rElRecoDn = 1. - denomfinal_elRecoDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
            rMuScaleUp = denomfinal_muScaleUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rMuScaleDn = 1. - denomfinal_muScaleDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
          }

          double rUp = std::sqrt( rvalSSFFup*rvalSSFFup + rvalOSFFup*rvalOSFFup + rHeepIdUp*rHeepIdUp
                                       + rPreCorrUp*rPreCorrUp + rElRecoUp*rElRecoUp + rMuScaleUp*rMuScaleUp + rnormSSFF*rnormSSFF + rnormOSFF*rnormOSFF );
          double rDn = std::sqrt( rvalSSFFdn*rvalSSFFdn + rvalOSFFdn*rvalOSFFdn + rHeepIdDn*rHeepIdDn
                                       + rPreCorrDn*rPreCorrDn + rElRecoDn*rElRecoDn + rMuScaleDn*rMuScaleDn + rnormSSFF*rnormSSFF + rnormOSFF*rnormOSFF );

          errRup.push_back( denomfinal_nominal->GetBinContent(idx) > 0. ? rUp : 0. );
          errRdn.push_back( denomfinal_nominal->GetBinContent(idx) > 0. ? rDn : 0. );
        }
      }

      padUp->cd();

      auto gr = new TGraphAsymmErrors(denomfinal_nominal->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+2);
      gr->SetLineColor(kGray+2);
      gr->SetFillStyle(3004);
      gr->Draw("2");

      if (padDn) {
        auto rgr = new TGraphAsymmErrors(dataHist_->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
        rgr->SetFillColor(kGray+2);
        rgr->SetLineColor(kGray+2);
        rgr->SetFillStyle(3004);
        padDn->cd();
        rgr->Draw("2");
      }

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
        dir_->WriteTObject(altMuScaleUp_.returnSS().release(),"SS_muMomentumScaleUp");
        dir_->WriteTObject(altMuScaleUp_.returnOS().release(),"OS_muMomentumScaleUp");
        dir_->WriteTObject(altMuScaleDn_.returnSS().release(),"SS_muMomentumScaleDown");
        dir_->WriteTObject(altMuScaleDn_.returnOS().release(),"OS_muMomentumScaleDown");
        dir_->WriteTObject(preCorrUp_.returnSS().release(),"SS_mergedEleEnCorrUp");
        dir_->WriteTObject(preCorrUp_.returnOS().release(),"OS_mergedEleEnCorrUp");
        dir_->WriteTObject(preCorrDn_.returnSS().release(),"SS_mergedEleEnCorrDown");
        dir_->WriteTObject(preCorrDn_.returnOS().release(),"OS_mergedEleEnCorrDown");
        dir_->WriteTObject(elRecoUp_.returnSS().release(),"SS_elRecoUp");
        dir_->WriteTObject(elRecoUp_.returnOS().release(),"OS_elRecoUp");
        dir_->WriteTObject(elRecoDn_.returnSS().release(),"SS_elRecoDown");
        dir_->WriteTObject(elRecoDn_.returnOS().release(),"OS_elRecoDown");

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

          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuScale").up_,sigFiles_.at(idx).name_+"_muMomentumScaleUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuScale").dn_,sigFiles_.at(idx).name_+"_muMomentumScaleDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuSmear").up_,sigFiles_.at(idx).name_+"_muMomentumSmearUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuSmear").dn_,sigFiles_.at(idx).name_+"_muMomentumSmearDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").up_,sigFiles_.at(idx).name_+"_highPtIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuId").dn_,sigFiles_.at(idx).name_+"_highPtIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").up_,sigFiles_.at(idx).name_+"_muLooseIsoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuIso").dn_,sigFiles_.at(idx).name_+"_muLooseIsoDown");
          // dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").up_,sigFiles_.at(idx).name_+"_muBoostIsoUp");
          // dir_->WriteTObject(sigSyst_.at(idx).at("sigMuBoostIso").dn_,sigFiles_.at(idx).name_+"_muBoostIsoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").up_,sigFiles_.at(idx).name_+"_muTrigUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigTrig").dn_,sigFiles_.at(idx).name_+"_muTrigDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").up_,sigFiles_.at(idx).name_+"_muRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigMuReco").dn_,sigFiles_.at(idx).name_+"_muRecoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").up_,sigFiles_.at(idx).name_+"_elRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").dn_,sigFiles_.at(idx).name_+"_elRecoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").up_,sigFiles_.at(idx).name_+"_PUrwgtUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").dn_,sigFiles_.at(idx).name_+"_PUrwgtDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").up_,sigFiles_.at(idx).name_+"_prefireUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").dn_,sigFiles_.at(idx).name_+"_prefireDown");
        }
      }

      denomSSfinal.release();
      denomOSfinal.release();
    }
  };

  canvas_1->cd();

  auto aloader3E = HistLoader2M(datafile,WZfile,ZZfile,sigsamples);
  auto aloader3E1 = HistLoader2M(datafile1,WZfile1,ZZfile1,sigsamples1);
  auto aloader3E2 = HistLoader2M(datafile2,WZfile2,ZZfile2,sigsamples2);
  auto aloader3E3 = HistLoader2M(datafile3,WZfile3,ZZfile3,sigsamples3);

  //canvas_1->SetLogy();
  aloader3E.load("2M_CRME_lll_invM","2M_antiME_lll_invM_CR",postfix.Data());
  aloader3E1.load("2M_CRME_lll_invM","2M_antiME_lll_invM_CR","20UL16APV");
  aloader3E2.load("2M_CRME_lll_invM","2M_antiME_lll_invM_CR","20UL17");
  aloader3E3.load("2M_CRME_lll_invM","2M_antiME_lll_invM_CR","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  //aloader3E.preparecard("MEMu2M_"+era+"_datacard.root","mergedEMu2M");
  aloader3E.compare(p1,1,p2); // 5
  SaveAs(canvas_2,"MEMuFF_2M_invM.pdf",p1);
  //aloader3E.close();
  canvas_1->SetLogy(0);

  return;

  aloader3E.load("2M_Et_CR_EB_CRME","2M_Et_CR_EB_antiME",postfix.Data());
  aloader3E1.load("2M_Et_CR_EB_CRME","2M_Et_CR_EB_antiME","20UL16APV");
  aloader3E2.load("2M_Et_CR_EB_CRME","2M_Et_CR_EB_antiME","20UL17");
  aloader3E3.load("2M_Et_CR_EB_CRME","2M_Et_CR_EB_antiME","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(p1,1,p2);
  SaveAs(canvas_2,"FF_2M_Et.pdf",p1);

  aloader3E.load("2M_eta_CR_EB_CRME","2M_antiME_eta",postfix.Data());
  aloader3E1.load("2M_eta_CR_EB_CRME","2M_antiME_eta","20UL16APV");
  aloader3E2.load("2M_eta_CR_EB_CRME","2M_antiME_eta","20UL17");
  aloader3E3.load("2M_eta_CR_EB_CRME","2M_antiME_eta","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(p1,2,p2);
  SaveAs(canvas_2,"FF_2M_eta.pdf",p1);

  /*aloader3E.load("2M_CRME_lll_mll","2M_antiME_lll_mll",postfix.Data());
  aloader3E1.load("2M_CRME_lll_mll","2M_antiME_lll_mll","20UL16APV");
  aloader3E2.load("2M_CRME_lll_mll","2M_antiME_lll_mll","20UL17");
  aloader3E3.load("2M_CRME_lll_mll","2M_antiME_lll_mll","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(p1,2,p2);
  SaveAs(canvas_2,"FF_2M_mll.png",p1);*/

/*  aloader3E.load("2M_CRME_lll_MET","2M_antiME_lll_MET",postfix.Data());
  aloader3E1.load("2M_CRME_lll_MET","2M_antiME_lll_MET","20UL16APV");
  aloader3E2.load("2M_CRME_lll_MET","2M_antiME_lll_MET","20UL17");
  aloader3E3.load("2M_CRME_lll_MET","2M_antiME_lll_MET","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(canvas_1,2);
  SaveAs(canvas_1,"FF_2M_MET.png");

  aloader3E.load("2M_CRME_lll_passConvVeto","2M_antiME_lll_passConvVeto",postfix.Data());
  aloader3E1.load("2M_CRME_lll_passConvVeto","2M_antiME_lll_passConvVeto","20UL16APV");
  aloader3E2.load("2M_CRME_lll_passConvVeto","2M_antiME_lll_passConvVeto","20UL17");
  aloader3E3.load("2M_CRME_lll_passConvVeto","2M_antiME_lll_passConvVeto","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2M_passConvVeto.png");*/

  canvas_2->cd();

  if (era=="run2")
    return;

  return;

  auto estimateCenter = [] (const std::vector<double>& vec) -> std::vector<double> {
    std::vector<double> out;
    out.push_back(100.);

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

  std::vector<double> xbins = {0,100,200,250,300,350,400,500,600,800,1000,1200,1500};
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
  numer->SetMaximum(2.5);
  numer->SetMinimum(0.0);
  numer->GetYaxis()->SetTitle("Fake factor");
  //numer->GetXaxis()->SetTitle("p_{T}");
  numer->SetLineColor(kBlack);
  numer->SetLineWidth(2);
  numer->SetStats(0);

  TF1* ffFunc = new TF1("MMFF","[0]*x+[1]",100,2000);
  ffFunc->SetLineColor(kRed);
  ffFunc->SetLineWidth(2);
  ffFunc->SetLineStyle(2);

  TFitResultPtr fitResult = numer->Fit(ffFunc,"RS");
  fitResult->SetName("fitResult");
  double ci[nbins+1];
  fitResult->GetConfidenceIntervals(nbins+1,1,0,&(xcen[0]),ci,0.95,false); // 0.6827
  std::vector<double> xbinw = estimateWidth(xbins);
  double ybin[nbins+1];

  for (unsigned idx = 0; idx < nbins+1; idx++) {
    ybin[idx] = ffFunc->Eval(xcen[idx]);
  }

  auto errGr = new TGraphErrors(nbins+1,&(xcen[0]),ybin,&(xbinw[0]),ci);
  errGr->SetFillColor(kRed);
  errGr->SetFillStyle(3005);

  canvas_2->cd();
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

  SaveAs(canvas_2,("MMFF_"+era+".pdf").Data());

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
