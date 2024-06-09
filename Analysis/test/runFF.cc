#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"
#include "TRandom3.h"

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

  TFile* datafile = new TFile("EleAnalyzer_"+fname+"_data.root","READ");
  TFile* WZfile = new TFile("EleAnalyzer_"+fname+"_WZFXFX.root","READ");
  TFile* ZZfile = new TFile("EleAnalyzer_"+fname+"_ZZ.root","READ");

  TFile* datafile1 = new TFile("EleAnalyzer_20UL16APV_data.root","READ");
  TFile* WZfile1 = new TFile("EleAnalyzer_20UL16APV_WZFXFX.root","READ");
  TFile* ZZfile1 = new TFile("EleAnalyzer_20UL16APV_ZZ.root","READ");

  TFile* datafile2 = new TFile("EleAnalyzer_20UL17_data.root","READ");
  TFile* WZfile2 = new TFile("EleAnalyzer_20UL17_WZFXFX.root","READ");
  TFile* ZZfile2 = new TFile("EleAnalyzer_20UL17_ZZ.root","READ");

  TFile* datafile3 = new TFile("EleAnalyzer_20UL18_data.root","READ");
  TFile* WZfile3 = new TFile("EleAnalyzer_20UL18_WZFXFX.root","READ");
  TFile* ZZfile3 = new TFile("EleAnalyzer_20UL18_ZZ.root","READ");

  //TFile* H250A1file = new TFile("MergedEleCR_"+era+"_H250A1.root","READ");
  //TFile* H750A1file = new TFile("MergedEleCR_"+era+"_H750A1.root","READ");
  //TFile* H2000A1file = new TFile("MergedEleCR_"+era+"_H2000A1.root","READ");
  //TFile* H250A10file = new TFile("MergedEleCR_"+era+"_H250A10.root","READ");
  //TFile* H750A10file = new TFile("MergedEleCR_"+era+"_H750A10.root","READ");
  //TFile* H2000A10file = new TFile("MergedEleCR_"+era+"_H2000A10.root","READ");

  class SigSample {
  public:
    SigSample(TFile* afile, TString name)
    : file_(afile), name_(name) {}

    ~SigSample()=default;

  public:
    TFile* file_;
    TString name_;
  };

  // auto H750A1sample = SigSample(new TFile("EleAnalyzer_"+era+"_H750A1.root","READ"),"H750A1");

  /*std::vector<SigSample> sigsamples = {
    SigSample(new TFile("EleAnalyzer_"+fname+"_H250A1.root","READ"),"H250A1"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H250A2.root","READ"),"H250A2"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H250A5.root","READ"),"H250A5"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H250A10.root","READ"),"H250A10"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H250A50.root","READ"),"H250A50"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H250A100.root","READ"),"H250A100"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A1.root","READ"),"H750A1"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A2.root","READ"),"H750A2"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A5.root","READ"),"H750A5"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A10.root","READ"),"H750A10"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A50.root","READ"),"H750A50"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A100.root","READ"),"H750A100"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H750A250.root","READ"),"H750A250"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A1.root","READ"),"H2000A1"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A2.root","READ"),"H2000A2"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A5.root","READ"),"H2000A5"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A10.root","READ"),"H2000A10"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A50.root","READ"),"H2000A50"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A100.root","READ"),"H2000A100"),
    SigSample(new TFile("EleAnalyzer_"+fname+"_H2000A750.root","READ"),"H2000A750")
  };

  std::vector<SigSample> sigsamples1 = {
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

  std::vector<SigSample> sigsamples2 = {
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

  std::vector<SigSample> sigsamples3 = {
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
  };*/

  auto sigsamples = std::vector<SigSample>{SigSample(new TFile("EleAnalyzer_"+fname+"_H750A1.root","READ"),"H750A1")};
  auto sigsamples1 = std::vector<SigSample>{SigSample(new TFile("EleAnalyzer_20UL16APV_H750A1.root","READ"),"H750A1")};
  auto sigsamples2 = std::vector<SigSample>{SigSample(new TFile("EleAnalyzer_20UL17_H750A1.root","READ"),"H750A1")};
  auto sigsamples3 = std::vector<SigSample>{SigSample(new TFile("EleAnalyzer_20UL18_H750A1.root","READ"),"H750A1")};

  std::vector<SigSample> sigsamples3E = {SigSample(new TFile("EleAnalyzer_"+fname+"_H750A10.root","READ"),"H750A10")};
  std::vector<SigSample> sigsamples3E1 = {SigSample(new TFile("EleAnalyzer_20UL16APV_H750A10.root","READ"),"H750A10")};
  std::vector<SigSample> sigsamples3E2 = {SigSample(new TFile("EleAnalyzer_20UL17_H750A10.root","READ"),"H750A10")};
  std::vector<SigSample> sigsamples3E3 = {SigSample(new TFile("EleAnalyzer_20UL18_H750A10.root","READ"),"H750A10")};

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

    void load(const std::string& nameNum, const std::string& name, int color, std::string anlyzrEra="", const std::string& denomtree="", bool isOS=true, const std::string& nameNum2="", const std::string& name2="", int color2=0) {
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

      auto rng = TRandom3(0);
      const double scale = 0.9964;
      const double scale2 = 1.032;
      const double smear = 0.0158;

      if (denomtree!="") {
        FFHist_ = new TH1D(("FFhist_"+denomtree).c_str(),";GeV;",FFHist_->GetNbinsX(),FFHist_->GetXaxis()->GetXmin(),FFHist_->GetXaxis()->GetXmax());
        TTree* atree = (TTree*)datafile_->Get( (std::string(anlyzrData+"/")+denomtree).c_str() );

        float pt1, eta1, phi1, pt2, eta2, phi2, wgt;
        int os;
        atree->SetBranchAddress("pt1",&pt1);
        atree->SetBranchAddress("eta1",&eta1);
        atree->SetBranchAddress("phi1",&phi1);
        atree->SetBranchAddress("pt2",&pt2);
        atree->SetBranchAddress("eta2",&eta2);
        atree->SetBranchAddress("phi2",&phi2);
        atree->SetBranchAddress("wgt",&wgt);
        atree->SetBranchAddress("isOS",&os);

        for (unsigned idx=0; idx<atree->GetEntries(); idx++) {
          atree->GetEntry(idx);

          if (static_cast<int>(isOS)!=os)
            continue;

          auto lvec1 = ROOT::Math::PtEtaPhiMVector(pt1,eta1,phi1,0.);
          auto lvec2 = ROOT::Math::PtEtaPhiMVector(pt2,eta2,phi2,0.);

          const double ptll = (lvec1+lvec2).Pt();

          if (ptll < 50.)
            lvec2 = ROOT::Math::PtEtaPhiMVector(pt2*scale2 + rng.Gaus(0.,pt2*scale2*smear),eta2,phi2,0.);
          else
            lvec2 = ROOT::Math::PtEtaPhiMVector(pt2*scale + rng.Gaus(0.,pt2*scale*smear),eta2,phi2,0.);

          double finalWgt = wgt;

          if (denomtree=="antiCRtree")
            finalWgt *= 0.5;

          FFHist_->Fill( (lvec1+lvec2).M(), finalWgt );

          if (denomtree=="antiCRtree") {
            if (ptll < 50.)
              lvec1 = ROOT::Math::PtEtaPhiMVector(pt1*scale2 + rng.Gaus(0.,pt1*scale2*smear),eta1,phi1,0.);
            else
              lvec1 = ROOT::Math::PtEtaPhiMVector(pt1*scale + rng.Gaus(0.,pt1*scale*smear),eta1,phi1,0.);

            lvec2 = ROOT::Math::PtEtaPhiMVector(pt2,eta2,phi2,0.);
            FFHist_->Fill( (lvec1+lvec2).M(), finalWgt );
          }
        }
      }

      FFHist_->SetLineWidth(0);
      FFHist_->SetFillColor(color);

      if ( TString(name).Contains("invM") || TString(name).Contains("eta") ) {
        TH1D* FFHistUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_up").c_str() )->Clone();
        TH1D* FFHistDn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_dn").c_str() )->Clone();
        syst_["FFHist"] = SystVariation(FFHistUp,FFHistDn);
      }

      if ( TString(name).Contains("invM") ) {
        TH1D* FFHistPreCorrUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name+"_preCorr").c_str() )->Clone();
        TH1D* FFHistPreCorrDn = variateDn(FFHist_,FFHistPreCorrUp);
        syst_["FFHistAbcdMergedEleEnCorr"] = SystVariation(FFHistPreCorrUp,FFHistPreCorrDn);
      }

      const double lumi = retrieveLumi(anlyzrEra);
      const double sigLumi = 0.01;
      sigHist_.clear();
      sigSyst_.clear();

      if (nameNum2!="") {
        dataHist2_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+nameNum2).c_str() )->Clone();
        FFHist2_ = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2).c_str() )->Clone();

        if (denomtree!="") {
          FFHist2_ = new TH1D(("FFhist2_"+denomtree).c_str(),";GeV;",FFHist2_->GetNbinsX(),FFHist2_->GetXaxis()->GetXmin(),FFHist2_->GetXaxis()->GetXmax());
          TTree* atree = (TTree*)datafile_->Get( (std::string(anlyzrData+"/")+denomtree).c_str() );

          float pt1, eta1, phi1, pt2, eta2, phi2, wgt;
          int os;
          atree->SetBranchAddress("pt1",&pt1);
          atree->SetBranchAddress("eta1",&eta1);
          atree->SetBranchAddress("phi1",&phi1);
          atree->SetBranchAddress("pt2",&pt2);
          atree->SetBranchAddress("eta2",&eta2);
          atree->SetBranchAddress("phi2",&phi2);
          atree->SetBranchAddress("wgt",&wgt);
          atree->SetBranchAddress("isOS",&os);

          for (unsigned idx=0; idx<atree->GetEntries(); idx++) {
            atree->GetEntry(idx);

            if (static_cast<int>(isOS)==os)
              continue;

            auto lvec1 = ROOT::Math::PtEtaPhiMVector(pt1,eta1,phi1,0.);
            auto lvec2 = ROOT::Math::PtEtaPhiMVector(pt2,eta2,phi2,0.);

            const double ptll = (lvec1+lvec2).Pt();

            if (ptll < 50.)
              lvec2 = ROOT::Math::PtEtaPhiMVector(pt2*scale2 + rng.Gaus(0.,pt2*scale2*smear),eta2,phi2,0.);
            else
              lvec2 = ROOT::Math::PtEtaPhiMVector(pt2*scale + rng.Gaus(0.,pt2*scale*smear),eta2,phi2,0.);

            double finalWgt = wgt;

            if (denomtree=="antiCRtree")
              finalWgt *= 0.5;

            FFHist2_->Fill( (lvec1+lvec2).M(), finalWgt );

            if (denomtree=="antiCRtree") {
              if (ptll < 50.)
                lvec1 = ROOT::Math::PtEtaPhiMVector(pt1*scale2 + rng.Gaus(0.,pt1*scale2*smear),eta1,phi1,0.);
              else
                lvec1 = ROOT::Math::PtEtaPhiMVector(pt1*scale + rng.Gaus(0.,pt1*scale*smear),eta1,phi1,0.);

              lvec2 = ROOT::Math::PtEtaPhiMVector(pt2,eta2,phi2,0.);
              FFHist2_->Fill( (lvec1+lvec2).M(), finalWgt );
            }
          }
        }

        FFHist2_->SetLineWidth(0);
        FFHist2_->SetFillColor(color2);

        if ( TString(name2).Contains("invM") ) {
          TH1D* FFHist2Up = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_up").c_str() )->Clone();
          TH1D* FFHist2Dn = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_dn").c_str() )->Clone();
          syst_["FFHist2"] = SystVariation(FFHist2Up,FFHist2Dn);

          TH1D* FFHist2preCorrUp = (TH1D*)datafile_->Get( (std::string(anlyzrData+"/")+name2+"_preCorr").c_str() )->Clone();
          TH1D* FFHist2preCorrDn = variateDn(FFHist2_,FFHist2preCorrUp);
          syst_["FFHist2abcdMergedEleEnCorr"] = SystVariation(FFHist2preCorrUp,FFHist2preCorrDn);
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
              TH1D* sigTrigUp = retrieveSigHist(sigFiles_.at(idx).file_,"_elTrigUp");
              TH1D* sigTrigDn = retrieveSigHist(sigFiles_.at(idx).file_,"_elTrigDn");
              init["sigElTrig"] = SystVariation(sigTrigUp,sigTrigDn);
              TH1D* sigRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoUp");
              TH1D* sigRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoDn");
              init["sigElReco"] = SystVariation(sigRecoUp,sigRecoDn);
              TH1D* sigPUrwgtUp = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtUp");
              TH1D* sigPUrwgtDn = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtDn");
              init["sigPUrwgt"] = SystVariation(sigPUrwgtUp,sigPUrwgtDn);
              TH1D* sigPrefireUp = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireUp");
              TH1D* sigPrefireDn = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireDn");
              init["sigPrefire"] = SystVariation(sigPrefireUp,sigPrefireDn);
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
      FFHist_->SetMarkerSize(0);

      if (dataHist2_) {
        dataHist2_->Rebin(rebin);
        FFHist2_->Rebin(rebin);
        FFHist2_->SetMarkerSize(0);
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
        FFHist2_->Add(FFadded);
        FFadded->Add(FFHist2_);
        auto* stack = new THStack("stack",";GeV;");
        stack->Add(FFadded);
        // stack->Add(FFHist2_);
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

      if ( TString(dataHist_->GetName()).Contains("invM") )
        dataHist->GetXaxis()->SetTitle("M_{ee} [GeV]");

      dataHist->SetLineWidth(2);
      dataHist->SetLineColor(kBlack);
      dataHist->Draw("E1");
      FFobj->Draw("hist&same");
      dataHist->Draw("E1&same");

      if (dataHist2_) {
        TLegend* legend = nullptr;

        if (TString(dataHist_->GetName()).Contains("CRME") && TString(dataHist_->GetName()).Contains("invM"))
          legend = new TLegend(0.58,0.65,0.95,0.93);
        else
          legend = new TLegend(0.58,0.7,0.95,0.93);

        legend->SetBorderSize(0);
        legend->AddEntry(dataHist,"Data");

//        legend->AddEntry(FFHist_,"Nonprompt bkg");
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

        TH1D* FFHistPreCorrUpAdded = nullptr;
        TH1D* FFHistPreCorrDnAdded = nullptr;

        if (syst_.count("FFHist2")) {
          FFHistUpAdded->Add(FFHist2_);
          FFHistDnAdded->Add(FFHist2_);

          FFHist2UpAdded = (TH1D*)syst_.at("FFHist2").up_->Clone();
          FFHist2DnAdded = (TH1D*)syst_.at("FFHist2").dn_->Clone();
          FFHist2UpAdded->Add(FFHist_);
          FFHist2DnAdded->Add(FFHist_);
        }

        if (syst_.count("FFHistAbcdMergedEleEnCorr")) {
          FFHistPreCorrUpAdded = (TH1D*)syst_.at("FFHistAbcdMergedEleEnCorr").up_->Clone();
          FFHistPreCorrDnAdded = (TH1D*)syst_.at("FFHistAbcdMergedEleEnCorr").dn_->Clone();
        }

        if (syst_.count("FFHist2abcdMergedEleEnCorr")) {
          FFHistPreCorrUpAdded->Add((TH1D*)syst_.at("FFHist2abcdMergedEleEnCorr").up_);
          FFHistPreCorrDnAdded->Add((TH1D*)syst_.at("FFHist2abcdMergedEleEnCorr").dn_);
        }

        for (unsigned idx = 1; idx <= FFadded->GetNbinsX(); idx++) {
          x0.push_back(FFadded->GetBinCenter(idx));
          y0.push_back(FFadded->GetBinContent(idx));
          errx.push_back(FFadded->GetBinWidth(idx)/2.);

          double valFFUp = syst_.at("FFHist").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
          double valFFDn = FFHist_->GetBinContent(idx) - syst_.at("FFHist").dn_->GetBinContent(idx);
          double valFF2Up = 0.;
          double valFF2Dn = 0.;
          double valFFpreCorrUp = 0.;
          double valFFpreCorrDn = 0.;
          double normFF = 0.2*FFHist_->GetBinContent(idx);
          double normFF2 = 0.;

          if (syst_.count("FFHistAbcdMergedEleEnCorr")) {
            valFFpreCorrUp = syst_.at("FFHistAbcdMergedEleEnCorr").up_->GetBinContent(idx) - FFHist_->GetBinContent(idx);
            valFFpreCorrDn = FFHist_->GetBinContent(idx) - syst_.at("FFHistAbcdMergedEleEnCorr").dn_->GetBinContent(idx);
          }

          if (syst_.count("FFHist2")) {
            valFF2Up = syst_.at("FFHist2").up_->GetBinContent(idx) - FFHist2_->GetBinContent(idx);
            valFF2Dn = FFHist2_->GetBinContent(idx) - syst_.at("FFHist2").dn_->GetBinContent(idx);
            normFF2 = 0.2*FFHist2_->GetBinContent(idx);

            if (syst_.count("FFHist2abcdMergedEleEnCorr")) {
              valFFpreCorrUp += syst_.at("FFHist2abcdMergedEleEnCorr").up_->GetBinContent(idx) - FFHist2_->GetBinContent(idx);
              valFFpreCorrDn += FFHist2_->GetBinContent(idx) - syst_.at("FFHist2abcdMergedEleEnCorr").dn_->GetBinContent(idx);
            }
          }

          erryUp.push_back( std::sqrt( (valFFUp+valFF2Up)*(valFFUp+valFF2Up) + valFFpreCorrUp*valFFpreCorrUp + (normFF+normFF2)*(normFF+normFF2) ) );
          erryDn.push_back( std::sqrt( (valFFDn+valFF2Dn)*(valFFDn+valFF2Dn) + valFFpreCorrDn*valFFpreCorrDn + (normFF+normFF2)*(normFF+normFF2) ) );

          if (ratio) {
            r0.push_back(1.);

            double rvalFFUp = FFHistUpAdded->GetBinContent(idx)/FFadded->GetBinContent(idx) - 1.;
            double rvalFFDn = 1. - FFHistDnAdded->GetBinContent(idx)/FFadded->GetBinContent(idx);
            double rvalFF2Up = 0.;
            double rvalFF2Dn = 0.;
            double rvalFFpreCorrUp = 0.;
            double rvalFFpreCorrDn = 0.;
            double rnormFF = normFF/FFadded->GetBinContent(idx);
            double rnormFF2 = 0.;

            if (syst_.count("FFHistAbcdMergedEleEnCorr")) {
              rvalFFpreCorrUp = FFHistPreCorrUpAdded->GetBinContent(idx)/FFadded->GetBinContent(idx) - 1.;
              rvalFFpreCorrDn = 1. - FFHistPreCorrDnAdded->GetBinContent(idx)/FFadded->GetBinContent(idx);
            }

            if (syst_.count("FFHist2")) {
              rvalFF2Up = FFHist2UpAdded->GetBinContent(idx)/FFadded->GetBinContent(idx) - 1.;
              rvalFF2Dn = 1. - FFHist2DnAdded->GetBinContent(idx)/FFadded->GetBinContent(idx);
              rnormFF2 = normFF2/FFadded->GetBinContent(idx);
            }

            double rUp = std::sqrt( (rvalFFUp+rvalFF2Up)*(rvalFFUp+rvalFF2Up) + rvalFFpreCorrUp*rvalFFpreCorrUp + (rnormFF+rnormFF2)*(rnormFF+rnormFF2) );
            double rDn = std::sqrt( (rvalFFDn+rvalFF2Dn)*(rvalFFDn+rvalFF2Dn) + rvalFFpreCorrDn*rvalFFpreCorrDn + (rnormFF+rnormFF2)*(rnormFF+rnormFF2) );

            errRup.push_back( FFadded->GetBinContent(idx) > 0. ? rUp : 0. );
            errRdn.push_back( FFadded->GetBinContent(idx) > 0. ? rDn : 0. );
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
          syst_.at("FFHist2").up_->Add(syst_.at("FFHist").up_);
          syst_.at("FFHist2").dn_->Add(syst_.at("FFHist").dn_);
          syst_.at("FFHist2abcdMergedEleEnCorr").up_->Add(syst_.at("FFHistAbcdMergedEleEnCorr").up_);
          syst_.at("FFHist2abcdMergedEleEnCorr").dn_->Add(syst_.at("FFHistAbcdMergedEleEnCorr").dn_);

          dir_->WriteTObject(dataHist_,"data_obs");
          // dir_->WriteTObject(FFHist_,"SS");
          dir_->WriteTObject(FFHist2_,"OS");

          // dir_->WriteTObject(syst_.at("FFHist").up_,"SS_mergedEleFakeFactorSSUp");
          // dir_->WriteTObject(syst_.at("FFHist").dn_,"SS_mergedEleFakeFactorSSDown");
          dir_->WriteTObject(syst_.at("FFHist2").up_,"OS_mergedEleFakeFactorOSUp");
          dir_->WriteTObject(syst_.at("FFHist2").dn_,"OS_mergedEleFakeFactorOSDown");
          // dir_->WriteTObject(syst_.at("FFHistAbcdMergedEleEnCorr").up_,"SS_mergedEleEnCorrUp");
          // dir_->WriteTObject(syst_.at("FFHistAbcdMergedEleEnCorr").dn_,"SS_mergedEleEnCorrDown");
          dir_->WriteTObject(syst_.at("FFHist2abcdMergedEleEnCorr").up_,"OS_mergedEleEnCorrUp");
          dir_->WriteTObject(syst_.at("FFHist2abcdMergedEleEnCorr").dn_,"OS_mergedEleEnCorrDown");

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
            dir_->WriteTObject(sigSyst_.at(idx).at("sigElTrig").up_,sigFiles_.at(idx).name_+"_elTrigUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigElTrig").dn_,sigFiles_.at(idx).name_+"_elTrigDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").up_,sigFiles_.at(idx).name_+"_elRecoUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").dn_,sigFiles_.at(idx).name_+"_elRecoDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").up_,sigFiles_.at(idx).name_+"_PUrwgtUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").dn_,sigFiles_.at(idx).name_+"_PUrwgtDown");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").up_,sigFiles_.at(idx).name_+"_prefireUp");
            dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").dn_,sigFiles_.at(idx).name_+"_prefireDown");
          }
        }
      }
    }
  };

  // invM
  auto aloader2E = HistLoader2E(datafile,sigsamples);
  auto aloader2E1 = HistLoader2E(datafile1,sigsamples1);
  auto aloader2E2 = HistLoader2E(datafile2,sigsamples2);
  auto aloader2E3 = HistLoader2E(datafile3,sigsamples3);

/*  aloader2E.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,postfix.Data(),"",false);
  aloader2E1.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL16APV");
  aloader2E2.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL17");
  aloader2E3.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",kCyan+1,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);

  aloader2E.compare(p1,true,p2,2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiSS.png",p1);*/

/*  aloader2E.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,postfix.Data(),"",true);
  aloader2E1.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,"20UL16APV");
  aloader2E2.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,"20UL17");
  aloader2E3.load("2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",kOrange,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_mixedOS.png",p1);*/

  p1->SetLogy();
/*  aloader2E.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,postfix.Data(),"",true);
  aloader2E1.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,"20UL16APV");
  aloader2E2.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,"20UL17");
  aloader2E3.load("2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",kOrange,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,false,p2,5);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiOS.png",p1);*/

  //p1->SetLogx();
  //p2->SetLogx();
  p1->SetLogy();
  aloader2E.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",TColor::GetColor("#f89c20"),postfix.Data(),"",false,"2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",TColor::GetColor("#f89c20"));
  aloader2E1.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",TColor::GetColor("#f89c20"),"20UL16APV","",false,"2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",TColor::GetColor("#f89c20"));
  aloader2E2.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",TColor::GetColor("#f89c20"),"20UL17","",false,"2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",TColor::GetColor("#f89c20"));
  aloader2E3.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",TColor::GetColor("#f89c20"),"20UL18","",false,"2E_CRME_OSll_invM","2E_mixedME_OSll_invM_xFF",TColor::GetColor("#f89c20"));
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  //aloader2E.preparecard("ME2E_"+era+"_datacard.root","mergedEle2E");
  aloader2E.compare(p1,false,p2,5); // 5
  SaveAs(canvas_2,"MEFF_2E_invM_denom_CR.pdf",p1);
  //aloader2E.close();

  p1->SetLogy();
  aloader2E.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",TColor::GetColor("#f89c20"),postfix.Data(),"",false,"2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",TColor::GetColor("#f89c20"));
  aloader2E1.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",TColor::GetColor("#f89c20"),"20UL16APV","",false,"2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",TColor::GetColor("#f89c20"));
  aloader2E2.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",TColor::GetColor("#f89c20"),"20UL17","",false,"2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",TColor::GetColor("#f89c20"));
  aloader2E3.load("2E_mixedME_SSll_invM","2E_antiME_SSll_invM_CR_xFF",TColor::GetColor("#f89c20"),"20UL18","",false,"2E_mixedME_OSll_invM","2E_antiME_OSll_invM_CR_xFF",TColor::GetColor("#f89c20"));
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,false,p2,5);
  SaveAs(canvas_2,"MEFF_2E_invM_mixed_CR.pdf",p1);

  //p1->SetLogx(0);
  //p2->SetLogx(0);

/*  aloader2E.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,postfix.Data(),"",false);
  aloader2E1.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL16APV");
  aloader2E2.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL17");
  aloader2E3.load("2E_CRME_SSll_invM","2E_mixedME_SSll_invM_xFF",kCyan+1,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedSS.png");*/

  // Et
/*  canvas_1->SetLogy();
  aloader2E.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,postfix.Data());
  aloader2E1.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,"20UL16APV");
  aloader2E2.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,"20UL17");
  aloader2E3.load("2E_Et_SSCR_EB_CRME","2E_Et_SSCR_EB_mixedAntiME_xFF",kCyan+1,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedSS.png");

  aloader2E.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,postfix.Data());
  aloader2E1.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,"20UL16APV");
  aloader2E2.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,"20UL17");
  aloader2E3.load("2E_Et_SSCR_EB_mixedME","2E_Et_SSCR_EB_antiME_xFF",kCyan+1,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiSS.png",p1);

  aloader2E.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,postfix.Data());
  aloader2E1.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,"20UL16APV");
  aloader2E2.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,"20UL17");
  aloader2E3.load("2E_Et_OSCR_EB_CRME","2E_Et_OSCR_EB_mixedAntiME_xFF",kOrange,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedOS.png");

  aloader2E.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,postfix.Data());
  aloader2E1.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,"20UL16APV");
  aloader2E2.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,"20UL17");
  aloader2E3.load("2E_Et_OSCR_EB_mixedME","2E_Et_OSCR_EB_antiME_xFF",kOrange,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,false,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiOS.png",p1);
  canvas_1->SetLogy(0);

  // eta
  p1->SetLogy(0);
  aloader2E.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,postfix.Data());
  aloader2E1.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,"20UL16APV");
  aloader2E2.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,"20UL17");
  aloader2E3.load("2E_CRME_SSll_eta","2E_mixedAntiME_SSll_eta_xFF",kCyan+1,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedSS.png");

  aloader2E.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,postfix.Data());
  aloader2E1.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,"20UL16APV");
  aloader2E2.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,"20UL17");
  aloader2E3.load("2E_mixedME_SSll_eta","2E_antiME_SSll_eta_xFF",kCyan+1,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,true,p2,2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiSS.png",p1);

  aloader2E.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,postfix.Data());
  aloader2E1.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,"20UL16APV");
  aloader2E2.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,"20UL17");
  aloader2E3.load("2E_CRME_OSll_eta","2E_mixedAntiME_OSll_eta_xFF",kOrange,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(canvas_1);
  SaveAs(canvas_1,"FF_2E_eta_denom_mixedOS.png");

  aloader2E.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,postfix.Data());
  aloader2E1.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,"20UL16APV");
  aloader2E2.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,"20UL17");
  aloader2E3.load("2E_mixedME_OSll_eta","2E_antiME_OSll_eta_xFF",kOrange,"20UL18");
  aloader2E.add(aloader2E1);
  aloader2E.add(aloader2E2);
  aloader2E.add(aloader2E3);
  aloader2E.compare(p1,true,p2);
  SaveAs(canvas_2,"FF_2E_eta_denom_antiOS.png",p1);*/

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

    public:
      TH1D* SSFFHist_ = nullptr;
      TH1D* SSWZHist_ = nullptr;
      TH1D* SSZZHist_ = nullptr;
      TH1D* OSWZHist_ = nullptr;
      TH1D* OSZZHist_ = nullptr;
      TH1D* OSZtagHist_ = nullptr;

    public:
      void add(const denomHists& other) {
        this->SSFFHist_->Add(other.SSFFHist_);
        this->SSWZHist_->Add(other.SSWZHist_);
        this->SSZZHist_->Add(other.SSZZHist_);
        this->OSWZHist_->Add(other.OSWZHist_);
        this->OSZZHist_->Add(other.OSZZHist_);
        this->OSZtagHist_->Add(other.OSZtagHist_);
      }

      void load(TString denomNameSS, TString denomNameOS, TString denomNameZtag, TFile* datafile, TFile* WZfile, TFile* ZZfile, std::string anlyzrEra="") {
        SSFFHist_ = (TH1D*)datafile->Get(anlyzrData+"/"+denomNameSS)->Clone();
        OSZtagHist_ = (TH1D*)datafile->Get(anlyzrData+"/"+denomNameZtag)->Clone();

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
        OSZtagHist_->Rebin(r);
      }

      TH1D* SSFFHist() { return SSFFHist_; }
      TH1D* SSWZHist() { return SSWZHist_; }
      TH1D* SSZZHist() { return SSZZHist_; }
      TH1D* OSWZHist() { return OSWZHist_; }
      TH1D* OSZZHist() { return OSZZHist_; }
      TH1D* OSZtagHist() { return OSZtagHist_; }

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
        denomOSfinal->Add(OSZtagHist_);
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
    denomHists preCorrUp_;
    denomHists preCorrDn_;
    denomHists SSFFup_;
    denomHists SSFFdn_;
    denomHists OSFFup_;
    denomHists OSFFdn_;
    denomHists heepIdUp_;
    denomHists heepIdDn_;
    denomHists elRecoUp_;
    denomHists elRecoDn_;
    denomHists elTrigUp_;
    denomHists elTrigDn_;

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
        this->preCorrUp_.add(other.preCorrUp_);
        this->preCorrDn_.add(other.preCorrDn_);
        this->elRecoUp_.add(other.elRecoUp_);
        this->elRecoDn_.add(other.elRecoDn_);
        this->elTrigUp_.add(other.elTrigUp_);
        this->elTrigDn_.add(other.elTrigDn_);
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
      nominal_.load(denomName+"_xSSFF",denomName+"_xOSFF",denomName+"_Ztag_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
      SSFFup_.load(denomName+"_xSSFF_up",denomName+"_xOSFF",denomName+"_Ztag_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
      SSFFdn_.load(denomName+"_xSSFF_dn",denomName+"_xOSFF",denomName+"_Ztag_xOSFF",datafile_,WZfile_,ZZfile_,anlyzrEra);
      OSFFup_.load(denomName+"_xSSFF",denomName+"_xOSFF_up",denomName+"_Ztag_xOSFF_up",datafile_,WZfile_,ZZfile_,anlyzrEra);
      OSFFdn_.load(denomName+"_xSSFF",denomName+"_xOSFF_dn",denomName+"_Ztag_xOSFF_dn",datafile_,WZfile_,ZZfile_,anlyzrEra);

      sigHist_.clear();
      sigSyst_.clear();

      if (numName.Contains("invM")) {
        const double lumi = retrieveLumi(anlyzrEra);
        const double sigLumi = 0.01;

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
          TH1D* sigElTrigUp = retrieveSigHist(sigFiles_.at(idx).file_,"_elTrigUp");
          TH1D* sigElTrigDn = retrieveSigHist(sigFiles_.at(idx).file_,"_elTrigDn");
          init["sigElTrig"] = SystVariation(sigElTrigUp,sigElTrigDn);
          TH1D* sigElRecoUp = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoUp");
          TH1D* sigElRecoDn = retrieveSigHist(sigFiles_.at(idx).file_,"_elRecoDn");
          init["sigElReco"] = SystVariation(sigElRecoUp,sigElRecoDn);
          TH1D* sigPUrwgtUp = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtUp");
          TH1D* sigPUrwgtDn = retrieveSigHist(sigFiles_.at(idx).file_,"_PUrwgtDn");
          init["sigPUrwgt"] = SystVariation(sigPUrwgtUp,sigPUrwgtDn);
          TH1D* sigPrefireUp = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireUp");
          TH1D* sigPrefireDn = retrieveSigHist(sigFiles_.at(idx).file_,"_prefireDn");
          init["sigPrefire"] = SystVariation(sigPrefireUp,sigPrefireDn);

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

        heepIdUp_.load(denomName+"_xSSFF_heepIdUp",denomName+"_xOSFF_heepIdUp",denomName+"_Ztag_xOSFF_heepIdUp",datafile_,WZfile_,ZZfile_,anlyzrEra);
        heepIdDn_.load(denomName+"_xSSFF_heepIdDn",denomName+"_xOSFF_heepIdDn",denomName+"_Ztag_xOSFF_heepIdDn",datafile_,WZfile_,ZZfile_,anlyzrEra);
        elRecoUp_.load(denomName+"_xSSFF_elRecoUp",denomName+"_xOSFF_elRecoUp",denomName+"_Ztag_xOSFF_elRecoUp",datafile_,WZfile_,ZZfile_,anlyzrEra);
        elRecoDn_.load(denomName+"_xSSFF_elRecoDn",denomName+"_xOSFF_elRecoDn",denomName+"_Ztag_xOSFF_elRecoDn",datafile_,WZfile_,ZZfile_,anlyzrEra);
        elTrigUp_.load(denomName+"_xSSFF_elTrigUp",denomName+"_xOSFF_elTrigUp",denomName+"_Ztag_xOSFF_elTrigUp",datafile_,WZfile_,ZZfile_,anlyzrEra);
        elTrigDn_.load(denomName+"_xSSFF_elTrigDn",denomName+"_xOSFF_elTrigDn",denomName+"_Ztag_xOSFF_elTrigDn",datafile_,WZfile_,ZZfile_,anlyzrEra);
        preCorrUp_.load(denomName+"_xSSFF_preCorr",denomName+"_xOSFF_preCorr",denomName+"_Ztag_xOSFF_preCorr",datafile_,WZfile_,ZZfile_,anlyzrEra);
        preCorrDn_.SSFFHist_ = variateDn(nominal_.SSFFHist(),preCorrUp_.SSFFHist());
        preCorrDn_.SSWZHist_ = variateDn(nominal_.SSWZHist(),preCorrUp_.SSWZHist());
        preCorrDn_.SSZZHist_ = variateDn(nominal_.SSZZHist(),preCorrUp_.SSZZHist());
        preCorrDn_.OSWZHist_ = variateDn(nominal_.OSWZHist(),preCorrUp_.OSWZHist());
        preCorrDn_.OSZZHist_ = variateDn(nominal_.OSZZHist(),preCorrUp_.OSZZHist());
        preCorrDn_.OSZtagHist_ = variateDn(nominal_.OSZtagHist(),preCorrUp_.OSZtagHist());
      }
    }

    void compare(TPad* padUp, TPad* padDn, int rebin=1) {
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

      auto denomSSfinal = nominal_.returnSS();
      auto denomOSfinal = nominal_.returnOS();

      auto denomfinal_nominal = nominal_.returnAdded();
      auto denomfinal_SSFFup = SSFFup_.returnAdded();
      auto denomfinal_SSFFdn = SSFFdn_.returnAdded();
      auto denomfinal_OSFFup = OSFFup_.returnAdded();
      auto denomfinal_OSFFdn = OSFFdn_.returnAdded();
      std::unique_ptr<TH1D> denomfinal_heepIdUp, denomfinal_heepIdDn, denomfinal_precorrUp, denomfinal_precorrDn;
      std::unique_ptr<TH1D> denomfinal_elRecoUp, denomfinal_elRecoDn, denomfinal_elTrigUp, denomfinal_elTrigDn;

      if ( TString(dataHist_->GetName()).Contains("invM") ) {
        heepIdUp_.rebin(heepIdUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        heepIdDn_.rebin(heepIdDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        elRecoUp_.rebin(elRecoUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        elRecoDn_.rebin(elRecoDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        elTrigUp_.rebin(elTrigUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        elTrigDn_.rebin(elTrigDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        preCorrUp_.rebin(preCorrUp_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        preCorrDn_.rebin(preCorrDn_.SSFFHist()->GetNbinsX()/dataHist_->GetNbinsX());
        denomfinal_heepIdUp = heepIdUp_.returnAdded();
        denomfinal_heepIdDn = heepIdDn_.returnAdded();
        denomfinal_elRecoUp = elRecoUp_.returnAdded();
        denomfinal_elRecoDn = elRecoDn_.returnAdded();
        denomfinal_elTrigUp = elTrigUp_.returnAdded();
        denomfinal_elTrigDn = elTrigDn_.returnAdded();
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
        dataHist_->GetXaxis()->SetTitle("M_{3e} [GeV]");
        //dataHist_->SetMaximum(5.*dataHist_->GetMaximum());
        //dataHist_->SetMinimum(0.2);
      } else if (TString(dataHist_->GetName()).Contains("Et")) {
        dataHist_->GetXaxis()->SetTitle("E_{T}^{5#times5} [GeV]");
        dataHist_->GetXaxis()->SetRangeUser(0.,250.);
      } else if (TString(dataHist_->GetName()).Contains("eta")) {
        dataHist_->GetXaxis()->SetTitle("#eta^{5#times5}");
        dataHist_->SetMaximum(1.1*dataHist_->GetMaximum());
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
          legend->AddEntry(sigHist_.at(0),"X750Y10");
          //legend_left->AddEntry(sigNums.at(isigDiv),"M_{A} = 10 GeV (#sigma = 10 fb)");
        }

        legend->Draw();
      }

      TH1D* ratio = nullptr;

      if (padDn) {
        ratio = (TH1D*)dataHist_->Clone();
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
        double normSSFF = 0.; // 0.2*nominal_.returnSS()->GetBinContent(idx);
        double normOSFF = 0.2*nominal_.returnOS()->GetBinContent(idx);

        double valHeepIdUp = 0., valHeepIdDn = 0., valPreCorrUp = 0., valPreCorrDn = 0.;
        double valElRecoUp = 0., valElRecoDn = 0., valElTrigUp = 0., valElTrigDn = 0.;

        if ( TString(dataHist_->GetName()).Contains("invM") ) {
          valHeepIdUp = denomfinal_heepIdUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valHeepIdDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_heepIdDn->GetBinContent(idx);
          valPreCorrUp = denomfinal_precorrUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valPreCorrDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_precorrDn->GetBinContent(idx);
          valElRecoUp = denomfinal_elRecoUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valElRecoDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_elRecoDn->GetBinContent(idx);
          valElTrigUp = denomfinal_elTrigUp->GetBinContent(idx) - denomfinal_nominal->GetBinContent(idx);
          valElTrigDn = denomfinal_nominal->GetBinContent(idx) - denomfinal_elTrigDn->GetBinContent(idx);
        }

        erryUp.push_back( std::sqrt( valSSFFup*valSSFFup + valOSFFup*valOSFFup + valHeepIdUp*valHeepIdUp + valPreCorrUp*valPreCorrUp + valElRecoUp*valElRecoUp + valElTrigUp*valElTrigUp
                                     + normSSFF*normSSFF + normOSFF*normOSFF ) );
        erryDn.push_back( std::sqrt( valSSFFdn*valSSFFdn + valOSFFdn*valOSFFdn + valHeepIdDn*valHeepIdDn + valPreCorrDn*valPreCorrDn + valElRecoDn*valElRecoDn + valElTrigDn*valElTrigDn
                                     + normSSFF*normSSFF + normOSFF*normOSFF ) );

        if (ratio) {
          r0.push_back(1.);

          double rvalSSFFup = denomfinal_SSFFup->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
          double rvalSSFFdn = 1. - denomfinal_SSFFdn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
          double rvalOSFFup = denomfinal_OSFFup->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
          double rvalOSFFdn = 1. - denomfinal_OSFFdn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
          double rnormSSFF = normSSFF/denomfinal_nominal->GetBinContent(idx);
          double rnormOSFF = normOSFF/denomfinal_nominal->GetBinContent(idx);

          double rHeepIdUp = 0., rHeepIdDn = 0., rPreCorrUp = 0., rPreCorrDn = 0.;
          double rElRecoUp = 0., rElRecoDn = 0., rElTrigUp = 0., rElTrigDn = 0.;

          if ( TString(dataHist_->GetName()).Contains("invM") ) {
            rHeepIdUp = denomfinal_heepIdUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rHeepIdDn = 1. - denomfinal_heepIdDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
            rPreCorrUp = denomfinal_precorrUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rPreCorrDn = 1. - denomfinal_precorrDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
            rElRecoUp = denomfinal_elRecoUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rElRecoDn = 1. - denomfinal_elRecoDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
            rElTrigUp = denomfinal_elTrigUp->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx) -1.;
            rElTrigDn = 1. - denomfinal_elTrigDn->GetBinContent(idx)/denomfinal_nominal->GetBinContent(idx);
          }

          double rUp = std::sqrt( rvalSSFFup*rvalSSFFup + rvalOSFFup*rvalOSFFup + rHeepIdUp*rHeepIdUp
                                       + rPreCorrUp*rPreCorrUp + rElRecoUp*rElRecoUp + rElTrigUp*rElTrigUp + rnormSSFF*rnormSSFF + rnormOSFF*rnormOSFF );
          double rDn = std::sqrt( rvalSSFFdn*rvalSSFFdn + rvalOSFFdn*rvalOSFFdn + rHeepIdDn*rHeepIdDn
                                       + rPreCorrDn*rPreCorrDn + rElRecoDn*rElRecoDn + rElTrigDn*rElTrigDn + rnormSSFF*rnormSSFF + rnormOSFF*rnormOSFF );

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
        dir_->WriteTObject(preCorrUp_.returnSS().release(),"SS_mergedEleEnCorrUp");
        dir_->WriteTObject(preCorrUp_.returnOS().release(),"OS_mergedEleEnCorrUp");
        dir_->WriteTObject(preCorrDn_.returnSS().release(),"SS_mergedEleEnCorrDown");
        dir_->WriteTObject(preCorrDn_.returnOS().release(),"OS_mergedEleEnCorrDown");
        dir_->WriteTObject(elRecoUp_.returnSS().release(),"SS_elRecoUp");
        dir_->WriteTObject(elRecoUp_.returnOS().release(),"OS_elRecoUp");
        dir_->WriteTObject(elRecoDn_.returnSS().release(),"SS_elRecoDown");
        dir_->WriteTObject(elRecoDn_.returnOS().release(),"OS_elRecoDown");
        dir_->WriteTObject(elTrigUp_.returnSS().release(),"SS_elTrigUp");
        dir_->WriteTObject(elTrigUp_.returnOS().release(),"OS_elTrigUp");
        dir_->WriteTObject(elTrigDn_.returnSS().release(),"SS_elTrigDown");
        dir_->WriteTObject(elTrigDn_.returnOS().release(),"OS_elTrigDown");

        for (unsigned idx=0; idx<sigHist_.size(); idx++) {
          dir_->WriteTObject(sigHist_.at(idx),sigFiles_.at(idx).name_);
          dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").up_,sigFiles_.at(idx).name_+"_modHeepIdUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigModHeep").dn_,sigFiles_.at(idx).name_+"_modHeepIdDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElTrig").up_,sigFiles_.at(idx).name_+"_elTrigUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElTrig").dn_,sigFiles_.at(idx).name_+"_elTrigDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").up_,sigFiles_.at(idx).name_+"_elRecoUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigElReco").dn_,sigFiles_.at(idx).name_+"_elRecoDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").up_,sigFiles_.at(idx).name_+"_PUrwgtUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPUrwgt").dn_,sigFiles_.at(idx).name_+"_PUrwgtDown");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").up_,sigFiles_.at(idx).name_+"_prefireUp");
          dir_->WriteTObject(sigSyst_.at(idx).at("sigPrefire").dn_,sigFiles_.at(idx).name_+"_prefireDown");
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

  canvas_2->cd();
  p1->SetLogy(0);

  auto aloader3E = HistLoader3E(datafile,WZfile,ZZfile,sigsamples3E);
  auto aloader3E1 = HistLoader3E(datafile1,WZfile1,ZZfile1,sigsamples3E1);
  auto aloader3E2 = HistLoader3E(datafile2,WZfile2,ZZfile2,sigsamples3E2);
  auto aloader3E3 = HistLoader3E(datafile3,WZfile3,ZZfile3,sigsamples3E3);

  //canvas_1->SetLogy();
  aloader3E.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR",postfix.Data());
  aloader3E1.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR","20UL16APV");
  aloader3E2.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR","20UL17");
  aloader3E3.load("3E_CRME_lll_invM","3E_antiME_lll_invM_CR","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  //aloader3E.preparecard("ME3E_"+era+"_datacard.root","mergedEle3E");
  aloader3E.compare(p1,p2,5); // 5
  SaveAs(canvas_2,"MEFF_3E_invM.pdf",p1);
  //aloader3E.close();
  //canvas_1->SetLogy(0);

  canvas_1->cd();

  canvas_1->SetLogy(0);
  aloader3E.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME",postfix.Data());
  aloader3E1.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME","20UL16APV");
  aloader3E2.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME","20UL17");
  aloader3E3.load("3E_Et_CR_EB_CRME","3E_Et_CR_EB_antiME","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(p1,p2);
  SaveAs(canvas_2,"MEFF_3E_Et.pdf",p1);
  canvas_1->SetLogy(0);

  aloader3E.load("3E_eta_CR_EB_CRME","3E_antiME_eta",postfix.Data());
  aloader3E1.load("3E_eta_CR_EB_CRME","3E_antiME_eta","20UL16APV");
  aloader3E2.load("3E_eta_CR_EB_CRME","3E_antiME_eta","20UL17");
  aloader3E3.load("3E_eta_CR_EB_CRME","3E_antiME_eta","20UL18");
  aloader3E.add(aloader3E1);
  aloader3E.add(aloader3E2);
  aloader3E.add(aloader3E3);
  aloader3E.compare(p1,p2,2);
  SaveAs(canvas_2,"MEFF_3E_eta.pdf",p1);

  return;
}
