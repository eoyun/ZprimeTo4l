#include "TROOT.h"
#include <cstdarg>
#include <cmath>
#include <algorithm>

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

static void fillHist(TTree* atree, TH1D* ahist, const TString& branch, TString multiply="", TF1* fit=nullptr, TFitResultPtr result=TFitResultPtr(-1), int updn=0) {
  float val = 0.;
  float wgt = 0.;
  float ffx = 0.;
  float e1pt, e1eta, e1phi, e2pt, e2eta, e2phi;
  float invm = 0.;

  if (branch!="e1e2Pt")
    atree->SetBranchAddress(branch.Data(),&val);

  if (multiply!="" && multiply!="e1e2Pt")
    atree->SetBranchAddress(multiply,&ffx);

  if (branch!="e1e2InvM")
    atree->SetBranchAddress("e1e2InvM",&invm);

  atree->SetBranchAddress("e1Pt",&e1pt);
  atree->SetBranchAddress("e1Eta",&e1eta);
  atree->SetBranchAddress("e1Phi",&e1phi);
  atree->SetBranchAddress("e2Pt",&e2pt);
  atree->SetBranchAddress("e2Eta",&e2eta);
  atree->SetBranchAddress("e2Phi",&e2phi);
  atree->SetBranchAddress("wgt",&wgt);

  for (int idx=0; idx < atree->GetEntries(); ++idx) {
    atree->GetEntry(idx);

    auto lvec1 = ROOT::Math::PtEtaPhiMVector(e1pt,e1eta,e1phi,0.);
    auto lvec2 = ROOT::Math::PtEtaPhiMVector(e2pt,e2eta,e2phi,0.);

    if (branch=="e1e2Pt")
      val = (lvec1+lvec2).Pt();

    if (multiply=="e1e2Pt")
      ffx = (lvec1+lvec2).Pt();

    if (branch=="e1Pt")
      val = e1pt;
    if (branch=="e2Pt")
      val = e2pt;

    double ff = 1.;

    if (branch=="e1e2InvM")
      invm = val;

    if (invm < 1.)
      continue;

    if (fit) {
      ff = fit->Eval( ffx );

      double xcen[1] = {ffx};
      double ci[1];
      result->GetConfidenceIntervals(1,1,0,xcen,ci,0.68,false);

      if (updn==1)
        ff = ff+ci[0];
      else if (updn==-1)
        ff = std::max(0.,ff-ci[0]);
    }

    ahist->Fill(val,wgt*ff);
  }

  atree->ResetBranchAddresses();

  return;
}

void estimateREFFdr03(TString era) {
  TH1::SetDefaultSumw2();

  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Internal";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

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
    lumi_13TeV = "137.6 fb^{-1}";
  } else {
    std::cout << "check era..." << std::endl;
  }

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

  auto canvas2 = std::make_unique<TCanvas>("canvas2","canvas2",800,800);

  canvas2->SetFillColor(0);
  canvas2->SetBorderMode(0);
  canvas2->SetFrameFillStyle(0);
  canvas2->SetFrameBorderMode(0);
  canvas2->SetLeftMargin( L/W );
  canvas2->SetRightMargin( R/W );
  canvas2->SetTopMargin( T/H );
  canvas2->SetBottomMargin( B/H );
  canvas2->SetTickx(0);
  canvas2->SetTicky(0);

  auto SaveAs = [&] (TCanvas* cv, const std::string& name, TPad* pad = nullptr) {
    cv->Update();

    // writing the lumi information and the CMS "logo"
    CMS_lumi( cv, iPeriod, iPos );

    if (pad) {
      pad->RedrawAxis();
      pad->GetFrame()->Draw();
    } else {
      cv->Update();
      cv->RedrawAxis();
      cv->GetFrame()->Draw();
    }

    cv->SaveAs(name.c_str());
  };

  class sample {
  public:
    TFile* file_;
    double xsec_;
    double lumi_;
    int color_;

    sample(TFile* afile, const double axsec, const double alumi, const int acolor) {
      file_ = afile;
      xsec_ = axsec;
      lumi_ = alumi;
      color_ = acolor;
    }
  };

  class samples {
  public:
    TString era_;
    double lumi_;
    std::vector<sample> samples_;
    std::vector<TH1D*> hists_;

    samples(TString aera, const double alumi) {
      era_ = aera;
      lumi_ = alumi;

      samples_ = {
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ZZ.root"),13.81,lumi_,TColor::GetColor("#7a21dd")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WZFXFX.root"),5.213,lumi_,TColor::GetColor("#9c9ca1")), // x0.65
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WW.root"),11.09	,lumi_,TColor::GetColor("#964a8b")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_tW_antitop.root"),79.3*0.5,lumi_,kGreen+2), // XSDB 32.51
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_tW_top.root"),79.3*0.5,lumi_,kGreen+2), // XSDB 32.45
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_t-ch_antitop.root"),80.0,lumi_,kGreen-3), // XSDB 71.75
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_t-ch_top.root"),134.2,lumi_,kGreen-3), // XSDB 119.7
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_ST_s-ch.root"),6.839,lumi_,kGreen-6), // XSDB 3.549 // NNLO 6.839
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_TTTo2L2Nu.root"),833.9*0.1062,lumi_,TColor::GetColor("#e42536")), // 88.29
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_TTsemi.root"),833.9*0.4394,lumi_,kRed-7),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_TThad.root"),833.9*0.4544,lumi_,kRed-9),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_DY_2J.root"),353.6,lumi_,TColor::GetColor("#f89c20")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_DY_1J.root"),983.5,lumi_,TColor::GetColor("#f89c20")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_DY_0J.root"),5090.0,lumi_,TColor::GetColor("#f89c20")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WJets_2J.root"),3276.0,lumi_,TColor::GetColor("#5790fc")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WJets_1J.root"),8832.0,lumi_,TColor::GetColor("#5790fc")),
        sample(TFile::Open("FakeEleAnalyzer_"+era_+"_WJets_0J.root"),52780.0	,lumi_,TColor::GetColor("#5790fc")),
      };
    }

    void retrieveHists(const TString& anlyzr, const TString& aname, const TString& branch, double xmin, double xmax, unsigned nbin=100, int rebin=1, TString multiply="", TF1* fit=nullptr, TFitResultPtr result=TFitResultPtr(-1), int updn=0) {
      for (auto ahist : hists_)
        delete ahist;

      hists_.clear();

      TString postfix = era_;

      for (const auto& el : samples_) {
        const double sumwgt = ((TH1D*)el.file_->Get("evtCounter/h_sumW"))->GetBinContent(1);

        TTree* atree = (TTree*)el.file_->Get(anlyzr+postfix+"/"+aname);
        TH1D* ahist = new TH1D( aname + branch + TString(el.file_->GetName())+TString(std::to_string(updn).c_str()), aname+branch, nbin, xmin, xmax );

        fillHist(atree,ahist,branch,multiply,fit,result,updn);

        hists_.push_back(ahist);
        hists_.back()->SetFillColor(el.color_);
        hists_.back()->SetLineColor(el.color_);
        hists_.back()->Scale(lumi_*1000.*el.xsec_/sumwgt);
        hists_.back()->Rebin(rebin);

        delete atree;
      }
    }

    void addHists(const samples& other) {
      for (unsigned idx=0; idx<hists_.size(); idx++)
        hists_.at(idx)->Add(other.hists_.at(idx));
    }

    std::unique_ptr<THStack> extractStack() {
      TString aformat;
      aformat.Form("%s;%s;%s",hists_.back()->GetTitle(),hists_.back()->GetXaxis()->GetTitle(),hists_.back()->GetYaxis()->GetTitle());
      auto astack = std::make_unique<THStack>("stack",aformat);

      for (auto* ahist : hists_)
        astack->Add(ahist);

      return std::move(astack);
    }

    std::unique_ptr<TH1D> extractNominal() {
      auto ahist = std::unique_ptr<TH1D>((TH1D*)hists_.front()->Clone());

      for (unsigned idx=1; idx<hists_.size(); idx++)
        ahist->Add(hists_.at(idx));

      return std::move(ahist);
    }
  };

  TFile* datafile1 = new TFile("FakeEleAnalyzer_20UL16APV_data.root","READ");
  TFile* datafile2 = new TFile("FakeEleAnalyzer_20UL16_data.root","READ");
  TFile* datafile3 = new TFile("FakeEleAnalyzer_20UL17_data.root","READ");
  TFile* datafile4 = new TFile("FakeEleAnalyzer_20UL18_data.root","READ");

  auto retrieveDataHist = [] (TFile* afile, const TString& anlyzr, const TString& aname, const TString& branch, double xmin, double xmax, unsigned nbin=100, int rebin=1, TString multiply="", TF1* fit=nullptr, TFitResultPtr result=TFitResultPtr(-1), int updn=0) {
    TTree* atree = (TTree*)afile->Get(anlyzr+"/"+aname);
    TH1D* ahist = new TH1D( aname + branch + TString(afile->GetName())+TString(std::to_string(updn).c_str()), aname+branch, nbin, xmin, xmax );

    fillHist(atree,ahist,branch,multiply,fit,result,updn);

    ahist->SetLineColor(kBlack);
    ahist->Rebin(rebin);

    delete atree;

    return ahist;
  };

  auto sample1 = samples("20UL16APV",19.5);
  auto sample2 = samples("20UL16",16.8);
  auto sample3 = samples("20UL17",41.48);
  auto sample4 = samples("20UL18",59.83);

  struct histo {
    TString name_;
    const double xmin_;
    const double xmax_;
    const unsigned nbin_;

    histo(TString name, double xmin, double xmax, unsigned nbin)
    : name_(name), xmin_(xmin), xmax_(xmax), nbin_(nbin) {}
  };

  auto subtractHist = [] (const TH1D* ahist, const TH1D* bhist) -> TH1D* {
    TH1D* result = (TH1D*)ahist->Clone();

    for (unsigned ibin = 0; ibin < ahist->GetNbinsX()+2; ibin++) {
      if (ahist->GetBinContent(ibin)==0.)
        continue;

      result->SetBinContent(ibin, std::max( ahist->GetBinContent(ibin) - bhist->GetBinContent(ibin), 0.) );
      result->SetBinError(ibin, result->GetBinContent(ibin)==0. ? 0. : std::hypot(ahist->GetBinError(ibin), bhist->GetBinError(ibin)) );
    }

    return result;
  };

  std::vector<double> binningPt = {0.,40.,50.,60.,70.,80.,100.,120.,300.};

  auto subtractedHist = [&] (const TString& atree, histo& ahist, bool doRebin=true, TString multiply="", TF1* fit=nullptr, TFitResultPtr result=TFitResultPtr(-1), int updn=0) -> TH1D* {
    TH1D* datahist = retrieveDataHist(datafile1,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn);
    datahist->Add( retrieveDataHist(datafile2,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn) );
    datahist->Add( retrieveDataHist(datafile3,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn) );
    datahist->Add( retrieveDataHist(datafile4,"resolvedFakeEleCRanalyzerData",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn) );
    sample1.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn);
    sample2.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn);
    sample3.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn);
    sample4.retrieveHists("resolvedFakeEleCRanalyzer",atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,1,multiply,fit,result,updn);
    sample1.addHists(sample2);
    sample1.addHists(sample3);
    sample1.addHists(sample4);

    std::vector<double> binning;

    if (ahist.name_=="e1e2Pt")
      binning = binningPt;
    else if (ahist.name_=="e1e2InvM")
      binning = {0.0,0.6,1.,1.4,2.,5.,10.};
    else if (ahist.name_=="e1e2DR")
      binning = {0.0,0.04,0.06,0.08,0.1,0.3};

    TH1D* subtracted = datahist;

    if (doRebin) {
      const int nbin = binning.size()-1;
      TH1D* rebin = (TH1D*)datahist->Rebin(nbin,TString(datahist->GetName())+"_rebin",&(binning[0]));
      delete datahist;

      for (unsigned idx=0; idx<sample1.hists_.size(); ++idx) {
        auto* temp = sample1.hists_.at(idx);
        sample1.hists_.at(idx) = (TH1D*)sample1.hists_.at(idx)->Rebin(nbin,TString(sample1.hists_.at(idx)->GetName())+"_rebin",&(binning[0]));
        delete temp;
      }

      subtracted = rebin;
    }

    std::vector<unsigned> range;

    if (atree.Contains("1m"))
      range = {0,1,2,3,4,5,6,7,8,11,12,13};
    else if (atree.Contains("2m"))
      range = {0,1};

    for (unsigned idx : range) {
      auto* temp = subtracted;
      subtracted = subtractHist(subtracted,sample1.hists_.at(idx));
      delete temp;
    }

    return subtracted;
  };

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

  std::vector<histo> varNames = {histo("e1e2Pt",0.,300.,300),}; // histo("e1e2InvM",0.,20.,200) histo("e1e2DR",0.,0.3,300)

  TH1D* selected = nullptr;
  TF1* func = nullptr;
  TFitResultPtr fitResult;

  canvas2->cd();

  for (auto ahist : varNames) {
    auto* numer = subtractedHist("numerTree1m",ahist);
    auto* denom = subtractedHist("denomTree1m",ahist);

    numer->Divide(denom);
    numer->SetMinimum(0.);
    numer->Draw("E1");

    if (ahist.name_=="e1e2Pt") {
      numer->GetXaxis()->SetTitle("p_{T}(ee) [GeV]");
      selected = numer;

      func = new TF1("REFFdr03","[0]+[1]*TMath::Exp(-[2]*x)",40.,300.);
      func->SetParameter(2,0.0001);
      func->SetLineColor(kRed);
      func->SetLineWidth(2);
      func->SetLineStyle(2);
      func->SetParLimits(0,0.,1.);
      fitResult = selected->Fit(func,"RS");

      canvas2->Update();
      TPaveStats* stats = (TPaveStats*)canvas2->GetPrimitive("stats");
      stats->SetTextColor(kRed);
      stats->SetX1NDC(.65);
      stats->SetX2NDC(.95);
      stats->SetY1NDC(.7);
      stats->SetY2NDC(.9);

      const unsigned nbinSel = selected->GetNbinsX();
      std::vector<double> xcen = estimateCenter(binningPt);
      std::vector<double> xbinw = estimateWidth(binningPt);
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
    }

    canvas2->Update();

    CMS_lumi( canvas2.get(), iPeriod, iPos );

    canvas2->SaveAs("REFFdr03"+ahist.name_+".pdf");

    delete denom;
  }

  std::cout << "write " << "REFFdr03_"+era+".root" << std::endl;

  TFile* outfile = new TFile("REFFdr03_"+era+".root","RECREATE");
  TH1D* output = (TH1D*)selected->Clone();
  output->SetName("REFFhist");
  output->Write();
  ((TF1*)func->Clone())->Write();
  outfile->Close();

  std::vector<histo> validateNames = {
    histo("e1scEta",-3.,3.,100),
    histo("e2scEta",-3.,3.,100),
    histo("e1e2Pt",0.,200.,100),
    histo("e1Pt",0.,200.,100),
    histo("e2Pt",0.,200.,100),
    histo("e1e2InvM",0.,20.,100),
    histo("e1e2DR",0.,0.4,100)
  };

  auto canvas = std::make_unique<TCanvas>("canvas","canvas",800,800);

  canvas->SetFillColor(0);
  canvas->SetBorderMode(0);
  canvas->SetFrameFillStyle(0);
  canvas->SetFrameBorderMode(0);
  canvas->SetLeftMargin( L/W );
  canvas->SetRightMargin( R/W );
  canvas->SetTopMargin( T/H );
  canvas->SetBottomMargin( B/H );
  canvas->SetTickx(0);
  canvas->SetTicky(0);

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
  p2->SetTopMargin(0.0);
  p2->SetBottomMargin(0.3);
  p2->SetLeftMargin( L/W );
  p2->SetRightMargin( R/W );
  p2->SetTickx(0);
  p2->SetTicky(0);
  p2->Draw();

  for (auto ahist : validateNames) {
    for (TString atree : {"2m","1m"}) {
      auto* denom = subtractedHist("denomTree"+atree,ahist,false,"e1e2Pt",func,fitResult,0);
      denom->SetFillColor(TColor::GetColor("#92dadd"));
      denom->SetLineColor(TColor::GetColor("#92dadd"));
      denom->Rebin(2);

      auto* denomUp = subtractedHist("denomTree"+atree,ahist,false,"e1e2Pt",func,fitResult,1);
      auto* denomDn = subtractedHist("denomTree"+atree,ahist,false,"e1e2Pt",func,fitResult,-1);
      denomUp->Rebin(2);
      denomDn->Rebin(2);

      auto* datahist = retrieveDataHist(datafile1,"resolvedFakeEleCRanalyzerData","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2);
      datahist->Add( retrieveDataHist(datafile2,"resolvedFakeEleCRanalyzerData","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2) );
      datahist->Add( retrieveDataHist(datafile3,"resolvedFakeEleCRanalyzerData","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2) );
      datahist->Add( retrieveDataHist(datafile4,"resolvedFakeEleCRanalyzerData","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2) );
      sample1.retrieveHists("resolvedFakeEleCRanalyzer","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2);
      sample2.retrieveHists("resolvedFakeEleCRanalyzer","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2);
      sample3.retrieveHists("resolvedFakeEleCRanalyzer","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2);
      sample4.retrieveHists("resolvedFakeEleCRanalyzer","numerTree"+atree,ahist.name_,ahist.xmin_,ahist.xmax_,ahist.nbin_,2);
      sample1.addHists(sample2);
      sample1.addHists(sample3);
      sample1.addHists(sample4);

      auto astack = std::make_unique<THStack>("stack","stack");

      std::vector<unsigned> range;

      if (atree=="1m")
        range = {0,1,2,3,4,5,6,7,8,11,12,13};
      else if (atree=="2m")
        range = {0,1};

      TH1D* added = (TH1D*)denom->Clone();

      for (unsigned idx : range) {
        astack->Add(sample1.hists_.at(idx));
        added->Add(sample1.hists_.at(idx));
      }

      astack->Add(denom);

      std::vector<double> x0, y0, errx, erryDn, erryUp;
      std::vector<double> r0, errRup, errRdn;

      for (int idx=1; idx<=added->GetNbinsX(); idx++) {
        x0.push_back(added->GetBinCenter(idx));
        y0.push_back(added->GetBinContent(idx));
        errx.push_back(added->GetBinWidth(idx)/2.);
        erryUp.push_back( denomUp->GetBinContent(idx)-denom->GetBinContent(idx) );
        erryDn.push_back( denom->GetBinContent(idx)-denomDn->GetBinContent(idx) );
      }

      datahist->SetLineColor(kBlack);
      datahist->SetLineWidth(2);

      datahist->SetMaximum( 1.5*datahist->GetMaximum() );
      datahist->SetMinimum( 0.01 );

      p1->cd();

      datahist->Draw("E1");

      astack->Draw("hist&same");
      datahist->Draw("E1&same");

      auto gr = new TGraphAsymmErrors(added->GetNbinsX(),&(x0[0]),&(y0[0]),&(errx[0]),&(errx[0]),&(erryDn[0]),&(erryUp[0]));
      gr->SetFillColor(kGray+3);
      gr->SetLineColor(kGray+3);
      gr->SetFillStyle(3004);
      gr->Draw("2");

      auto legend = std::make_unique<TLegend>(0.6,0.55,0.95,0.9);
      legend->SetBorderSize(0);
      legend->SetNColumns(2);
      legend->AddEntry(datahist,"Data");

      if (atree=="1m") {
        legend->AddEntry(sample1.hists_.at(11),"DY");
        legend->AddEntry(sample1.hists_.at(8),"TTlep");
        legend->AddEntry(sample1.hists_.at(7),"ST_s-ch");
        legend->AddEntry(sample1.hists_.at(5),"ST_t-ch");
        legend->AddEntry(sample1.hists_.at(3),"ST_tW");
        legend->AddEntry(sample1.hists_.at(2),"WW");
      }

      legend->AddEntry(sample1.hists_.at(1),"WZ");
      legend->AddEntry(sample1.hists_.at(0),"ZZ");
      legend->AddEntry(denom,"nonprompt");
      legend->Draw();

      p2->cd();

      TH1* ratio = (TH1*)datahist->Clone("ratio");
      ratio->SetStats(0);
      ratio->SetTitle("");
      ratio->Divide(added);
      ratio->GetYaxis()->SetTitle("Obs/Exp");
      ratio->GetYaxis()->SetTitleSize(0.1);
      ratio->GetYaxis()->SetTitleOffset(0.4);
      ratio->GetXaxis()->SetLabelSize(0.1);
      ratio->GetYaxis()->SetLabelSize(0.1);
      ratio->GetXaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetRangeUser(0.01,1.99);
      ratio->GetYaxis()->SetNdivisions(505);
      ratio->GetXaxis()->SetTitleSize(0.12);
      ratio->GetXaxis()->SetTitleOffset(1.0);
      ratio->SetLineColor(kBlack);

      if (ahist.name_=="e1e2Pt")
        ratio->GetXaxis()->SetTitle("p_{T}(ee) [GeV]");
      else if (ahist.name_=="e1e2DR")
        ratio->GetXaxis()->SetTitle("#Delta R(ee)");

      ratio->Draw("E1");

      for (int idx=1; idx<=added->GetNbinsX(); idx++) {
        r0.push_back(1.);
        errRup.push_back( added->GetBinContent(idx)!=0. ? erryUp.at(idx-1)/added->GetBinContent(idx) : 0. );
        errRdn.push_back( added->GetBinContent(idx)!=0. ? erryDn.at(idx-1)/added->GetBinContent(idx) : 0. );
      }

      auto rgr = new TGraphAsymmErrors(added->GetNbinsX(),&(x0[0]),&(r0[0]),&(errx[0]),&(errx[0]),&(errRdn[0]),&(errRup[0]));
      rgr->SetFillColor(kGray+2);
      rgr->SetLineColor(kGray+2);
      rgr->SetFillStyle(3004);
      rgr->Draw("2");

      canvas->Update();

      CMS_lumi( canvas.get(), iPeriod, iPos );

      p1->RedrawAxis();
      p1->GetFrame()->Draw();

      canvas->SaveAs(atree+"_"+ahist.name_+".png");

      if (ahist.name_=="e1e2Pt" || ahist.name_=="e1e2DR")  
        SaveAs(canvas.get(),TString(atree+"_"+ahist.name_+".pdf").Data(),p1); 
    }
  }
}
