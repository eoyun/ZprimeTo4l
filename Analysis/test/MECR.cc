#include "TROOT.h"
#include <cstdarg>
#include <cmath>
#include <algorithm>

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void MECR(TString era) {
  TH1::SetDefaultSumw2();

  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  double valLumi = 0.;
  TString postfix = era;

  if (era=="20UL16APV") {
    lumi_13TeV = "19.5 fb^{-1}";
    valLumi = 19.5;
  } else if (era=="20UL16") {
    lumi_13TeV = "16.8 fb^{-1}";
    valLumi = 16.8;
    postfix = "";
  } else if (era=="20UL17") {
    lumi_13TeV = "41.48 fb^{-1}";
    valLumi = 41.48;
  } else if (era=="20UL18") {
    lumi_13TeV = "59.83 fb^{-1}";
    valLumi = 59.83;
  } else {
    std::cout << "check params!" << std::endl;
    throw;
  }


  int iPeriod = 4;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
  int iPos = 0;

  if( iPos==0 )
    relPosX = 0.12;

  int W = 800;
  int H = 600;

  int H_ref = 600;
  int W_ref = 800;

  // references for T, B, L, R
  float T = 0.08*H_ref;
  float B = 0.12*H_ref;
  float L = 0.12*W_ref;
  float R = 0.04*W_ref;

  class sample {
  public:
    TFile* file;
    double xsec;
    int color;

    sample(TFile* afile, const double axsec, const int acolor) {
      file = afile;
      xsec = axsec;
      color = acolor;
    }
  };

  std::vector<sample> bkglist = {
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT50to100.root"),187700000.0,kCyan-10),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT100to200.root"),23640000.0,kCyan-9),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT200to300.root"),1546000.0,kCyan-7),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT300to500.root"),321600.0,kCyan-6),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT500to700.root"),30980.0,kCyan-3),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT700to1000.root"),6364.0,kCyan-2),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT1000to1500.root"),1117.0,kCyan+2),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT1500to2000.root"),108.4,kCyan+3),
    sample(TFile::Open("MergedEleCR_"+era+"_QCD_HT2000toInf.root"),22.36,kCyan+4),
    sample(TFile::Open("MergedEleCR_"+era+"_WJets_0J.root"),53330.0,kGreen-3),
    sample(TFile::Open("MergedEleCR_"+era+"_WJets_1J.root"),8875.0,kGreen+2),
    sample(TFile::Open("MergedEleCR_"+era+"_WJets_2J.root"),3338.0,kGreen-1),
    sample(TFile::Open("MergedEleCR_"+era+"_DY_0J.root"),5090.0*6077.22/6424.0,kOrange), // reweight to NNLO xsec
    sample(TFile::Open("MergedEleCR_"+era+"_DY_1J.root"),983.5*6077.22/6424.0,kOrange+2),
    sample(TFile::Open("MergedEleCR_"+era+"_DY_2J.root"),353.6*6077.22/6424.0,kOrange+3),
    sample(TFile::Open("MergedEleCR_"+era+"_TT.root"),831.76,kRed+1), // FXFX inclusive
    sample(TFile::Open("MergedEleCR_"+era+"_WZ.root"),5.213,kViolet+1),
    sample(TFile::Open("MergedEleCR_"+era+"_ZZ.root"),1.325,kBlue+1),
    sample(TFile::Open("MergedEleCR_"+era+"_ZG.root"),51.1,kAzure-4)
  };

  auto retrieveHists = [](const TString aname, std::vector<sample> list) {
    std::vector<TH1D*> result;

    for (const auto& element : list)
      result.push_back((TH1D*)(element.file)->Get(aname));

    return result;
  };

  TFile* file_data = TFile::Open("MergedEleCR_"+era+"_data.root");

  auto setAlphabetLabel = [](TH1* h, const char* name[]) {
    for (int idx = 1; idx <= h->GetNbinsX(); idx++)
      h->GetXaxis()->SetBinLabel(idx,name[idx-1]);

    h->GetXaxis()->SetLabelSize(0.15);
  };

  auto GetHists = [&](TFile* file, bool isMC=true) -> std::vector<TH1D*> {
    std::vector<TH1D*> result;

    if (!isMC)
      postfix = "Data";

    TDirectory* dir = dynamic_cast<TDirectory*>(file->Get("mergedEleCRanalyzer"+postfix));

    for ( auto&& keyobj : *dir->GetListOfKeys() ) {
      auto key = (TKey*) keyobj;

      if ( TString(key->GetName()).Contains("totWeightedSum") )
        continue;

      if ( TString(key->GetName()).Contains("FF") )
        continue;

      if ( TString(key->GetName()).Contains("mva") )
        continue;

      if (key->GetClassName()==TString("TH1D")) {
        result.push_back( static_cast<TH1D*>(dir->Get(key->GetName())) );

        if (TString(key->GetName()).Contains("3E_antiME_lll_invM_CR"))
          result.back()->Rebin(2);

        if (TString(key->GetName()).Contains("2E_Et") && TString(key->GetName()).Contains("mixed"))
          result.back()->Rebin(2);

        if (TString(key->GetName()).Contains("2E_antiME") && TString(key->GetName()).Contains("invM") && TString(key->GetName()).Contains("CR"))
          result.back()->GetXaxis()->SetRangeUser(50.,200.);

        if (TString(key->GetName()).Contains("2E_mixedME") && TString(key->GetName()).Contains("invM"))
          result.back()->GetXaxis()->SetRangeUser(50.,200.);
      }
    }

    return std::move(result);
  };

  auto scaling = [&](std::vector<TH1D*>& hists, int color, float xsec=0., double sumwgt=0., bool isSig=false) {
    for (unsigned idx = 0; idx < hists.size(); idx++) {
      hists.at(idx)->Scale(valLumi*1000.*xsec/sumwgt); // 20UL16 16.8 // 14.964 for v8 v9

      if ( !isSig )
        hists.at(idx)->SetFillColor(color);

      hists.at(idx)->SetLineColor(color);
    }
  };

  auto stacking = [&](std::vector<std::unique_ptr<THStack>>& stacks, std::initializer_list<std::vector<TH1D*>> aBkgHists) {
    for (unsigned idx = 0; idx < stacks.size(); idx++) {
      for (auto& element : aBkgHists)
        stacks.at(idx)->Add(element.at(idx));
    }
  };

  auto sumwgt_bkgs = retrieveHists("mergedEleCRanalyzer"+postfix+"/totWeightedSum", bkglist);

  std::vector<std::vector<TH1D*>> histowner;
  histowner.reserve(bkglist.size());

  for (unsigned idx = 0; idx < bkglist.size(); idx++) {
    std::string numstr = std::to_string(idx);
    auto hists = GetHists(bkglist.at(idx).file);
    histowner.emplace_back(std::move(hists));
  }

  auto dataHists_0 = GetHists(file_data,false);

  for (unsigned idx = 0; idx < bkglist.size(); idx++)
    scaling(histowner.at(idx),
            bkglist.at(idx).color,
            bkglist.at(idx).xsec,
            sumwgt_bkgs.at(idx)->GetBinContent(1));

  std::vector<std::unique_ptr<THStack>> stacks;

  for (unsigned idx = 0; idx < histowner.front().size(); idx++)
    stacks.emplace_back(std::make_unique<THStack>(TString(histowner.front().at(idx)->GetName())+"_stack",histowner.front().at(idx)->GetTitle()));

  stacking(stacks,{
    histowner.at(17), // zz
    histowner.at(16), // zz
    histowner.at(18), // zz'
    histowner.at(15), // zz

    histowner.at(14),
    histowner.at(13),
    histowner.at(12), // dy

    histowner.at(11),
    histowner.at(10),
    histowner.at(9), // w

    histowner.at(8),
    histowner.at(7),
    histowner.at(6),
    histowner.at(5),
    histowner.at(4),
    histowner.at(3),
    histowner.at(2),
    histowner.at(1),
    histowner.at(0)
  });

  auto draw = [&](std::vector<TH1D*>& dataHists,
                 std::vector<std::unique_ptr<THStack>>& bkgHists) {
    auto canvas = std::make_unique<TCanvas>("canvas","canvas",800,600);

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

    auto legend = std::make_unique<TLegend>(0.7,0.75,0.9,0.9);
    legend->SetNColumns(2);
    legend->SetBorderSize(0);
    legend->AddEntry(dataHists.at(0),era);
    legend->AddEntry(histowner.at(4).at(0),"QCD");
    legend->AddEntry(histowner.at(10).at(0),"WJets");
    legend->AddEntry(histowner.at(12).at(0),"DY");
    legend->AddEntry(histowner.at(15).at(0),"TT");
    legend->AddEntry(histowner.at(16).at(0),"WZ");
    legend->AddEntry(histowner.at(17).at(0),"ZZ");
    legend->AddEntry(histowner.at(18).at(0),"ZG");

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
    p2->SetTopMargin(0.03);
    p2->SetBottomMargin(0.3);
    p2->SetLeftMargin( L/W );
    p2->SetRightMargin( R/W );
    p2->SetTickx(0);
    p2->SetTicky(0);
    p2->Draw();

    for (unsigned idx = 0; idx < bkgHists.size(); idx++) {
      p1->cd();
      p1->SetLogy();
      double theMax = dataHists.at(idx)->GetMaximum();

      dataHists.at(idx)->SetLineColor(kBlack);
      dataHists.at(idx)->SetLineWidth(2);

      if (TString(dataHists.at(idx)->GetName()).Contains("N-1")) {
        p1->SetLogy(0);
        dataHists.at(idx)->SetMaximum(1.5*theMax);
      } else {
        dataHists.at(idx)->SetMaximum(10.*theMax);
        dataHists.at(idx)->SetMinimum( 0.2 );
      }

      dataHists.at(idx)->Draw("E1");

      bkgHists.at(idx)->Draw("hist&same");
      theMax = std::max(bkgHists.at(idx)->GetMaximum(),theMax);
      dataHists.at(idx)->Draw("E1&same");
      legend->Draw();

      p2->cd();

      TList* stackHists = bkgHists.at(idx)->GetHists();

      TH1* tmpHist = (TH1*)stackHists->At(0)->Clone();

      for (int idx = 1; idx < stackHists->GetSize(); ++idx)
        tmpHist->Add((TH1*)stackHists->At(idx));

      TH1* ratio = (TH1*)dataHists.at(idx)->Clone("ratio");
      ratio->SetStats(0);
      ratio->SetTitle("");
      ratio->Divide(tmpHist);
      ratio->GetYaxis()->SetTitle("Data/MC");
      ratio->GetYaxis()->SetTitleSize(0.1);
      ratio->GetYaxis()->SetTitleOffset(0.4);
      ratio->GetXaxis()->SetLabelSize(0.1);
      ratio->GetYaxis()->SetLabelSize(0.1);
      ratio->GetXaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetLabelOffset(0.01);
      ratio->GetYaxis()->SetRangeUser(0.5,1.5);
      ratio->GetXaxis()->SetTitleSize(0.12);
      ratio->GetXaxis()->SetTitleOffset(0.75);
      ratio->SetLineColor(kBlack);
      // ratio->SetLineWidth(5);
      // ratio->Draw("Hist");

      if (TString(dataHists.at(idx)->GetName()).Contains("N-1")) {
        const char* nmLabel[] = {"","min E_{T}","#eta_{SC}","#Delta#eta_{in}(seed)","#Delta#phi_{in}(SC)","full5x5 #sigma_{i#eta i#eta}","E2x5/E5x5","H/E","Trk iso","Ecal+HcalD1 iso","dxy","missing hits","is Ecal driven","","",""};
        setAlphabetLabel(ratio,nmLabel);
      }

      ratio->Draw("E1");

      // TString histname = bkgHists.at(idx)->GetName();
      // if ( histname.Contains("pt") || histname.Contains("eta") || histname.Contains("phi") ) {
      //   canvas->SetLogy(0);
      //   dataHists.at(idx)->SetMaximum(1.5*theMax);
      //   dataHists.at(idx)->SetMinimum(0.);
      // }

      canvas->Update();

      // writing the lumi information and the CMS "logo"
      CMS_lumi( canvas.get(), iPeriod, iPos );

      p1->RedrawAxis();
      p1->GetFrame()->Draw();

      canvas->SaveAs(TString("")+bkgHists.at(idx)->GetName()+".png");

      delete tmpHist;
      delete ratio;
    } // histograms
  };

  draw(dataHists_0,stacks);

  return;
}
