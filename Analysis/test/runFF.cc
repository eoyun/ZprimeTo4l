#include "TFile.h"
#include "TH1D.h"
#include "THStack.h"
#include "TCanvas.h"

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void runFF() {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)
  lumi_13TeV = "59.83 fb^{-1}";

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

  TFile* datafile = new TFile("MergedEleCR_20UL18_data.root","READ");
  TFile* WZfile = new TFile("MergedEleCR_20UL18_WZ.root","READ");
  TFile* ZZfile = new TFile("MergedEleCR_20UL18_ZZ.root","READ");
  TFile* FFfile = new TFile("MEFF_20UL18.root","READ");

  TF1* sslow = static_cast<TF1*>(FFfile->FindObjectAny("sslow"));
  TF1* sshigh = static_cast<TF1*>(FFfile->FindObjectAny("sshigh"));
  TF1* oslow = static_cast<TF1*>(FFfile->FindObjectAny("oslow"));
  TF1* oshigh = static_cast<TF1*>(FFfile->FindObjectAny("oshigh"));

  const double ssBoundary = 200.;
  const double osBoundary = 500.;

  TH1D* SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_CRME")->Clone();
  TH1D* SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_mixedAntiME")->Clone();
  TH1D* SSnumAnti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_mixedME")->Clone();
  TH1D* SSanti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_SSCR_antiME")->Clone();

  TH1D* OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_CRME")->Clone();
  TH1D* OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_mixedAntiME")->Clone();
  TH1D* OSnumAnti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_mixedME")->Clone();
  TH1D* OSanti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_Et_OSCR_antiME")->Clone();

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

  auto compareEt = [] (const TF1* low,
                       const TF1* high,
                       TH1D* num,
                       TH1D* denom,
                       const double boundary,
                       const int color,
                       TPad* pad) -> TH1D* {
    auto* cloned = (TH1D*)denom->Clone();
    for (unsigned idx = 0; idx < cloned->GetNbinsX()+2; idx++) {
      const double con = cloned->GetBinContent(idx);

      if ( con==0. )
        continue;

      const double cen = cloned->GetBinCenter(idx);
      const double ff = (cen > boundary) ? high->Eval(cen) : low->Eval(cen);

      cloned->SetBinContent(idx,ff*con);
      cloned->SetBinError(idx,ff*denom->GetBinError(idx));
    }

    std::cout << num->Integral() <<" "<< cloned->Integral()<< std::endl;

    cloned->Rebin(2);
    num->Rebin(2);

    cloned->SetFillColor(color);
    num->SetLineWidth(2);
    num->SetLineColor(kBlack);
    num->GetXaxis()->SetRangeUser(0.,500.);
    pad->cd();
    num->Draw("E1");
    cloned->SetLineWidth(0);
    cloned->Draw("hist&same");
    num->Draw("E1&same");

    return cloned;
  };

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

  compareEt(sslow,sshigh,SSnum,SSdenom,200.,kCyan+1,canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedSS.png");
  compareEt(oslow,oshigh,OSnum,OSdenom,500.,kOrange,canvas_1);
  SaveAs(canvas_1,"FF_2E_Et_denom_mixedOS.png");
  auto* denom_EtSS = compareEt(sslow,sshigh,SSnumAnti,SSanti,200.,kCyan+1,p1);
  drawRatio(SSnumAnti,denom_EtSS,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiSS.png",p1);
  auto* denom_EtOS = compareEt(oslow,oshigh,OSnumAnti,OSanti,500.,kOrange,p1);
  drawRatio(OSnumAnti,denom_EtOS,p2);
  SaveAs(canvas_2,"FF_2E_Et_denom_antiOS.png",p1);

  TH1D* invM_SSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_CRME_SSll_invM")->Clone();
  TH1D* invM_SSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_SSll_invM_xFF")->Clone();
  TH1D* invM_SSmixed = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_SSll_invM")->Clone();
  TH1D* invM_SSanti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_SSll_invM_CR_xFF")->Clone();

  TH1D* invM_OSnum = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_CRME_OSll_invM")->Clone();
  TH1D* invM_OSdenom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_OSll_invM_xFF")->Clone();
  TH1D* invM_OSmixed = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_mixedME_OSll_invM")->Clone();
  TH1D* invM_OSanti = (TH1D*)datafile->Get("mergedEleCRanalyzerData/2E_antiME_OSll_invM_CR_xFF")->Clone();

  std::cout << invM_OSnum->Integral() << " and " << invM_OSdenom->Integral() << std::endl;

  auto compare = [] (TH1D* num, TH1D* denom, const int color, TPad* pad) {
    if ( denom->GetNbinsX() % num->GetNbinsX()!=0 ) {
      // hardcode (rebin to GCD)
      denom->Rebin( 5 );
      num->Rebin( 2 );
    } else {
      denom->Rebin( denom->GetNbinsX()/num->GetNbinsX() );
    }

    denom->SetFillColor(color);
    pad->cd();
    num->SetMaximum( 1.2*std::max(num->GetMaximum(),denom->GetMaximum()) );
    num->SetLineWidth(2);
    num->SetLineColor(kBlack);
    num->Draw("E1");
    denom->SetLineWidth(0);
    denom->Draw("hist&same");
    num->Draw("E1&same");
  };

  compare(invM_SSnum,invM_SSdenom,kCyan+1,canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedSS.png");
  compare(invM_SSmixed,invM_SSanti,kCyan+1,p1);
  drawRatio(invM_SSmixed,invM_SSanti,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiSS.png",p1);
  compare(invM_OSnum,invM_OSdenom,kOrange,canvas_1);
  SaveAs(canvas_1,"FF_2E_invM_denom_mixedOS.png");
  compare(invM_OSmixed,invM_OSanti,kOrange,p1);
  drawRatio(invM_OSmixed,invM_OSanti,p2);
  SaveAs(canvas_2,"FF_2E_invM_denom_antiOS.png",p1);

  TH1D* invM_3E_denomSS = (TH1D*)datafile->Get("mergedEleCRanalyzerData/3E_antiME_lll_invM_CR_xSSFF")->Clone();
  TH1D* invM_3E_num = (TH1D*)datafile->Get("mergedEleCRanalyzerData/3E_CRME_lll_invM")->Clone();
  TH1D* Et_3E_denom = (TH1D*)datafile->Get("mergedEleCRanalyzerData/3E_Et_CR_antiME")->Clone();
  TH1D* Et_3E_num = (TH1D*)datafile->Get("mergedEleCRanalyzerData/3E_Et_CR_CRME")->Clone();

  TString MCanalyzerName = "mergedEleCRanalyzer20UL18";
  double intlumi = 59.83;
  double WZxsec = 5.213;
  double ZZxsec = 1.325;
  double WZsumwgt = ((TH1D*)WZfile->Get(MCanalyzerName+"/totWeightedSum"))->GetBinContent(1);
  double ZZsumwgt = ((TH1D*)ZZfile->Get(MCanalyzerName+"/totWeightedSum"))->GetBinContent(1);

  TH1D* invM_3E_denomSS_WZ = (TH1D*)WZfile->Get(MCanalyzerName+"/3E_antiME_lll_invM_CR_xSSFF")->Clone();
  invM_3E_denomSS_WZ->Scale( intlumi*1000.*WZxsec/WZsumwgt );
  TH1D* invM_3E_denomSS_ZZ = (TH1D*)ZZfile->Get(MCanalyzerName+"/3E_antiME_lll_invM_CR_xSSFF")->Clone();
  invM_3E_denomSS_ZZ->Scale( intlumi*1000.*ZZxsec/ZZsumwgt );
  TH1D* invM_3E_denomOS_WZ = (TH1D*)WZfile->Get(MCanalyzerName+"/3E_antiME_lll_invM_CR_xOSFF")->Clone();
  invM_3E_denomOS_WZ->Scale( intlumi*1000.*WZxsec/WZsumwgt );
  TH1D* invM_3E_denomOS_ZZ = (TH1D*)ZZfile->Get(MCanalyzerName+"/3E_antiME_lll_invM_CR_xOSFF")->Clone();
  invM_3E_denomOS_ZZ->Scale( intlumi*1000.*ZZxsec/ZZsumwgt );
  TH1D* Et_3E_denom_WZ = (TH1D*)WZfile->Get(MCanalyzerName+"/3E_Et_CR_antiME")->Clone();
  Et_3E_denom_WZ->Scale( intlumi*1000.*WZxsec/WZsumwgt );
  TH1D* Et_3E_denom_ZZ = (TH1D*)ZZfile->Get(MCanalyzerName+"/3E_Et_CR_antiME")->Clone();
  Et_3E_denom_ZZ->Scale( intlumi*1000.*ZZxsec/ZZsumwgt );

  auto square = [](double x) { return x*x; };

  auto subtractHist = [&square](const TH1D* denom, const TH1D* denom_prompt) -> TH1D* {
    auto* cloned = (TH1D*)denom->Clone();

    for (unsigned idx = 0; idx < cloned->GetNbinsX()+2; idx++) {
      cloned->SetBinContent(idx,cloned->GetBinContent(idx)-denom_prompt->GetBinContent(idx));
      cloned->SetBinError(idx,std::sqrt(square(cloned->GetBinError(idx))+square(denom_prompt->GetBinError(idx))));
    }

    return cloned;
  };

  canvas_2->cd();

  TH1D* invM_3E_denomSSfinal = subtractHist( subtractHist(invM_3E_denomSS,invM_3E_denomSS_WZ), invM_3E_denomSS_ZZ);
  invM_3E_denomSSfinal->SetFillColor(kCyan+1);
  invM_3E_denomSSfinal->Rebin( invM_3E_denomSSfinal->GetNbinsX()/invM_3E_num->GetNbinsX() );
  TH1D* invM_3E_denomOSfinal = (TH1D*)invM_3E_denomOS_WZ->Clone();
  invM_3E_denomOSfinal->Add(invM_3E_denomOS_ZZ);
  invM_3E_denomOSfinal->SetFillColor(kOrange);
  invM_3E_denomOSfinal->Rebin( invM_3E_denomOSfinal->GetNbinsX()/invM_3E_num->GetNbinsX() );

  THStack* invM_3E_denomFinal = new THStack("invM_3E","invM_3E;GeV;");
  invM_3E_denomOSfinal->SetLineWidth(0);
  invM_3E_denomSSfinal->SetLineWidth(0);
  invM_3E_denomFinal->Add(invM_3E_denomOSfinal);
  invM_3E_denomFinal->Add(invM_3E_denomSSfinal);

  invM_3E_num->SetLineWidth(2);
  invM_3E_num->SetLineColor(kBlack);

  invM_3E_num->Draw("E1");
  invM_3E_denomFinal->Draw("hist&same");
  invM_3E_num->Draw("E1&same");
  std::cout<<invM_3E_num->Integral()<< " " <<invM_3E_denomSSfinal->Integral()<<std::endl;

  TLegend* legend_left = new TLegend(0.15,0.7,0.4,0.9);
  legend_left->SetBorderSize(0);
  legend_left->AddEntry(invM_3E_num,"Data");
  legend_left->AddEntry(invM_3E_denomSSfinal,"Data-driven (X #rightarrow ME)");
  legend_left->AddEntry(invM_3E_denomOSfinal,"Data-driven (e #rightarrow ME)");
  legend_left->Draw();

  SaveAs(canvas_2,"FF_3E_invM.png");

  TH1D* Et_3E_denom_OS = (TH1D*)Et_3E_denom_WZ->Clone();
  Et_3E_denom_OS->Add(Et_3E_denom_ZZ);

  for (unsigned idx = 0; idx < Et_3E_denom->GetNbinsX()+2; idx++) {
    const double cen = Et_3E_denom->GetBinCenter(idx);
    const double ssff = (cen > ssBoundary) ? sshigh->Eval(cen) : sslow->Eval(cen);
    const double osff = (cen > osBoundary) ? oshigh->Eval(cen) : oslow->Eval(cen);
    Et_3E_denom->SetBinContent( idx, ssff*(Et_3E_denom->GetBinContent(idx) - Et_3E_denom_WZ->GetBinContent(idx) - Et_3E_denom_ZZ->GetBinContent(idx)) );
    Et_3E_denom->SetBinError( idx, ssff*std::sqrt( square(Et_3E_denom->GetBinError(idx)) + square(Et_3E_denom_WZ->GetBinError(idx)) + square(Et_3E_denom_ZZ->GetBinError(idx)) ) );

    Et_3E_denom_OS->SetBinContent( idx, osff*Et_3E_denom_OS->GetBinContent(idx) );
    Et_3E_denom_OS->SetBinError( idx, osff*Et_3E_denom_OS->GetBinError(idx) );
  }

  Et_3E_denom->Rebin(2);
  Et_3E_denom_OS->Rebin(2);
  Et_3E_num->Rebin(2);

  THStack* Et_3E_denomFinal = new THStack("Et_3E","Et_3E;E_{T} [GeV];");
  Et_3E_denom_OS->SetFillColor(kOrange);
  Et_3E_denom->SetFillColor(kCyan+1);
  Et_3E_denom_OS->SetLineWidth(0);
  Et_3E_denom->SetLineWidth(0);
  Et_3E_denomFinal->Add(Et_3E_denom_OS);
  Et_3E_denomFinal->Add(Et_3E_denom);

  Et_3E_num->SetLineWidth(2);
  Et_3E_num->SetLineColor(kBlack);
  Et_3E_num->GetXaxis()->SetRangeUser(0.,500.);
  Et_3E_num->SetMaximum(1.2*std::max(Et_3E_num->GetMaximum(),Et_3E_denomFinal->GetMaximum()));

  Et_3E_num->Draw("E1");
  Et_3E_denomFinal->Draw("hist&same");
  Et_3E_num->Draw("E1&same");

  std::cout<<Et_3E_num->Integral()<<" "<<Et_3E_denom->Integral()<<std::endl;

  TLegend* legend_right = new TLegend(0.68,0.7,0.93,0.9);
  legend_right->SetBorderSize(0);
  legend_right->AddEntry(Et_3E_num,"Data");
  legend_right->AddEntry(Et_3E_denom,"Data-driven (X #rightarrow ME)");
  legend_right->AddEntry(Et_3E_denom_OS,"Data-driven (e #rightarrow ME)");
  legend_right->Draw();

  SaveAs(canvas_2,"FF_3E_Et.png");

  return;
}
