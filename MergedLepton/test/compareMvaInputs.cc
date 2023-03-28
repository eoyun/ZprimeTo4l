#include "TROOT.h"
#include <cstdarg>
#include <cmath>
#include <initializer_list>

#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/tdrstyle.C"
#include "/u/user/sako/ModHEEP/CMSSW_10_6_29/work/CMS_lumi.C"

typedef struct {
  float weight;
  float pt, eta, phi, en, et;
  int charge;
  float enSC, etSC, etaSC, phiSC, etaSCWidth, phiSCWidth;
  float full5x5_sigmaIetaIeta, full5x5_sigmaIphiIphi;
  float full5x5_E1x5, full5x5_E2x5, full5x5_E5x5, full5x5_hOverE, full5x5_r9;
  float dEtaIn, dPhiIn, dPhiSeed, dEtaEle, dPhiEle, dEtaSeed;
  float EseedOverP, EOverP;
  float ecalEn, ecalErr, trkErr, combErr, PFcombErr;
  float dr03TkSumPtHEEP, dr03EcalRecHitSumEt, dr03HcalDepth1TowerSumEt;
  float modTrkIso, modEcalIso;
  int lostHits, nValidHits, nValidPixelHits, GsfHits;
  float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz;
  float Gsfpt, Gsfeta, Gsfphi;
  int expectedMissingInnerHits;
  float convVtxFitProb, convVtxChi2, convDist, convDcot, convRadius;
  int passConversionVeto, nbrem;
  float fbrem, fbremSC;
  float full5x5_e2x5Left, full5x5_e2x5Right, full5x5_e2x5Top, full5x5_e2x5Bottom, full5x5_eLeft, full5x5_eRight, full5x5_eTop, full5x5_eBottom;
  float full5x5_eMax, full5x5_e2nd;
} ElectronStruct;

typedef struct {
  float weight;
  float Gsfpt, Gsfeta, Gsfphi;
  int lostHits, nValidHits, nValidPixelHits;
  float chi2, d0, d0Err, dxyErr, vz, dzErr, dxy, dz, ptErr;
} AddGsfStruct;

constexpr double pi = 3.14159265358979323846;

double deltaPhi(double phi1, double phi2) {
  double dphi = phi2 - phi1;

  if (dphi > pi)
    dphi -= 2.*pi;
  else if (dphi <= -pi)
    dphi += 2.*pi;

  return dphi;
}

double deltaR2(double eta1, double phi1,double eta2,double phi2) {
  return std::abs(eta1-eta2)*std::abs(eta1-eta2) + std::abs(deltaPhi(phi1,phi2))*std::abs(deltaPhi(phi1,phi2));
}

void compareMvaInputs(const std::string ang, const std::string etThres, const std::string EBEE, const std::string era) {
  TH1::SetDefaultSumw2();

  setTDRStyle();

  writeExtraText = true;            // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";            // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;                  // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
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

  TFile* file_sig;
  std::string signame = "";
  int nbins = 50;

  if (era!="20UL16APV" && era!="20UL16" && era!="20UL17" && era!="20UL18") {
    std::cout << "check params!" << std::endl;
    throw;
  }

  if (ang=="DR2"&&etThres=="Et1") {
    file_sig = TFile::Open(("MergedEleMva_"+era+"_H200A1.root").c_str());
    signame = "H200A1";
  } else if (ang=="DR2"&&etThres=="Et2") {
    file_sig = TFile::Open(("MergedEleMva_"+era+"_H800A10.root").c_str());
    signame = "H800A10";
  } else if (ang=="DR1"&&etThres=="Et2") {
    file_sig = TFile::Open(("MergedEleMva_"+era+"_H800A1.root").c_str()); // no EE
    signame = "H800A1 (w/ GSF)";
  } else if (ang==""&&etThres=="Et2") {
    file_sig = TFile::Open(("MergedEleMva_"+era+"_H800A1.root").c_str()); // no EE
    signame = "H800A1 (w/o GSF)";
    nbins = 100;
  } else {
    std::cout << "check params!" << std::endl;
    throw;
  }

  std::vector<sample> bkglist = {
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT50to100.root"),187700000.0,kCyan-10),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT100to200.root"),23640000.0,kCyan-9),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT200to300.root"),1546000.0,kCyan-7),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT300to500.root"),321600.0,kCyan-6),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT500to700.root"),30980.0,kCyan-3),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT700to1000.root"),6364.0,kCyan-2),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT1000to1500.root"),1117.0,kCyan+2),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT1500to2000.root"),108.4,kCyan+2),
    // sample(TFile::Open("MergedEleMva_20UL16_QCD_HT2000toInf.root"),22.36,kCyan+4),

    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-0To70.root").c_str()),53870.0,kGreen),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-70To100.root").c_str()),1283.0,kGreen-7),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-100To200.root").c_str()),1244.0,kGreen+1),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-200To400.root").c_str()),337.8,kGreen-6),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-400To600.root").c_str()),44.93,kGreen-8),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-600To800.root").c_str()),11.19,kGreen-5),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-800To1200.root").c_str()),4.926,kGreen+2),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-1200To2500.root").c_str()),1.152,kGreen+3),
    sample(TFile::Open(("MergedEleMva_"+era+"_WJets_HT-2500ToInf.root").c_str()),0.02646,kGreen+4),

    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-0To70.root").c_str()),5379.0,kOrange),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-70To100.root").c_str()),140.0,kOrange-4),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-100To200.root").c_str()),139.2,kOrange+6),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-200To400.root").c_str()),38.4,kOrange-3),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-400To600.root").c_str()),5.174,kOrange+7),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-600To800.root").c_str()),1.258,kOrange-5),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-800To1200.root").c_str()),0.5598,kOrange+5),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-1200To2500.root").c_str()),0.1305,kOrange-6),
    sample(TFile::Open(("MergedEleMva_"+era+"_DY_HT-2500ToInf.root").c_str()),0.002997,kOrange+4),
    sample(TFile::Open(("MergedEleMva_"+era+"_TT.root").c_str()),471.7,kRed+1) // 831.76
  };

  auto retrieveTrees = [](const TString aname, std::vector<sample> list) {
    std::vector<TTree*> result;

    for (const auto& element : list)
      result.push_back((TTree*)(element.file)->Get(aname));

    return result;
  };

  auto retrieveHists = [](const TString aname, std::vector<sample> list) {
    std::vector<TH1D*> result;

    for (const auto& element : list)
      result.push_back((TH1D*)(element.file)->Get(aname));

    return result;
  };

  TTree* mergedElTree = nullptr;
  TTree* mergedElGsfTree = nullptr;
  std::vector<TTree*> fakeTrees;
  std::vector<TTree*> fakeGsfTrees;

  if (!ang.empty()) {
    mergedElTree = (TTree*)file_sig->Get("mergedEleSigMvaInput/mergedEl1_elTree");
    mergedElGsfTree = (TTree*)file_sig->Get("mergedEleSigMvaInput/mergedEl1Gsf_addGsfTree");
    fakeTrees = retrieveTrees("mergedEleBkgMvaInput/fake_elTree", bkglist);
    fakeGsfTrees = retrieveTrees("mergedEleBkgMvaInput/fakeGsf_addGsfTree", bkglist);
  } else {
    mergedElTree = (TTree*)file_sig->Get("mergedEleSigMvaInput/mergedEl2_elTree");
    fakeTrees = retrieveTrees("mergedEleBkgMvaInput/bkg_elTree", bkglist);
  }

  auto fillElectrons = [&](TTree* tr, TTree* gtr, TString postfix) -> std::vector<std::shared_ptr<TH1D>> {
    auto h_1x5o5x5 = std::make_shared<TH1D>("E1x5o5x5_"+postfix,"E1x5/5x5;E1x5/5x5;nEle",nbins,0.,1.2);
    auto h_r9 = std::make_shared<TH1D>("R9_"+postfix,"R9;;nEle",nbins,0.,1.);
    auto h_dphiIn = std::make_shared<TH1D>("dPhiIn_"+postfix,"#Delta#phi_{in};#Delta#phi_{in};nEle",nbins,-0.1,0.1);
    auto h_EOP = std::make_shared<TH1D>("EoverP_"+postfix,"E/p;E/p;nEle",nbins,0.,10.);
    auto h_fbrem = std::make_shared<TH1D>("fbrem_"+postfix,"f_{brem};f_{brem};nEle",nbins,-1.,1.);

    auto h_sigieie = std::make_shared<TH1D>("sigieie_"+postfix,"#sigma_{i#eta i#eta};;nEle",nbins,0.,0.05);
    auto h_sigipip = std::make_shared<TH1D>("sigipip_"+postfix,"#sigma_{i#phi i#phi};;nEle",nbins,0.,0.08);
    auto h_HoE = std::make_shared<TH1D>("HoE_"+postfix,"h/E;;nEle",nbins,0.,0.1);
    auto h_HcalD1iso = std::make_shared<TH1D>("HcalD1iso_"+postfix,"HCAL depth1 iso;;nEle",nbins,0.,20.);

    auto h_detaInSeed = std::make_shared<TH1D>("dEtaInSeed_"+postfix,"#Delta#eta_{in, seed};#Delta#eta_{in, seed};nEle",nbins,-0.05,0.05);

    ElectronStruct elstruct;
    tr->SetBranchAddress("ElectronStruct",&elstruct);

    AddGsfStruct addgsfstruct;

    if (gtr)
      gtr->SetBranchAddress("AddGsfStruct",&addgsfstruct);

    for (int ievt = 0; ievt < tr->GetEntries(); ievt++) {
      tr->GetEntry(ievt);

      if (gtr)
        gtr->GetEntry(ievt);

      double weight = elstruct.weight/std::abs(elstruct.weight);

      if (EBEE=="EB") {
        if (std::abs(elstruct.etaSC) > 1.5)
          continue;
      } else {
        if (std::abs(elstruct.etaSC) < 1.5)
          continue;
      }

      bool isEB = (std::abs(elstruct.etaSC) < 1.5);

      double dr2 = 0.;

      if (gtr)
        dr2 = deltaR2(addgsfstruct.Gsfeta,addgsfstruct.Gsfphi,elstruct.Gsfeta,elstruct.Gsfphi);

      if (ang=="DR1") {
        if (isEB) {
          if ( dr2 > 0.0174*0.0174 )
            continue;
        } else {
          if ( dr2 > (0.00864*std::abs(std::sinh(elstruct.eta)))*(0.00864*std::abs(std::sinh(elstruct.eta))) )
            continue;
        }
      } else if (ang=="DR2") {
        if (isEB) {
          if ( dr2 < 0.0174*0.0174 )
            continue;
        } else {
          if ( dr2 < (0.00864*std::abs(std::sinh(elstruct.eta)))*(0.00864*std::abs(std::sinh(elstruct.eta))) )
            continue;
        }
      } else if (ang=="") {

      } else {
        std::cout << "check opening angle param!" << std::endl;
        throw;
      }

      if (elstruct.et < 50)
        continue;

      // if (etThres=="Et1") {
      //   if (isEB) {
      //     if (elstruct.et > 200 || elstruct.et < 50)
      //       continue;
      //   } else {
      //     if (elstruct.et > 150 || elstruct.et < 50)
      //       continue;
      //   }
      // } else if (etThres=="Et2") {
      //   if (isEB) {
      //     if (elstruct.et < 200)
      //       continue;
      //   } else {
      //     if (elstruct.et < 150)
      //       continue;
      //   }
      // } else {
      //   std::cout << "check Et threshold param!" << std::endl;
      //   throw;
      // }

      h_1x5o5x5->Fill(elstruct.full5x5_E1x5/elstruct.full5x5_E5x5,weight);
      h_r9->Fill(elstruct.full5x5_r9,weight);
      h_dphiIn->Fill(elstruct.dPhiIn,weight);
      h_EOP->Fill(elstruct.EOverP,weight);
      h_fbrem->Fill(elstruct.fbrem,weight);
      h_sigieie->Fill(elstruct.full5x5_sigmaIetaIeta,weight);
      h_sigipip->Fill(elstruct.full5x5_sigmaIphiIphi,weight);
      h_HoE->Fill(elstruct.full5x5_hOverE,weight);
      h_HcalD1iso->Fill(elstruct.dr03HcalDepth1TowerSumEt,weight);
      h_detaInSeed->Fill(elstruct.dEtaSeed,weight);
    }

    std::vector<std::shared_ptr<TH1D>> hists;
    hists.push_back(std::move(h_1x5o5x5));
    hists.push_back(std::move(h_r9));
    hists.push_back(std::move(h_dphiIn));
    hists.push_back(std::move(h_EOP));
    hists.push_back(std::move(h_fbrem));
    hists.push_back(std::move(h_sigieie));
    hists.push_back(std::move(h_sigipip));
    hists.push_back(std::move(h_HoE));
    hists.push_back(std::move(h_HcalD1iso));
    hists.push_back(std::move(h_detaInSeed));

    return std::move(hists);
  };

  auto scaling = [&](TTree* tr, std::vector<std::shared_ptr<TH1D>>& hists, int color, float xsec=0., double sumwgt=0.) {
    int nEntries = tr->GetEntries();

    for (unsigned idx = 0; idx < hists.size(); idx++) {
      hists.at(idx)->Scale(xsec/sumwgt);
      hists.at(idx)->SetFillColor(color);
      hists.at(idx)->SetLineColor(color);
      hists.at(idx)->SetMarkerSize(0);
    }
  };

  auto stacking = [&](std::vector<std::unique_ptr<THStack>>& stacks, std::initializer_list<std::vector<std::shared_ptr<TH1D>>> aBkgHists) {
    for (unsigned idx = 0; idx < stacks.size(); idx++) {
      for (auto& element : aBkgHists)
        stacks.at(idx)->Add(element.at(idx).get());
    }
  };

  auto sumwgt_bkgs = retrieveHists("mergedEleBkgMvaInput/totWeightedSum", bkglist);

  double sumwgt_sig = ((TH1D*)(file_sig->Get("mergedEleSigMvaInput/totWeightedSum")))->GetBinContent(1);

  auto sigHists = fillElectrons(mergedElTree,mergedElGsfTree,"sig");

  std::vector<std::vector<std::shared_ptr<TH1D>>> histowner;
  histowner.reserve(bkglist.size());

  for (unsigned idx = 0; idx < bkglist.size(); idx++) {
    std::string numstr = std::to_string(idx);
    std::vector<std::shared_ptr<TH1D>> hists;

    if (!ang.empty())
      hists = fillElectrons(fakeTrees.at(idx),fakeGsfTrees.at(idx),static_cast<std::string>("fake")+numstr);
    else
      hists = fillElectrons(fakeTrees.at(idx),nullptr,static_cast<std::string>("bkg")+numstr);

    histowner.emplace_back(std::move(hists));
  }

  for (unsigned idx = 0; idx < bkglist.size(); idx++)
    scaling(fakeTrees.at(idx),
            histowner.at(idx),
            bkglist.at(idx).color,
            bkglist.at(idx).xsec,
            sumwgt_bkgs.at(idx)->GetBinContent(1));

  std::vector<std::unique_ptr<THStack>> stacks;

  for (unsigned idx = 0; idx < histowner.front().size(); idx++)
    stacks.emplace_back(std::make_unique<THStack>(
      TString(histowner.front().at(idx)->GetName())+"_stack",
      TString(histowner.front().at(idx)->GetTitle())+";"+TString(histowner.front().at(idx)->GetXaxis()->GetTitle())+";"
    ));

  // set this 0 and uncomment below lines if you want to include QCD MC
  const int offset = 9;

  stacking(stacks,{
    histowner.at(27-offset),

    histowner.at(26-offset),
    histowner.at(25-offset),
    histowner.at(24-offset),
    histowner.at(23-offset),
    histowner.at(22-offset),
    histowner.at(21-offset),
    histowner.at(20-offset),
    histowner.at(19-offset),
    histowner.at(18-offset),

    histowner.at(17-offset),
    histowner.at(16-offset),
    histowner.at(15-offset),
    histowner.at(14-offset),
    histowner.at(13-offset),
    histowner.at(12-offset),
    histowner.at(11-offset),
    histowner.at(10-offset),
    histowner.at(9-offset)

    // histowner.at(8),
    // histowner.at(7),
    // histowner.at(6),
    // histowner.at(5),
    // histowner.at(4),
    // histowner.at(3),
    // histowner.at(2),
    // histowner.at(1),
    // histowner.at(0)
  });

  auto draw = [&](std::vector<std::shared_ptr<TH1D>>& sigHists,
                  std::vector<std::unique_ptr<THStack>>& bkgHists) {
    auto canvas = std::make_unique<TCanvas>("canvas","canvas",50,50,W,H);
    canvas->SetLogy();

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

    auto legend = std::make_unique<TLegend>(0.68,0.75,0.9,0.9);
    legend->SetNColumns(2);
    legend->SetBorderSize(0);
    legend->AddEntry(sigHists.at(0).get(),signame.c_str());
    // legend->AddEntry(histowner.at(4).at(0).get(),"QCD");
    legend->AddEntry(histowner.at(13-offset).at(0).get(),"W");
    legend->AddEntry(histowner.at(22-offset).at(0).get(),"DY");
    legend->AddEntry(histowner.at(27-offset).at(0).get(),"TT");

    for (unsigned idx = 0; idx < sigHists.size(); idx++) {
      double theMax = 0.;

      bkgHists.at(idx)->Draw("hist");
      theMax = std::max(bkgHists.at(idx)->GetMaximum(),theMax);

      auto* tmpHist = (TH1D*)bkgHists.at(idx)->GetHists()->At(0)->Clone();

      for (size_t jdx = 1; jdx < bkgHists.at(idx)->GetHists()->GetSize(); jdx++)
        tmpHist->Add((TH1D*)bkgHists.at(idx)->GetHists()->At(jdx));

      double integral = tmpHist->Integral();
      sigHists.at(idx)->Scale(integral/sigHists.at(idx)->Integral());
      sigHists.at(idx)->SetLineColor(kRed);
      sigHists.at(idx)->SetLineWidth(2);
      sigHists.at(idx)->SetMarkerSize(0);

      sigHists.at(idx)->Draw("hists&same");
      theMax = std::max(sigHists.at(idx)->GetMaximum(),theMax);

      bkgHists.at(idx)->SetMaximum(10.*theMax);
      bkgHists.at(idx)->SetMinimum(10e-5*theMax);

      legend->Draw();

      canvas->Update();

      // writing the lumi information and the CMS "logo"
      CMS_lumi( canvas.get(), iPeriod, iPos );

      canvas->Update();
      canvas->RedrawAxis();
      canvas->GetFrame()->Draw();

      canvas->SaveAs(TString("MvaInputs/")+TString(ang)+TString(etThres)+TString(EBEE)+"_"+sigHists.at(idx)->GetName()+".png");
    } // histograms
  };

  draw(sigHists,stacks);

  return;
}
