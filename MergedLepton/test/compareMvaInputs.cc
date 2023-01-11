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
  float EoverP_1st, EoverP_2nd, dEtaIn_trk1xtal1, dPhiIn_trk1xtal1, dEtaIn_trk2xtal2, dPhiIn_trk2xtal2;
  float E5x5, E1x1;
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

void compareMvaInputs(const std::string ang, const std::string etThres, const std::string EBEE) {
  TH1::SetDefaultSumw2();

  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "Work in progress";  // default extra text is "Preliminary"
  lumi_sqrtS = "13 TeV";       // used with iPeriod = 0, e.g. for simulation-only plots (default is an empty string)

  int iPeriod = 0;    // 1=7TeV, 2=8TeV, 3=7+8TeV, 7=7+8+13TeV, 0=free form (uses lumi_sqrtS)
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

  // const std::string ang = "DR1";
  // const std::string etThres = "Et1";
  // const std::string EBEE = "EB";

  TFile* file_sig;
  std::string signame = "";
  int nbins = 50;

  if (ang=="DR2"&&etThres=="Et1") {
    file_sig = TFile::Open("MergedEleMva_20UL16_H200A1.root");
    signame = "H200A1";
  } else if (ang=="DR2"&&etThres=="Et2") {
    file_sig = TFile::Open("MergedEleMva_20UL16_H800A10.root");
    signame = "H800A10";
  } else if (ang=="DR1"&&etThres=="Et2") {
    file_sig = TFile::Open("MergedEleMva_20UL16_H800A1.root"); // no EE
    signame = "H800A1 (w/ GSF)";
  } else if (ang==""&&etThres=="Et2") {
    file_sig = TFile::Open("MergedEleMva_20UL16_H800A1.root"); // no EE
    signame = "H800A1 (w/o GSF)";
    nbins = 100;
  } else {
    std::cout << "check params!" << std::endl;
    throw;
  }

  std::vector<sample> bkglist = {
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-15to20_EM.root"),1324000.0,kCyan-10),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-20to30_EM.root"),4896000.0,kCyan-9),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-30to50_EM.root"),6447000.0,kCyan-7),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-50to80_EM.root"),1988000.0,kCyan-6),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-80to120_EM.root"),367500.0,kCyan-3),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-120to170_EM.root"),66590.0,kCyan-2),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-170to300_EM.root"),16620.0,kCyan+2),
    sample(TFile::Open("MergedEleMva_20UL16_QCD_Pt-300toInf_EM.root"),1104.0,kCyan+3),
    sample(TFile::Open("MergedEleMva_20UL16_WJets.root"),53870.0,kGreen-2),
    sample(TFile::Open("MergedEleMva_20UL16_DY.root"),6077.22,kOrange),
    sample(TFile::Open("MergedEleMva_20UL16_TT.root"),54.17,kRed+1) // 831.76
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
    auto h_etaSCwidth = std::make_shared<TH1D>("etaSCWidth_"+postfix,"SC #eta width;#Delta#eta_{SC};nEle",nbins,0.,0.03);
    auto h_phiSCwidth = std::make_shared<TH1D>("phiSCWidth_"+postfix,"SC #phi width;#Delta#phi_{SC};nEle",nbins,0.,0.1);
    auto h_1x5o5x5 = std::make_shared<TH1D>("E1x5o5x5_"+postfix,"E1x5/5x5;E1x5/5x5;nEle",nbins,0.,1.);
    auto h_r9 = std::make_shared<TH1D>("R9_"+postfix,"R9;;nEle",nbins,0.,1.);
    auto h_detaIn = std::make_shared<TH1D>("dEtaIn_"+postfix,"#Delta#eta_{in};#Delta#eta_{in};nEle",nbins,-0.05,0.05);
    auto h_dphiIn = std::make_shared<TH1D>("dPhiIn_"+postfix,"#Delta#phi_{in};#Delta#phi_{in};nEle",nbins,-0.1,0.1);
    auto h_EOP = std::make_shared<TH1D>("EoverP_"+postfix,"E/p;E/p;nEle",nbins,0.,10.);
    auto h_fbrem = std::make_shared<TH1D>("fbrem_"+postfix,"f_{brem};f_{brem};nEle",nbins,-1.,1.);

    auto h_EoP1 = std::make_shared<TH1D>("EoP1st_"+postfix,"E(1st xtal)/p(1st GSF);E/p;nEle",nbins,0.,10.);
    auto h_EoP2 = std::make_shared<TH1D>("EoP2nd_"+postfix,"E(2nd xtal)/p(2nd GSF);E/p;nEle",nbins,0.,10.);
    auto h_detaIn_trk1xtal1 = std::make_shared<TH1D>("dEtaIn_trk1xtal1_"+postfix,"#Delta#eta(1st xtal, 1st GSF);#Delta#eta_{in};nEle",nbins,-0.15,0.15);
    auto h_dphiIn_trk1xtal1 = std::make_shared<TH1D>("dPhiIn_trk1xtal1_"+postfix,"#Delta#phi(1st xtal, 1st GSF);#Delta#phi_{in};nEle",nbins,-0.3,0.3);
    auto h_detaIn_trk2xtal2 = std::make_shared<TH1D>("dEtaIn_trk2xtal2_"+postfix,"#Delta#eta(2nd xtal, 2nd GSF);#Delta#eta_{in};nEle",nbins,-0.15,0.15);
    auto h_dphiIn_trk2xtal2 = std::make_shared<TH1D>("dPhiIn_trk2xtal2_"+postfix,"#Delta#phi(2nd xtal, 2nd GSF);#Delta#phi_{in};nEle",nbins,-0.3,0.3);
    auto h_r1 = std::make_shared<TH1D>("R1_"+postfix,"R1;;nEle",nbins,0.,1.);
    auto h_E1x1oE3x3 = std::make_shared<TH1D>("E1x1oE3x3_"+postfix,"E1x1/E3x3;;nEle",nbins,0.,1.);

    auto h_sigieie = std::make_shared<TH1D>("sigieie_"+postfix,"#sigma_{i#eta i#eta};;nEle",nbins,0.,0.05);
    auto h_sigipip = std::make_shared<TH1D>("sigipip_"+postfix,"#sigma_{i#phi i#phi};;nEle",nbins,0.,0.08);
    auto h_HoE = std::make_shared<TH1D>("HoE_"+postfix,"h/E;;nEle",nbins,0.,0.1);
    auto h_HcalD1iso = std::make_shared<TH1D>("HcalD1iso_"+postfix,"HCAL depth1 iso;;nEle",nbins,0.,20.);

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
        if (std::abs(elstruct.eta) > 1.5)
          continue;
      } else {
        if (std::abs(elstruct.eta) < 1.5)
          continue;
      }

      bool isEB = (std::abs(elstruct.eta) < 1.5);

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

      h_etaSCwidth->Fill(elstruct.etaSCWidth,weight);
      h_phiSCwidth->Fill(elstruct.phiSCWidth,weight);
      h_1x5o5x5->Fill(elstruct.full5x5_E1x5/elstruct.full5x5_E5x5,weight);
      h_r9->Fill(elstruct.full5x5_r9,weight);
      h_detaIn->Fill(elstruct.dEtaIn,weight);
      h_dphiIn->Fill(elstruct.dPhiIn,weight);
      h_EOP->Fill(elstruct.EOverP,weight);
      h_fbrem->Fill(elstruct.fbrem,weight);
      h_EoP1->Fill(elstruct.EoverP_1st,weight);
      h_EoP2->Fill(elstruct.EoverP_2nd,weight);
      h_detaIn_trk1xtal1->Fill(elstruct.dEtaIn_trk1xtal1,weight);
      h_dphiIn_trk1xtal1->Fill(elstruct.dPhiIn_trk1xtal1,weight);
      h_detaIn_trk2xtal2->Fill(elstruct.dEtaIn_trk2xtal2,weight);
      h_dphiIn_trk2xtal2->Fill(elstruct.dPhiIn_trk2xtal2,weight);
      h_r1->Fill(elstruct.E1x1/elstruct.E5x5,weight);
      h_E1x1oE3x3->Fill(elstruct.E1x1/(elstruct.E5x5*elstruct.full5x5_r9),weight);
      h_sigieie->Fill(elstruct.full5x5_sigmaIetaIeta,weight);
      h_sigipip->Fill(elstruct.full5x5_sigmaIphiIphi,weight);
      h_HoE->Fill(elstruct.full5x5_hOverE,weight);
      h_HcalD1iso->Fill(elstruct.dr03HcalDepth1TowerSumEt,weight);
    }

    std::vector<std::shared_ptr<TH1D>> hists;
    hists.push_back(std::move(h_etaSCwidth));
    hists.push_back(std::move(h_phiSCwidth));
    hists.push_back(std::move(h_1x5o5x5));
    hists.push_back(std::move(h_r9));
    hists.push_back(std::move(h_detaIn));
    hists.push_back(std::move(h_dphiIn));
    hists.push_back(std::move(h_EOP));
    hists.push_back(std::move(h_fbrem));
    hists.push_back(std::move(h_EoP1));
    hists.push_back(std::move(h_EoP2));
    hists.push_back(std::move(h_detaIn_trk1xtal1));
    hists.push_back(std::move(h_dphiIn_trk1xtal1));
    hists.push_back(std::move(h_detaIn_trk2xtal2));
    hists.push_back(std::move(h_dphiIn_trk2xtal2));
    hists.push_back(std::move(h_r1));
    hists.push_back(std::move(h_E1x1oE3x3));
    hists.push_back(std::move(h_sigieie));
    hists.push_back(std::move(h_sigipip));
    hists.push_back(std::move(h_HoE));
    hists.push_back(std::move(h_HcalD1iso));

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

  stacking(stacks,{
    histowner.at(10),
    histowner.at(9),
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
    legend->AddEntry(histowner.at(3).at(0).get(),"QCD");
    legend->AddEntry(histowner.at(8).at(0).get(),"W");
    legend->AddEntry(histowner.at(9).at(0).get(),"DY");
    legend->AddEntry(histowner.at(10).at(0).get(),"TT");

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
