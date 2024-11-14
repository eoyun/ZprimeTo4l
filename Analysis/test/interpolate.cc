#include "TFile.h"
#include "TH1D.h"

#include "RooWorkspace.h"

#include <regex>

#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/tdrstyle.C"
#include "/home/ko/Desktop/Study/Zprime/ZprimeTo4l/work/CMS_lumi.C"

void interpolate(TString filename) {
  setTDRStyle();

  writeExtraText = true;       // if extra text
  extraText  = "   Simulation Internal";  // default extra text is "Preliminary"

  lumi_sqrtS = "13 TeV";
  lumi_13TeV = "";

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

  TFile* afile = TFile::Open(filename,"READ");
  auto keyList = afile->GetListOfKeys();

  TFile* bfile = new TFile(filename.ReplaceAll(".root","_rebin.root"),"RECREATE");

  for (auto keyObj : *keyList) {
    auto key = (TKey*)keyObj;

    if (TString(key->GetClassName())==TString("TDirectoryFile")) {
      TDirectory* adir = (TDirectory*)key->ReadObj();
      TString dirname = adir->GetName();
      auto keyList2 = adir->GetListOfKeys();

      bfile->cd();
      TDirectory* bdir = bfile->mkdir(dirname);
      bdir->cd();

      std::vector<TString> listSyst;
      std::vector<std::pair<int,int>> masspoints;

      std::vector<double> binEdges;
      std::string rooVarRange = "x[0,2500]";

      if (dirname==TString("mergedEMu2M") || dirname==TString("mergedEle3E")) {
        binEdges = {200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,300.,320.,340.,360.,380.,400.,420.,440.,460.,480.,500.,
                    550.,600.,700.,800.,1000.,2500.};
      } else if (dirname.Contains("resolved") || dirname==TString("mergedEle2E")) {
        binEdges = {200.,210.,220.,230.,240.,250.,260.,270.,280.,290.,
                    300.,310.,320.,330.,340.,350.,360.,370.,380.,390.,
                    400.,410.,420.,430.,440.,450.,460.,470.,480.,490.,
                    500.,520.,540.,560.,580.,600.,620.,640.,660.,680.,700.,720.,740.,760.,780.,
                    800.,850.,900.,950.,
                    1000.,1100.,1200.,1500.,2500.};
      } else if (dirname==TString("mergedEMu1M")) {
        binEdges = {250.,275.,300.,325.,350.,375.,400.,425.,450.,475.,500.,
                    550.,600.,700.,800.,1000.,2500.};
      } else if (dirname==TString("mergedMu3M")) {
        binEdges = {250.,300.,350.,400.,450.,500.,
                    600.,700.,800.,2500.};
      }

      for (auto keyObj2 : *keyList2) {
        auto key2 = (TKey*)keyObj2;

        if (TString(key2->GetClassName()).Contains("TH1")) {
          std::string keyname = std::string(key2->GetName());
          std::regex sigPattern("H([0-9]+)A([0-9]+).*");
          std::regex sigSystPattern("H[0-9]+A[0-9]+_(.*)");
          std::smatch XYmass, systName;

          bool isSig = std::regex_match(keyname, XYmass, sigPattern);
          bool isSyst = std::regex_match(keyname, systName, sigSystPattern);

          int mx = 0, my = 0;

          if (isSig) {
            mx = std::stoi(XYmass[1].str());
            my = std::stoi(XYmass[2].str());
          }

          if (isSig && !isSyst)
            masspoints.push_back( std::make_pair(mx,my) );

          if ( isSyst && std::find(listSyst.begin(), listSyst.end(),TString(systName[1].str()))==listSyst.end() )
            listSyst.push_back(TString(systName[1].str()));

          if (!isSig && !isSyst) {
            TH1* ahist = (TH1*)key2->ReadObj();
            auto temp = std::unique_ptr<TH1>((TH1*)ahist->Clone());
            TH1* bhist = temp->Rebin(binEdges.size()-1,key2->GetName(),&(binEdges[0]));
            bdir->WriteTObject(bhist,key2->GetName());
          }
        }
      }

      // std::vector<double> masses = {250., 750., 2000.}; // {1.,2.,5.,10.,50.,100.};
      std::vector<double> masses = {1.,2.,5.,10.,50.,100.};
      listSyst.push_back("");

      for (const auto syst : listSyst) {

      for (const double mass : masses) {

        auto canvas = std::make_unique<TCanvas>("canvas","canvas",800,800,W,H);
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

        auto work = std::make_unique<RooWorkspace>("w");
        work->factory(rooVarRange.c_str());
        work->factory("m[250,2000]");
//        work->factory("m[1,100]");

        unsigned massLen = 0;

        for (auto apoint : masspoints) {
          if (static_cast<double>(apoint.second)==mass) //
            massLen++;
        }

        auto pdfs = RooArgList();
        auto paramVec = TVectorD(massLen);

        unsigned idxy = 0;
        unsigned nbins = 0;

        for (auto keyObj2 : *keyList2) {
          auto key2 = (TKey*)keyObj2;

          if (TString(key2->GetClassName()).Contains("TH1")) {
            std::string keyname = std::string(key2->GetName());
            std::regex sigPattern("H([0-9]+)A([0-9]+).*");
            std::regex sigSystPattern("H[0-9]+A[0-9]+_(.*)");
            std::smatch XYmass, systName;

            bool isSig = std::regex_match(keyname, XYmass, sigPattern);
            bool isSyst = std::regex_match(keyname, systName, sigSystPattern);

            int mx = 0, my = 0;

            if (isSig) {
              mx = std::stoi(XYmass[1].str());
              my = std::stoi(XYmass[2].str());
            }

            bool pickPoints = static_cast<double>(my)==mass; //
            bool pickSysts = (isSyst && TString(systName[1])==syst) || (!isSyst && syst=="");

            if (!pickPoints)
              continue;

            if (!pickSysts)
              continue;

            TH1* ahist = (TH1*)key2->ReadObj();
            TH1* bhist = (TH1*)ahist->Clone();

            if (nbins==0)
              nbins=ahist->GetNbinsX();
            else if (nbins!=ahist->GetNbinsX())
              throw std::invalid_argument("uneven histogram binning!");

            TH1* temp = (TH1*)bhist->Clone();
            TH1* chist = temp->Rebin(binEdges.size()-1,key2->GetName(),&(binEdges[0]));
            // TH1* chist = (TH1*)ahist->Clone();
            bdir->WriteTObject(chist,key2->GetName());

            if (isSig) {
              bhist->Scale(1000.);
              bhist->SetMaximum(40.*bhist->GetMaximum());

              if (idxy==0)
                bhist->Draw("hist");
              else
                bhist->Draw("hist&same");

              RooDataHist rooHist(key2->GetName(),key2->GetName(),*work->var("x"),bhist);
              work->import(rooHist);

              work->factory((std::string("HistPdf::")+key2->GetName()+"_pdf(x,"+key2->GetName()+")").c_str());
              auto* pdf = work->pdf((std::string(key2->GetName())+"_pdf").c_str());
              pdfs.add(*pdf);
              paramVec[idxy] = static_cast<double>(mx); //
              idxy++;
            }
          }
        } // loop hists

        RooArgList varlist;
        varlist.add(*work->var("x"));

        // std::vector<TString> newmassList = {"1p2","1p4","1p6","1p8","3","4","6","7","8","9","20","30","40"};
        std::vector<TString> newmassList = {"300","350","400","450","500","600","850","1000","1200","1400","1600","1800"};

        for (auto newmass : newmassList) {
          TString decimal = TString(newmass).ReplaceAll("p",".");
          const double valmass = std::stod(decimal.Data());
          TString newName = TString("H") + newmass + "A" + TString(std::to_string( static_cast<int>(mass) ).c_str());
          // TString newName = TString("H") + TString(std::to_string( static_cast<int>(mass) ).c_str()) + "A" + newmass;

          if (syst!="")
            newName += (TString("_") + syst);

          auto* alpha = work->var("m");
          alpha->setVal(valmass);

          RooMomentMorph morph(newName+"_morph",newName+"_morph",*alpha,varlist,pdfs,paramVec,RooMomentMorph::Linear); // NonLinearPosFractions

          TH1* morphedHist = morph.createHistogram(newName+"_createHist",*work->var("x"),RooFit::Binning(nbins,0.,2500.),RooFit::ConditionalObservables({*work->var("x")}));
          morphedHist->SetLineColor(kBlue);
          morphedHist->SetLineWidth(2);
          morphedHist->Sumw2();
          morphedHist->Draw("hist&same");

          TH1* morphedHist_clone = (TH1*)morphedHist->Clone();
          morphedHist_clone->Scale(1./1000.);

          TH1* morphedHist_rebin = morphedHist_clone->Rebin(binEdges.size()-1,newName+"_rebin",&(binEdges[0]));
          bdir->WriteTObject(morphedHist_rebin,newName);
        }

        if (syst=="") {
          canvas->Update();

          CMS_lumi( canvas.get(), iPeriod, iPos );

          canvas->Update();
          canvas->RedrawAxis();
          canvas->GetFrame()->Draw();

          // canvas->SaveAs(TString(std::to_string(static_cast<int>(mass)))+".pdf");
        }
      } // loop masses

      for (auto heavyY : {250.,750.}) {
        unsigned nbins = 0;

        for (auto keyObj2 : *keyList2) {
          auto key2 = (TKey*)keyObj2;

          if (TString(key2->GetClassName()).Contains("TH1")) {
            std::string keyname = std::string(key2->GetName());
            std::regex sigPattern("H([0-9]+)A([0-9]+).*");
            std::regex sigSystPattern("H[0-9]+A[0-9]+_(.*)");
            std::smatch XYmass, systName;

            bool isSig = std::regex_match(keyname, XYmass, sigPattern);
            bool isSyst = std::regex_match(keyname, systName, sigSystPattern);

            int mx = 0, my = 0;

            if (isSig) {
              mx = std::stoi(XYmass[1].str());
              my = std::stoi(XYmass[2].str());
            }

            bool pickPoints = static_cast<double>(my)==heavyY; //
            bool pickSysts = (isSyst && TString(systName[1])==syst) || (!isSyst && syst=="");

            if (!pickPoints)
              continue;

            if (!pickSysts)
              continue;

            TH1* ahist = (TH1*)key2->ReadObj();
            TH1* bhist = (TH1*)ahist->Clone();

            if (nbins==0)
              nbins=ahist->GetNbinsX();
            else if (nbins!=ahist->GetNbinsX())
              throw std::invalid_argument("uneven histogram binning!");

            TH1* temp = (TH1*)bhist->Clone();
            TH1* chist = temp->Rebin(binEdges.size()-1,key2->GetName(),&(binEdges[0]));
            // TH1* chist = (TH1*)ahist->Clone();
            bdir->WriteTObject(chist,key2->GetName());
          }
        } // loop hists
      } // loop masses

      } // loop systs
      
    } // if TDirectoryFile

  } // loop TDirectory

  afile->Close();
  bfile->Close();

  return;
}
