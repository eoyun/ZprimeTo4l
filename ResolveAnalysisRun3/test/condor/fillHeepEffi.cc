// Condor worker: process one sample's ntuple files -> 6 dR_gen histograms.
// Usage: fillHeepEffi <out.root> <filelist.txt>
//        fillHeepEffi <out.root> file1.root [file2.root ...]
// filelist.txt = one ntuple path per line (absolute /pnfs paths work on worker nodes).
#include "HeepEffiCore.h"
#include <TChain.h>
#include <TFile.h>
#include <fstream>
#include <string>
#include <cstdio>
#include <cstdlib>

int main(int argc, char** argv){
  if(argc<3){ printf("usage: %s <out.root> <filelist.txt | file1.root [file2.root ...]>\n", argv[0]); return 2; }
  std::string out = argv[1];

  TChain ch("resolvedNtuplizer/Events");
  std::string a2 = argv[2];
  int nadded=0;
  if(argc==3 && a2.size()>4 && a2.substr(a2.size()-4)==".txt"){
    std::ifstream fin(a2);
    if(!fin){ printf("cannot open filelist %s\n", a2.c_str()); return 1; }
    std::string line;
    while(std::getline(fin,line)){ if(line.empty()||line[0]=='#') continue; ch.Add(line.c_str()); nadded++; }
  } else {
    for(int i=2;i<argc;i++){ ch.Add(argv[i]); nadded++; }
  }
  printf("added %d files, entries=%lld\n", nadded, ch.GetEntries());

  bool req2 = (std::getenv("REQUIRE_2NDTRK") != nullptr);   // denominator = electrons with a 2nd (additional) track
  if(req2) printf("[require 2nd track: denominator restricted to electrons with valid dPerpIn]\n");
  heep::Hists h = heep::makeHists();
  heep::Cutflow cf; heep::processTree(&ch, h, cf, req2);
  heep::printCutflow(cf, h, out.c_str());

  TFile* fo = TFile::Open(out.c_str(),"RECREATE");
  for(TH1D* x : {h.denEB,h.befEB,h.aftEB,h.denEE,h.befEE,h.aftEE,
                 h.denEBpt,h.befEBpt,h.aftEBpt,h.denEEpt,h.befEEpt,h.aftEEpt,
                 h.r3EB,h.r3EE,h.r3EBpt,h.r3EEpt,
                 h.cfEB,h.cfEE,h.cfV70EB,h.cfV70EE,h.cfR3EB,h.cfR3EE,
                 h.cfHiEB,h.cfHiEE,h.cfHiV70EB,h.cfHiV70EE,h.cfHiR3EB,h.cfHiR3EE}) x->Write();
  for(int c=0;c<12;c++) for(int d=0;d<3;d++) for(int rr=0;rr<2;rr++){ h.cut[c][d][rr]->Write(); h.cutpt[c][d][rr]->Write(); }
  for(int rr=0;rr<2;rr++){
    for(int v=0;v<2;v++){ h.isoTrk[v][rr]->Write(); h.drTrk[v][rr]->Write(); }
    for(int v=0;v<3;v++) h.isoEcal[v][rr]->Write();
    h.isoHcal[rr]->Write(); }
  fo->Close();
  printf("wrote %s\n", out.c_str());
  return 0;
}
