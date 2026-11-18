#!/usr/bin/env python3
"""DY/W+jets: prompt efficiency & fake rate for 3 IDs (V70, modHEEP, modHEEP Run3).

Run-3-dedicated, self-contained (does not touch the existing Study-2 scripts).
Two SEPARATE measurements, EB and EE split:
  - prompt efficiency: reco e matched to a hard-process final-state gen e (dR<0.001),
                       passing rate of each ID.
  - fake rate:         reco e with NO eligible gen e (dR>=0.1), passing rate of each ID.
Each ID's "pass" = its 12-bit VID bitmap == 4095. Unweighted; DY components summed.
NOT PF-isolated; no W/Z ancestry (prompt/fake are gen-match proxies); minLeptonsSkim=2 bias.

  source env.sh
  python3 -u dyw_run3.py --max-files 100     # or 0 for all files
  python3 dyw_run3.py --plots-only
Input production: ../26.10.08_22preEE (schema v6, has ModifiedHEEPRun3ID). Baseline reco pT>20.
"""
import argparse, csv, glob, re
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True); ROOT.gStyle.SetOptStat(0); ROOT.TH1.AddDirectory(False)
BASE = Path(__file__).resolve().parent
PNFS = '/pnfs/knu.ac.kr/data/cms/store/user/yeo'
ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument('--prod', default=str(BASE/'..'/'26.10.08_22preEE'))
ap.add_argument('--max-files', type=int, default=100)
ap.add_argument('--plots-only', action='store_true')
args = ap.parse_args()
tag = 'all' if args.max_files == 0 else f'{args.max_files}files'
OUT = BASE/'out'/f'dyw_run3_{tag}'; IMG = BASE/'imgs'/f'dyw_run3_{tag}'
OUT.mkdir(parents=True, exist_ok=True); IMG.mkdir(parents=True, exist_ok=True)

# name -> primary dataset dir, from the production inputdataMC.dat
components = {'DY': ['DY_Jets','DY_Jets_ext1','DY_Jets_ext2'], 'W': ['W_Jets']}
IDS = [('v70','electron_heepV70Bitmap','HEEP V70',ROOT.kBlack,20),
       ('mod','electron_modHeepBitmap','modified HEEP',ROOT.kBlue+1,24),
       ('run3','electron_bitmapModifiedHEEPRun3ID','modified HEEP Run3',ROOT.kRed+1,21)]
POPS = [('prompt',0.001,'false','prompt efficiency'), ('fake',0.1,'true','fake rate')]
cut_labels = ['Baseline','ET>20','SC #eta','#Delta#eta seed','#Delta#phi','#sigma_{i#etai#eta}',
              'E ratios','H/E','Track iso','ECAL+HCAL iso','d_{xy}','Missing hits','ECAL driven']

ROOT.gInterpreter.Declare(r'''
#include <ROOT/RVec.hxx>
#include <cmath>
#include <algorithm>
using ROOT::VecOps::RVec;
// region: 0 reject, 1 EB, 2 EE ; baked-in reco pT>20 and prompt/fake gen matching.
RVec<int> regionPop(const RVec<float>& pt,const RVec<float>& eta,const RVec<float>& phi,
 const RVec<float>& sc,const RVec<float>& ge,const RVec<float>& gf,const RVec<int>& pid,
 double drCut,bool unmatched){
  RVec<int> out(pt.size(),0);
  for(size_t i=0;i<pt.size();++i){
    if(pt[i]<=20) continue;
    double best=1e9;
    for(size_t j=0;j<pid.size();++j) if(std::abs(pid[j])==11 && std::abs(ge[j])<2.5)
      best=std::min(best,std::hypot(double(eta[i]-ge[j]),std::remainder(double(phi[i]-gf[j]),2*M_PI)));
    if(unmatched ? best<drCut : best>=drCut) continue;
    double a=std::abs(sc[i]);
    if(a<1.4442) out[i]=1; else if(a>1.566 && a<2.5) out[i]=2;
  }
  return out;
}
RVec<float> basePts(const RVec<float>& pt,const RVec<int>& reg,int r){
  RVec<float> o; for(size_t i=0;i<pt.size();++i) if(reg[i]==r) o.push_back(pt[i]); return o; }
RVec<float> passPts(const RVec<float>& pt,const RVec<int>& reg,const RVec<int>& bits,int r){
  RVec<float> o; for(size_t i=0;i<pt.size();++i) if(reg[i]==r && (bits[i]&4095)==4095) o.push_back(pt[i]); return o; }
// cumulative cutflow stages (0.5 baseline, then pass-through of cuts 0..c)
RVec<float> cfStages(const RVec<int>& reg,const RVec<int>& bits,int r){
  RVec<float> s; for(size_t i=0;i<reg.size();++i) if(reg[i]==r){
    s.push_back(.5f); int b=bits[i]&4095;
    for(int c=0;c<12;++c){ int m=(1<<(c+1))-1; if((b&m)!=m) break; s.push_back(c+1.5f); } }
  return s; }
''')

_m = re.search(r'(\d{2})\.(\d{2})\.(\d{2})', Path(args.prod).name)
TSPREFIX = ''.join(_m.groups()) if _m else '261010'   # production timestamp prefix, e.g. 261010

def resolve(name, dataset):
    primary = dataset.split('/')[1]
    ts = sorted(glob.glob(f'{PNFS}/{primary}/muonJpsi_{name}/{TSPREFIX}*'))
    if not ts: raise RuntimeError(f'no 261008 dir for {name}')
    files = sorted(glob.glob(f'{ts[-1]}/**/Events*.root', recursive=True))
    if args.max_files and len(files) > args.max_files:
        idx = [round(i*(len(files)-1)/max(1,args.max_files-1)) for i in range(args.max_files)]
        files = [files[i] for i in idx]
    return files

pt_edges = ROOT.std.vector('double')([20,25,30,35,40,50,70,100,150,200,300,500,1000])

if not args.plots_only:
    datasets = {}
    for line in (Path(args.prod)/'inputdataMC.dat').read_text().splitlines():
        p = line.split()
        if len(p) >= 2: datasets[p[0]] = p[1]
    for group, names in components.items():
        for name in names:
            files = resolve(name, datasets[name])
            (OUT/f'{name}_files.txt').write_text('\n'.join(files)+'\n')
            print(f'{name}: {len(files)} files', flush=True)
            vfiles = ROOT.std.vector('string')()
            for f in files: vfiles.push_back(f)
            df = ROOT.RDataFrame('resolvedNtuplizer/Events', vfiles)
            jobs = []
            for pop, dr, unm, _ in POPS:
                df = df.Define(f'reg_{pop}',
                   f'regionPop(electron_ptRaw,electron_etaRaw,electron_phiRaw,electron_etaSC,gen_eta,gen_phi,gen_pdgId,{dr},{unm})')
                for r, region in [(1,'EB'),(2,'EE')]:
                    base = f'base_{pop}_{region}'
                    df = df.Define(base, f'basePts(electron_ptRaw,reg_{pop},{r})')
                    jobs.append(df.Histo1D((base,';p_{T} [GeV];e',pt_edges.size()-1,pt_edges.data()),base))
                    for idn, br, *_ in IDS:
                        num = f'num_{pop}_{region}_{idn}'
                        df = df.Define(num, f'passPts(electron_ptRaw,reg_{pop},{br},{r})')
                        jobs.append(df.Histo1D((num,';p_{T} [GeV];e',pt_edges.size()-1,pt_edges.data()),num))
                        cf = f'cf_{pop}_{region}_{idn}'
                        df = df.Define(cf, f'cfStages(reg_{pop},{br},{r})')
                        jobs.append(df.Histo1D((cf,'',13,0,13),cf))
            fo = ROOT.TFile(str(OUT/f'{name}.root'),'RECREATE')
            for h in jobs: h.GetValue().Write()
            fo.Close(); print(f'{name}: done', flush=True)

# sum components into DY / W
hist = {}
for group, names in components.items():
    for name in names:
        f = ROOT.TFile.Open(str(OUT/f'{name}.root'))
        if not f or f.IsZombie(): raise RuntimeError(f'missing {name}.root')
        for k in f.GetListOfKeys():
            h = k.ReadObj(); key = (group, k.GetName())
            if key not in hist: hist[key] = h.Clone(group+'_'+k.GetName()); hist[key].SetDirectory(0)
            else: hist[key].Add(h)
        f.Close()

def save(c, name):
    for e in ('png','pdf'): c.SaveAs(str(IMG/f'{name}.{e}'))
    c.Close()

keep = []; rate_rows = []
for group in components:
    for pop, dr, unm, poplabel in POPS:
        for region in ('EB','EE'):
            den = hist[(group, f'base_{pop}_{region}')]
            # rate vs pT
            c = ROOT.TCanvas('r','',900,700); c.SetLogx(); c.SetGridx(); c.SetGridy()
            leg = ROOT.TLegend(.5,.76,.9,.9); leg.SetBorderSize(0); leg.SetFillStyle(0)
            graphs = []
            for idn, br, idlabel, col, mk in IDS:
                num = hist[(group, f'num_{pop}_{region}_{idn}')]
                if not ROOT.TEfficiency.CheckConsistency(num, den): raise RuntimeError('bad counts')
                g = ROOT.TGraphAsymmErrors(num, den, 'cl=0.682689492 cp')
                g.SetLineColor(col); g.SetMarkerColor(col); g.SetMarkerStyle(mk); g.SetMarkerSize(.9)
                g.SetTitle(f'{group}+jets | {region} | {poplabel};Electron raw p_{{T}} [GeV];{poplabel}')
                graphs.append(g); leg.AddEntry(g, idlabel, 'lp')
            graphs[0].SetMinimum(0); graphs[0].SetMaximum(1.05)
            graphs[0].Draw('AP'); graphs[0].GetXaxis().SetLimits(20,1000)
            for g in graphs[1:]: g.Draw('P SAME')
            leg.Draw(); keep += graphs+[leg]
            save(c, f'{group}_{pop}_{region}_rate_vs_pt')
            # cutflow (3 IDs)
            c = ROOT.TCanvas('cf','',1000,700); c.SetGridy(); c.SetBottomMargin(.2)
            leg = ROOT.TLegend(.55,.74,.9,.89); leg.SetBorderSize(0); leg.SetFillStyle(0)
            first = True
            for idn, br, idlabel, col, mk in IDS:
                h = hist[(group, f'cf_{pop}_{region}_{idn}')].Clone('d_'+idn+group+pop+region)
                n = h.GetBinContent(1)
                for b in range(1, 14):
                    v = h.GetBinContent(b)
                    h.SetBinContent(b, v/n if n else 0); h.SetBinError(b, (v/n*(1-v/n)/n)**.5 if n and 0<=v/n<=1 else 0)
                for b, lab in enumerate(cut_labels, 1): h.GetXaxis().SetBinLabel(b, lab)
                h.GetXaxis().LabelsOption('v'); h.GetYaxis().SetTitle('Cumulative passing fraction')
                h.SetMinimum(0); h.SetMaximum(1.08); h.SetLineColor(col); h.SetMarkerColor(col)
                h.SetMarkerStyle(mk); h.SetMarkerSize(.8); h.SetLineWidth(2)
                h.SetTitle(f'{group}+jets | {region} | {poplabel} cutflow')
                h.Draw('P' if first else 'P SAME'); first = False
                leg.AddEntry(h, f'{idlabel} ({int(hist[(group,f"cf_{pop}_{region}_{idn}")].GetBinContent(13))}/{int(n)})', 'lp')
                keep.append(h)
            leg.Draw(); keep.append(leg)
            save(c, f'{group}_{pop}_{region}_cutflow')
            # integrated rates CSV (pT>20 and pT>35)
            for ptmin in (20, 35):
                lo = den.FindBin(ptmin+1e-6)
                d = den.Integral(lo, den.GetNbinsX()+1)
                for idn, br, idlabel, *_ in IDS:
                    num = hist[(group, f'num_{pop}_{region}_{idn}')]
                    nn = num.Integral(lo, num.GetNbinsX()+1)
                    rate = nn/d if d else ''
                    lohi = (ROOT.TEfficiency.ClopperPearson(int(d),int(nn),.682689492,False),
                            ROOT.TEfficiency.ClopperPearson(int(d),int(nn),.682689492,True)) if d else ('','')
                    rate_rows.append([group,pop,region,ptmin,idn,int(d),int(nn),rate,lohi[0],lohi[1]])

with (OUT/'rates.csv').open('w') as f:
    w = csv.writer(f)
    w.writerow(['group','population','region','pt_min','ID','denominator','numerator','rate','lo68','hi68'])
    w.writerows(rate_rows)
print(f'Plots: {IMG}\nCounts: {OUT}', flush=True)
