#!/usr/bin/env python3
"""Run after `source env.sh`; use --plots-only to redraw cached histograms.

Uses the HeepEffiCore.h denominator: greedy nearest opposite-charge gen pairs,
both |eta|<2.5, leading gen pt>20; reco pt>10 and nearest selected gen dR<0.001.
No reconstructed electron ID requirement. Unit event weights; signed matching
variables. Invalid values are counted separately, never folded into the plots.
"""
import argparse
import csv
import re
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
BASE = Path(__file__).resolve().parent
OUT = BASE / 'out/id_inputs'
IMG = BASE / 'imgs/id_inputs'
for p in (OUT, IMG):
    p.mkdir(parents=True, exist_ok=True)
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--plots-only', action='store_true')
args = parser.parse_args()

ROOT.gInterpreter.Declare(r'''
#include <ROOT/RVec.hxx>
#include <cmath>
#include <algorithm>
using ROOT::VecOps::RVec;
double inputDR(double e1,double p1,double e2,double p2) {
  return std::hypot(e1-e2,std::remainder(p1-p2,2*M_PI));
}
RVec<int> inputRegions(const RVec<float>& gp,const RVec<float>& ge,
 const RVec<float>& gf,const RVec<int>& pid,const RVec<float>& ep,
 const RVec<float>& ee,const RVec<float>& ef,const RVec<float>& sc) {
  std::vector<int> pos,neg;
  for(size_t i=0;i<pid.size();++i) {
    if(pid[i]==-11) pos.push_back(i);
    if(pid[i]==11) neg.push_back(i);
  }
  std::vector<bool> usedp(pos.size()),usedn(neg.size()),selected(pid.size());
  while(true) {
    double best=1e9; int bp=-1,bn=-1;
    for(size_t a=0;a<pos.size();++a) if(!usedp[a])
      for(size_t b=0;b<neg.size();++b) if(!usedn[b]) {
        double d=inputDR(ge[pos[a]],gf[pos[a]],ge[neg[b]],gf[neg[b]]);
        if(d<best) {best=d;bp=a;bn=b;}
      }
    if(bp<0) break;
    usedp[bp]=true;usedn[bn]=true;
    int a=pos[bp],b=neg[bn];
    if(std::abs(ge[a])<2.5 && std::abs(ge[b])<2.5 && std::max(gp[a],gp[b])>20)
      selected[a]=selected[b]=true;
  }
  RVec<int> region(ep.size(),0);
  for(size_t i=0;i<ep.size();++i) {
    if(ep[i]<=10) continue;
    double best=1e9;
    for(size_t j=0;j<pid.size();++j) if(selected[j])
      best=std::min(best,inputDR(ee[i],ef[i],ge[j],gf[j]));
    if(best>=0.001) continue;
    double eta=std::abs(sc[i]);
    if(eta<1.4442) region[i]=1;
    else if(eta>1.566 && eta<2.5) region[i]=2;
  }
  return region;
}
RVec<float> inputValues(const RVec<float>& values,const RVec<int>& region,int r) {
  RVec<float> result;
  if(values.size()!=region.size()) throw std::runtime_error("electron vector size mismatch");
  for(size_t i=0;i<values.size();++i)
    if(region[i]==r && std::isfinite(values[i]) && values[i]!=-999.f && std::abs(values[i])<1e30)
      result.push_back(values[i]);
  return result;
}
''')

variables = {
    'dPerpIn': ('electron_dPerpIn', '#Delta u_{5#times5}^{in} (signed)', (-0.03, 0.03), (-0.03, 0.03)),
    'dEtaInSeed': ('electron_dEtaInSeed', '#Delta#eta_{in}^{seed} (signed)', (-0.03, 0.03), (-0.03, 0.03)),
    'sigmaIetaIeta': ('electron_sigmaIeta', 'full5x5 #sigma_{i#etai#eta}', (0., 0.03), (0., 0.10)),
}
samples = []
for fl in BASE.glob('filelists/HToAATo4L_H*A*_22.txt'):
    match = re.fullmatch(r'HToAATo4L_H(\d+)A([\dp]+)_22', fl.stem)
    if match:
        samples.append((int(match[1]), float(match[2].replace('p', '.')), fl))
samples.sort()
if not samples:
    raise RuntimeError('No signal filelists found')

if not args.plots_only:
    rows = []
    for mh, ma, fl in samples:
        paths = [s.strip() for s in fl.read_text().splitlines() if s.strip() and not s.startswith('#')]
        print(f'Processing {fl.stem}: {len(paths)} files', flush=True)
        # Check every file; fail rather than silently dropping missing branches/files.
        for path in paths:
            f = ROOT.TFile.Open(path)
            if not f or f.IsZombie():
                raise RuntimeError(f'Cannot open {path}')
            t = f.Get('resolvedNtuplizer/Events')
            needed = [v[0] for v in variables.values()] + ['gen_pt','gen_eta','gen_phi','gen_pdgId',
                      'electron_ptRaw','electron_etaRaw','electron_phiRaw','electron_etaSC']
            if not t or any(not t.GetBranch(n) for n in needed):
                raise RuntimeError(f'Missing tree/required branches: {path}')
            f.Close()
        files = ROOT.std.vector('string')()
        for path in paths:
            files.push_back(path)
        df = ROOT.RDataFrame('resolvedNtuplizer/Events', files).Define('region',
            'inputRegions(gen_pt,gen_eta,gen_phi,gen_pdgId,electron_ptRaw,electron_etaRaw,electron_phiRaw,electron_etaSC)')
        count = df.Count()
        jobs = []
        for r, region in enumerate(('EB','EE'), 1):
            node = df.Define(f'n_{region}', f'(ULong64_t)Sum(region=={r})')
            total = node.Sum('n_'+region)
            for name, (branch, title, eb, ee) in variables.items():
                col = name+'_'+region
                node = node.Define(col, f'inputValues({branch},region,{r})')
                lo, hi = eb if r==1 else ee
                h = node.Histo1D((col, ';'+title+';Electrons', 120, lo, hi), col)
                jobs.append((name,region,h,total))
        n_events = count.GetValue()  # all booked actions run in one event loop
        output = ROOT.TFile(str(OUT / (fl.stem+'.root')), 'RECREATE')
        for name,region,h,total in jobs:
            hist = h.GetValue()
            hist.Write()
            valid = int(hist.GetEntries())
            rows.append([fl.stem,mh,ma,n_events,name,region,int(total.GetValue()),valid,
                         int(total.GetValue())-valid,int(hist.GetBinContent(0)),
                         int(hist.GetBinContent(hist.GetNbinsX()+1))])
        output.Close()
        print(f'Finished {fl.stem}: {n_events} events', flush=True)
    with (OUT/'counts.csv').open('w') as f:
        writer = csv.writer(f)
        writer.writerow(['sample','mH_GeV','mA_GeV','events','variable','region','selected','valid','invalid','underflow','overflow'])
        writer.writerows(rows)

histograms = {}
for mh,ma,fl in samples:
    f = ROOT.TFile.Open(str(OUT/(fl.stem+'.root')))
    if not f or f.IsZombie():
        raise RuntimeError(f'Missing histogram file for {fl.stem}')
    for name in variables:
        for region in ('EB','EE'):
            h = f.Get(name+'_'+region)
            if not h:
                raise RuntimeError(f'Missing histogram {name}_{region}')
            histograms[mh,ma,name,region] = h.Clone(fl.stem+'_'+name+'_'+region)
    f.Close()

colors = [ROOT.kBlue+1,ROOT.kRed+1,ROOT.kGreen+2,ROOT.kMagenta+1,ROOT.kOrange+7]
def panel(pad, keys, title, labels):
    pad.cd()
    pad.SetLeftMargin(.15); pad.SetBottomMargin(.14); pad.SetTopMargin(.11)
    pad.SetLogy()
    hs=[]
    for i,key in enumerate(keys):
        h=histograms[key].Clone('display_'+str(key))
        # Normalize to all valid electrons, including under/overflow.
        n=h.Integral(0,h.GetNbinsX()+1)
        if n: h.Scale(1/n)
        h.SetLineColor(colors[i%len(colors)]); h.SetLineWidth(2)
        hs.append(h)
    ymax=max([h.GetMaximum() for h in hs]+[0.001])
    legend=ROOT.TLegend(.53,.64,.94,.88)
    legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.033)
    for i,h in enumerate(hs):
        h.SetTitle(title)
        h.GetYaxis().SetTitle('Fraction of valid electrons / bin')
        h.GetYaxis().SetTitleOffset(1.65)
        h.SetMinimum(1e-5);h.SetMaximum(ymax*12)
        h.Draw('HIST' if i==0 else 'HIST SAME')
        legend.AddEntry(h,labels[i], 'l')
    legend.Draw()
    return hs,legend

keep=[]
for mh,ma,fl in samples:
    c=ROOT.TCanvas(fl.stem,fl.stem,1600,900); c.Divide(3,2)
    for row,region in enumerate(('EB','EE')):
        for col,name in enumerate(variables):
            key=(mh,ma,name,region)
            h=histograms[key]
            label=f'N valid = {int(h.GetEntries())}'
            keep.append(panel(c.cd(row*3+col+1),[key],f'm_{{H}}={mh}, m_{{A}}={ma:g} GeV | {region}',[label]))
    for ext in ('png','pdf'): c.SaveAs(str(IMG/(fl.stem+'.'+ext)))
    c.Close()

# Compare A masses at fixed H, and H masses at fixed A.
for fixed in ('H','A'):
    values=sorted({s[0 if fixed=='H' else 1] for s in samples})
    for mass in values:
        group=[s for s in samples if s[0 if fixed=='H' else 1]==mass]
        for name in variables:
            c=ROOT.TCanvas('comparison','comparison',1400,600);c.Divide(2,1)
            for i,region in enumerate(('EB','EE')):
                keys=[(h,a,name,region) for h,a,_ in group]
                labels=[f'm_{{{"A" if fixed=="H" else "H"}}} = {a if fixed=="H" else h:g} GeV' for h,a,_ in group]
                keep.append(panel(c.cd(i+1),keys,f'm_{{{fixed}}}={mass:g} GeV | {region}',labels))
            stem=f'compare_{fixed}{mass:g}_{name}'.replace('.','p')
            for ext in ('png','pdf'): c.SaveAs(str(IMG/(stem+'.'+ext)))
            c.Close()
print(f'Plots: {IMG}\nHistograms and counts: {OUT}', flush=True)
