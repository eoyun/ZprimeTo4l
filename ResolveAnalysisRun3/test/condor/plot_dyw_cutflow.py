#!/usr/bin/env python3
"""Only cumulative VID cutflows for DY/W matched and unmatched electrons.
source env.sh; python3 plot_dyw_cutflow.py [--plots-only]
Uses the exact 400-file manifests of the previous diagnostic.
"""
import argparse
import csv
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--plots-only',action='store_true')
args=parser.parse_args()
base=Path(__file__).resolve().parent
out=base/'out/dyw_cutflow_only';img=base/'imgs/dyw_cutflow_only'
out.mkdir(parents=True,exist_ok=True);img.mkdir(parents=True,exist_ok=True)
ROOT.gInterpreter.Declare(r'''
#include <ROOT/RVec.hxx>
#include <cmath>
#include <algorithm>
using ROOT::VecOps::RVec;
// Sign encodes matching: positive=matched, negative=unmatched; magnitude EB/EE.
RVec<int> cfCategory(const RVec<float>& pt,const RVec<float>& eta,
 const RVec<float>& phi,const RVec<float>& sc,const RVec<float>& ge,
 const RVec<float>& gp,const RVec<int>& id) {
 RVec<int> result(pt.size(),0);
 for(size_t i=0;i<pt.size();++i) {
  if(pt[i]<=20) continue;
  double a=std::abs(sc[i]);int r=a<1.4442?1:(a>1.566 && a<2.5?2:0);
  if(!r) continue;
  double dr=1e9;
  for(size_t j=0;j<id.size();++j) if(std::abs(id[j])==11 && std::abs(ge[j])<2.5)
   dr=std::min(dr,std::hypot(double(eta[i]-ge[j]),std::remainder(double(phi[i]-gp[j]),2*M_PI)));
  result[i]=dr<.1?r:-r;
 }
 return result;
}
RVec<float> cfStages(const RVec<int>& category,const RVec<float>& pt,
 const RVec<int>& bits,const RVec<int>& pass,int region,int threshold) {
 RVec<float> stages;
 for(size_t i=0;i<pt.size();++i) if(category[i]==region && pt[i]>threshold) {
  if(((bits[i]&4095)==4095)!=bool(pass[i])) throw std::runtime_error("VID flag/bitmap mismatch");
  stages.push_back(.5f);
  for(int c=0;c<12;++c) {
   int mask=(1<<(c+1))-1;
   if((bits[i]&mask)!=mask) break;
   stages.push_back(c+1.5f);
  }
 }
 return stages;
}
''')
samples=['DY_Jets','DY_Jets_ext1','DY_Jets_ext2','W_Jets']
labels=['Baseline','p_{T}','SC #eta','#Delta#eta seed','#Delta#phi','#sigma_{i#etai#eta}',
        'E ratios','H/E','Track iso','ECAL+HCAL iso','d_{xy}','Missing hits','ECAL driven']
if not args.plots_only:
 for sample in samples:
  manifest=base/'out/prompt_heep_100files'/f'{sample}_files.txt'
  paths=manifest.read_text().splitlines()
  (out/f'{sample}_files.txt').write_text(manifest.read_text())
  print(f'{sample}: {len(paths)} files',flush=True)
  files=ROOT.std.vector('string')()
  for path in paths: files.push_back(path)
  df=ROOT.RDataFrame('resolvedNtuplizer/Events',files).Define('cf_category',
    'cfCategory(electron_ptRaw,electron_etaRaw,electron_phiRaw,electron_etaSC,gen_eta,gen_phi,gen_pdgId)')
  jobs=[];n_events=df.Count()
  for pop,sign in [('prompt',1),('unmatched',-1)]:
   for region,r in [('EB',1),('EE',2)]:
    for pt in [20,35]:
     for ident,bits,flag in [('nominal','electron_heepV70Bitmap','electron_passHeepV70'),
                             ('modified','electron_modHeepBitmap','electron_passModHeep')]:
      col=f'{pop}_{region}_pt{pt}_{ident}'
      node=df.Define(col,f'cfStages(cf_category,electron_ptRaw,{bits},{flag},{sign*r},{pt})')
      jobs.append(node.Histo1D((col,'',13,0,13),col))
  events=n_events.GetValue()
  f=ROOT.TFile(str(out/f'{sample}.root'),'RECREATE')
  for h in jobs: h.GetValue().Write()
  ROOT.TNamed('events',str(events)).Write();f.Close()
  print(f'{sample}: {events} events complete',flush=True)

hist={}
for sample in samples:
 f=ROOT.TFile.Open(str(out/f'{sample}.root'))
 if not f or f.IsZombie(): raise RuntimeError(f'Missing {sample} cutflows')
 group='DY' if sample.startswith('DY') else 'W'
 for key in f.GetListOfKeys():
  h=key.ReadObj()
  if not h.InheritsFrom('TH1'): continue
  name=key.GetName()
  if (group,name) not in hist: hist[group,name]=h.Clone(group+'_'+name)
  else: hist[group,name].Add(h)
 f.Close()

# Cross-check independent cumulative implementation against earlier unmatched counters.
old=list(csv.DictReader((base/'out/heep_unmatched_dr0p1_100files/cuts.csv').open()))
for r in old:
 for ident in ['nominal','modified']:
  h=hist[r['sample'],f'unmatched_{r["region"]}_pt{r["pt_min"]}_{ident}']
  assert h.GetBinContent(1)==int(r['selected'])
  assert h.GetBinContent(int(r['bit'])+2)==int(r[ident+'_cumulative'])

rows=[];keep=[]
for group in ['DY','W']:
 for pt in [20,35]:
  canvas=ROOT.TCanvas(f'{group}{pt}','',1600,1100);canvas.Divide(2,2)
  for row,pop in enumerate(['prompt','unmatched']):
   for col,region in enumerate(['EB','EE']):
    pad=canvas.cd(row*2+col+1);pad.SetBottomMargin(.27);pad.SetLeftMargin(.12);pad.SetTopMargin(.12)
    legend=ROOT.TLegend(.15,.68,.67,.81);legend.SetFillStyle(0);legend.SetBorderSize(0);legend.SetTextSize(.032)
    den=None
    for i,ident in enumerate(['nominal','modified']):
     raw=hist[group,f'{pop}_{region}_pt{pt}_{ident}']
     n=raw.GetBinContent(1)
     if den is not None: assert den==n
     den=n
     vals=[int(raw.GetBinContent(b)) for b in range(1,14)]
     assert all(a>=b for a,b in zip(vals,vals[1:]))
     for stage,count in enumerate(vals):
      rows.append([group,pop,region,pt,ident,stage,labels[stage],int(n),count,
                   count/n if n else '',(vals[stage-1]-count) if stage else 0])
     h=raw.Clone('draw_'+raw.GetName())
     if n: h.Scale(1/n)
     title=('Prompt-matched (#DeltaR<0.1)' if pop=='prompt' else 'Unmatched (#DeltaR#geq0.1)')
     h.SetTitle(f'{group}+jets | {title} | {region}, p_{{T}}>{pt} GeV')
     h.GetYaxis().SetTitle('Cumulative passing fraction');h.GetYaxis().SetTitleOffset(1.25)
     for b,label in enumerate(labels,1): h.GetXaxis().SetBinLabel(b,label)
     h.GetXaxis().LabelsOption('v');h.GetXaxis().SetLabelSize(.04)
     h.SetMinimum(0);h.SetMaximum(1.08)
     h.SetLineColor(ROOT.kBlue+1 if i==0 else ROOT.kRed+1);h.SetLineWidth(3)
     h.SetLineStyle(2 if i==0 else 1);h.SetMarkerStyle(24 if i==0 else 20)
     h.SetMarkerColor(h.GetLineColor());h.SetMarkerSize(.65)
     h.Draw('HIST' if i==0 else 'HIST SAME')
     legend.AddEntry(h,('Nominal HEEP V70' if i==0 else 'Modified HEEP')+f' ({vals[-1]:,}/{int(n):,})','l')
     keep.append(h)
    legend.Draw();keep.append(legend)
  for ext in ['png','pdf']: canvas.SaveAs(str(img/f'{group}_cutflow_pt{pt}.{ext}'))
  canvas.Close()
with (out/'cutflow.csv').open('w') as f:
 w=csv.writer(f);w.writerow(['sample','population','region','pt_min','ID','stage','cut','baseline','passing','fraction','lost_at_step']);w.writerows(rows)
print('PASS: unmatched counts reproduce previous study; all cutflows monotonic, paired baselines identical.',flush=True)
print(f'Plots: {img}',flush=True)
