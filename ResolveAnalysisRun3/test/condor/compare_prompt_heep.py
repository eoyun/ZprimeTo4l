#!/usr/bin/env python3
"""DY/W prompt gen-matched electron diagnostic; source env.sh first.

--max-files 100 selects evenly spaced files per production; 0 uses all files.
This is not a PF-isolated sample. Uses stored VID bits, not reimplemented cuts.
"""
import argparse
import csv
from pathlib import Path
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.TH1.AddDirectory(False)
BASE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--max-files', type=int, default=100)
parser.add_argument('--plots-only', action='store_true')
parser.add_argument('--population', choices=['matched','unmatched'], default='matched')
parser.add_argument('--match-dr', type=float, default=None)
args = parser.parse_args()
if args.max_files < 0: parser.error('--max-files must be nonnegative')
match_dr = args.match_dr if args.match_dr is not None else (.1 if args.population=='unmatched' else .001)
if match_dr <= 0: parser.error('--match-dr must be positive')
tag = 'all' if args.max_files == 0 else f'{args.max_files}files'
prefix = 'prompt_heep' if args.population=='matched' and match_dr==.001 else f'heep_{args.population}_dr{match_dr:g}'.replace('.','p')
OUT = BASE / 'out' / (prefix+'_'+tag)
IMG = BASE / 'imgs' / (prefix+'_'+tag)
population_label = 'gen matched' if args.population=='matched' else f'gen unmatched (#DeltaR#geq{match_dr:g})'
OUT.mkdir(parents=True, exist_ok=True)
IMG.mkdir(parents=True, exist_ok=True)
ROOT.gInterpreter.Declare(r'''
#include <ROOT/RVec.hxx>
#include <cmath>
#include <algorithm>
#include <stdexcept>
using ROOT::VecOps::RVec;
RVec<int> promptRegion(const RVec<float>& ep,const RVec<float>& ee,
 const RVec<float>& ef,const RVec<float>& sc,const RVec<float>& ge,
 const RVec<float>& gf,const RVec<int>& pid,double drCut,bool unmatched) {
  RVec<int> result(ep.size(),0);
  for(size_t i=0;i<ep.size();++i) {
    if(ep[i]<=20) continue;
    double best=1e9;
    for(size_t j=0;j<pid.size();++j) if(std::abs(pid[j])==11 && std::abs(ge[j])<2.5)
      best=std::min(best,std::hypot(double(ee[i]-ge[j]),std::remainder(double(ef[i]-gf[j]),2*M_PI)));
    if(unmatched ? best<drCut : best>=drCut) continue;
    double eta=std::abs(sc[i]);
    if(eta<1.4442) result[i]=1;
    else if(eta>1.566 && eta<2.5) result[i]=2;
  }
  return result;
}
bool promptValid(float x) {return std::isfinite(x) && x!=-999.f && std::abs(x)<1e30;}
RVec<float> ratePts(const RVec<float>& pt,const RVec<int>& reg,
 const RVec<int>& nom,const RVec<int>& mod,int r,int group) {
  RVec<float> out;
  for(size_t i=0;i<pt.size();++i) if(reg[i]==r && pt[i]>20) {
    bool n=(nom[i]&4095)==4095,m=(mod[i]&4095)==4095;
    if(group==0 || (group==1 && n) || (group==2 && m) || (group==3 && !n && m)) out.push_back(pt[i]);
  }
  return out;
}
RVec<float> promptValues(const RVec<float>& x,const RVec<int>& reg,
 const RVec<float>& pt,const RVec<int>& nom,const RVec<int>& mod,int r,int group) {
  if(x.size()!=reg.size()) throw std::runtime_error("electron vector length mismatch");
  RVec<float> out;
  for(size_t i=0;i<x.size();++i) if(reg[i]==r && pt[i]>35 && promptValid(x[i])) {
    bool n=(nom[i]&4095)==4095,m=(mod[i]&4095)==4095;
    if(group==0 || (group==1 && n) || (group==2 && m) ||
      (group==3 && !(nom[i]&4) && (mod[i]&4))) out.push_back(x[i]);
  }
  return out;
}
RVec<float> promptMatching(const RVec<float>& eta,const RVec<float>& perp,
 const RVec<int>& reg,const RVec<float>& pt,int r,bool usePerp) {
  RVec<float> out;
  for(size_t i=0;i<eta.size();++i) if(reg[i]==r && pt[i]>35 && promptValid(eta[i]) && promptValid(perp[i]))
    out.push_back(std::abs(usePerp ? perp[i] : eta[i]));
  return out;
}
RVec<float> promptStats(const RVec<int>& reg,const RVec<float>& pt,
 const RVec<int>& nom,const RVec<int>& mod,const RVec<int>& passn,
 const RVec<int>& passm,const RVec<float>& perp,int r,int threshold) {
  RVec<float> out;
  for(size_t i=0;i<reg.size();++i) if(reg[i]==r && pt[i]>threshold) {
    int n=nom[i]&4095,m=mod[i]&4095;
    auto add=[&](int k) {out.push_back(k+.5f);};
    add(0);
    if(n==4095) add(1);
    if(m==4095) add(2);
    if(n!=4095 && m==4095) add(3);
    if(n==4095 && m!=4095) add(4);
    if(promptValid(perp[i])) add(5);
    if((n==4095)!=bool(passn[i]) || (m==4095)!=bool(passm[i])) add(6);
    for(int c=0;c<12;++c) {
      int b=1<<c,seq=(1<<(c+1))-1,other=4095^b;
      if(n&b) add(10+c);
      if(m&b) add(22+c);
      if((n&seq)==seq) add(34+c);
      if((m&seq)==seq) add(46+c);
      if((n&other)==other) add(58+c);
      if((m&other)==other) add(70+c);
      if(n==other) add(82+c);
      if(m==other) add(94+c);
      if(!(n&b) && (m&b)) add(106+c);
      if((n&b) && !(m&b)) add(118+c);
    }
  }
  return out;
}
''')
cuts = ['pT','SC eta','dEta seed','dPhi','sigma ieta','E ratios','H/E','track iso','ECAL+HCAL iso','dxy','missing hits','ECAL driven']
variables = {
 'dEtaInSeed':('electron_dEtaInSeed','#Delta#eta_{in}^{seed} (signed)',-.03,.03),
 'dPerpIn':('electron_dPerpIn','#Delta u_{5#times5}^{in} (signed)',-.03,.03),
 'sigmaIetaIeta':('electron_sigmaIeta','full5x5 #sigma_{i#etai#eta}',0.,.10),
 'dPhiIn':('electron_dPhiIn','#Delta#phi_{in}',-.15,.15),
 'hOverE':('electron_hOverE','H/E',0.,.20),
 'modTrkIso':('electron_modTrkIso','Modified track isolation [GeV]',0.,20.),
}
components = ['DY_Jets','DY_Jets_ext1','DY_Jets_ext2','W_Jets']
audit=[]
if not args.plots_only:
 for sample in components:
    paths=(BASE/'filelists'/f'{sample}.txt').read_text().splitlines()
    available=len(paths)
    if args.max_files and len(paths)>args.max_files:
        indices=[round(i*(len(paths)-1)/max(1,args.max_files-1)) for i in range(args.max_files)]
        paths=[paths[i] for i in indices]
    (OUT/f'{sample}_files.txt').write_text('\n'.join(paths)+'\n')
    print(f'{sample}: checking {len(paths)}/{available} files',flush=True)
    required=[v[0] for v in variables.values()]+['electron_ptRaw','electron_etaRaw','electron_phiRaw',
        'electron_etaSC','gen_eta','gen_phi','gen_pdgId','electron_heepV70Bitmap',
        'electron_modHeepBitmap','electron_passHeepV70','electron_passModHeep']
    entries=0
    for path in paths:
        f=ROOT.TFile.Open(path)
        if not f or f.IsZombie(): raise RuntimeError(f'Unreadable file {path}')
        t=f.Get('resolvedNtuplizer/Events')
        if not t or any(not t.GetBranch(n) for n in required):
            raise RuntimeError(f'Missing required tree/branch in {path}')
        entries+=t.GetEntries(); f.Close()
    print(f'{sample}: processing {entries} events',flush=True)
    files=ROOT.std.vector('string')()
    for path in paths: files.push_back(path)
    df=ROOT.RDataFrame('resolvedNtuplizer/Events',files).Define('prompt_region',
      f'promptRegion(electron_ptRaw,electron_etaRaw,electron_phiRaw,electron_etaSC,gen_eta,gen_phi,gen_pdgId,{match_dr},{str(args.population=="unmatched").lower()})')
    total=df.Count(); jobs=[]
    for r,region in enumerate(['EB','EE'],1):
        from array import array
        pt_edges=array('d',[20,25,30,35,40,50,70,100,150,200,300,500,1000])
        for g in range(4):
            col=f'ratePt_{region}_g{g}'
            node=df.Define(col,f'ratePts(electron_ptRaw,prompt_region,electron_heepV70Bitmap,electron_modHeepBitmap,{r},{g})')
            jobs.append(node.Histo1D((col,';Electron raw p_{T} [GeV];Electrons',len(pt_edges)-1,pt_edges),col))
        for threshold in [20,35]:
            col=f'stats_{region}_pt{threshold}'
            node=df.Define(col,f'promptStats(prompt_region,electron_ptRaw,electron_heepV70Bitmap,electron_modHeepBitmap,electron_passHeepV70,electron_passModHeep,electron_dPerpIn,{r},{threshold})')
            jobs.append(node.Histo1D((col,'',130,0,130),col))
        for name,(branch,title,lo,hi) in variables.items():
            if name=='sigmaIetaIeta' and r==1: hi=.03
            for g in range(4):
                col=f'{name}_{region}_g{g}'
                node=df.Define(col,f'promptValues({branch},prompt_region,electron_ptRaw,electron_heepV70Bitmap,electron_modHeepBitmap,{r},{g})')
                jobs.append(node.Histo1D((col,';'+title+';Electrons',120,lo,hi),col))
        for usePerp,name in [(False,'absEta_sameSubset'),(True,'absPerp_sameSubset')]:
            col=f'{name}_{region}'
            node=df.Define(col,f'promptMatching(electron_dEtaInSeed,electron_dPerpIn,prompt_region,electron_ptRaw,{r},{str(usePerp).lower()})')
            jobs.append(node.Histo1D((col,';Absolute matching residual;Electrons',100,0,.03),col))
    actual=total.GetValue()
    if actual!=entries: raise RuntimeError('Event count mismatch')
    output=ROOT.TFile(str(OUT/f'{sample}.root'),'RECREATE')
    for h in jobs: h.GetValue().Write()
    output.Close()
    audit.append([sample,available,len(paths),actual])
    print(f'{sample}: complete',flush=True)
 with (OUT/'input_audit.csv').open('w') as f:
    writer=csv.writer(f);writer.writerow(['sample','available_files','processed_files','events']);writer.writerows(audit)

# Sum the DY productions with unit electron weights; these are diagnostic
# efficiencies in the processed sample, not cross-section weighted predictions.
histograms={}
for sample in components:
    f=ROOT.TFile.Open(str(OUT/f'{sample}.root'))
    if not f or f.IsZombie(): raise RuntimeError(f'Missing output for {sample}')
    group='DY' if sample.startswith('DY') else 'W'
    for key in f.GetListOfKeys():
        name=key.GetName();h=key.ReadObj()
        if (group,name) not in histograms: histograms[group,name]=h.Clone(group+'_'+name)
        else: histograms[group,name].Add(h)
    f.Close()

def write_csv(name,header,rows):
    with (OUT/name).open('w') as f:
        writer=csv.writer(f);writer.writerow(header);writer.writerows(rows)

cutrows=[];summary=[];flows=[]
for group in ['DY','W']:
 for region in ['EB','EE']:
  for pt in [20,35]:
    h=histograms[group,f'stats_{region}_pt{pt}']
    v=lambda k:int(h.GetBinContent(k+1))
    if v(6): raise RuntimeError('Stored ID boolean and VID bitmap disagree')
    assert v(2)-v(1)==v(3)-v(4)
    summary.append([group,region,pt]+[v(k) for k in range(7)])
    for c,label in enumerate(cuts):
        assert v(22+c)-v(10+c)==v(106+c)-v(118+c)
        cutrows.append([group,region,pt,c,label,v(0)]+[v(offset+c) for offset in [10,22,34,46,58,70,82,94,106,118]])
 for (g,name),h in histograms.items():
    if g==group and not name.startswith('stats'):
        flows.append([g,name,int(h.GetEntries()),h.GetBinContent(0),h.GetBinContent(h.GetNbinsX()+1)])
write_csv('summary.csv',['sample','region','pt_min','selected','nominal_pass','modified_pass','recovered','lost','valid_dPerpIn','bitmap_mismatch'],summary)
write_csv('cuts.csv',['sample','region','pt_min','bit','cut','selected','nominal_individual','modified_individual','nominal_cumulative','modified_cumulative','nominal_Nminus1_den','modified_Nminus1_den','nominal_only_failure','modified_only_failure','cut_recovered','cut_lost'],cutrows)
write_csv('flows.csv',['sample','histogram','valid_entries','underflow','overflow'],flows)

keep=[]
def save(c,name):
    for ext in ['png','pdf']: c.SaveAs(str(IMG/(name+'.'+ext)))
    c.Close()

for group in ['DY','W']:
 for pt in [20,35]:
  c=ROOT.TCanvas('cutimpact','cutimpact',1800,950);c.Divide(3,2)
  for row,region in enumerate(['EB','EE']):
    source=histograms[group,f'stats_{region}_pt{pt}']
    den=source.GetBinContent(1)
    for col,(mode,offsets) in enumerate([('Individual cut pass',[10,22]),('Cumulative pass',[34,46]),('Failure after other 11 cuts',[82,94])]):
        pad=c.cd(row*3+col+1);pad.SetBottomMargin(.26);pad.SetLeftMargin(.13)
        legend=ROOT.TLegend(.15,.79,.56,.89);legend.SetBorderSize(0);legend.SetFillStyle(0)
        if col != 2:
            legend.SetY1NDC(.20);legend.SetY2NDC(.30)
        rates=[]
        for j,offset in enumerate(offsets):
            for bit in range(12):
                d=source.GetBinContent((58 if j==0 else 70)+bit+1) if col==2 else den
                if d: rates.append(source.GetBinContent(offset+bit+1)/d)
        for j,offset in enumerate(offsets):
            h=ROOT.TH1D(f'{group}_{region}_{pt}_{col}_{j}',f'{group}+jets | {region} | p_{{T}}>{pt} GeV',12,0,12)
            for bit,label in enumerate(cuts):
                n=source.GetBinContent(offset+bit+1)
                d=source.GetBinContent((58 if j==0 else 70)+bit+1) if col==2 else den
                h.SetBinContent(bit+1,n/d if d else 0)
                h.GetXaxis().SetBinLabel(bit+1,label)
            h.GetXaxis().LabelsOption('v');h.GetYaxis().SetTitle(mode+' fraction')
            h.SetMinimum(0 if col==2 else max(0,min(rates+[1])-.08))
            h.SetMaximum(max(rates+[.001])*1.4 if col==2 else 1.02);h.SetLineWidth(2)
            h.SetLineColor(ROOT.kBlue+1 if j==0 else ROOT.kRed+1)
            h.Draw('HIST' if j==0 else 'HIST SAME')
            legend.AddEntry(h,'Nominal HEEP V70' if j==0 else 'Modified HEEP','l');keep.append(h)
        legend.Draw();keep.append(legend)
  save(c,f'{group}_cut_impact_pt{pt}')
 for name in variables:
  c=ROOT.TCanvas('inputs','inputs',1400,650);c.Divide(2,1)
  for row,region in enumerate(['EB','EE']):
    pad=c.cd(row+1);pad.SetLogy();pad.SetLeftMargin(.15);pad.SetBottomMargin(.14)
    legend=ROOT.TLegend(.51,.67,.95,.89);legend.SetBorderSize(0);legend.SetFillStyle(0);legend.SetTextSize(.029)
    labels=['All '+args.population,'Pass nominal','Pass modified','dEta cut recovered']
    colors=[ROOT.kBlack,ROOT.kBlue+1,ROOT.kRed+1,ROOT.kGreen+2]
    for g in range(4):
        raw=histograms[group,f'{name}_{region}_g{g}'];h=raw.Clone('draw_'+raw.GetName())
        n=h.Integral(0,h.GetNbinsX()+1)
        if n: h.Scale(1/n)
        h.SetTitle(f'{group}+jets | {region} | p_{{T}}>35 GeV, {population_label}')
        h.SetLineColor(colors[g]);h.SetLineWidth(2);h.SetMinimum(1e-5);h.SetMaximum(2)
        h.GetYaxis().SetTitle('Fraction of valid electrons / bin');h.GetYaxis().SetTitleOffset(1.6)
        h.Draw('HIST' if g==0 else 'HIST SAME')
        legend.AddEntry(h,f'{labels[g]} (N={int(raw.GetEntries())})','l');keep.append(h)
    legend.Draw();keep.append(legend)
  save(c,f'{group}_{name}')
 c=ROOT.TCanvas('same_subset','same_subset',1400,650);c.Divide(2,1)
 for row,region in enumerate(['EB','EE']):
    pad=c.cd(row+1);pad.SetLogy();pad.SetLeftMargin(.15);pad.SetBottomMargin(.14)
    legend=ROOT.TLegend(.50,.74,.94,.89);legend.SetBorderSize(0);legend.SetFillStyle(0)
    for j,name in enumerate(['absEta_sameSubset','absPerp_sameSubset']):
        raw=histograms[group,f'{name}_{region}'];h=raw.Clone('draw_'+raw.GetName())
        n=h.Integral(0,h.GetNbinsX()+1)
        if n: h.Scale(1/n)
        h.SetTitle(f'{group}+jets | {region} | p_{{T}}>35, both variables valid')
        h.SetLineColor(ROOT.kBlue+1 if j==0 else ROOT.kRed+1);h.SetLineWidth(2)
        h.SetMinimum(1e-5);h.SetMaximum(2);h.GetYaxis().SetTitle('Fraction / bin')
        h.Draw('HIST' if j==0 else 'HIST SAME');keep.append(h)
        legend.AddEntry(h,'|#Delta#eta_{in}^{seed}|' if j==0 else '|#Delta u_{5#times5}^{in}|','l')
    legend.Draw();keep.append(legend)
 save(c,f'{group}_matching_same_subset')
# A zoom on the actual changes: small differences can be invisible in an
# efficiency plot whose full vertical scale is 0..1.
for group in ['DY','W']:
 c=ROOT.TCanvas('transitions','transitions',1400,650);c.Divide(2,1)
 for row,region in enumerate(['EB','EE']):
    pad=c.cd(row+1);pad.SetBottomMargin(.25);pad.SetLeftMargin(.15)
    source=histograms[group,f'stats_{region}_pt35']
    denominator=source.GetBinContent(1)
    legend=ROOT.TLegend(.17,.77,.87,.89);legend.SetBorderSize(0);legend.SetFillStyle(0)
    maximum=max([source.GetBinContent(o+b+1) for o in [106,118] for b in range(12)]+[1])
    for j,offset in enumerate([106,118]):
        h=ROOT.TH1D(f'changes_{group}_{region}_{j}',f'{group}+jets | {region} | p_{{T}}>35 GeV',12,0,12)
        for bit,label in enumerate(cuts):
            h.SetBinContent(bit+1,100*source.GetBinContent(offset+bit+1)/denominator if denominator else 0)
            h.GetXaxis().SetBinLabel(bit+1,label)
        h.GetXaxis().LabelsOption('v');h.GetYaxis().SetTitle('Changed electrons / denominator [%]')
        h.SetMinimum(0);h.SetMaximum(100*maximum/max(1,denominator)*1.4)
        h.SetFillColor(ROOT.kBlue-7 if j==0 else ROOT.kRed-7)
        h.SetBarWidth(.4);h.SetBarOffset(.05 if j==0 else .5)
        h.Draw('BAR' if j==0 else 'BAR SAME');keep.append(h)
        legend.AddEntry(h,'Nominal fail #rightarrow modified pass' if j==0 else 'Nominal pass #rightarrow modified fail','f')
    legend.Draw();keep.append(legend)
 save(c,f'{group}_cut_changes_pt35')
rate_rows=[]
for row in summary:
    group,region,pt,den,nom,mod,*_=row
    for label,num in [('nominal',nom),('modified',mod)]:
        rate=num/den if den else ''
        lo=ROOT.TEfficiency.ClopperPearson(den,num,.682689492,False) if den else ''
        hi=ROOT.TEfficiency.ClopperPearson(den,num,.682689492,True) if den else ''
        rate_rows.append([group,region,pt,label,den,num,rate,lo,hi])
write_csv('rates.csv',['sample','region','pt_min','ID','denominator','numerator','rate','lower_68pct','upper_68pct'],rate_rows)
for group in ['DY','W']:
 if (group,'ratePt_EB_g0') not in histograms: continue  # old cached outputs
 c=ROOT.TCanvas('rates_pt','rates_pt',1400,650);c.Divide(2,1)
 for row,region in enumerate(['EB','EE']):
    pad=c.cd(row+1);pad.SetLogx();pad.SetLeftMargin(.15);pad.SetBottomMargin(.14)
    den=histograms[group,f'ratePt_{region}_g0']
    graphs=[]
    for j in [1,2]:
        num=histograms[group,f'ratePt_{region}_g{j}']
        if not ROOT.TEfficiency.CheckConsistency(num,den): raise RuntimeError('Invalid efficiency counts')
        graph=ROOT.TGraphAsymmErrors(num,den,'cl=0.682689492 cp')
        graph.SetLineColor(ROOT.kBlue+1 if j==1 else ROOT.kRed+1)
        graph.SetMarkerColor(graph.GetLineColor());graph.SetMarkerStyle(24 if j==1 else 20)
        graph.SetTitle(f'{group}+jets | {region} | {population_label};Electron raw p_{{T}} [GeV];ID passing fraction')
        graphs.append(graph)
    ymax=max([graphs[j].GetPointY(k)+graphs[j].GetErrorYhigh(k) for j in range(2) for k in range(graphs[j].GetN())]+[.01])
    graphs[0].SetMinimum(0);graphs[0].SetMaximum(min(1.05,ymax*1.3))
    graphs[0].Draw('AP');graphs[1].Draw('P SAME');graphs[0].GetXaxis().SetLimits(20,1000)
    legend=ROOT.TLegend(.49,.76,.94,.88);legend.SetBorderSize(0);legend.SetFillStyle(0)
    legend.AddEntry(graphs[0],'Nominal HEEP V70','lp');legend.AddEntry(graphs[1],'Modified HEEP','lp')
    legend.Draw();keep.extend(graphs+[legend])
 save(c,f'{group}_passing_rate_vs_pt')
print(f'Completed. Plots: {IMG}; counts: {OUT}',flush=True)
