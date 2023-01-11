#!/usr/bin/env python3

import ROOT
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
import hyperopt
import math
import argparse
import os

import xgbPreprocess as zprep
import xgbVis as zvis

parser = argparse.ArgumentParser(description="BDT training for merged electrons")
parser.add_argument("--det",dest="det",type=str,help="detector (EB/EE)")
parser.add_argument("--ang",dest="angle",type=str,help="opening angle (DR1/DR2/None)")
parser.add_argument("--et",dest="et",type=str,help="Et range (Et1/Et2)")
parser.add_argument("--opt",dest="opt",type=str,help="perform training")
parser.add_argument("--era",dest="era",type=str,help="era (20UL16APV/20UL16/20UL17/20UL18)")
parser.add_argument("--model",dest="model",type=str,help="path to model",default="")

args = parser.parse_args()

# choose the right signal model for relevant phase spaces
aFile = None
mergedTree = None
mergedGsfTree = None

if args.angle=="DR2" and args.et=="Et1":
    aFile = ROOT.TFile.Open("data/MergedEleMva_"+args.era+"_H200A1.root")
elif args.angle=="DR2" and args.et=="Et2":
    aFile = ROOT.TFile.Open("data/MergedEleMva_"+args.era+"_H800A10.root")
elif args.angle=="DR1" and args.et=="Et2":
    aFile = ROOT.TFile.Open("data/MergedEleMva_"+args.era+"_H800A1.root")
elif args.angle=="None" and args.et=="Et2":
    aFile = ROOT.TFile.Open("data/MergedEleMva_"+args.era+"_H800A1.root")
else:
    raise NameError('check dr/Et argument [DR1/DR2/None][Et1/Et2]')

if args.angle=="None":
    mergedTree = aFile.Get("mergedEleSigMvaInput/mergedEl2_elTree")
else:
    mergedTree = aFile.Get("mergedEleSigMvaInput/mergedEl1_elTree")
    mergedGsfTree = aFile.Get("mergedEleSigMvaInput/mergedEl1Gsf_addGsfTree")

mergedWgtHist= aFile.Get("mergedEleSigMvaInput/totWeightedSum")
mergedTotWgtSum = mergedWgtHist.GetBinContent(1)

workflowname = args.angle+args.et+'_'+args.det+'_'+args.era
dfProducer = zprep.DataframeInitializer(args.det,args.angle,args.et)

bkglist = [
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-15to20_EM.root",1324000.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-20to30_EM.root",4896000.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-30to50_EM.root",6447000.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-50to80_EM.root",1988000.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-80to120_EM.root",367500.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-120to170_EM.root",66590.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-170to300_EM.root",16620.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_QCD_Pt-300toInf_EM.root",1104.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_WJets.root",53870.0),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_DY.root",6077.22),
    zprep.SampleProcessor("data/MergedEleMva_"+args.era+"_TT.root",831.76)
]

for aproc in bkglist:
    aproc.read(dfProducer)

df_mergedEl, wgts_mergedEl, pts_mergedEl, etas_mergedEl, invMs_mergedEl = dfProducer.fill_arr(mergedTree,mergedGsfTree)
col_names = dfProducer.col_names()

# equalize sumwgt(bkg) to sumwgt(sig) to avoid class imbalance
sum_sig = np.sum(wgts_mergedEl)
sum_bkg = 0.

for aproc in bkglist:
    sum_bkg += np.sum(aproc.wgts_)

for aproc in bkglist:
    aproc.wgts_ = aproc.wgts_*(sum_sig/sum_bkg)

# concat bkg
df_bkg = pd.concat([ aproc.df_ for aproc in bkglist ], axis=0, ignore_index=True)
wgts_bkg = np.concatenate([ aproc.wgts_ for aproc in bkglist ], axis=0)
pts_bkg = np.concatenate([ aproc.pts_ for aproc in bkglist ], axis=0)

# prepare labels
y_mergedEl = np.ones(shape=(df_mergedEl.shape[0],))
y_bkgEl = np.zeros(shape=(df_bkg.shape[0],))

# concat sig+bkg
df_total = pd.concat([df_mergedEl,df_bkg], axis=0, ignore_index=True)
y_total = np.concatenate((y_mergedEl,y_bkgEl), axis=0)
wgts_total = np.concatenate((wgts_mergedEl,wgts_bkg), axis=0)

transformer = sklearn.preprocessing.StandardScaler()
x_total = transformer.fit_transform(df_total)
np_meanstd = np.array([transformer.mean_,transformer.scale_])

print('total number of merged objects = %d' % y_mergedEl.shape[0])
print('total number of bkg objects = %d' % y_bkgEl.shape[0])

dtotal = xgb.DMatrix(x_total, weight=wgts_total, label=y_total, feature_names=col_names)

# automatic hyper-optimization
def objective(params,dmat):
    param = {
        'max_depth': int(params['max_depth']),
        'eta': params['eta'],
        'gamma': params['gamma'],
        'lambda': params['lambda'],
        'min_child_weight':params['min_child_weight'],
        'objective':'binary:logistic',
        'subsample':0.5,
        'eval_metric':'logloss',
        'tree_method':'gpu_hist',
        'nthread':4
    }

    xgb_cv = xgb.cv(dtrain=dmat,nfold=5,num_boost_round=200,metrics='logloss',early_stopping_rounds=20,params=param)

    return xgb_cv['test-logloss-mean'].min()

if args.opt=="True": # hyper-optimization
    weightSum = np.sum(wgts_total)

    param_space = {
        'max_depth': hyperopt.hp.quniform('max_depth',1,5,1), # lower depth, more robust # 4 for EE dr2 dr3
        'eta': hyperopt.hp.loguniform('eta',-3.,-1.), # from exp(-3) to exp(1)
        'gamma': hyperopt.hp.uniform('gamma',0.,10.),
        'lambda': hyperopt.hp.uniform('lambda',0.,2.),
        'min_child_weight': hyperopt.hp.loguniform('min_child_weight',math.log(weightSum/1000.),math.log(weightSum/10.))
    }

    trials = hyperopt.Trials()
    objec = lambda x: objective(x, dtotal)
    best = hyperopt.fmin(fn=objec, space=param_space, max_evals=256, algo=hyperopt.tpe.suggest, trials=trials)

    print('the best params are')
    print(best)

    # save hyper parameters
    np.savez('npz/'+workflowname+'.npz', meanstd=np_meanstd, bestparams=np.array(best))

# load hyper parameters
npz_params = np.array([])

if args.model is not "" and os.path.splitext(args.model)[1]=='.model':
    npz_params = np.load('npz/'+os.path.splitext(os.path.basename(args.model))[0]+'.npz',allow_pickle=True)
elif args.opt=="True":
    npz_params = np.load('npz/'+workflowname+'.npz',allow_pickle=True)
else:
    raise NameError('Neither loading existing model nor running the training, please check opt/model arguments')

param = npz_params['bestparams'][()] # numpy stores dict as 0-d array
param['max_depth'] = int(param['max_depth'])
param['objective'] = 'binary:logistic'
param['subsample'] = 0.5
param['eval_metric'] = 'logloss'
param['tree_method'] = 'exact'
param['nthread'] = 4

# run k-fold training, then pick the model with the smallest overtraining i.e. minimum |AUC(train)-AUC(test)|
kf = sklearn.model_selection.KFold(n_splits=5,shuffle=True)
modelPerforms_list = []
model_list = []

if args.opt=='True':
    for train, test in kf.split(x_total):
        x_train = x_total[train]
        y_train = y_total[train]
        wgts_train = wgts_total[train]

        x_test = x_total[test]
        y_test = y_total[test]
        wgts_test = wgts_total[test]

        dtrain = xgb.DMatrix(x_train, weight=wgts_train, label=y_train, feature_names=col_names)
        dtest = xgb.DMatrix(x_test, weight=wgts_test, label=y_test, feature_names=col_names)

        evallist = [(dtest, 'eval'), (dtrain, 'train')]
        num_round = 200

        bst = xgb.Booster(param)
        early_stop = xgb.callback.EarlyStopping(rounds=20,metric_name='logloss',data_name='eval')

        if args.opt=='True':
            bst = xgb.train(param, dtrain, num_round, evallist, callbacks=[early_stop])
            model_list.append(bst)

        dTrainPredict    = bst.predict(dtrain)
        dTestPredict     = bst.predict(dtest)

        modelPerforms_list.append( zvis.calROC(dTrainPredict,dTestPredict,y_train,y_test) )

targetFpr = 0.2
nbins = 50
idx_max = -1
targetThres = -1
trainThres = -1
bst = None
plotname = workflowname

# calculate score threshold at a certain FPR, which differs by w/ or w/o GSF scenarios
if args.angle=="None":
    targetFpr = 0.05
    nbins = 100
elif args.angle=="DR2" and args.et=="Et2" and args.det=="EB":
    targetFpr = 0.1
elif args.angle=="DR2" and args.et=="Et2" and args.det=="EE":
    targetFpr = 0.1
else:
    pass

# append the model name when not training
if args.model is not "" and os.path.splitext(args.model)[1]=='.model':
    plotname += ('_'+os.path.splitext(os.path.basename(args.model))[0])

if args.opt=='True':
    # draw ROC curve and calculate score threshold
    idx_max = zvis.drawROC(modelPerforms_list, plotname+'_roc')
    targetThres, trainThres = zvis.drawThr(modelPerforms_list[idx_max], targetFpr, plotname+'_thr')

    # use train threshold instead when test set lacks of statistics
    if args.angle=="DR2" and args.et=="Et2" and args.det=="EB":
        targetThres = trainThres
    elif args.angle=="DR2" and args.et=="Et2" and args.det=="EE":
        targetThres = trainThres
    else:
        pass

    # save model & importance plot
    bst = model_list[idx_max]
    bst.save_model('model/'+workflowname+'.model')
    gain = bst.get_score(importance_type='gain')
    cover = bst.get_score(importance_type='cover')
    zvis.drawImportance(gain,cover,col_names,plotname+'_importance')
else:
    bst = xgb.Booster(param)
    targetThres = 0.791 # set manually for now
    bst.load_model('model/'+workflowname+'.model')

# calculate scores by physics processes
def bkg_predict(abst,atransformer,adf,awgt,colname):
    if adf.shape[0]==0:
        return np.zeros(shape=(adf.shape[0],))

    atrans= atransformer.transform(adf)
    admat = xgb.DMatrix(atrans, weight=awgt, label=np.zeros(shape=(adf.shape[0],)), feature_names=colname)

    return abst.predict(admat)

np_meanstd = npz_params['meanstd'][()]
transformer.mean_ = np_meanstd[0]
transformer.scale_ = np_meanstd[1]
trans_mergedEl = transformer.transform(df_mergedEl)
dSig = xgb.DMatrix(trans_mergedEl, weight=wgts_mergedEl, label=y_mergedEl, feature_names=col_names)
dSigPredict = bst.predict(dSig)

score_list = []

for aproc in bkglist:
    score_list.append( bkg_predict(bst,transformer,aproc.df_,aproc.wgts_,col_names) )

scores_bkg = np.concatenate(score_list, axis=0)

zvis.drawScoreByProcess(
    dSigPredict,
    wgts_mergedEl,
    list(np.array([aproc for aproc in reversed(score_list)],dtype=object)),
    list(np.array([aproc.wgts_ for aproc in reversed(bkglist)],dtype=object)),
    ['TT','DY','WJets','_','_','_','QCD','_','_','_','_'],
    ['tomato','wheat','green','steelblue','cadetblue','darkcyan','darkturquoise','cyan','turquoise','paleturquoise','lightcyan'],
    nbins,
    plotname+'_scoreProcess'
)

zvis.drawScoreByProcess(
    None,
    None,
    # hail python
    list(np.array([aproc.etas_[score_list[idx] > targetThres] for idx, aproc in reversed(list(enumerate(bkglist)))],dtype=object)),
    list(np.array([aproc.wgts_[score_list[idx] > targetThres] for idx, aproc in reversed(list(enumerate(bkglist)))],dtype=object)),
    ['TT','DY','WJets','_','_','_','QCD','_','_','_','_'],
    ['tomato','wheat','green','steelblue','cadetblue','darkcyan','darkturquoise','cyan','turquoise','paleturquoise','lightcyan'],
    100,
    plotname+'_etaSC'
)

if args.angle!="None":
    zvis.drawScoreByProcess(
        None,
        None,
        list(np.array([aproc.invMs_[score_list[idx] > targetThres] for idx, aproc in reversed(list(enumerate(bkglist)))],dtype=object)),
        list(np.array([aproc.wgts_[score_list[idx] > targetThres] for idx, aproc in reversed(list(enumerate(bkglist)))],dtype=object)),
        ['TT','DY','WJets','_','_','_','QCD','_','_','_','_'],
        ['tomato','wheat','green','steelblue','cadetblue','darkcyan','darkturquoise','cyan','turquoise','paleturquoise','lightcyan'],
        100,
        plotname+'_invM'
    )

# efficiency as a function of Et at threshold
etThresLow = 0.
etThresHigh = 1000.

if args.et=="Et1":
    if args.det=="EB":
        etThresHigh = 200.
    else:
        etThresHigh = 150.
elif args.et=="Et2":
    if args.det=="EB":
        etThresLow = 200.
    else:
        etThresLow = 150.
else:
    raise NameError('check Et argument [Et1/Et2]')

effplot_sig = ROOT.TEfficiency("eff_sig",";GeV;Eff",20,etThresLow,etThresHigh)
effplot_bkg = ROOT.TEfficiency("eff_bkg",";GeV;Eff",20,etThresLow,etThresHigh)

for idx in range(len(wgts_mergedEl)):
    effplot_sig.FillWeighted(dSigPredict[idx] > targetThres, wgts_mergedEl[idx], pts_mergedEl[idx])

for idx in range(len(wgts_bkg)):
    effplot_bkg.FillWeighted(scores_bkg[idx] > targetThres, wgts_bkg[idx], pts_bkg[idx])

# some serious plotting
import CMS_lumi, tdrstyle

#set the tdr style
tdrstyle.setTDRStyle()

#change the CMS_lumi variables (see CMS_lumi.py)
CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Work in progress"
CMS_lumi.lumi_sqrtS = "13 TeV" # used with iPeriod = 0

iPos = 0

if( iPos==0 ):
    CMS_lumi.relPosX = 0.12

H_ref = 600;
W_ref = 800;
W = W_ref
H = H_ref

iPeriod = 0

# references for T, B, L, R
T = 0.08*H_ref
B = 0.12*H_ref
L = 0.12*W_ref
R = 0.04*W_ref

canvas = ROOT.TCanvas("c2","c2",50,50,W,H)
canvas.SetFillColor(0)
canvas.SetBorderMode(0)
canvas.SetFrameFillStyle(0)
canvas.SetFrameBorderMode(0)
canvas.SetLeftMargin( L/W )
canvas.SetRightMargin( R/W )
canvas.SetTopMargin( T/H )
canvas.SetBottomMargin( B/H )
canvas.SetTickx(0)
canvas.SetTicky(0)

effplot_bkg.SetLineWidth(2)
effplot_sig.SetLineWidth(2)
effplot_bkg.SetLineColor(ROOT.kRed)
effplot_sig.SetLineColor(ROOT.kBlue)
effplot_bkg.Draw()
canvas.Update()
agraph = effplot_bkg.GetPaintedGraph()
agraph.GetYaxis().SetTitleOffset(1)
agraph.SetMinimum(0.);
agraph.SetMaximum(1.2);
canvas.Update()
effplot_sig.Draw("sames")

legend = ROOT.TLegend(0.85,0.8,0.98,0.9)
legend.SetBorderSize(0);
legend.AddEntry(effplot_sig,"sig")
legend.AddEntry(effplot_bkg,"bkg")
legend.Draw()

#draw the lumi text on the canvas
CMS_lumi.CMS_lumi(canvas, iPeriod, iPos)
canvas.cd()
canvas.Update()
canvas.RedrawAxis()
frame = canvas.GetFrame()
frame.Draw()

canvas.SaveAs("plot/"+plotname+"_eff.png")
