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

import ZprimeTo4l.MergedLepton.xgbPreprocess as zprep
import ZprimeTo4l.MergedLepton.xgbVis as zvis

parser = argparse.ArgumentParser(description="BDT training for merged electrons")
parser.add_argument("-s","--sig",dest="sigFile",type=str,help="signal sample")
parser.add_argument("--dy",dest="dyFile",type=str,help="DY sample")
parser.add_argument("--tt",dest="ttFile",type=str,help="ttbar sample")
parser.add_argument("--det",dest="det",type=str,help="detector (EB/EE)")
parser.add_argument("--ang",dest="angle",type=str,help="opening angle (dr1/dr2/dr3)")
parser.add_argument("--opt",dest="opt",type=str,help="perform training")
parser.add_argument("--load",dest="load",type=str,help="load model",default="False")
parser.add_argument("--model",dest="model",type=str,help="path to model",default="")

args = parser.parse_args()

aFile = ROOT.TFile.Open(args.sigFile)
mergedTree = aFile.Get("mergedLeptonAnalyzer/mergedEl1_elTree")
mergedGsfTree = aFile.Get("mergedLeptonAnalyzer/mergedEl1Gsf_addGsfTree")
mergedWgtHist= aFile.Get("mergedLeptonAnalyzer/totWeightedSum")
mergedTotWgtSum = mergedWgtHist.GetBinContent(1)

dyFile = ROOT.TFile.Open(args.dyFile)
dyTree = dyFile.Get("mergedFakeAnalyzer/fake_elTree")
dyGsfTree = dyFile.Get("mergedFakeAnalyzer/fakeGsf_addGsfTree")
dyWgtHist = dyFile.Get("mergedFakeAnalyzer/totWeightedSum")
dyTotWgtSum = dyWgtHist.GetBinContent(1)

ttFile = ROOT.TFile.Open(args.ttFile)
ttTree = ttFile.Get("mergedFakeAnalyzer/fake_elTree")
ttGsfTree = ttFile.Get("mergedFakeAnalyzer/fakeGsf_addGsfTree")
ttWgtHist = ttFile.Get("mergedFakeAnalyzer/totWeightedSum")
ttTotWgtSum = ttWgtHist.GetBinContent(1)

dfProducer = zprep.DataframeInitializer(args.det,args.angle)
df_mergedEl, wgts_mergedEl = dfProducer.fill_arr(mergedTree,mergedGsfTree)
df_dyEl, wgts_dyEl = dfProducer.fill_arr(dyTree,dyGsfTree)
tt_dyEl, wgts_ttEl = dfProducer.fill_arr(ttTree,ttGsfTree)
col_names = dfProducer.col_names()

# normalize bkg by xsec
wgts_dyEl = wgts_dyEl*(6077.22*1000./dyTotWgtSum)
wgts_ttEl = wgts_ttEl*(831.76*1000./ttTotWgtSum)

# equalize sumwgt(bkg) to sumwgt(sig) to avoid class imbalance
sum_sig = np.sum(wgts_mergedEl)
sum_bkg = np.sum(wgts_dyEl) + np.sum(wgts_ttEl)

wgts_dyEl = wgts_dyEl*(sum_sig/sum_bkg)
wgts_ttEl = wgts_ttEl*(sum_sig/sum_bkg)

# concat bkg
df_bkg = pd.concat([df_dyEl,tt_dyEl], axis=0, ignore_index=True)
wgts_bkg = np.concatenate((wgts_dyEl,wgts_ttEl), axis=0)

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

x_train, x_test, y_train, y_test, wgts_train, wgts_test = sklearn.model_selection.train_test_split(x_total,y_total,wgts_total,test_size=0.5)

print('total number of merged objects = %d' % y_mergedEl.shape[0])
print('total number of bkg objects = %d' % y_bkgEl.shape[0])

print('bkg weights are')
print(wgts_dyEl)
print(wgts_ttEl)

dtrain = xgb.DMatrix(x_train, weight=wgts_train, label=y_train, feature_names=col_names)
dtest = xgb.DMatrix(x_test, weight=wgts_test, label=y_test, feature_names=col_names)
dtotal = xgb.DMatrix(x_total, weight=wgts_total, label=y_total, feature_names=col_names)

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
        'min_child_weight': hyperopt.hp.loguniform('min_child_weight',math.log(weightSum/1000.),math.log(weightSum/10.)) # 1/200 for EB dr1 & EE
    }

    trials = hyperopt.Trials()
    objec = lambda x: objective(x, dtotal)
    best = hyperopt.fmin(fn=objec, space=param_space, max_evals=256, algo=hyperopt.tpe.suggest, trials=trials)

    print('the best params are')
    print(best)

    np.savez('npz/'+args.angle+'_'+args.det+'.npz', meanstd=np_meanstd, bestparams=np.array(best))

npz_params = np.array([])

if args.load=="True" and os.path.splitext(args.model)[1]=='.model':
   npz_params = np.load('npz/'+os.path.splitext(os.path.basename(args.model))[0]+'.npz',allow_pickle=True)
else:
   npz_params = np.load('npz/'+args.angle+'_'+args.det+'.npz',allow_pickle=True)

param = npz_params['bestparams'][()] # numpy stores dict as 0-d array
param['max_depth'] = int(param['max_depth'])
param['objective'] = 'binary:logistic'
param['subsample'] = 0.5
param['eval_metric'] = 'logloss'
param['tree_method'] = 'exact'
param['nthread'] = 4

evallist = [(dtest, 'eval'), (dtrain, 'train')]
num_round = 200

bst = xgb.Booster(param)
early_stop = xgb.callback.EarlyStopping(rounds=20,metric_name='logloss',data_name='eval')

if args.opt=='True':
    bst = xgb.train(param, dtrain, num_round, evallist, callbacks=[early_stop])
    bst.save_model('model/'+args.angle+'_'+args.det+'.model')
else:
    if args.load=="True" and os.path.splitext(args.model)[1]=='.model':
        print('load model %s' % args.model)
        bst.load_model(args.model)
    else:
        bst.load_model('model/'+args.angle+'_'+args.det+'.model')

dTrainPredict    = bst.predict(dtrain)
dTestPredict     = bst.predict(dtest)

scoreMerged_train = dTrainPredict[ y_train == 1 ]
scoreMerged_test = dTestPredict[ y_test == 1 ]
scoreBkg_train = dTrainPredict[ y_train == 0 ]
scoreBkg_test = dTestPredict[ y_test == 0 ]

plotname = args.angle+'_'+args.det
if args.load=="True" and os.path.splitext(args.model)[1]=='.model':
    plotname += ('_'+os.path.splitext(os.path.basename(args.model))[0])

zvis.drawScoreOverlay(scoreMerged_train,scoreBkg_train,scoreMerged_test,scoreBkg_test,plotname)

if args.opt=='True':
    gain = bst.get_score( importance_type='gain')
    cover = bst.get_score(importance_type='cover')
    zvis.drawImportance(gain,cover,col_names,args.angle+'_'+args.det+'_importance')

np_meanstd = npz_params['meanstd'][()]
transformer.mean_ = np_meanstd[0]
transformer.scale_ = np_meanstd[1]
trans_mergedEl = transformer.transform(df_mergedEl)
trans_DY = transformer.transform(df_dyEl)
trans_TT = transformer.transform(tt_dyEl)
dSig = xgb.DMatrix(trans_mergedEl, weight=wgts_mergedEl, label=y_mergedEl, feature_names=col_names)
dDy = xgb.DMatrix(trans_DY, weight=wgts_dyEl, label=np.zeros(shape=(df_dyEl.shape[0],)), feature_names=col_names)
dTT = xgb.DMatrix(trans_TT, weight=wgts_ttEl, label=np.zeros(shape=(tt_dyEl.shape[0],)), feature_names=col_names)

dSigPredict = bst.predict(dSig)
dDyPredict = bst.predict(dDy)
dTTPredict = bst.predict(dTT)

zvis.drawScoreByProcess(dSigPredict, wgts_mergedEl, dDyPredict, wgts_dyEl, dTTPredict, wgts_ttEl, args.angle+'_'+args.det+'_scoreProcess')
