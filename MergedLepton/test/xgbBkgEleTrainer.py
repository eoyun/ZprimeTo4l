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
parser.add_argument("--opt",dest="opt",type=str,help="perform training")
parser.add_argument("--load",dest="load",type=str,help="load model",default="False")
parser.add_argument("--model",dest="model",type=str,help="path to model",default="")

args = parser.parse_args()

aFile = ROOT.TFile.Open(args.sigFile)
mergedTree = aFile.Get("mergedLeptonAnalyzer/mergedEl2_elTree")
mergedWgtHist= aFile.Get("mergedLeptonAnalyzer/totWeightedSum")
mergedTotWgtSum = mergedWgtHist.GetBinContent(1)

dfProducer = zprep.DataframeInitializer(args.det,"dr1")

dyProcessor = zprep.SampleProcessor(args.dyFile,6077.22)
df_dyEl, wgts_dyEl = dyProcessor.read_bkg(dfProducer)

# ttProcessor = zprep.SampleProcessor(args.ttFile,831.76) # FXFX inclusive
ttProcessor = zprep.SampleProcessor(args.ttFile,54.17) # MLM dilepton
df_ttEl, wgts_ttEl = ttProcessor.read_bkg(dfProducer)

wjetsProcessor = zprep.SampleProcessor("WJets_bkg.root",53870.0)
df_wjetsEl, wgts_wjetsEl = wjetsProcessor.read_bkg(dfProducer)

qcdProcessor_1 = zprep.SampleProcessor("QCD_Pt-15to20_bkg.root",1324000.0)
df_qcdEl_1, wgts_qcdEl_1 = qcdProcessor_1.read_bkg(dfProducer)

qcdProcessor_2 = zprep.SampleProcessor("QCD_Pt-20to30_bkg.root",4896000.0)
df_qcdEl_2, wgts_qcdEl_2 = qcdProcessor_2.read_bkg(dfProducer)

qcdProcessor_3 = zprep.SampleProcessor("QCD_Pt-30to50_bkg.root",6447000.0)
df_qcdEl_3, wgts_qcdEl_3 = qcdProcessor_3.read_bkg(dfProducer)

qcdProcessor_4 = zprep.SampleProcessor("QCD_Pt-50to80_bkg.root",1988000.0)
df_qcdEl_4, wgts_qcdEl_4 = qcdProcessor_4.read_bkg(dfProducer)

qcdProcessor_5 = zprep.SampleProcessor("QCD_Pt-80to120_bkg.root",367500.0)
df_qcdEl_5, wgts_qcdEl_5 = qcdProcessor_5.read_bkg(dfProducer)

qcdProcessor_6 = zprep.SampleProcessor("QCD_Pt-120to170_bkg.root",66590.0)
df_qcdEl_6, wgts_qcdEl_6 = qcdProcessor_6.read_bkg(dfProducer)

qcdProcessor_7 = zprep.SampleProcessor("QCD_Pt-170to300_bkg.root",16620.0)
df_qcdEl_7, wgts_qcdEl_7 = qcdProcessor_7.read_bkg(dfProducer)

qcdProcessor_8 = zprep.SampleProcessor("QCD_Pt-300toInf_bkg.root",1104.0)
df_qcdEl_8, wgts_qcdEl_8 = qcdProcessor_8.read_bkg(dfProducer)

df_mergedEl, wgts_mergedEl = dfProducer.fill_bkg(mergedTree)
col_names = dfProducer.col_names()

# equalize sumwgt(bkg) to sumwgt(sig) to avoid class imbalance
sum_sig = np.sum(wgts_mergedEl)
sum_bkg = (np.sum(wgts_dyEl) + np.sum(wgts_ttEl) + np.sum(wgts_wjetsEl)
          + np.sum(wgts_qcdEl_1) + np.sum(wgts_qcdEl_2) + np.sum(wgts_qcdEl_3) + np.sum(wgts_qcdEl_4)
          + np.sum(wgts_qcdEl_5) + np.sum(wgts_qcdEl_6) + np.sum(wgts_qcdEl_7) + np.sum(wgts_qcdEl_8))

wgts_dyEl = wgts_dyEl*(sum_sig/sum_bkg)
wgts_ttEl = wgts_ttEl*(sum_sig/sum_bkg)
wgts_wjetsEl = wgts_wjetsEl*(sum_sig/sum_bkg)
wgts_qcdEl_1 = wgts_qcdEl_1*(sum_sig/sum_bkg)
wgts_qcdEl_2 = wgts_qcdEl_2*(sum_sig/sum_bkg)
wgts_qcdEl_3 = wgts_qcdEl_3*(sum_sig/sum_bkg)
wgts_qcdEl_4 = wgts_qcdEl_4*(sum_sig/sum_bkg)
wgts_qcdEl_5 = wgts_qcdEl_5*(sum_sig/sum_bkg)
wgts_qcdEl_6 = wgts_qcdEl_6*(sum_sig/sum_bkg)
wgts_qcdEl_7 = wgts_qcdEl_7*(sum_sig/sum_bkg)
wgts_qcdEl_8 = wgts_qcdEl_8*(sum_sig/sum_bkg)

# concat bkg
df_bkg = pd.concat([df_dyEl,df_ttEl,df_wjetsEl,
                    df_qcdEl_1,
                    df_qcdEl_2,
                    df_qcdEl_3,
                    df_qcdEl_4,
                    df_qcdEl_5,
                    df_qcdEl_6,
                    df_qcdEl_7,
                    df_qcdEl_8], axis=0, ignore_index=True)
wgts_bkg = np.concatenate((wgts_dyEl,wgts_ttEl,wgts_wjetsEl,
                           wgts_qcdEl_1,
                           wgts_qcdEl_2,
                           wgts_qcdEl_3,
                           wgts_qcdEl_4,
                           wgts_qcdEl_5,
                           wgts_qcdEl_6,
                           wgts_qcdEl_7,
                           wgts_qcdEl_8), axis=0)

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
    best = hyperopt.fmin(fn=objec, space=param_space, max_evals=64, algo=hyperopt.tpe.suggest, trials=trials)

    print('the best params are')
    print(best)

    np.savez('npz/'+'bkg'+'_'+args.det+'.npz', meanstd=np_meanstd, bestparams=np.array(best))

npz_params = np.array([])

if args.load=="True" and os.path.splitext(args.model)[1]=='.model':
   npz_params = np.load('npz/'+os.path.splitext(os.path.basename(args.model))[0]+'.npz',allow_pickle=True)
else:
   npz_params = np.load('npz/'+'bkg'+'_'+args.det+'.npz',allow_pickle=True)

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
    bst.save_model('model/'+'bkg'+'_'+args.det+'.model')
else:
    if args.load=="True" and os.path.splitext(args.model)[1]=='.model':
        print('load model %s' % args.model)
        bst.load_model(args.model)
    else:
        bst.load_model('model/'+'bkg'+'_'+args.det+'.model')

dTrainPredict    = bst.predict(dtrain)
dTestPredict     = bst.predict(dtest)

scoreMerged_train = dTrainPredict[ y_train == 1 ]
scoreMerged_test = dTestPredict[ y_test == 1 ]
scoreBkg_train = dTrainPredict[ y_train == 0 ]
scoreBkg_test = dTestPredict[ y_test == 0 ]
wgtsMerged_train = wgts_train[ y_train == 1 ]
wgtsMerged_test = wgts_test[ y_test == 1 ]
wgtsBkg_train = wgts_train[ y_train == 0 ]
wgtsBkg_test = wgts_test[ y_test == 0 ]

plotname = 'bkg'+'_'+args.det
if args.load=="True" and os.path.splitext(args.model)[1]=='.model':
    plotname += ('_'+os.path.splitext(os.path.basename(args.model))[0])

zvis.drawScoreOverlay(scoreMerged_train,
                      wgtsMerged_train,
                      scoreBkg_train,
                      wgtsBkg_train,
                      scoreMerged_test,
                      wgtsMerged_test,
                      scoreBkg_test,
                      wgtsBkg_test,
                      plotname)

if args.opt=='True':
    gain = bst.get_score( importance_type='gain')
    cover = bst.get_score(importance_type='cover')
    zvis.drawImportance(gain,cover,col_names,'bkg'+'_'+args.det+'_importance')

np_meanstd = npz_params['meanstd'][()]
transformer.mean_ = np_meanstd[0]
transformer.scale_ = np_meanstd[1]
trans_mergedEl = transformer.transform(df_mergedEl)
trans_DY = transformer.transform(df_dyEl)
trans_TT = transformer.transform(df_ttEl)
trans_WJ = transformer.transform(df_wjetsEl)
trans_QCD1 = transformer.transform(df_qcdEl_1)
trans_QCD2 = transformer.transform(df_qcdEl_2)
trans_QCD3 = transformer.transform(df_qcdEl_3)
trans_QCD4 = transformer.transform(df_qcdEl_4)
trans_QCD5 = transformer.transform(df_qcdEl_5)
trans_QCD6 = transformer.transform(df_qcdEl_6)
trans_QCD7 = transformer.transform(df_qcdEl_7)
trans_QCD8 = transformer.transform(df_qcdEl_8)
dSig = xgb.DMatrix(trans_mergedEl, weight=wgts_mergedEl, label=y_mergedEl, feature_names=col_names)
dDy = xgb.DMatrix(trans_DY, weight=wgts_dyEl, label=np.zeros(shape=(df_dyEl.shape[0],)), feature_names=col_names)
dTT = xgb.DMatrix(trans_TT, weight=wgts_ttEl, label=np.zeros(shape=(df_ttEl.shape[0],)), feature_names=col_names)
dWJ = xgb.DMatrix(trans_WJ, weight=wgts_wjetsEl, label=np.zeros(shape=(df_wjetsEl.shape[0],)), feature_names=col_names)
dQCD1 = xgb.DMatrix(trans_QCD1, weight=wgts_qcdEl_1, label=np.zeros(shape=(df_qcdEl_1.shape[0],)), feature_names=col_names)
dQCD2 = xgb.DMatrix(trans_QCD2, weight=wgts_qcdEl_2, label=np.zeros(shape=(df_qcdEl_2.shape[0],)), feature_names=col_names)
dQCD3 = xgb.DMatrix(trans_QCD3, weight=wgts_qcdEl_3, label=np.zeros(shape=(df_qcdEl_3.shape[0],)), feature_names=col_names)
dQCD4 = xgb.DMatrix(trans_QCD4, weight=wgts_qcdEl_4, label=np.zeros(shape=(df_qcdEl_4.shape[0],)), feature_names=col_names)
dQCD5 = xgb.DMatrix(trans_QCD5, weight=wgts_qcdEl_5, label=np.zeros(shape=(df_qcdEl_5.shape[0],)), feature_names=col_names)
dQCD6 = xgb.DMatrix(trans_QCD6, weight=wgts_qcdEl_6, label=np.zeros(shape=(df_qcdEl_6.shape[0],)), feature_names=col_names)
dQCD7 = xgb.DMatrix(trans_QCD7, weight=wgts_qcdEl_7, label=np.zeros(shape=(df_qcdEl_7.shape[0],)), feature_names=col_names)
dQCD8 = xgb.DMatrix(trans_QCD8, weight=wgts_qcdEl_8, label=np.zeros(shape=(df_qcdEl_8.shape[0],)), feature_names=col_names)

dSigPredict = bst.predict(dSig)
dDyPredict = bst.predict(dDy)
dTTPredict = bst.predict(dTT)
dWJPredict = bst.predict(dWJ)
dQCD1Predict = bst.predict(dQCD1)
dQCD2Predict = bst.predict(dQCD2)
dQCD3Predict = bst.predict(dQCD3)
dQCD4Predict = bst.predict(dQCD4)
dQCD5Predict = bst.predict(dQCD5)
dQCD6Predict = bst.predict(dQCD6)
dQCD7Predict = bst.predict(dQCD7)
dQCD8Predict = bst.predict(dQCD8)

zvis.drawScoreByProcess(
    dSigPredict,
    wgts_mergedEl,
    list(np.array([dTTPredict,dDyPredict,dWJPredict,
    dQCD1Predict,dQCD2Predict,dQCD3Predict,dQCD4Predict,dQCD5Predict,dQCD6Predict,dQCD7Predict,dQCD8Predict],dtype=object)),
    list(np.array([wgts_ttEl,wgts_dyEl,wgts_wjetsEl,
    wgts_qcdEl_1,wgts_qcdEl_2,wgts_qcdEl_3,wgts_qcdEl_4,wgts_qcdEl_5,wgts_qcdEl_6,wgts_qcdEl_7,wgts_qcdEl_8],dtype=object)),
    ['TT','DY','WJets','QCD_Pt-15to20','QCD_Pt-20to30','QCD_Pt-30to50','QCD_Pt-50to80','QCD_Pt-80to120',
    'QCD_Pt-120to170','QCD_Pt-170to300','QCD_Pt-300toInf'],
    ['tomato','wheat','green','lightcyan','paleturquoise','turquoise','cyan','darkturquoise','darkcyan','cadetblue','steelblue'],
    'bkg'+'_'+args.det+'_scoreProcess'
)
