import ROOT
import numpy as np
import pandas as pd
import sklearn
import xgboost as xgb
import hyperopt
import math
import sys
import os

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# usage
# train: python3 mergedElePreprocess.py HxAx_el.root True True dr1 EB
# test : python3 mergedElePreprocess.py H2000A10_el.root False False dr2 EB model/HxAx_el_dr2_EB.model

argv = sys.argv
filename = argv[1]

aFile = ROOT.TFile.Open(filename)
resolvedTree1 = aFile.Get("mergedLeptonAnalyzer/heep1_elTree")
resolvedTree2 = aFile.Get("mergedLeptonAnalyzer/heep2_elTree")
mergedTree = aFile.Get("mergedLeptonAnalyzer/mergedEl1_elTree")

resolvedGsfTree1 = aFile.Get("mergedLeptonAnalyzer/heep1Gsf_addGsfTree")
resolvedGsfTree2 = aFile.Get("mergedLeptonAnalyzer/heep2Gsf_addGsfTree")
mergedGsfTree = aFile.Get("mergedLeptonAnalyzer/mergedEl1Gsf_addGsfTree")

nFeature = 17

arr_mergedEl = np.zeros(shape=(mergedTree.GetEntries(),nFeature))
arr_resolvedEl = np.zeros(shape=(resolvedTree1.GetEntries()+resolvedTree2.GetEntries(),nFeature))

def deltaPhi(phi1, phi2):
    dphi = phi2 - phi1
    if dphi > math.pi:
        dphi -= 2.*math.pi
    elif dphi <= -math.pi:
        dphi += 2.*math.pi

    return dphi

def deltaR2(eta1,phi1,eta2,phi2):
    return abs(eta1-eta2)**2 + abs(deltaPhi(phi1,phi2))**2

def fill_arr(aTree,aGsfTree,anArr,startIdx):
    for iEl, el in enumerate(aTree):
        if argv[5]=='EB':
            if abs(el.eta) > 1.5:
                continue
        elif argv[5]=='EE':
            if abs(el.eta) < 1.5:
                continue
        else:
            NameError('Please check EB/EE argument, the current argument is %s' % argv[5])

        aGsfTree.GetEntry(iEl)

        if argv[4]=='ss':
            if argv[5]=='EB':
                if not abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) > 2*0.0174:
                    continue
            else:
                if not abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) > 2*0.00864*abs(math.sinh(el.eta)):
                    continue
        elif argv[4]=='trk':
            if argv[5]=='EB':
                if not abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) < 0.0174:
                    continue
            else:
                if not abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) < 0.00864*abs(math.sinh(el.eta)):
                    continue
        elif argv[4]=='int':
            if argv[5]=='EB':
                if not ( abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) > 0.0174 and
                         abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) < 2.*0.0174 ):
                    continue
            else:
                if not ( abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) > 0.00864*abs(math.sinh(el.eta)) and
                         abs(deltaPhi(aGsfTree.Gsfphi,el.phi)) < 2.*0.00864*abs(math.sinh(el.eta)) ):
                    continue
        elif argv[4]=='dr1':
            if argv[5]=='EB':
                if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < 0.0174**2:
                    continue
            else:
                if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < (0.00864*abs(math.sinh(el.eta)))**2:
                    continue
        elif argv[4]=='dr2':
            if argv[5]=='EB':
                if not ( deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > 0.0174**2 and
                         deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < (2*0.0174)**2 ):
                    continue
            else:
                if not ( deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > (0.00864*abs(math.sinh(el.eta)))**2 and
                         deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) < (2*0.00864*abs(math.sinh(el.eta)))**2 ):
                    continue
        elif argv[4]=='dr3':
            if argv[5]=='EB':
                if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > (2*0.0174)**2:
                    continue
            else:
                if not deltaR2(aGsfTree.Gsfeta,aGsfTree.Gsfphi,el.eta,el.phi) > (2*0.00864*abs(math.sinh(el.eta)))**2:
                    continue
        elif argv[4]=='all':
            pass
        else:
            NameError('Please check dPhi(trk/int/ss) argument, the current argument is %s' % argv[4])

        idx = iEl + startIdx
        anArr[idx,0] = float(el.charge)
        anArr[idx,1] = el.etaSCWidth
        anArr[idx,2] = el.phiSCWidth
        anArr[idx,3] = el.full5x5_sigmaIetaIeta
        anArr[idx,4] = el.full5x5_sigmaIphiIphi
        anArr[idx,5] = el.full5x5_E1x5/el.full5x5_E5x5
        anArr[idx,6] = el.full5x5_E2x5/el.full5x5_E5x5
        anArr[idx,7] = el.full5x5_r9
        anArr[idx,8] = el.dEtaIn
        anArr[idx,9] = el.dPhiIn
        anArr[idx,10] = el.dPhiSeed
        anArr[idx,11] = el.dEtaEle
        anArr[idx,12] = el.dPhiEle
        anArr[idx,13] = el.dEtaSeed
        anArr[idx,14] = el.EseedOverP
        anArr[idx,15] = el.EOverP
        anArr[idx,16] = el.fbremSC

fill_arr(mergedTree,mergedGsfTree,arr_mergedEl,0)
fill_arr(resolvedTree1,resolvedGsfTree1,arr_resolvedEl,0)
fill_arr(resolvedTree2,resolvedGsfTree2,arr_resolvedEl,resolvedTree1.GetEntries())

col_names = [
'charge',
'etaSCWidth',
'phiSCWidth',
'full5x5_sigmaIetaIeta',
'full5x5_sigmaIphiIphi',
'full5x5_E1x5/E5x5',
'full5x5_E2x5/E5x5',
'full5x5_r9',
'dEtaIn',
'dPhiIn',
'dPhiSeed',
'dEtaEle',
'dPhiEle',
'dEtaSeed',
'EseedOverP',
'EOverP',
'fbremSC'
]

# remove rows with only zeros
arr_mergedEl = arr_mergedEl[~np.all(arr_mergedEl==0.0,axis=1)]
arr_resolvedEl = arr_resolvedEl[~np.all(arr_resolvedEl==0.0,axis=1)]

df_mergedEl = pd.DataFrame(arr_mergedEl,columns=col_names)
df_resolvedEl = pd.DataFrame(arr_resolvedEl,columns=col_names)
df_total = pd.concat([df_mergedEl,df_resolvedEl], axis=0, ignore_index=True)

y_mergedEl = np.ones(shape=(df_mergedEl.shape[0],))
y_resolvedEl = np.zeros(shape=(df_resolvedEl.shape[0],))
y_total = np.concatenate((y_mergedEl,y_resolvedEl), axis=0)

transformer = sklearn.preprocessing.StandardScaler()
x_total = transformer.fit_transform(df_total)
np_meanstd = np.array([transformer.mean_,transformer.scale_])

x_train, x_test, y_train, y_test = sklearn.model_selection.train_test_split(x_total,y_total,test_size=0.5)

wgts_val = sklearn.utils.class_weight.compute_class_weight('balanced',np.unique(y_total),y_total)
wgts_train = np.full(y_train.shape[0],1.)
wgts_test = np.full(y_test.shape[0],1.)
wgts_total = np.full(y_total.shape[0],1.)

print('total number of merged objects = %d' % y_mergedEl.shape[0])
print('total number of resolved objects = %d' % y_resolvedEl.shape[0])

print('weights are')
print(wgts_val)

for idx,val in enumerate(wgts_val):
    wgts_train = np.multiply(wgts_train,np.where(y_train==idx,val,1.))
    wgts_test = np.multiply(wgts_test,np.where(y_test==idx,val,1.))
    wgts_total = np.multiply(wgts_total,np.where(y_total==idx,val,1.))

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

if argv[2]=="True": # hyper-optimization
    weightSum = np.sum(wgts_total)

    param_space = {
        'max_depth': hyperopt.hp.quniform('max_depth',1,5,1), # lower depth, more robust
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

    np.savez('npz/'+os.path.splitext(filename)[0]+'_'+argv[4]+'_'+argv[5]+'.npz', meanstd=np_meanstd, bestparams=np.array(best))

npz_params = np.array([])

if len(argv) > 6 and os.path.splitext(argv[6])[1]=='.model':
   npz_params = np.load('npz/'+os.path.splitext(os.path.basename(argv[6]))[0]+'.npz',allow_pickle=True)
else:
   npz_params = np.load('npz/'+os.path.splitext(filename)[0]+'_'+argv[4]+'_'+argv[5]+'.npz',allow_pickle=True)

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

if argv[3]=='True':
    bst = xgb.train(param, dtrain, num_round, evallist, callbacks=[early_stop])
    bst.save_model('model/'+os.path.splitext(filename)[0]+'_'+argv[4]+'_'+argv[5]+'.model')
else:
    if len(argv) > 6 and os.path.splitext(argv[6])[1]=='.model':
        print('load model %s' % argv[6])
        bst.load_model(argv[6])
    else:
        bst.load_model('model/'+os.path.splitext(filename)[0]+'_'+argv[4]+'_'+argv[5]+'.model')

dTrainPredict    = bst.predict(dtrain)
dTestPredict     = bst.predict(dtest)

scoreMerged_train = dTrainPredict[ y_train == 1 ]
scoreMerged_test = dTestPredict[ y_test == 1 ]
scoreResolved_train = dTrainPredict[ y_train == 0 ]
scoreResolved_test = dTestPredict[ y_test == 0 ]

def drawScoreOverlay(dSigPredictTrain, dBkgPredictTrain, dSigPredictTest, dBkgPredictTest, plotname, dirname="plot"):
    plt.figure(figsize=(6,4))
    plt.rc('font', size=12)
    plt.hist(dSigPredictTrain, 100, density=True, histtype=u'step', label='Train:Sig', range=(0,1), color='navy')
    plt.hist(dBkgPredictTrain, 100, density=True, histtype=u'step', label='Train:Bkg', range=(0,1), color='maroon')
    plt.hist(dSigPredictTest, 100, density=True, alpha=0.4, label='Test:Sig', range=(0,1), color='blue')
    plt.hist(dBkgPredictTest, 100, density=True, alpha=0.4, label='Test:Bkg', range=(0,1), color='red')
    plt.grid()
    plt.yscale('log')
    plt.xlim([0,1])
    plt.xlabel('Score')
    plt.ylabel('a.u.')
    plt.legend(loc=(0.65,0.65))

    plt.savefig(dirname+'/'+plotname+'.png',dpi=300, bbox_inches='tight')
    plt.close()

    return

def drawImportance(gain, cover, colname_full, plotname, dirname="plot"):
    colname = [ col for col in colname_full if col in gain.keys() ]
    valGain     = np.asarray( [ gain[x]  for x in colname ] )
    sortedCover = np.asarray( [ cover[x] for x in colname ] )

    plt.figure(figsize=(6,4))
    barwidth = 0.4
    b1 = plt.barh(np.arange(len(gain)) -barwidth/2., 100.*valGain/np.sum(valGain),         barwidth, color='r', label='gain')
    b2 = plt.barh(np.arange(len(cover))+barwidth/2., 100.*sortedCover/np.sum(sortedCover), barwidth, color='b', label='cover')
    plt.yticks(range(len(gain)), colname, fontsize=5)
    plt.legend( (b1[0],b2[0]), ('gain','cover'), fontsize=5 )

    plt.savefig(dirname+'/'+plotname+'.png',dpi=300, bbox_inches='tight')
    plt.close()

    return

plotname = os.path.splitext(filename)[0]+'_'+argv[4]+'_'+argv[5]
if len(argv) > 6 and os.path.splitext(argv[6])[1]=='.model':
    plotname += ('_'+os.path.splitext(os.path.basename(argv[6]))[0])

drawScoreOverlay(scoreMerged_train,scoreResolved_train,scoreMerged_test,scoreResolved_test,plotname)

if argv[3]=='True':
    gain = bst.get_score( importance_type='gain')
    cover = bst.get_score(importance_type='cover')
    drawImportance(gain,cover,col_names,os.path.splitext(filename)[0]+'_'+argv[4]+'_'+argv[5]+'_importance')
