import numpy as np
import sklearn
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

class modelPerformance(object):
    """docstring for modelPerformance."""
    def __init__(self,fprTr,tprTr,thrTr,aucTr,fprTe,tprTe,thrTe,aucTe):

        super(modelPerformance, self).__init__()
        self.fprTrain_ = fprTr
        self.tprTrain_ = tprTr
        self.thrTrain_ = thrTr
        self.aucTrain_ = aucTr
        self.fprTest_ = fprTe
        self.tprTest_ = tprTe
        self.thrTest_ = thrTe
        self.aucTest_ = aucTe

def calROC(dTrainPredict, dTestPredict, y_train, y_test):
    fprTr, tprTr, thrTr = sklearn.metrics.roc_curve(y_train, dTrainPredict)
    aucTr = sklearn.metrics.roc_auc_score(y_train, dTrainPredict)

    fprTe, tprTe, thrTe = sklearn.metrics.roc_curve(y_test, dTestPredict)
    aucTe = sklearn.metrics.roc_auc_score(y_test, dTestPredict)

    return modelPerformance(fprTr,tprTr,thrTr,aucTr,fprTe,tprTe,thrTe,aucTe)

def drawScoreByProcess(dSigPredict, sigWgt, dBkgPredicts, bkgWgts, labellist, colorlist, nbins, plotname, dirname="plot"):
    plt.figure(figsize=(6,4))
    plt.rc('font', size=12)
    plt.hist(dBkgPredicts, nbins, stacked=True, weights=bkgWgts, label=labellist, range=(0,1), color=colorlist)
    plt.hist(dSigPredict, nbins, weights=sigWgt, histtype=u'step', range=(0,1), color='navy')
    plt.title(r"$\bf{CMS}$"+"$\it{\;Simulation,\;work\;in\;progress}$", loc='left',fontsize=12)
    plt.grid()
    plt.yscale('log')
    plt.xlim([0,1])
    ymin, ymax = plt.gca().get_ylim()
    plt.ylim([ymin,10.*ymax])
    plt.xlabel('Score')
    plt.ylabel('a.u.')
    plt.legend(loc=(0.8,0.73), fontsize=8)

    plt.savefig(dirname+'/'+plotname+'.png',dpi=300, bbox_inches='tight')
    plt.close()

    return

def drawImportance(gain, cover, colname_full, plotname, dirname="plot"):
    colname = [ col for col in colname_full if col in gain.keys() ]
    valGain     = np.asarray( [ gain[x]  for x in colname ] )
    sortedCover = np.asarray( [ cover[x] for x in colname ] )

    plt.figure(figsize=(6,4))
    barwidth = 0.3
    b1 = plt.barh(np.arange(len(gain)) -barwidth/2., 100.*valGain/np.sum(valGain),         barwidth, color='r', label='gain')
    b2 = plt.barh(np.arange(len(cover))+barwidth/2., 100.*sortedCover/np.sum(sortedCover), barwidth, color='b', label='cover')
    plt.yticks(range(len(gain)), colname, fontsize=5)
    plt.title(r"$\bf{CMS}$"+"$\it{\;Simulation,\;work\;in\;progress}$", loc='left')
    plt.legend( (b1[0],b2[0]), ('gain','cover'), fontsize=8 )

    plt.savefig(dirname+'/'+plotname+'.png',dpi=300, bbox_inches='tight')
    plt.close()

    return

def drawROC(modelPerforms_list, plotname, dirname="plot"):
    plt.figure(figsize=(6,4))

    for idx in range(len(modelPerforms_list)):
        plt.plot(modelPerforms_list[idx].fprTrain_, modelPerforms_list[idx].tprTrain_, color='lightsalmon')
        plt.plot(modelPerforms_list[idx].fprTest_,  modelPerforms_list[idx].tprTest_,  color='lightblue')

    # pick the fold with the lowest abs( test AUC - train AUC )
    auclist = [ aperform.aucTest_ for aperform in modelPerforms_list ]
    maxauc = max(auclist)
    idx_max = auclist.index(maxauc)

    plt.plot(modelPerforms_list[idx_max].fprTrain_, modelPerforms_list[idx_max].tprTrain_, color='r', label='Train ROC (AUC = %.4f)' % modelPerforms_list[idx_max].aucTrain_)
    plt.plot(modelPerforms_list[idx_max].fprTest_,  modelPerforms_list[idx_max].tprTest_,  color='b', label='Test ROC (AUC = %.4f)' % modelPerforms_list[idx_max].aucTest_)
    plt.title(r"$\bf{CMS}$"+"$\it{\;Simulation,\;work\;in\;progress}$", loc='left')

    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.legend(loc='lower right')
    plt.grid()

    plt.savefig(dirname+'/'+plotname+'.png',dpi=300, bbox_inches='tight')
    plt.close()

    return idx_max

def drawThr(amodelperform, targetFpr, plotname, dirname="plot"):
    targetThres = np.interp(targetFpr,amodelperform.fprTest_,amodelperform.thrTest_)

    plt.figure(figsize=(6,4))
    plt.plot(amodelperform.thrTrain_, amodelperform.fprTrain_, color='r', label='Train thr')
    plt.plot(amodelperform.thrTest_,  amodelperform.fprTest_,  color='b', label='Test thr (%.3f at FPR=%.2f)' % (targetThres,targetFpr))
    plt.title(r"$\bf{CMS}$"+"$\it{\;Simulation,\;work\;in\;progress}$", loc='left')
    plt.xlabel('Threshold')
    plt.ylabel('False Positive Rate')
    plt.xlim(xmin=0.0, xmax=1.0)
    plt.legend(loc='upper right')
    plt.grid()

    plt.savefig(dirname+'/'+plotname+'.png',dpi=300, bbox_inches='tight')
    plt.close()

    return targetThres
