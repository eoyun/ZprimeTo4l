import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

def drawScoreByProcess(dSigPredict, sigWgt, dBkgPredicts, bkgWgts, labellist, colorlist, plotname, dirname="plot"):
    plt.figure(figsize=(6,4))
    plt.rc('font', size=12)
    plt.hist(dBkgPredicts, 100, stacked=True, weights=bkgWgts, label=labellist, range=(0,1), color=colorlist)
    plt.hist(dSigPredict, 100, weights=sigWgt, histtype=u'step', label='Sig', range=(0,1), color='navy')
    plt.grid()
    plt.yscale('log')
    plt.xlim([0,1])
    plt.xlabel('Score')
    plt.ylabel('a.u.')
    plt.legend(loc=(0.75,0.73))

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
