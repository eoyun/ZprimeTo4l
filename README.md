# ZprimeTo4l
repository for Z'to4l study

```sh
cmsrel CMSSW_10_6_39
cd CMSSW_10_6_39/src
# checkout ZprimeTo4l anlaysis package
git cms-merge-topic SanghyunKo:Zprime_from-CMSSW_10_6_39

# checkout RunII UL S&S energy correction files
git clone -b ULSSfiles_correctScaleSysMC https://github.com/jainshilpi/EgammaAnalysis-ElectronTools.git EgammaAnalysis/ElectronTools/data/
git cms-addpkg EgammaAnalysis/ElectronTools

# checkout EgammaPostRecoTools
git clone git@github.com:cms-egamma/EgammaPostRecoTools.git EgammaUser/EgammaPostRecoTools

# checkout TnPTreeProducer (needed for SimpleEventCounter.cc)
git clone -b RunIIfinal https://github.com/cms-egamma/EgammaAnalysis-TnPTreeProducer.git EgammaAnalysis/TnPTreeProducer

# compile
scram b -j8
```

Example cfg files and analysis scripts can be found in the `Analysis/test` or `MergedLepton/test` folders
