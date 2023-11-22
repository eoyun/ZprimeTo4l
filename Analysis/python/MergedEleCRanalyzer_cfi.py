import FWCore.ParameterSet.Config as cms

mergedEleCRanalyzer = cms.EDAnalyzer("MergedEleCRanalyzer",
  isMC = cms.bool(True),
  selectHT = cms.bool(False),
  selectZpT = cms.bool(False),
  maxHT = cms.double(70.),
  maxVpT = cms.double(50.),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  generator = cms.InputTag("generator"),
  srcGenPtc=cms.InputTag("prunedGenParticles"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  lheEvent = cms.InputTag("externalLHEProducer"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  addPackedCandMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddPackedCand"),
  union5x5dEtaIn = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5dEtaIn"),
  union5x5Energy = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","union5x5Energy"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  trigList = cms.vstring(
  # https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter"
  ),
  etThres1 = cms.double(35.),
  mergedEleSFmuHasTrk = cms.double(1.012),
  mergedEleSFmuNoTrk = cms.double(1.004),
  mergedEleSFcl95HasTrk = cms.double(0.07007),
  mergedEleSFcl95NoTrk = cms.double(0.2448),
  mergedEleSFupperHasTrk = cms.double(1.083),
  mergedEleSFupperNoTrk = cms.double(1.785),
  mergedElePolHasTrkStr = cms.string("0.0007899*x+0.970"),
  mergedElePolNoTrkStr = cms.string("0.004586*x+0.919")
)

mergedEleCRanalyzer20UL16APV = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL16APV.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16APV.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.)
)

mergedEleCRanalyzer20UL17 = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL17.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(35.),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2017.root"),
  modHeepSFmuEB1 = cms.double(1.015),
  modHeepSFmuEB2 = cms.double(1.002),
  modHeepSFmuEE = cms.double(1.003),
  modHeepSFcl95EB1 = cms.double(0.03528),
  modHeepSFcl95EB2 = cms.double(0.02765),
  modHeepSFcl95EE = cms.double(0.02814),
  modHeepSFupperEB1 = cms.double(1.118),
  modHeepSFupperEB2 = cms.double(1.093),
  modHeepSFupperEE = cms.double(1.144),
  modHeepPolEB1str = cms.string("1.119*0.00001*x+1.012"),
  modHeepPolEB2str = cms.string("2.415*0.00001*x+0.995"),
  modHeepPolEEstr = cms.string("2.153*0.00001*x+1.001")
)

mergedEleCRanalyzer20UL18 = mergedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/MEFF_20UL18.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle25_CaloIdL_MW_v*"
  ),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "Flag_EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  etThres1 = cms.double(28.),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2018.root"),
  modHeepSFmuEB1 = cms.double(0.998),
  modHeepSFmuEB2 = cms.double(0.991),
  modHeepSFmuEE = cms.double(1.005),
  modHeepSFcl95EB1 = cms.double(0.02489),
  modHeepSFcl95EB2 = cms.double(0.03772),
  modHeepSFcl95EE = cms.double(0.04466),
  modHeepSFupperEB1 = cms.double(1.110),
  modHeepSFupperEB2 = cms.double(1.128),
  modHeepSFupperEE = cms.double(1.125),
  modHeepPolEB1str = cms.string("-8.873*0.000001*x+1.000"),
  modHeepPolEB2str = cms.string("2.09*0.000001*x+0.990"),
  modHeepPolEEstr = cms.string("2.354*0.00001*x+0.999")
)
