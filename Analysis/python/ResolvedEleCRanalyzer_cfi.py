import FWCore.ParameterSet.Config as cms

resolvedEleCRanalyzer = cms.EDAnalyzer("ResolvedEleCRanalyzer",
  isMC = cms.bool(True),
  generator = cms.InputTag("generator"),
  pileupSummary = cms.InputTag("slimmedAddPileupInfo"),
  triggerResults = cms.InputTag("TriggerResults","","HLT"),
  srcEle = cms.InputTag("slimmedElectrons"),
  srcPv = cms.InputTag("offlineSlimmedPrimaryVertices"),
  addGsfTrkMap = cms.InputTag("modifiedHEEPIDVarValueMaps2nd","eleAddGsfTrk"),
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016postVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL16.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  METfilters = cms.InputTag("TriggerResults","","PAT"),
  triggerObjects = cms.InputTag("slimmedPatTrigger"),
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
  trigList = cms.vstring(
  # https://docs.google.com/spreadsheets/d/1Yy1VYIp-__pVUDFUs6JDNUh4XGlL2SlOsqXjdCsNiGU/edit#gid=663660886
  # https://twiki.cern.ch/twiki/bin/view/CMS/EgHLTRunIISummary
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.),
  etThres2 = cms.double(20.),
  ffSystCL95 = cms.double(0.5),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016postVFP.root"),
  modHeepSFmuEB1 = cms.double(0.992),
  modHeepSFmuEB2 = cms.double(1.010),
  modHeepSFmuEE = cms.double(1.035),
  modHeepSFcl95EB1 = cms.double(0.03516),
  modHeepSFcl95EB2 = cms.double(0.03548),
  modHeepSFcl95EE = cms.double(0.06211),
  modHeepSFupperEB1 = cms.double(1.098),
  modHeepSFupperEB2 = cms.double(1.142),
  modHeepSFupperEE = cms.double(1.344),
  modHeepPolEB1str = cms.string("-7.294*0.000001*x+0.995"),
  modHeepPolEB2str = cms.string("-6.237*0.000001*x+1.012"),
  modHeepPolEEstr = cms.string("7.755*0.00001*x+1.026")
)

resolvedEleCRanalyzer20UL16APV = resolvedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2016preVFP.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL16APV.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL16.root"),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*",
    "HLT_DoubleEle33_CaloIdL_GsfTrkIdVL_v*"
  ),
  etThres1 = cms.double(35.),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2016preVFP.root"),
  modHeepSFmuEB1 = cms.double(0.995),
  modHeepSFmuEB2 = cms.double(1.003),
  modHeepSFmuEE = cms.double(1.007),
  modHeepSFcl95EB1 = cms.double(0.03915),
  modHeepSFcl95EB2 = cms.double(0.03809),
  modHeepSFcl95EE = cms.double(0.05239),
  modHeepSFupperEB1 = cms.double(1.116),
  modHeepSFupperEB2 = cms.double(1.112),
  modHeepSFupperEE = cms.double(1.227),
  modHeepPolEB1str = cms.string("-1.039*0.00001*x+0.998"),
  modHeepPolEB2str = cms.string("-3.185*0.00001*x+1.012"),
  modHeepPolEEstr = cms.string("7.773*0.000001*x+1.005")
)

resolvedEleCRanalyzer20UL17 = resolvedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2017.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL17.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL17.root"),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  trigList = cms.vstring(
    "HLT_DoubleEle33_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(35.),
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

resolvedEleCRanalyzer20UL18 = resolvedEleCRanalyzer.clone(
  recoSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_ptAbove20_EGM2D_UL2018.root"),
  modHeepSFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/egammaEffi_modHeep_EGM2D_UL2018.root"),
  FFpath = cms.FileInPath("ZprimeTo4l/Analysis/data/REFF_20UL18.root"),
  PUrwgt = cms.FileInPath("ZprimeTo4l/Analysis/data/PUrwgt20UL18.root"),
  METfilterList = cms.vstring(
    "Flag_goodVertices",
    "Flag_globalSuperTightHalo2016Filter",
    "Flag_HBHENoiseFilter",
    "Flag_HBHENoiseIsoFilter",
    "EcalDeadCellTriggerPrimitiveFilter",
    "Flag_BadPFMuonFilter",
    "Flag_BadPFMuonDzFilter",
    "Flag_eeBadScFilter",
    "Flag_hfNoisyHitsFilter",
    "Flag_ecalBadCalibFilter"
  ),
  trigList = cms.vstring(
    "HLT_DoubleEle25_CaloIdL_MW_v*"
  ),
  etThres1 = cms.double(28.),
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
