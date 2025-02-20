from CRABClient.UserUtilities import config

config = config()

config.General.requestName = "WtoLNu-2Jets_0JEE"
config.General.workArea = "crab_projects"
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "runMergedLeptonIDImageMC_run3_cfg.py"
#config.JobType.maxMemoryMB = 4000
#config.JobType.numCores = 8

config.Data.inputDataset = "/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EEMiniAODv4-130X_mcRun3_2022_realistic_postEE_v6-v3/MINIAODSIM" 
config.Data.outLFNDirBase = "/store/user/yeo/"# % (getUsernameFromCRIC())
config.Data.outputDatasetTag = "WtoLNu-2Jets_0JEE"
config.Data.inputDBS = "global"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True

config.Site.storageSite = "T3_KR_KNU"
config.Site.whitelist = ["T2_*","T3_*"]
