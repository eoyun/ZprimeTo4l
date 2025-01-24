from CRABClient.UserUtilities import config

config = config()

config.General.requestName = "HToAATo4L_H250A10_22EE"
config.General.workArea = "crab_projects"
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "runMergedLeptonIDImageMC_run3_cfg.py"
#config.JobType.maxMemoryMB = 4000
#config.JobType.numCores = 8

config.Data.inputDataset = "/HToAATo4L_H250A10_TuneCP5_13p6TeV-pythia8/yeo-Run3Summer22EE_MiniAODv4-8d987a117c9c64189fa569ace8eb11af/USER" 
config.Data.outLFNDirBase = "/store/user/yeo/"# % (getUsernameFromCRIC())
config.Data.outputDatasetTag = "HToAATo4L_H250A10_22EE"
config.Data.inputDBS = "phys03"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True

config.Site.storageSite = "T3_KR_KNU"
config.Site.whitelist = ["T2_*","T3_*"]
