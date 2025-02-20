from CRABClient.UserUtilities import config

config = config()

config.General.requestName = "DYto2L_2J_MLMEE"
config.General.workArea = "crab_projects"
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "runMergedLeptonIDImageMC_run3_cfg.py"
#config.JobType.maxMemoryMB = 4000
#config.JobType.numCores = 8

config.Data.inputDataset = "dataset" 
config.Data.outLFNDirBase = "/store/user/yeo/"# % (getUsernameFromCRIC())
config.Data.outputDatasetTag = "DYto2L_2J_MLMEE"
config.Data.inputDBS = "global"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True

config.Site.storageSite = "T3_KR_KNU"
config.Site.whitelist = ["T2_*","T3_*"]
