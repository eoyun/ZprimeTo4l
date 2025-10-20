from CRABClient.UserUtilities import config

config = config()

config.General.requestName = "ParkingSingleMuon4_2024I"
config.General.workArea = "crab_projects"
config.General.transferLogs = True
config.General.transferOutputs = True

config.JobType.pluginName = "Analysis"
config.JobType.psetName = "runMergedLeptonIDJpsiAnalyzerData_run3_2024_cfg.py"
#config.JobType.maxMemoryMB = 4000
#config.JobType.numCores = 8

config.Data.inputDataset = "/ParkingSingleMuon4/Run2024I-PromptReco-v2/MINIAOD" 
config.Data.outLFNDirBase = "/store/user/yeo/"# % (getUsernameFromCRIC())
config.Data.outputDatasetTag = "ParkingSingleMuon4_2024I"
config.Data.inputDBS = "global"
config.Data.splitting = "FileBased"
config.Data.unitsPerJob = 1
config.Data.publication = False
config.Data.ignoreLocality = True

config.Site.storageSite = "T3_KR_KNU"
config.Site.whitelist = ["T2_*","T3_*"]
