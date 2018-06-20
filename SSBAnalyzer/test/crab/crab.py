from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'test'

config.section_('JobType')
config.JobType.pluginName  = 'Analysis'
config.JobType.psetName    = '../ssbanalyzer_cfg.py'
config.JobType.outputFiles = ['SSBTree.root']
#config.JobType.pyCfgParams = ['label=HWW', 'id=12345', 'scale=1', 'outputFile=stepB_MC_ggHww.root', 'doNoFilter=True', 'doMuonIsoId=True', 'doGen=True', 'doLHE=True', 'runPUPPISequence=True']
config.JobType.allowUndistributedCMSSW = True

config.section_('Data')
config.Data.inputDataset = '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM'
config.Data.splitting    = 'FileBased'  #'LumiBased'
config.Data.unitsPerJob  = 1  # Since files based, 10 files per job
config.Data.inputDBS     = 'global'
config.Data.outLFNDirBase = '/store/user/spak/CP_Violation/MC/test'

config.section_('Site')
config.Site.storageSite = 'T2_KR_KNU'
config.Site.blacklist = ['T2_US_Vanderbilt','T2_RU_IHEP','T1_RU_JINR']
