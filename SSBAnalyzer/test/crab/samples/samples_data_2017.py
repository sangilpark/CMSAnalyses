#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# DAS query: https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+dataset%3D%2F*%2FRun2017*31Mar2018*%2FMINIAOD
#
# For GT and more, see https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookMiniAOD
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

samples['DoubleEG_Run2017B-31Mar2018-v1']       = ['/DoubleEG/Run2017B-31Mar2018-v1/MINIAOD',       ['label=DoubleEG']]
samples['DoubleMuon_Run2017B-31Mar2018-v1']     = ['/DoubleMuon/Run2017B-31Mar2018-v1/MINIAOD',     ['label=DoubleMuon']]
samples['MuonEG_Run2017B-31Mar2018-v1']         = ['/MuonEG/Run2017B-31Mar2018-v1/MINIAOD',         ['label=MuonEG']]
samples['SingleElectron_Run2017B-31Mar2018-v1'] = ['/SingleElectron/Run2017B-31Mar2018-v1/MINIAOD', ['label=SingleElectron']]
samples['SingleMuon_Run2017B-31Mar2018-v1']     = ['/SingleMuon/Run2017B-31Mar2018-v1/MINIAOD',     ['label=SingleMuon']]

samples['DoubleEG_Run2017C-31Mar2018-v1']       = ['/DoubleEG/Run2017C-31Mar2018-v1/MINIAOD',       ['label=DoubleEG']]
samples['DoubleMuon_Run2017C-31Mar2018-v1']     = ['/DoubleMuon/Run2017C-31Mar2018-v1/MINIAOD',     ['label=DoubleMuon']]
samples['MuonEG_Run2017C-31Mar2018-v1']         = ['/MuonEG/Run2017C-31Mar2018-v1/MINIAOD',         ['label=MuonEG']]
samples['SingleElectron_Run2017C-31Mar2018-v1'] = ['/SingleElectron/Run2017C-31Mar2018-v1/MINIAOD', ['label=SingleElectron']]
samples['SingleMuon_Run2017C-31Mar2018-v1']     = ['/SingleMuon/Run2017C-31Mar2018-v1/MINIAOD',     ['label=SingleMuon']]

samples['DoubleEG_Run2017D-31Mar2018-v1']       = ['/DoubleEG/Run2017D-31Mar2018-v1/MINIAOD',       ['label=DoubleEG']]
samples['DoubleMuon_Run2017D-31Mar2018-v1']     = ['/DoubleMuon/Run2017D-31Mar2018-v1/MINIAOD',     ['label=DoubleMuon']]
samples['MuonEG_Run2017D-31Mar2018-v1']         = ['/MuonEG/Run2017D-31Mar2018-v1/MINIAOD',         ['label=MuonEG']]
samples['SingleElectron_Run2017D-31Mar2018-v1'] = ['/SingleElectron/Run2017D-31Mar2018-v1/MINIAOD', ['label=SingleElectron']]
samples['SingleMuon_Run2017D-31Mar2018-v1']     = ['/SingleMuon/Run2017D-31Mar2018-v1/MINIAOD',     ['label=SingleMuon']]

samples['DoubleEG_Run2017E-31Mar2018-v1']       = ['/DoubleEG/Run2017E-31Mar2018-v1/MINIAOD',       ['label=DoubleEG']]
samples['DoubleMuon_Run2017E-31Mar2018-v1']     = ['/DoubleMuon/Run2017E-31Mar2018-v1/MINIAOD',     ['label=DoubleMuon']]
samples['MuonEG_Run2017E-31Mar2018-v1']         = ['/MuonEG/Run2017E-31Mar2018-v1/MINIAOD',         ['label=MuonEG']]
samples['SingleElectron_Run2017E-31Mar2018-v1'] = ['/SingleElectron/Run2017E-31Mar2018-v1/MINIAOD', ['label=SingleElectron']]
samples['SingleMuon_Run2017E-31Mar2018-v1']     = ['/SingleMuon/Run2017E-31Mar2018-v1/MINIAOD',     ['label=SingleMuon']]

samples['DoubleEG_Run2017F-31Mar2018-v1']       = ['/DoubleEG/Run2017F-31Mar2018-v1/MINIAOD',       ['label=DoubleEG']]
samples['DoubleMuon_Run2017F-31Mar2018-v1']     = ['/DoubleMuon/Run2017F-31Mar2018-v1/MINIAOD',     ['label=DoubleMuon']]
samples['MuonEG_Run2017F-31Mar2018-v1']         = ['/MuonEG/Run2017F-31Mar2018-v1/MINIAOD',         ['label=MuonEG']]
samples['SingleElectron_Run2017F-31Mar2018-v1'] = ['/SingleElectron/Run2017F-31Mar2018-v1/MINIAOD', ['label=SingleElectron']]
samples['SingleMuon_Run2017F-31Mar2018-v1']     = ['/SingleMuon/Run2017F-31Mar2018-v1/MINIAOD',     ['label=SingleMuon']]


## additional python configuration param
pyCfgParams.append('globalTag=94X_dataRun2_v6')

# Luminosity: 41.29 /fb
config.Data.lumiMask       = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions17/13TeV/Final/Cert_294927-306462_13TeV_PromptReco_Collisions17_JSON.txt'
config.Data.splitting      = 'LumiBased'
config.Data.unitsPerJob    = 15
config.Data.outLFNDirBase  = '/store/user/spak/CP_Violation/data/test'

#   ------------------------------
#     dataset | from run | to run
#   ----------+----------+--------
#    Run2016B |   297046 | 299329
#    Run2016C |   299368 | 302029
#    Run2016D |   302030 | 303434
#    Run2016E |   303824 | 304797
#    Run2016F |   305040 | 306462
#   ------------------------------
#
# https://twiki.cern.ch/twiki/bin/view/CMS/PdmV2017Analysis
