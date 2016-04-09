from WMCore.Configuration import Configuration

prodTag = "07Apr16"
sampleN = 0   # 0 - 8
jobName = "Run2015D"
jobTag  = jobName + '_' + prodTag


datasets = ['/HLTPhysics3/Run2015D-v1/RAW',
            '/SingleMuon/Run2015D-v1/RAW',
            '/SingleElectron/Run2015D-v1/RAW',
            '/ZeroBias1/Run2015D-v1/RAW',
            '/ZeroBias2/Run2015D-v1/RAW',
            '/ZeroBias3/Run2015D-v1/RAW',
            '/ZeroBias4/Run2015D-v1/RAW',
            ]
labels   = ['HLTPhysics3',
            'SingleMuon',
            'SingleElectron',
            'ZeroBias1',
            'ZeroBias2',
            'ZeroBias3',
            'ZeroBias4',
            ]

dataset = datasets[ sampleN ]
label   = labels[ sampleN ]

config = Configuration()
config.section_('General')
config.General.workArea        = jobTag
config.General.transferOutputs = True
config.General.requestName     = label
config.section_('JobType')
config.JobType.psetName = '../test/hlt_data_alphaT.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['hltReRunResults.root']
config.section_('Data')
config.Data.outLFNDirBase   = '/store/user/mbaber/' + jobTag
config.Data.inputDBS        = 'global'
config.Data.inputDataset = dataset
config.Data.publication  = False
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
config.Data.lumiMask    = 'Selected_Trigger_Runs.JSON'
#config.Data.totalUnits = 1
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_UK_London_IC'
