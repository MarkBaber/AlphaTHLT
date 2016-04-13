from WMCore.Configuration import Configuration

prodTag = "13Apr16"
sampleN = 0   # 0 - 8
jobName = "MC_76X_RECOTest"
jobTag  = jobName + '_' + prodTag


datasets = ['/SMS-T2tt_mStop-500_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW',
            ]
labels   = ['T2tt_mStop_500_mLSP_325',
            ]

dataset = datasets[ sampleN ]
label   = labels[ sampleN ]

config = Configuration()
config.section_('General')
config.General.workArea        = jobTag
config.General.transferOutputs = True
config.General.requestName     = label
config.section_('JobType')
config.JobType.psetName = '../test/hlt_mc_alphaT.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['hltReRunResults.root']
config.section_('Data')
config.Data.outLFNDirBase   = '/store/user/mbaber/' + jobTag
config.Data.inputDBS        = 'global'
#config.Data.inputDataset          = dataset
config.Data.inputDataset          = '/SMS-T2tt_mStop-500_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/AODSIM'
config.Data.secondaryInputDataset = '/SMS-T2tt_mStop-500_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/RunIIFall15DR76-25nsFlat10to50NzshcalRaw_76X_mcRun2_asymptotic_v12-v1/GEN-SIM-RAW'
config.Data.publication  = False
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 1
# DEBUGGING
config.Data.totalUnits  = 1
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_UK_London_IC'
