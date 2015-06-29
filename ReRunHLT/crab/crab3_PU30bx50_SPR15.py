from WMCore.Configuration import Configuration

prodTag = "13Jun15"
sampleN = 5   # 0 - 8
jobName = "SPR1574_STV1_752_PU30bx50"


datasets = ['/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_30_BX_50ns_tsg_MCRUN2_74_V6-v1/GEN-SIM-RAW',
            ]
labels   = ['QCD30to50',
            'QCD50to80',
            'QCD80to120',
            'QCD120to170',
            'QCD170to300',
            'QCD300to470',
            'QCD470to600',
            'QCD600to800',
            'QCD800to1000',
            ]

dataset = datasets[ sampleN ]
label   = labels[ sampleN ]

config = Configuration()
config.section_('General')
config.General.workArea = jobName + '_' + prodTag
config.General.transferOutputs = True
config.General.requestName = label
config.section_('JobType')
config.JobType.psetName = '../test/AlphaT_HLT_745_50ns_SPRING15.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['hltReRunResults.root']
config.section_('Data')
config.Data.outLFNDirBase   = '/store/user/mbaber/' + jobName + '_' + prodTag
config.Data.inputDBS        = 'global'
config.Data.inputDataset = dataset
config.Data.publication  = False
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 1
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_UK_London_IC'
