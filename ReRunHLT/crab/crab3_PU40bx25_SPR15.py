from WMCore.Configuration import Configuration

prodTag = "25Jul15"
sampleN = 0   # 0 - 8
jobName = "74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3"


datasets = ['/QCD_Pt_30to50_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',    #4,970,341
            '/QCD_Pt_50to80_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v2/GEN-SIM-RAW',    #4,978,005
            '/QCD_Pt_80to120_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',   #996,507
            '/QCD_Pt_120to170_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v2/GEN-SIM-RAW',  #981,327
            '/QCD_Pt_170to300_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',  #495,618
            '/QCD_Pt_300to470_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',  #492,994
            '/QCD_Pt_470to600_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW',  #496,179
            '/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v2/GEN-SIM-RAW',  #199,880
            '/QCD_Pt_800to1000_TuneCUETP8M1_13TeV_pythia8/RunIISpring15Digi74-AVE_40_BX_25ns_tsg_MCRUN2_74_V7-v1/GEN-SIM-RAW', #199,872
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
config.JobType.psetName = '../test/AlphaT_HLT_742V2_25ns.py'
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
