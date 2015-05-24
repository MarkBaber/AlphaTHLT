from WMCore.Configuration import Configuration

prodTag = "24May15"
sampleN = 0   # 0 - 10
jobName = "FALL1374_50V0_742_PU40bx50"


datasets = ['/TT_Tune4C_13TeV-pythia8-tauola/Fall13dr-tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/DYJetsToLL_M-50_13TeV-madgraph-pythia8/Fall13dr-tsg_PU40bx50_POSTLS162_V1-v1/GEN-SIM-RAW',
            '/QCD_Pt-30to50_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-50to80_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-80to120_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-120to170_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-170to300_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-300to470_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            '/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Fall13dr-castor_tsg_PU40bx50_POSTLS162_V2-v1/GEN-SIM-RAW',
            ]
labels   = ['TT',
            'DYJetsToLL',
            'QCD30to50',
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
config.JobType.psetName = '../test/AlphaT_HLT_742V3_50ns.py'
config.JobType.pluginName = 'Analysis'
config.JobType.outputFiles = ['hltReRunResults.root']
config.section_('Data')
config.Data.outLFN   = '/store/user/mbaber/' + jobName + '_' + prodTag
config.Data.inputDBS = 'global'
config.Data.inputDataset = dataset
config.Data.publication  = False
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.totalUnits = 1
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_UK_London_IC'
