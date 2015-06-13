from WMCore.Configuration import Configuration

prodTag = "26May15v2"
sampleN =  6  # 0 - 8
jobName = "PHY1474_STV4_745_PU30bx50"
signal = True

datasets = {}
if signal:
    datasets = {'T1bbbb_2J_mGl_1000_mLSP_900':'/SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'T1qqqq_2J_mGl_1000_mLSP_800':'/SMS-T1qqqq_2J_mGl-1000_mLSP-800_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v3/GEN-SIM-RAW',
                'T1tttt_2J_mGl_1200_mLSP_800':'/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/GEN-SIM-RAW',
                'T1tttt_2J_mGl_1500_mLSP_100':'/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'T2tt_2J_mStop_650_mLSP_325' :'/SMS-T2tt_2J_mStop-650_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/GEN-SIM-RAW',
                'T2tt_2J_mStop_500_mLSP_325' :'/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v2/GEN-SIM-RAW',
                'T2tt_2J_mStop_425_mLSP_325' :'/SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola/Phys14DR-AVE30BX50_tsg_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                }
else:
    datasets = {'QCD30to50'   :'/QCD_Pt-30to50_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v2/GEN-SIM-RAW',
                'QCD50to80'   :'/QCD_Pt-50to80_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD80to120'  :'/QCD_Pt-80to120_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD120to170' :'/QCD_Pt-120to170_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD170to300' :'/QCD_Pt-170to300_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD300to470' :'/QCD_Pt-300to470_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD470to600' :'/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD600to800' :'/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                'QCD800to1000':'/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Phys14DR-AVE30BX50_tsg_castor_PHYS14_ST_V1-v1/GEN-SIM-RAW',
                }

dataset = datasets.values()[ sampleN ]
label   = datasets.keys()  [ sampleN ]

# ------------------------------------------------------------------------------------------------------------------------

config = Configuration()
config.section_('General')
config.General.workArea = jobName + '_' + prodTag
config.General.transferOutputs = True
config.General.requestName = label
config.section_('JobType')
config.JobType.psetName = '../test/AlphaT_HLT_745_PHYS14_50ns.py'
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
