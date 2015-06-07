from WMCore.Configuration import Configuration
from collections import OrderedDict

prodTag = "07Jun15"
sampleN =  0  # 0 - 8
jobName = "PHY1474_STV4_742_PU20bx25"
signal  = True
datasets= {}

if signal == True:
        datasets = {'WToMuNu'   :'/WToMuNu_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/GEN-SIM-RAW',
                    'DYToMuMu'  :'/DYToMuMu_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/GEN-SIM-RAW',
                    # 'WToENu'   :'/WToENu_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/GEN-SIM-RAW',
                    # 'DYToEE'  :'/DYToEE_Tune4C_13TeV-pythia8/Phys14DR-AVE20BX25_tsg_PHYS14_25_V3-v1/GEN-SIM-RAW',
                    }

else:

    datasets = {'QCD30to50'   :'/QCD_Pt-30to50_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v2/GEN-SIM-RAW',  
                'QCD50to80'   :'/QCD_Pt-50to80_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',  
                'QCD80to120'  :'/QCD_Pt-80to120_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW', 
                'QCD120to170' :'/QCD_Pt-120to170_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',
                'QCD170to300' :'/QCD_Pt-170to300_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',
                'QCD300to470' :'/QCD_Pt-300to470_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',
                'QCD470to600' :'/QCD_Pt-470to600_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',
                'QCD600to800' :'/QCD_Pt-600to800_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',
                'QCD800to1000':'/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/Phys14DR-AVE20BX25_tsg_castor_PHYS14_25_V3-v1/GEN-SIM-RAW',
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
config.JobType.psetName = '../test/AlphaT_HLT_742V2_25ns_PHYS14.py'
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
