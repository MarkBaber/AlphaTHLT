from WMCore.Configuration import Configuration

prodTag = "24May15"
sampleN = 13   # 0 - 13
jobName = "PHY1474_25V4_742_PU20BX25_HCAL3"

# datasets = [
#  '/SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/AODSIM',
#  '/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/AODSIM',
#  '/SMS-T2tt_2J_mStop-850_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',
#  '/SMS-T2tt_2J_mStop-650_mLSP-325_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',
#  '/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',
#  '/SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',
#  '/SMS-T2bb_2J_mStop-600_mLSP-580_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',
#  '/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/AODSIM',
#             ]
# labels   = ['T1bbbb_2J_mGl_1000_mLSP_900',
#             'T1tttt_2J_mGl_1200_mLSP_800',
#             'T2tt_2J_mStop_850_mLSP_100',
#             'T2tt_2J_mStop_650_mLSP_325',
#             'T2tt_2J_mStop_500_mLSP_325',
#             'T2tt_2J_mStop_425_mLSP_325',
#             'T2tt_2J_mStop_600_mLSP_580',
#             'T2tt_2J_mStop_600_mLSP_550',
#             ]


datasets = [
 '/SMS-T1bbbb_2J_mGl-1000_mLSP-900_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/GEN-SIM-RAW',
 '/SMS-T1bbbb_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/GEN-SIM-RAW',
 '/SMS-T1qqqq_2J_mGl-1000_mLSP-800_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/GEN-SIM-RAW',
 '/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/GEN-SIM-RAW',
 '/SMS-T1tttt_2J_mGl-1200_mLSP-800_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/GEN-SIM-RAW',
 '/SMS-T1tttt_2J_mGl-1500_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v2/GEN-SIM-RAW',
 '/SMS-T2bb_2J_mStop-600_mLSP-580_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2bb_2J_mStop-900_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2qq_2J_mStop-1200_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2tt_2J_mStop-425_mLSP-325_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2tt_2J_mStop-650_mLSP-325_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
 '/SMS-T2tt_2J_mStop-850_mLSP-100_Tune4C_13TeV-madgraph-tauola/Spring14dr-PU20bx25_POSTLS170_V5-v1/GEN-SIM-RAW',
]
labels   = ['T1bbbb_2J_mGl_1000_mLSP_900',
            'T1bbbb_2J_mGl_1500_mLSP_100',
            'T1qqqq_2J_mGl_1000_mLSP_800',
            'T1qqqq_2J_mGl_1400_mLSP_100',
            'T1tttt_2J_mGl_1200_mLSP_800',
            'T1tttt_2J_mGl_1500_mLSP_100',
            'T2bb_2J_mStop_600_mLSP_580',
            'T2bb_2J_mStop_900_mLSP_100',
            'T2qq_2J_mStop_1200_mLSP_100',
            'T2qq_2J_mStop_600_mLSP_550',
            'T2tt_2J_mStop_425_mLSP_325',
            'T2tt_2J_mStop_500_mLSP_325',
            'T2tt_2J_mStop_650_mLSP_325',
            'T2tt_2J_mStop_850_mLSP_100',

            ]


dataset = datasets[ sampleN ]
label   = labels[ sampleN ]

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
config.Data.outLFNDirBase = '/store/user/mbaber/' + jobName + '_' + prodTag
#config.Data.useParent    = True
config.Data.inputDBS = 'global'
config.Data.inputDataset = dataset
config.Data.publication  = False
config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 5
#config.Data.ignoreLocality = True  # Use AAA

#config.Data.totalUnits = 1
config.section_('User')
config.section_('Site')
config.Site.storageSite = 'T2_UK_London_IC'
#config.Site.whitelist = ["T2_UK*"] # To ensure AAA is in same region (samples at RAL)
