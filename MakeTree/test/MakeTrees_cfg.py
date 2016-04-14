import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("TreeMaker")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
   input = cms.untracked.int32(-1)
#   input = cms.untracked.int32(1000)
)

# --------------------------------------------------------------------------------
# Select bx to process
#bx = "25ns" #bx = "50ns" #bx = "20PU25ns" #bx = "40PU25ns" #bx = "AVE30BX50"
#bx = "MC"
bx = "Data"
#bx = "ZB"

if   (bx == "25ns"):      from AlphaTHLT.MakeTree.samples.FALL1374_25V4_742_PU40bx25_HCAL3_24May15_cfi import * 
elif (bx == "50ns"):      from AlphaTHLT.MakeTree.samples.FALL1374_50V0_742_PU40bx50_24May15_cfi import * 
elif (bx == "AVE30BX50"): from AlphaTHLT.MakeTree.samples.PHY1474_STV4_742_PU30bx50_26May15v2_cfi import *
elif (bx == "20PU25ns"):  from AlphaTHLT.MakeTree.samples.PHY1474_STV4_742_PU20bx25_28May15_cfi import * 
elif (bx == "40PU25ns"):  from AlphaTHLT.MakeTree.samples._74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15_cfi import *
elif (bx == "MC"):        from AlphaTHLT.MakeTree.samples.MC_76X_Full_14Apr16_cfi    import *
elif (bx == "Data"):      from AlphaTHLT.MakeTree.samples.Run2015D_11Apr16_cfi       import *
elif (bx == "ZB"):        from AlphaTHLT.MakeTree.samples.Run2015D_Loose_12Apr16_cfi import *
else:  
    print "Error: Bunch spacing '", bx, "' not recognised\n"
    exit(0)

# Samples: Specify which samples to ntuplise as defined in: AlphaT/MakeTree/python/samples
# ----------------------------------------
samples = [] # Container for sample PSets

if (bx == "AVE30BX50"):
    samples = [QCD30to50,    # 0
               QCD50to80,    # 1
               QCD80to120,   # 2
               QCD120to170,  # 3
               QCD170to300,  # 4
               QCD300to470,  # 5 
               QCD470to600,  # 6
               QCD600to800,  # 7
               QCD800to1000, # 8
               # T2tt_2J_mStop_425_mLSP_325, # 1M events
               # T1tttt_2J_mGl_1500_mLSP_100,
               # T2tt_2J_mStop_650_mLSP_325,
               # T1tttt_2J_mGl_1200_mLSP_800,
               ]
elif (bx == "20PU25ns" or bx == "40PU25ns"):
    samples = [QCD30to50,    # 0
               QCD50to80,    # 1
               QCD80to120,   # 2
               QCD120to170,  # 3
               QCD170to300,  # 4
               QCD300to470,  # 5 
               QCD470to600,  # 6
               QCD600to800,  # 7
               QCD800to1000, # 8
               # T2tt_2J_mStop_425_mLSP_325, # 1M events
               # T1tttt_2J_mGl_1500_mLSP_100,
               # T2tt_2J_mStop_650_mLSP_325,
               # T1tttt_2J_mGl_1200_mLSP_800,
               ]

    # samples = [
    # T2tt_2J_mStop_425_mLSP_325 , # 1M events - 0
    # T1bbbb_2J_mGl_1000_mLSP_900,
    # #T1tttt_2J_mGl_1200_mLSP_800,
    # T2qq_2J_mStop_1200_mLSP_100,
    # T2tt_2J_mStop_650_mLSP_325 ,
    # T2qq_2J_mStop_600_mLSP_550 ,
    # T1qqqq_2J_mGl_1400_mLSP_100,
    # T2bb_2J_mStop_900_mLSP_100 ,
    # T2bb_2J_mStop_600_mLSP_580 ,
    # T2tt_2J_mStop_500_mLSP_325 ,
    # T1tttt_2J_mGl_1500_mLSP_100,
    # T2tt_2J_mStop_850_mLSP_100 ,
    # T1bbbb_2J_mGl_1500_mLSP_100,
    # T1qqqq_2J_mGl_1000_mLSP_800, # 13
    # ]

elif (bx == "25ns"):
    samples = [ QCD30to50,    # 0
                QCD50to80,    # 1
                QCD80to120,   # 2
                QCD120to170,  # 3
                QCD170to300,  # 4
                QCD300to470,  # 5 
                QCD470to600,  # 6
                QCD600to800,  # 7
                QCD800to1000, # 8
                TTbar,        # 9
                DYJets ]      # 10
elif (bx == "50ns"):
    samples = [QCD30to50,    # 0
               QCD50to80,    # 1
               QCD80to120,   # 2
               QCD120to170,  # 3
               QCD170to300,  # 4
               QCD300to470,  # 5 
               QCD470to600,  # 6
               QCD600to800,  # 7
               QCD800to1000, # 8
               TTbar,        # 9
               DYJets ]      # 10
elif (bx == "MC"):
    samples = [T2tt_2J_mStop_500_mLSP_325,
               T2tt_2J_mStop_425_mLSP_325,
               T1tttt_2J_mGluino_1200_mLSP_800,
               T1tttt_2J_mGluino_1500_mLSP_100,
               T1bbbb_2J_mGluino_1500_mLSP_100,
               T1bbbb_2J_mGluino_1000_mLSP_900,
               ]
elif (bx == "Data"):
    samples = [HLTPhysics1,
               HLTPhysics2,
               HLTPhysics3,
               HLTPhysics4
               ]
elif (bx == "ZB"):
    samples = [ZeroBias1,
               ZeroBias2,
               ZeroBias3,
               ZeroBias4,
               ]

selectedSample = samples[3]


# --------------------------------------------------------------------------------

process.source = cms.Source ("PoolSource",
#                             fileNames = cms.untracked.vstring( 'file:hltReRunResults.root' ),  # Test on local file
                             fileNames = selectedSample.files,                                  # Execute file from imported PSet
#                             fileNames = cms.untracked.vstring( 'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/MCRUN2_72_V3A_74X_PU40bx25/TT_Tune4C_13TeV-pythia8-tauola/crab_TT/150322_224937/0000/hltReRunResults_1.root')

)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
                                   fileName = selectedSample.name 
) 

# --------------------------------------------------------------------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# # override the GlobalTag, connection string and pfnPrefix                                                                                  
# GT = ""

# if (bx == "25ns" or bx == "20PU25ns" or bx == "40PU25ns"):
#     GT = 'MCRUN2_72_V3A::All' #GT = 'PHY1474_25V4'
# elif (bx == "50ns" or bx == "50nsPU1"):
#     GT = 'MCRUN2_72_V4A::All' #GT = 'FALL1374_50V0'
# elif (bx == "AVE30BX50"):
#     GT = 'MCRUN2_72_V4A::All' #GT = 'PHYS14_50_V1' #GT = 'PHY1474_STV4'

# print "\nProcessing sample :\t", selectedSample.name, "\nBeam BX scenario  :\t", bx, "\nGlobaltag         :\t", GT, "\n\n"
# if 'GlobalTag' in process.__dict__:
#     from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
#     process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = GT)
#     process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
#     process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
#     for pset in process.GlobalTag.toGet.value():
#         pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
#     # fix for multi-run processing                                 
#     process.GlobalTag.RefreshEachRun   = cms.untracked.bool( False )
#     process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )
# --------------------------------------------------------------------------------



process.load("AlphaTHLT.MakeTree.MakeTrees_cfi")

process.HLTJetProducer = cms.EDProducer('HLTJetProducer',
                                        HLTResults = cms.untracked.InputTag("TriggerResults::HLT2"),
                                        JetMinPT   = cms.double( 30. ),
                                        HLTAk4Calojets = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed"),
)

# Output module
process.out = cms.OutputModule("PoolOutputModule",
                               fileName        = cms.untracked.string('file:patReRunResults.root'),
                               SelectEvents    = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands  = cms.untracked.vstring('drop *')
                              )

process.p = cms.Path()
from PhysicsTools.PatAlgos.tools.trigTools import *
# switchOnTrigger( process, path = 'p', hltProcess = 'HLT2') #,outputModule = 'out')

#process.RemovePileUpDominatedEvents = cms.EDFilter("RemovePileUpDominatedEvents")


process.content = cms.EDAnalyzer("EventContentAnalyzer")

process.p1 = cms.Path(
      #process.RemovePileUpDominatedEvents*
      #process.content *
      process.MakeTrees
      )
