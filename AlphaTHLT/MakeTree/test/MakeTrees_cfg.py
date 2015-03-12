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
bx = "25ns"
#bx = "50nsPU1"
#bx = "50ns"
#bx = "20PU25ns"


if (bx == "25ns"):
#    from AlphaTHLT.MakeTree.samples.PU40bx25_MfHCALM2_10Mar15.TTbar.TTbar_cfi import *
#    from AlphaTHLT.MakeTree.samples.MCRUN2_72_V3A_PU40bx25_MfHCALM2.QCD600to800_cfi import *
    from AlphaTHLT.MakeTree.samples.MCRUN2_72_V3A_PU40bx25_MfHCALM2_cfi import *
#    from AlphaTHLT.MakeTree.samples.PU40bx25_MfHCALM2_10Mar15.10Mar15_PU40bx25_MfHCALM2_cfi import *
    # from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V1A_PU40bx25_12Nov14_cfi import * # 25ns, HCAL fix
    # from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V1A_PU20bx25_17Nov14_cfi import * # 25ns, HCAL fix, signal
elif (bx == "50ns"):
    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V2A_PU40bx50_12Nov14_cfi import * # 50ns, HCAL fix
elif (bx == "50nsPU1"):
#    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V4A_PU1bx50_29Jan15_cfi  import * # PU1
    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V4A_PU1bx50_16Feb15_cfi  import * # PU1
elif (bx == "20PU25ns"):
    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V1A_PU20bx25_27Nov14_cfi import * # 25ns, PU20
    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V3A_PU20bx25_13Jan15_cfi import * # 25ns, PU20, RECO
else:
    print "Error: Bunch spacing '", bx, "' not recognised\n"
    exit(0)

# Samples:
# ----------------------------------------
samples = []

if (bx == "20PU25ns"):
    # samples = [QCD30to50,    # 0
    #            QCD50to80,    # 1
    #            QCD80to120,   # 2
    #            QCD120to170,  # 3
    #            QCD170to300,  # 4
    #            QCD300to470,  # 5 
    #            QCD470to600,  # 6
    #            QCD600to800,  # 7
    #            #QCD800to1000, # 8
    #            TTbar,        # 8
    #            DYJets ]      # 9
    samples = [T2tt_2J_mStop_425_mLSP_325, # 1M events
               T2bb_2J_mStop_600_mLSP_580,
               T2tt_2J_mStop_650_mLSP_325,
               T2qq_2J_mStop_600_mLSP_550,
               T2tt_2J_mStop_500_mLSP_325,
               T2tt_2J_mStop_850_mLSP_100
               ]
elif (bx == "25ns"):
    samples = [
                QCD30to50,    # 0
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
elif (bx == "50nsPU1"):
    samples = [QCD30to50,    # 0
               # QCD50to80,    # 1
               # QCD80to120,   # 2
               # QCD120to170,  # 3
               QCD170to300,  # 4
               QCD300to470,  # 5 
               # QCD470to600,  # 6
               # QCD600to800,  # 7
               #QCD800to1000] # 8
               ]
selectedSample = samples[10]


# --------------------------------------------------------------------------------

process.source = cms.Source ("PoolSource",

                             fileNames = selectedSample.files,

#                             fileNames = cms.untracked.vstring( 'file:/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_3_0/src/AlphaTHLT/ReRunHLT/test/hltReRunResults.root' ), 
                             #fileNames = cms.untracked.vstring( 'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/AlphaTHLT_PRE_LS172_V16_PU40bx25_10Oct14/QCD_Pt-30to50_Tune4C_13TeV_pythia8/hltReRunResults_302_1_UtZ.root' )

)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
                                   fileName = selectedSample.name 
) 



# --------------------------------------------------------------------------------
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# override the GlobalTag, connection string and pfnPrefix                                                                                  
GT = ""

if (bx == "25ns" or bx == "20PU25ns"):
    GT = 'MCRUN2_72_V3A::All'
elif (bx == "50ns" or bx == "50nsPU1"):
    GT = 'MCRUN2_72_V4A::All'

print "\nProcessing sample :\t", selectedSample.name, "\nBeam BX scenario  :\t", bx, "\nGlobaltag         :\t", GT, "\n\n"
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = GT)
    process.GlobalTag.connect = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing                                 
    process.GlobalTag.RefreshEachRun   = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

# --------------------------------------------------------------------------------



process.load("AlphaTHLT.MakeTree.MakeTrees_cfi")

process.HLTJetProducer = cms.EDProducer('HLTJetProducer',
                                        HLTResults = cms.untracked.InputTag("TriggerResults::HLT2"),
                                        JetMinPT   = cms.double( 20. ),
                                        HLTAk4Calojets = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed"),
)



# Output module
process.out = cms.OutputModule("PoolOutputModule",
                               fileName = cms.untracked.string('file:patReRunResults.root'),
                               SelectEvents    = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
                               outputCommands  = cms.untracked.vstring('drop *')
                              )

process.p = cms.Path()
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTrigger( process, path = 'p', hltProcess = 'HLT2') #,outputModule = 'out')



process.p1 = cms.Path(
      process.MakeTrees
      )
