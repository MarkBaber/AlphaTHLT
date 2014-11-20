import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("TreeMaker")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 10
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(1000)
)

# --------------------------------------------------------------------------------

# Select bx to process
bx = "50ns"



if (bx == "25ns"):
    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V1A_PU40bx25_12Nov14_cfi import * # 25ns, HCAL fix
elif (bx == "50ns"):
    from AlphaTHLT.MakeTree.samples.AlphaTHLT_MCRUN2_72_V2A_PU40bx50_12Nov14_cfi import * # 50ns, HCAL fix
else:
    print "Error: Bunch spacing '", bx, "' not recognised\n"
    exit(0)

# Samples:
samples = [QCD30to50,  QCD50to80, QCD80to120, QCD120to170, QCD170to300, QCD300to470,  QCD470to600, QCD600to800, QCD800to1000, 
           TTbar, DYJets ]
           #T2cc_250_210, T2tt_500_250, T2tt_300_200 ]
#samples = [T1bbbb_2J_mGl_1000_mLSP_900, T2tt_2J_mStop_850_mLSP_100, T2tt_2J_mStop_425_mLSP_325, T2tt_2J_mStop_500_mLSP_325, T1tttt_2J_mGl_1200_mLSP_800 ]

selectedSample = samples[10]


# --------------------------------------------------------------------------------

process.source = cms.Source ("PoolSource",


                             fileNames = selectedSample.files,

#                             fileNames = cms.untracked.vstring( 'file:/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_2_1_patch2/src/AlphaTHLT/ReRunHLT/test/hltReRunResults.root' ), 
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

if bx   == "25ns":
    GT = 'MCRUN2_72_V1A::All'
elif bx == "50ns":
    GT = 'MCRUN2_72_V2A::All'

print "BX = ", bx, ", using GT:", GT
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
