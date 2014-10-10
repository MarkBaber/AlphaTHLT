import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("TreeMaker")

process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
#    input = cms.untracked.int32(10)
    )


# --------------------------------------------------------------------------------

from AlphaTHLT.MakeTree.samples.samples_PU40_10Oct14_cfi import *

# TTbar, DYJets, 
# QCD30to50, QCD50to80, QCD80to120, QCD120to170, QCD170to300, QCD300to470,  QCD470to600, QCD600to800, QCD800to1000
# NuGun
# T2cc_250_210, T2tt_500_250
selectedSample = QCD30to50

# --------------------------------------------------------------------------------

process.source = cms.Source ("PoolSource",

                             fileNames = selectedSample.files,
                             #fileNames = cms.untracked.vstring( 'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/AlphaTHLT_PRE_LS172_V16_PU40bx25_10Oct14/QCD_Pt-30to50_Tune4C_13TeV_pythia8/hltReRunResults_302_1_UtZ.root' )

)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
                                   fileName = selectedSample.name 
) 

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
#process.GlobalTag.globaltag = 'POSTLS161_V2::All'
process.GlobalTag.globaltag = 'PRE_LS172_V16::All'



process.load("AlphaTHLT.MakeTree.MakeTrees_cfi")
# process.MakeTrees = cms.EDAnalyzer("MakeTrees",
# )

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

#process.p += process.HLTJetProducer
#process.e = cms.EndPath(process.out)


