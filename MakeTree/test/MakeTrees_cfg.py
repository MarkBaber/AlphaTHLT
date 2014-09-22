import FWCore.ParameterSet.Config as cms
import os

from FWCore.ParameterSet.VarParsing import VarParsing
process = cms.Process("TreeMaker")


process.load("FWCore.MessageLogger.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
    )


# --------------------------------------------------------------------------------

from AlphaTHLT.MakeTree.samples.samples_21Sep14_cfi import *

# TTbar       = DataSets_cfi.TTbar
# NeutrinoGun = DataSets_cfi.NeutrinoGun
# Muon        = DataSets_cfi.Muon
# MonoJet     = DataSets_cfi.MonoJet
# MonoJet1000 = DataSets_cfi.MonoJet1000

selectedSample = TTbar

# --------------------------------------------------------------------------------

process.source = cms.Source ("PoolSource",
                             #fileNames = selectedSample.files
                             fileNames = cms.untracked.vstring(

        'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/AlphaTHLT_POSTLS162_V2_PU40bx25_18Sep14/TT_Tune4C_13TeV-pythia8-tauola/hltReRunResults_159_1_hJl.root',

        # 'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mcitron/skims/140917_ngun50_alphatL1/l1_skim_2014-09-17_1466_1_SGk.root',
        #        'file:///home/hep/mc3909/AlphaT/CMSSW_7_1_0_pre8/src/AlphaTL1/MakeL1Trees/test/l1_skim_2014-09-17.root'

        )
)
process.source.duplicateCheckMode = cms.untracked.string('noDuplicateCheck')

process.TFileService = cms.Service("TFileService",
                                   fileName = selectedSample.name 
) 

process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.GlobalTag.globaltag = 'POSTLS161_V2::All'


process.load("AlphaTHLT.MakeTree.MakeTrees_cfi")



# process.MakeTrees = cms.EDAnalyzer("MakeTrees",

# )
# process.HLTJetProducer = cms.EDProducer('HLTJetProducer',
#                                         HLTResults = cms.untracked.InputTag("TriggerResults::HLT2"),
#                                         JetMinPT   = cms.double( 20. ),
# )




process.p1 = cms.Path(
      process.MakeTrees
      )
