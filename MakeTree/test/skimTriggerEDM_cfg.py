import FWCore.ParameterSet.Config as cms

process = cms.Process("SKIM")
# initialize  MessageLogger
process.load("FWCore.MessageService.MessageLogger_cfi")


from AlphaTHLT.MakeTree.samples.PHY1474_STV4_745_PU30bx50_26Jun15_cfi import *

selectedSample = DYJets

process.source = cms.Source("PoolSource",
    fileNames = selectedSample.files
)


process.load('Configuration.StandardSequences.EndOfProcess_cff')

process.MINIAODSIMoutput = cms.OutputModule("PoolOutputModule",
    compressionLevel = cms.untracked.int32(4),
    compressionAlgorithm = cms.untracked.string('LZMA'),
    eventAutoFlushCompressedSize = cms.untracked.int32(15728640),
#    outputCommands = process.MINIAODSIMEventContent.outputCommands,
    outputCommands  = cms.untracked.vstring('drop *', 'keep edmTriggerResults_*_*_*',),
    fileName        = selectedSample.name,
    dataset         = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('')
    ),
    dropMetaData = cms.untracked.string('ALL'),
    fastCloning = cms.untracked.bool(False),
    overrideInputFileSplitLevels = cms.untracked.bool(True)
)


# Path and EndPath definitions 
process.endjob_step = cms.EndPath(process.endOfProcess)
process.endpath = cms.EndPath(process.MINIAODSIMoutput)
