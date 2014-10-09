import FWCore.ParameterSet.Config as cms

process = cms.Process('GenSkim')

################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(5)
)

# process.options = cms.untracked.PSet(
#     wantSummary = cms.untracked.bool( True )
# )
################################################################

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')



# Input source
process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('root://xrootd.unl.edu//store/mc/Fall13dr/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/0016CB18-D9A6-E311-BC1A-003048C68A9C.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications')
)




# ================================================================================
# Additional reconstruction
# ================================================================================

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *  ", # this is the default
#    "drop pdgId = {Z0} & status = 2"
    "keep++ pdgId = {e+}   & status = 1",
    "keep++ pdgId = {e-}   & status = 1",
    "keep++ pdgId = {mu+}  & status = 1",
    "keep++ pdgId = {mu-}  & status = 1",
    # "keep++ pdgId = {tau+} & status = 1",
    # "keep++ pdgId = {tau-} & status = 1",

    )
)

process.GenParticleFilter = cms.EDFilter(
    "GenParticleFilter",
    particleCands = cms.InputTag("prunedGenParticles"),
    ptMin         = cms.double(5.0),
)

process.GenParticleSkimmer = cms.Path(
    process.prunedGenParticles
#    *process.GenParticleFilter
)

process.GenParticleSkimSchedule = cms.Schedule( process.GenParticleSkimmer )


# ================================================================================


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring( 
        'drop *',
        'keep *_prunedGenParticles_*_GenSkim',
#        'keep *_*_GenParticleFilter_GenSkim',
      ),
    fileName = cms.untracked.string('file:hltReRunResults.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('AODSIM')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.GenParticleSkimSchedule)
process.schedule.extend([process.endjob_step,process.output_step])

