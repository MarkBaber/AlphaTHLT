# Auto generated configuration file
# using: 
# Revision: 1.20 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --conditions POSTLS162_V2::All -n 10 --eventcontent HLTDEBUG -s DIGI,L1,HLT --datatier GEN-SIM-DIGI-RAW --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry Extended2015 --magField 38T_PostLS1 --no_exec --filein /store/mc/Fall13dr/VBF_HToInv_M-125_13TeV_powheg-pythia6/GEN-SIM-RAW/tsg_PU20bx25_POSTLS162_V2-v1/00000/001A6D13-E16C-E311-BE56-00266CFAE7C4.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT2')

################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(10)
)

# process.options = cms.untracked.PSet(
#     wantSummary = cms.untracked.bool( True )
# )
################################################################


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('AlphaTHLT.ReRunHLT.HLT_AlphaT_cff')

# For producing L1 Extra objects 
process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cff")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(200)
)

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/243D919C-8C77-E311-9364-002618943810.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/2E02D87C-8C77-E311-9B7E-0030486790B8.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/3886AA7F-5677-E311-8A87-0025905A611C.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/621BA364-8C77-E311-8DCC-003048678FF6.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/72CF499A-8C77-E311-B856-0030486790B8.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/7476B7ED-5E77-E311-8413-003048678A78.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/8403B991-4E77-E311-8244-003048FFCB9E.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/A80DFA9E-5E77-E311-B5C2-0025905A48E4.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/C471CC8E-5D77-E311-9654-003048678A78.root',
      '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/F2AFED81-5A77-E311-99E1-0025905964BE.root',

),
    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/0012EEAC-F577-E311-9AC6-003048FFCB8C.root'),
#    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/VBF_HToInv_M-125_13TeV_powheg-pythia6/GEN-SIM-RAW/tsg_PU20bx25_POSTLS162_V2-v1/00000/001A6D13-E16C-E311-BE56-00266CFAE7C4.root')
)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.20 $'),
    annotation = cms.untracked.string('step2 nevts:10'),
    name = cms.untracked.string('Applications')
)

# ================================================================================
# Additional reconstruction
# ================================================================================


# ------------------------------------------------------------
# RECOJets                                             
# ------------------------------------------------------------

# process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")

# from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
# from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
# from JetMETCorrections.Configuration.DefaultJEC_cff import *

# ak5CaloL1Fastjet.srcRho = 'kt6CaloJetsCentral:rho'
# ak5CaloL1FastL2L3 = ak5CaloL1Fastjet.clone()

# process.ak5CaloJetsL1FastL2L3 = cms.EDProducer(
#     'CaloJetCorrectionProducer',
#     src         = cms.InputTag('ak5CaloJets'),
#     correctors  = cms.vstring('ak5CaloL1FastL2L3')
#     )


# process.ak5PFJetsL1FastL2L3 = cms.EDProducer(
#     'PFJetCorrectionProducer',
#     src         = cms.InputTag('ak5PFJets'),
#     correctors  = cms.vstring('ak5PFL1FastL2L3')
#     )


# ------------------------------------------------------------
# GenJets                                             
# ------------------------------------------------------------
process.load('RecoJets.Configuration.GenJetParticles_cff')
process.load("RecoJets.JetProducers.ak5GenJets_cfi")

# GenJets - NoMuNoNu                                                                                                                       
process.ak5NoMuNoNuGenJets = process.ak5GenJets.clone(
    src = cms.InputTag("genParticlesForJetsNoMuNoNu")
    )
process.ak4NoMuNoNuGenJets = process.ak5GenJets.clone(
    src = cms.InputTag("genParticlesForJetsNoMuNoNu"),
    rParam = cms.double(0.4)
    )

# GenJets - NoNu                                   
process.ak5NoNuGenJets =  process.ak5GenJets.clone(
    src = cms.InputTag("genParticlesForJetsNoNu"),
    )
process.ak4NoNuGenJets = process.ak5GenJets.clone(
    src = cms.InputTag("genParticlesForJetsNoNu"),
    rParam = cms.double(0.4)
    )
process.ak4GenJets = process.ak5GenJets.clone(
    rParam = cms.double(0.4)
    )


# ================================================================================


process.p = cms.Path(
  process.l1extraParticles

  # RECO jets
  #*process.ak5CaloJetsL1FastL2L3
  #*process.ak5PFJetsL1FastL2L3
  
  # GenJets                                                                                                                               
  *process.genParticlesForJetsNoMuNoNu
  *process.ak4NoMuNoNuGenJets
  *process.ak5NoMuNoNuGenJets
  
  *process.genParticlesForJetsNoNu
  *process.ak4NoNuGenJets
  *process.ak5NoNuGenJets
  
  *process.genParticlesForJets
  *process.ak4GenJets
  *process.ak5GenJets

)

# Output definition
process.HLTDEBUGoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring( 
      'drop *',
      'keep *_l1extraParticles_*_*',
      'keep *_genMetCalo_*_*',
      'keep *_genMetCaloAndNonPrompt_*_*',
      'keep *_genMetTrue_*_*',

      'keep *_*GenJets_*_HLT2',
      'keep *_ak5PFJets*_*_*',
      'keep *_ak5CaloJets*_*_*',
      'keep *_*_Rho_*',
      'keep PileupSummaryInfos_*_*_*',


      'keep *_hltL1GtObjectMap_*_HLT2',
      'keep FEDRawDataCollection_rawDataCollector_*_HLT2',
      'keep FEDRawDataCollection_source_*_HLT2',
      'keep edmTriggerResults_*_*_HLT2',
      'keep triggerTriggerEvent_*_*_HLT2' ),
    fileName = cms.untracked.string('file:hltReRunResults.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS162_V2::All', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup_GRun', '')


# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.HLTDEBUGoutput_step = cms.EndPath(process.HLTDEBUGoutput)

# Schedule definition
process.schedule = cms.Schedule(process.digitisation_step,process.L1simulation_step,process.p)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.HLTDEBUGoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

#call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
process = customizeHLTforMC(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

#call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1(process)

# End of customisation functions
