# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: AlphaT_ReRunHLT --conditions PRE_LS172_V16::All -n 10 --eventcontent FEVTDEBUGHLT -s HLT --datatier GEN-SIM-DIGI --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry Extended2015 --magField 38T_PostLS1 --no_exec --filein file:step1.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

# Load menu
#process.load('HLTrigger.Configuration.HLT_GRun_cff')
#process.load('AlphaTHLT.ReRunHLT.hlt_stage1_newL1Menu_Slim_cff')
#process.load('AlphaTHLT.ReRunHLT.hlt_stage1_cff')
process.load('AlphaTHLT.ReRunHLT.hlt_stage1_AlphaT_cff')

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(3)
)

# Input source
process.source = cms.Source("PoolSource",


#    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/4AAC6F5B-9979-E311-9DED-0025905A48F0.root')
    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx25_POSTLS162_V2-v1/00000/0002521B-7C9D-E311-874A-003048FFCB8C.root')


)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('AlphaT_ReRunHLT nevts:10'),
    name = cms.untracked.string('Applications')
)

# # Output definition

# process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
#     splitLevel = cms.untracked.int32(0),
#     eventAutoFlushCompressedSize = cms.untracked.int32(1048576),
#     outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
#     fileName = cms.untracked.string('file:step2.root'),
#     dataset = cms.untracked.PSet(
#         filterName = cms.untracked.string(''),
#         dataTier = cms.untracked.string('GEN-SIM-DIGI')
#     )
# )



# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring( 
       'drop *',

        # NVTX
        'keep PileupSummaryInfos_*_*_*',

        # GEN
        'keep GenEventInfoProduct_*_*_*',
        'keep *_prunedGenParticles_*_*',
        #        'keep *_*_GenParticleFilter_*', 
        'keep recoGenJets_*GenJets_*_HLT2',

        'keep *_genMetCalo_*_*',
        'keep *_genMetCaloAndNonPrompt_*_*',
        'keep *_genMetTrue_*_*',

        

       #   'keep *_hlt*_*_HLT2',


       #  'keep *_*_*_HLT2',
       # 'drop *_*_*Digi*_',
       # 'drop *_mix_*_*',

# gen particle skimmer

# BXVector<l1t::Jet>                    "caloStage1FinalDigis"      ""                "HLT2"

        # UCT
       'keep *BXVector*_*_*_*',
       'keep *_caloStage1FinalDigis_*_*',


      # UCT
      'keep l1extra*_*_*_*',
      # HLT
      'keep recoPFJets_hltAK4PFJetsCorrected*_*_*',
      'keep recoCaloJets_hltAK4CaloJetsCorrected*_*_*',


        'keep triggerTriggerFilterObjectWithRefs_*_*_HLT2',
        'keep recoMETs_*_*_HLT2',
        'keep recoCaloMETs_*_*_HLT2',


# #       # HLT jet triggers
#       'keep recoCaloJets_hltCaloJetL1FastJetCorrected_*_HLT2',
#       'keep recoPFJets_hltAntiKT4PFJets*_*_HLT2',
#       'keep *_hltSingleJet20_*_*',
#       'keep *_hlt1PFJet20_*_*',
#       'keep *_hlt1PFJetNoPU20_*_*',

#       # Offline jets
#       'keep *_ak5*Jets*_*_*',
#       'drop *_ak5*Jets*_rho_*',
#       'drop *_ak5*Jets*_sigma*_*',
#       'keep *_kt6CaloJetsCentral_rho_*',
# #      'keep *_*Jets*_*_*',


      'keep *_ak4PFJets*_*_*',
      'keep *_ak4CaloJets*_*_*',

      # 'drop *_hlt*_rho_*',
      # 'drop *_*Gen*_rho_*',



      'keep *_hltL1GtObjectMap_*_HLT2',
      'keep FEDRawDataCollection_rawDataCollector_*_HLT2',
      'keep FEDRawDataCollection_source_*_HLT2',
      'keep edmTriggerResults_*_*_HLT2',
      'keep triggerTriggerEvent_*_*_HLT2' 



      ),
    fileName = cms.untracked.string('file:hltReRunResults.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        #        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    )
)


# ================================================================================
# Additional reconstruction
# ================================================================================

# ------------------------------------------------------------
# UCT
# ------------------------------------------------------------

process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')

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

process.JetProducer = cms.Path(
#  process.l1extraParticles

  # # GCT
  # process.RawToDigi
  # *process.valRctDigis
  # *process.valGctDigis
  # *process.L1Extra

  # UCT
  #  *process.L1TCaloStage1_PPFromRaw
  
  process.L1TCaloStage1_PPFromRaw
  *process.l1ExtraLayer2

  # GenJets                                                                                                                               
  *process.genParticlesForJetsNoMuNoNu
  *process.ak4NoMuNoNuGenJets
  *process.ak5NoMuNoNuGenJets
  
  *process.genParticlesForJetsNoNu
  *process.ak4NoNuGenJets
  *process.ak5NoNuGenJets
  
  # *process.genParticlesForJets
  # *process.ak4GenJets
  # *process.ak5GenJets

)

process.JetProducerSchedule = cms.Schedule( process.JetProducer )

# ------------------------------------------------------------
# GenParticles
# ------------------------------------------------------------

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")
process.prunedGenParticles = cms.EDProducer(
    "GenParticlePruner",
    src = cms.InputTag("genParticles"),
    select = cms.vstring(
    "drop  *  ", 
    "keep++ pdgId = {e+}   & status = 1",
    "keep++ pdgId = {e-}   & status = 1",
    "keep++ pdgId = {mu+}  & status = 1",
    "keep++ pdgId = {mu-}  & status = 1",
    )
)

# process.GenParticleFilter = cms.EDFilter(
#     "GenParticleFilter",
#     particleCands = cms.InputTag("prunedGenParticles"),
#     ptMin         = cms.double(5.0),
# )

process.GenParticleSkimmer = cms.Path(
    process.prunedGenParticles
#    *process.GenParticleFilter
)

process.GenParticleSkimSchedule = cms.Schedule( process.GenParticleSkimmer )

# ================================================================================






# Additional output definition






# Other statements
from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
#process.GlobalTag = GlobalTag(process.GlobalTag, 'PRE_LS172_V16::All', '')
process.GlobalTag = GlobalTag(process.GlobalTag, 'PRE_LS172_V16', '')
#process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.GenParticleSkimSchedule)
#process.schedule.extend([process.digitisation_step]) #,process.L1simulation_step])
process.schedule.extend(process.JetProducerSchedule)
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.output_step])

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

