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
process.load('AlphaTHLT.ReRunHLT.HLT_AlphaT_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
#    input = cms.untracked.int32(10)
    input = cms.untracked.int32(-1)
)


#from AlphaTHLT.ReRunHLT.samples.samples_Signal_MC_cfi import *
#selectedSample =  T2cc_250_210 #T2tt_500_250 #T2cc_250_210


# Input source
process.source = cms.Source("PoolSource",


    fileNames = cms.untracked.vstring("root://xrootd.unl.edu//store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/02B79593-F47F-E311-8FF6-003048FFD796.root")


#     fileNames = selectedSample.files,


#    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/4AAC6F5B-9979-E311-9DED-0025905A48F0.root')
#    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx25_POSTLS162_V2-v1/00000/0002521B-7C9D-E311-874A-003048FFCB8C.root')


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

# load PostLS1 customisation 
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
process = customisePostLS1(process)


process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')

# GT 
from L1Trigger.Configuration.SimL1Emulator_cff import simGtDigis
process.simGtDigis = simGtDigis.clone()
process.simGtDigis.GmtInputTag = 'simGmtDigis'
process.simGtDigis.GctInputTag = 'caloStage1LegacyFormatDigis'
process.simGtDigis.TechnicalTriggersInputTags = cms.VInputTag( )


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

  # # UCT
  # process.L1TCaloStage1_PPFromRaw
  # +process.simGtDigis
  # +process.l1ExtraLayer2



  # GenJets                                                                                                                               
  process.genParticlesForJetsNoMuNoNu
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






# # Other statements
# from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
# #process.GlobalTag = GlobalTag(process.GlobalTag, 'PRE_LS172_V16::All', '')
# process.GlobalTag = GlobalTag(process.GlobalTag, 'PRE_LS172_V16', '')
# #process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc_GRun', '')

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


# --------------------------------------------------------------------------------

# customise the HLT menu for running on MC
from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC
process = customizeHLTforMC(process)

# load PostLS1 customisation
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1
process = customisePostLS1(process)

# override the GlobalTag, connection string and pfnPrefix 
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'auto:run2_mc_GRun')
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing 
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )

# customize the L1 emulator to run customiseL1EmulatorFromRaw with HLT to switchToSimStage1Digis 
process.load( 'Configuration.StandardSequences.RawToDigi_cff' )
process.load( 'Configuration.StandardSequences.SimL1Emulator_cff' )
import L1Trigger.Configuration.L1Trigger_custom 
# 

# 2015 Run2 emulator
import L1Trigger.L1TCalorimeter.L1TCaloStage1_customForHLT
process = L1Trigger.L1TCalorimeter.L1TCaloStage1_customForHLT.customiseL1EmulatorFromRaw( process )

#
process = L1Trigger.Configuration.L1Trigger_custom.customiseResetPrescalesAndMasks( process )
# customize the HLT to use the emulated results
import HLTrigger.Configuration.customizeHLTforL1Emulator
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimStage1Digis( process )



# --------------------------------------------------------------------------------


