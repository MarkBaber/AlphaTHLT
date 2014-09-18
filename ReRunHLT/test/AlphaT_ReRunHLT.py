# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --conditions POSTLS162_V2::All -n 10 --eventcontent FEVTDEBUGHLT -s HLT --datatier GEN-SIM-RAW --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry Extended2015 --magField 38T_PostLS1 --no_exec --filein file:step1.root --fileout file:step2.root
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
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.Geometry.GeometryExtended2015Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

#process.load('HLTrigger.Configuration.HLT_GRun_cff')
#process.load('AlphaTHLT.ReRunHLT.HLT_AlphaT_cff')
#process.load('AlphaTHLT.ReRunHLT.HLT_GRun_cff')
process.load('AlphaTHLT.ReRunHLT.HLT_GRun_slim_cff')





# Input source
process.source = cms.Source("PoolSource",
#      secondaryFileNames = cms.untracked.vstring(
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/4AAC6F5B-9979-E311-9DED-0025905A48F0.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/4E2A7324-B979-E311-8123-00248C55CC40.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/4ED2613D-A079-E311-9270-001BFCDBD1BA.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/747CE870-B679-E311-9A50-00261894391F.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/78B3D2C6-9A79-E311-B2BB-0026189438CB.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/7A7D1028-A379-E311-8468-0025905A60F8.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/86B5449A-B879-E311-BCA1-0026189438E1.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/A649988B-BA79-E311-9B73-00261894388A.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/C4F45D34-B079-E311-92A0-003048FFCC18.root',
# # '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/D09B073E-8B79-E311-98B7-0025905A608C.root',
#                             'root://xrootd.unl.edu//store/mc/Fall13dr/TT_Tune4C_13TeV-pythia8-tauola/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/DEA553B3-B575-E311-9A10-00304867BFF2.root',
#  ),
#     fileNames = cms.untracked.vstring(
# '/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/006AD633-407A-E311-B6DB-0025905A60A0.root',
#        
# ),

    # Test locally
    # --------------------
    # secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:Neutrino_Pt-2to20_gun_tsg_PU40bx25.root')
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


# ------------------------------------------------------------
# L1Jets                                             
# ------------------------------------------------------------

# For producing L1 Extra objects 
#process.load("L1Trigger.L1ExtraFromDigis.l1extraParticles_cff")



# GCT modifications for 2015
# --------------------------------------------------------------------------------

process.load("Configuration.StandardSequences.RawToDigi_cff")

process.load('L1Trigger.Configuration.L1Extra_cff')
process.load("L1TriggerConfig.GctConfigProducers.L1GctConfig_cff")
process.load("Configuration.StandardSequences.L1HwVal_cff")

# Employ new 10 GeV seed
process.L1GctConfigProducers.JetFinderCentralJetSeed = 10.0
process.L1GctConfigProducers.JetFinderForwardJetSeed = 10.0

process.es_prefer_gct = cms.ESPrefer("L1GctConfigProducers")

process.l1extraParticles.isolatedEmSource      = cms.InputTag('valGctDigis', 'isoEm')
process.l1extraParticles.nonIsolatedEmSource   = cms.InputTag('valGctDigis', 'nonIsoEm')
process.l1extraParticles.centralJetSource      = cms.InputTag('valGctDigis', 'cenJets')
process.l1extraParticles.forwardJetSource      = cms.InputTag('valGctDigis', 'forJets')
process.l1extraParticles.tauJetSource          = cms.InputTag('valGctDigis', 'tauJets')
process.l1extraParticles.etTotalSource         = cms.InputTag('valGctDigis')
process.l1extraParticles.etHadSource           = cms.InputTag('valGctDigis')
process.l1extraParticles.etMissSource          = cms.InputTag('valGctDigis')
process.l1extraParticles.htMissSource          = cms.InputTag("valGctDigis")
process.l1extraParticles.hfRingEtSumsSource    = cms.InputTag("valGctDigis")
process.l1extraParticles.hfRingBitCountsSource = cms.InputTag("valGctDigis")


# ------------------------------------------------------------
# UCT
# ------------------------------------------------------------

process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')

# # ------------------------------------------------------------
# # RECOJets                                             
# # ------------------------------------------------------------

# process.load("Configuration.StandardSequences.RawToDigi_cff")
# process.load("Configuration.StandardSequences.Reconstruction_Data_cff")
# process.load("RecoParticleFlow.PFProducer.particleFlow_cff")

# process.load("RecoJets.Configuration.RecoJets_cff")

# from RecoJets.JetProducers.ak4CaloJets_cfi import ak4CaloJets
# from RecoJets.JetProducers.ak4PFJets_cfi   import ak4PFJets
# #from RecoJets.JetProducers.kt6CaloJets_cfi import kt6CaloJets

# # ak4 
# # -------------------- 
# process.ak4CaloJets.doPVCorrection = cms.bool( False )
# process.ak4CaloJets.doAreaFastjet  = cms.bool( True )
# process.ak4CaloJets.doRhoFastjet   = cms.bool( True )
# process.ak4CaloJets.doPUOffsetCorr = cms.bool( False )
# # NOTE: CURRENTLY RESTRICTING RHO DETERMINATION TO |eta| < 2.6 
# #process.ak4CaloJets.Rho_EtaMax     = cms.double( 2.6 )
# process.ak4CaloJets.jetPtMin       = cms.double( 0.00001 )


# # kt6 (For rho estimation) 
# # -------------------- 
# process.kt6CaloJets.doPVCorrection = cms.bool( False )
# process.kt6CaloJets.rParam         = cms.double(0.6)
# process.kt6CaloJets.doRhoFastjet   = cms.bool( True )
# process.kt6CaloJets.doAreaFastjet  = cms.bool( True )
# process.kt6CaloJets.doPUOffsetCorr = cms.bool( False )
# # NOTE: CURRENTLY RESTRICTING RHO DETERMINATION TO |eta| < 2.4 
# process.kt6CaloJets.Rho_EtaMax     = cms.double( 2.4 ) 



# process.load("JetMETCorrections.Configuration.JetCorrectionServices_cff")

# from JetMETCorrections.Configuration.JetCorrectionProducersAllAlgos_cff import *
# from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *
# from JetMETCorrections.Configuration.DefaultJEC_cff import *

# # ak5CaloL1Fastjet.srcRho = 'kt6CaloJetsCentral:rho'
# # ak5CaloL1FastL2L3 = ak5CaloL1Fastjet.clone()

# # process.ak5CaloJetsL1FastL2L3 = cms.EDProducer(
# #     'CaloJetCorrectionProducer',
# #     src         = cms.InputTag('ak5CaloJets'),
# #     correctors  = cms.vstring('ak5CaloL1FastL2L3')
# #     )


# process.ak4PFJetsL1FastL2L3 = cms.EDProducer(
#     'PFJetCorrectionProducer',
#     src         = cms.InputTag('ak4PFJets'),
#     correctors  = cms.vstring('ak4PFL1FastL2L3')
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




process.JetProducer = cms.Path(
#  process.l1extraParticles

  # GCT
  process.RawToDigi
  *process.valRctDigis
  *process.valGctDigis
  *process.L1Extra

  # UCT
  *process.L1TCaloStage1_PPFromRaw
  
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

process.JetProducerSchedule = cms.Schedule( process.JetProducer )

# ================================================================================


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring( 
      'drop *',

      # GCT
      'keep *_l1extraParticles_*_*',
      # UCT
      'keep *_caloStage1LegacyFormatDigis_*_*',

      'keep *_genMetCalo_*_*',
      'keep *_genMetCaloAndNonPrompt_*_*',
      'keep *_genMetTrue_*_*',

      # HLT jet triggers
      'keep recoCaloJets_hltCaloJetL1FastJetCorrected_*_HLT2',
      'keep recoPFJets_hltAntiKT4PFJets*_*_HLT2',
      'keep *_hltSingleJet20_*_*',
      'keep *_hlt1PFJet20_*_*',
      'keep *_hlt1PFJetNoPU20_*_*',

      # Offline jets
      'keep *_ak5*Jets*_*_*',
      'drop *_ak5*Jets*_rho_*',
      'drop *_ak5*Jets*_sigma*_*',
      'keep *_kt6CaloJetsCentral_rho_*',
#      'keep *_*Jets*_*_*',

      'keep recoGenJets_*GenJets_*_HLT2',
      'keep *_ak4PFJets*_*_*',
      'keep *_ak4CaloJets*_*_*',

      # 'drop *_hlt*_rho_*',
      # 'drop *_*Gen*_rho_*',
      'keep PileupSummaryInfos_*_*_*',

      'drop *_*_*Digi*_',
      'keep *_hltL1GtObjectMap_*_HLT2',
      'keep FEDRawDataCollection_rawDataCollector_*_HLT2',
      'keep FEDRawDataCollection_source_*_HLT2',
      'keep edmTriggerResults_*_*_HLT2',
      'keep triggerTriggerEvent_*_*_HLT2' ),
    fileName = cms.untracked.string('file:hltReRunResults.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        #        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS162_V2::All', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend([process.digitisation_step,process.L1simulation_step])
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
