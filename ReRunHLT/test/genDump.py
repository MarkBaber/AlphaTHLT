# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: step2 --conditions POSTLS162_V2::All -n 10 --eventcontent FEVTDEBUGHLT -s HLT --datatier GEN-SIM-RAW --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry Extended2015 --magField 38T_PostLS1 --no_exec --filein file:step1.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('GenDump')

################################################################
process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
#process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(1000)
#  input = cms.untracked.int32(-1)
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

# from AlphaTHLT.ReRunHLT.samples.samples_Signal_MC_cfi import *
# selectedSample = T2tt_500_250 #T2cc_250_210

# Input source
process.source = cms.Source("PoolSource",


                            fileNames = cms.untracked.vstring('/store/mc/Fall13dr/QCD_Pt-30to50_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx25_POSTLS162_V2-v1/00000/0000AAF3-F8A6-E311-B13C-0025905964BA.root')

#                            fileNames = cms.untracked.vstring('/store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00001/4AAC6F5B-9979-E311-9DED-0025905A48F0.root')


#                            fileNames = selectedSample.files,


    #  secondaryFileNames = cms.untracked.vstring(
    #     'root://xrootd.unl.edu//store/mc/Fall13dr/TT_Tune4C_13TeV-pythia8-tauola/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00000/DEA553B3-B575-E311-9A10-00304867BFF2.root', 
    #     ),
    # fileNames = cms.untracked.vstring(
    #     'root://xrootd.unl.edu//store/mc/Fall13dr/TT_Tune4C_13TeV-pythia8-tauola/AODSIM/tsg_PU40bx25_POSTLS162_V2-v1/00000/A27B8942-5676-E311-A7C4-00261894398D.root', 
    #     ),

#     # Test locally
#     # --------------------
# #    fileNames = cms.untracked.vstring('file:Neutrino_Pt-2to20_gun_tsg_PU40bx25.root')






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

# process.GenParticleSkimmer = cms.Path(
#     process.prunedGenParticles
# #    *process.GenParticleFilter
# )


#process.GenParticleSkimSchedule = cms.Schedule( process.GenParticleSkimmer )


# ================================================================================


# Output definition
process.output = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = cms.untracked.vstring( 
       'drop *',


# gen particle skimmer
#        'keep *_prunedGenParticles_*_*',
#        'keep *_*_GenParticleFilter_*', 


#       'keep *_*_*_SIM',

       'keep GenEventInfoProduct_*_*_*',


      'keep *_genMetCalo_*_*',
      'keep *_genMetCaloAndNonPrompt_*_*',
      'keep *_genMetTrue_*_*',

      'keep recoGenJets_*GenJets_*_*',
      'keep *_ak4PFJets*_*_*',
      'keep *_ak4CaloJets*_*_*',
      'keep PileupSummaryInfos_*_*_*',


      ),
    fileName = cms.untracked.string('file:genDump.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        #        dataTier = cms.untracked.string('GEN-SIM-DIGI-RAW')
        dataTier = cms.untracked.string('GEN-SIM-RAW')
    )
)

# Additional output definition

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS1', '')

# Path and EndPath definitions
process.digitisation_step = cms.Path(process.pdigi)
#process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.output_step = cms.EndPath(process.output)

# Schedule definition
process.schedule = cms.Schedule()
#process.schedule.extend(process.GenParticleSkimSchedule)
process.schedule.extend([process.digitisation_step]) #,process.L1simulation_step])
process.schedule.extend(process.JetProducerSchedule)
#process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.output_step])

# ----------------------------------------------------------------------------------------------------
# customisation of the process.


# # customize the L1 emulator to run customiseL1EmulatorFromRaw with HLT to switchToSimStage1Digis 
# process.load( 'Configuration.StandardSequences.RawToDigi_cff' )
# process.load( 'Configuration.StandardSequences.SimL1Emulator_cff' )

# #process.load('L1Trigger.L1TCalorimeter.L1TCaloStage1_PPFromRaw_cff')
# process.load('L1Trigger/L1TCalorimeter/caloStage1RegionSF_cfi')
# import L1Trigger.L1TCalorimeter.L1TCaloStage1_customForHLT
# import L1Trigger.Configuration.L1Trigger_custom
# process = L1Trigger.L1TCalorimeter.L1TCaloStage1_customForHLT.customiseL1EmulatorFromRaw( process )
# process = L1Trigger.Configuration.L1Trigger_custom.customiseResetPrescalesAndMasks( process )
# ## ccla Add in additional stuff from postLS1Customs 
# from SLHCUpgradeSimulations.Configuration.postLS1Customs import *
# process = customise_HLT( process )

# # customize the HLT to use the emulated results 
# import HLTrigger.Configuration.customizeHLTforL1Emulator
# process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
# process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimStage1Digis( process )

# # --------------------




# # Automatic addition of the customisation function from HLTrigger.Configuration.customizeHLTforMC
# from HLTrigger.Configuration.customizeHLTforMC import customizeHLTforMC 

# #call to customisation function customizeHLTforMC imported from HLTrigger.Configuration.customizeHLTforMC
# process = customizeHLTforMC(process)

# # Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
# from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1 

# #call to customisation function customisePostLS1 imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
# process = customisePostLS1(process)

# # End of customisation functions
