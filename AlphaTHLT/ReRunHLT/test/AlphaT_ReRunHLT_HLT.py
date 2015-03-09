# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: AlphaT_ReRunHLT --conditions PRE_LS172_V16::All -n 10 --eventcontent FEVTDEBUGHLT -s HLT --datatier GEN-SIM-DIGI --customise SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --geometry Extended2015 --magField 38T_PostLS1 --no_exec --filein file:step1.root --fileout file:step2.root
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT2')


# ================================================================================
# Select sample scenario
# ================================================================================
bx = "25ns"
#bx = "50ns"

#era = "Fall13"
era = "Spring14_70X" # Currently the GT/customisation is not working: No "L1RCTNoisyChannelMaskRcd" record found in the EventSetup.

# ================================================================================


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
#process.load('AlphaTHLT.ReRunHLT.HLT_AlphaT_cff')
process.load('AlphaTHLT.ReRunHLT.hlt_stage1_New_cff')


process.load("FWCore.MessageService.MessageLogger_cfi")
#process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
#    input = cms.untracked.int32(-1)
)


from AlphaTHLT.ReRunHLT.samples.samples_Signal_MC_cfi import *
selectedSample = T2tt_500_250 #T2tt_300_200 #T2tt_500_250 #T2cc_250_210


# Input source
process.source = cms.Source("PoolSource",

                            # Run on private MC samples (Run 'voms-proxy-init -out ~/myproxy -voms cms' first)
                            # fileNames = selectedSample.files,

                            # Fall13 Samples
                            # ------------------
                            
                            # 50ns
                            # fileNames = cms.untracked.vstring('/store/mc/Fall13dr/QCD_Pt-30to50_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx50_POSTLS162_V2-v1/00000/004C06A8-7FB9-E311-91B6-003048FEAED4.root')
                            # 25ns
                            # fileNames = cms.untracked.vstring('/store/mc/Fall13dr/QCD_Pt-30to50_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx25_POSTLS162_V2-v1/00000/0000AAF3-F8A6-E311-B13C-0025905964BA.root')

                            #    fileNames = cms.untracked.vstring('/store/mc/Fall13dr/QCD_Pt-800to1000_Tune4C_13TeV_pythia8/GEN-SIM-RAW/castor_tsg_PU40bx25_POSTLS162_V2-v1/00000/0002521B-7C9D-E311-874A-003048FFCB8C.root')
                            #    fileNames = cms.untracked.vstring("root://xrootd.unl.edu//store/mc/Fall13dr/Neutrino_Pt-2to20_gun/GEN-SIM-RAW/tsg_PU40bx25_POSTLS162_V2-v1/00005/02B79593-F47F-E311-8FF6-003048FFD796.root")
                            

                            # Spring14 Samples
                            # ------------------
                            # 25ns20PU
 fileNames = cms.untracked.vstring("/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/6893CF6A-5788-E311-A6DD-0002C90A3698.root")


 #                            secondaryFileNames = cms.untracked.vstring(
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/70666666-DD3E-E411-B4E1-0025B31E3C28.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/0E53A1A7-E53E-E411-8DC1-002590200868.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/2CA2EA99-E23E-E411-8EF3-001E67396874.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/7882CC0D-ED3E-E411-AD6F-001E67396FA9.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/844A3443-E03E-E411-9500-002590200B0C.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/88D536EB-F13E-E411-88D2-002481E14E58.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/BAC9A47D-E93E-E411-BE76-002590A80DE0.root",
#         "/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/GEN-SIM-RAW/PU20bx25_POSTLS170_V5-v1/00000/EAFF646A-EA3E-E411-BEE2-002590A80DE0.root",
#         ),
# 
                            # fileNames = cms.untracked.vstring("/store/mc/Spring14dr/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/745A841F-0E3F-E411-BED1-002590200828.root")


# secondaryFileNames = cms.untracked.vstring(
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/887CC49F-278A-E311-9039-002590E2F9D4.root",

# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/FC43D02E-C188-E311-BBC5-002590D60036.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/FA4B7C63-5888-E311-B8F2-0002C90A370A.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/E63DBC91-C088-E311-991C-00215E2EB6E2.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/BC3FE595-C088-E311-996C-90B11C050AD4.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/AE211795-C088-E311-8658-6C3BE5B59210.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/A2FAC29F-278A-E311-B3F7-002590E2F9D4.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/8C77B7DE-2F88-E311-B820-D48564447752.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/6E26EB88-C088-E311-B02F-001F296544A8.root",
# "/store/mc/Fall13/DYJetsToLL_M-50_13TeV-madgraph-pythia8/GEN-SIM/POSTLS162_V1-v2/10000/6893CF6A-5788-E311-A6DD-0002C90A3698.root",
# ),

#  fileNames = cms.untracked.vstring("/store/mc/Spring14dr/DYJetsToLL_M-50_13TeV-madgraph-pythia8/AODSIM/PU20bx25_POSTLS170_V5-v1/00000/00350711-84CB-E311-BDBF-0025901D4C44.root")


)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('AlphaT_ReRunHLT nevts:10'),
    name = cms.untracked.string('Applications')
)


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


       # 'keep *_hlt*_*_HLT2',
       # 'keep *_*_*_HLT2',
       # 'drop *_*_*Digi*_',
       # 'drop *_mix_*_*',

        # UCT
       'keep *BXVector*_*_*_*',
       'keep *_caloStage1FinalDigis_*_*',
       'keep l1extra*_*_*_*',
       # HLT
       'keep recoPFJets_hlt*_*_*',
       'keep recoCaloJets_hltAK4CaloJets*_*_*',
       # RECO
       'keep *_ak4PFJets*_*_*',
       'keep *_ak4CaloJets*_*_*',
       'keep *_fixedGridRho*_*_*',


       'keep triggerTriggerFilterObjectWithRefs_*_*_HLT2',
       'keep recoMETs_*_*_HLT2',
       'keep recoCaloMETs_*_*_HLT2',
       

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


# --------------------------------------------------------------------------------
# Get the correct GlobalTag

GT = ""
if era == "Fall13":
    if bx   == "25ns":
        GT = 'MCRUN2_72_V3A::All'
    elif bx == "50ns":
        GT = 'MCRUN2_72_V4A::All'
    else:
        print "Error: Bunch spacing '", bx, "' not recognised.\n"
        exit(0)
    pass
elif era == "Spring14_70X":
    if bx   == "25ns":
        GT = 'PHYS14_25_V1'
    elif bx == "50ns":
        GT = 'PHYS14_50_V1'
    else:
        print "Error: Bunch spacing '", bx, "' not recognised.\n"
        exit(0)
    pass
else:
    print "Error: Era '", era, "' not recognised.\n"
    exit(0)



print "Era = '", era, "', BX = '", bx, "', using GT:", GT
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

if era == "Spring14_70X":
    process.hltCsc2DRecHits.wireDigiTag  = cms.InputTag("simMuonCSCDigis","MuonCSCWireDigi")
    process.hltCsc2DRecHits.stripDigiTag = cms.InputTag("simMuonCSCDigis","MuonCSCStripDigi")


# # FIX incorrect HCAL tower configurations
# process.GlobalTag.toGet.append(
#     cms.PSet(
#         record  = cms.string( 'HcalRecoParamsRcd' ),
#         label   = cms.untracked.string( '' ),
#         connect = cms.untracked.string( 'frontier://FrontierProd/CMS_COND_44X_HCAL' ),
#         tag     = cms.string( 'HcalRecoParams_v8.0_mc' )
#     ) 
# )

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


