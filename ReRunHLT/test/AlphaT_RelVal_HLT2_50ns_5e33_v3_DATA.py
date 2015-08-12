# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: RelVal --step=HLT:50ns_5e33_v3 --conditions=auto:run2_hlt_50ns_5e33_v3 --filein=file:RelVal_HLT_50ns_5e33_v3_DATA.root --custom_conditions= --fileout=RelVal_HLT2_50ns_5e33_v3_DATA.root --number=100 --data --no_exec --datatier SIM-DIGI-RAW-HLTDEBUG --eventcontent=RAW --customise=HLTrigger/Configuration/CustomConfigs.L1THLT --customise=SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1_50ns --customise= --scenario=pp --python_filename=RelVal_HLT2_50ns_5e33_v3_DATA.py --processName=HLT2
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('HLTrigger.Configuration.HLT_50ns_5e33_v3_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

# Input source
process.source = cms.Source("PoolSource",
#    fileNames = cms.untracked.vstring('file:RelVal_HLT_50ns_5e33_v3_DATA.root'),
    fileNames = cms.untracked.vstring('/store/data/Run2015B/SingleMu/RAW/v1/000/251/025/00000/6A3297AE-8524-E511-B952-02163E0136DD.root'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RAWoutput = cms.OutputModule("PoolOutputModule",
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string('SIM-DIGI-RAW-HLTDEBUG'),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    #fileName = cms.untracked.string('RelVal_HLT2_50ns_5e33_v3_DATA.root'),
    fileName = cms.untracked.string('Run2015B_SingleMu.root'),
    outputCommands = process.RAWEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
from HLTrigger.Configuration.CustomConfigs import ProcessName
process = ProcessName(process)

from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_hlt_50ns_5e33_v3', '')

# Path and EndPath definitions
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWoutput_step = cms.EndPath(process.RAWoutput)

# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
process.schedule.extend([process.endjob_step,process.RAWoutput_step])

# customisation of the process.

# Automatic addition of the customisation function from HLTrigger.Configuration.CustomConfigs
from HLTrigger.Configuration.CustomConfigs import L1THLT 

#call to customisation function L1THLT imported from HLTrigger.Configuration.CustomConfigs
process = L1THLT(process)

# Automatic addition of the customisation function from SLHCUpgradeSimulations.Configuration.postLS1Customs
from SLHCUpgradeSimulations.Configuration.postLS1Customs import customisePostLS1_50ns 

#call to customisation function customisePostLS1_50ns imported from SLHCUpgradeSimulations.Configuration.postLS1Customs
process = customisePostLS1_50ns(process)

# End of customisation functions

