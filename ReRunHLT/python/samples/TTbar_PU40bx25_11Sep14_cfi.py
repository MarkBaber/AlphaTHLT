import FWCore.ParameterSet.Config as cms

#maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring() 
source = cms.Source ("PoolSource",fileNames = readFiles, secondaryFileNames = secFiles)
readFiles.extend( [
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_116_1_yVG.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_92_1_zX7.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_87_1_Ysn.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_74_1_3OZ.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_88_1_kyA.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_121_1_NDo.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_115_1_KVR.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_43_1_FwH.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_113_1_KET.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_120_1_jvy.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_119_1_nsx.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_90_1_kZP.root',
'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/TTbar_PU40bx25_11Sep14/hltReRunResults_41_1_tof.root',
    ])

secFiles.extend( [
               ] )
