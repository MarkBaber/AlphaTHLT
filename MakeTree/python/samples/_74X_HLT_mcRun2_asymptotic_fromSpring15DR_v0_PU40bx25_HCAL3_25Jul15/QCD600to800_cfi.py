import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# QCD600to800
# ---------------------------------------------------------------------------------------------------------------------------------------
QCD600to800 = cms.PSet(
	name  = cms.string("QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8.root"),
	files = cms.untracked.vstring()
)
QCD600to800.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_29.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_3.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_1.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_21.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_8.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_5.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_4.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_7.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_6.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_18.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_22.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_9.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_13.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_20.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_23.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_14.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_19.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_15.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_2.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_11.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_10.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_24.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_27.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_25.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_17.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_16.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_26.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/74X_HLT_mcRun2_asymptotic_fromSpring15DR_v0_PU40bx25_HCAL3_25Jul15/QCD_Pt_600to800_TuneCUETP8M1_13TeV_pythia8/crab_QCD600to800/150725_221939/0000/hltReRunResults_28.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------