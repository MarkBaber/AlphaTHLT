import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T1bbbb_2J_mGluino_1000_mLSP_900
# ---------------------------------------------------------------------------------------------------------------------------------------
T1bbbb_2J_mGluino_1000_mLSP_900 = cms.PSet(
	name  = cms.string("SMS-T1bbbb_mGluino-1000_mLSP-900_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
	files = cms.untracked.vstring()
)
T1bbbb_2J_mGluino_1000_mLSP_900.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/MC_76X_Full_13Apr16/SMS-T1bbbb_mGluino-1000_mLSP-900_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T1bbbb_mGluino_1000_mLSP_900/160413_185635/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------