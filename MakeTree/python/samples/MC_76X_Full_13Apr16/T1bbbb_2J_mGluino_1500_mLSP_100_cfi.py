import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T1bbbb_2J_mGluino_1500_mLSP_100
# ---------------------------------------------------------------------------------------------------------------------------------------
T1bbbb_2J_mGluino_1500_mLSP_100 = cms.PSet(
	name  = cms.string("SMS-T1bbbb_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
	files = cms.untracked.vstring()
)
T1bbbb_2J_mGluino_1500_mLSP_100.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/MC_76X_Full_13Apr16/SMS-T1bbbb_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T1bbbb_mGluino_1500_mLSP_100/160413_185517/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------