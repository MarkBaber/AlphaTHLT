import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T1tttt_2J_mGluino_1500_mLSP_100
# ---------------------------------------------------------------------------------------------------------------------------------------
T1tttt_2J_mGluino_1500_mLSP_100 = cms.PSet(
	name  = cms.string("SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
	files = cms.untracked.vstring()
)
T1tttt_2J_mGluino_1500_mLSP_100.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/MC_76X_Full_13Apr16/SMS-T1tttt_mGluino-1500_mLSP-100_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T1tttt_mGluino_1500_mLSP_100/160413_185304/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------