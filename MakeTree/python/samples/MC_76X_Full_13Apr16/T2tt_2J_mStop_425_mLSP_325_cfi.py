import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T2tt_2J_mStop_425_mLSP_325
# ---------------------------------------------------------------------------------------------------------------------------------------
T2tt_2J_mStop_425_mLSP_325 = cms.PSet(
	name  = cms.string("SMS-T2tt_mStop-425_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8.root"),
	files = cms.untracked.vstring()
)
T2tt_2J_mStop_425_mLSP_325.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/MC_76X_Full_13Apr16/SMS-T2tt_mStop-425_mLSP-325_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/crab_T2tt_mStop_425_mLSP_325/160413_185142/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------