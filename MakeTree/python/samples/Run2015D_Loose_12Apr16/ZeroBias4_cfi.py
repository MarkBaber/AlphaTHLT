import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# ZeroBias4
# ---------------------------------------------------------------------------------------------------------------------------------------
ZeroBias4 = cms.PSet(
	name  = cms.string("ZeroBias4.root"),
	files = cms.untracked.vstring()
)
ZeroBias4.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_37.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_72.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_34.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_35.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_45.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_20.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_27.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_23.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_22.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_24.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_19.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_18.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_17.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_13.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_11.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_16.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_9.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_7.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_15.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_8.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_6.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_14.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias4/crab_ZeroBias4/160412_125313/0000/hltReRunResults_10.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------