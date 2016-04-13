import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# ZeroBias3
# ---------------------------------------------------------------------------------------------------------------------------------------
ZeroBias3 = cms.PSet(
	name  = cms.string("ZeroBias3.root"),
	files = cms.untracked.vstring()
)
ZeroBias3.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_85.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_167.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_51.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_18.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_29.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_35.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_28.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_41.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_27.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_42.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_10.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_19.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_9.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_14.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_15.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_6.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_11.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_7.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_8.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_13.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------