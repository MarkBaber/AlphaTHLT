import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# ZeroBias3
# ---------------------------------------------------------------------------------------------------------------------------------------
ZeroBias3 = cms.PSet(
	name  = cms.string("ZeroBias3.root"),
	files = cms.untracked.vstring()
)
ZeroBias3.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_76.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_145.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_16.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_68.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_121.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_94.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_128.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_75.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_123.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_38.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_71.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_108.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_122.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_59.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_53.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_139.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_26.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_91.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_47.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_105.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_27.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_77.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_135.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_141.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_60.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_120.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_109.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_146.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_49.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_130.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_137.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_112.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_56.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_95.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_82.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_127.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_70.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_99.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_80.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_40.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_114.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_107.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_30.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_34.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_37.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_39.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_78.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_89.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_24.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_92.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_33.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_55.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_28.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_41.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_62.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_61.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_125.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_69.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_65.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_44.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_140.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_126.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_86.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_157.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_172.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_149.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_22.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_148.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_131.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_144.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_156.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_83.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_136.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_58.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_170.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_160.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_169.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_132.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_25.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_17.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_167.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_142.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_87.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_153.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_147.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_81.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_117.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_159.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_97.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_119.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_50.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_106.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_152.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_181.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_164.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_48.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_154.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_103.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_102.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_52.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_31.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_133.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_138.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_155.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_151.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_143.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_174.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_45.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_46.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_184.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_163.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_74.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_93.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_177.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_179.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_43.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_19.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_178.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_90.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_64.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_63.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_110.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_98.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_84.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_116.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_104.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_111.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_113.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_66.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_96.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_171.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_54.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_101.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_175.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_79.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_124.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_36.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_20.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_57.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_73.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_166.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_85.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_51.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_18.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_29.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_35.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_42.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_Loose_12Apr16/ZeroBias3/crab_ZeroBias3/160412_124956/0000/hltReRunResults_10.root',
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