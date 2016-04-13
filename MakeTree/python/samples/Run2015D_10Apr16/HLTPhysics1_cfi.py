import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# HLTPhysics1
# ---------------------------------------------------------------------------------------------------------------------------------------
HLTPhysics1 = cms.PSet(
	name  = cms.string("HLTPhysics1.root"),
	files = cms.untracked.vstring()
)
HLTPhysics1.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_10Apr16/HLTPhysics1/crab_HLTPhysics1/160410_215126/0000/hltReRunResults_18.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/Run2015D_10Apr16/HLTPhysics1/crab_HLTPhysics1/160410_215126/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------