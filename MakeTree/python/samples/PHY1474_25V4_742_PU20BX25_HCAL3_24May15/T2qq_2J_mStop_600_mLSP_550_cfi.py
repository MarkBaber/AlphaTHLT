import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T2qq_2J_mStop_600_mLSP_550
# ---------------------------------------------------------------------------------------------------------------------------------------
T2qq_2J_mStop_600_mLSP_550 = cms.PSet(
	name  = cms.string("SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola.root"),
	files = cms.untracked.vstring()
)
T2qq_2J_mStop_600_mLSP_550.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_11.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_7.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_13.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_10.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_6.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_9.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_8.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_4.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_3.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_2.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_14.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_5.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2qq_2J_mStop-600_mLSP-550_Tune4C_13TeV-madgraph-tauola/crab_T2qq_2J_mStop_600_mLSP_550/150524_185327/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------