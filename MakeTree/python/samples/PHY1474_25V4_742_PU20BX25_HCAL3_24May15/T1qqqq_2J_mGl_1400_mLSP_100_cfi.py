import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T1qqqq_2J_mGl_1400_mLSP_100
# ---------------------------------------------------------------------------------------------------------------------------------------
T1qqqq_2J_mGl_1400_mLSP_100 = cms.PSet(
	name  = cms.string("SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola.root"),
	files = cms.untracked.vstring()
)
T1qqqq_2J_mGl_1400_mLSP_100.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_26.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_20.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_21.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_19.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_7.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_22.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_24.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_4.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_23.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_16.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_25.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_14.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_17.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_11.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/failed',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_15.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_8.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_9.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_5.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_6.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_3.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_10.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_2.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T1qqqq_2J_mGl-1400_mLSP-100_Tune4C_13TeV-madgraph-tauola/crab_T1qqqq_2J_mGl_1400_mLSP_100/150524_184916/0000/hltReRunResults_1.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------