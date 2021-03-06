import FWCore.ParameterSet.Config as cms
# ---------------------------------------------------------------------------------------------------------------------------------------
# T2tt_2J_mStop_500_mLSP_325
# ---------------------------------------------------------------------------------------------------------------------------------------
T2tt_2J_mStop_500_mLSP_325 = cms.PSet(
	name  = cms.string("SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola.root"),
	files = cms.untracked.vstring()
)
T2tt_2J_mStop_500_mLSP_325.files.extend([
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_13.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_16.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_9.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_14.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_7.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_11.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_15.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_8.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_6.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_10.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_12.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_17.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_3.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_4.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_1.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_5.root',
	'root://gfe02.grid.hep.ph.ic.ac.uk/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/PHY1474_25V4_742_PU20BX25_HCAL3_24May15/SMS-T2tt_2J_mStop-500_mLSP-325_Tune4C_13TeV-madgraph-tauola/crab_T2tt_2J_mStop_500_mLSP_325/150524_185455/0000/hltReRunResults_2.root',
])

# ---------------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------------