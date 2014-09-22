import FWCore.ParameterSet.Config as cms

import datetime
now = datetime.datetime.now()
date = now.strftime("%Y-%m-%d")

#
# TTbar 
#
TTbar = cms.PSet(
    name  = cms.string("TT_Tune4C_13TeV-pythia8-tauola_" + date + ".root"),
    files = cms.untracked.vstring()
)
TTbar.files.extend([
        '/pnfs/hep.ph.ic.ac.uk/data/cms/store/user/mbaber/AlphaTHLT_POSTLS162_V2_PU40bx25_18Sep14/TT_Tune4C_13TeV-pythia8-tauola/hltReRunResults_55_1_ivn.root',
        
        ])





# #
# # ttbar 
# #
# ttbar = cms.PSet(
#     name  = cms.string("TT_Tune4C_13TeV-pythia8-tauola_" + date + ".root"),
#     files = cms.untracked.vstring()
# )
# ttbar.files.extend([
# ])
# ttbar = cms.PSet(
#     name  = cms.string("TT_Tune4C_13TeV-pythia8-tauola_" + date + ".root"),
#     files = cms.untracked.vstring()
#     files.extend([
#             ])
# )
