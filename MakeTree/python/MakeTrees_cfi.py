import FWCore.ParameterSet.Config as cms

# General purpose tree maker
MakeTrees = cms.EDAnalyzer("MakeTrees",


    #UCT inputs  
    #srcUctMht = cms.InputTag("caloStage1FinalDigis","", "L1TEMULATION"),
    srcUctMht = cms.InputTag("caloStage1FinalDigis"),
    srcUctMet = cms.InputTag("caloStage1FinalDigis"),
#    srcUctJet = cms.InputTag("caloStage1FinalDigis"),
    srcUctJet = cms.InputTag("caloStage1LegacyFormatDigis"),


    srcS2Met = cms.InputTag("Stage2JetProducer", "l1Stage2Met"),
    #Stage 2 inputs
    srcS2DonutMht = cms.InputTag("Stage2JetProducer","l1Stage2DonutPUSMht"),
    srcS2DonutJetCentral = cms.VInputTag(cms.InputTag("Stage2JetProducer","l1Stage2JetsDonutPUS")),
    #Stage 2 inputs
    srcS2NopusMht = cms.InputTag("Stage2JetProducer","l1Stage2NoPUSMht"),
    srcS2NopusJetCentral = cms.VInputTag(cms.InputTag("Stage2JetProducer","l1Stage2JetsNoPUS")),

    #Stage 2 inputs
    srcS2GlobalMht = cms.InputTag("Stage2JetProducer","l1Stage2GlobalPUSMht"),
    srcS2GlobalJetCentral = cms.VInputTag(cms.InputTag("Stage2JetProducer","l1Stage2JetsGlobalPUS")),

    #GCT inputs
    # srcGctMht = cms.InputTag("l1extraParticles", "MHT"),
    # srcGctMet = cms.InputTag("l1extraParticles", "MET"),
    srcGctJetCentral = cms.VInputTag(cms.InputTag("l1extraParticles", "Central")),
    srcGctJetAll = cms.VInputTag(cms.InputTag("l1extraParticles", "Central"),cms.InputTag("l1extraParticles", "Forward")),

    #Gen inputs
    srcGen4Jet = cms.VInputTag(cms.InputTag("ak4NoMuNoNuGenJets","","")),
    srcGen5Jet = cms.VInputTag(cms.InputTag("ak5NoMuNoNuGenJets","","")),

    # HLT jets
    srcHLTAk4PF     = cms.VInputTag(cms.InputTag("hltAntiKT4PFJets","","")),                           
    srcHLTAk4PFNoPU = cms.VInputTag(cms.InputTag("hltAntiKT4PFJetsNoPU","","")),                           
    srcHLTAk4Calo   = cms.VInputTag(cms.InputTag("hltCaloJetL1FastJetCorrected","","")),                           


    #Reco inputs
#     srcCaloJet = cms.VInputTag(cms.InputTag("ak5CaloJets")),
#     srcPfJet = cms.VInputTag(cms.InputTag("ak5PFJets")),
#     srcCaloJet = cms.VInputTag(cms.InputTag("PUsubAK5CaloJetProducer")),
#     srcPfJet   = cms.VInputTag(cms.InputTag("PUsubAK5PFJetProducer")),

    #srcCaloJet = cms.VInputTag(cms.InputTag("ak5CaloJetsL1FastL2L3")),
    #srcPfJet   = cms.VInputTag(cms.InputTag("ak5PFJetsL1FastL2L3")),
    srcCaloJet = cms.VInputTag(),
    srcPfJet   = cms.VInputTag(),
    

    #Parameters for HT
    #htThreshold = cms.double(30.0),
    maxjet      = cms.uint32(10),
    usePU       = cms.bool(True),

    # Cuts for jet skims
    jetMinPt    = cms.double(20.0),
    jetMaxEta   = cms.double(3.0),

    srcRegionEt = cms.InputTag("UCT2015Producer","regionEt"),
    srcRegionEta = cms.InputTag("UCT2015Producer","regionEta"),
    srcRegionPhi = cms.InputTag("UCT2015Producer","regionPhi"),
    )

