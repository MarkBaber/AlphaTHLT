import FWCore.ParameterSet.Config as cms

# General purpose tree maker
MakeTrees = cms.EDAnalyzer("MakeTrees",

    HLTResults = cms.untracked.InputTag("TriggerResults::HLT2"),


    #UCT inputs  
    #srcUctMht = cms.InputTag("caloStage1FinalDigis","", "L1TEMULATION"),
    srcUctMht = cms.InputTag("caloStage1FinalDigis"),
    srcUctMet = cms.InputTag("caloStage1FinalDigis"),
    srcUctJet = cms.InputTag("caloStage1FinalDigis"),
#    srcUctJet = cms.InputTag("caloStage1LegacyFormatDigis"),


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
    srcGctMht = cms.InputTag("hltL1extraParticles", "MHT"),
    srcGctMet = cms.InputTag("hltL1extraParticles", "MET"),

    srcGctJetCentral = cms.VInputTag(cms.InputTag("hltL1extraParticles", "Central")),
    srcGctJetForward = cms.VInputTag(cms.InputTag("hltL1extraParticles", "Forward")),

    srcGctJetAll = cms.VInputTag(cms.InputTag("hltL1extraParticles", "Central"),cms.InputTag("hltL1extraParticles", "Forward")),

    #Gen inputs
    srcGen4Jet = cms.VInputTag(cms.InputTag("ak4NoNuGenJets","","")),
#    srcGen4Jet = cms.VInputTag(cms.InputTag("ak4NoMuNoNuGenJets","","")),
#    srcGen5Jet = cms.VInputTag(cms.InputTag("ak5NoMuNoNuGenJets","","")),

    srcGenMetCalo             = cms.InputTag( "genMetCalo" ),
    srcGenMetCaloAndNonPrompt = cms.InputTag( "genMetCaloAndNonPrompt" ),
    srcGenMetTrue             = cms.InputTag( "genMetTrue" ),

                           
    srcHLTMetCalo                = cms.InputTag( "hltMet" ), 
    srcHLTMetCleanCalo           = cms.InputTag( "hltMetClean" ), 
    srcHLTMetCleanUsingJetIDCalo = cms.InputTag( "hltMetCleanUsingJetID" ), 
    srcHLTMetPF                  = cms.InputTag( "hltPFMETProducer" ),

    srcHLTMhtCalo             = cms.InputTag( "hltHtMht" ), 
    srcHLTMhtPF               = cms.InputTag( "hltPFHT" ), 


    # Gen particles                           
    MakeGenParticles          = cms.bool( True ),
    srcGenParticles           = cms.InputTag( "prunedGenParticles"),
    genElectronMinPt          = cms.double( 20. ), 
    genElectronMaxEta         = cms.double( 2.5 ),                           
    genMuonMinPt              = cms.double( 10. ),                           
    genMuonMaxEta             = cms.double( 2.5 ),                           
    genPhotonMinPt            = cms.double( 25. ),                           
    genPhotonMaxEta           = cms.double( 2.5 ),                           



    # HLT jets
    srcHLTAk4PF            = cms.VInputTag(cms.InputTag("hltAK4PFJetsCorrected","","")),                           
    srcHLTAk4Calo          = cms.VInputTag(cms.InputTag("hltAK4CaloJetsCorrected","","")), 
    srcHLTAk4CaloID        = cms.VInputTag(cms.InputTag("hltAK4CaloJetsCorrectedIDPassed","","")), 
    srcHLTAk4CaloNoFastJet = cms.VInputTag(cms.InputTag("hltAK4CaloJetsCorrectedIDPassedNoFastJet","","")),                           
#    srcHLTAk4PF     = cms.VInputTag(cms.InputTag("hltAntiKT4PFJets","","")),                           
#    srcHLTAk4PFNoPU = cms.VInputTag(cms.InputTag("hltAntiKT4PFJetsNoPU","","")),                           
#    srcHLTAk4Calo   = cms.VInputTag(cms.InputTag("hltCaloJetL1FastJetCorrected","","")),                           


    #Reco inputs
#     srcCaloJet = cms.VInputTag(cms.InputTag("ak5CaloJets")),
#     srcPfJet = cms.VInputTag(cms.InputTag("ak5PFJets")),
#     srcCaloJet = cms.VInputTag(cms.InputTag("PUsubAK5CaloJetProducer")),
#     srcPfJet   = cms.VInputTag(cms.InputTag("PUsubAK5PFJetProducer")),
                           

    # srcAk4Calo = cms.VInputTag(cms.InputTag("ak4CaloJetsL1FastL2L3")),
    # srcAk4PF   = cms.VInputTag(cms.InputTag("ak4PFJetsCHSL1FastL2L3")),
    srcAk4Calo = cms.VInputTag(),
    srcAk4PF   = cms.VInputTag(),
    

    #Parameters for HT
    #htThreshold = cms.double(30.0),
    maxjet      = cms.uint32(20),
    usePU       = cms.bool(True),

    # Cuts for jet skims
    jetMinPt    = cms.double(30.0),

    cenJetMinEta = cms.double(0.0),
    cenJetMaxEta = cms.double(3.0),
    forJetMinEta = cms.double(3.0),
    forJetMaxEta = cms.double(5.0),


    # Dynamic HT, AlphaT low-jet pt threshold
    dynamicJetThreshold = cms.double(30.0),


    )

