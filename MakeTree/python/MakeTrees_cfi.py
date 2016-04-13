import FWCore.ParameterSet.Config as cms

# General purpose tree maker
MakeTrees = cms.EDAnalyzer("MakeTrees",
                       
    # ----------------------------------------
    # Level -1
    # ----------------------------------------
    srcGctMht = cms.InputTag("hltL1extraParticles", "MHT"),
    srcGctMet = cms.InputTag("hltL1extraParticles", "MET"),
    srcGctJetCentral = cms.VInputTag(cms.InputTag("hltL1extraParticles", "Central")),
    srcGctJetForward = cms.VInputTag(cms.InputTag("hltL1extraParticles", "Forward")),
    srcGctJetAll = cms.VInputTag(cms.InputTag("hltL1extraParticles", "Central"),cms.InputTag("hltL1extraParticles", "Forward")),

    # MET
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

    # Jets
    hltCaloSrc          = cms.InputTag("hltAK4CaloJetsCorrected"),
    hltCaloIDSrc        = cms.InputTag("hltAK4CaloJetsCorrectedIDPassed"),
    hltPFSrc            = cms.InputTag("hltAK4PFJetsCorrected"),
    genSrc              = cms.InputTag("ak4GenJetsNoNu"),
    

    HLTResults     = cms.InputTag("TriggerResults"),          # Remulated
    HLTResultsData = cms.InputTag("TriggerResults","","HLT"), # Recorded (In ntuple)

    #hltCaloMetSrc = cms.InputTag( "hltMet","","HLT2"),
    hltCaloMetSrc = cms.InputTag( "hltMet"),
    hltPFMetSrc   = cms.InputTag( "hltPFMETProducer"),


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

