import FWCore.ParameterSet.Config as cms

# General purpose tree maker
MakeTrees = cms.EDAnalyzer("MakeTrees",

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
    
    # MET
    genMetCaloSrc                = cms.InputTag( "genMetCalo" ),
    #genMetCaloAndNonPromptSrc   = cms.InputTag( "genMetCaloAndNonPrompt" ),
    genMetTrueSrc                = cms.InputTag( "genMetTrue" ),

    hltCaloMetSrc                = cms.InputTag( "hltMet"),
    hltCaloMetCleanSrc           = cms.InputTag( "hltMetClean"),
    hltCaloMetCleanUsingJetIDSrc = cms.InputTag( "hltMetCleanUsingJetID"),
    hltPFMetSrc                  = cms.InputTag( "hltPFMETProducer"),
    # srcHLTMhtCalo             = cms.InputTag( "hltHtMht" ), 
    # srcHLTMhtPF               = cms.InputTag( "hltPFHT" ), 


    # Trigger                           
    HLTResults     = cms.InputTag("TriggerResults"),          # Remulated
    HLTResultsData = cms.InputTag("TriggerResults","","HLT"), # Recorded (In ntuple)

)

