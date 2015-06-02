
// **************************************************
// Switches
// **************************************************

// UNCOMMENT TO RUN ON RECO
#define RECO
// Remove isolated leptons from gen and HLT jets
//#define LEPTON_XCLEANING
// **************************************************

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "TTree.h"
// L1
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"
#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"
// Jets
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/METReco/interface/CaloMET.h"
// NVTX
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
// PAT trigger
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

// Generator information
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "GeneratorInterface/LHEInterface/interface/LHEEvent.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"

#include "DataFormats/METReco/interface/GenMET.h" 


// Custom modules
#include "AlphaTHLT/Core/interface/AlphaT.h"
#include "AlphaTHLT/Core/interface/JetMatch.h"

TLorentzVector P4toTLV (reco::Particle::LorentzVector a){
  return TLorentzVector( a.px(), a.py(), a.pz(), a.energy() );
}




typedef std::vector<edm::InputTag> VInputTag;


class MakeTrees : public edm::EDAnalyzer {
  public:
    MakeTrees(const edm::ParameterSet& pset);
    virtual ~MakeTrees();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);


    void storeJet( TString jetCollName, const std::vector<const reco::Candidate*>& jetColl, 
                   std::map<TString,std::vector<Float_t>* >& pt,
		   std::map<TString,std::vector<Float_t>* >& px,
		   std::map<TString,std::vector<Float_t>* >& py,
		   std::map<TString,std::vector<Float_t>* >& eta,
		   std::map<TString,std::vector<Float_t>* >& phi);

    std::vector< const reco::Candidate*> skimJets(const std::vector<const reco::Candidate*>& inJets, 
						  double minPt, double minEta, double maxEta);



    float leadL1GenDeltaR( const std::vector<const reco::Candidate*>& gctCen, const std::vector<const reco::Candidate*>& gctFor, 
			   const std::vector<const reco::Candidate*>& genAk4);
    float leadHLTGenDeltaR( const std::vector<const reco::Candidate*>& hltAk4, const std::vector<const reco::Candidate*>& genAk4);


    void   lepJetDeltaR( std::vector<Float_t>* leptonEta, std::vector<Float_t>* leptonPhi, const std::vector<const reco::Candidate*>& jet,
                         std::vector<int>   &matchedJetIndexDeltaR, std::vector<float> &lepJetMinDeltaR );

  // void crosscleanIsolatedLeptons( std::vector<Float_t>* genLepPt, std::vector<Float_t>* genLepEta, std::vector<Float_t>* genLepPhi,
  // 				  std::vector<const reco::Candidate*>& genJet, std::vector<const reco::Candidate*>& hltJet, 
  // 				  int& nIsoLeptons, double maxDeltaR, double minPt);
  void crosscleanIsolatedLeptons( std::vector<Float_t>* genLepPt, std::vector<Float_t>* genLepEta, std::vector<Float_t>* genLepPhi,
				  std::vector<const reco::Candidate*>& genJet,
				  std::vector<const reco::Candidate*>& hltJet, 
				  std::vector<const reco::Candidate*>& hltCaloJet, 
				  int& nIsoLeptons, double maxDeltaR, double minPt,
				  std::vector<const reco::PFJet*>& hltPFJet,
				  std::vector<float>& genLeptonPt ,
				  std::vector<float>& genLeptonMatchedGenJetPt ,
				  std::vector<float>& genLeptonMatchedHLTPFJetPt ,
				  std::vector<float>& genLeptonMatchedHLTPFJetMuonEF ,
				  std::vector<float>& genLeptonMatchedHLTPFJetElectronEF );


  void removeLeptonsFromJets( std::vector <const reco::Candidate*>& inputJets, std::vector<int> cleaningIndex );
  
  


  private:
    //Get the UCT stuff
    edm::InputTag srcUctMET_;
    edm::InputTag srcUctMht_;
    edm::InputTag srcUctJet_;

    //Get the Stage 2 stuff
    edm::InputTag srcS2Met_;
    edm::InputTag srcS2DonutMht_;
    VInputTag srcS2DonutJetCentral_;

    edm::InputTag srcS2GlobalMht_;
    VInputTag srcS2GlobalJetCentral_;

    edm::InputTag srcS2NopusMht_;
    VInputTag srcS2NopusJetCentral_;

    //Get the GCT stuff
    edm::InputTag srcGctMht_;
    edm::InputTag srcGctMet_;
    VInputTag srcGctJetCentral_;
    VInputTag srcGctJetForward_;
    VInputTag srcGctJetAll_;

    //Get the Gen stuff
    VInputTag srcGen4Jet_;
    //VInputTag srcGen5Jet_;

    edm::InputTag srcGenMetCalo_;
    edm::InputTag srcGenMetCaloAndNonPrompt_;
    edm::InputTag srcGenMetTrue_;

    edm::InputTag srcHLTMetCalo_;
    edm::InputTag srcHLTMetCleanCalo_;
    edm::InputTag srcHLTMetCleanUsingJetIDCalo_;
    edm::InputTag srcHLTMetPF_;
    edm::InputTag srcHLTMhtCalo_;
    edm::InputTag srcHLTMhtPF_;



    edm::InputTag srcGenParticles_;
    bool          makeGenParticles;
    double        genElectronMinPt;
    double        genElectronMaxEta;
    double        genMuonMinPt;
    double        genMuonMaxEta;
    double        genPhotonMinPt;
    double        genPhotonMaxEta;

    // L1 seeds
    bool L1HTT175;
    bool L1ETM70;
    bool L1HTT175OrETM70;
    float L1Jet_DPhi;


    VInputTag srcHLTAk4Calo;
    VInputTag srcHLTAk4CaloNoFastJet;
    VInputTag srcHLTAk4PF;
    VInputTag srcHLTAk4PFNoPU;

    VInputTag srcAk4Calo;
    VInputTag srcAk4PF;


    //Define the levels
    std::vector<TString> lvl_;
    TTree* tree;


    std::map<TString,Float_t> mhtPt_;
    std::map<TString,Float_t> alphaT_;
    std::map<TString,Float_t> mhtPhi_;
    std::map<TString,Float_t> ht_;
    std::map<TString,Float_t> dht_;
    std::map<TString,Float_t> mhtDivHt_;

    std::map<TString,Float_t> metPt_;
    std::map<TString,Float_t> metPhi_;
    std::map<TString,Float_t> et_;

    std::map<TString,std::vector<Float_t>* > jetPt;
    std::map<TString,std::vector<Float_t>* > jetPx;
    std::map<TString,std::vector<Float_t>* > jetPy;
    std::map<TString,std::vector<Float_t>* > jetMuons_;
    std::map<TString,std::vector<Float_t>* > jetPhi;
    std::map<TString,std::vector<Float_t>* > jetEta;
    std::map<TString,std::vector<Float_t>* > jetPtsAll_;
    std::map<TString,std::vector<Float_t>* > jetPhisAll_;
    std::map<TString,std::vector<Float_t>* > jetEtasAll_;
    std::map<TString,Float_t> jetDPhi12_;
    UInt_t Muon_;
    UInt_t NVTX;

    // Gen leptons
    bool genLeptonVeto;
    bool genElectronVeto;
    bool genMuonVeto;
    bool genPhotonVeto;
    bool genIsoLeptonVeto;
    bool genIsoElectronVeto;
    std::vector<Float_t>* genElectronPt;
    std::vector<Float_t>* genElectronEta;
    std::vector<Float_t>* genElectronPhi;
    std::vector<Float_t>* genMuonPt;
    std::vector<Float_t>* genMuonEta;
    std::vector<Float_t>* genMuonPhi;
    std::vector<Float_t>* genPhotonPt;
    std::vector<Float_t>* genPhotonEta;
    std::vector<Float_t>* genPhotonPhi;

    // HLT paths
    std::vector<TString>    hltPathNames;
    std::map<TString, bool> hltPathFired;


    std::pair<float,float> genAk4AlphaTHT40;
    std::pair<float,float> hltAk4PFAlphaTHT40;
    std::pair<float,float> hltAk4CaloAlphaTHT40;
  //    std::pair<float,float> hltAk4CaloNoFastJetAlphaTHT40;
    std::pair<float,float> recoAk4PFAlphaTHT40;
    std::pair<float,float> recoAk4CaloAlphaTHT40;

    std::pair<float,float> genAk4AlphaTHT50;
    std::pair<float,float> hltAk4PFAlphaTHT50;
    std::pair<float,float> hltAk4CaloAlphaTHT50;
    std::pair<float,float> hltAk4CaloNoFastJetAlphaTHT50;
    std::pair<float,float> recoAk4PFAlphaTHT50;
    std::pair<float,float> recoAk4CaloAlphaTHT50;

  float genAk4_AlphaTPrime40;
  float hltAk4Calo_AlphaTPrime40;
  float hltAk4PF_AlphaTPrime40;  
  float recoAk4Calo_AlphaTPrime40;
  float recoAk4PF_AlphaTPrime40; 

    // Dynamic HT, AlphaT
    float dynamicJetThreshold;
    std::vector<std::pair<float,float> > genAk4DynamicAlphaTHT40;
    std::vector<std::pair<float,float> > hltAk4PFDynamicAlphaTHT40;
    std::vector<std::pair<float,float> > recoAk4PFDynamicAlphaTHT40;


    std::pair<float,float> genAk4MHT40;
    std::pair<float,float> hltAk4PFMHT40;
    std::pair<float,float> hltAk4CaloMHT40;
    std::pair<float,float> recoAk4PFMHT40;
    std::pair<float,float> recoAk4CaloMHT40;

    std::pair<float,float> genAk4ForMHT40;
    std::pair<float,float> hltAk4PFForMHT40;
    std::pair<float,float> recoAk4PFForMHT40;

  // DeltaR between leading jets of collections
  float L1GenDeltaR;
  float HLTGenDeltaR;
  UInt_t hpuVeto;

  float hltMetCaloPFMht40_DeltaPhi;
  float genMetCaloMht40_DeltaPhi;

  UInt_t genAk4NJet40;
  UInt_t hltAk4PFNJet40;
  UInt_t hltAk4CaloNJet40;
  UInt_t recoAk4PFNJet40;
  UInt_t recoAk4CaloNJet40;
  
  UInt_t genAk4NJet50;
  UInt_t hltAk4PFNJet50;
  UInt_t hltAk4CaloNJet50;
  UInt_t recoAk4PFNJet50;
  UInt_t recoAk4CaloNJet50;

  Int_t genAk4NJetBin40;
  Int_t hltAk4PFNJetBin40;
  Int_t hltAk4CaloNJetBin40;
  Int_t recoAk4PFNJetBin40;
  Int_t recoAk4CaloNJetBin40;

  Int_t genAk4HTBin40;
  Int_t hltAk4PFHTBin40;
  Int_t hltAk4CaloHTBin40;
  Int_t recoAk4PFHTBin40;
  Int_t recoAk4CaloHTBin40;


    UInt_t maxjet_;
    bool usePU_; 
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;

    // Jet skim cuts
    double minPt;
    double minEtaCen;
    double maxEtaCen;
    double minEtaFor;
    double maxEtaFor;


    edm::InputTag HLTResultsTag;
  
    float PThat;

  float genAk4ForMaxPt;
  float recoAk4PFForMaxPt;
  float hltAk4PFForMaxPt;

  // Lead jet
  float genAk4LeadJetPt;
  float recoAk4PFLeadJetPt;
  float hltAk4PFLeadJetPt;
  float hltAk4CaloLeadJetPt;
  // Second jet
  float genAk4SecondJetPt;
  float recoAk4PFSecondJetPt;
  float hltAk4PFSecondJetPt;
  float hltAk4CaloSecondJetPt;
  // Dijet avg
  float genAk4DijetAvgPt;
  float recoAk4PFDijetAvgPt;
  float hltAk4PFDijetAvgPt;
  float hltAk4CaloDijetAvgPt;


  // Biased deltaPhi
  float genAk4BiasedDPhi;   
  float recoAk4PFBiasedDPhi; 
  float hltAk4PFBiasedDPhi;  
  // Tag jet index
  Int_t genAk4BiasedDPhiIndex;   
  Int_t recoAk4PFBiasedDPhiIndex; 
  Int_t hltAk4PFBiasedDPhiIndex;  
  // AlphaT - Vector
  float genAk4VecAlphaT40;
  float recoAk4PFVecAlphaT40;
  float hltAk4PFVecAlphaT40;
  // BetaT - Scalar
  float genAk4ScaBetaT40;
  float recoAk4PFScaBetaT40;
  float hltAk4PFScaBetaT40;
  // BetaT - Vector
  float genAk4VecBetaT40;
  float recoAk4PFVecBetaT40;
  float hltAk4PFVecBetaT40;


  // Matched indices
  std::vector<int> genMatchedAk4HLTPF;

  // Znunu emulation
  //std::vector<int> genMuonMatchedGenMuonIndex;
  std::vector<float> genMuonMatchedGenMuonPt;
  std::vector<float> genMuonMatchedGenJetPt;
  std::vector<float> genMuonMatchedHLTPFJetPt;
  std::vector<float> genMuonMatchedHLTPFJetMuonEF;
  std::vector<float> genMuonMatchedHLTPFJetElectronEF;

  //std::vector<int> genMuonMatchedGenMuonIndex;
  std::vector<float> genElectronMatchedGenElectronPt;
  std::vector<float> genElectronMatchedGenJetPt;
  std::vector<float> genElectronMatchedHLTPFJetPt;
  std::vector<float> genElectronMatchedHLTPFJetMuonEF;
  std::vector<float> genElectronMatchedHLTPFJetElectronEF;

  // Gen-isolated leptons
  int nIsoElectrons;
  int nIsoMuons;


};

MakeTrees::MakeTrees(const edm::ParameterSet& pset){

    // Initialize the ntuple builder
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("Ntuple", "Ntuple");



    lvl_.push_back("genAk4");
    lvl_.push_back("genAk4For");

    lvl_.push_back("hltAk4PF");
    lvl_.push_back("hltAk4PFFor");
    lvl_.push_back("hltAk4Calo");
    lvl_.push_back("hltAk4CaloFor");
    lvl_.push_back("recoAk4PF");
    lvl_.push_back("recoAk4PFFor");
    lvl_.push_back("recoAk4Calo");
    lvl_.push_back("recoAk4CaloFor");

    lvl_.push_back("gctCen");
    lvl_.push_back("gctFor");
    // lvl_.push_back("hltAk4PFNoPU");



    // ********************************************************************************
    // Set the tree branches
    // ********************************************************************************

    NVTX  = 0;
    PThat = 0;
    genElectronPt  = new std::vector<Float_t>();
    genElectronEta = new std::vector<Float_t>();
    genElectronPhi = new std::vector<Float_t>();
    genMuonPt      = new std::vector<Float_t>();
    genMuonEta     = new std::vector<Float_t>();
    genMuonPhi     = new std::vector<Float_t>();
    genPhotonPt    = new std::vector<Float_t>();
    genPhotonEta   = new std::vector<Float_t>();
    genPhotonPhi   = new std::vector<Float_t>();

    // Forward jet max Pt
    genAk4ForMaxPt    = 0;
    recoAk4PFForMaxPt = 0;
    hltAk4PFForMaxPt  = 0;

    // Lead jet 
    genAk4LeadJetPt    = 0;
    recoAk4PFLeadJetPt = 0;
    hltAk4PFLeadJetPt  = 0;
    hltAk4CaloLeadJetPt= 0;
    // Second jet 
    genAk4SecondJetPt    = 0;
    recoAk4PFSecondJetPt = 0;
    hltAk4PFSecondJetPt  = 0;
    hltAk4CaloSecondJetPt= 0;
    // Dijet avg 
    genAk4DijetAvgPt    = 0;
    recoAk4PFDijetAvgPt = 0;
    hltAk4PFDijetAvgPt  = 0;
    hltAk4CaloDijetAvgPt= 0;


    genMuonMatchedGenMuonPt          = std::vector<float>();
    genMuonMatchedGenJetPt           = std::vector<float>();
    genMuonMatchedHLTPFJetPt         = std::vector<float>();
    genMuonMatchedHLTPFJetMuonEF     = std::vector<float>();
    genMuonMatchedHLTPFJetElectronEF = std::vector<float>();

    genElectronMatchedGenElectronPt      = std::vector<float>();
    genElectronMatchedGenJetPt           = std::vector<float>();
    genElectronMatchedHLTPFJetPt         = std::vector<float>();
    genElectronMatchedHLTPFJetMuonEF     = std::vector<float>();
    genElectronMatchedHLTPFJetElectronEF = std::vector<float>();



    // Event
    tree->Branch("NVTX",  &NVTX,  "NVTX/i");
    tree->Branch("PThat", &PThat, "PThat/f");


    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){

	mhtPt_[*iLvl]        = 0.0;
	mhtPhi_[*iLvl]       = 0.0;
	ht_[*iLvl]           = 0.0;
	dht_[*iLvl]          = 0.0;
	alphaT_[*iLvl]       = 0.0;
	mhtDivHt_[*iLvl]     = 0.0;

	metPt_[*iLvl]        = 0.0;
	metPhi_[*iLvl]       = 0.0;
	et_[*iLvl]           = 0.0;
	jetPt[*iLvl]       = new std::vector<Float_t>();
	jetPx[*iLvl]       = new std::vector<Float_t>();
	jetPy[*iLvl]       = new std::vector<Float_t>();
	jetPhi[*iLvl]      = new std::vector<Float_t>();
	jetEta[*iLvl]      = new std::vector<Float_t>();

	tree->Branch(*iLvl + "_Pt",  "std::vector<float>", &jetPt[*iLvl]);
	tree->Branch(*iLvl + "_Px",  "std::vector<float>", &jetPx[*iLvl]);
	tree->Branch(*iLvl + "_Py",  "std::vector<float>", &jetPy[*iLvl]);

	tree->Branch(*iLvl + "_Phi", "std::vector<float>", &jetPhi[*iLvl]);
	tree->Branch(*iLvl + "_Eta", "std::vector<float>", &jetEta[*iLvl]);

    } // End ilvl loop

    // Maximum forward jet PT
    tree->Branch("genAk4For_MaxPt",      &genAk4ForMaxPt,               "genAk4For_MaxPt/f");
    tree->Branch("recoAk4PFFor_MaxPt",   &recoAk4PFForMaxPt,            "recoAk4PFFor_MaxPt/f");
    tree->Branch("hltAk4PFFor_MaxPt",    &hltAk4PFForMaxPt,             "hltAk4PFFor_MaxPt/f");


    // Lead jet 
    tree->Branch("genAk4Lead_Pt",      &genAk4LeadJetPt,               "genAk4Lead_Pt/f");
    tree->Branch("recoAk4PFLead_Pt",   &recoAk4PFLeadJetPt,            "recoAk4PFLead_Pt/f");
    tree->Branch("hltAk4PFLead_Pt",    &hltAk4PFLeadJetPt,             "hltAk4PFLead_Pt/f");
    tree->Branch("hltAk4CaloLead_Pt",  &hltAk4CaloLeadJetPt,           "hltAk4CaloLead_Pt/f");
    // Second jet 
    tree->Branch("genAk4Second_Pt",      &genAk4SecondJetPt,               "genAk4Second_Pt/f");
    tree->Branch("recoAk4PFSecond_Pt",   &recoAk4PFSecondJetPt,            "recoAk4PFSecond_Pt/f");
    tree->Branch("hltAk4PFSecond_Pt",    &hltAk4PFSecondJetPt,             "hltAk4PFSecond_Pt/f");
    tree->Branch("hltAk4CaloSecond_Pt",  &hltAk4CaloSecondJetPt,           "hltAk4CaloSecond_Pt/f");
    // Dijet avg 
    tree->Branch("genAk4DijetAvg_Pt",      &genAk4DijetAvgPt,               "genAk4DijetAvg_Pt/f");
    tree->Branch("recoAk4PFDijetAvg_Pt",   &recoAk4PFDijetAvgPt,            "recoAk4PFDijetAvg_Pt/f");
    tree->Branch("hltAk4PFDijetAvg_Pt",    &hltAk4PFDijetAvgPt,             "hltAk4PFDijetAvg_Pt/f");
    tree->Branch("hltAk4CaloDijetAvg_Pt",  &hltAk4CaloDijetAvgPt,           "hltAk4CaloDijetAvg_Pt/f");


    tree->Branch("hpuVeto",                  &hpuVeto,                      "hpuVeto/b");
    tree->Branch("leadL1GenDeltaR",          &L1GenDeltaR,                  "leadL1GenDeltaR/f");
    tree->Branch("leadHLTGenDeltaR",         &HLTGenDeltaR,                 "leadHLTGenDeltaR/f");


    // Biased deltaPhi
    tree->Branch("genAk4_BiasedDPhi",         &genAk4BiasedDPhi,         "genAk4_BiasedDPhi/f");
    tree->Branch("recoAk4PF_BiasedDPhi",      &recoAk4PFBiasedDPhi,      "recoAk4PF_BiasedDPhi/f");
    tree->Branch("hltAk4PF_BiasedDPhi",       &hltAk4PFBiasedDPhi,       "hltAk4PF_BiasedDPhi/f");
    tree->Branch("genAk4_BiasedDPhiIndex",    &genAk4BiasedDPhiIndex,    "genAk4_BiasedDPhiIndex/I");
    tree->Branch("recoAk4PF_BiasedDPhiIndex", &recoAk4PFBiasedDPhiIndex, "recoAk4PF_BiasedDPhiIndex/I");
    tree->Branch("hltAk4PF_BiasedDPhiIndex",  &hltAk4PFBiasedDPhiIndex,  "hltAk4PF_BiasedDPhiIndex/I");

    // AlphaT - Vector 
    tree->Branch("genAk4_VecAlphaT40",    &genAk4VecAlphaT40,    "genAk4_VecAlphaT40/f"); 
    tree->Branch("recoAk4PF_VecAlphaT40", &recoAk4PFVecAlphaT40, "recoAk4PF_VecAlphaT40/f"); 
    tree->Branch("hltAk4PF_VecAlphaT40",  &hltAk4PFVecAlphaT40,  "hltAk4PF_VecAlphaT40/f"); 
    // BetaT - Scalar 
    tree->Branch("genAk4_ScaBetaT40",     &genAk4ScaBetaT40,     "genAk4_ScaBetaT40/f"); 
    tree->Branch("recoAk4PF_ScaBetaT40",  &recoAk4PFScaBetaT40,  "recoAk4PF_ScaBetaT40/f"); 
    tree->Branch("hltAk4PF_ScaBetaT40",   &hltAk4PFScaBetaT40,   "hltAk4PF_ScaBetaT40/f"); 
    // BetaT - Vector 
    tree->Branch("genAk4_VecBetaT40",     &genAk4VecBetaT40,     "genAk4_VecBetaT40/f"); 
    tree->Branch("recoAk4PF_VecBetaT40",  &recoAk4PFVecBetaT40,  "recoAk4PF_VecBetaT40/f"); 
    tree->Branch("hltAk4PF_VecBetaT40",   &hltAk4PFVecBetaT40,   "hltAk4PF_VecBetaT40/f"); 



    // AlphaT, HT  
    tree->Branch("genAk4_AlphaT40",      &genAk4AlphaTHT40.first,       "genAk4_AlphaT40/f");
    tree->Branch("genAk4_HT40",          &genAk4AlphaTHT40.second,      "genAk4_HT40/f");
    tree->Branch("hltAk4PF_AlphaT40",    &hltAk4PFAlphaTHT40.first,     "hltAk4PF_AlphaT40/f");
    tree->Branch("hltAk4PF_HT40",        &hltAk4PFAlphaTHT40.second,    "hltAk4PF_HT40/f");
    tree->Branch("hltAk4Calo_AlphaT40",  &hltAk4CaloAlphaTHT40.first,   "hltAk4Calo_AlphaT40/f");
    tree->Branch("hltAk4Calo_HT40",      &hltAk4CaloAlphaTHT40.second,  "hltAk4Calo_HT40/f");
    tree->Branch("recoAk4PF_AlphaT40",   &recoAk4PFAlphaTHT40.first,    "recoAk4PF_AlphaT40/f");
    tree->Branch("recoAk4PF_HT40",       &recoAk4PFAlphaTHT40.second,   "recoAk4PF_HT40/f");
    tree->Branch("recoAk4Calo_AlphaT40", &recoAk4CaloAlphaTHT40.first,  "recoAk4Calo_AlphaT40/f");
    tree->Branch("recoAk4Calo_HT40",     &recoAk4CaloAlphaTHT40.second, "recoAk4Calo_HT40/f");

    tree->Branch("genAk4_AlphaT50",      &genAk4AlphaTHT50.first,       "genAk4_AlphaT50/f");
    tree->Branch("genAk4_HT50",          &genAk4AlphaTHT50.second,      "genAk4_HT50/f");
    tree->Branch("hltAk4PF_AlphaT50",    &hltAk4PFAlphaTHT50.first,     "hltAk4PF_AlphaT50/f");
    tree->Branch("hltAk4PF_HT50",        &hltAk4PFAlphaTHT50.second,    "hltAk4PF_HT50/f");
    tree->Branch("hltAk4Calo_AlphaT50",  &hltAk4CaloAlphaTHT50.first,   "hltAk4Calo_AlphaT50/f");
    tree->Branch("hltAk4Calo_HT50",      &hltAk4CaloAlphaTHT50.second,  "hltAk4Calo_HT50/f");
    tree->Branch("recoAk4PF_AlphaT50",   &recoAk4PFAlphaTHT50.first,    "recoAk4PF_AlphaT50/f");
    tree->Branch("recoAk4PF_HT50",       &recoAk4PFAlphaTHT50.second,   "recoAk4PF_HT50/f");
    tree->Branch("recoAk4Calo_AlphaT50", &recoAk4CaloAlphaTHT50.first,  "recoAk4Calo_AlphaT50/f");
    tree->Branch("recoAk4Calo_HT50",     &recoAk4CaloAlphaTHT50.second, "recoAk4Calo_HT50/f");


    tree->Branch("genAk4_AlphaTPrime40",      &genAk4_AlphaTPrime40,      "genAk4_AlphaTPrime40/f");
    tree->Branch("hltAk4Calo_AlphaTPrime40",  &hltAk4Calo_AlphaTPrime40,  "hltAk4Calo_AlphaTPrime40/f");
    tree->Branch("hltAk4PF_AlphaTPrime40",    &hltAk4PF_AlphaTPrime40,    "hltAk4PF_AlphaTPrime40/f");
    tree->Branch("recoAk4Calo_AlphaTPrime40", &recoAk4Calo_AlphaTPrime40, "recoAk4Calo_AlphaTPrime40/f");
    tree->Branch("recoAk4PF_AlphaTPrime40",   &recoAk4PF_AlphaTPrime40,   "recoAk4PF_AlphaTPrime40/f");


    // Dynamic AlphaT, HT
    tree->Branch("genAk4_DynamicAlphaTHT40",    "std::vector<std::pair<float,float>>", &genAk4DynamicAlphaTHT40);
    tree->Branch("hltAk4PF_DynamicAlphaTHT40",  "std::vector<std::pair<float,float>>", &hltAk4PFDynamicAlphaTHT40);
    tree->Branch("recoAk4PF_DynamicAlphaTHT40", "std::vector<std::pair<float,float>>", &recoAk4PFDynamicAlphaTHT40);


    // MHT
    tree->Branch("genAk4_MhtPT40",       &genAk4MHT40.first,       "genAk4_MhtPT40/f");
    tree->Branch("genAk4_MhtPhi40",      &genAk4MHT40.second,      "genAk4_MhtPhi40/f");
    tree->Branch("hltAk4PF_MhtPT40",     &hltAk4PFMHT40.first,     "hltAk4PF_MhtPT40/f");
    tree->Branch("hltAk4PF_MhtPhi40",    &hltAk4PFMHT40.second,    "hltAk4PF_MhtPhi40/f");
    tree->Branch("hltAk4Calo_MhtPT40",   &hltAk4CaloMHT40.first,   "hltAk4Calo_MhtPT40/f");
    tree->Branch("hltAk4Calo_MhtPhi40",  &hltAk4CaloMHT40.second,  "hltAk4Calo_MhtPhi40/f");
    tree->Branch("recoAk4PF_MhtPT40",    &recoAk4PFMHT40.first,    "recoAk4PF_MhtPT40/f");
    tree->Branch("recoAk4PF_MhtPhi40",   &recoAk4PFMHT40.second,   "recoAk4PF_MhtPhi40/f");
    tree->Branch("recoAk4Calo_MhtPT40",  &recoAk4CaloMHT40.first,  "recoAk4Calo_MhtPT40/f");
    tree->Branch("recoAk4Calo_MhtPhi40", &recoAk4CaloMHT40.second, "recoAk4Calo_MhtPhi40/f");

    tree->Branch("genAk4For_MhtPT40",       &genAk4ForMHT40.first,       "genAk4For_MhtPT40/f");
    tree->Branch("genAk4For_MhtPhi40",      &genAk4ForMHT40.second,      "genAk4For_MhtPhi40/f");
    tree->Branch("hltAk4PFFor_MhtPT40",     &hltAk4PFForMHT40.first,     "hltAk4PFFor_MhtPT40/f");
    tree->Branch("hltAk4PFFor_MhtPhi40",    &hltAk4PFForMHT40.second,    "hltAk4PFFor_MhtPhi40/f");
    tree->Branch("recoAk4PFFor_MhtPT40",    &recoAk4PFForMHT40.first,    "recoAk4PFFor_MhtPT40/f");
    tree->Branch("recoAk4PFFor_MhtPhi40",   &recoAk4PFForMHT40.second,   "recoAk4PFFor_MhtPhi40/f");


    tree->Branch("hltMetCaloPFMht40_DeltaPhi", &hltMetCaloPFMht40_DeltaPhi, "hltMetCaloPFMht40_DeltaPhi/f");
    tree->Branch("genMetCaloMht40_DeltaPhi",   &genMetCaloMht40_DeltaPhi,   "genMetCaloMht40_DeltaPhi/f");

    tree->Branch("genAk4_NJet40",      &genAk4NJet40,       "genAk4_NJet40/i");
    tree->Branch("hltAk4PF_NJet40",    &hltAk4PFNJet40,     "hltAk4PF_NJet40/i");
    tree->Branch("hltAk4Calo_NJet40",  &hltAk4CaloNJet40,   "hltAk4Calo_NJet40/i");
    tree->Branch("recoAk4PF_NJet40",   &recoAk4PFNJet40,    "recoAk4PF_NJet40/i");
    tree->Branch("recoAk4Calo_NJet40", &recoAk4CaloNJet40,  "recoAk4Calo_NJet40/i");

    tree->Branch("genAk4_NJet50",      &genAk4NJet50,       "genAk4_NJet50/i");
    tree->Branch("hltAk4PF_NJet50",    &hltAk4PFNJet50,     "hltAk4PF_NJet50/i");
    tree->Branch("hltAk4Calo_NJet50",  &hltAk4CaloNJet50,   "hltAk4Calo_NJet50/i");
    tree->Branch("recoAk4PF_NJet50",   &recoAk4PFNJet50,    "recoAk4PF_NJet50/i");
    tree->Branch("recoAk4Calo_NJet50", &recoAk4CaloNJet50,  "recoAk4Calo_NJet50/i");

    tree->Branch("genAk4_NJetBin40",      &genAk4NJetBin40,       "genAk4_NJetBin40/I");
    tree->Branch("hltAk4PF_NJetBin40",    &hltAk4PFNJetBin40,     "hltAk4PF_NJetBin40/I");
    tree->Branch("hltAk4Calo_NJetBin40",  &hltAk4CaloNJetBin40,   "hltAk4Calo_NJetBin40/I");
    tree->Branch("recoAk4PF_NJetBin40",   &recoAk4PFNJetBin40,    "recoAk4PF_NJetBin40/I");
    tree->Branch("recoAk4Calo_NJetBin40", &recoAk4CaloNJetBin40,  "recoAk4Calo_NJetBin40/I");

    tree->Branch("genAk4_HTBin40",      &genAk4HTBin40,       "genAk4_HTBin40/I");
    tree->Branch("hltAk4PF_HTBin40",    &hltAk4PFHTBin40,     "hltAk4PF_HTBin40/I");
    tree->Branch("hltAk4Calo_HTBin40",  &hltAk4CaloHTBin40,   "hltAk4Calo_HTBin40/I");
    tree->Branch("recoAk4PF_HTBin40",   &recoAk4PFHTBin40,    "recoAk4PF_HTBin40/I");
    tree->Branch("recoAk4Calo_HTBin40", &recoAk4CaloHTBin40,  "recoAk4Calo_HTBin40/I");

    // Energy sums
    tree->Branch("gct_Ht",       &ht_["gct"],      "gct_Ht/f");
    tree->Branch("gct_MhtPt",    &mhtPt_["gct"],   "gct_MhtPt/f");
    tree->Branch("gct_MhtPhi",   &mhtPhi_["gct"],  "gct_MhtPhi/f");
    tree->Branch("gct_MhtDivHt", &mhtDivHt_["gct"],"gct_MhtDivHt/f"); 

    tree->Branch("gct_Et",     &et_["gct"],     "gct_Et/f");
    tree->Branch("gct_MetPt",  &metPt_["gct"],  "gct_MetPt/f");
    tree->Branch("gct_MetPhi", &metPhi_["gct"], "gct_MetPhi/f");

    // GEN MET
    tree->Branch("genMetCalo_MetPt",              &metPt_["genMetCalo"],              "genMetCalo_MetPt/f");
    tree->Branch("genMetCaloAndNonPrompt_MetPt",  &metPt_["genMetCaloAndNonPrompt"],  "genMetCaloAndNonPrompt_MetPt/f");
    tree->Branch("genMetTrue_MetPt",              &metPt_["genMetTrue"],              "genMetTrue_MetPt/f");
    tree->Branch("genMetCalo_MetPhi",             &metPhi_["genMetCalo"],             "genMetCalo_MetPhi/f");
    tree->Branch("genMetCaloAndNonPrompt_MetPhi", &metPhi_["genMetCaloAndNonPrompt"], "genMetCaloAndNonPrompt_MetPhi/f");
    tree->Branch("genMetTrue_MetPhi",             &metPhi_["genMetTrue"],             "genMetTrue_MetPhi/f");

    // RECO MET
    tree->Branch("hltMetCalo_MetPT",                 &metPt_["hltMetCalo"],                 "hltMetCalo_MetPt/f");
    tree->Branch("hltMetCalo_MetPhi",                &metPhi_["hltMetCalo"],                "hltMetCalo_MetPhi/f");
    tree->Branch("hltMetCleanCalo_MetPT",            &metPt_["hltMetCleanCalo"],            "hltMetCleanCalo_MetPt/f");
    tree->Branch("hltMetCleanCalo_MetPhi",           &metPhi_["hltMetCleanCalo"],           "hltMetCleanCalo_MetPhi/f");
    tree->Branch("hltMetCleanUsingJetIDCalo_MetPT",  &metPt_["hltMetCleanUsingJetIDCalo"],  "hltMetCleanUsingJetIDCalo_MetPt/f");
    tree->Branch("hltMetCleanUsingJetIDCalo_MetPhi", &metPhi_["hltMetCleanUsingJetIDCalo"], "hltMetCleanUsingJetIDCalo_MetPhi/f");
    tree->Branch("hltMetPF_MetPT",                   &metPt_["hltMetPF"],                   "hltMetPF_MetPt/f");
    tree->Branch("hltMetPF_MetPhi",                  &metPhi_["hltMetPF"],                  "hltMetPF_MetPhi/f");
    // MHT
    tree->Branch("hltMhtCalo_MhtPT",              &mhtPt_["hltMhtCalo"],             "hltMhtCalo_MetPt/f");
    tree->Branch("hltMhtPF_MhtPT",                &mhtPt_["hltMhtPF"],               "hltMhtPF_MetPt/f");

    // Gen leptons
    tree->Branch("genLeptonVeto",    &genLeptonVeto,       "genLeptonVeto/b");
    tree->Branch("genElectronVeto",  &genElectronVeto,     "genElectronVeto/b");
    tree->Branch("genMuonVeto",      &genMuonVeto,         "genMuonVeto/b");
    tree->Branch("genPhotonVeto",    &genPhotonVeto,       "genPhotonVeto/b");
    tree->Branch("genElectron_Pt",   "std::vector<float>", &genElectronPt);
    tree->Branch("genElectron_Eta",  "std::vector<float>", &genElectronEta);
    tree->Branch("genElectron_Phi",  "std::vector<float>", &genElectronPhi);		
    tree->Branch("genMuon_Pt",       "std::vector<float>", &genMuonPt);
    tree->Branch("genMuon_Eta",      "std::vector<float>", &genMuonEta);
    tree->Branch("genMuon_Phi",      "std::vector<float>", &genMuonPhi);		
    tree->Branch("genPhoton_Pt",     "std::vector<float>", &genPhotonPt);
    tree->Branch("genPhoton_Eta",    "std::vector<float>", &genPhotonEta);
    tree->Branch("genPhoton_Phi",    "std::vector<float>", &genPhotonPhi);		


    tree->Branch("genMatchedAk4HLTPF", "std::vector<int>", &genMatchedAk4HLTPF);


    tree->Branch("nIsoMuons",     &nIsoMuons,     "nIsoMuons/i");
    tree->Branch("nIsoElectrons", &nIsoElectrons, "nIsoMuons/i");

    tree->Branch("genMuonMatchedGenMuon_Pt",          "std::vector<float>", &genMuonMatchedGenMuonPt);
    tree->Branch("genMuonMatchedGenJet_Pt",           "std::vector<float>", &genMuonMatchedGenJetPt);
    tree->Branch("genMuonMatchedHLTPFJet_Pt",         "std::vector<float>", &genMuonMatchedHLTPFJetPt);
    tree->Branch("genMuonMatchedHLTPFJet_MuonEF",     "std::vector<float>", &genMuonMatchedHLTPFJetMuonEF);
    tree->Branch("genMuonMatchedHLTPFJet_ElectronEF", "std::vector<float>", &genMuonMatchedHLTPFJetElectronEF);

    tree->Branch("genElectronMatchedGenElectron_Pt",          "std::vector<float>", &genElectronMatchedGenElectronPt);
    tree->Branch("genElectronMatchedGenJet_Pt",           "std::vector<float>",     &genElectronMatchedGenJetPt);
    tree->Branch("genElectronMatchedHLTPFJet_Pt",         "std::vector<float>",     &genElectronMatchedHLTPFJetPt);
    tree->Branch("genElectronMatchedHLTPFJet_MuonEF",     "std::vector<float>",     &genElectronMatchedHLTPFJetMuonEF);
    tree->Branch("genElectronMatchedHLTPFJet_ElectronEF", "std::vector<float>",     &genElectronMatchedHLTPFJetElectronEF);



    // Store L1 seeds
    // ------------------------------------------------------------

    tree->Branch("L1HTT175",        &L1HTT175,         "L1HTT175/b");
    tree->Branch("L1ETM70",         &L1ETM70,          "L1ETM70/b");
    tree->Branch("L1HTT175OrETM70", &L1HTT175OrETM70,  "L1HTT175OrETM70/b");
    tree->Branch("L1Jet_DPhi",      &L1Jet_DPhi,       "L1Jet_DPhi/f");


    // Store HLT paths
    // ------------------------------------------------------------


     hltPathNames.push_back("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57_v1"); 
     hltPathNames.push_back("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_v1"); 
     hltPathNames.push_back("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53_v1"); 
     hltPathNames.push_back("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52_v1"); 
     hltPathNames.push_back("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51_v1"); 

     hltPathNames.push_back("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v1"); 
     hltPathNames.push_back("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_v1");
     hltPathNames.push_back("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_v1"); 
     hltPathNames.push_back("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_v1"); 
     hltPathNames.push_back("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_v1");

                                    
     hltPathNames.push_back("HLT_PFHT200_v1");
     hltPathNames.push_back("HLT_PFHT250_v1");
     hltPathNames.push_back("HLT_PFHT300_v1");
     hltPathNames.push_back("HLT_PFHT350_v2");
     hltPathNames.push_back("HLT_PFHT400_v1");
     hltPathNames.push_back("HLT_PFHT475_v1"); 
                                    
     hltPathNames.push_back("HLT_PFHT800_v1");
     hltPathNames.push_back("HLT_PFHT350_PFMET100_NoiseCleaned_v1");
     hltPathNames.push_back("HLT_PFMET170_NoiseCleaned_v2");
     hltPathNames.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v2");
                                  


     hltPathNames.push_back("HLT_Rsq0p25_v1"); 
     hltPathNames.push_back("HLT_Rsq0p30_v1"); 
     hltPathNames.push_back("HLT_RsqMR240_Rsq0p09_MR200_v1");
     hltPathNames.push_back("HLT_RsqMR240_Rsq0p09_MR200_4jet_v1");
     hltPathNames.push_back("HLT_RsqMR270_Rsq0p09_MR200_v1");
     hltPathNames.push_back("HLT_RsqMR270_Rsq0p09_MR200_4jet_v1");
                                 

     hltPathNames.push_back("HLT_PFHT600_v2");
     hltPathNames.push_back("HLT_PFHT650_v2");








    
    // Trigger bits
    for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
      TString path = hltPathNames[ iPath ];
      hltPathFired[ path ] = false;
      tree->Branch( path, &hltPathFired[ path ], path + "/b" );
    }






    HLTResultsTag = pset.getUntrackedParameter("HLTResults", edm::InputTag("TriggerResults","HLT"));

    srcUctMET_ = pset.getParameter<edm::InputTag>("srcUctMet");
    srcUctMht_ = pset.getParameter<edm::InputTag>("srcUctMht");
    srcUctJet_ = pset.getParameter<edm::InputTag>("srcUctJet");

    srcGenMetCalo_                = pset.getParameter<edm::InputTag>("srcGenMetCalo");
    srcGenMetCaloAndNonPrompt_    = pset.getParameter<edm::InputTag>("srcGenMetCaloAndNonPrompt");
    srcGenMetTrue_                = pset.getParameter<edm::InputTag>("srcGenMetTrue");

    srcHLTMetCalo_                = pset.getParameter<edm::InputTag>("srcHLTMetCalo");
    srcHLTMetCleanCalo_           = pset.getParameter<edm::InputTag>("srcHLTMetCleanCalo");   
    srcHLTMetCleanUsingJetIDCalo_ = pset.getParameter<edm::InputTag>("srcHLTMetCleanUsingJetIDCalo");
    srcHLTMetPF_                  = pset.getParameter<edm::InputTag>("srcHLTMetPF");
    srcHLTMhtCalo_                = pset.getParameter<edm::InputTag>("srcHLTMhtCalo");
    srcHLTMhtPF_                  = pset.getParameter<edm::InputTag>("srcHLTMhtPF");


    // Gen particles
    makeGenParticles           = pset.getParameter<bool>("MakeGenParticles");
    srcGenParticles_           = pset.getParameter<edm::InputTag>("srcGenParticles");
    genElectronMinPt           = pset.getParameter<double>("genElectronMinPt");
    genElectronMaxEta          = pset.getParameter<double>("genElectronMaxEta");
    genMuonMinPt               = pset.getParameter<double>("genMuonMinPt");
    genMuonMaxEta              = pset.getParameter<double>("genMuonMaxEta");
    genPhotonMinPt             = pset.getParameter<double>("genPhotonMinPt");
    genPhotonMaxEta            = pset.getParameter<double>("genPhotonMaxEta");


    srcGctMht_        = pset.getParameter<edm::InputTag>("srcGctMht");
    srcGctMet_        = pset.getParameter<edm::InputTag>("srcGctMet");
    srcGctJetCentral_ = pset.getParameter<VInputTag>("srcGctJetCentral");
    srcGctJetForward_ = pset.getParameter<VInputTag>("srcGctJetForward");
    // srcGctJetAll_     = pset.getParameter<VInputTag>("srcGctJetAll");

    srcGen4Jet_        = pset.getParameter<VInputTag>("srcGen4Jet");
    //srcGen5Jet_        = pset.getParameter<VInputTag>("srcGen5Jet");

    srcHLTAk4Calo          = pset.getParameter<VInputTag>("srcHLTAk4Calo");
    srcHLTAk4CaloNoFastJet = pset.getParameter<VInputTag>("srcHLTAk4CaloNoFastJet");
    srcHLTAk4PF            = pset.getParameter<VInputTag>("srcHLTAk4PF");
    // srcHLTAk4PFNoPU    = pset.getParameter<VInputTag>("srcHLTAk4PFNoPU");

#ifdef RECO    
    srcAk4Calo       = pset.getParameter<VInputTag>("srcAk4Calo");
    srcAk4PF         = pset.getParameter<VInputTag>("srcAk4PF");
#endif
    //Parameters for the HSums
    //htThreshold_      = pset.getParameter<double>("htThreshold");
    maxjet_           = pset.getParameter<unsigned int>("maxjet");
    usePU_            = pset.getParameter<bool>("usePU");

    // Jet skim cuts
    minPt             = pset.getParameter<double>("jetMinPt");
    minEtaCen         = pset.getParameter<double>("cenJetMinEta");
    maxEtaCen         = pset.getParameter<double>("cenJetMaxEta");
    minEtaFor         = pset.getParameter<double>("forJetMinEta");
    maxEtaFor         = pset.getParameter<double>("forJetMaxEta");

    // Dynamic HT, AlphaT low-jet threshold
    dynamicJetThreshold = pset.getParameter<double>("dynamicJetThreshold");

}

MakeTrees::~MakeTrees() {
}


namespace {

    // Predicate to sort candidates by descending pt
    class CandPtSorter {
	public:
	    bool operator()(const reco::Candidate* candA, const reco::Candidate* candB)
		const {
		    return candA->pt() > candB->pt();
		}
    };

  // Turn a set of InputTags into a collection of candidate pointers.
  std::vector<const reco::Candidate*> getCollections(const edm::Event& evt, const VInputTag& collections) {
    std::vector<const reco::Candidate*> output;
    // Loop over collections
    for (size_t i = 0; i < collections.size(); ++i) {
      edm::Handle<edm::View<reco::Candidate> > handle;
      evt.getByLabel(collections[i], handle);
      // Loop over objects in current collection
      for (size_t j = 0; j < handle->size(); ++j) {
	const reco::Candidate& object = handle->at(j);
	output.push_back(&object);
      }
    }
    return output;
  }

  // Turn a set of InputTags into a collection of candidate pointers.
  std::vector<const reco::PFJet*> getPFCollections(const edm::Event& evt, const VInputTag& collections) {
    std::vector<const reco::PFJet*> output;
    // Loop over collections
    for (size_t i = 0; i < collections.size(); ++i) {
      edm::Handle<edm::View<reco::PFJet> > handle;
      evt.getByLabel(collections[i], handle);
      // Loop over objects in current collection
      for (size_t j = 0; j < handle->size(); ++j) {
	const reco::PFJet& object = handle->at(j);
	output.push_back(&object);
      }
    }
    return output;
  }

  void getValue(const edm::Event& evt, const edm::InputTag& tag, Float_t& pt, Float_t& phi) {
	edm::Handle<edm::View<reco::Candidate> > handle;
	evt.getByLabel(tag, handle);
	if(!handle.isValid()){
	  //std::cout << "MissingProduct " << tag.label() << std::endl;
	  pt  = 0;
	  phi = 0;
	}
	else{
	  pt  = handle->at(0).pt();
	  phi = handle->at(0).phi();
	}
    }

    void getSumEtL1(const edm::Event& evt, const edm::InputTag& tag, Float_t& sumet,bool upgrade) {
	if(!upgrade) {
	    edm::Handle<l1extra::L1EtMissParticleCollection> handle;
	    evt.getByLabel(tag, handle);
	    sumet = handle->at(0).etTotal();
	} else{
	    edm::Handle<edm::View<reco::Candidate> > handle;
	    evt.getByLabel(tag, handle);
	    sumet = handle->at(0).pt();
	}
    }

    const double PI = 3.14159265359;


    // inline float deltaPhi( float phi1, float phi2 ){

    // 	float const  PI        = ROOT::Math::Pi();
    // 	float const  TWOPI     = 2.*PI;
    // 	float dPhi = (phi1 - phi2);
    // 	while (dPhi >= PI) dPhi -= TWOPI;
    // 	while (dPhi < -PI) dPhi += TWOPI;
    // 	return dPhi;

    // }


  UInt_t calculateNJet(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
    
    UInt_t NJets(0);
    for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){
      if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
      NJets++;
    }
    return NJets;
  }
  Int_t calculateNJetBin(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
    // Returns the following:
    // 0 = No bin
    //-1 = eq2a,-2 = eq3a,-3 = ge4a
    // 1 = eq2j, 2 = eq3j, 3 = ge4j
    Int_t NJetBin(-1); 
    for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){
      if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
      NJetBin++;
      if (NJetBin == 3){ break; }
    }
    if ( NJetBin == -1){ NJetBin = 0; } // No jets

    // Asymmetric bin
    if ( NJetBin > 0 ){ // Require two jets
      if ( jets.at(1)->pt() < 100. ){ NJetBin *= -1; } // Bin is asymmetric
    }

    return NJetBin;
  }

  Int_t calculateHTBin( float ht, float alphaT ){

    Int_t HTBin(-1); 
    if ( ht >= 200){
      if (ht < 250)     { // 200 < HT < 250
	if (alphaT >= 0.65){ HTBin = 0; }
      }
      else if (ht < 300){ // 250 < HT < 300
	if (alphaT >= 0.60){ HTBin = 1; }
      }
      else if (ht < 350){ // 300 < HT < 350
	if (alphaT >= 0.55){ HTBin = 2; }
      }
      else if (ht < 400){ // 350 < HT < 400
	if (alphaT >= 0.53){ HTBin = 3; }
      }
      else if (ht < 500){ // 400 < HT < 500
	if (alphaT >= 0.52){ HTBin = 4; }
      }
      else{               // HT > 500
	if (alphaT >= 0.52){ HTBin = 5; }
      }

    }

    return HTBin;
  }


  std::pair<float, float> calculateMHT(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
    TLorentzVector mhtVec;

    for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){
      
      if ( jets.at(iJet)->pt() < jetThreshold ){ break; }

      TLorentzVector jet;
      jet.SetPtEtaPhiM( jets.at(iJet)->pt(), jets.at(iJet)->eta(), jets.at(iJet)->phi(), jets.at(iJet)->mass() );
      mhtVec += jet;
      
    }
    // Rotate by 180 degrees
    mhtVec = -mhtVec;

    return std::make_pair( mhtVec.Pt(), mhtVec.Phi() );
  }


//   std::pair<float, float> calculateAlphaTHT(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
     
//       // Momentum sums in transverse plane
//       float sum_et(0), sum_px(0), sum_py(0);

//       // check the size of the input collection
//       if (jets.size() <= 1){
// 	return std::make_pair(0., sum_et);
//       }


//       // Jet collection restricted to jet threshold
//       std::vector<float> jetPTNew;

//       for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){

// 	if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
// 	jetPTNew.push_back( jets.at(iJet)->pt() );

// 	sum_et += jets.at(iJet)->pt();
// 	sum_px += jets.at(iJet)->px();
// 	sum_py += jets.at(iJet)->py();

//       }
//       // check the size of the new input collection 
//       if (jetPTNew.size() <= 1){
// 	// empty jet collection, return AlphaT = 0 
// 	return std::make_pair( 0., sum_et);
//       }

//       // Minimum Delta Et for two pseudo-jets 
//       double min_delta_sum_et = sum_et;

//       for (unsigned int i = 0; i < (1U << (jetPTNew.size() - 1)); i++) { //@@ iterate through different combinations 
// 	double delta_sum_et = 0.;

// 	for (unsigned int j = 0; j < jetPTNew.size(); ++j) { //@@ iterate through jets 
// 	  if (i & (1U << j))
// 	    delta_sum_et -= jetPTNew[j];
// 	  else
// 	    delta_sum_et += jetPTNew[j];
// 	}
// 	delta_sum_et = std::abs(delta_sum_et);
// 	if (delta_sum_et < min_delta_sum_et) {
// 	  min_delta_sum_et = delta_sum_et;
// 	}

//       }

//       // Return a large value of alphaT 
//       if ( (sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py)) <= 0 ){
// 	return std::make_pair(11.,  sum_et);
//       }

//       // Alpha_T 
//       float alphaT = 0.5 * (sum_et - min_delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) );
//       if (alphaT > 10){
// 	alphaT = 10;
//       }
//      return std::make_pair( alphaT, sum_et );


//     }

}


void MakeTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& es) {

  
  // Get NVTX from simulation
  edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  NVTX = 0;
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    int BX = PVI->getBunchCrossing();
    if(BX == 0) {
      NVTX = PVI->getPU_NumInteractions();
      break;
    }
  }

  // Get gen info
  edm::Handle<GenEventInfoProduct> geninfo;  iEvent.getByLabel("generator",geninfo);
  std::auto_ptr<bool>   genInfoValid ( new bool( geninfo.isValid() && !geninfo->binningValues().empty()));
  std::auto_ptr<double> pthat (new double(*genInfoValid ? geninfo->binningValues()[0] : -1.));
  PThat = *pthat;






  // ------------------------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------------------------

  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByLabel(HLTResultsTag, hltresults);
  
  // Get the PAT TriggerEvent
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( "patTriggerEvent", triggerEvent );

  


  // Get a vector of all HLT paths
  std::vector<pat::TriggerPath> const* paths = triggerEvent->paths();
  
  // Find the full label of the chosen HLT path (i.e. with the version number)
  std::string full_name;


  // Iterate through paths run
  // ------------------------------------------------------------
  for (unsigned i = 0; i < paths->size(); ++i) {
    std::string name = paths->at(i).name();
    //std::cout << name << "\n";


    if ( hltPathFired.find( name ) != hltPathFired.end() ){
      //      std::cout << name << " found. Fired = " << paths->at(i).wasAccept() << "\n";
      hltPathFired[ name ] = paths->at(i).wasAccept();
    }

  }


  // ------------------------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------------------------
  
    // Input jets without eta or pT requirements
    // --------------------------------------------------------------------------------
    std::vector<const reco::Candidate*> gctCenUnskimmed                 = getCollections( iEvent, srcGctJetCentral_);
    std::vector<const reco::Candidate*> genJet4Unskimmed                = getCollections( iEvent, srcGen4Jet_);
    std::vector<const reco::Candidate*> hltAk4CaloUnskimmed             = getCollections( iEvent, srcHLTAk4Calo );
    std::vector<const reco::Candidate*> hltAk4PFUnskimmed               = getCollections( iEvent, srcHLTAk4PF );
    std::vector<const reco::Candidate*> recoAk4CaloUnskimmed            = getCollections( iEvent, srcAk4Calo );
    std::vector<const reco::Candidate*> recoAk4PFUnskimmed              = getCollections( iEvent, srcAk4PF   );

    std::vector<const reco::PFJet*>     hltAk4PFJetUnskimmed            = getPFCollections( iEvent, srcHLTAk4PF );

    std::vector<const reco::Candidate*> gctForUnskimmed                 = getCollections( iEvent, srcGctJetForward_);
    std::vector<const reco::Candidate*> genJet4ForUnskimmed             = getCollections( iEvent, srcGen4Jet_);
    std::vector<const reco::Candidate*> hltAk4CaloForUnskimmed          = getCollections( iEvent, srcHLTAk4Calo );
    std::vector<const reco::Candidate*> hltAk4PFForUnskimmed            = getCollections( iEvent, srcHLTAk4PF );
    std::vector<const reco::Candidate*> recoAk4CaloForUnskimmed         = getCollections( iEvent, srcAk4Calo );
    std::vector<const reco::Candidate*> recoAk4PFForUnskimmed           = getCollections( iEvent, srcAk4PF   );


 
    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
	jetPt[*iLvl] ->clear(); 
	jetPx[*iLvl] ->clear(); 
	jetPy[*iLvl] ->clear(); 
	jetPhi[*iLvl]->clear();
	jetEta[*iLvl]->clear();
    }


    // Clear previous event's objects
    // --------------------------------------------------------------------------------

    genElectronPt->clear();
    genElectronEta->clear();
    genElectronPhi->clear();
    genMuonPt ->clear();
    genMuonEta->clear();
    genMuonPhi->clear();
    genPhotonPt ->clear();
    genPhotonEta->clear();
    genPhotonPhi->clear();

    genMuonMatchedGenMuonPt          .clear();
    genMuonMatchedGenJetPt           .clear();
    genMuonMatchedHLTPFJetPt         .clear();
    genMuonMatchedHLTPFJetMuonEF     .clear();
    genMuonMatchedHLTPFJetElectronEF .clear();

    genElectronMatchedGenElectronPt      .clear();
    genElectronMatchedGenJetPt           .clear();
    genElectronMatchedHLTPFJetPt         .clear();
    genElectronMatchedHLTPFJetMuonEF     .clear();
    genElectronMatchedHLTPFJetElectronEF .clear();

    nIsoElectrons = 0;
    nIsoMuons     = 0;

    

    // Gen particles
    // --------------------------------------------------------------------------------

    genLeptonVeto      = false;
    genElectronVeto    = false;
    genMuonVeto        = false;
    genPhotonVeto      = false;
    genIsoLeptonVeto   = false;
    genIsoElectronVeto = false;


    std::vector<int> selLeptonIndices;
    int particleIndex(0);

    if (makeGenParticles){
      edm::Handle< std::vector<reco::GenParticle> > genParticles;
      iEvent.getByLabel(srcGenParticles_, genParticles);
      for(std::vector<reco::GenParticle>::const_iterator iter = genParticles->begin(); iter != genParticles->end(); ++iter){
    	const reco::GenParticle& genParticle = *iter;

    	double genParticlePt  = genParticle.p4().pt();
    	double genParticleEta = genParticle.p4().eta();
    	double genParticlePhi = genParticle.p4().phi();

    	// Electron
    	if(TMath::Abs(genParticle.pdgId()) == 11){
    	  if ( (genParticlePt >= genElectronMinPt) && (TMath::Abs(genParticleEta) <= genElectronMaxEta) ){
    	    genElectronPt ->push_back( genParticlePt  );
    	    genElectronEta->push_back( genParticleEta );
    	    genElectronPhi->push_back( genParticlePhi );
    	    genLeptonVeto   = true;
    	    genElectronVeto = true;
	    selLeptonIndices.push_back( particleIndex );
    	  }
    	} // End lepton requirement
    	// Muon
    	if(TMath::Abs(genParticle.pdgId()) == 13){
    	  if ( (genParticlePt >= genMuonMinPt) && (TMath::Abs(genParticleEta) <= genMuonMaxEta) ){
    	    genMuonPt ->push_back( genParticlePt  );
   	    genMuonEta->push_back( genParticleEta );
    	    genMuonPhi->push_back( genParticlePhi );
    	    genLeptonVeto = true;
	    genMuonVeto   = true;
	    selLeptonIndices.push_back( particleIndex );
    	  }
    	} // End Muon requirement

    	// Photon
    	if(TMath::Abs(genParticle.pdgId()) == 22){
    	  if ( (genParticlePt >= genPhotonMinPt) && (TMath::Abs(genParticleEta) <= genPhotonMaxEta) ){
    	    genPhotonPt ->push_back( genParticlePt  );
   	    genPhotonEta->push_back( genParticleEta );
    	    genPhotonPhi->push_back( genParticlePhi );
    	    genPhotonVeto = true;
    	  }
    	} // End Photon requirement

	particleIndex++;
      } // End loop
    } // End gen particles




    // Jet lepton cleaning
    // --------------------------------------------------------------------------------

    // std::cout << "\nBefore cleaning: " << genJet4Unskimmed.size() << "\t" << hltAk4PFUnskimmed.size() 
    // 	      << "\t" << nIsoMuons << "\t" << nIsoElectrons << "\n";
#ifdef LEPTON_XCLEANING
    if (genMuonVeto){
      //crosscleanIsolatedLeptons( genMuonPt,  genMuonEta,  genMuonPhi, genJet4Unskimmed, hltAk4PFUnskimmed, nIsoMuons, 0.3, 40. );
      crosscleanIsolatedLeptons( genMuonPt,  genMuonEta,  genMuonPhi, genJet4Unskimmed, hltAk4PFUnskimmed, hltAk4CaloUnskimmed, 
				 nIsoMuons, 0.3, 40.,
				 hltAk4PFJetUnskimmed,
				 genMuonMatchedGenMuonPt,
				 genMuonMatchedGenJetPt,
				 genMuonMatchedHLTPFJetPt,
				 genMuonMatchedHLTPFJetMuonEF,
				 genMuonMatchedHLTPFJetElectronEF);

    }
    if (genElectronVeto){
      // crosscleanIsolatedLeptons( genElectronPt,  genElectronEta,  genElectronPhi, genJet4Unskimmed, hltAk4PFUnskimmed, 
      // 				 nIsoElectrons, 0.3, 40. );
      crosscleanIsolatedLeptons( genElectronPt,  genElectronEta,  genElectronPhi, genJet4Unskimmed, hltAk4PFUnskimmed,  hltAk4CaloUnskimmed,
				 nIsoElectrons, 0.3, 40. ,
				 hltAk4PFJetUnskimmed,
				 genElectronMatchedGenElectronPt,
                                 genElectronMatchedGenJetPt,
                                 genElectronMatchedHLTPFJetPt,
                                 genElectronMatchedHLTPFJetMuonEF,
                                 genElectronMatchedHLTPFJetElectronEF);

    }
    // std::cout << "\nAfter cleaning: " << genJet4Unskimmed.size() << "\t" << hltAk4PFUnskimmed.size() 
    // 	      << "\t" << nIsoMuons << "\t" << nIsoElectrons << "\n";
#endif


    // std::vector<int>   matchedGenJetIndexDeltaR, matchedHLTJetIndexDeltaR;
    // std::vector<float> lepGenJetMinDeltaR,    lepHLTJetMinDeltaR;

    // if (genMuonVeto){

    //   lepJetDeltaR( genMuonEta, genMuonPhi, genJet4Unskimmed,   matchedGenJetIndexDeltaR, lepGenJetMinDeltaR );
    //   lepJetDeltaR( genMuonEta, genMuonPhi, hltAk4PFUnskimmed, matchedHLTJetIndexDeltaR, lepHLTJetMinDeltaR );

    //   for (uint i = 0; i < matchedGenJetIndexDeltaR.size(); ++i){

    // 	int genIndex = matchedGenJetIndexDeltaR[i];
    // 	if (genIndex == -1){ continue; }

    // 	// std::cout << "( " << genAk4.at(index)->eta() << ", " << genAk4.at(index)->phi() << " )\t" 
    // 	// 	  << "( " << (*genMuonEta)[i]        << ", " << (*genMuonPhi)[i]        << " )\t"
    // 	// 	  << lepJetMinDeltaR[i] << "\n";
	
    // 	float deltaR = lepGenJetMinDeltaR[i];
    // 	if ( deltaR < 0.3 ){ // Lepton matched to genJet

    // 	  float deltaPtGen = genJet4Unskimmed.at(genIndex)->pt() - (*genMuonPt)[i]; 

    // 	  // reco::Candidate*               jet;
    // 	  // reco::Candidate::LorentzVector totalP4 = jet->p4();
  
    // 	  // // Get muon from genParticles with index in vector
    // 	  // jetP4 -= (*muon)->p4();
    // 	  // jet->setP4(totalP4);

    //       if (deltaPtGen < 40.){ 
    //         // Remove genjet 
    // 	    //
    // 	    //	    
  
    // 	    if (i < matchedHLTJetIndexDeltaR.size() ){ // Check whether HLT reconstructed lepton as jet
    // 	      int hltIndex = matchedHLTJetIndexDeltaR[i];
    // 	      if (hltIndex == -1){ continue; }	    

    // 	      //genMuonMatchedGenMuonIndex.push_back( i );
    // 	      genMuonMatchedGenMuonPt.push_back(  (*genMuonPt)[i] );
    // 	      genMuonMatchedGenJetPt.push_back( genJet4Unskimmed.at(genIndex)->pt() );
    // 	      genMuonMatchedHLTPFJetPt.push_back( hltAk4PFJetUnskimmed.at(hltIndex)->pt() );
    // 	      genMuonMatchedHLTPFJetMuonEF.push_back( hltAk4PFJetUnskimmed.at(hltIndex)->muonEnergyFraction() );
    // 	      genMuonMatchedHLTPFJetElectronEF.push_back(hltAk4PFJetUnskimmed.at(hltIndex)->electronEnergyFraction() );
	      
    // 	      // //	      float deltaPtHLT = hltAk4PFUnskimmed.at(hltIndex)->pt() - (*genMuonPt)[i];
    // 	      // float deltaPtHLT = hltAk4PFJetUnskimmed.at(hltIndex)->pt() - (*genMuonPt)[i];
	      
    // 	      // 	std::cout << "GEN:\n \tgenJetPt = " << genJet4Unskimmed.at(genIndex)->pt() 
    // 	      // 		  << "\tGenMuonPt = " << (*genMuonPt)[i] 
    // 	      // 		  << "\tdeltaPt = " << deltaPtGen 
    // 	      // 		  << "\n";
		
    // 	      // 	std::cout << "\tGenJet( " << genJet4Unskimmed.at(genIndex)->eta() << ", " << genJet4Unskimmed.at(genIndex)->phi() << " )\t" 
    // 	      // 		  << "GenMuon( " << (*genMuonEta)[i]                     << ", " << (*genMuonPhi)[i]        << " )\tdeltaR = "
    // 	      // 		  << lepGenJetMinDeltaR[i] << "\n";
		
    // 	      // 	std::cout << "HLT:\n \thltJetPt = " << hltAk4PFJetUnskimmed.at(hltIndex)->pt() 
    // 	      // 		  << "\tGenMuonPt = " << (*genMuonPt)[i] 
    // 	      // 		  << "\tdeltaPt = " << deltaPtHLT << "\n"
    // 	      // 		  << "\thltJet( " << hltAk4PFJetUnskimmed.at(hltIndex)->eta() 
    // 	      // 		  << ", " << hltAk4PFJetUnskimmed.at(hltIndex)->phi() << " )\t" 
    // 	      // 		  << "GenMuon( " << (*genMuonEta)[i]        << ", " << (*genMuonPhi)[i]        << " )\tdeltaR = "
    // 	      // 		  << lepHLTJetMinDeltaR[i] << "\n"
    // 	      // 		  << "\thltMuonPt = "   << hltAk4PFJetUnskimmed.at(hltIndex)->muonEnergy()*TMath::Sin(hltAk4PFJetUnskimmed.at(hltIndex)->theta())
    // 	      // 		  << "\n\n"

		
    // 	      // 		  << "\thltJetEn  = "   << hltAk4PFJetUnskimmed.at(hltIndex)->energy() << "\n"
    // 	      // 	          << "\thltMuonEn = "   << hltAk4PFJetUnskimmed.at(hltIndex)->muonEnergy()
    // 	      // 		  << "\thltMuonEF = "   << hltAk4PFJetUnskimmed.at(hltIndex)->muonEnergyFraction() 
    // 	      // 		  << "\thltMuonMult = " << hltAk4PFJetUnskimmed.at(hltIndex)->muonMultiplicity() << "\n"
    // 	      // 	          << "\thltEleEn  = "   << hltAk4PFJetUnskimmed.at(hltIndex)->electronEnergy()
    // 	      // 		  << "\thltEleEF  = "   << hltAk4PFJetUnskimmed.at(hltIndex)->electronEnergyFraction() 
    // 	      // 		  << "\thltEleMult  = " << hltAk4PFJetUnskimmed.at(hltIndex)->electronMultiplicity() 
    // 	      // 	          << "\n\thltCHEn   = "   << hltAk4PFJetUnskimmed.at(hltIndex)->chargedHadronEnergy()
    // 	      // 		  << "\thltCHEF   = "     << hltAk4PFJetUnskimmed.at(hltIndex)->chargedHadronEnergyFraction() 
    // 	      // 		  << "\thltCHMult   = "   << hltAk4PFJetUnskimmed.at(hltIndex)->chargedHadronMultiplicity() 
    // 	      // 	          // << "\n\thltCEMEn  = "   << hltAk4PFJetUnskimmed.at(hltIndex)->chargedEmEnergy()
    // 	      // 		  // << "\thltCEMEF  = "     << hltAk4PFJetUnskimmed.at(hltIndex)->chargedEmEnergyFraction() 
    // 	      // 	          << "\n\thltNHEn   = "   << hltAk4PFJetUnskimmed.at(hltIndex)->neutralHadronEnergy()
    // 	      // 		  << "\thltNHEF   = "     << hltAk4PFJetUnskimmed.at(hltIndex)->neutralHadronEnergyFraction() 
    // 	      // 		  << "\thltNHMult   = "   << hltAk4PFJetUnskimmed.at(hltIndex)->neutralHadronMultiplicity() 
    // 	      // 	          << "\n\thltNEMEn  = "   << hltAk4PFJetUnskimmed.at(hltIndex)->neutralEmEnergy()
    // 	      // 		  << "\thltNEMEF  = "     << hltAk4PFJetUnskimmed.at(hltIndex)->neutralEmEnergyFraction() 
    // 	      // 		  << "\n";
		
    // 	      // 	std::cout << "\n\n";

    // 	      // if (deltaPtHLT > 40.){
    // 	      // 	// REMOVE ENERGY
		
    
    // 	      // }
    // 	      // else{
    // 	      // 	// REMOVE JET
    // 	      // }

    // 	    }
	    
    // 	  }

    // 	} // Matched lepton to jet

    //   } // DeltaR loop

    // }




    // Skim jet collections
    // ----------------------------------------

    // Central jets
    // --------------------
    std::vector<const reco::Candidate*> gctCen              = skimJets(gctCenUnskimmed,              minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> genAk4              = skimJets(genJet4Unskimmed,             minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltAk4Calo          = skimJets(hltAk4CaloUnskimmed,          minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltAk4PF            = skimJets(hltAk4PFUnskimmed,            minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> recoAk4Calo         = skimJets(recoAk4CaloUnskimmed,         minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> recoAk4PF           = skimJets(recoAk4PFUnskimmed,           minPt, minEtaCen, maxEtaCen );


    // Forward jets
    // --------------------
    std::vector<const reco::Candidate*> gctFor                 = skimJets(gctForUnskimmed,                 minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> genAk4For              = skimJets(genJet4ForUnskimmed,             minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltAk4CaloFor          = skimJets(hltAk4CaloForUnskimmed,          minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltAk4PFFor            = skimJets(hltAk4PFForUnskimmed,            minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> recoAk4CaloFor         = skimJets(recoAk4CaloForUnskimmed,         minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> recoAk4PFFor           = skimJets(recoAk4PFForUnskimmed,           minPt, minEtaFor, maxEtaFor );


    // Jets
    // --------------------------------------------------------------------------------


    // Forward jet
    genAk4ForMaxPt    = 0;
    recoAk4PFForMaxPt = 0;
    hltAk4PFForMaxPt  = 0;
    if ( genAk4For.size()   > 0){    genAk4ForMaxPt =    genAk4For.at(0)->pt(); }
    if (recoAk4PFFor.size() > 0){ recoAk4PFForMaxPt = recoAk4PFFor.at(0)->pt(); }
    if ( hltAk4PFFor.size() > 0){  hltAk4PFForMaxPt =  hltAk4PFFor.at(0)->pt(); }

    // Lead jet 
    genAk4LeadJetPt= 0;
    recoAk4PFLeadJetPt= 0;
    hltAk4PFLeadJetPt= 0;
    hltAk4CaloLeadJetPt= 0;
    if ( genAk4.size()     > 0){ genAk4LeadJetPt      =     genAk4.at(0)->pt(); }
    if (recoAk4PF.size()   > 0){ recoAk4PFLeadJetPt   =  recoAk4PF.at(0)->pt(); }
    if ( hltAk4PF.size()   > 0){ hltAk4PFLeadJetPt    =   hltAk4PF.at(0)->pt(); }
    if ( hltAk4Calo.size() > 0){ hltAk4CaloLeadJetPt  = hltAk4Calo.at(0)->pt(); }

    // Second jet 
    genAk4SecondJetPt= 0;
    recoAk4PFSecondJetPt= 0;
    hltAk4PFSecondJetPt= 0;
    hltAk4CaloSecondJetPt= 0;
    if ( genAk4.size()     > 1){ genAk4SecondJetPt      =     genAk4.at(1)->pt(); }
    if (recoAk4PF.size()   > 1){ recoAk4PFSecondJetPt   =  recoAk4PF.at(1)->pt(); }
    if ( hltAk4PF.size()   > 1){ hltAk4PFSecondJetPt    =   hltAk4PF.at(1)->pt(); }
    if ( hltAk4Calo.size() > 1){ hltAk4CaloSecondJetPt  = hltAk4Calo.at(1)->pt(); }

    // Dijet avg 
    genAk4DijetAvgPt= 0;
    recoAk4PFDijetAvgPt= 0;
    hltAk4PFDijetAvgPt= 0;
    hltAk4CaloDijetAvgPt= 0;
    if ( genAk4.size()     > 1){ genAk4DijetAvgPt      = 0.5*(genAk4.at(0)->pt()     + genAk4.at(1)->pt()); }
    if (recoAk4PF.size()   > 1){ recoAk4PFDijetAvgPt   = 0.5*(recoAk4PF.at(0)->pt()  + recoAk4PF.at(1)->pt()); }
    if ( hltAk4PF.size()   > 1){ hltAk4PFDijetAvgPt    = 0.5*(hltAk4PF.at(0)->pt()   + hltAk4PF.at(1)->pt()); }
    if ( hltAk4Calo.size() > 1){ hltAk4CaloDijetAvgPt  = 0.5*(hltAk4Calo.at(0)->pt() + hltAk4Calo.at(1)->pt()); }



    // ********************************************************************************
    // *                                  Energy sums                                 *
    // ********************************************************************************
    getValue(iEvent,   srcGctMht_, mhtPt_["gct"], mhtPhi_["gct"]);
    getValue(iEvent,   srcGctMet_, metPt_["gct"], metPhi_["gct"]);
    getSumEtL1(iEvent, srcGctMht_, ht_["gct"], false);
    getSumEtL1(iEvent, srcGctMet_, et_["gct"], false); 
    if( ht_["gct"] > 0.){
      mhtDivHt_["gct"] = mhtPt_["gct"]/ht_["gct"];
    }
    else{
      mhtDivHt_["gct"] = 0;
    }

    // genMet
    getValue(iEvent, srcGenMetCalo_,                metPt_["genMetCalo"],                metPhi_["genMetCalo"]);
    getValue(iEvent, srcGenMetCaloAndNonPrompt_,    metPt_["genMetCaloAndNonPrompt"],    metPhi_["genMetCaloAndNonPrompt"]);
    getValue(iEvent, srcGenMetTrue_,                metPt_["genMetTrue"],                metPhi_["genMetTrue"]);
    // hltMet
    getValue(iEvent, srcHLTMetCalo_,                metPt_["hltMetCalo"],                metPhi_["hltMetCalo"]);
    getValue(iEvent, srcHLTMetCleanCalo_,           metPt_["hltMetCleanCalo"],           metPhi_["hltMetCleanCalo"]);
    getValue(iEvent, srcHLTMetCleanUsingJetIDCalo_, metPt_["hltMetCleanUsingJetIDCalo"], metPhi_["hltMetCleanUsingJetIDCalo"]);
    getValue(iEvent, srcHLTMetPF_,                  metPt_["hltMetPF"],                  metPhi_["hltMetPF"]);
    // hltMHT
    getValue(iEvent, srcHLTMhtCalo_,             mhtPt_["hltMhtCalo"],             mhtPhi_["hltMhtCalo"]);
    getValue(iEvent, srcHLTMhtPF_,               mhtPt_["hltMhtPF"],               mhtPhi_["hltMhtPF"]);


    // ********************************************************************************
    // *                           Loop over Jet collections                          *
    // ********************************************************************************

    // HT and AlphaT
    genAk4AlphaTHT40              = calculateAlphaTHT( genAk4,      40.);
    hltAk4CaloAlphaTHT40          = calculateAlphaTHT( hltAk4Calo,  40.);
    hltAk4PFAlphaTHT40            = calculateAlphaTHT( hltAk4PF,    40.);
    recoAk4CaloAlphaTHT40         = calculateAlphaTHT( recoAk4Calo, 40.);
    recoAk4PFAlphaTHT40           = calculateAlphaTHT( recoAk4PF,   40.);

    genAk4AlphaTHT50              = calculateAlphaTHT( genAk4,      50.);
    hltAk4CaloAlphaTHT50          = calculateAlphaTHT( hltAk4Calo,  50.);
    hltAk4PFAlphaTHT50            = calculateAlphaTHT( hltAk4PF,    50.);
    recoAk4CaloAlphaTHT50         = calculateAlphaTHT( recoAk4Calo, 50.);
    recoAk4PFAlphaTHT50           = calculateAlphaTHT( recoAk4PF,   50.);

    // Dynamic HT and AlphaT 
    genAk4DynamicAlphaTHT40    = calculateDynamicAlphaTPairs( genAk4,    dynamicJetThreshold );
    hltAk4PFDynamicAlphaTHT40  = calculateDynamicAlphaTPairs( hltAk4PF,  dynamicJetThreshold );
    recoAk4PFDynamicAlphaTHT40 = calculateDynamicAlphaTPairs( recoAk4PF, dynamicJetThreshold );


    // ********************************************************************************
    // Calculate: Biased deltaPhi, tagJet index, AlphaTVector, BetaTScalar, BetaTVector
    // ********************************************************************************
    // genAk4BiasedDPhi    = calculateBiasedDeltaPhi( genAk4,    genAk4BiasedDPhiIndex);
    // recoAk4PFBiasedDPhi = calculateBiasedDeltaPhi( recoAk4PF, recoAk4PFBiasedDPhiIndex);
    // hltAk4PFBiasedDPhi  = calculateBiasedDeltaPhi( hltAk4PF,  hltAk4PFBiasedDPhiIndex);

    calculateBDPhiAlphaBetaT( genAk4,    40., genAk4BiasedDPhi,  genAk4BiasedDPhiIndex, 
			      genAk4VecAlphaT40,    genAk4ScaBetaT40,    genAk4VecBetaT40);
    calculateBDPhiAlphaBetaT( recoAk4PF, 40., recoAk4PFBiasedDPhi, recoAk4PFBiasedDPhiIndex,
			      recoAk4PFVecAlphaT40, recoAk4PFScaBetaT40, recoAk4PFVecBetaT40);
    calculateBDPhiAlphaBetaT( hltAk4PF,  40., hltAk4PFBiasedDPhi,  hltAk4PFBiasedDPhiIndex, 
			      hltAk4PFVecAlphaT40,  hltAk4PFScaBetaT40,  hltAk4PFVecBetaT40);


    // ******************************************************************************** 
    // Calculate: DeltaR of leading jets L1-Gen, HLT-Gen
    // ******************************************************************************** 
    L1GenDeltaR  = leadL1GenDeltaR(    gctCenUnskimmed,  gctForUnskimmed, genJet4Unskimmed );
    HLTGenDeltaR = leadHLTGenDeltaR( hltAk4PFUnskimmed, genJet4Unskimmed );
    hpuVeto      = (L1GenDeltaR < 0.5);


    // MHT
    genAk4MHT40              = calculateMHT( genAk4,      40.);
    hltAk4CaloMHT40          = calculateMHT( hltAk4Calo,  40.);
    hltAk4PFMHT40            = calculateMHT( hltAk4PF,    40.);
    recoAk4CaloMHT40         = calculateMHT( recoAk4Calo, 40.);
    recoAk4PFMHT40           = calculateMHT( recoAk4PF,   40.);

    genAk4ForMHT40           = calculateMHT( genAk4For,      40.);
    hltAk4PFForMHT40         = calculateMHT( hltAk4PFFor,    40.);
    recoAk4PFForMHT40        = calculateMHT( recoAk4PFFor,   40.);

    // AlphaT prime
    genAk4_AlphaTPrime40        = calculateAlphaTPrime( genAk4MHT40.first,      genAk4AlphaTHT40.second );
    hltAk4Calo_AlphaTPrime40    = calculateAlphaTPrime( hltAk4CaloMHT40.first,  hltAk4CaloAlphaTHT40.second );
    hltAk4PF_AlphaTPrime40      = calculateAlphaTPrime( hltAk4PFMHT40.first,    hltAk4PFAlphaTHT40.second );
    recoAk4Calo_AlphaTPrime40   = calculateAlphaTPrime( recoAk4CaloMHT40.first, recoAk4CaloAlphaTHT40.second );
    recoAk4PF_AlphaTPrime40     = calculateAlphaTPrime( recoAk4PFMHT40.first,   recoAk4PFAlphaTHT40.second );


    // ****************************************
    // MHT-MET deltaPhi
    // ****************************************
    hltMetCaloPFMht40_DeltaPhi = fabs( deltaPhi( metPhi_["hltMetCalo"], hltAk4PFMHT40.second ) );
    genMetCaloMht40_DeltaPhi   = fabs( deltaPhi( metPhi_["genMetCalo"], genAk4MHT40.second )   );



    // Count analysis jet multiplicity
    genAk4NJet40              = calculateNJet( genAk4,      40.);
    hltAk4CaloNJet40          = calculateNJet( hltAk4Calo,  40.);
    hltAk4PFNJet40            = calculateNJet( hltAk4PF,    40.);
    recoAk4CaloNJet40         = calculateNJet( recoAk4Calo, 40.);
    recoAk4PFNJet40           = calculateNJet( recoAk4PF,   40.);

    genAk4NJet50              = calculateNJet( genAk4,      50.);
    hltAk4CaloNJet50          = calculateNJet( hltAk4Calo,  50.);
    hltAk4PFNJet50            = calculateNJet( hltAk4PF,    50.);
    recoAk4CaloNJet50         = calculateNJet( recoAk4Calo, 50.);
    recoAk4PFNJet50           = calculateNJet( recoAk4PF,   50.);

    // Analysis bins
    // ----------------------------------------
    
    genAk4NJetBin40           = calculateNJetBin( genAk4,      40.);
    hltAk4CaloNJetBin40       = calculateNJetBin( hltAk4Calo,  40.);
    hltAk4PFNJetBin40         = calculateNJetBin( hltAk4PF,    40.);
    recoAk4CaloNJetBin40      = calculateNJetBin( recoAk4Calo, 40.);
    recoAk4PFNJetBin40        = calculateNJetBin( recoAk4PF,   40.);


    genAk4HTBin40           = calculateHTBin( genAk4AlphaTHT40.second,      genAk4AlphaTHT40.first );
    hltAk4CaloHTBin40       = calculateHTBin( hltAk4CaloAlphaTHT40.second,  hltAk4CaloAlphaTHT40.first );
    hltAk4PFHTBin40         = calculateHTBin( hltAk4PFAlphaTHT40.second,    hltAk4PFAlphaTHT40.first );
    recoAk4CaloHTBin40      = calculateHTBin( recoAk4CaloAlphaTHT40.second, recoAk4CaloAlphaTHT40.first );
    recoAk4PFHTBin40        = calculateHTBin( recoAk4PFAlphaTHT40.second,   recoAk4PFAlphaTHT40.first );



    // Store jet collections
    // ----------------------------------------

    storeJet( "gctCen",                 gctCen,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk4",                 genAk4,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4Calo",             hltAk4Calo,             jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4PF",               hltAk4PF,               jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4Calo",            recoAk4Calo,            jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4PF",              recoAk4PF,              jetPt, jetPx, jetPy, jetEta, jetPhi );

    storeJet( "gctFor",                 gctFor,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk4For",              genAk4For,              jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4CaloFor",          hltAk4CaloFor,          jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4PFFor",            hltAk4PFFor,            jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4CaloFor",         recoAk4CaloFor,         jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4PFFor",           recoAk4PFFor,           jetPt, jetPx, jetPy, jetEta, jetPhi );

    
    // L1 seeds
    // ----------------------------------------

    L1HTT175        = ( ht_["gct"]    >= 175);
    L1ETM70         = ( metPt_["gct"] >= 70);
    L1HTT175OrETM70 = (L1HTT175 || L1ETM70);

    L1Jet_DPhi  = -1;
    if (jetPhi["gctCen"]->size() > 1){
      L1Jet_DPhi      = fabs( deltaPhi( (*jetPhi["gctCen"])[0], (*jetPhi["gctCen"])[1] ) );
    }

    
    // Perform jet matching
    // ----------------------------------------

    genMatchedAk4HLTPF  = matchJetCollections( jetPt["genAk4"],   jetEta["genAk4"],   jetPhi["genAk4"],
    					       jetPt["hltAk4PF"], jetEta["hltAk4PF"], jetPhi["hltAk4PF"], 20., 0.25 );
    
    // genMatchedAk4RECOPF = matchJetCollections( jetPt["genAk4"],    jetEta["genAk4"],    jetPhi["genAk4"],
    // 					       jetPt["recoAk4PF"], jetEta["recoAk4PF"], jetPhi["recoAk4PF"], 20., 0.25 );
    


    // Fill the trees
    tree->Fill();
}






std::vector< const reco::Candidate*> 
MakeTrees::skimJets(const std::vector<const reco::Candidate*>& inputJets, double minPt, double minEta, double maxEta){

    std::vector <const reco::Candidate*> skimmedJets; 
    for ( unsigned int iJet = 0; iJet < inputJets.size(); iJet++ ){
      if ( (inputJets.at(iJet)->pt() >= minPt) 
	   && ( fabs(inputJets.at(iJet)->eta()) >= minEta) 
	   && ( fabs(inputJets.at(iJet)->eta()) <  maxEta) ){
	    skimmedJets.push_back(inputJets.at(iJet));
	}
    }

    return skimmedJets;
}




// std::vector <int>
// MakeTrees::correctLeptonHLTJets(const std::vector<const reco::Candidate*>& inputJets,
// 				std::vector <int> skimmedJetIndices,
// 				std::vector<int> matchedGenJetIndexDeltaR, std::vector<float> lepGenJetMinDeltaR,
// 				double maxDeltaR, double minPt, int& isoLepCount ){

	  
// 	    if (i < matchedHLTJetIndexDeltaR.size() ){ // Check whether HLT reconstructed lepton as jet
// 	      int hltIndex = matchedHLTJetIndexDeltaR[i];
// 	      if (hltIndex == -1){ continue; }	    
	      
// 	      float deltaPtHLT = hltAk4PFUnskimmed.at(hltIndex)->pt() - (*genMuonPt)[i];

// 	      if (deltaPtHLT > minPt){
// 		// REMOVE ENERGY
		
// 		std::cout << "GEN: " << genJet4Unskimmed.at(genIndex)->pt() << "\t" << (*genMuonPt)[i] << "\t" << deltaPtGen << "\n";
		
// 		std::cout << "( " << genJet4Unskimmed.at(genIndex)->eta() << ", " << genJet4Unskimmed.at(genIndex)->phi() << " )\t" 
// 			  << "( " << (*genMuonEta)[i]                     << ", " << (*genMuonPhi)[i]        << " )\t"
// 			  << lepGenJetMinDeltaR[i] << "\n";
		
// 		std::cout << "HLT: " << hltAk4PFUnskimmed.at(hltIndex)->pt() << "\t" 
// 			  << (*genMuonPt)[i] << "\t" << deltaPtHLT << "\n";
		
// 		std::cout << "( " << hltAk4PFUnskimmed.at(hltIndex)->eta() << ", " << hltAk4PFUnskimmed.at(hltIndex)->phi() << " )\t" 
// 			  << "( " << (*genMuonEta)[i]        << ", " << (*genMuonPhi)[i]        << " )\t"
// 			  << lepHLTJetMinDeltaR[i] << "\n";
// 		std::cout << "\n\n";
    
// 	      }
// 	      else{
// 		// REMOVE JET
// 	      }
// 	    }
// 	  }
//       }
//     }



 




// }







void
MakeTrees::storeJet( TString jetCollName, const std::vector<const reco::Candidate*>& jetColl,
		     std::map<TString,std::vector<Float_t>* >& pt,
		     std::map<TString,std::vector<Float_t>* >& px,
		     std::map<TString,std::vector<Float_t>* >& py,
		     std::map<TString,std::vector<Float_t>* >& eta,
		     std::map<TString,std::vector<Float_t>* >& phi){

  for ( std::vector<const reco::Candidate*>::const_iterator itr = jetColl.begin(); itr != jetColl.end(); ++itr ){
    pt [jetCollName]->push_back( (*itr)->pt() );
    px [jetCollName]->push_back( (*itr)->px() );
    py [jetCollName]->push_back( (*itr)->py() );
    eta[jetCollName]->push_back( (*itr)->eta() );
    phi[jetCollName]->push_back( (*itr)->phi() );
  }

}


float
MakeTrees::leadL1GenDeltaR( const std::vector<const reco::Candidate*>& gctCen, const std::vector<const reco::Candidate*>& gctFor, const std::vector<const reco::Candidate*>& genAk4){

  const reco::Candidate* leadL1Jet  = NULL;
  const reco::Candidate* leadGenJet = NULL;
  float deltaR(20.);

  // Get lead genJet in acceptance
  for ( std::vector<const reco::Candidate*>::const_iterator itr = genAk4.begin(); itr != genAk4.end(); ++itr ){
    if ( fabs((*itr)->eta()) < 5. ){ leadGenJet = (*itr); break; }
  }
  if ( leadGenJet == NULL )    { return deltaR; }
  if ( leadGenJet->pt() < 10. ){ return deltaR; }

  // Determine lead L1 jet
  // ----------------------------------------
  float cenMaxPT(-1), forMaxPT(-1);
  if (gctCen.size() > 0){ cenMaxPT = gctCen.at(0)->pt(); }
  if (gctFor.size() > 0){ forMaxPT = gctFor.at(0)->pt(); }
  if (( cenMaxPT == -1 ) && ( forMaxPT == -1 )){ return deltaR; }
  if (cenMaxPT > forMaxPT){
    leadL1Jet = gctCen.at(0);
    deltaR    = reco::deltaR( leadL1Jet->eta(), leadL1Jet->phi(), leadGenJet->eta(), leadGenJet->phi() );
  }
  else if (cenMaxPT < forMaxPT){
    leadL1Jet = gctFor.at(0);
    deltaR    = reco::deltaR( leadL1Jet->eta(), leadL1Jet->phi(), leadGenJet->eta(), leadGenJet->phi() );
  }
  else{
    leadL1Jet     = gctCen.at(0);
    deltaR        = reco::deltaR( leadL1Jet->eta(), leadL1Jet->phi(), leadGenJet->eta(), leadGenJet->phi() );
    leadL1Jet     = gctFor.at(0);
    float deltaR2 = reco::deltaR( leadL1Jet->eta(), leadL1Jet->phi(), leadGenJet->eta(), leadGenJet->phi() );
    if (deltaR2 < deltaR){ deltaR = deltaR2; }
  }

  return deltaR;

}

float
MakeTrees::leadHLTGenDeltaR( const std::vector<const reco::Candidate*>& hltAk4,
			     const std::vector<const reco::Candidate*>& genAk4){

  const reco::Candidate* leadHLTJet = NULL;
  const reco::Candidate* leadGenJet = NULL;
  float deltaR(20.);

  // Get lead genJet in acceptance
  for ( std::vector<const reco::Candidate*>::const_iterator itr = genAk4.begin(); itr != genAk4.end(); ++itr ){
    if ( fabs((*itr)->eta()) < 5. ){
      leadGenJet = (*itr); break;
    }
  }
  if ( leadGenJet == NULL )    { return deltaR; }
  if ( leadGenJet->pt() < 10. ){ return deltaR; }

  // Get lead hltJet
  if (hltAk4.size() > 0){ leadHLTJet = hltAk4.at(0); }
  else{ return deltaR; }
  deltaR    = reco::deltaR( leadHLTJet->eta(), leadHLTJet->phi(), leadGenJet->eta(), leadGenJet->phi() );

  return deltaR;

}







// Not one-to-one matching
void
MakeTrees::lepJetDeltaR( std::vector<Float_t>* leptonEta, std::vector<Float_t>* leptonPhi, const std::vector<const reco::Candidate*>& jet,
			 std::vector<int>   &matchedJetIndexDeltaR, std::vector<float> &lepJetMinDeltaR ){

  // No objects
  if ( (leptonEta->size() == 0) || (jet.size() == 0)){ return; }
  float deltaRMin(20.), deltaR(20.);
  int jetIndex(0);
  int matchedJet(-1);

  // Lepton loop
  for (uint iLep = 0; iLep < leptonEta->size(); ++iLep){
    //for ( std::vector<const reco::Candidate*>::const_iterator lepItr = lepton.begin(); lepItr != lepton.end(); ++lepItr ){
    // Reset variables
    jetIndex   = 0;
    deltaRMin  = 20;
    matchedJet = -1;

    // Jet loop
    for ( std::vector<const reco::Candidate*>::const_iterator jetItr = jet.begin(); jetItr != jet.end(); ++jetItr ){
      deltaR = reco::deltaR( (*leptonEta)[iLep], (*leptonPhi)[iLep], (*jetItr)->eta(), (*jetItr)->phi() );

      if (deltaR < deltaRMin){
	matchedJet = jetIndex;
	deltaRMin  = deltaR;
      }


      // std::cout << "( " << (*jetItr)->eta() << ", " << (*jetItr)->phi() << " )\t"
      // 		<< "( " << (*leptonEta)[iLep]      << ", " << (*leptonPhi)[iLep]      << " )\t"
      // 		<< deltaR << "\n";

      jetIndex++;
    }


    // std::cout << "================================================================================\n"
    // 	      << "( " << jet.at(matchedJet)->eta() << ", " << jet.at(matchedJet)->phi() << " )\t" 
    // 	      << "( " << (*leptonEta)[iLep]      << ", " << (*leptonPhi)[iLep]      << " )\t" 
    // 	      << deltaRMin << "\n"
    // 	      << "================================================================================\n";

    matchedJetIndexDeltaR.push_back( matchedJet );
    lepJetMinDeltaR.push_back( deltaRMin );

  }

}


void 
MakeTrees::crosscleanIsolatedLeptons( std::vector<Float_t>* genLepPt, std::vector<Float_t>* genLepEta, std::vector<Float_t>* genLepPhi,
				      std::vector<const reco::Candidate*>& genJet, 
				      std::vector<const reco::Candidate*>& hltJet, 
				      std::vector<const reco::Candidate*>& hltCaloJet, 
				      int& nIsoLeptons, double maxDeltaR, double minPt,
				      std::vector<const reco::PFJet*>& hltPFJet,
				      std::vector<float>& genLeptonPt ,
				      std::vector<float>& genLeptonMatchedGenJetPt ,
				      std::vector<float>& genLeptonMatchedHLTPFJetPt ,
				      std::vector<float>& genLeptonMatchedHLTPFJetMuonEF ,
				      std::vector<float>& genLeptonMatchedHLTPFJetElectronEF ){




  // Clear number of isolated leptons
  nIsoLeptons = 0;

  std::vector<int>   genJetCleaningIndex,      hltJetCleaningIndex,      hltCaloJetCleaningIndex;
  std::vector<int>   matchedGenJetIndexDeltaR, matchedHLTJetIndexDeltaR, matchedHLTCaloJetIndexDeltaR;
  std::vector<float> lepGenJetMinDeltaR,       lepHLTJetMinDeltaR,       lepHLTCaloJetMinDeltaR;
  
  
  // Determine deltaR between jets and lepton
  lepJetDeltaR( genLepEta, genLepPhi, genJet,      matchedGenJetIndexDeltaR,     lepGenJetMinDeltaR );
  lepJetDeltaR( genLepEta, genLepPhi, hltJet,      matchedHLTJetIndexDeltaR,     lepHLTJetMinDeltaR );
  lepJetDeltaR( genLepEta, genLepPhi, hltCaloJet,  matchedHLTCaloJetIndexDeltaR, lepHLTCaloJetMinDeltaR );


    for (uint iLep = 0; iLep < matchedGenJetIndexDeltaR.size(); ++iLep){
      float lepPt = (*genLepPt)[iLep];

      int genIndex     = matchedGenJetIndexDeltaR[iLep];
      int hltIndex     = matchedHLTJetIndexDeltaR[iLep];      
      int hltCaloIndex = matchedHLTCaloJetIndexDeltaR[iLep];      
      if (genIndex     == -1){ continue; }
      if (hltIndex     == -1){ continue; }	    
      if (hltCaloIndex == -1){ continue; }	    
      
      float deltaRGen      = lepGenJetMinDeltaR[iLep];
      float deltaRHLT      = lepHLTJetMinDeltaR[iLep];
      float deltaRCaloHLT  = lepHLTJetMinDeltaR[iLep];
      float deltaPtGen = genJet.at(genIndex)->pt() - lepPt; 
      //float deltaPtHLT = hltJet.at(hltIndex)->pt() - lepPt;

      if ( (deltaRGen < maxDeltaR)&&(deltaPtGen < minPt) ){ // Lepton matched to genJet and deemed to be gen-isolated
	nIsoLeptons++;
	genJetCleaningIndex.push_back(     genIndex );
	// Temporary - for diagnostics
	genLeptonPt                       .push_back( lepPt );
	genLeptonMatchedGenJetPt          .push_back( genJet.at(genIndex)->pt() );
	    
	// Remove HLTJet if can be matched to GenLepton
	if ( (deltaRHLT < maxDeltaR) ){                     // Lepton matched to HLTJet
	  hltJetCleaningIndex.push_back(     hltIndex );	
	  genLeptonMatchedHLTPFJetPt        .push_back( hltJet.at(hltIndex)->pt() );
	  genLeptonMatchedHLTPFJetMuonEF    .push_back( hltPFJet.at(hltIndex)->muonEnergyFraction() );
	  genLeptonMatchedHLTPFJetElectronEF.push_back( hltPFJet.at(hltIndex)->electronEnergyFraction() );
	}
	if ( (deltaRCaloHLT < maxDeltaR) ){                     // Lepton matched to HLTJet
	  hltCaloJetCleaningIndex.push_back( hltCaloIndex );
	}

      } // End genJet matching requirements
    } // End lepton loop

    // Clean jets associated to isolated leptons
    removeLeptonsFromJets( genJet,     genJetCleaningIndex );
    removeLeptonsFromJets( hltJet,     hltJetCleaningIndex );
    removeLeptonsFromJets( hltCaloJet, hltCaloJetCleaningIndex );

    return;
}


void 
MakeTrees::removeLeptonsFromJets( std::vector <const reco::Candidate*>& inputJets, std::vector<int> cleaningIndices ){
  
  
  std::vector <const reco::Candidate*> skimmedJets; 
  for ( int iJet = 0; iJet < int(inputJets.size()); iJet++ ){ 
    bool skim = false; 
    for (auto cleaningIndex: cleaningIndices){ if ( cleaningIndex == iJet ){ skim = true; break; } } 
    if (!skim){ skimmedJets.push_back(inputJets.at(iJet)); }
  } 
  
  inputJets = skimmedJets;
}


 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakeTrees);
