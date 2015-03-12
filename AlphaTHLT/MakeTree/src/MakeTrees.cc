// UNCOMMENT TO RUN ON RECO
#define RECO
 
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
    bool genIsoLeptonVeto;
    bool genIsoElectronVeto;
    std::vector<Float_t>* genElectronPt;
    std::vector<Float_t>* genElectronEta;
    std::vector<Float_t>* genElectronPhi;
    std::vector<Float_t>* genMuonPt;
    std::vector<Float_t>* genMuonEta;
    std::vector<Float_t>* genMuonPhi;

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

};

MakeTrees::MakeTrees(const edm::ParameterSet& pset){

    // Initialize the ntuple builder
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("Ntuple", "Ntuple");


    lvl_.push_back("gctCen");
    lvl_.push_back("genAk4");
    //lvl_.push_back("genAk5");
    lvl_.push_back("hltAk4Calo");
    //lvl_.push_back("hltAk4CaloNoFastJet");
    lvl_.push_back("hltAk4PF");
    lvl_.push_back("recoAk4Calo");
    lvl_.push_back("recoAk4PF");


    lvl_.push_back("gctFor");
    lvl_.push_back("genAk4For");
    //lvl_.push_back("genAk5For");
    lvl_.push_back("hltAk4CaloFor");
    //lvl_.push_back("hltAk4CaloNoFastJetFor");
    lvl_.push_back("hltAk4PFFor");
    lvl_.push_back("recoAk4CaloFor");
    lvl_.push_back("recoAk4PFFor");

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

    // Forward jet max Pt
    genAk4ForMaxPt    = 0;
    recoAk4PFForMaxPt = 0;
    hltAk4PFForMaxPt  = 0;


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
    tree->Branch("genElectron_Pt",   "std::vector<float>", &genElectronPt);
    tree->Branch("genElectron_Eta",  "std::vector<float>", &genElectronEta);
    tree->Branch("genElectron_Phi",  "std::vector<float>", &genElectronPhi);		
    tree->Branch("genMuon_Pt",       "std::vector<float>", &genMuonPt);
    tree->Branch("genMuon_Eta",      "std::vector<float>", &genMuonEta);
    tree->Branch("genMuon_Phi",      "std::vector<float>", &genMuonPhi);		


    tree->Branch("genMatchedAk4HLTPF", "std::vector<int>", &genMatchedAk4HLTPF);


    // Store HLT paths
    // ------------------------------------------------------------

    // hltPathNames.push_back("HLT_CaloJet20_v1");
    // hltPathNames.push_back("HLT_PFJet20_v1");
    // hltPathNames.push_back("HLT_HT100_v1");
    // hltPathNames.push_back("HLT_PFHT100_v1");
    
    // hltPathNames.push_back("HLT_HT200_AlphaT0p57_NoL1_v1");
    // hltPathNames.push_back("HLT_HT250_AlphaT0p55_NoL1_v1");
    // hltPathNames.push_back("HLT_HT300_AlphaT0p53_NoL1_v1");
    // hltPathNames.push_back("HLT_HT350_AlphaT0p52_NoL1_v1");
    // hltPathNames.push_back("HLT_HT400_AlphaT0p51_NoL1_v1");
    // hltPathNames.push_back("HLT_HT200_AlphaT0p57_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT250_AlphaT0p55_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT300_AlphaT0p53_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT350_AlphaT0p52_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT400_AlphaT0p51_L1HTT175OrETM70_v1");
    
    // hltPathNames.push_back("HLT_HT200_AlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT250_AlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT300_AlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT350_AlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT200_AlphaT0p5_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT250_AlphaT0p5_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT300_AlphaT0p5_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT350_AlphaT0p5_L1HTT175OrETM70_v1");

    // hltPathNames.push_back("HLT_HT200_PFAlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT250_PFAlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT300_PFAlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT350_PFAlphaT0p5_NoL1_v1");
    // hltPathNames.push_back("HLT_HT200_PFAlphaT0p5_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT250_PFAlphaT0p5_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT300_PFAlphaT0p5_L1HTT175OrETM70_v1");
    // hltPathNames.push_back("HLT_HT350_PFAlphaT0p5_L1HTT175OrETM70_v1");
          
    // hltPathNames.push_back("HLT_PFHT350_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT600_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT350_v1");
    // hltPathNames.push_back("HLT_PFHT600_v1");
    
    // hltPathNames.push_back("HLT_PFHT900_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_NoL1_v1");
    // hltPathNames.push_back("HLT_PFMET170_NoiseCleaned_NoL1_v1");
    // hltPathNames.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT900_v1");
    // hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_v1");
    // hltPathNames.push_back("HLT_PFMET170_NoiseCleaned_v1");
    // hltPathNames.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v1");
    
    // hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_NoL1_v1");
    // hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_4jet_NoL1_v1");
    // hltPathNames.push_back("HLT_Rsq0p36_NoL1_v1");
    // hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_v1");
    // hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_4jet_v1");
    // hltPathNames.push_back("HLT_Rsq0p36_v1");
    
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
    // std::vector<const reco::Candidate*> gctJetAll                      = getCollections( iEvent, srcGctJetAll_);
    std::vector<const reco::Candidate*> genJet4Unskimmed                = getCollections( iEvent, srcGen4Jet_);
    //std::vector<const reco::Candidate*> genJet5Unskimmed                = getCollections( iEvent, srcGen5Jet_);
    std::vector<const reco::Candidate*> hltAk4CaloUnskimmed             = getCollections( iEvent, srcHLTAk4Calo );
    // std::vector<const reco::Candidate*> hltAk4CaloNoFastJetUnskimmed    = getCollections( iEvent, srcHLTAk4CaloNoFastJet);
    std::vector<const reco::Candidate*> hltAk4PFUnskimmed               = getCollections( iEvent, srcHLTAk4PF );
    std::vector<const reco::Candidate*> recoAk4CaloUnskimmed            = getCollections( iEvent, srcAk4Calo );
    std::vector<const reco::Candidate*> recoAk4PFUnskimmed              = getCollections( iEvent, srcAk4PF   );


    std::vector<const reco::Candidate*> gctForUnskimmed                 = getCollections( iEvent, srcGctJetForward_);
    std::vector<const reco::Candidate*> genJet4ForUnskimmed             = getCollections( iEvent, srcGen4Jet_);
    //std::vector<const reco::Candidate*> genJet5ForUnskimmed             = getCollections( iEvent, srcGen5Jet_);
    std::vector<const reco::Candidate*> hltAk4CaloForUnskimmed          = getCollections( iEvent, srcHLTAk4Calo );
    //std::vector<const reco::Candidate*> hltAk4CaloNoFastJetForUnskimmed = getCollections( iEvent, srcHLTAk4CaloNoFastJet);
    std::vector<const reco::Candidate*> hltAk4PFForUnskimmed            = getCollections( iEvent, srcHLTAk4PF );
    std::vector<const reco::Candidate*> recoAk4CaloForUnskimmed         = getCollections( iEvent, srcAk4Calo );
    std::vector<const reco::Candidate*> recoAk4PFForUnskimmed           = getCollections( iEvent, srcAk4PF   );
 
 

    // Skim jet collections
    // ----------------------------------------

    // Central jets
    // --------------------
    std::vector<const reco::Candidate*> gctCen              = skimJets(gctCenUnskimmed,              minPt, minEtaCen, maxEtaCen );
    //    std::vector<const reco::Candidate*> uctCen  = skimJets(uctCenUnskimmed,  minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> genAk4              = skimJets(genJet4Unskimmed,             minPt, minEtaCen, maxEtaCen );
    //std::vector<const reco::Candidate*> genAk5              = skimJets(genJet5Unskimmed,             minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltAk4Calo          = skimJets(hltAk4CaloUnskimmed,          minPt, minEtaCen, maxEtaCen );
    //std::vector<const reco::Candidate*> hltAk4CaloNoFastJet = skimJets(hltAk4CaloNoFastJetUnskimmed, minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltAk4PF            = skimJets(hltAk4PFUnskimmed,            minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> recoAk4Calo         = skimJets(recoAk4CaloUnskimmed,         minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> recoAk4PF           = skimJets(recoAk4PFUnskimmed,           minPt, minEtaCen, maxEtaCen );


    // Forward jets
    // --------------------
    std::vector<const reco::Candidate*> gctFor                 = skimJets(gctForUnskimmed,                 minPt, minEtaFor, maxEtaFor );
    //    std::vector<const reco::Candidate*> uctCenFor  = skimJets(uctForUnskimmed,  minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> genAk4For              = skimJets(genJet4ForUnskimmed,             minPt, minEtaFor, maxEtaFor );
    //std::vector<const reco::Candidate*> genAk5For              = skimJets(genJet5ForUnskimmed,             minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltAk4CaloFor          = skimJets(hltAk4CaloForUnskimmed,          minPt, minEtaFor, maxEtaFor );
    //std::vector<const reco::Candidate*> hltAk4CaloNoFastJetFor = skimJets(hltAk4CaloNoFastJetForUnskimmed, minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltAk4PFFor            = skimJets(hltAk4PFForUnskimmed,            minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> recoAk4CaloFor         = skimJets(recoAk4CaloForUnskimmed,         minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> recoAk4PFFor           = skimJets(recoAk4PFForUnskimmed,           minPt, minEtaFor, maxEtaFor );


    // Clear previous event's objects
    // --------------------------------------------------------------------------------

    genElectronPt->clear();
    genElectronEta->clear();
    genElectronPhi->clear();
    genMuonPt ->clear();
    genMuonEta->clear();
    genMuonPhi->clear();

    genAk4ForMaxPt    = 0;
    recoAk4PFForMaxPt = 0;
    hltAk4PFForMaxPt  = 0;
    if ( genAk4For.size()   > 0){    genAk4ForMaxPt =    genAk4For.at(0)->pt(); }
    if (recoAk4PFFor.size() > 0){ recoAk4PFForMaxPt = recoAk4PFFor.at(0)->pt(); }
    if ( hltAk4PFFor.size() > 0){  hltAk4PFForMaxPt =  hltAk4PFFor.at(0)->pt(); }

    
    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
	jetPt[*iLvl] ->clear(); 
	jetPx[*iLvl] ->clear(); 
	jetPy[*iLvl] ->clear(); 
	jetPhi[*iLvl]->clear();
	jetEta[*iLvl]->clear();
    }


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
    //hltAk4CaloNoFastJetAlphaTHT40 = calculateAlphaTHT( hltAk4CaloNoFastJet, 40.);
    hltAk4PFAlphaTHT40            = calculateAlphaTHT( hltAk4PF,    40.);
    recoAk4CaloAlphaTHT40         = calculateAlphaTHT( recoAk4Calo, 40.);
    recoAk4PFAlphaTHT40           = calculateAlphaTHT( recoAk4PF,   40.);

    genAk4AlphaTHT50              = calculateAlphaTHT( genAk4,      50.);
    hltAk4CaloAlphaTHT50          = calculateAlphaTHT( hltAk4Calo,  50.);
    //    hltAk4CaloNoFastJetAlphaTHT50 = calculateAlphaTHT( hltAk4CaloNoFastJet, 50.);
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


    // Store jet collections
    // ----------------------------------------

    storeJet( "gctCen",                 gctCen,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk4",                 genAk4,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    //storeJet( "genAk5",                 genAk5,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4Calo",             hltAk4Calo,             jetPt, jetPx, jetPy, jetEta, jetPhi );
    //storeJet( "hltAk4CaloNoFastJet",    hltAk4CaloNoFastJet,    jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4PF",               hltAk4PF,               jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4Calo",            recoAk4Calo,            jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4PF",              recoAk4PF,              jetPt, jetPx, jetPy, jetEta, jetPhi );

    storeJet( "gctFor",                 gctFor,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk4For",              genAk4For,              jetPt, jetPx, jetPy, jetEta, jetPhi );
    //storeJet( "genAk5For",              genAk5For,              jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4CaloFor",          hltAk4CaloFor,          jetPt, jetPx, jetPy, jetEta, jetPhi );
    //storeJet( "hltAk4CaloNoFastJetFor", hltAk4CaloNoFastJetFor, jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4PFFor",            hltAk4PFFor,            jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4CaloFor",         recoAk4CaloFor,         jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "recoAk4PFFor",           recoAk4PFFor,           jetPt, jetPx, jetPy, jetEta, jetPhi );

    
    
    // Perform jet matching
    // ----------------------------------------

    // genMatchedAk4HLTPF  = matchJetCollections( jetPt["genAk4"],   jetEta["genAk4"],   jetPhi["genAk4"],
    // 					       jetPt["hltAk4PF"], jetEta["hltAk4PF"], jetPhi["hltAk4PF"], 20., 0.25 );
    
    // genMatchedAk4RECOPF = matchJetCollections( jetPt["genAk4"],    jetEta["genAk4"],    jetPhi["genAk4"],
    // 					       jetPt["recoAk4PF"], jetEta["recoAk4PF"], jetPhi["recoAk4PF"], 20., 0.25 );
    

    // Gen particles
    // --------------------------------------------------------------------------------

    genLeptonVeto      = false;
    genElectronVeto    = false;
    genIsoLeptonVeto   = false;
    genIsoElectronVeto = false;

    if (makeGenParticles){
      edm::Handle< std::vector<reco::GenParticle> > genParticles;
      iEvent.getByLabel(srcGenParticles_, genParticles);
      for(std::vector<reco::GenParticle>::const_iterator iter = genParticles->begin(); iter != genParticles->end(); ++iter){
    	const reco::GenParticle& genLep = *iter;

    	double genLepPt  = genLep.p4().pt();
    	double genLepEta = genLep.p4().eta();
    	double genLepPhi = genLep.p4().phi();

    	// Electron
    	if(TMath::Abs(genLep.pdgId()) == 11){
    	  if ( (genLepPt >= genElectronMinPt) && (TMath::Abs(genLepEta) <= genElectronMaxEta) ){
    	    genElectronPt ->push_back( genLepPt  );
    	    genElectronEta->push_back( genLepEta );
    	    genElectronPhi->push_back( genLepPhi );
    	    genLeptonVeto   = true;
    	    genElectronVeto = true;
    	  }
    	} // End lepton requirement
    	// Muon
    	if(TMath::Abs(genLep.pdgId()) == 13){
    	  if ( (genLepPt >= genMuonMinPt) && (TMath::Abs(genLepEta) <= genMuonMaxEta) ){
    	    genMuonPt ->push_back( genLepPt  );
   	    genMuonEta->push_back( genLepEta );
    	    genMuonPhi->push_back( genLepPhi );
    	    genLeptonVeto = true;
    	  }
    	} // End Muon requirement

      }
    } // End gen particles

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


 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakeTrees);
