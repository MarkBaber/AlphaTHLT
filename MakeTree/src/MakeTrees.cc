// **************************************************
// Switches
// **************************************************
//#define DATA
#define SIMULATION
//#define L1
//#define RECO
// Remove isolated leptons from gen and HLT jets
//#define LEPTON_XCLEANING
// **************************************************
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
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
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/METReco/interface/GenMET.h"
#include "DataFormats/METReco/interface/GenMETCollection.h"
#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/CaloMETCollection.h"

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

TLorentzVector P4toTLV (reco::Particle::LorentzVector a){ return TLorentzVector( a.px(), a.py(), a.pz(), a.energy() ); }


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

    std::vector< const reco::Candidate*> skimJets(const std::vector<const reco::Candidate*>& inJets, double minPt, double minEta, double maxEta);

    float leadL1GenDeltaR( const std::vector<const reco::Candidate*>& gctCen, const std::vector<const reco::Candidate*>& gctFor, 
			   const std::vector<const reco::Candidate*>& gen);
    float leadHLTGenDeltaR( const std::vector<const reco::Candidate*>& hlt, const std::vector<const reco::Candidate*>& gen);

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

    edm::InputTag srcGctMht_;
    edm::InputTag srcGctMet_;
    VInputTag srcGctJetCentral_;
    VInputTag srcGctJetForward_;
    VInputTag srcGctJetAll_;

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


    std::pair<float,float> genAlphaTHT40;
    std::pair<float,float> hltPFAlphaTHT40;
    std::pair<float,float> hltCaloAlphaTHT40;
    std::pair<float,float> hltCaloIDAlphaTHT40;
  //    std::pair<float,float> hltCaloNoFastJetAlphaTHT40;
  // std::pair<float,float> recoPFAlphaTHT40;
  // std::pair<float,float> recoCaloAlphaTHT40;


  float gen_AlphaTPrime40;
  float hltCalo_AlphaTPrime40;
  float hltCaloID_AlphaTPrime40;
  float hltPF_AlphaTPrime40;  
  // float recoCalo_AlphaTPrime40;
  // float recoPF_AlphaTPrime40; 

    // Dynamic HT, AlphaT
    float dynamicJetThreshold;
    std::vector<std::pair<float,float> > genDynamicAlphaTHT40;
    std::vector<std::pair<float,float> > hltPFDynamicAlphaTHT40;
    std::vector<float> genDynamicAlphaT40;
    std::vector<float> genDynamicHT40;
    std::vector<float> hltPFDynamicAlphaT40;
    std::vector<float> hltPFDynamicHT40;

  //    std::vector<std::pair<float,float> > recoPFDynamicAlphaTHT40;


    std::pair<float,float> genMHT40;
    std::pair<float,float> hltPFMHT40;
    std::pair<float,float> hltCaloMHT40;
    std::pair<float,float> hltCaloIDMHT40;
    // std::pair<float,float> recoPFMHT40;
    // std::pair<float,float> recoCaloMHT40;

    std::pair<float,float> genForMHT40;
    std::pair<float,float> hltCaloForMHT40;
    std::pair<float,float> hltPFForMHT40;
  //std::pair<float,float> recoPFForMHT40;

  // DeltaR between leading jets of collections
  float L1GenDeltaR;
  float HLTGenDeltaR;
  UInt_t hpuVeto;

  float hltMetCaloPFMht40_DeltaPhi;
  float genMetCaloMht40_DeltaPhi;

  UInt_t genNJet40;
  UInt_t hltPFNJet40;
  UInt_t hltCaloNJet40;
  UInt_t hltCaloIDNJet40;
  // UInt_t recoPFNJet40;
  // UInt_t recoCaloNJet40;
  
  Int_t genNJetBin40;
  Int_t hltPFNJetBin40;
  Int_t hltCaloNJetBin40;
  Int_t hltCaloIDNJetBin40;
  // Int_t recoPFNJetBin40;
  // Int_t recoCaloNJetBin40;

  Int_t genHTBin40;
  Int_t hltPFHTBin40;
  Int_t hltCaloHTBin40;
  Int_t hltCaloIDHTBin40;
  // Int_t recoPFHTBin40;
  // Int_t recoCaloHTBin40;


    UInt_t maxjet_;
    bool usePU_; 
    unsigned int run;
    unsigned int lumi;
    unsigned int event;

    // Jet skim cuts
    double minPt;
    double minEtaCen;
    double maxEtaCen;
    double minEtaFor;
    double maxEtaFor;


    edm::InputTag HLTResultsTag;
  
    float PThat;

  float genForMaxPt;
  //  float recoPFForMaxPt;
  float hltPFForMaxPt;

  // Lead jet
  float genLeadJetPt;
  //  float recoPFLeadJetPt;
  float hltPFLeadJetPt;
  float hltCaloLeadJetPt;
  float hltCaloIDLeadJetPt;
  // Second jet
  float genSecondJetPt;
  //  float recoPFSecondJetPt;
  float hltPFSecondJetPt;
  float hltCaloSecondJetPt;
  float hltCaloIDSecondJetPt;
  // Dijet avg
  float genDijetAvePt;
  //  float recoPFDijetAvePt;
  float hltPFDijetAvePt;
  float hltCaloDijetAvePt;
  float hltCaloIDDijetAvePt;


  // Biased deltaPhi
  float genBiasedDPhi;   
  //  float recoPFBiasedDPhi; 
  float hltPFBiasedDPhi;  
  // Tag jet index
  Int_t genBiasedDPhiIndex;   
  //  Int_t recoPFBiasedDPhiIndex; 
  Int_t hltPFBiasedDPhiIndex;  


  // Matched indices
  std::vector<int> genMatchedHLTPF;

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

  edm::InputTag hltCaloMetTag_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBits_;
  edm::EDGetTokenT<edm::TriggerResults> triggerBitsData_;

  // Jets
  edm::EDGetTokenT<reco::CaloJetCollection> hltCaloJetToken_;
  edm::EDGetTokenT<reco::CaloJetCollection> hltCaloJetIDToken_;
  edm::EDGetTokenT<reco::PFJetCollection>   hltPFJetToken_;
  edm::EDGetTokenT<reco::GenJetCollection>  genJetToken_;

  // Met
  edm::EDGetTokenT<reco::CaloMETCollection> hltCaloMetToken_;
  edm::EDGetTokenT<reco::METCollection>     hltPFMetToken_;
  edm::EDGetTokenT<reco::GenMETCollection>  genMetCaloToken_;
  edm::EDGetTokenT<reco::GenMETCollection>  genMetTrueToken_;


};

MakeTrees::MakeTrees(const edm::ParameterSet& pset){


    // Initialize the ntuple builder
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("tree", "tree");


#ifdef SIMULATION
    lvl_.push_back("gen");
    lvl_.push_back("genFor");
#endif
    lvl_.push_back("hltPF");
    lvl_.push_back("hltPFFor");
    lvl_.push_back("hltCalo");
    lvl_.push_back("hltCaloFor");
    lvl_.push_back("hltCaloID");
    lvl_.push_back("hltCaloIDFor");

    // lvl_.push_back("recoPF");
    // lvl_.push_back("recoPFFor");
    // lvl_.push_back("recoCalo");
    // lvl_.push_back("recoCaloFor");
#ifdef L1
    lvl_.push_back("gctCen");
    lvl_.push_back("gctFor");
#endif

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
    genForMaxPt    = 0;
    //    recoPFForMaxPt = 0;
    hltPFForMaxPt  = 0;

    // Lead jet 
    genLeadJetPt    = 0;
    //    recoPFLeadJetPt = 0;
    hltPFLeadJetPt  = 0;
    hltCaloLeadJetPt= 0;
    hltCaloIDLeadJetPt= 0;
    // Second jet 
    genSecondJetPt    = 0;
    //    recoPFSecondJetPt = 0;
    hltPFSecondJetPt  = 0;
    hltCaloSecondJetPt= 0;
    hltCaloIDSecondJetPt= 0;
    // Dijet avg 
    genDijetAvePt    = 0;
    //    recoPFDijetAvePt = 0;
    hltPFDijetAvePt  = 0;
    hltCaloDijetAvePt= 0;
    hltCaloIDDijetAvePt= 0;


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

    tree->Branch("run",  &run,  "run/i");
    tree->Branch("lumi", &lumi, "lumi/i");
    tree->Branch("evt",  &event,"evt/i");

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




#ifdef SIMULATION
    tree->Branch("genFor_MaxPt",      &genForMaxPt,               "genFor_MaxPt/f");
    tree->Branch("genLead_Pt",      &genLeadJetPt,               "genLead_Pt/f");
    tree->Branch("genDijetAve_Pt",      &genDijetAvePt,               "genDijetAve_Pt/f");
    tree->Branch("genSecond_Pt",      &genSecondJetPt,               "genSecond_Pt/f");
#endif

    // Maximum forward jet PT
    //tree->Branch("recoPFFor_MaxPt",   &recoPFForMaxPt,            "recoPFFor_MaxPt/f");
    tree->Branch("hltPFFor_MaxPt",    &hltPFForMaxPt,             "hltPFFor_MaxPt/f");


    // Lead jet 
    //    tree->Branch("recoPFLead_Pt",   &recoPFLeadJetPt,            "recoPFLead_Pt/f");
    tree->Branch("hltPFLead_Pt",    &hltPFLeadJetPt,             "hltPFLead_Pt/f");
    tree->Branch("hltCaloLead_Pt",  &hltCaloLeadJetPt,           "hltCaloLead_Pt/f");
    tree->Branch("hltCaloIDLead_Pt",&hltCaloIDLeadJetPt,         "hltCaloIDLead_Pt/f");
    // Second jet 
    //tree->Branch("recoPFSecond_Pt",   &recoPFSecondJetPt,            "recoPFSecond_Pt/f");
    tree->Branch("hltPFSecond_Pt",    &hltPFSecondJetPt,             "hltPFSecond_Pt/f");
    tree->Branch("hltCaloSecond_Pt",  &hltCaloSecondJetPt,           "hltCaloSecond_Pt/f");
    tree->Branch("hltCaloIDSecond_Pt",&hltCaloIDSecondJetPt,         "hltCaloIDSecond_Pt/f");
    // Dijet avg 

    //tree->Branch("recoPFDijetAve_Pt",   &recoPFDijetAvePt,            "recoPFDijetAve_Pt/f");
    tree->Branch("hltPFDijetAve_Pt",    &hltPFDijetAvePt,             "hltPFDijetAve_Pt/f");
    tree->Branch("hltCaloDijetAve_Pt",  &hltCaloDijetAvePt,           "hltCaloDijetAve_Pt/f");
    tree->Branch("hltCaloIDDijetAve_Pt",&hltCaloIDDijetAvePt,         "hltCaloIDDijetAve_Pt/f");

#ifdef SIMULATION 
    // Used for rate studies only
    // tree->Branch("hpuVeto",                  &hpuVeto,                      "hpuVeto/b");
    // tree->Branch("leadL1GenDeltaR",          &L1GenDeltaR,                  "leadL1GenDeltaR/f");
    // tree->Branch("leadHLTGenDeltaR",         &HLTGenDeltaR,                 "leadHLTGenDeltaR/f");
#endif

    // Biased deltaPhi
    tree->Branch("gen_BiasedDPhi",         &genBiasedDPhi,         "gen_BiasedDPhi/f");
    //    tree->Branch("recoPF_BiasedDPhi",      &recoPFBiasedDPhi,      "recoPF_BiasedDPhi/f");
    tree->Branch("hltPF_BiasedDPhi",       &hltPFBiasedDPhi,       "hltPF_BiasedDPhi/f");
    tree->Branch("gen_BiasedDPhiIndex",    &genBiasedDPhiIndex,    "gen_BiasedDPhiIndex/I");
    //tree->Branch("recoPF_BiasedDPhiIndex", &recoPFBiasedDPhiIndex, "recoPF_BiasedDPhiIndex/I");
    tree->Branch("hltPF_BiasedDPhiIndex",  &hltPFBiasedDPhiIndex,  "hltPF_BiasedDPhiIndex/I");


#ifdef SIMULATION
    tree->Branch("gen_AlphaT40",      &genAlphaTHT40.first,       "gen_AlphaT40/f");
    tree->Branch("gen_HT40",          &genAlphaTHT40.second,      "gen_HT40/f");
    tree->Branch("gen_AlphaTPrime40",      &gen_AlphaTPrime40,      "gen_AlphaTPrime40/f");
    tree->Branch("gen_DynamicAlphaTHT40",    "std::vector<std::pair<float,float>>", &genDynamicAlphaTHT40);
    tree->Branch("gen_DynamicAlphaT40",   "std::vector<float>", &genDynamicAlphaT40);
    tree->Branch("gen_DynamicHT40",       "std::vector<float>", &genDynamicHT40);

    tree->Branch("gen_MhtPT40",       &genMHT40.first,       "gen_MhtPT40/f");
    tree->Branch("gen_MhtPhi40",      &genMHT40.second,      "gen_MhtPhi40/f");
    tree->Branch("genFor_MhtPT40",       &genForMHT40.first,       "genFor_MhtPT40/f");
    tree->Branch("genFor_MhtPhi40",      &genForMHT40.second,      "genFor_MhtPhi40/f");
#endif

    // AlphaT, HT  
    tree->Branch("hltPF_AlphaT40",    &hltPFAlphaTHT40.first,     "hltPF_AlphaT40/f");
    tree->Branch("hltPF_HT40",        &hltPFAlphaTHT40.second,    "hltPF_HT40/f");
    tree->Branch("hltCalo_AlphaT40",  &hltCaloAlphaTHT40.first,   "hltCalo_AlphaT40/f");
    tree->Branch("hltCalo_HT40",      &hltCaloAlphaTHT40.second,  "hltCalo_HT40/f");
    tree->Branch("hltCaloID_AlphaT40",  &hltCaloIDAlphaTHT40.first,   "hltCaloID_AlphaT40/f");
    tree->Branch("hltCaloID_HT40",      &hltCaloIDAlphaTHT40.second,  "hltCaloID_HT40/f");

    // tree->Branch("recoPF_AlphaT40",   &recoPFAlphaTHT40.first,    "recoPF_AlphaT40/f");
    // tree->Branch("recoPF_HT40",       &recoPFAlphaTHT40.second,   "recoPF_HT40/f");
    // tree->Branch("recoCalo_AlphaT40", &recoCaloAlphaTHT40.first,  "recoCalo_AlphaT40/f");
    // tree->Branch("recoCalo_HT40",     &recoCaloAlphaTHT40.second, "recoCalo_HT40/f");

    tree->Branch("hltCalo_AlphaTPrime40",  &hltCalo_AlphaTPrime40,  "hltCalo_AlphaTPrime40/f");
    tree->Branch("hltCaloID_AlphaTPrime40",  &hltCaloID_AlphaTPrime40,  "hltCaloID_AlphaTPrime40/f");
    tree->Branch("hltPF_AlphaTPrime40",    &hltPF_AlphaTPrime40,    "hltPF_AlphaTPrime40/f");
    // tree->Branch("recoCalo_AlphaTPrime40", &recoCalo_AlphaTPrime40, "recoCalo_AlphaTPrime40/f");
    // tree->Branch("recoPF_AlphaTPrime40",   &recoPF_AlphaTPrime40,   "recoPF_AlphaTPrime40/f");


    // Dynamic AlphaT, HT
    tree->Branch("hltPF_DynamicAlphaTHT40",  "std::vector<std::pair<float,float>>", &hltPFDynamicAlphaTHT40);
    tree->Branch("hltPF_DynamicAlphaT40", "std::vector<float>", &hltPFDynamicAlphaT40);
    tree->Branch("hltPF_DynamicHT40",     "std::vector<float>", &hltPFDynamicHT40);

    //tree->Branch("recoPF_DynamicAlphaTHT40", "std::vector<std::pair<float,float>>", &recoPFDynamicAlphaTHT40);


    // MHT
    tree->Branch("hltPF_MhtPT40",     &hltPFMHT40.first,     "hltPF_MhtPT40/f");
    tree->Branch("hltPF_MhtPhi40",    &hltPFMHT40.second,    "hltPF_MhtPhi40/f");
    tree->Branch("hltCalo_MhtPT40",   &hltCaloMHT40.first,   "hltCalo_MhtPT40/f");
    tree->Branch("hltCalo_MhtPhi40",  &hltCaloMHT40.second,  "hltCalo_MhtPhi40/f");
    tree->Branch("hltCaloID_MhtPT40",   &hltCaloIDMHT40.first,   "hltCaloID_MhtPT40/f");
    tree->Branch("hltCaloID_MhtPhi40",  &hltCaloIDMHT40.second,  "hltCaloID_MhtPhi40/f");

    // tree->Branch("recoPF_MhtPT40",    &recoPFMHT40.first,    "recoPF_MhtPT40/f");
    // tree->Branch("recoPF_MhtPhi40",   &recoPFMHT40.second,   "recoPF_MhtPhi40/f");
    // tree->Branch("recoCalo_MhtPT40",  &recoCaloMHT40.first,  "recoCalo_MhtPT40/f");
    // tree->Branch("recoCalo_MhtPhi40", &recoCaloMHT40.second, "recoCalo_MhtPhi40/f");

    tree->Branch("hltCaloFor_MhtPT40",   &hltCaloForMHT40.first,   "hltCaloFor_MhtPT40/f");
    tree->Branch("hltCaloFor_MhtPhi40",  &hltCaloForMHT40.second,  "hltCaloFor_MhtPhi40/f");
    tree->Branch("hltPFFor_MhtPT40",     &hltPFForMHT40.first,     "hltPFFor_MhtPT40/f");
    tree->Branch("hltPFFor_MhtPhi40",    &hltPFForMHT40.second,    "hltPFFor_MhtPhi40/f");
    // tree->Branch("recoPFFor_MhtPT40",    &recoPFForMHT40.first,    "recoPFFor_MhtPT40/f");
    // tree->Branch("recoPFFor_MhtPhi40",   &recoPFForMHT40.second,   "recoPFFor_MhtPhi40/f");


    tree->Branch("hltMetCaloPFMht40_DeltaPhi", &hltMetCaloPFMht40_DeltaPhi, "hltMetCaloPFMht40_DeltaPhi/f");
    //tree->Branch("genMetCaloMht40_DeltaPhi",   &genMetCaloMht40_DeltaPhi,   "genMetCaloMht40_DeltaPhi/f");

    //tree->Branch("gen_NJet40",      &genNJet40,       "gen_NJet40/i");
    tree->Branch("hltPF_NJet40",    &hltPFNJet40,     "hltPF_NJet40/i");
    tree->Branch("hltCalo_NJet40",  &hltCaloNJet40,   "hltCalo_NJet40/i");
    tree->Branch("hltCaloID_NJet40",  &hltCaloIDNJet40,   "hltCaloID_NJet40/i");
    // tree->Branch("recoPF_NJet40",   &recoPFNJet40,    "recoPF_NJet40/i");
    // tree->Branch("recoCalo_NJet40", &recoCaloNJet40,  "recoCalo_NJet40/i");

    //tree->Branch("gen_NJetBin40",      &genNJetBin40,       "gen_NJetBin40/I");
    tree->Branch("hltPF_NJetBin40",    &hltPFNJetBin40,     "hltPF_NJetBin40/I");
    tree->Branch("hltCalo_NJetBin40",  &hltCaloNJetBin40,   "hltCalo_NJetBin40/I");
    tree->Branch("hltCaloID_NJetBin40",  &hltCaloIDNJetBin40,   "hltCaloID_NJetBin40/I");
    // tree->Branch("recoPF_NJetBin40",   &recoPFNJetBin40,    "recoPF_NJetBin40/I");
    // tree->Branch("recoCalo_NJetBin40", &recoCaloNJetBin40,  "recoCalo_NJetBin40/I");

    //tree->Branch("gen_HTBin40",      &genHTBin40,       "gen_HTBin40/I");
    tree->Branch("hltPF_HTBin40",    &hltPFHTBin40,     "hltPF_HTBin40/I");
    tree->Branch("hltCalo_HTBin40",  &hltCaloHTBin40,   "hltCalo_HTBin40/I");
    tree->Branch("hltCaloID_HTBin40",  &hltCaloIDHTBin40,   "hltCaloID_HTBin40/I");
    // tree->Branch("recoPF_HTBin40",   &recoPFHTBin40,    "recoPF_HTBin40/I");
    // tree->Branch("recoCalo_HTBin40", &recoCaloHTBin40,  "recoCalo_HTBin40/I");

#ifdef L1
    // Energy sums
    tree->Branch("gct_Ht",       &ht_["gct"],      "gct_Ht/f");
    tree->Branch("gct_MhtPt",    &mhtPt_["gct"],   "gct_MhtPt/f");
    tree->Branch("gct_MhtPhi",   &mhtPhi_["gct"],  "gct_MhtPhi/f");
    tree->Branch("gct_MhtDivHt", &mhtDivHt_["gct"],"gct_MhtDivHt/f"); 

    tree->Branch("gct_Et",     &et_["gct"],     "gct_Et/f");
    tree->Branch("gct_MetPt",  &metPt_["gct"],  "gct_MetPt/f");
    tree->Branch("gct_MetPhi", &metPhi_["gct"], "gct_MetPhi/f");
#endif

#ifdef SIMULATION
    // GEN MET
    tree->Branch("genMetCalo_MetPt",              &metPt_["genMetCalo"],              "genMetCalo_MetPt/f");
    tree->Branch("genMetCaloAndNonPrompt_MetPt",  &metPt_["genMetCaloAndNonPrompt"],  "genMetCaloAndNonPrompt_MetPt/f");
    tree->Branch("genMetTrue_MetPt",              &metPt_["genMetTrue"],              "genMetTrue_MetPt/f");
    tree->Branch("genMetCalo_MetPhi",             &metPhi_["genMetCalo"],             "genMetCalo_MetPhi/f");
    tree->Branch("genMetCaloAndNonPrompt_MetPhi", &metPhi_["genMetCaloAndNonPrompt"], "genMetCaloAndNonPrompt_MetPhi/f");
    tree->Branch("genMetTrue_MetPhi",             &metPhi_["genMetTrue"],             "genMetTrue_MetPhi/f");
#endif

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

#ifdef SIMULATION
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

    tree->Branch("genMatchedHLTPF", "std::vector<int>", &genMatchedHLTPF);

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
#endif

    tree->Branch("nIsoMuons",     &nIsoMuons,     "nIsoMuons/i");
    tree->Branch("nIsoElectrons", &nIsoElectrons, "nIsoMuons/i");

    // Store L1 seeds
    // ------------------------------------------------------------
#ifdef L1
    tree->Branch("L1HTT175",        &L1HTT175,         "L1HTT175/b");
    tree->Branch("L1ETM70",         &L1ETM70,          "L1ETM70/b");
    tree->Branch("L1HTT175OrETM70", &L1HTT175OrETM70,  "L1HTT175OrETM70/b");
    tree->Branch("L1Jet_DPhi",      &L1Jet_DPhi,       "L1Jet_DPhi/f");
#endif

    // Store HLT paths
    // ------------------------------------------------------------

    hltPathNames.push_back("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p57_v2"); 
    hltPathNames.push_back("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p55_v2"); 
    hltPathNames.push_back("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p53_v2"); 
    hltPathNames.push_back("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p52_v2"); 
    hltPathNames.push_back("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p51_v2"); 

     hltPathNames.push_back("HLT_PFHT200_DiPFJetAve90_PFAlphaT0p63_v2"); 
     hltPathNames.push_back("HLT_PFHT250_DiPFJetAve90_PFAlphaT0p58_v2");
     hltPathNames.push_back("HLT_PFHT300_DiPFJetAve90_PFAlphaT0p54_v2"); 
     hltPathNames.push_back("HLT_PFHT350_DiPFJetAve90_PFAlphaT0p53_v2"); 
     hltPathNames.push_back("HLT_PFHT400_DiPFJetAve90_PFAlphaT0p52_v2");
                                 
     hltPathNames.push_back("HLT_PFHT200_PFAlphaT0p51_v2");
   
     hltPathNames.push_back("HLT_PFHT200_v2");
     hltPathNames.push_back("HLT_PFHT250_v2");
     hltPathNames.push_back("HLT_PFHT300_v2");
     hltPathNames.push_back("HLT_PFHT350_v3");
     hltPathNames.push_back("HLT_PFHT400_v2");
     hltPathNames.push_back("HLT_PFHT475_v2"); 
     hltPathNames.push_back("HLT_PFHT600_v3");
     hltPathNames.push_back("HLT_PFHT650_v3");
     hltPathNames.push_back("HLT_PFHT800_v2");

     hltPathNames.push_back("HLT_PFMET90_PFMHT90_IDTight_v2");
     hltPathNames.push_back("HLT_PFMET100_PFMHT100_IDTight_v2");
     hltPathNames.push_back("HLT_PFMET110_PFMHT110_IDTight_v2");
     hltPathNames.push_back("HLT_PFMET120_PFMHT120_IDTight_v2");

     hltPathNames.push_back("HLT_IsoMu20_v3");
    
    // Trigger bits
    triggerBits_     = consumes<edm::TriggerResults> (pset.getParameter<edm::InputTag>("HLTResults"));
    triggerBitsData_ = consumes<edm::TriggerResults> (pset.getParameter<edm::InputTag>("HLTResultsData"));
    for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
      TString path = hltPathNames[ iPath ];
      hltPathFired[ path ] = false;
      tree->Branch( path, &hltPathFired[ path ], path + "/b" );
#ifdef DATA
      TString dataPath = "Data_"+path;
      hltPathFired[ dataPath ] = false;
      tree->Branch( dataPath, &hltPathFired[ dataPath ], dataPath + "/b" );
#endif
    }



    // Gen 
    srcGenMetCalo_                = pset.getParameter<edm::InputTag>("srcGenMetCalo");
    srcGenMetCaloAndNonPrompt_    = pset.getParameter<edm::InputTag>("srcGenMetCaloAndNonPrompt");
    srcGenMetTrue_                = pset.getParameter<edm::InputTag>("srcGenMetTrue");
    makeGenParticles           = pset.getParameter<bool>("MakeGenParticles");
    srcGenParticles_           = pset.getParameter<edm::InputTag>("srcGenParticles");
    genElectronMinPt           = pset.getParameter<double>("genElectronMinPt");
    genElectronMaxEta          = pset.getParameter<double>("genElectronMaxEta");
    genMuonMinPt               = pset.getParameter<double>("genMuonMinPt");
    genMuonMaxEta              = pset.getParameter<double>("genMuonMaxEta");
    genPhotonMinPt             = pset.getParameter<double>("genPhotonMinPt");
    genPhotonMaxEta            = pset.getParameter<double>("genPhotonMaxEta");

#ifdef L1
    srcUctMET_ = pset.getParameter<edm::InputTag>("srcUctMet");
    srcUctMht_ = pset.getParameter<edm::InputTag>("srcUctMht");
    srcUctJet_ = pset.getParameter<edm::InputTag>("srcUctJet");

    srcGctMht_        = pset.getParameter<edm::InputTag>("srcGctMht");
    srcGctMet_        = pset.getParameter<edm::InputTag>("srcGctMet");
    srcGctJetCentral_ = pset.getParameter<VInputTag>("srcGctJetCentral");
    srcGctJetForward_ = pset.getParameter<VInputTag>("srcGctJetForward");
#endif


    hltCaloJetToken_   = consumes<reco::CaloJetCollection>(pset.getParameter<edm::InputTag>("hltCaloSrc"));
    hltCaloJetIDToken_ = consumes<reco::CaloJetCollection>(pset.getParameter<edm::InputTag>("hltCaloIDSrc"));
    hltPFJetToken_     = consumes<reco::PFJetCollection>  (pset.getParameter<edm::InputTag>("hltPFSrc"));
    genJetToken_       = consumes<reco::GenJetCollection> (pset.getParameter<edm::InputTag>("genSrc"));

    hltCaloMetToken_   = consumes<reco::CaloMETCollection>(pset.getParameter<edm::InputTag>("hltCaloMetSrc"));
    hltPFMetToken_     = consumes<reco::METCollection>    (pset.getParameter<edm::InputTag>("hltPFMetSrc"));
    genMetCaloToken_   = consumes<reco::GenMETCollection> (pset.getParameter<edm::InputTag>("genMetCaloSrc")); 
    genMetTrueToken_   = consumes<reco::GenMETCollection> (pset.getParameter<edm::InputTag>("genMetTrueSrc")); 


#ifdef RECO    
    srcCalo       = pset.getParameter<VInputTag>("srcCalo");
    srcPF         = pset.getParameter<VInputTag>("srcPF");
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
		const { return candA->pt() > candB->pt(); }
    };


  // Turn a set of InputTags into a collection of candidate pointers.
  //edm::Handle<std::vector<reco::CaloJet> >&
  std::vector<const reco::Candidate*> getJetCollection(const edm::Event& evt, const edm::Handle<std::vector<reco::CaloJet> >& handle, float ptMin=20) {
    //std::vector<const reco::Candidate*> getJetCollection(const edm::Event& evt, const edm::EDGetTokenT<reco::CandidateCollection>& token, float ptMin=20){
    std::vector<const reco::Candidate*> output;
    // edm::Handle<reco::CandidateCollection> handle;
    // evt.getByToken(token, handle);

    if (handle.isValid()){
      for (size_t j = 0; j < handle->size(); ++j) {
	const reco::Candidate& object = handle->at(j);
	output.push_back(&object);
      }
    }
    return output;
  }
  
  std::vector<const reco::Candidate*> getJetCollection(const edm::Event& evt, const edm::Handle<std::vector<reco::PFJet> >& handle, float ptMin=20) {
    std::vector<const reco::Candidate*> output;
    for (size_t j = 0; j < handle->size(); ++j) {
      const reco::Candidate& object = handle->at(j);
      output.push_back(&object);
    }
    return output;
  }
  std::vector<const reco::Candidate*> getJetCollection(const edm::Event& evt, const edm::Handle<std::vector<reco::GenJet> >& handle, float ptMin=20){
    std::vector<const reco::Candidate*> output;
    for (size_t j = 0; j < handle->size(); ++j) {
      const reco::Candidate& object = handle->at(j);
      output.push_back(&object);
    }
    return output;
  }

  void getValue(const edm::Event& evt, edm::Handle<std::vector<reco::CaloMET> >& handle, Float_t& pt, Float_t& phi) {
	if(!handle.isValid()){ pt  = 0; phi = 0; }
	else{ pt  = handle->at(0).pt(); phi = handle->at(0).phi(); }
  }
  void getValue(const edm::Event& evt, edm::Handle<std::vector<reco::MET> >& handle, Float_t& pt, Float_t& phi) {
	if(!handle.isValid()){ pt  = 0; phi = 0; }
	else{ pt  = handle->at(0).pt(); phi = handle->at(0).phi(); }
  }
  void getValue(const edm::Event& evt, edm::Handle<std::vector<reco::GenMET> >& handle, Float_t& pt, Float_t& phi) {
	if(!handle.isValid()){ pt  = 0; phi = 0; }
	else{ pt  = handle->at(0).pt(); phi = handle->at(0).phi(); }
  }

  //   void getSumEtL1(const edm::Event& evt, const edm::InputTag& tag, Float_t& sumet,bool upgrade) {
  // 	if(!upgrade) {
  // 	    edm::Handle<l1extra::L1EtMissParticleCollection> handle;
  // 	    evt.getByLabel(tag, handle);
  // 	    sumet = handle->at(0).etTotal();
  // 	} else{
  // 	    edm::Handle<edm::View<reco::Candidate> > handle;
  // 	    evt.getByLabel(tag, handle);
  // 	    sumet = handle->at(0).pt();
  // 	}
  //   }

  // UInt_t calculateNJet(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
    
  //   UInt_t NJets(0);
  //   for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){
  //     if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
  //     NJets++;
  //   }
  //   return NJets;
  // }
  // Int_t calculateNJetBin(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
  //   // Returns the following:
  //   // 0 = No bin
  //   //-1 = eq2a,-2 = eq3a,-3 = ge4a
  //   // 1 = eq2j, 2 = eq3j, 3 = ge4j
  //   Int_t NJetBin(-1); 
  //   for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){
  //     if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
  //     NJetBin++;
  //     if (NJetBin == 3){ break; }
  //   }
  //   if ( NJetBin == -1){ NJetBin = 0; } // No jets

  //   // Asymmetric bin
  //   if ( NJetBin > 0 ){ // Require two jets
  //     if ( jets.at(1)->pt() < 100. ){ NJetBin *= -1; } // Bin is asymmetric
  //   }

  //   return NJetBin;
  // }

  // Int_t calculateHTBin( float ht, float alphaT ){

  //   Int_t HTBin(-1); 
  //   if ( ht >= 200){
  //     if (ht < 250)     { // 200 < HT < 250
  // 	if (alphaT >= 0.65){ HTBin = 0; }
  //     }
  //     else if (ht < 300){ // 250 < HT < 300
  // 	if (alphaT >= 0.60){ HTBin = 1; }
  //     }
  //     else if (ht < 350){ // 300 < HT < 350
  // 	if (alphaT >= 0.55){ HTBin = 2; }
  //     }
  //     else if (ht < 400){ // 350 < HT < 400
  // 	if (alphaT >= 0.53){ HTBin = 3; }
  //     }
  //     else if (ht < 500){ // 400 < HT < 500
  // 	if (alphaT >= 0.52){ HTBin = 4; }
  //     }
  //     else{               // HT > 500
  // 	if (alphaT >= 0.52){ HTBin = 5; }
  //     }

  //   }

  //   return HTBin;
  // }


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

}


void MakeTrees::analyze(const edm::Event& iEvent, const edm::EventSetup& es) {

  // Event information
  run   = (unsigned int) iEvent.id().run();
  lumi  = (unsigned int) iEvent.id().luminosityBlock();
  event = (unsigned int) iEvent.id().event();


  NVTX = 0;
#ifdef SIMULATION 
  // // Get NVTX from simulation
  // edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
  // iEvent.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);

  // std::vector<PileupSummaryInfo>::const_iterator PVI;
  // for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
  //   int BX = PVI->getBunchCrossing();
  //   if(BX == 0) {
  //     NVTX = PVI->getPU_NumInteractions();
  //     break;
  //   }
  // }

  // // Get gen info
  // edm::Handle<GenEventInfoProduct> geninfo;  iEvent.getByLabel("generator",geninfo);
  // std::auto_ptr<bool>   genInfoValid ( new bool( geninfo.isValid() && !geninfo->binningValues().empty()));
  // std::auto_ptr<double> pthat (new double(*genInfoValid ? geninfo->binningValues()[0] : -1.));
  // PThat = *pthat;
#endif



  // ------------------------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------------------------

    edm::Handle<edm::TriggerResults> triggerBits;
    iEvent.getByToken(triggerBits_, triggerBits);
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i) { 
      if ( hltPathFired.find( names.triggerName(i) ) != hltPathFired.end() ){ hltPathFired[ names.triggerName(i) ] = triggerBits->accept(i); } 
    }

#ifdef DATA
    edm::Handle<edm::TriggerResults> triggerBitsData;
    iEvent.getByToken(triggerBitsData_, triggerBitsData);
    const edm::TriggerNames &namesData = iEvent.triggerNames(*triggerBitsData);
    for (unsigned int i = 0, n = triggerBitsData->size(); i < n; ++i) {
      if ( hltPathFired.find( namesData.triggerName(i) ) != hltPathFired.end() ){ hltPathFired[ "Data_"+namesData.triggerName(i) ] = triggerBitsData->accept(i); }
    }
#endif

  // ------------------------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------------------------

    // --------------------------------------------------------------------------------
    // Clear  objects
    // --------------------------------------------------------------------------------
  
    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
	jetPt[*iLvl] ->clear(); 
	jetPx[*iLvl] ->clear(); 
	jetPy[*iLvl] ->clear(); 
	jetPhi[*iLvl]->clear();
	jetEta[*iLvl]->clear();
    }

//     genElectronPt ->clear();
//     genElectronEta->clear();
//     genElectronPhi->clear();
//     genMuonPt ->clear();
//     genMuonEta->clear();
//     genMuonPhi->clear();
//     genPhotonPt ->clear();
//     genPhotonEta->clear();
//     genPhotonPhi->clear();

//     genMuonMatchedGenMuonPt          .clear();
//     genMuonMatchedGenJetPt           .clear();
//     genMuonMatchedHLTPFJetPt         .clear();
//     genMuonMatchedHLTPFJetMuonEF     .clear();
//     genMuonMatchedHLTPFJetElectronEF .clear();

//     genElectronMatchedGenElectronPt      .clear();
//     genElectronMatchedGenJetPt           .clear();
//     genElectronMatchedHLTPFJetPt         .clear();
//     genElectronMatchedHLTPFJetMuonEF     .clear();
//     genElectronMatchedHLTPFJetElectronEF .clear();

//     nIsoElectrons = 0;
//     nIsoMuons     = 0;


    // Input jets without eta or pT requirements
    // --------------------------------------------------------------------------------
#ifdef L1
    std::vector<const reco::Candidate*> gctCenUnskimmed                 = getJetCollection( iEvent, srcGctJetCentral_);
    std::vector<const reco::Candidate*> gctForUnskimmed                 = getJetCollection( iEvent, srcGctJetForward_);
#endif
#ifdef SIMULATION
    std::vector<const reco::Candidate*> genJetUnskimmed;
    std::vector<const reco::Candidate*> genJetForUnskimmed;
#endif

  std::vector<const reco::Candidate*> hltCaloJetUnskimmed, hltCaloJetForUnskimmed;
  std::vector<const reco::Candidate*> hltCaloJetIDUnskimmed, hltCaloJetIDForUnskimmed;
  std::vector<const reco::Candidate*> hltPFJetUnskimmed,   hltPFJetForUnskimmed;

  edm::Handle<reco::CaloJetCollection> hltCaloJetHandle;   iEvent.getByToken(hltCaloJetToken_, hltCaloJetHandle);
  if (hltCaloJetHandle.isValid()){ hltCaloJetUnskimmed     = getJetCollection( iEvent, hltCaloJetHandle); }
  edm::Handle<reco::CaloJetCollection> hltCaloJetIDHandle; iEvent.getByToken(hltCaloJetIDToken_, hltCaloJetIDHandle);
  if (hltCaloJetIDHandle.isValid()){ hltCaloJetIDUnskimmed = getJetCollection( iEvent, hltCaloJetIDHandle); }
  edm::Handle<reco::PFJetCollection> hltPFJetHandle;       iEvent.getByToken(hltPFJetToken_, hltPFJetHandle);
  if (hltPFJetHandle.isValid())  { hltPFJetUnskimmed       = getJetCollection( iEvent, hltPFJetHandle); }
  edm::Handle<reco::GenJetCollection> genJetHandle;        iEvent.getByToken(genJetToken_, genJetHandle);
  if (genJetHandle.isValid())  { genJetUnskimmed           = getJetCollection( iEvent, genJetHandle); }

  // Forward jets
  hltCaloJetForUnskimmed   = hltCaloJetUnskimmed;
  hltCaloJetIDForUnskimmed = hltCaloJetIDUnskimmed;
  hltPFJetForUnskimmed     = hltPFJetUnskimmed;
  genJetForUnskimmed       = genJetUnskimmed;
  

//     // std::vector<const reco::Candidate*> recoCaloUnskimmed            = getJetCollection( iEvent, srcCalo );
//     // std::vector<const reco::Candidate*> recoPFUnskimmed              = getJetCollection( iEvent, srcPF   );
//     // std::vector<const reco::Candidate*> recoCaloForUnskimmed         = getJetCollection( iEvent, srcCalo );
//     // std::vector<const reco::Candidate*> recoPFForUnskimmed           = getJetCollection( iEvent, srcPF   );

    

//     // Gen particles
//     // --------------------------------------------------------------------------------

//     genLeptonVeto      = false;
//     genElectronVeto    = false;
//     genMuonVeto        = false;
//     genPhotonVeto      = false;
//     genIsoLeptonVeto   = false;
//     genIsoElectronVeto = false;

//     std::vector<int> selLeptonIndices;
//     int particleIndex(0);

//     if (makeGenParticles){
//       edm::Handle< std::vector<reco::GenParticle> > genParticles;
//       iEvent.getByLabel(srcGenParticles_, genParticles);
//       for(std::vector<reco::GenParticle>::const_iterator iter = genParticles->begin(); iter != genParticles->end(); ++iter){
//     	const reco::GenParticle& genParticle = *iter;

//     	double genParticlePt  = genParticle.p4().pt();
//     	double genParticleEta = genParticle.p4().eta();
//     	double genParticlePhi = genParticle.p4().phi();


//     	if(TMath::Abs(genParticle.pdgId()) == 11){    	// Electron
//     	  if ( (genParticlePt >= genElectronMinPt) && (TMath::Abs(genParticleEta) <= genElectronMaxEta) ){
//     	    genElectronPt ->push_back( genParticlePt  );
//     	    genElectronEta->push_back( genParticleEta );
//     	    genElectronPhi->push_back( genParticlePhi );
//     	    genLeptonVeto   = true;
//     	    genElectronVeto = true;
// 	    selLeptonIndices.push_back( particleIndex );
//     	  }
//     	} // End lepton requirement
//     	else if(TMath::Abs(genParticle.pdgId()) == 13){ // Muon
//     	  if ( (genParticlePt >= genMuonMinPt) && (TMath::Abs(genParticleEta) <= genMuonMaxEta) ){
//     	    genMuonPt ->push_back( genParticlePt  );
//    	    genMuonEta->push_back( genParticleEta );
//     	    genMuonPhi->push_back( genParticlePhi );
//     	    genLeptonVeto = true;
// 	    genMuonVeto   = true;
// 	    selLeptonIndices.push_back( particleIndex );
//     	  }
//     	} // End Muon requirement
//     	else if(TMath::Abs(genParticle.pdgId()) == 22){ // Photon
//     	  if ( (genParticlePt >= genPhotonMinPt) && (TMath::Abs(genParticleEta) <= genPhotonMaxEta) ){
//     	    genPhotonPt ->push_back( genParticlePt  );
//    	    genPhotonEta->push_back( genParticleEta );
//     	    genPhotonPhi->push_back( genParticlePhi );
//     	    genPhotonVeto = true;
//     	  }
//     	} // End Photon requirement

// 	particleIndex++;
//       } // End loop
//     } // End gen particles




//     // Jet lepton cleaning
//     // --------------------------------------------------------------------------------

// #ifdef LEPTON_XCLEANING
//     if (genMuonVeto){
//       //crosscleanIsolatedLeptons( genMuonPt,  genMuonEta,  genMuonPhi, genJetUnskimmed, hltPFUnskimmed, nIsoMuons, 0.3, 40. );
//       crosscleanIsolatedLeptons( genMuonPt,  genMuonEta,  genMuonPhi, genJetUnskimmed, hltPFUnskimmed, hltCaloUnskimmed, 
// 				 nIsoMuons, 0.3, 40.,
// 				 hltPFJetUnskimmed,
// 				 genMuonMatchedGenMuonPt,
// 				 genMuonMatchedGenJetPt,
// 				 genMuonMatchedHLTPFJetPt,
// 				 genMuonMatchedHLTPFJetMuonEF,
// 				 genMuonMatchedHLTPFJetElectronEF);

//     }
//     if (genElectronVeto){
//       // crosscleanIsolatedLeptons( genElectronPt,  genElectronEta,  genElectronPhi, genJetUnskimmed, hltPFUnskimmed, 
//       // 				 nIsoElectrons, 0.3, 40. );
//       crosscleanIsolatedLeptons( genElectronPt,  genElectronEta,  genElectronPhi, genJetUnskimmed, hltPFUnskimmed,  hltCaloUnskimmed,
// 				 nIsoElectrons, 0.3, 40. ,
// 				 hltPFJetUnskimmed,
// 				 genElectronMatchedGenElectronPt,
//                                  genElectronMatchedGenJetPt,
//                                  genElectronMatchedHLTPFJetPt,
//                                  genElectronMatchedHLTPFJetMuonEF,
//                                  genElectronMatchedHLTPFJetElectronEF);

//     }
//     // std::cout << "\nAfter cleaning: " << genJetUnskimmed.size() << "\t" << hltPFUnskimmed.size() 
//     // 	      << "\t" << nIsoMuons << "\t" << nIsoElectrons << "\n";
// #endif




    // Skim jet collections
    // ----------------------------------------

#ifdef L1
    std::vector<const reco::Candidate*> gctFor              = skimJets(gctForUnskimmed,                 minPt, minEtaFor, maxEtaFor ); 
    std::vector<const reco::Candidate*> gctCen              = skimJets(gctCenUnskimmed,              minPt, minEtaCen, maxEtaCen );
#endif
#ifdef SIMULATION 
    std::vector<const reco::Candidate*> gen              = skimJets(genJetUnskimmed,                minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> genFor           = skimJets(genJetForUnskimmed,             minPt, minEtaFor, maxEtaFor );
#endif


    // Central jets
    // --------------------
    std::vector<const reco::Candidate*> hltCalo          = skimJets(hltCaloJetUnskimmed,          minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltCaloID        = skimJets(hltCaloJetIDUnskimmed,        minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltPF            = skimJets(hltPFJetUnskimmed,            minPt, minEtaCen, maxEtaCen );
    // std::vector<const reco::Candidate*> recoCalo         = skimJets(recoCaloUnskimmed,         minPt, minEtaCen, maxEtaCen );
    // std::vector<const reco::Candidate*> recoPF           = skimJets(recoPFUnskimmed,           minPt, minEtaCen, maxEtaCen );

    // Forward jets
    // --------------------
    std::vector<const reco::Candidate*> hltCaloFor       = skimJets(hltCaloJetForUnskimmed,          minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltCaloIDFor     = skimJets(hltCaloJetIDForUnskimmed,        minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltPFFor         = skimJets(hltPFJetForUnskimmed,            minPt, minEtaFor, maxEtaFor );
    // std::vector<const reco::Candidate*> recoCaloFor         = skimJets(recoCaloForUnskimmed,         minPt, minEtaFor, maxEtaFor );
    // std::vector<const reco::Candidate*> recoPFFor           = skimJets(recoPFForUnskimmed,           minPt, minEtaFor, maxEtaFor );


    // Jets
    // --------------------------------------------------------------------------------


    // Forward jet
#ifdef SIMULATION 
    genForMaxPt    = 0;
    if ( genFor.size()   > 0){    genForMaxPt =    genFor.at(0)->pt(); }
#endif
    //recoPFForMaxPt = 0;
    hltPFForMaxPt  = 0;
    //if (recoPFFor.size() > 0){ recoPFForMaxPt = recoPFFor.at(0)->pt(); }
    if ( hltPFFor.size() > 0){  hltPFForMaxPt =  hltPFFor.at(0)->pt(); }

    // Lead jet 
#ifdef SIMULATION 
    genLeadJetPt= 0;
    if ( gen.size()     > 0){ genLeadJetPt      =     gen.at(0)->pt(); }
#endif
    //recoPFLeadJetPt= 0;
    hltPFLeadJetPt= 0;
    hltCaloLeadJetPt= 0;
    hltCaloIDLeadJetPt= 0;
    //if (recoPF.size()   > 0){ recoPFLeadJetPt   =  recoPF.at(0)->pt(); }
    if ( hltPF.size()   > 0){ hltPFLeadJetPt    =   hltPF.at(0)->pt(); }
    if ( hltCalo.size() > 0){ hltCaloLeadJetPt  = hltCalo.at(0)->pt(); }
    if ( hltCaloID.size() > 0){ hltCaloIDLeadJetPt  = hltCaloID.at(0)->pt(); }

    // Second jet
#ifdef SIMULATION 
    genSecondJetPt= 0;
    if ( gen.size()     > 1){ genSecondJetPt      =     gen.at(1)->pt(); }
#endif
    //    recoPFSecondJetPt= 0;
    hltPFSecondJetPt= 0;
    hltCaloSecondJetPt= 0;
    hltCaloIDSecondJetPt= 0;
    //    if (recoPF.size()   > 1){ recoPFSecondJetPt   =  recoPF.at(1)->pt(); }
    if ( hltPF.size()   > 1){ hltPFSecondJetPt    =   hltPF.at(1)->pt(); }
    if ( hltCalo.size() > 1){ hltCaloSecondJetPt  = hltCalo.at(1)->pt(); }
    if ( hltCaloID.size() > 1){ hltCaloIDSecondJetPt  = hltCaloID.at(1)->pt(); }

    // Dijet avg 
#ifdef SIMULATION
    genDijetAvePt= 0;
    if ( gen.size()     > 1){ genDijetAvePt      = 0.5*(gen.at(0)->pt()     + gen.at(1)->pt()); }
#endif
    //    recoPFDijetAvePt= 0;
    hltPFDijetAvePt    = 0;
    hltCaloDijetAvePt  = 0;
    hltCaloIDDijetAvePt= 0;
    //    if (recoPF.size()   > 1){ recoPFDijetAvePt   = 0.5*(recoPF.at(0)->pt()  + recoPF.at(1)->pt()); }
    if ( hltPF.size()   > 1){ hltPFDijetAvePt    = 0.5*(hltPF.at(0)->pt()   + hltPF.at(1)->pt()); }
    if ( hltCalo.size() > 1){ hltCaloDijetAvePt  = 0.5*(hltCalo.at(0)->pt() + hltCalo.at(1)->pt()); }
    if ( hltCaloID.size() > 1){ hltCaloIDDijetAvePt  = 0.5*(hltCaloID.at(0)->pt() + hltCaloID.at(1)->pt()); }

    // ********************************************************************************
    // *                                  Energy sums                                 *
    // ********************************************************************************

    // hltMet
    // getValue(iEvent, srcHLTMetCalo_,                metPt_["hltMetCalo"],                metPhi_["hltMetCalo"]);
    // getValue(iEvent, srcHLTMetCleanCalo_,           metPt_["hltMetCleanCalo"],           metPhi_["hltMetCleanCalo"]);
    // getValue(iEvent, srcHLTMetCleanUsingJetIDCalo_, metPt_["hltMetCleanUsingJetIDCalo"], metPhi_["hltMetCleanUsingJetIDCalo"]);

    edm::Handle<reco::CaloMETCollection> hltCaloMetHandle;  iEvent.getByToken(hltCaloMetToken_, hltCaloMetHandle);    
    getValue(iEvent, hltCaloMetHandle,              metPt_["hltMetCalo"], metPhi_["hltMetCalo"]);

    edm::Handle<reco::METCollection> hltPFMetHandle;  iEvent.getByToken(hltPFMetToken_, hltPFMetHandle);    
    getValue(iEvent, hltPFMetHandle,                metPt_["hltMetPF"],   metPhi_["hltMetPF"]);

#ifdef SIMULATION
    edm::Handle<reco::GenMETCollection> genMetCaloHandle;  iEvent.getByToken(genMetCaloToken_, genMetCaloHandle);    
    getValue(iEvent, genMetCaloHandle,              metPt_["genMetCalo"], metPhi_["genMetCalo"]);
    edm::Handle<reco::GenMETCollection> genMetTrueHandle;  iEvent.getByToken(genMetTrueToken_, genMetTrueHandle);    
    getValue(iEvent, genMetTrueHandle,              metPt_["genMetTrue"], metPhi_["genMetTrue"]);
#endif



    // ********************************************************************************
    //                            Produce jet event quantities
    // ********************************************************************************

    // HT and AlphaT
#ifdef SIMULATION
    genAlphaTHT40              = calculateAlphaTHT( gen,      40.);
#endif
    hltCaloAlphaTHT40          = calculateAlphaTHT( hltCalo,  40.);
    hltCaloIDAlphaTHT40        = calculateAlphaTHT( hltCaloID,  40.);
    hltPFAlphaTHT40            = calculateAlphaTHT( hltPF,    40.);
    // recoCaloAlphaTHT40         = calculateAlphaTHT( recoCalo, 40.);
    // recoPFAlphaTHT40           = calculateAlphaTHT( recoPF,   40.);

//     // Dynamic HT and AlphaT 
//     // genDynamicAlphaTHT40    = calculateDynamicAlphaTPairs( gen,    dynamicJetThreshold );
//     // hltPFDynamicAlphaTHT40  = calculateDynamicAlphaTPairs( hltPF,  dynamicJetThreshold );
// #ifdef SIMULATION
//     calculateDynamicAlphaTPairs( gen, dynamicJetThreshold, genDynamicAlphaT40, genDynamicHT40 );
// #endif
//     calculateDynamicAlphaTPairs( hltPF,  dynamicJetThreshold, hltPFDynamicAlphaT40, hltPFDynamicHT40 );

//     //    recoPFDynamicAlphaTHT40 = calculateDynamicAlphaTPairs( recoPF, dynamicJetThreshold );

//     // ******************************************************************************** 
//     // Calculate: DeltaR of leading jets L1-Gen, HLT-Gen
//     // ******************************************************************************** 

#ifdef SIMULATION
//     L1GenDeltaR  = leadL1GenDeltaR(    gctCenUnskimmed,  gctForUnskimmed, genJetUnskimmed );
//     HLTGenDeltaR = leadHLTGenDeltaR( hltPFUnskimmed, genJetUnskimmed );
//     hpuVeto      = (L1GenDeltaR < 0.5);
    gen_AlphaTPrime40     = calculateAlphaTPrime( genMHT40.first,      genAlphaTHT40.second );
    genMHT40              = calculateMHT( gen,      40.);
    genForMHT40           = calculateMHT( genFor,      40.);
#endif

    // MHT
    hltCaloMHT40          = calculateMHT( hltCalo,    40.);
    hltCaloForMHT40       = calculateMHT( hltCaloFor, 40.);
    hltCaloIDMHT40        = calculateMHT( hltCaloID,  40.);
    hltPFMHT40            = calculateMHT( hltPF,      40.);
    hltPFForMHT40         = calculateMHT( hltPFFor,   40.);
    //recoPFMHT40         = calculateMHT( recoPF,     40.);
    //recoPFForMHT40      = calculateMHT( recoPFFor,  40.);

    // AlphaT prime
    hltCalo_AlphaTPrime40    = calculateAlphaTPrime( hltCaloMHT40.first,  hltCaloAlphaTHT40.second );
    hltCaloID_AlphaTPrime40  = calculateAlphaTPrime( hltCaloIDMHT40.first,  hltCaloIDAlphaTHT40.second );
    hltPF_AlphaTPrime40      = calculateAlphaTPrime( hltPFMHT40.first,    hltPFAlphaTHT40.second );
    // recoCalo_AlphaTPrime40   = calculateAlphaTPrime( recoCaloMHT40.first, recoCaloAlphaTHT40.second );
    // recoPF_AlphaTPrime40     = calculateAlphaTPrime( recoPFMHT40.first,   recoPFAlphaTHT40.second );


    // ****************************************
    // MHT-MET deltaPhi
    // ****************************************
    hltMetCaloPFMht40_DeltaPhi = fabs( deltaPhi( metPhi_["hltMetCalo"], hltPFMHT40.second ) );
    //genMetCaloMht40_DeltaPhi   = fabs( deltaPhi( metPhi_["genMetCalo"], genMHT40.second )   );

//     // ----------------------------------------
//     // Analysis binning
//     // ----------------------------------------
    
// #ifdef SIMULATION
//     genNJet40              = calculateNJet( gen,      40.);
//     genNJetBin40           = calculateNJetBin( gen,      40.);
//     genHTBin40           = calculateHTBin( genAlphaTHT40.second,      genAlphaTHT40.first );
// #endif
//     hltCaloIDNJet40        = calculateNJet( hltCaloID,  40.);
//     hltCaloNJet40          = calculateNJet( hltCalo,  40.);
//     hltPFNJet40            = calculateNJet( hltPF,    40.);
//     // recoCaloNJet40         = calculateNJet( recoCalo, 40.);
//     // recoPFNJet40           = calculateNJet( recoPF,   40.);

//     hltCaloNJetBin40       = calculateNJetBin( hltCalo,  40.);
//     hltCaloIDNJetBin40     = calculateNJetBin( hltCaloID,  40.);
//     hltPFNJetBin40         = calculateNJetBin( hltPF,    40.);
//     // recoCaloNJetBin40      = calculateNJetBin( recoCalo, 40.);
//     // recoPFNJetBin40        = calculateNJetBin( recoPF,   40.);

//     hltCaloHTBin40       = calculateHTBin( hltCaloAlphaTHT40.second,  hltCaloAlphaTHT40.first );
//     hltCaloIDHTBin40     = calculateHTBin( hltCaloIDAlphaTHT40.second,  hltCaloIDAlphaTHT40.first );
//     hltPFHTBin40         = calculateHTBin( hltPFAlphaTHT40.second,    hltPFAlphaTHT40.first );
//     // recoCaloHTBin40      = calculateHTBin( recoCaloAlphaTHT40.second, recoCaloAlphaTHT40.first );
//     // recoPFHTBin40        = calculateHTBin( recoPFAlphaTHT40.second,   recoPFAlphaTHT40.first );


    // ----------------------------------------
    // Store jet collections
    // ----------------------------------------

#ifdef SIMULATION
    storeJet( "gen",                 gen,                 jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genFor",              genFor,              jetPt, jetPx, jetPy, jetEta, jetPhi );
#endif

    storeJet( "hltCalo",             hltCalo,             jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltCaloID",           hltCaloID,           jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltPF",               hltPF,               jetPt, jetPx, jetPy, jetEta, jetPhi );
    // storeJet( "recoCalo",            recoCalo,            jetPt, jetPx, jetPy, jetEta, jetPhi );
    // storeJet( "recoPF",              recoPF,              jetPt, jetPx, jetPy, jetEta, jetPhi );

    storeJet( "hltCaloFor",          hltCaloFor,          jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltCaloIDFor",        hltCaloIDFor,          jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltPFFor",            hltPFFor,            jetPt, jetPx, jetPy, jetEta, jetPhi );
    // storeJet( "recoCaloFor",         recoCaloFor,         jetPt, jetPx, jetPy, jetEta, jetPhi );
    // storeJet( "recoPFFor",           recoPFFor,           jetPt, jetPx, jetPy, jetEta, jetPhi );

    
//     // ----------------------------------------
//     // Perform jet matching
//     // ----------------------------------------
// #ifdef SIMULATION
//     genMatchedHLTPF  = matchJetCollections( jetPt["gen"],   jetEta["gen"],   jetPhi["gen"],
//     					       jetPt["hltPF"], jetEta["hltPF"], jetPhi["hltPF"], 20., 0.25 );
// #endif
//     // genMatchedRECOPF = matchJetCollections( jetPt["gen"],    jetEta["gen"],    jetPhi["gen"],
//     // 					       jetPt["recoPF"], jetEta["recoPF"], jetPhi["recoPF"], 20., 0.25 );
    

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
MakeTrees::leadL1GenDeltaR( const std::vector<const reco::Candidate*>& gctCen, const std::vector<const reco::Candidate*>& gctFor, const std::vector<const reco::Candidate*>& gen){

  const reco::Candidate* leadL1Jet  = NULL;
  const reco::Candidate* leadGenJet = NULL;
  float deltaR(20.);

  // Get lead genJet in acceptance
  for ( std::vector<const reco::Candidate*>::const_iterator itr = gen.begin(); itr != gen.end(); ++itr ){
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
MakeTrees::leadHLTGenDeltaR( const std::vector<const reco::Candidate*>& hlt,
			     const std::vector<const reco::Candidate*>& gen){

  const reco::Candidate* leadHLTJet = NULL;
  const reco::Candidate* leadGenJet = NULL;
  float deltaR(20.);

  // Get lead genJet in acceptance
  for ( std::vector<const reco::Candidate*>::const_iterator itr = gen.begin(); itr != gen.end(); ++itr ){
    if ( fabs((*itr)->eta()) < 5. ){
      leadGenJet = (*itr); break;
    }
  }
  if ( leadGenJet == NULL )    { return deltaR; }
  if ( leadGenJet->pt() < 10. ){ return deltaR; }

  // Get lead hltJet
  if (hlt.size() > 0){ leadHLTJet = hlt.at(0); }
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


 

DEFINE_FWK_MODULE(MakeTrees);
