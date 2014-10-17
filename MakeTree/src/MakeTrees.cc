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
    VInputTag srcGen5Jet_;

    edm::InputTag srcGenMetCalo_;
    edm::InputTag srcGenMetCaloAndNonPrompt_;
    edm::InputTag srcGenMetTrue_;

    edm::InputTag srcGenParticles_;
    bool          makeGenParticles;
    double        genElectronMinPt;
    double        genElectronMaxEta;
    double        genMuonMinPt;
    double        genMuonMaxEta;



    VInputTag srcHLTAk4Calo;
    VInputTag srcHLTAk4PF;
    VInputTag srcHLTAk4PFNoPU;


    VInputTag srcCaloJet_;
    VInputTag srcPfJet_;


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

    std::pair<float,float> genAk4AlphaTHT50;
    std::pair<float,float> hltAk4PFAlphaTHT50;
    std::pair<float,float> hltAk4CaloAlphaTHT50;

    UInt_t maxjet_;
    bool usePU_; 
    UInt_t run_;
    UInt_t lumi_;
    ULong64_t event_;
    bool jetVeto50_;
    bool jetVeto30_;
    //Parameters for the HSums
    //double htThreshold_;

    // Jet skim cuts
    double  minPt;
    double minEtaCen;
    double maxEtaCen;
    double minEtaFor;
    double maxEtaFor;


    edm::InputTag HLTResultsTag;
  
    float PThat;



  // /// module config parameters
  // std::string   processName_;
  // std::string   triggerName_;

  // edm::InputTag                                   triggerResultsTag_;
  // edm::EDGetTokenT<edm::TriggerResults>           triggerResultsToken_;
  // edm::InputTag                                   triggerEventWithRefsTag_;
  // edm::EDGetTokenT<trigger::TriggerEventWithRefs> triggerEventWithRefsToken_;

  // edm::Handle<edm::TriggerResults>           triggerResultsHandle_;
  // edm::Handle<trigger::TriggerEventWithRefs> triggerEventWithRefsHandle_;
  // trigger::Vids        jetIds_;

  // HLTConfigProvider hltConfig_;

};

MakeTrees::MakeTrees(const edm::ParameterSet& pset){
// : 
//   processName_(pset.getParameter<std::string>("processName")),
//   triggerName_(pset.getParameter<std::string>("triggerName")),
//   triggerResultsTag_(pset.getParameter<edm::InputTag>("triggerResults")),
//   triggerResultsToken_(consumes<edm::TriggerResults>(triggerResultsTag_)),
//   triggerEventWithRefsTag_(pset.getParameter<edm::InputTag>("triggerEventWithRefs")),
//   triggerEventWithRefsToken_(consumes<trigger::TriggerEventWithRefs>(triggerEventWithRefsTag_)) {

    // Initialize the ntuple builder
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("Ntuple", "Ntuple");




    lvl_.push_back("gctCen");
    lvl_.push_back("genAk4");
    lvl_.push_back("genAk5");
    lvl_.push_back("hltAk4Calo");
    lvl_.push_back("hltAk4PF");


    lvl_.push_back("gctFor");
    lvl_.push_back("genAk4For");
    lvl_.push_back("genAk5For");
    lvl_.push_back("hltAk4CaloFor");
    lvl_.push_back("hltAk4PFFor");

    // lvl_.push_back("hltAk4PFNoPU");
    // lvl_.push_back("Calo");
    // lvl_.push_back("Pf");



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

    // Event
    tree->Branch("NVTX",  &NVTX,  "NVTX/i");
    tree->Branch("PThat", &PThat, "PThat/f");


    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){

	mhtPt_[*iLvl]        = -10.0;
	mhtPhi_[*iLvl]       = -10.0;
	ht_[*iLvl]           = -10.0;
	dht_[*iLvl]           = -10.0;
	alphaT_[*iLvl]       = 0.0;
	mhtDivHt_[*iLvl]     = 0.0;

	metPt_[*iLvl]        = -10.0;
	metPhi_[*iLvl]       = -10.0;
	et_[*iLvl]       = -10.0;
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



    // AlphaT, HT  
    tree->Branch("genAk4_AlphaT40",     &genAk4AlphaTHT40.first,      "genAk4_AlphaT40/f");
    tree->Branch("genAk4_HT40",         &genAk4AlphaTHT40.second,     "genAk4_HT40/f");
    tree->Branch("hltAk4PF_AlphaT40",   &hltAk4PFAlphaTHT40.first,    "hltAk4PF_AlphaT40/f");
    tree->Branch("hltAk4PF_HT40",       &hltAk4PFAlphaTHT40.second,   "hltAk4PF_HT40/f");
    tree->Branch("hltAk4Calo_AlphaT40", &hltAk4CaloAlphaTHT40.first,  "hltAk4Calo_AlphaT40/f");
    tree->Branch("hltAk4Calo_HT40",     &hltAk4CaloAlphaTHT40.second, "hltAk4Calo_HT40/f");

    tree->Branch("genAk4_AlphaT50",     &genAk4AlphaTHT50.first,      "genAk4_AlphaT50/f");
    tree->Branch("genAk4_HT50",         &genAk4AlphaTHT50.second,     "genAk4_HT50/f");
    tree->Branch("hltAk4PF_AlphaT50",   &hltAk4PFAlphaTHT50.first,    "hltAk4PF_AlphaT50/f");
    tree->Branch("hltAk4PF_HT50",       &hltAk4PFAlphaTHT50.second,   "hltAk4PF_HT50/f");
    tree->Branch("hltAk4Calo_AlphaT50", &hltAk4CaloAlphaTHT50.first,  "hltAk4Calo_AlphaT50/f");
    tree->Branch("hltAk4Calo_HT50",     &hltAk4CaloAlphaTHT50.second, "hltAk4Calo_HT50/f");

    // Energy sums
    tree->Branch("gct_Ht",       &ht_["gct"],      "gct_Ht/f");
    tree->Branch("gct_MhtPt",    &mhtPt_["gct"],   "gct_MhtPt/f");
    tree->Branch("gct_MhtPhi",   &mhtPhi_["gct"],  "gct_MhtPhi/f");
    tree->Branch("gct_MhtDivHt", &mhtDivHt_["gct"],"GCT Mht divided by Ht/f"); 

    tree->Branch("gct_Et",     &et_["gct"],     "gct_Et/f");
    tree->Branch("gct_MetPt",  &metPt_["gct"],  "gct_MetPt/f");
    tree->Branch("gct_MetPhi", &metPhi_["gct"], "gct_MetPhi/f");

    tree->Branch("genMetCalo_MetPt",              &metPt_["genMetCalo"],             "genMetCalo_MetPt/f");
    tree->Branch("genMetCaloAndNonPrompt_MetPt",  &metPt_["genMetCaloAndNonPrompt"], "genMetCaloAndNonPrompt_MetPt/f");
    tree->Branch("genMetTrue_MetPt",              &metPt_["genMetTrue"],             "genMetTrue_MetPt/f");


    // recoCaloMETs_hltMet
    tree->Branch("hltMetCalo_MetPT",              &metPt_["hltMetCalo"],             "hltMetCalo_MetPt/f");
    // recoMETs_hltHtMht
    tree->Branch("hltMhtCalo_MetPT",              &mhtPt_["hltMhtCalo"],             "hltMhtCalo_MetPt/f");
    // recoMETs_hltPFHT
    tree->Branch("hltMhtPF_MetPT",                &mhtPt_["hltMhtPF"],               "hltMhtPF_MetPt/f");

    // NOTE: These are not always filled due to threshold requirements
    // recoCaloMETs_hltMetClean
    tree->Branch("hltMetCleanCalo_MetPT",         &metPt_["hltMetCleanCalo"],        "hltMetCleanCalo_MetPt/f");


    // Gen leptons
    tree->Branch("genLeptonVeto",   &genLeptonVeto,       "genLeptonVeto/b");
    tree->Branch("genElectron_Pt",  "std::vector<float>", &genElectronPt);
    tree->Branch("genElectron_Eta", "std::vector<float>", &genElectronEta);
    tree->Branch("genElectron_Phi", "std::vector<float>", &genElectronPhi);		
    tree->Branch("genMuon_Pt",      "std::vector<float>", &genMuonPt);
    tree->Branch("genMuon_Eta",     "std::vector<float>", &genMuonEta);
    tree->Branch("genMuon_Phi",     "std::vector<float>", &genMuonPhi);		



    // Store HLT paths
    // hltPathNames.push_back("HLT_CaloJet20_v1");
    // hltPathNames.push_back("HLT_PFJet20_v1");
    // hltPathNames.push_back("HLT_HT100_v1");
    // hltPathNames.push_back("HLT_PFHT100_v1");
    // hltPathNames.push_back("HLT_HT200_AlphaT0p57_v1");
    // hltPathNames.push_back("HLT_HT250_AlphaT0p55_v1");
    // hltPathNames.push_back("HLT_HT300_AlphaT0p53_v1");
    // hltPathNames.push_back("HLT_HT350_AlphaT0p52_v1");
    // hltPathNames.push_back("HLT_HT400_AlphaT0p51_v1");
    // hltPathNames.push_back("HLT_HT200_AlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT250_AlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT300_AlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT350_AlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT200_PFAlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT250_PFAlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT300_PFAlphaT0p5_v1");
    // hltPathNames.push_back("HLT_HT350_PFAlphaT0p5_v1");
    // hltPathNames.push_back("HLT_PFHT350_v1");
    // hltPathNames.push_back("HLT_PFHT600_v1");
    // hltPathNames.push_back("HLT_PFHT900_v1");
    // hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_v1");
    // hltPathNames.push_back("HLT_PFHT350_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT600_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT900_NoL1_v1");
    // hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_NoL1_v1");

    hltPathNames.push_back("HLT_CaloJet20_v1");
    hltPathNames.push_back("HLT_PFJet20_v1");
    hltPathNames.push_back("HLT_HT100_v1");
    hltPathNames.push_back("HLT_PFHT100_v1");
    
    hltPathNames.push_back("HLT_HT200_AlphaT0p57_NoL1_v1");
    hltPathNames.push_back("HLT_HT250_AlphaT0p55_NoL1_v1");
    hltPathNames.push_back("HLT_HT300_AlphaT0p53_NoL1_v1");
    hltPathNames.push_back("HLT_HT350_AlphaT0p52_NoL1_v1");
    hltPathNames.push_back("HLT_HT400_AlphaT0p51_NoL1_v1");
    hltPathNames.push_back("HLT_HT200_AlphaT0p57_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT250_AlphaT0p55_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT300_AlphaT0p53_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT350_AlphaT0p52_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT400_AlphaT0p51_L1HTT175OrETM70_v1");
    
    hltPathNames.push_back("HLT_HT200_AlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT250_AlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT300_AlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT350_AlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT200_AlphaT0p5_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT250_AlphaT0p5_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT300_AlphaT0p5_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT350_AlphaT0p5_L1HTT175OrETM70_v1");

    hltPathNames.push_back("HLT_HT200_PFAlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT250_PFAlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT300_PFAlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT350_PFAlphaT0p5_NoL1_v1");
    hltPathNames.push_back("HLT_HT200_PFAlphaT0p5_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT250_PFAlphaT0p5_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT300_PFAlphaT0p5_L1HTT175OrETM70_v1");
    hltPathNames.push_back("HLT_HT350_PFAlphaT0p5_L1HTT175OrETM70_v1");
          
    hltPathNames.push_back("HLT_PFHT350_NoL1_v1");
    hltPathNames.push_back("HLT_PFHT600_NoL1_v1");
    hltPathNames.push_back("HLT_PFHT350_v1");
    hltPathNames.push_back("HLT_PFHT600_v1");
    
    hltPathNames.push_back("HLT_PFHT900_NoL1_v1");
    hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_NoL1_v1");
    hltPathNames.push_back("HLT_PFMET170_NoiseCleaned_NoL1_v1");
    hltPathNames.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_NoL1_v1");
    hltPathNames.push_back("HLT_PFHT900_v1");
    hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_v1");
    hltPathNames.push_back("HLT_PFMET170_NoiseCleaned_v1");
    hltPathNames.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v1");
    
    hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_NoL1_v1");
    hltPathNames.push_back("HLT_RsqMR260_Rsq0p09_MR200_4jet_NoL1_v1");
    hltPathNames.push_back("HLT_Rsq0p36_NoL1_v1");
    hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_v1");
    hltPathNames.push_back("HLT_RsqMR260_Rsq0p09_MR200_4jet_v1");
    hltPathNames.push_back("HLT_Rsq0p36_v1");
    



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

    srcGenMetCalo_             = pset.getParameter<edm::InputTag>("srcGenMetCalo");
    srcGenMetCaloAndNonPrompt_ = pset.getParameter<edm::InputTag>("srcGenMetCaloAndNonPrompt");
    srcGenMetTrue_             = pset.getParameter<edm::InputTag>("srcGenMetTrue");

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
    srcGen5Jet_        = pset.getParameter<VInputTag>("srcGen5Jet");

    srcHLTAk4Calo      = pset.getParameter<VInputTag>("srcHLTAk4Calo");
    srcHLTAk4PF        = pset.getParameter<VInputTag>("srcHLTAk4PF");
    //    srcHLTAk4PFNoPU    = pset.getParameter<VInputTag>("srcHLTAk4PFNoPU");


    srcCaloJet_       = pset.getParameter<VInputTag>("srcCaloJet");
    srcPfJet_         = pset.getParameter<VInputTag>("srcPfJet");

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

    void getValue(const edm::Event& evt, const edm::InputTag& tag, Float_t& et, Float_t& phi) {
	edm::Handle<edm::View<reco::Candidate> > handle;
	evt.getByLabel(tag, handle);
	et = handle->at(0).pt();
	phi = handle->at(0).phi();
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


    inline float deltaPhi( float phi1, float phi2 ){

	float const  PI        = ROOT::Math::Pi();
	float const  TWOPI     = 2.*PI;
	float dPhi = (phi1 - phi2);
	while (dPhi >= PI) dPhi -= TWOPI;
	while (dPhi < -PI) dPhi += TWOPI;
	return dPhi;

    }



  std::pair<float, float> calculateAlphaTHT(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
     
      // Momentum sums in transverse plane
      float sum_et(0), sum_px(0), sum_py(0);

      // check the size of the input collection
      if (jets.size() <= 1){
	return std::make_pair(0., sum_et);
      }


      // Jet collection restricted to jet threshold
      std::vector<float> jetPTNew;

      for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){

	if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
	jetPTNew.push_back( jets.at(iJet)->pt() );

	sum_et += jets.at(iJet)->pt();
	sum_px += jets.at(iJet)->px();
	sum_py += jets.at(iJet)->py();

      }
      // check the size of the new input collection 
      if (jetPTNew.size() <= 1){
	// empty jet collection, return AlphaT = 0 
	return std::make_pair( 0., sum_et);
      }

      // Minimum Delta Et for two pseudo-jets 
      double min_delta_sum_et = sum_et;

      for (unsigned int i = 0; i < (1U << (jetPTNew.size() - 1)); i++) { //@@ iterate through different combinations 
	double delta_sum_et = 0.;

	for (unsigned int j = 0; j < jetPTNew.size(); ++j) { //@@ iterate through jets 
	  if (i & (1U << j))
	    delta_sum_et -= jetPTNew[j];
	  else
	    delta_sum_et += jetPTNew[j];
	}
	delta_sum_et = std::abs(delta_sum_et);
	if (delta_sum_et < min_delta_sum_et) {
	  min_delta_sum_et = delta_sum_et;
	}

      }

      // Return a large value of alphaT 
      if ( (sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py)) <= 0 ){
	return std::make_pair(10.,  sum_et);
      }

      // Alpha_T 
     return std::make_pair( (0.5 * (sum_et - min_delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) )), sum_et );


    }

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



  // // Get UCT jets
  // edm::Handle< BXVector<l1t::Jet>> uctCalibJets;
  // iEvent.getByLabel(srcUctJet_, uctCalibJets);

  // std::vector<const reco::Candidate*> uctCenUnskimmed;
  // // std::vector<const reco::Candidate*> uctJetAll;

  // std::cout << "UCTBX:\n";
  // for ( auto itr = uctCalibJets->begin(0); itr != uctCalibJets->end(0); ++itr ) { 
  //   if(fabs(itr->eta())<=3.0){ 

  //     std::cout <<itr->pt() << "\t" <<  itr->eta() << "\t" << itr->phi() << "\n";
  //     //      uctCenUnskimmed.push_back(&(*itr)); 
  //   } 
  //   //    uctJetAll.push_back(&(*itr)); 
  // } 


  // edm::Handle< BXVector<l1t::Jet>> uctCalibJets;
  //   srcUctMET_ = pset.getParameter<edm::InputTag>("srcUctMet");
  //   srcUctMht_ = pset.getParameter<edm::InputTag>("srcUctMht");



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
  for (unsigned i = 0; i < paths->size(); ++i) {
    std::string name = paths->at(i).name();
    //std::cout << name << "\n";


    if ( hltPathFired.find( name ) != hltPathFired.end() ){
      //      std::cout << name << " found\n";
      hltPathFired[ name ] = paths->at(i).wasAccept();
    }

    // if ( name.find("HLT_PFHT900_v1")                             != name.npos ){
    //   HLT_PFHT900_v1 = paths->at(i).wasAccept();
    // }else if ( name.find("HLT_PFHT350_PFMET120_NoiseCleaned_v1") != name.npos ){
    //   HLT_PFHT350_PFMET120_NoiseCleaned_v1 = paths->at(i).wasAccept();
    // }else if ( name.find("HLT_HT200_AlphaT0p4_v1")               != name.npos ){
    //   HLT_HT200_AlphaT0p4_v1 = paths->at(i).wasAccept();
    // }else if ( name.find("HLT_HT200_PFAlphaT0p4_v1")             != name.npos ){
    //   HLT_HT200_PFAlphaT0p4_v1 = paths->at(i).wasAccept();
    // }

  }

   // std::cout << "Tree: " << HLT_PFHT900_v1 << "\t" << HLT_PFHT350_PFMET120_NoiseCleaned_v1 << "\t" 
   // 	    << HLT_HT200_AlphaT0p4_v1 << "\t" << HLT_HT200_PFAlphaT0p4_v1    << "\n";




  // Extracting triggerRefs

  // // get event products
  // iEvent.getByToken(triggerResultsToken_,triggerResultsHandle_);
  // if (!triggerResultsHandle_.isValid()) {
  //   std::cout << "HLTEventAnalyzerRAW::analyze: Error in getting TriggerResults product from Event!" << std::endl;
  //   return;
  // }
  // iEvent.getByToken(triggerEventWithRefsToken_,triggerEventWithRefsHandle_);
  // if (!triggerEventWithRefsHandle_.isValid()) {
  //   std::cout << "HLTEventAnalyzerRAW::analyze: Error in getting TriggerEventWithRefs product from Event!" << std::endl;
  //   return;
  // }



  // const unsigned int n(hltConfig_.size());
  // const unsigned int triggerIndex(hltConfig_.triggerIndex(triggerName));
  // assert(triggerIndex==iEvent.triggerNames(*triggerResultsHandle_).triggerIndex(triggerName));

  // // abort on invalid trigger name
  // if (triggerIndex>=n) {
  //   std::cout << "HLTEventAnalyzerRAW::analyzeTrigger: path "
  // 	      << triggerName << " - not found!" << std::endl;
  //   return;
  // }
  
  // std::cout << "HLTEventAnalyzerRAW::analyzeTrigger: path "
  // 	    << triggerName << " [" << triggerIndex << "]" << std::endl;
  // // modules on this trigger path
  // const unsigned int m(hltConfig_.size(triggerIndex));
  // const vector<string>& moduleLabels(hltConfig_.moduleLabels(triggerIndex));

  // // Results from TriggerResults product
  // std::cout << " Trigger path status:"
  //      << " WasRun=" << triggerResultsHandle_->wasrun(triggerIndex)
  //      << " Accept=" << triggerResultsHandle_->accept(triggerIndex)
  //      << " Error =" << triggerResultsHandle_->error(triggerIndex)
  // 	    << std::endl;


  // const unsigned int moduleIndex(triggerResultsHandle_->index(triggerIndex));
  // std::cout << " Last active module - label/type: "
  //      << moduleLabels[moduleIndex] << "/" << hltConfig_.moduleType(moduleLabels[moduleIndex])
  //      << " [" << moduleIndex << " out of 0-" << (m-1) << " on this path]"
  //      << std::endl;

  // for (unsigned int j=0; j<=moduleIndex; ++j) {
  //   const string& moduleLabel(moduleLabels[j]);
  //   const string  moduleType(hltConfig_.moduleType(moduleLabel));
  //   // check whether the module is packed up in TriggerEventWithRef product
  //   const unsigned int filterIndex(triggerEventWithRefsHandle_->filterIndex(InputTag(moduleLabel,"",processName_)));
  //   if (filterIndex<triggerEventWithRefsHandle_->size()) {
  //     std::cout << " Filter in slot " << j << " - label/type " << moduleLabel << "/" << moduleType << endl;
  //     std::cout << " Filter packed up at: " << filterIndex << endl;
  //     std::cout << "  Accepted objects:" << endl;



  //     jetIds_.clear();
  //     triggerEventWithRefsHandle_->getObjects(filterIndex,jetIds_,jetRefs_);
  //     const unsigned int nJets(jetIds_.size());
  //     if (nJets>0) {
  // 	std::cout << "   Jets: " << nJets << "  - the objects: # id pt" << std::endl;
  // 	for (unsigned int i=0; i!=nJets; ++i) {
  // 	  std::cout << "   " << i << " " << jetIds_[i]
  // 	       << " " << jetRefs_[i]->pt()
  // 	       << std::endl;
  // 	}
  //     }
      
      
      
  //   }
  // }




  // ------------------------------------------------------------------------------------------------------------------------
  // ------------------------------------------------------------------------------------------------------------------------
    
    // Input jets without eta or pT requirements
    // --------------------------------------------------------------------------------
    // std::vector<const reco::Candidate*> s2NopusJetCentral    = getCollections( iEvent, srcS2NopusJetCentral_);
    // std::vector<const reco::Candidate*> s2GlobalJetCentral    = getCollections( iEvent, srcS2GlobalJetCentral_);
    // std::vector<const reco::Candidate*> s2DonutJetCentral    = getCollections( iEvent, srcS2DonutJetCentral_);
    std::vector<const reco::Candidate*> gctCenUnskimmed       = getCollections( iEvent, srcGctJetCentral_);
    // std::vector<const reco::Candidate*> gctJetAll        = getCollections( iEvent, srcGctJetAll_);
    std::vector<const reco::Candidate*> genJet4Unskimmed      = getCollections( iEvent, srcGen4Jet_);
    std::vector<const reco::Candidate*> genJet5Unskimmed      = getCollections( iEvent, srcGen5Jet_);

    std::vector<const reco::Candidate*> hltAk4CaloUnskimmed   = getCollections( iEvent, srcHLTAk4Calo );
    std::vector<const reco::Candidate*> hltAk4PFUnskimmed     = getCollections( iEvent, srcHLTAk4PF );


    std::vector<const reco::Candidate*> gctForUnskimmed          = getCollections( iEvent, srcGctJetForward_);
    std::vector<const reco::Candidate*> genJet4ForUnskimmed      = getCollections( iEvent, srcGen4Jet_);
    std::vector<const reco::Candidate*> genJet5ForUnskimmed      = getCollections( iEvent, srcGen5Jet_);
    std::vector<const reco::Candidate*> hltAk4CaloForUnskimmed   = getCollections( iEvent, srcHLTAk4Calo );
    std::vector<const reco::Candidate*> hltAk4PFForUnskimmed     = getCollections( iEvent, srcHLTAk4PF );

    //    std::vector<const reco::Candidate*> hltAk4PFNoPUUnskimmed = getCollections( iEvent, srcHLTAk4PFNoPU );

    // std::vector<const reco::Candidate*> caloJetUnskimmed = getCollections( iEvent, srcCaloJet_);
    // std::vector<const reco::Candidate*> pfJetUnskimmed   = getCollections( iEvent, srcPfJet_);
    /*
       edm::Handle<std::vector<double> > regionEtHandle;
       edm::Handle<std::vector<Int_t> > regionEtaHandle;
       edm::Handle<std::vector<Int_t> > regionPhiHandle;

       iEvent.getByLabel(srcRegionEt_,regionEtHandle);
       iEvent.getByLabel(srcRegionEta_,regionEtaHandle);
       iEvent.getByLabel(srcRegionPhi_,regionPhiHandle);

       std::vector<double> regionEtTemp =  (*regionEtHandle); 
       std::vector<Int_t> regionEtaTemp = (*regionEtaHandle); 
       std::vector<Int_t> regionPhiTemp = (*regionPhiHandle); 


       std::vector<regionEtStruct> regionEts;
       for  (unsigned int i = 0; i < regionEtTemp.size();i++)
       {
       regionEtStruct temp;
       temp.et = regionEtTemp.at(i);
       double phi = (regionPhiTemp.at(i) > 9) ? (regionPhiTemp.at(i)*20.-420.):regionPhiTemp.at(i)*20.;  
       temp.phi = (phi * 3.1415927)/180.;
       double eta = (regionEtaTemp.at(i)/2.1)-5;
       temp.eta = eta;
       regionEts.push_back(temp);
       }

       std::sort(regionEts.begin(),regionEts.end(),compareByEt);

       regionEtTemp.clear(); 
       regionEtaTemp.clear();
       regionPhiTemp.clear();

       for (unsigned int i = 0; i < regionEts.size(); i++)
       {
       regionEtTemp.push_back(regionEts.at(i).et);
       regionEtaTemp.push_back(regionEts.at(i).eta);
       regionPhiTemp.push_back(regionEts.at(i).phi);
       }
       regionEt_ = &(regionEtTemp);
       regionEta_ = &(regionEtaTemp);
       regionPhi_ = &(regionPhiTemp);
       */
    // jetVeto30_ = false;
    // jetVeto50_ = false;
    // for (size_t i = 0; i < pfJetUnskimmed.size(); ++i) {
    // 	if(pfJetUnskimmed[i]->pt() >= 30 && abs(pfJetUnskimmed[i]->eta()) > 3) {jetVeto30_ = true;}
    // 	if(pfJetUnskimmed[i]->pt() >= 50 && abs(pfJetUnskimmed[i]->eta()) > 3) {jetVeto50_ = true;}
    // 	if (pfJetUnskimmed[i]->pt() < 30 ) {break;}
    // }


    // Skim jet collections
    // ----------------------------------------

    // Central jets
    // --------------------
    std::vector<const reco::Candidate*> gctCen     = skimJets(gctCenUnskimmed,     minPt, minEtaCen, maxEtaCen );

    //    std::cout << gctCen.size() << "\n";
    //    std::cout << "UCTHLT:\n";
    // for ( std::vector<const reco::Candidate*>::const_iterator itr = gctCen.begin(); itr != gctCen.end(); ++itr ){
    //   std::cout << (*itr)->pt() << "\t" <<  (*itr)->eta() << "\t" << (*itr)->phi() << "\n";
    // }
    //std::cout << "\n";

    //    std::vector<const reco::Candidate*> uctCen  = skimJets(uctCenUnskimmed,  minPt, minEtaCen, maxEtaCen );

    std::vector<const reco::Candidate*> genAk4     = skimJets(genJet4Unskimmed,    minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> genAk5     = skimJets(genJet5Unskimmed,    minPt, minEtaCen, maxEtaCen );

    std::vector<const reco::Candidate*> hltAk4Calo = skimJets(hltAk4CaloUnskimmed, minPt, minEtaCen, maxEtaCen );
    std::vector<const reco::Candidate*> hltAk4PF   = skimJets(hltAk4PFUnskimmed,   minPt, minEtaCen, maxEtaCen );
    //    std::vector<const reco::Candidate*> hltAk4PFNoPU = skimJets(hltAk4PFNoPUUnskimmed, minPt, maxEta );
    // std::vector<const reco::Candidate*> caloJet  = skimJets(caloJetUnskimmed,  minPt, maxEta );
    // std::vector<const reco::Candidate*> pfJet    = skimJets(pfJetUnskimmed,    minPt, maxEta );


    // Forward jets
    // --------------------
    std::vector<const reco::Candidate*> gctFor        = skimJets(gctForUnskimmed,     minPt, minEtaFor, maxEtaFor );
    //    std::vector<const reco::Candidate*> uctCenFor  = skimJets(uctForUnskimmed,  minPt, minEtaFor, maxEtaFor );

    std::vector<const reco::Candidate*> genAk4For     = skimJets(genJet4ForUnskimmed,    minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> genAk5For     = skimJets(genJet5ForUnskimmed,    minPt, minEtaFor, maxEtaFor );

    std::vector<const reco::Candidate*> hltAk4CaloFor = skimJets(hltAk4CaloForUnskimmed, minPt, minEtaFor, maxEtaFor );
    std::vector<const reco::Candidate*> hltAk4PFFor   = skimJets(hltAk4PFForUnskimmed,   minPt, minEtaFor, maxEtaFor );



    // Clear previous event's objects
    // --------------------------------------------------------------------------------

    genElectronPt->clear();
    genElectronEta->clear();
    genElectronPhi->clear();
    genMuonPt ->clear();
    genMuonEta->clear();
    genMuonPhi->clear();


    
    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
	jetPt[*iLvl] ->clear(); 
	jetPx[*iLvl] ->clear(); 
	jetPy[*iLvl] ->clear(); 
	jetPhi[*iLvl]->clear();
	jetEta[*iLvl]->clear();
	// if(*iLvl=="Uct" || *iLvl=="Gct" )
	// {
	//     jetPtsAll_[*iLvl]->clear();
	//     jetPhisAll_[*iLvl]->clear();
	//     jetEtasAll_[*iLvl]->clear();
	// }
    }

    // // Setup meta info
    // run_ = iEvent.id().run();
    // lumi_ = iEvent.id().luminosityBlock();
    // event_ = iEvent.id().event();

    // ********************************************************************************
    // *                           Loop over Jet collections                          *
    // ********************************************************************************

    genAk4AlphaTHT40     = calculateAlphaTHT( genAk4,     40.);
    hltAk4CaloAlphaTHT40 = calculateAlphaTHT( hltAk4Calo, 40.);
    hltAk4PFAlphaTHT40   = calculateAlphaTHT( hltAk4PF,   40.);

    genAk4AlphaTHT50     = calculateAlphaTHT( genAk4,     50.);
    hltAk4CaloAlphaTHT50 = calculateAlphaTHT( hltAk4Calo, 50.);
    hltAk4PFAlphaTHT50   = calculateAlphaTHT( hltAk4PF,   50.);

    //    std::cout << alphaTHT.first << "\t" << alphaTHT.second << "\n";

    // Store jet collections
    // ----------------------------------------

    storeJet( "gctCen", gctCen, jetPt, jetPx, jetPy, jetEta, jetPhi );
    //    storeJet( "uctCen", uctCen, jetPt, jetPx, jetPy, jetEta, jetPhi );

    storeJet( "genAk4", genAk4, jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk5", genAk5, jetPt, jetPx, jetPy, jetEta, jetPhi );

    storeJet( "hltAk4Calo",   hltAk4Calo,   jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4PF",     hltAk4PF,     jetPt, jetPx, jetPy, jetEta, jetPhi );
    //    storeJet( "hltAk4PFNoPU", hltAk4PFNoPU, jetPt, jetPx, jetPy, jetEta, jetPhi );


    storeJet( "gctFor",        gctFor,        jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk4For",     genAk4For,     jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "genAk5For",     genAk5For,     jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4CaloFor", hltAk4CaloFor, jetPt, jetPx, jetPy, jetEta, jetPhi );
    storeJet( "hltAk4PFFor",   hltAk4PFFor,   jetPt, jetPx, jetPy, jetEta, jetPhi );



    // Energy sums
    // --------------------
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

    getValue(iEvent, srcGenMetCalo_,             metPt_["genMetCalo"],             metPhi_["genMetCalo"]);
    getValue(iEvent, srcGenMetCaloAndNonPrompt_, metPt_["genMetCaloAndNonPrompt"], metPhi_["genMetCaloAndNonPrompt"]);
    getValue(iEvent, srcGenMetTrue_,             metPt_["genMetTrue"],             metPhi_["genMetTrue"]);



    // Gen particles
    // --------------------------------------------------------------------------------

    genLeptonVeto = false;
    if (makeGenParticles){
      edm::Handle< std::vector<reco::GenParticle> > genParticles;
      iEvent.getByLabel(srcGenParticles_, genParticles);
      for(std::vector<reco::GenParticle>::const_iterator iter = genParticles->begin(); iter != genParticles->end(); ++iter){
    	const reco::GenParticle& gen = *iter;

    	double genPt  = gen.p4().pt();
    	double genEta = gen.p4().eta();
    	double genPhi = gen.p4().phi();

    	// Electron
    	if(TMath::Abs(gen.pdgId()) == 11){
    	  if ( (genPt >= genElectronMinPt) && (TMath::Abs(genEta) <= genElectronMaxEta) ){
    	    genElectronPt ->push_back(  genPt );
    	    genElectronEta->push_back( genEta );
    	    genElectronPhi->push_back( genPhi );
    	    genLeptonVeto = true;
    	  }
    	}
    	// Muon
    	if(TMath::Abs(gen.pdgId()) == 13){
    	  if ( (genPt >= genMuonMinPt) && (TMath::Abs(genEta) <= genMuonMaxEta) ){
    	    genMuonPt ->push_back(  genPt );
    	    genMuonEta->push_back( genEta );
    	    genMuonPhi->push_back( genPhi );
    	    genLeptonVeto = true;
    	  }
    	}

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


 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakeTrees);
