/*
 * =====================================================================================
 *
 *       Filename:  MakeTrees.cc
 *
 *    Description:  Produces a flat tree for trigger studies
 *
 *        Initial Author:  Evan Friis, evan.friis@cern.ch
 *        Company:  UW Madison
 *
 *        Modified by: Adam Elwood, adam.elwood09@imperial.ac.uk
 *
 * =====================================================================================
 */

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticle.h"
#include "DataFormats/L1Trigger/interface/L1EtMissParticleFwd.h"

#include "DataFormats/L1Trigger/interface/Jet.h"
#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "TTree.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/GenJet.h"




struct alphatSums{
  double ht;
  double phi;
  double mhtx;
  double mhty;
  double mht;
  double dht;
  double alphat;
};
struct regionEtStruct{
double et;
double phi;
double eta;
};
bool compareByEt(const regionEtStruct &a, const regionEtStruct &b)
{
return a.et > b.et;
}

typedef std::vector<edm::InputTag> VInputTag;

class MakeTrees : public edm::EDAnalyzer {
  public:
    MakeTrees(const edm::ParameterSet& pset);
    virtual ~MakeTrees();
    void analyze(const edm::Event& evt, const edm::EventSetup& es);


    void storeJet( TString jetCollName, const std::vector<const reco::Candidate*>& jetColl,                                                         
                   std::map<TString,std::vector<Float_t>* >& pt,
		   std::map<TString,std::vector<Float_t>* >& eta,
		   std::map<TString,std::vector<Float_t>* >& phi);

		   //	 std::vector<Float_t>* pt, std::vector<Float_t>* eta, std::vector<Float_t>* phi);
    std::vector< const reco::Candidate*> skimJets(const std::vector<const reco::Candidate*>& inJets, double minPt, double maxEta);


  private:
    //Get the UCT stuff
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
    VInputTag srcGctJetAll_;

    //Get the Gen stuff
    VInputTag srcGen4Jet_;
    VInputTag srcGen5Jet_;
    VInputTag srcHLTAk4Calo;
    VInputTag srcHLTAk4PF;
    VInputTag srcHLTAk4PFNoPU;


    VInputTag srcCaloJet_;
    VInputTag srcPfJet_;




//Get the Regions;
/*    edm::InputTag srcRegionEt_;
    edm::InputTag srcRegionEta_;
    edm::InputTag srcRegionPhi_;
*/
    //Define the levels
    std::vector<TString> lvl_;
    TTree* tree;

    std::map<TString,Float_t> mhtPt30_;
    std::map<TString,Float_t> alphaT30_;
    std::map<TString,Float_t> mhtPhi30_;
    std::map<TString,Float_t> ht30_;
    std::map<TString,Float_t> dht30_;
    std::map<TString,Float_t> mhtDivHt30_;

    std::map<TString,Float_t> mhtPt_;
    std::map<TString,Float_t> alphaT_;
    std::map<TString,Float_t> mhtPhi_;
    std::map<TString,Float_t> ht_;
    std::map<TString,Float_t> dht_;
    std::map<TString,Float_t> mhtDivHt_;

    std::map<TString,Float_t> metPt_;
    std::map<TString,Float_t> metPhi_;
    std::map<TString,Float_t> et_;

    // std::map<TString,int> multiplicity_;
    // std::map<TString,int> multiplicity30_;
    // std::map<TString,int> multiplicity50_; 
  //    std::vector<Float_t> * invEnergy_;
    std::map<TString,std::vector<Float_t>* > jetPt;
    std::map<TString,std::vector<Float_t>* > jetMuons_;
    std::map<TString,std::vector<Float_t>* > jetPhi;
    std::map<TString,std::vector<Float_t>* > jetEta;
    std::map<TString,std::vector<Float_t>* > jetPtsAll_;
    std::map<TString,std::vector<Float_t>* > jetPhisAll_;
    std::map<TString,std::vector<Float_t>* > jetEtasAll_;
    std::map<TString,Float_t> jetDPhi12_;
    UInt_t Muon_;
    UInt_t NVTX;

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
    double maxEta;
    /*std::vector<double> * regionEt_;
    std::vector<Int_t> * regionEta_;
    std::vector<Int_t> * regionPhi_;*/
};

MakeTrees::MakeTrees(const edm::ParameterSet& pset) {

    // Initialize the ntuple builder
    edm::Service<TFileService> fs;
    tree = fs->make<TTree>("Ntuple", "Ntuple");

    // lvl_.push_back("Uct");
    // lvl_.push_back("S2Nopus");
    // lvl_.push_back("S2Donut");
    // lvl_.push_back("S2Global");
    // lvl_.push_back("Gct");

    //    lvl_.push_back("uctCen");
    lvl_.push_back("gctCen");

    lvl_.push_back("genAk4");
    lvl_.push_back("genAk5");

    lvl_.push_back("hltAk4Calo");
    lvl_.push_back("hltAk4PF");
    lvl_.push_back("hltAk4PFNoPU");


    // lvl_.push_back("Calo");
    // lvl_.push_back("Pf");
    //    invEnergy_ = new std::vector<Float_t>;
    //regionEt_ = new std::vector<double>;
    //regionEta_ = new std::vector<Int_t>;
    //regionPhi_ = new std::vector<Int_t>;


    //    NVTX = 0;
    //    tree->Branch("NVTX", &NVTX, "NVTX/i");


    //Set the tree branches
    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
	mhtPt30_[*iLvl]        = -10.0;
	mhtPhi30_[*iLvl]       = -10.0;
	ht30_[*iLvl]           = -10.0;
	dht30_[*iLvl]           = -10.0;
	alphaT30_[*iLvl]       = 0.0;
	mhtDivHt30_[*iLvl]     = 0.0;

	mhtPt_[*iLvl]        = -10.0;
	mhtPhi_[*iLvl]       = -10.0;
	ht_[*iLvl]           = -10.0;
	dht_[*iLvl]           = -10.0;
	alphaT_[*iLvl]       = 0.0;
	mhtDivHt_[*iLvl]     = 0.0;

	metPt_[*iLvl]        = -10.0;
	metPhi_[*iLvl]       = -10.0;
	et_[*iLvl]       = -10.0;
	// multiplicity_[*iLvl] = -1;
	// multiplicity30_[*iLvl] = -1;
	// multiplicity50_[*iLvl] = -1;
	jetPt[*iLvl]       = new std::vector<Float_t>();
	jetPhi[*iLvl]      = new std::vector<Float_t>();
	jetEta[*iLvl]      = new std::vector<Float_t>();
	//	jetDPhi12_[*iLvl]    = -1.0;
	// if(*iLvl=="Uct" || *iLvl=="Gct"){
	//     jetPtsAll_[*iLvl]  = new std::vector<Float_t>();
	//     jetPhisAll_[*iLvl] = new std::vector<Float_t>();
	//     jetEtasAll_[*iLvl] = new std::vector<Float_t>();
	// }

	// if(*iLvl=="Gen4" || *iLvl=="Gen5" || *iLvl=="Calo" || *iLvl=="Pf"){ 
	//     tree->Branch("mhtPt30"+*iLvl, &mhtPt30_[*iLvl], "mhtPt30"+*iLvl+"/f");
	//     tree->Branch("mhtPhi30"+*iLvl, &mhtPhi30_[*iLvl], "mhtPhi30"+*iLvl+"/f");
	//     tree->Branch("ht30"+*iLvl, &ht30_[*iLvl], "ht30"+*iLvl+"/f");
	//     tree->Branch("dht30"+*iLvl, &dht30_[*iLvl], "dht30"+*iLvl+"/f");
	//     tree->Branch("alphaT30"+*iLvl, &alphaT30_[*iLvl], "alphaT30"+*iLvl+"/f");
	//     tree->Branch("mhtDivHt30"+*iLvl, &mhtDivHt30_[*iLvl], "mhtDivHt30"+*iLvl+"/f");
	// tree->Branch("multiplicity30"+*iLvl, &multiplicity30_[*iLvl], "multiplicity30"+*iLvl+"/i");
	// tree->Branch("multiplicity50"+*iLvl, &multiplicity50_[*iLvl], "multiplicity50"+*iLvl+"/i");
	// }
	// tree->Branch("mhtPt"+*iLvl, &mhtPt_[*iLvl], "mhtPt"+*iLvl+"/f");
	// tree->Branch("mhtPhi"+*iLvl, &mhtPhi_[*iLvl], "mhtPhi"+*iLvl+"/f");
	// tree->Branch("ht"+*iLvl, &ht_[*iLvl], "ht"+*iLvl+"/f");
	// tree->Branch("dht"+*iLvl, &dht_[*iLvl], "dht"+*iLvl+"/f");
	// tree->Branch("alphaT"+*iLvl, &alphaT_[*iLvl], "alphaT"+*iLvl+"/f");
	// tree->Branch("mhtDivHt"+*iLvl, &mhtDivHt_[*iLvl], "mhtDivHt"+*iLvl+"/f");

	// if(*iLvl!="S2Global" && *iLvl!="S2Nopus" && *iLvl!="S2Donut"){ 
	// tree->Branch("metPt"+*iLvl, &metPt_[*iLvl], "metPt"+*iLvl+"/f");
	// tree->Branch("metPhi"+*iLvl, &metPhi_[*iLvl], "metPhi"+*iLvl+"/f");
	// tree->Branch("et"+*iLvl, &et_[*iLvl], "et"+*iLvl+"/f");
	// }

	// tree->Branch("multiplicity"+*iLvl, &multiplicity_[*iLvl], "multiplicity"+*iLvl+"/i");
	// tree->Branch("jetDPhi12"+*iLvl, &jetDPhi12_[*iLvl], "jetDPhi12"+*iLvl+"/f");

	tree->Branch(*iLvl + "_Pt",  "std::vector<float>", &jetPt[*iLvl]);
	tree->Branch(*iLvl + "_Phi", "std::vector<float>", &jetPhi[*iLvl]);
	tree->Branch(*iLvl + "_Eta", "std::vector<float>", &jetEta[*iLvl]);

	// if(*iLvl=="Uct" || *iLvl=="Gct"){
	//     tree->Branch("jetPtsAll"+*iLvl, "std::vector<float>", &jetPtsAll_[*iLvl]);
	//     tree->Branch("jetPhisAll"+*iLvl, "std::vector<float>", &jetPhisAll_[*iLvl]);
	//     tree->Branch("jetEtasAll"+*iLvl, "std::vector<float>", &jetEtasAll_[*iLvl]);
	// }
    }
    //    tree->Branch("invEnergy", "std::vector<float>", &invEnergy_);

    //tree->Branch("regionEt", "std::vector<double>", &regionEt_);
    //tree->Branch("regionEta", "std::vector<int>", &regionEta_);
    //tree->Branch("regionPhi", "std::vector<int>", &regionPhi_);

    // tree->Branch("metPtS2", &metPt_["S2"], "metPtS2/f");
    // tree->Branch("metPhiS2", &metPhi_["S2"], "metPhiS2/f");
    // tree->Branch("etS2", &et_["S2"], "etS2/f");

    // tree->Branch("jetVeto50",   &jetVeto50_,   "jetVeto50/O");
    // tree->Branch("jetVeto30",   &jetVeto30_,   "jetVeto30/O");

    // jetMuons_["Pf"]       = new std::vector<Float_t>();
    // tree->Branch("jetMuonsPf", "std::vector<float>", &jetMuons_["Pf"]);

    //    tree->Branch("MuonNumber",   &Muon_,   "NMuons/i");
    // tree->Branch("run",   &run_,   "run/i");
    // tree->Branch("lumi",  &lumi_,  "run/i");
    // tree->Branch("event", &event_, "run/l");


    // srcUctMht_        = pset.getParameter<edm::InputTag>("srcUctMht");
    srcUctJet_ = pset.getParameter<edm::InputTag>("srcUctJet");

    // srcS2Met_        = pset.getParameter<edm::InputTag>("srcS2Met");
    // srcS2DonutMht_        = pset.getParameter<edm::InputTag>("srcS2DonutMht");
    // srcS2DonutJetCentral_ = pset.getParameter<VInputTag>("srcS2DonutJetCentral");
    // srcS2NopusMht_        = pset.getParameter<edm::InputTag>("srcS2NopusMht");
    // srcS2NopusJetCentral_ = pset.getParameter<VInputTag>("srcS2NopusJetCentral");
    // srcS2GlobalMht_        = pset.getParameter<edm::InputTag>("srcS2GlobalMht");
    // srcS2GlobalJetCentral_ = pset.getParameter<VInputTag>("srcS2GlobalJetCentral");

    // srcGctMht_        = pset.getParameter<edm::InputTag>("srcGctMht");
    // srcGctMet_        = pset.getParameter<edm::InputTag>("srcGctMet");
    srcGctJetCentral_ = pset.getParameter<VInputTag>("srcGctJetCentral");
    // srcGctJetAll_     = pset.getParameter<VInputTag>("srcGctJetAll");

    srcGen4Jet_        = pset.getParameter<VInputTag>("srcGen4Jet");
    srcGen5Jet_        = pset.getParameter<VInputTag>("srcGen5Jet");

    srcHLTAk4Calo      = pset.getParameter<VInputTag>("srcHLTAk4Calo");
    srcHLTAk4PF        = pset.getParameter<VInputTag>("srcHLTAk4PF");
    srcHLTAk4PFNoPU    = pset.getParameter<VInputTag>("srcHLTAk4PFNoPU");


    srcCaloJet_       = pset.getParameter<VInputTag>("srcCaloJet");
    srcPfJet_         = pset.getParameter<VInputTag>("srcPfJet");

    //Parameters for the HSums
    //htThreshold_      = pset.getParameter<double>("htThreshold");
    maxjet_           = pset.getParameter<unsigned int>("maxjet");
    usePU_            = pset.getParameter<bool>("usePU");

    // Jet skim cuts
    minPt             = pset.getParameter<double>("jetMinPt");
    maxEta            = pset.getParameter<double>("jetMaxEta");
/*
    srcRegionEt_ = pset.getParameter<edm::InputTag>("srcRegionEt");
    srcRegionEta_ = pset.getParameter<edm::InputTag>("srcRegionEta");
    srcRegionPhi_ = pset.getParameter<edm::InputTag>("srcRegionPhi");
*/

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
    std::vector<const reco::Candidate*> getCollections(
	    const edm::Event& evt, const VInputTag& collections) {
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

    //Return the jet energy sums in a vector: ht, mhtPt, mhtPhi, mhtPx, mhtPy
    alphatSums getHSum(const std::vector<const reco::Candidate*>& jets, double htThreshold,unsigned int maxJet){

	double ht   = 0.0;
	double mhtX = 0.0;
	double mhtY = 0.0;

	int jetsAboveThresh(0);
	for (size_t i = 0; i < jets.size(); ++i) {

	    if ( jets[i]->pt() > htThreshold){

		jetsAboveThresh++;
		ht   += jets[i]->pt();
		mhtX += jets[i]->px();
		mhtY += jets[i]->py();
	    }
	    else{break;}
	}
	double mht = sqrt(mhtX*mhtX + mhtY*mhtY);
	alphatSums hSumVec;
	//
	double phi = atan2(mhtY,mhtX) + PI;

	//Make sure Phi is in range -pi to pi
	phi = (phi>PI) ? phi-2.0*PI : phi;

	// Minimum Delta Et for two pseudo-jets
	//  if ( (jets.size() > 1) && jets.size() <= maxjet) // bit harsh - removes events with large jet multiplicities
	double alphaT = 0.;
	double min_delta_sum_et = -1.;
	if (jetsAboveThresh > 1)
	{

	    unsigned int jetLimit = jetsAboveThresh;
	    if ( jetLimit > maxJet){ jetLimit = maxJet; }

	    //      for ( unsigned i=0; i < unsigned(1<<(jets.size()-1)); i++ )  //@@ iterate through different combinations
	    for ( unsigned i=0; i < unsigned(1<<( jetLimit - 1)); i++ ) { //@@ iterate through different combinations
		double delta_sum_et = 0.;
		for ( unsigned j=0; j < jetLimit; j++ ) { //@@ iterate through jets
		    delta_sum_et += jets.at(j)->pt() * ( 1 - 2 * (int(i>>j)&1) ); 
		}
		if ( ( fabs(delta_sum_et) < min_delta_sum_et || min_delta_sum_et < 0. ) ) {
		    min_delta_sum_et = fabs(delta_sum_et);
		}
	    }
	    if ( min_delta_sum_et < 0. ) { alphaT=0.; }

	    // Alpha_T
	    alphaT = 0.5 * ( ht - min_delta_sum_et ) / sqrt( ht*ht - (mht*mht) );
	    if (alphaT > 10){ alphaT = 10; }
	    //    return ( 0.5 * ( sum_et - min_delta_sum_et ) / sqrt( sum_et*sum_et - (sum_met*sum_met) ) );
	}
	else
	{
	    alphaT = 0.;
	}

	hSumVec.ht=ht;
	hSumVec.dht=min_delta_sum_et;
	hSumVec.mht = mht;
	hSumVec.phi =phi;
	hSumVec.mhtx = mhtX;
	hSumVec.mhty = mhtY;
	hSumVec.alphat = alphaT;

	return hSumVec;
    }


    inline float deltaPhi( float phi1, float phi2 ){

	float const  PI        = ROOT::Math::Pi();
	float const  TWOPI     = 2.*PI;
	float dPhi = (phi1 - phi2);
	while (dPhi >= PI) dPhi -= TWOPI;
	while (dPhi < -PI) dPhi += TWOPI;
	return dPhi;

    }
}


void MakeTrees::analyze(const edm::Event& evt, const edm::EventSetup& es) {


    // // get nvertices from simulation
    // edm::Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    // evt.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
    // if(usePU_){
    //   NVTX = 0;

    // 	std::vector<PileupSummaryInfo>::const_iterator PVI;
    // 	for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    // 	    int BX = PVI->getBunchCrossing();
    // 	    //	    std::cout << "bx = " << BX << std::endl;
    // 	    if(BX == 0) {
    // 		NVTX = PVI->getPU_NumInteractions();
    // 		//		std::cout << "NPV = " << NVTX << std::endl;
    // 		break;
    // 	    }
    // 	}
    // }


  //   // Get UCT jets
  // //    edm::Handle< BXVector<l1t::Jet>> uctCalibJets;
  // edm::Handle< std::vector<L1GctJetCand> > uctCalibJets;
  //   evt.getByLabel(srcUctJet_, uctCalibJets);

  //   std::vector<const reco::Candidate*> uctCenUnskimmed;
  //   // std::vector<const reco::Candidate*> uctJetAll;

  //   for ( std::vector<L1GctJetCand>::const_iterator itr = uctCalibJets->begin(); itr != uctCalibJets->end(); ++itr ){
  //     //    for ( auto itr = uctCalibJets->begin(0); itr != uctCalibJets->end(0); ++itr ) {
  //     //    	if(fabs(itr->eta())<=3.0){
  // 	  std::cout << (itr)->etaIndex() << "\n";
  // 	  //    	    uctCenUnskimmed.push_back(&(*itr));
  // 	  //    	}
  // 	//    	uctJetAll.push_back(&(*itr));
  //   }

    
    // Input jets without eta or pT requirements
    // --------------------------------------------------------------------------------
    // std::vector<const reco::Candidate*> s2NopusJetCentral    = getCollections( evt, srcS2NopusJetCentral_);
    // std::vector<const reco::Candidate*> s2GlobalJetCentral    = getCollections( evt, srcS2GlobalJetCentral_);
    // std::vector<const reco::Candidate*> s2DonutJetCentral    = getCollections( evt, srcS2DonutJetCentral_);
    std::vector<const reco::Candidate*> gctCenUnskimmed       = getCollections( evt, srcGctJetCentral_);
    // std::vector<const reco::Candidate*> gctJetAll        = getCollections( evt, srcGctJetAll_);
    std::vector<const reco::Candidate*> genJet4Unskimmed      = getCollections( evt, srcGen4Jet_);
    std::vector<const reco::Candidate*> genJet5Unskimmed      = getCollections( evt, srcGen5Jet_);

    std::vector<const reco::Candidate*> hltAk4CaloUnskimmed   = getCollections( evt, srcHLTAk4Calo );
    std::vector<const reco::Candidate*> hltAk4PFUnskimmed     = getCollections( evt, srcHLTAk4PF );
    std::vector<const reco::Candidate*> hltAk4PFNoPUUnskimmed = getCollections( evt, srcHLTAk4PFNoPU );

    // std::vector<const reco::Candidate*> caloJetUnskimmed = getCollections( evt, srcCaloJet_);
    // std::vector<const reco::Candidate*> pfJetUnskimmed   = getCollections( evt, srcPfJet_);
    /*
       edm::Handle<std::vector<double> > regionEtHandle;
       edm::Handle<std::vector<Int_t> > regionEtaHandle;
       edm::Handle<std::vector<Int_t> > regionPhiHandle;

       evt.getByLabel(srcRegionEt_,regionEtHandle);
       evt.getByLabel(srcRegionEta_,regionEtaHandle);
       evt.getByLabel(srcRegionPhi_,regionPhiHandle);

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
    std::vector<const reco::Candidate*> gctCen  = skimJets(gctCenUnskimmed,   minPt, maxEta );
    //    std::vector<const reco::Candidate*> uctCen  = skimJets(uctCenUnskimmed,     minPt, maxEta ); 

    std::vector<const reco::Candidate*> genAk4  = skimJets(genJet4Unskimmed,  minPt, maxEta );
    std::vector<const reco::Candidate*> genAk5  = skimJets(genJet5Unskimmed,  minPt, maxEta );

    std::vector<const reco::Candidate*> hltAk4Calo   = skimJets(hltAk4CaloUnskimmed,   minPt, maxEta );
    std::vector<const reco::Candidate*> hltAk4PF     = skimJets(hltAk4PFUnskimmed,     minPt, maxEta );
    std::vector<const reco::Candidate*> hltAk4PFNoPU = skimJets(hltAk4PFNoPUUnskimmed, minPt, maxEta );


    // std::vector<const reco::Candidate*> caloJet  = skimJets(caloJetUnskimmed,  minPt, maxEta );
    // std::vector<const reco::Candidate*> pfJet    = skimJets(pfJetUnskimmed,    minPt, maxEta );
    // Clear previous event's objects
    //    jetMuons_["Pf"]->clear();
    for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
	jetPt[*iLvl] ->clear(); 
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
    // run_ = evt.id().run();
    // lumi_ = evt.id().luminosityBlock();
    // event_ = evt.id().event();

    // ********************************************************************************
    // *                           Loop over Jet collections                          *
    // ********************************************************************************
    //First L1 Collections
    // for (size_t i = 0; i < uctJetCentral.size(); ++i) {
    // 	jetPts_["Uct"]->push_back(uctJetCentral[i]->pt());
    // 	jetEtas_["Uct"]->push_back(uctJetCentral[i]->eta());
    // 	jetPhis_["Uct"]->push_back(uctJetCentral[i]->phi());
    // }
    //if(uctJetCentral.size()>1) jetDPhi12_["Uct"] = deltaPhi(uctJetCentral[0]->phi(), uctJetCentral[1]->phi());

    // for (size_t i = 0; i < s2NopusJetCentral.size(); ++i) {
    // 	jetPts_["S2Nopus"]->push_back(s2NopusJetCentral[i]->pt());
    // 	jetEtas_["S2Nopus"]->push_back(s2NopusJetCentral[i]->eta());
    // 	jetPhis_["S2Nopus"]->push_back(s2NopusJetCentral[i]->phi());
    // }
    // for (size_t i = 0; i < s2GlobalJetCentral.size(); ++i) {
    // 	jetPts_["S2Global"]->push_back(s2GlobalJetCentral[i]->pt());
    // 	jetEtas_["S2Global"]->push_back(s2GlobalJetCentral[i]->eta());
    // 	jetPhis_["S2Global"]->push_back(s2GlobalJetCentral[i]->phi());
    // }
    // for (size_t i = 0; i < s2DonutJetCentral.size(); ++i) {
    // 	jetPts_["S2Donut"]->push_back(s2DonutJetCentral[i]->pt());
    // 	jetEtas_["S2Donut"]->push_back(s2DonutJetCentral[i]->eta());
    // 	jetPhis_["S2Donut"]->push_back(s2DonutJetCentral[i]->phi());
    // }


    // for (size_t i = 0; i < uctJetAll.size(); ++i) {
    // 	jetPtsAll_["Uct"]->push_back(uctJetAll[i]->pt());
    // 	jetEtasAll_["Uct"]->push_back(uctJetAll[i]->eta());
    // 	jetPhisAll_["Uct"]->push_back(uctJetAll[i]->phi());
    // }
    // for (size_t i = 0; i < gctJetCentral.size(); ++i) {
    // 	jetPts_["Gct"]->push_back(gctJetCentral[i]->pt());
    // 	jetEtas_["Gct"]->push_back(gctJetCentral[i]->eta());
    // 	jetPhis_["Gct"]->push_back(gctJetCentral[i]->phi());
    // }
    // if(gctJetCentral.size()>1) jetDPhi12_["Gct"] = deltaPhi(gctJetCentral[0]->phi(), gctJetCentral[1]->phi());

    // for (size_t i = 0; i < gctJetAll.size(); ++i) {
    // 	jetPtsAll_["Gct"]->push_back(gctJetAll[i]->pt());
    // 	jetEtasAll_["Gct"]->push_back(gctJetAll[i]->eta());
    // 	jetPhisAll_["Gct"]->push_back(gctJetAll[i]->phi());
    // }
    //Now Gen Collections
    //    invEnergy_->clear();


    // Store jet collections
    // ----------------------------------------

    storeJet( "gctCen", gctCen, jetPt, jetEta, jetPhi );
    //    storeJet( "uctCen", uctCen, jetPt, jetEta, jetPhi );


    storeJet( "genAk4", genAk4, jetPt, jetEta, jetPhi );
    storeJet( "genAk5", genAk5, jetPt, jetEta, jetPhi );

    storeJet( "hltAk4Calo",   hltAk4Calo,   jetPt, jetEta, jetPhi );
    storeJet( "hltAk4PF",     hltAk4PF,     jetPt, jetEta, jetPhi );
    storeJet( "hltAk4PFNoPU", hltAk4PFNoPU, jetPt, jetEta, jetPhi );




    // for (size_t i = 0; i < genJet4.size(); ++i) {
    // 	jetPts_["genak4"] ->push_back(genJet4[i]->pt());
    // 	jetEtas_["genak4"]->push_back(genJet4[i]->eta());
    // 	jetPhis_["genak4"]->push_back(genJet4[i]->phi());
    // 	invEnergy_->push_back((static_cast<const reco::GenJet *>(genJet4[i]))->invisibleEnergy());
    // }
    // //    if(genJet4.size()>1) jetDPhi12_["Gen4"] = deltaPhi(genJet4[0]->phi(), genJet4[1]->phi());
    // for (size_t i = 0; i < genJet5.size(); ++i) {
    // 	jetPts_["genak5"] ->push_back(genJet5[i]->pt());
    // 	jetEtas_["genak5"]->push_back(genJet5[i]->eta());
    // 	jetPhis_["genak5"]->push_back(genJet5[i]->phi());
    // }
    //    if(genJet5.size()>1) jetDPhi12_["genak5"] = deltaPhi(genJet5[0]->phi(), genJet5[1]->phi());

    //Now Reco Collections
    // for (size_t i = 0; i < caloJet.size(); ++i) {
    // 	jetPts_["Calo"] ->push_back(caloJet[i]->pt());
    // 	jetEtas_["Calo"]->push_back(caloJet[i]->eta());
    // 	jetPhis_["Calo"]->push_back(caloJet[i]->phi());
    // }
    // if(caloJet.size()>1) jetDPhi12_["Calo"] = deltaPhi(caloJet[0]->phi(), caloJet[1]->phi());
    // Muon_=0;
    // for (size_t i = 0; i < pfJet.size(); ++i) {
    // 	jetPts_["Pf"] ->push_back(pfJet[i]->pt());
    // 	jetEtas_["Pf"]->push_back(pfJet[i]->eta());
    // 	jetPhis_["Pf"]->push_back(pfJet[i]->phi());
    // 	jetMuons_["Pf"]->push_back((static_cast<const reco::PFJet *>(pfJet[i]))->muonEnergy());
    // 	Muon_ += (static_cast<const reco::PFJet *>(pfJet[i]))->muonMultiplicity();
    // }
    // if(pfJet.size()>1) jetDPhi12_["Pf"] = deltaPhi(pfJet[0]->phi(), pfJet[1]->phi());

    // //Set the Mht and ht sums
    // edm::Handle< BXVector<l1t::EtSum> > uctHtHandle;
    // evt.getByLabel(srcUctMht_,uctHtHandle);

    // for(auto itr = uctHtHandle->begin(0); itr != uctHtHandle->end(0); itr++){

    // 	if(itr->getType() == l1t::EtSum::EtSumType::kMissingHt){
    // 	    mhtPt_["Uct"]=itr->pt();
    // 	    mhtPhi_["Uct"]=itr->phi();
    // 	}
    // 	if(itr->getType() == l1t::EtSum::EtSumType::kMissingEt){
    // 	    metPt_["Uct"]=itr->pt();
    // 	    metPhi_["Uct"]=itr->phi();
    // 	}
    // 	if(itr->getType() == l1t::EtSum::EtSumType::kTotalHt){
    // 	    ht_["Uct"]=itr->pt();
    // 	}
    // 	if(itr->getType() == l1t::EtSum::EtSumType::kTotalEt){
    // 	    et_["Uct"]=itr->pt();
    // 	}
    // }






    //    getValue(evt, srcGctMht_,mhtPt_["Gct"],mhtPhi_["Gct"]);
    // getValue(evt, srcS2NopusMht_,mhtPt_["S2Nopus"],mhtPhi_["S2Nopus"]);
    // getValue(evt, srcS2GlobalMht_,mhtPt_["S2Global"],mhtPhi_["S2Global"]);
    // getValue(evt, srcS2DonutMht_,mhtPt_["S2Donut"],mhtPhi_["S2Donut"]);
    // getSumEtL1(evt, srcS2NopusMht_, ht_["S2Nopus"], false);
    // getSumEtL1(evt, srcS2GlobalMht_, ht_["S2Global"], false);
    // getSumEtL1(evt, srcS2DonutMht_, ht_["S2Donut"], false);
    // getSumEtL1(evt, srcGctMht_, ht_["Gct"], false);
    // getValue(evt, srcGctMet_,metPt_["Gct"],metPhi_["Gct"]);
    // getValue(evt, srcS2Met_,metPt_["S2"],metPhi_["S2"]);
    // getSumEtL1(evt, srcGctMet_, et_["Gct"], false);
    // getSumEtL1(evt, srcS2Met_, et_["S2"], false);

    // alphatSums hSumVec = getHSum(gctJetCentral, 50,3);
    // alphaT_["Gct"]       = hSumVec.alphat;
    // dht_["Gct"]       = hSumVec.dht;
    // multiplicity_["Gct"] = gctJetCentral.size();

    // hSumVec = getHSum(s2NopusJetCentral, 50,3);
    // alphaT_["S2Nopus"]       = hSumVec.alphat;
    // dht_["S2Nopus"]       = hSumVec.dht;
    // multiplicity_["S2Nopus"] = s2NopusJetCentral.size();

    // hSumVec = getHSum(s2GlobalJetCentral, 50,3);
    // alphaT_["S2Global"]       = hSumVec.alphat;
    // dht_["S2Global"]       = hSumVec.dht;
    // multiplicity_["S2Global"] = s2GlobalJetCentral.size();

    // hSumVec = getHSum(s2DonutJetCentral, 50,3);
    // alphaT_["S2Donut"]       = hSumVec.alphat;
    // dht_["S2Donut"]       = hSumVec.dht;
    // multiplicity_["S2Donut"] = s2DonutJetCentral.size();

    // hSumVec              = getHSum(uctJetCentral, 50,3);
    // alphaT_["Uct"]       = hSumVec.alphat;
    // dht_["Uct"]       = hSumVec.dht;
    // multiplicity_["Uct"] = uctJetCentral.size();

    // hSumVec = getHSum(genJet4, 50, maxjet_);
    // ht_["Gen4"]           = hSumVec.ht;
    // dht_["Gen4"]           = hSumVec.dht;
    // mhtPt_["Gen4"]        = hSumVec.mht;
    // mhtPhi_["Gen4"]       = hSumVec.phi;
    // alphaT_["Gen4"]       = hSumVec.alphat;


    //    unsigned jetCount = 0;
    // while (genJet4.size() != jetCount && genJet4.at(jetCount)->pt() > 50) jetCount++;
    // multiplicity50_["Gen4"] = jetCount;

    // hSumVec = getHSum(genJet4, 30, maxjet_);
    // ht30_["Gen4"]           = hSumVec.ht;
    // dht30_["Gen4"]           = hSumVec.dht;
    // mhtPt30_["Gen4"]        = hSumVec.mht;
    // mhtPhi30_["Gen4"]       = hSumVec.phi;
    // alphaT30_["Gen4"]       = hSumVec.alphat;
    // //std::cout << "713" << std::endl;
    // jetCount = 0;
    // while (genJet4.size() != jetCount && genJet4.at(jetCount)->pt() > 30) jetCount++;
    // multiplicity30_["Gen4"] = jetCount;
    // multiplicity_["Gen4"] = genJet4.size();

    // hSumVec = getHSum(genJet5, 50, maxjet_);
    // ht_["Gen5"]           = hSumVec.ht;
    // dht_["Gen5"]           = hSumVec.dht;
    // mhtPt_["Gen5"]        = hSumVec.mht;
    // mhtPhi_["Gen5"]       = hSumVec.phi;
    // alphaT_["Gen5"]       = hSumVec.alphat;

    // //std::cout << "726" << std::endl;
    // jetCount = 0;
    // while (genJet5.size() != jetCount && genJet5.at(jetCount)->pt() > 50) jetCount++;
    // multiplicity50_["Gen5"] = jetCount;

    // hSumVec = getHSum(genJet5, 30, maxjet_);
    // ht30_["Gen5"]           = hSumVec.ht;
    // dht30_["Gen5"]           = hSumVec.dht;
    // mhtPt30_["Gen5"]        = hSumVec.mht;
    // mhtPhi30_["Gen5"]       = hSumVec.phi;
    // alphaT30_["Gen5"]       = hSumVec.alphat;

    // jetCount = 0;
    // while (genJet5.size() != jetCount && genJet5.at(jetCount)->pt() > 30) jetCount++;
    // multiplicity30_["Gen5"] = jetCount;
    // multiplicity_["Gen5"] = genJet5.size();
    // //std::cout << "742" << std::endl;

    // //std::cout << alphaT_["Gen"] <<std::endl;
    // //hSumVec.clear();
    // hSumVec               = getHSum(caloJet, 50,maxjet_);
    // ht_["Calo"]           = hSumVec.ht;
    // dht_["Calo"]           = hSumVec.dht;
    // mhtPt_["Calo"]        = hSumVec.mht;
    // mhtPhi_["Calo"]       = hSumVec.phi;
    // alphaT_["Calo"]       = hSumVec.alphat;
    // multiplicity_["Calo"] = caloJet.size();

    // jetCount = 0;
    // while (caloJet.size() != jetCount && caloJet.at(jetCount)->pt() > 50) jetCount++;
    // multiplicity50_["Calo"] = jetCount;

    // hSumVec               = getHSum(caloJet, 30,maxjet_);
    // ht30_["Calo"]           = hSumVec.ht;
    // dht30_["Calo"]           = hSumVec.dht;
    // mhtPt30_["Calo"]        = hSumVec.mht;
    // mhtPhi30_["Calo"]       = hSumVec.phi;
    // alphaT30_["Calo"]       = hSumVec.alphat;

    // jetCount = 0;
    // while (caloJet.size() != jetCount && caloJet.at(jetCount)->pt() > 30) jetCount++;
    // multiplicity30_["Calo"] = jetCount;

    // //std::cout << "769" << std::endl;
    // //hSumVec.clear();
    // hSumVec             = getHSum(pfJet, 50,maxjet_);
    // ht_["Pf"]           = hSumVec.ht;
    // dht_["Pf"]           = hSumVec.dht;
    // mhtPt_["Pf"]        = hSumVec.mht;
    // mhtPhi_["Pf"]       = hSumVec.phi;
    // alphaT_["Pf"]       = hSumVec.alphat;
    // multiplicity_["Pf"] = caloJet.size();
    // jetCount = 0;
    // while (pfJet.size() != jetCount && pfJet.at(jetCount)->pt() > 50) jetCount++;
    // multiplicity50_["Pf"] = jetCount;

    // hSumVec             = getHSum(pfJet, 30,maxjet_);
    // ht30_["Pf"]           = hSumVec.ht;
    // dht30_["Pf"]           = hSumVec.dht;
    // mhtPt30_["Pf"]        = hSumVec.mht;
    // mhtPhi30_["Pf"]       = hSumVec.phi;
    // alphaT30_["Pf"]       = hSumVec.alphat;

    // jetCount = 0;
    // while (pfJet.size() != jetCount && pfJet.at(jetCount)->pt() > 30) jetCount++;
    // multiplicity30_["Pf"] = jetCount;
    // //multiplicity_["Pf"] = pfJet.size();

    // for(std::vector<TString>::const_iterator iLvl=lvl_.begin(); iLvl!=lvl_.end(); iLvl++){
    // 	if(ht_[*iLvl]>1.) mhtDivHt_[*iLvl]=mhtPt_[*iLvl]/ht_[*iLvl];
    // }

    // Fill the trees

    tree->Fill();
}






std::vector< const reco::Candidate*> 
MakeTrees::skimJets(const std::vector<const reco::Candidate*>& inputJets, double minPt, double maxEta){

    std::vector <const reco::Candidate*> skimmedJets; 
    for ( unsigned int iJet = 0; iJet < inputJets.size(); iJet++ ){
	if ( (inputJets.at(iJet)->pt() > minPt) && ( fabs(inputJets.at(iJet)->eta()) < maxEta) ){
	    skimmedJets.push_back(inputJets.at(iJet));
	}
    }

    return skimmedJets;
}



void
MakeTrees::storeJet( TString jetCollName, const std::vector<const reco::Candidate*>& jetColl,
		     std::map<TString,std::vector<Float_t>* >& pt,
		     std::map<TString,std::vector<Float_t>* >& eta,
		     std::map<TString,std::vector<Float_t>* >& phi){

  for ( std::vector<const reco::Candidate*>::const_iterator itr = jetColl.begin(); itr != jetColl.end(); ++itr ){

    pt [jetCollName]->push_back( (*itr)->pt() );
    eta[jetCollName]->push_back( (*itr)->eta() );
    phi[jetCollName]->push_back( (*itr)->phi() );
    
  }

}


 
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(MakeTrees);
