#define TEST
#define SIGNAL
//#define NEUTRINO
// Switch between PF and Calojets at HLT
#define HLT_CALOJET

const int N_50NS_BUNCHES = 1368;
const int N_25NS_BUNCHES = 2508;
const int N_MAX_BUNCHES  = 3564;
const double CROSSING_TIME  = 25E-9;
const int LHC_FREQUENCY  = 11246;
const int ZB_XSECTION    = LHC_FREQUENCY*N_25NS_BUNCHES;

const double INST_LUMI_25NS = 1.4E34;
const double PB_TO_CM2      = 1e-36;



#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TStyle.h>
#include <TCut.h>
#include <TMath.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>

#include <vector>
#include <iostream>
#include <math.h>

// Function prototypes
//TH2D * makeCumu(TH2D * input);
void reverseCumulative( TH1* histogram, TH1* rCumulHist, double scale);
void reverseCumulative2D( TH2* histogram, TH2* rCumulHist, double scale);
void fillUniform2D( TH2* histogram, double value);
void clearTriangle( TH1* histogram, int sign, float extra );
void clearRectangleX( TH1* histogram, float cut );

void fillUniformY( TH2* histogram, TH2* rCumulHist);
void reverseCumulativeY( TH2* histogram, TH2* rCumulHist, double scale);

std::vector< std::pair< int, int> > getEfficiencyBins( TH2* histogram, double maxRate );
bool passesAlphaTSelection( std::vector<float> *jetPT, int &jetsAbove50 );



float calculateAlphaT( std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold );
float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold, uint maxJets, 
			     float dynamicJetThreshold, float dynamicAlphaTThreshold);
std::pair<float,float> calculateDynamicAlphaTHT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold, uint maxJets,
                                                float dynamicJetThreshold, float dynamicAlphaTThreshold, float dynamicHTThreshold);

std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold, uint maxJets, float dynamicJetThreshold);


float const  PI        = TMath::Pi();
float const  TWOPI     = 2.*PI;
float const  RADTODEG  = 180/PI;

// Returns smallest angle between two azimuths in the range [-PI, PI]
inline float deltaPhi( float phi1, float phi2 ){

  float dPhi = (phi1 - phi2);
  while (dPhi >= PI) dPhi -= TWOPI;
  while (dPhi < -PI) dPhi += TWOPI;
  return dPhi;

}


inline float deltaPhiFrom180( float phi1, float phi2 ){
  float dijetDeltaPhi  = fabs(deltaPhi( phi1, phi2 ));
  dijetDeltaPhi *= RADTODEG;
  dijetDeltaPhi = 180 - dijetDeltaPhi;
  return dijetDeltaPhi; 
}


struct triggerFiredBin{
  int xBin;
  int yBin;
  triggerFiredBin():xBin(0),yBin(0){};
  triggerFiredBin(int aXBin, int aYBin):xBin(aXBin),yBin(aYBin){};
};


class dynamicRate{


private:
  // Fired, < xBin, yBin >
  std::vector< triggerFiredBin >     firedBin;
  std::vector< std::vector< bool > > cumulativeFired;



  void integrate(){

    //    std::cout << "1\n";
    // Store fired bins in array
    for (uint iFired = 0; iFired < firedBin.size(); ++iFired){
      //      std::cout << "Fired bin = " << firedBin[iFired].xBin << "\t" << firedBin[iFired].yBin << "\n";
      cumulativeFired[ firedBin[iFired].xBin ][ firedBin[iFired].yBin ] = true;
    }

//      std::cout << "Input:\n";
//      print();

    // Integrate fired bins over y
    for (int x = xBins - 1;x > -1;--x){
      for (int y = yBins - 2;y > -1;--y){
	if ( cumulativeFired[ x ][ y + 1 ] == true ){
	  cumulativeFired[ x ][ y ] = cumulativeFired[ x ][ y + 1 ] ;
	}
      }
    }

//      std::cout << "Int x:\n";
//      print();
    
    // Integrate fired bins over x
    for (int y = yBins - 1;y > -1;--y){
      for (int x = xBins - 2;x > -1;--x){
	if ( cumulativeFired[ x + 1 ][ y ] == true ){
	  cumulativeFired[ x ][ y ] = cumulativeFired[ x + 1 ][ y ] ;
	}
      }
    }

//      std::cout << "Int y:\n";
//      print();

  }
  void addToSum(){
    
    for (int y = 0;y < yBins;++y){
      for (int x = 0;x < xBins;++x){
	
	if ( cumulativeFired[x][y] == true ){ cumulativeSum[x][y]++; }
	else{ break; }
	
      }
    }
    
  }
  void reset(){
    firedBin.clear();
    cumulativeFired.clear();
    cumulativeFired.resize( xBins, std::vector<bool>( yBins, false ) );

  }

public:


  int xBins;
  float xMin;
  float xMax;
  int yBins;
  float yMin;
  float yMax;


  std::vector< std::vector< int > >  cumulativeSum;


  dynamicRate():xBins(0),xMin(0),xMax(0),yBins(0),yMin(0),yMax(0){};
  dynamicRate(int aXBins, float aXMin, float aXMax, int aYBins, float aYMin, float aYMax):xBins(aXBins),xMin(aXMin),xMax(aXMax),yBins(aYBins),yMin(aYMin),yMax(aYMax){
    cumulativeFired.resize( xBins, std::vector<bool>( yBins, false ) );
    cumulativeSum.resize(   xBins, std::vector<int> ( yBins, 0 ) );
  };
  
  void triggerFired( float x, float y ){

    int xBin = floor( (x - xMin)/(xMax - xMin)*xBins);
    int yBin = floor( (y - yMin)/(yMax - yMin)*yBins);

    //    std::cout << "HT = " << x << "\tAlphaT = " << y << "\n";
    
    // Discard underflow
    if ( (xBin < 0) || (yBin < 0) ){ return; }
    // Keep overflow
    if ( xBin >= xBins ){ xBin = xBins - 1; }
    if ( yBin >= yBins ){ yBin = yBins - 1; }

    //    std::cout << "xBin = " << xBin << "\tyBin = " << yBin << "\n\n";
    firedBin.push_back( triggerFiredBin( xBin, yBin ) );
  }
  void endEvent(){
    integrate();
    addToSum();

    reset();
  }
      
  void print(){

//     for (uint y = yBins - 1;y > 0;--y){
//       for (uint x = 0;x < xBins;++x){
// 	std::cout << cumulativeFired[x][y] << " ";
// 	//	std::cout << "(" << x << ", " << y << ") = " << cumulativeFired[x][y] << "\n";
//       }
//       std::cout << "\n";
//     }
    
    for (int y = yBins - 1;y > 0;--y){
      for (int x = 0;x < xBins;++x){
	std::cout << cumulativeSum[x][y] << " ";
	//	std::cout << "(" << x << ", " << y << ") = " << cumulativeFired[x][y] << "\n";
      }
      std::cout << "\n";
    }
    
  }
  
  
};



// ********************************************************************************
// *                                  Selections                                  *
// ********************************************************************************


const int MAX_JETS = 4;


  float secondJetThreshold = 100;

  float caloJetThreshold = 50; //50;
  float pfJetThreshold   = 50; //50;

  // Dynamic alphaT variables
  int maxCaloJet              = 15;
  float caloJetDynThreshold   = 50;  
  float caloJetHTThreshold    = 200;
  float caloJetAlphaThreshold = 0.55;  
  
  // offline selection for measuring calojet efficiency 
  float pfJetHTThreshold     = 200;
  float pfJetAlphaTThreshold = 0.65;
  float pfJet2PTThreshold = 100;

// ********************************************************************************



int    nHTBINS = 3;
double HTBINS[] = {200,275,325,375};

// int    nHTBINS  = 10;
// double HTBINS[] = {200, 275, 325, 375, 475, 575, 675, 775, 875, 975, 1075};

int    nJetBINS = 3;
double JetBINS[] = {1.5,2.5,3.5,4.5};
int    nAlphaTBINS = 10;
double AlphaTBINS[] = {0.3,0.35,0.4,0.45,0.5,0.55,0.60,0.65,0.70,0.75,0.8};





unsigned int maxEvents(100000);


struct sample{
  TString Name;
  TChain* Chain;
  double  CrossSection;

  sample(TString name, TString branch, TString files, double xs):Name(name),CrossSection(xs){
    Chain = new TChain( branch );
    Chain->Add( files );
  }
};

// struct sampleCollection{
//   std::map<TString, sample> Sample;
//   void AddSample( TString name, TString branch, TString files, double xs ){
//     Sample[ name ] = sample( branch, files, xs );
//   }
// }



void makeSUSYHLTAlphaT(){

  TChain* rateChain = new TChain("MakeTrees/Ntuple");
  TChain* sigChain  = new TChain("MakeTrees/Ntuple");


  // ------------------------------------------------------------------------------------------------------------------------
  // Samples
  // ------------------------------------------------------------------------------------------------------------------------
  TString branch    = "MakeTrees/Ntuple";
  TString sampleDir = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/03Oct14/";  // Now includes forward jets, leptons, lepton veto
  
  //QCD_Pt-30to50    = 161500000 pb
  //QCD_Pt-50to80    = 22110000  pb
  //QCD_Pt-80to120   = 3000114.3 pb
  //QCD_Pt-120to170  = 493200    pb
  //QCD_Pt-170to300  = 120300    pb
  //QCD_Pt-300to470  = 7475      pb
  //QCD_Pt-470to600  = 587.1     pb
  //QCD_Pt-600to800  = 167       pb
  //QCD_Pt-800to1000 = 28.25     pb

  sample TTBar        = sample("TTBar",        branch, sampleDir + 
			       "TT_Tune4C_13TeV-pythia8-tauola/TT_Tune4C_13TeV-pythia8-tauola_*.root",                 491.2);
  sample DYJets       = sample("DYJets",       branch, sampleDir + 
			       "DYJetsToLL_M-50_13TeV-madgraph-pythia8/DYJetsToLL_M-50_13TeV-madgraph-pythia8_*.root", 4746);
  sample NuGun        = sample("NuGun",        branch, sampleDir + 
			       "Neutrino_Pt-2to20_gun/PU40bx25_Neutrino_Pt-2to20_gun_*.root", 
			       double(ZB_XSECTION)/(PB_TO_CM2*INST_LUMI_25NS));
  sample QCD30to50    = sample("QCD30to50",    branch, sampleDir + 
			       "QCD_Pt-30to50_Tune4C_13TeV_pythia8/QCD_Pt-30to50_Tune4C_13TeV_pythia8_*.root",       161500000);
  sample QCD50to80    = sample("QCD50to80",    branch, sampleDir + 
			       "QCD_Pt-50to80_Tune4C_13TeV_pythia8/QCD_Pt-50to80_Tune4C_13TeV_pythia8_*.root",       22110000);
  sample QCD80to120   = sample("QCD80to120",   branch, sampleDir + 
			       "QCD_Pt-80to120_Tune4C_13TeV_pythia8/QCD_Pt-80to120_Tune4C_13TeV_pythia8_*.root",     3000114.3);
  sample QCD120to170  = sample("QCD120to170",  branch, sampleDir + 
			       "QCD_Pt-120to170_Tune4C_13TeV_pythia8/QCD_Pt-120to170_Tune4C_13TeV_pythia8_*.root",   493200);
  sample QCD170to300  = sample("QCD170to300",  branch, sampleDir + 
			       "QCD_Pt-170to300_Tune4C_13TeV_pythia8/QCD_Pt-170to300_Tune4C_13TeV_pythia8_*.root",   120300);
  sample QCD300to470  = sample("QCD300to470",  branch, sampleDir + 
			       "QCD_Pt-300to470_Tune4C_13TeV_pythia8/QCD_Pt-300to470_Tune4C_13TeV_pythia8_*.root",   7475);
  sample QCD470to600  = sample("QCD470to600",  branch, sampleDir + 
			       "QCD_Pt-470to600_Tune4C_13TeV_pythia8/QCD_Pt-470to600_Tune4C_13TeV_pythia8_*.root",   587.1);
  sample QCD600to800  = sample("QCD600to800",  branch, sampleDir + 
			       "QCD_Pt-600to800_Tune4C_13TeV_pythia8/QCD_Pt-600to800_Tune4C_13TeV_pythia8_*.root",   167);
  sample QCD800to1000 = sample("QCD800to1000", branch, sampleDir + 
			       "QCD_Pt-800to1000_Tune4C_13TeV_pythia8/QCD_Pt-800to1000_Tune4C_13TeV_pythia8_*.root", 28.25);

  sample T2cc_250_210 = sample("T2cc_250_210", branch, sampleDir + 
			       "T2cc_250_210_Fall13/T2cc_250_210_Fall13*.root", 1);
  sample T2tt_500_250 = sample("T2tt_500_250", branch, sampleDir + 
			       "T2tt_500_250_Fall13/T2tt_500_250_Fall13*.root", 1);

  // ------------------------------------------------------------------------------------------------------------------------
  sample selectedSample = TTBar; //QCD30to50; // T2tt_500_250; //T2cc_250_210; //DYJets; //NuGun; //DYJets; //TTBar; //DYJets;

  TString ntupleType = "";
  #ifdef SIGNAL
  sigChain  = selectedSample.Chain;
  ntupleType = "Eff";
  #endif
  #ifdef NEUTRINO  
  rateChain = selectedSample.Chain;
  ntupleType = "Rate";
  #endif

  // Sample cross section (cm^2)
  double sampleXS = PB_TO_CM2 * selectedSample.CrossSection;
  // ------------------------------------------------------------------------------------------------------------------------


  // rate_mc = collrate * (1 - math.exp(-1* (xs*ilumi*counts/nevt)/collrate))
  //     collrate = (nfillb/mfillb)/xtime
  //     nfillb = the number of filled bunches (usually 2662 for 25ns bunch spacing, 1331 for 50ns bunch spacing)
  //     mfillb = the maximum number of bunches, which is 3564
  //     xtime = the spacing between the bunches (25e-9 for 25ns bunch spacing, 50e-9 for 50ns bunch spacing)
  //     xs = the cross section of your sample multiplied by 1e-36 to convert from pb to cm^2
  //     ilumi = the luminosity you want to find the rate for. For example, the biggest luminosity we can 
  //             achieve in 2015 with 25ns bunch spacing is estimated to be 1.4e34
  //     counts = the number of events that pass your trigger
  //     nevt = the total number of events of your sample


  // ------------------------------------------------------------------------------------------------------------------------


  //  rateChain ->Add("/vols/ssd00/cms/mbaber/AlphaT/Trigger/24Jul14/NeutrinoGun/NeutrinoGun_2014-07-24_*.root" );
  // rateChain ->Add("/vols/ssd00/cms/mbaber/AlphaT/Trigger/07Aug14/NeutrinoGun/NeutrinoGun_2014-08-07_*.root" );
  // sigChain->Add("/vols/ssd00/cms/mbaber/AlphaT/Trigger/24Jul14/TTBar/TTbar_2014-07-24_*.root");
  //    sigChain->Add("/vols/ssd00/cms/mbaber/AlphaT/Trigger/24Jul14/DYJetsToLL/DYJetsToLL_2014-07-24_*.root" );

  unsigned int nuNEvents  = (unsigned int)rateChain ->GetEntries(); 
  unsigned int sigNEvents = (unsigned int)sigChain  ->GetEntries();
  //  double  rateScaleFactor = double(ZB_XSECTION)/nuNEvents;

  // New scale factor calculation for QCD samples
  double rateScaleFactor = (sampleXS * INST_LUMI_25NS)/nuNEvents;

    // rate_mc_simple = (xs * ilumi * counts)/nevt
    // The equation for the error of your rate is:
    // rateerr_mc = xs * ilumi * ((math.sqrt(counts + ((counts)**2)/nevt))/nevt)
    //   rateerr_mc_simple = ((xs * ilumi)/nevt)*math.sqrt(counts) 
    //   So if you use several samples (like QCD binned in Pt), you 



// #ifdef TEST
//   //  rateScaleFactor = double(ZB_XSECTION)/maxEvents;
//   rateScaleFactor = (sampleXS * INST_LUMI_25NS)/maxEvents;
// #endif

  // Output files
  TFile* fOut = new TFile("output/" + selectedSample.Name + "_" + ntupleType + ".root","recreate");

  // Containers
  std::map< TString, TH1*> hist1D;
  std::map< TString, TH2*> hist2D;
  std::map< TString, TH1*> hist1DRate;
  std::map< TString, TH2*> hist2DRate;
  std::map< TString, TEfficiency*> hist1DEff;
  std::map< TString, TEfficiency*> hist2DEff;
  // Efficiency given online selections
  std::map< TString, TEfficiency*> hist2DOnEff;

  // Efficiency given hlt selections
  std::map< TString, TEfficiency*> hist2DHLTEff;
  std::map< TString, TH1*>         histHLTRate;

  std::vector<float> *pfJetPT  = new std::vector<float>();
  std::vector<float> *pfJetPx  = new std::vector<float>();
  std::vector<float> *pfJetPy  = new std::vector<float>();
  std::vector<float> *pfJetEta = new std::vector<float>();
  //  std::vector<float> *pfJetPhi = new std::vector<float>();
  std::vector<float> *pfJetForPT  = new std::vector<float>();
  std::vector<float> *pfJetForPx  = new std::vector<float>();
  std::vector<float> *pfJetForPy  = new std::vector<float>();
  std::vector<float> *pfJetForEta = new std::vector<float>();
  //  std::vector<float> *pfJetForPhi = new std::vector<float>();

  float pfHT(0), pfAlphaT(0);
  std::vector<float> *caloJetPT  = new std::vector<float>();
  std::vector<float> *caloJetPx  = new std::vector<float>();
  std::vector<float> *caloJetPy  = new std::vector<float>();
  std::vector<float> *caloJetEta = new std::vector<float>();
  //  std::vector<float> *caloJetPhi = new std::vector<float>();
  float caloHT(0), caloAlphaT(0);


  std::vector<float> *hltPFJetPT  = new std::vector<float>();
  std::vector<float> *hltPFJetPx  = new std::vector<float>();
  std::vector<float> *hltPFJetPy  = new std::vector<float>();
  std::vector<float> *hltPFJetEta = new std::vector<float>();
  //  std::vector<float> *hltPFJetPhi = new std::vector<float>();
  float hltPFHT(0), hltPFAlphaT(0);


  bool genLeptonVeto(false);

  // ------------------------------------------------------------
  // UCT
  // ------------------------------------------------------------
  std::vector<float> *uctJetPT  = new std::vector<float>();
  std::vector<float> *uctJetPhi = new std::vector<float>();
  float uctHT(0), uctMHT(0), uctMHToverHT(0), uctMET(0);

  //  float uctDijetDeltaPhi(0);
//   float uctJet2PT(0);

  // // ------------------------------------------------------------
  // // HLT
  // // ------------------------------------------------------------
  // // HLT CaloJet
  // std::vector<float> *hltCaloJetPT  = new std::vector<float>();
  // std::vector<float> *hltCaloJetPx  = new std::vector<float>();
  // std::vector<float> *hltCaloJetPy  = new std::vector<float>();
  // std::vector<float> *hltCaloJetEta = new std::vector<float>();
  // //  std::vector<float> *hltCaloJetPhi = new std::vector<float>();
  // float hltCaloJetHT(0), hltCaloJetAlphaT(0);
  // // HLT PF
  // std::vector<float> *hltPFJetPT  = new std::vector<float>();
  // std::vector<float> *hltPFJetPx  = new std::vector<float>();
  // std::vector<float> *hltPFJetPy  = new std::vector<float>();
  // std::vector<float> *hltPFJetEta = new std::vector<float>();
  // //  std::vector<float> *hltPFJetPhi = new std::vector<float>();
  // float hltPFJetHT(0), hltPFJetAlphaT(0);
  // // HLT PFNoPU
  // // std::vector<float> *hltPFNoPUPT  = new std::vector<float>();
  // // std::vector<float> *hltPFNoPUPx  = new std::vector<float>();
  // // std::vector<float> *hltPFNoPUPy  = new std::vector<float>();
  // // std::vector<float> *hltPFNoPUEta = new std::vector<float>();
  // // //  std::vector<float> *hltPFNoPUPhi = new std::vector<float>();
  // // float hltPFNoPUHT(0), hltPFNoPUAlphaT(0);

  // // ------------------------------------------------------------
  // // GEN
  // // ------------------------------------------------------------
  // std::vector<float> *genJetPT  = new std::vector<float>();
  // std::vector<float> *genJetPx  = new std::vector<float>();
  // std::vector<float> *genJetPy  = new std::vector<float>();
  // std::vector<float> *genJetEta = new std::vector<float>();
  // //  std::vector<float> *genJetPhi = new std::vector<float>();
  // float genHT(0), genAlphaT(0), genMET(0);




  // ********************************************************************************
  // *                                Book histograms                               *
  // ********************************************************************************

  std::vector<TString> alphaTTypes;
  alphaTTypes.push_back("alphaT");
  alphaTTypes.push_back("alphaTDyn");
  alphaTTypes.push_back("alphaTDyn2");

  std::vector<float> alphaTCuts;
  alphaTCuts.push_back(0.6);
  alphaTCuts.push_back(0.65);
  alphaTCuts.push_back(0.7);
  alphaTCuts.push_back(0.75);
  alphaTCuts.push_back(0.8);

  std::vector<float> jet2PTCuts;
  jet2PTCuts.push_back(50.);
  jet2PTCuts.push_back(60.);
  jet2PTCuts.push_back(70.);
  jet2PTCuts.push_back(80.);
  jet2PTCuts.push_back(90.);
  jet2PTCuts.push_back(100.);


  std::vector< std::pair<float, float> > anaJetBins;
  anaJetBins.push_back( std::make_pair( 2, 3 ) );
  anaJetBins.push_back( std::make_pair( 3, 4 ) );
  anaJetBins.push_back( std::make_pair( 4, 9 ) );

  std::vector< std::pair<float, float> > anaHtBins;
  // anaHtBins.push_back( std::make_pair( 200, 275 ) );
  // anaHtBins.push_back( std::make_pair( 275, 325 ) );
  // anaHtBins.push_back( std::make_pair( 325, 375 ) );
  // anaHtBins.push_back( std::make_pair( 375, 9999 ) );
  anaHtBins.push_back( std::make_pair( 200, 300 ) );
  anaHtBins.push_back( std::make_pair( 300, 400 ) );
  anaHtBins.push_back( std::make_pair( 400, 500 ) );
  anaHtBins.push_back( std::make_pair( 500, 9999 ) );

  std::vector< float > anaAlphaTBins;
  anaAlphaTBins.push_back( 0.65 );
  anaAlphaTBins.push_back( 0.6 );
  anaAlphaTBins.push_back( 0.55 );
  anaAlphaTBins.push_back( 0.55 );


  
  int   HTBins(50);
  float HTMin(0);
  float HTMax(1000);

  int   MHTBins(50);
  float MHTMin(0);
  float MHTMax(250);

  int alphaTBins(60);
  float alphaTMin(0.4);
  float alphaTMax(1.0);

  int   PTBins(50);
  float PTMin(5);
  float PTMax(305);

  int   HT12Bins(100);
  float HT12Min(5);
  float HT12Max(505);


  hist1DEff["HT_UCT200TurnOn"] = new TEfficiency("HT_UCT200TurnOn","H_{T} turn on (UCT H_{T} > 200);Efficiency",  20,25,1025);
  
  // HT Correlations
  hist2D["CaloHT_vs_CaloHTDyn2"]  = new TH2D("CaloHT_vs_CaloHTDyn2", "Dynamic2 calo H_{T} vs calo H_{T};H_{T}^{Dynamic Calo} (GeV);H_{T}^{Calo} (GeV)",    HTBins,HTMin,HTMax, HTBins,HTMin,HTMax);
  hist2D["PFHT_vs_CaloHT"]        = new TH2D("PFHT_vs_CaloHT",   "PF H_{T} vs calo H_{T};H_{T}^{Calo} (GeV);H_{T}^{PF} (GeV)",    HTBins,HTMin,HTMax, HTBins,HTMin,HTMax);
  hist2D["PFHT_vs_CaloHTDyn2"]    = new TH2D("PFHT_vs_CaloHTDyn2","PF H_{T} vs dynamic2 calo H_{T};H_{T}^{Dynamic Calo} (GeV);H_{T}^{PF} (GeV)",    HTBins,HTMin,HTMax, HTBins,HTMin,HTMax);

  // AlphaT Correlations
  hist2D["CaloAlphaT_vs_CaloAlphaTDyn"]   = new TH2D("CaloAlphaT_vs_CaloAlphaTDyn",   "Calo #alpha_{T} vs dynamic calo #alpha_{T};#alpha_{T}^{Dynamic Calo};#alpha_{T}^{Calo}",    alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
  hist2D["PFAlphaT_vs_CaloAlphaT"]        = new TH2D("PFAlphaT_vs_CaloAlphaT",   "PF #alpha_{T} vs calo #alpha_{T};#alpha_{T}^{Calo};#alpha_{T}^{PF}",    alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
  hist2D["PFAlphaT_vs_CaloAlphaTDyn"]     = new TH2D("PFAlphaT_vs_CaloAlphaTDyn","PF #alpha_{T} vs dynamic calo #alpha_{T};#alpha_{T}^{Dynamic Calo};#alpha_{T}^{PF}",    alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);



  hist2D["UCTMET_vs_UCTMHT"]  = new TH2D("UCTMET_vs_UCTMHT", "UCT #slash{E}_{T} vs #slash{H}_{T};#slash{H}_{T} (GeV);#slash{E}_{T} (GeV)",    MHTBins,MHTMin,MHTMax, MHTBins,MHTMin,MHTMax);
  hist2D["UCTMET_vs_UCTHT"]  = new TH2D("UCTMET_vs_UCTHT", "UCT #slash{E}_{T} vs {H}_{T};{H}_{T} (GeV);#slash{E}_{T} (GeV)",  HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);
  hist2D["UCTMHT_vs_UCTHT"]  = new TH2D("UCTMHT_vs_UCTHT", "UCT #slash{H}_{T} vs {H}_{T};{H}_{T} (GeV);#slash{H}_{T} (GeV)",  HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);

  // ********************************************************************************
  // *                                     Rates                                    *
  // ********************************************************************************

  histHLTRate["HLTCaloHTRate"]                    = new TH1D("HLTCaloHTRate","HLTCaloHTRate;HLT Calo HT (GeV);Rate (Hz)", 100, 0, 2000.); 

  histHLTRate["L1HTT200_HLTCaloHTRate"]           = new TH1D("L1HTT200_HLTCaloHTRate",
							     "L1HTT200_HLTCaloHTRate;HLT Calo HT (GeV);Rate (Hz)", 
							     100, 0, 2000.); 

  histHLTRate["L1HTT200_HLTCaloHT200_HLTPFHTRate"]= new TH1D("L1HTT200_HLTCaloHT200_HLTPFHTRate",
							     "L1HTT200_HLTCaloHT200_HLTPFHTRate;HLT PF HT (GeV);Rate (Hz)", 
							     100, 0, 2000.); 




  hist2DRate["UCTMHT_vs_UCTHT_Rate"] = new TH2D("UCTMHT_vs_UCTHT_Rate","UCT #slash{H}_{T} vs H_{T};H_{T} (GeV);#slash{H}_{T} (GeV)",  HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);
  hist2DRate["UCTMET_vs_UCTHT_Rate"] = new TH2D("UCTMET_vs_UCTHT_Rate","UCT #slash{E}_{T} vs H_{T};H_{T} (GeV);#slash{E}_{T} (GeV)",  HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);
  hist2DRate["UCTMET_vs_UCTHT_MHTdivHT0p3_Rate"] = new TH2D("UCTMET_vs_UCTHT_MHTdivHT0p3_Rate","UCT #slash{E}_{T} vs H_{T} (#slash{E}_{T}/H_{T} > 0.3);H_{T} (GeV);#slash{E}_{T} (GeV)",  HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);

  // MHT vs HT
  hist2DRate["MHT_vs_HT_Rate"] = new TH2D("MHT_vs_HT_Rate","#slash{H}_{T} vs H_{T};H_{T} (GeV);#slash{H}_{T} (GeV)",  HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);

  // AlphaT vs HT

  hist2DRate["AlphaTStd_vs_HT_Rate"] = new TH2D("AlphaTStd_vs_HT_Rate","#alpha_{T} vs H_{T};H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

  hist2DRate["AlphaTStd_vs_HT_Jet2PT100_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2PT100_Rate","#alpha_{T} vs H_{T} p_{T}^{j_{2}} > 100 GeV;H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  hist2DRate["AlphaTStd_vs_HT_Jet2PT90_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2PT90_Rate","#alpha_{T} vs H_{T} p_{T}^{j_{2}} > 90 GeV;H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  hist2DRate["AlphaTStd_vs_HT_Jet2PT80_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2PT80_Rate","#alpha_{T} vs H_{T} p_{T}^{j_{2}} > 80 GeV;H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  hist2DRate["AlphaTStd_vs_HT_Jet2PT70_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2PT70_Rate","#alpha_{T} vs H_{T} p_{T}^{j_{2}} > 70 GeV;H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  hist2DRate["AlphaTStd_vs_HT_Jet2PT60_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2PT60_Rate","#alpha_{T} vs H_{T} p_{T}^{j_{2}} > 60 GeV;H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  hist2DRate["AlphaTStd_vs_HT_Jet2PT50_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2PT50_Rate","#alpha_{T} vs H_{T} p_{T}^{j_{2}} > 50 GeV;H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);




  hist2DRate["AlphaTStd_vs_HT_DynJet2PT_Rate"] = new TH2D("AlphaTStd_vs_HT_DynJet2PT_Rate","#alpha_{T} vs H_{T} dynamic p_{T}^{j_{2}};H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

//   hist2DRate["AlphaTDyn_vs_HT_Rate"] = new TH2D("AlphaTDyn_vs_HT_Rate","#alpha_{T}^{Dynamic} vs H_{T};H_{T} (GeV);#alpha_{T}^{Dynamic}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//   hist2DRate["AlphaTDy2_vs_HT_Rate"] = new TH2D("AlphaTDy2_vs_HT_Rate","#alpha_{T}^{Dynamic} vs H_{T}^{Dynamic};H_{T}^{Dynamic} (GeV);#alpha_{T}^{Dynamic}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

  TH2D *AlphaTStd_vs_HT_Rate = new TH2D("AlphaTStd_vs_HT_Rate_Check","#alpha_{T} vs H_{T};H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  TH2D *AlphaTDyn_vs_HT_Rate = new TH2D("AlphaTDyn_vs_HT_Rate","#alpha_{T}^{Dynamic} vs H_{T};H_{T} (GeV);#alpha_{T}^{Dynamic}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
  TH2D *AlphaTDy2_vs_HT_Rate = new TH2D("AlphaTDy2_vs_HT_Rate","#alpha_{T}^{Dynamic} vs H_{T}^{Dynamic};H_{T}^{Dynamic} (GeV);#alpha_{T}^{Dynamic}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

  hist1DRate["MonoJetPT_Rate"] = new TH1D("MonoJetPT_Rate","Monojet rate;Lead jet p_{T} (GeV);Rate (Hz)",             PTBins,PTMin,PTMax);
  hist1DEff["MonoJetPT200_Eff"]   = new TEfficiency("MonoJetPT200_Eff","Monojet efficiency (Calo p_{T} > 200 (GeV);Lead jet p_{T} (GeV);Rate (Hz)", PTBins,PTMin,PTMax);
  hist1DEff["MonoJetPT210_Eff"]   = new TEfficiency("MonoJetPT210_Eff","Monojet efficiency (Calo p_{T} > 210 (GeV);Lead jet p_{T} (GeV);Rate (Hz)", PTBins,PTMin,PTMax);
  hist1DEff["MonoJetPT220_Eff"]   = new TEfficiency("MonoJetPT220_Eff","Monojet efficiency (Calo p_{T} > 220 (GeV);Lead jet p_{T} (GeV);Rate (Hz)", PTBins,PTMin,PTMax);


  // Validate rates                                                                                                                       
  // --------------------------------------------------------------------------------                                                     
  hist1DRate["Trig1_UCTHT_Rate"] = new TH1D("Trig1_UCTHT_Rate","Trigger 1 H_{T} Rate;H_{T} (GeV);Rate (Hz)", HTBins,HTMin,HTMax);
  hist1DRate["Trig2_UCTHT_Rate"] = new TH1D("Trig2_UCTHT_Rate","Trigger 2 H_{T} Rate;H_{T} (GeV);Rate (Hz)", HTBins,HTMin,HTMax);
  hist1DRate["Trig3_UCTHT_Rate"] = new TH1D("Trig3_UCTHT_Rate","Trigger 3 H_{T} Rate;H_{T} (GeV);Rate (Hz)", HTBins,HTMin,HTMax);
  hist1DRate["Trig4_UCTHT_Rate"] = new TH1D("Trig4_UCTHT_Rate","Trigger 4 H_{T} Rate;H_{T} (GeV);Rate (Hz)", HTBins,HTMin,HTMax);
  hist1DRate["Trig5_UCTHT_Rate"] = new TH1D("Trig5_UCTHT_Rate","Trigger 5 H_{T} Rate;H_{T} (GeV);Rate (Hz)", HTBins,HTMin,HTMax);
  hist1DRate["Trig6_UCTHT_Rate"] = new TH1D("Trig6_UCTHT_Rate","Trigger 6 H_{T} Rate;H_{T} (GeV);Rate (Hz)", HTBins,HTMin,HTMax);
  

  hist1DRate["TrigX_UCTHT_Rate"] = new TH1D("TrigX_UCTHT_Rate","Trigger X Rate;#slash{E}_{T} (GeV);Rate (Hz)", 200,0,200.);
  hist1DRate["TrigY_UCTHT_Rate"] = new TH1D("TrigY_UCTHT_Rate","Trigger Y Rate;#slash{E}_{T} (GeV);Rate (Hz)", 200,0,200.);
  hist1DRate["SingleJet_Rate"]   = new TH1D("SingleJet_Rate",  "SingleJet Rate;p_{T} (GeV);Rate (Hz)", 300,0,300.);
  hist1DRate["DoubleJet_Rate"]   = new TH1D("DoubleJet_Rate",  "DoubleJet Rate;p_{T} (GeV);Rate (Hz)", 300,0,300.);
  hist1DRate["QuadJet_Rate"]     = new TH1D("QuadJet_Rate",  "QuadJet Rate;p_{T} (GeV);Rate (Hz)", 300,0,300.);

  // ********************************************************************************
  // *                                Distributions                                 *
  // ********************************************************************************

  hist1D["CaloJet1PT"]  = new TH1D("CaloJet1PT","Lead calojet p_{T};p_{T} (GeV);Entries",  250,0,250);
  hist1D["CaloJet2PT"]  = new TH1D("CaloJet2PT","Second calojet p_{T};p_{T} (GeV);Entries",250,0,250);
  hist1D["CaloJet3PT"]  = new TH1D("CaloJet3PT","Third calojet p_{T};p_{T} (GeV);Entries", 250,0,250);
  hist1D["CaloJet4PT"]  = new TH1D("CaloJet4PT","Fourth calojet p_{T};p_{T} (GeV);Entries",250,0,250);

  hist1D["CaloJet1Eta"]  = new TH1D("CaloJet1Eta","Lead calojet #eta;#eta;Entries",  62,-3.1,3.1);
  hist1D["CaloJet2Eta"]  = new TH1D("CaloJet2Eta","Second calojet #eta;#eta;Entries",62,-3.1,3.1);
  hist1D["CaloJet3Eta"]  = new TH1D("CaloJet3Eta","Third calojet #eta;#eta;Entries", 62,-3.1,3.1);
  hist1D["CaloJet4Eta"]  = new TH1D("CaloJet4Eta","Fourth calojet #eta;#eta;Entries",62,-3.1,3.1);

  hist1D["CaloJetAlphaT"]     = new TH1D("CaloJetAlphaT",   "Calojet #alpha_{T};#alpha_{T};Entries", 60,0.4,1.);
  hist1D["CaloJetAlphaTDyn"]  = new TH1D("CaloJetAlphaTDyn","Calojet dynamic #alpha_{T};#alpha_{T};Entries", 60,0.4,1.);

  hist1D["PFJet1PT"]  = new TH1D("PFJet1PT","Lead PFjet p_{T};p_{T} (GeV);Entries",  250,0,250);
  hist1D["PFJet2PT"]  = new TH1D("PFJet2PT","Second PFjet p_{T};p_{T} (GeV);Entries",250,0,250);
  hist1D["PFJet3PT"]  = new TH1D("PFJet3PT","Third PFjet p_{T};p_{T} (GeV);Entries", 250,0,250);
  hist1D["PFJet4PT"]  = new TH1D("PFJet4PT","Fourth PFjet p_{T};p_{T} (GeV);Entries",250,0,250);

  hist1D["PFJet1Eta"]  = new TH1D("PFJet1Eta","Lead PFjet #eta;#eta;Entries",  62,-3.1,3.1);
  hist1D["PFJet2Eta"]  = new TH1D("PFJet2Eta","Second PFjet #eta;#eta;Entries",62,-3.1,3.1);
  hist1D["PFJet3Eta"]  = new TH1D("PFJet3Eta","Third PFjet #eta;#eta;Entries", 62,-3.1,3.1);
  hist1D["PFJet4Eta"]  = new TH1D("PFJet4Eta","Fourth PFjet #eta;#eta;Entries",62,-3.1,3.1);

  hist1D["PFJetAlphaT"]     = new TH1D("PFJetAlphaT",   "PFjet #alpha_{T};#alpha_{T};Entries", 60,0.4,1.);


  hist1D["AlphaTRes"]     = new TH1D("AlphaTRes", "#alpha_{T} resolution;(#alpha_{T}^{PF} - #alpha_{T}^{Calo})/#alpha_{T}^{Calo};Entries", 40,-2,2.);



  // Study different determinations of alphaT and HT
  for (uint iType = 0; iType < alphaTTypes.size(); ++iType){

    TString aTType = alphaTTypes[ iType ];

    for (uint iCut = 0; iCut < alphaTCuts.size(); ++iCut){
      
      float alphaTCut = alphaTCuts[ iCut ];
      TString aTLab = Form("%1.2f", alphaTCut );
      TString aTStr = aTLab;
      aTStr.ReplaceAll(".","p");
      
      hist2DRate["HT12_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate"] = new TH2D("HT12_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate","H_{T}^{1,2} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ");Lead jet p_{T} (GeV);H_{T}^{1,2} (GeV)",  PTBins,PTMin,PTMax, HT12Bins,HT12Min,HT12Max);
      hist2DRate["Jet2PT_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate"] = new TH2D("Jet2PT_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate","Second jet p_{T} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ");Lead jet p_{T} (GeV);Second jet p_{T} (GeV)",  PTBins,PTMin,PTMax, PTBins,PTMin,PTMax);
      
      hist2DRate["PFHT12_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate"] = new TH2D("PFHT12_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate","H_{T}^{1,2} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ");Lead jet p_{T} (GeV);H_{T}^{1,2} (GeV)",  PTBins,PTMin,PTMax, HT12Bins,HT12Min,HT12Max);
      hist2DRate["PFJet2PT_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate"] = new TH2D("PFJet2PT_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate","Second jet p_{T} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ");Lead jet p_{T} (GeV);Second jet p_{T} (GeV)",  PTBins,PTMin,PTMax, PTBins,PTMin,PTMax);


      // WW scattering excess 
      hist2DRate["Jet4PT_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate"] = new TH2D("Jet4PT_vs_Jet1PT_" + aTType + "gt" + aTStr + "_Rate","Fourth jet p_{T} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ");Lead jet p_{T} (GeV);Fourth jet p_{T} (GeV)",  PTBins,PTMin,PTMax, PTBins,PTMin,PTMax);
      hist2DEff["Jet4PT_vs_Jet1PT_" + aTType + "gt" + aTStr]  = new TEfficiency("Jet4PT_vs_Jet1PT_" + aTType + "gt" + aTStr,"Fourth jet p_{T} vs lead jet p_{T} efficiency (#alpha_{T} > " + aTLab + ");Lead jet p_{T} (GeV);Fourth jet p_{T} (GeV)", PTBins,PTMin,PTMax, PTBins,PTMin,PTMax);


      for (int jetMult = 2; jetMult <= MAX_JETS; ++jetMult){
	
	TString jetBinStr = TString(Form( "%d", jetMult)) + TString("Jets");
	TString suffix =  aTType + "gt" + aTStr + "_" + jetBinStr + "_Rate";

	hist2DRate["HT12_vs_Jet1PT_" + suffix] = new TH2D("HT12_vs_Jet1PT_" + suffix,"H_{T}^{1,2} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ", " + jetBinStr + ");Lead jet p_{T} (GeV);H_{T}^{1,2} (GeV)",  PTBins,PTMin,PTMax, HT12Bins,HT12Min,HT12Max);
	hist2DRate["Jet2PT_vs_Jet1PT_" + suffix] = new TH2D("Jet2PT_vs_Jet1PT_" + suffix,"Second jet p_{T} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ", " + jetBinStr + ");Lead jet p_{T} (GeV);Second jet p_{T} (GeV)",  PTBins,PTMin,PTMax, PTBins,PTMin,PTMax);
	

	hist2DRate["PFHT12_vs_Jet1PT_" + suffix] = new TH2D("PFHT12_vs_Jet1PT_" + suffix,"H_{T}^{1,2} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ", " + jetBinStr + ");Lead jet p_{T} (GeV);H_{T}^{1,2} (GeV)",  PTBins,PTMin,PTMax, HT12Bins,HT12Min,HT12Max);
	hist2DRate["PFJet2PT_vs_Jet1PT_" + suffix] = new TH2D("PFJet2PT_vs_Jet1PT_" + suffix,"Second jet p_{T} vs lead jet p_{T} rate (#alpha_{T} > " + aTLab + ", " + jetBinStr + ");Lead jet p_{T} (GeV);Second jet p_{T} (GeV)",  PTBins,PTMin,PTMax, PTBins,PTMin,PTMax);


	
	hist2D["CaloAlphaT_vs_CaloAlphaTDyn_" + jetBinStr] = new TH2D("CaloAlphaT_vs_CaloAlphaTDyn_" + jetBinStr, "Calo #alpha_{T} vs dynamic calo #alpha_{T} (" + jetBinStr + ");#alpha_{T}^{Dynamic Calo};#alpha_{T}^{Calo}", alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
	hist2D["PFAlphaT_vs_CaloAlphaT_" + jetBinStr]      = new TH2D("PFAlphaT_vs_CaloAlphaT_" + jetBinStr,   "PF #alpha_{T} vs calo #alpha_{T};#alpha_{T}^{Calo};#alpha_{T}^{PF}",    alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
	hist2D["PFAlphaT_vs_CaloAlphaTDyn_" + jetBinStr]   = new TH2D("PFAlphaT_vs_CaloAlphaTDyn_" + jetBinStr,"PF #alpha_{T} vs dynamic calo #alpha_{T};#alpha_{T}^{Dynamic Calo};#alpha_{T}^{PF}",    alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
	
	

      }

    }
  
  }


  // ********************************************************************************
  // Calojet efficiency for offline selection
  // ********************************************************************************
  
  // Vary the second jet PT cut
  for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
    float jet2PTCut = jet2PTCuts[ iJet2Cut ];
    TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
    TString jet2PTCutLab = "p_{T}^{j2} > " + TString(Form("%1.0f", jet2PTCut ));

    TString stdStr  = "On_AlphaTStd_vs_HT_"  + jet2PTCutStr;
    TString dynStr  = "On_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
    TString dyn2Str = "On_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
    TString dyn3Str = "On_AlphaTDyn3_vs_HT_" + jet2PTCutStr;

    // Inclusive
    hist2DOnEff[stdStr  + "_Inclusive"] = new TEfficiency(stdStr  + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DOnEff[dynStr  + "_Inclusive"] = new TEfficiency(dynStr  + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DOnEff[dyn2Str + "_Inclusive"] = new TEfficiency(dyn2Str + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DOnEff[dyn3Str + "_Inclusive"] = new TEfficiency(dyn3Str + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

    // N_Jet binned
    for (int jetMult = 2; jetMult <= MAX_JETS; ++jetMult){
      TString jetBinStr = TString(Form( "%d", jetMult)) + TString("Jets");
      
      hist2DOnEff[stdStr  + "_" + jetBinStr] = new TEfficiency(stdStr  + "_" + jetBinStr,"Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DOnEff[dynStr  + "_" + jetBinStr] = new TEfficiency(dynStr  + "_" + jetBinStr,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DOnEff[dyn2Str + "_" + jetBinStr] = new TEfficiency(dyn2Str + "_" + jetBinStr,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DOnEff[dyn3Str + "_" + jetBinStr] = new TEfficiency(dyn3Str + "_" + jetBinStr,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      


      // Make candidate HLTtrigger plots                                                                                                            
      for (uint trigN = 1; trigN <= 5; ++trigN ){
	TString trigStr = TString("HLT") + Form("%d", trigN);

	if ( (jet2PTCut == 100.) || (jet2PTCut == 50.) ){
	  hist2DHLTEff[trigStr + stdStr  + "_" + jetBinStr] = new TEfficiency(trigStr + stdStr  + "_" + jetBinStr,"Gen #alpha_{T}^{static} vs H_{T}^{static} efficiency "  + jetBinStr + " (" + jet2PTCutLab + ");Gen H_{T} (GeV);Gen #alpha_{T}",   20.,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	}
      }



    // Analysis binned
    for (uint iHtBin = 0; iHtBin < anaHtBins.size(); ++iHtBin){

      // HT bin
      float htLow  = anaHtBins[iHtBin].first;
      float htHigh = anaHtBins[iHtBin].second;
      TString anaHtStr = TString("HT") + Form("%1.0f", htLow ) + TString("to") + Form("%1.0f", htHigh );
      TString anaHtLab = Form("%1.0f", htLow ) + TString(" #leq H_{T} < ") + Form("%1.0f", htHigh );
      // AlphaT
      float alphaT         = anaAlphaTBins[iHtBin];
      TString anaAlphaTLab = TString("#alpha_{T} > ") + Form("%1.2f", alphaT );

      for (uint iJetBin = 0; iJetBin < anaJetBins.size(); ++iJetBin){
	
	int jetLow  = anaJetBins[ iJetBin ].first;
	int jetHigh = anaJetBins[ iJetBin ].second;
	TString anaJetStr = Form("%d", jetLow) + TString("to") + Form("%d", jetHigh);

	TString suffix = TString("_") + anaHtStr + TString("_") + anaJetStr;
	TString label  = anaHtLab + TString(", ") + anaAlphaTLab + TString(", ") + anaJetStr;

	hist2DOnEff[stdStr  + suffix] = new TEfficiency(stdStr  + suffix,"Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff[dynStr  + suffix] = new TEfficiency(dynStr  + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff[dyn2Str + suffix] = new TEfficiency(dyn2Str + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff[dyn3Str + suffix] = new TEfficiency(dyn3Str + suffix,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


	// Make candidate L1 trigger plots
	for (uint trigN = 1; trigN <= 7; ++trigN ){
	  TString trigStr = TString("Trig") + Form("%d", trigN);

	hist2DOnEff[trigStr + stdStr  + suffix] = new TEfficiency(trigStr + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff[trigStr + dynStr  + suffix] = new TEfficiency(trigStr + dynStr  + suffix,trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff[trigStr + dyn2Str + suffix] = new TEfficiency(trigStr + dyn2Str + suffix,trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff[trigStr + dyn3Str + suffix] = new TEfficiency(trigStr + dyn3Str + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	}
	hist2DOnEff["NoL1" + stdStr  + suffix] = new TEfficiency("NoL1" + stdStr  + suffix,"No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff["NoL1" + dynStr  + suffix] = new TEfficiency("NoL1" + dynStr  + suffix,"No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff["NoL1" + dyn2Str + suffix] = new TEfficiency("NoL1" + dyn2Str + suffix,"No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DOnEff["NoL1" + dyn3Str + suffix] = new TEfficiency("NoL1" + dyn3Str + suffix,"No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);



// 	// Rates
// 	stdStr  = "OnRate_AlphaTStd_vs_HT_"  + jet2PTCutStr;
// 	dynStr  = "OnRate_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
// 	dyn2Str = "OnRate_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
// 	dyn3Str = "OnRate_AlphaTDyn3_vs_HT_" + jet2PTCutStr;

// 	hist2DRate[stdStr  + suffix] = new TH2D(stdStr  + suffix,"Calojet #alpha_{T}^{static} vs H_{T}^{static} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
// 	hist2DRate[dynStr  + suffix] = new TH2D(dynStr  + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
// 	hist2DRate[dyn2Str + suffix] = new TH2D(dyn2Str + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
// 	hist2DRate[dyn3Str + suffix] = new TH2D(dyn3Str + suffix,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	
	
      }
    }


  }
  // ********************************************************************************

  // ********************************************************************************
  // Calojet rate
  // ********************************************************************************
  
  // Vary the second jet PT cut
  for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
    float jet2PTCut = jet2PTCuts[ iJet2Cut ];
    TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
    TString jet2PTCutLab = "p_{T}^{j2} > " + TString(Form("%1.0f", jet2PTCut ));

    TString stdStr  = "OnRate_AlphaTStd_vs_HT_"  + jet2PTCutStr;
    TString dynStr  = "OnRate_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
    TString dyn2Str = "OnRate_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
    TString dyn3Str = "OnRate_AlphaTDyn3_vs_HT_" + jet2PTCutStr;
  
    // Inclusive
    hist2DRate[stdStr  + "_Inclusive"] = new TH2D(stdStr  + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DRate[dynStr  + "_Inclusive"] = new TH2D(dynStr  + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DRate[dyn2Str + "_Inclusive"] = new TH2D(dyn2Str + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DRate[dyn3Str + "_Inclusive"] = new TH2D(dyn3Str + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


    
    // Make candidate L1 trigger plots
    for (uint trigN = 1; trigN <= 7; ++trigN ){
      TString trigStr = TString("Trig") + Form("%d", trigN);

      hist2DRate[trigStr + stdStr  + "_Inclusive"] = new TH2D(trigStr + stdStr  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[trigStr + dynStr  + "_Inclusive"] = new TH2D(trigStr + dynStr  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[trigStr + dyn2Str + "_Inclusive"] = new TH2D(trigStr + dyn2Str + "_Inclusive",trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[trigStr + dyn3Str + "_Inclusive"] = new TH2D(trigStr + dyn3Str + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    }

      hist2DRate["NoL1" + stdStr  + "_Inclusive"] = new TH2D("NoL1" + stdStr  + "_Inclusive","No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate["NoL1" + dynStr  + "_Inclusive"] = new TH2D("NoL1" + dynStr  + "_Inclusive","No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate["NoL1" + dyn2Str + "_Inclusive"] = new TH2D("NoL1" + dyn2Str + "_Inclusive","No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate["NoL1" + dyn3Str + "_Inclusive"] = new TH2D("NoL1" + dyn3Str + "_Inclusive","No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

  
    // N_Jet binned
    for (int jetMult = 2; jetMult <= MAX_JETS; ++jetMult){
      TString jetBinStr = TString(Form( "%d", jetMult)) + TString("Jets");
      
      hist2DRate[stdStr  + "_" + jetBinStr] = new TH2D(stdStr  + "_" + jetBinStr,"Calojet #alpha_{T}^{static} vs H_{T}^{static} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[dynStr  + "_" + jetBinStr] = new TH2D(dynStr  + "_" + jetBinStr,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[dyn2Str + "_" + jetBinStr] = new TH2D(dyn2Str + "_" + jetBinStr,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[dyn3Str + "_" + jetBinStr] = new TH2D(dyn3Str + "_" + jetBinStr,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate " + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      
    }
  }
  }
  // ********************************************************************************


    // N_Jet binned
    for (int jetMult = 2; jetMult <= MAX_JETS; ++jetMult){

      TString jetBinStr = TString(Form( "%d", jetMult)) + TString("Jets");

      hist1DEff["AlphaTStd_TurnOn_" + jetBinStr] = new TEfficiency("AlphaTStd_TurnOn_" + jetBinStr,"#alpha_{T} turn on (" + jetBinStr + ");PF #alpha_{T};Efficiency",  alphaTBins,alphaTMin,alphaTMax);
      hist1DEff["AlphaTDyn_TurnOn_" + jetBinStr] = new TEfficiency("AlphaTDyn_TurnOn_" + jetBinStr,"#alpha_{T} turn on (" + jetBinStr + ");PF #alpha_{T};Efficiency",  alphaTBins,alphaTMin,alphaTMax);
      hist1DEff["AlphaTDy2_TurnOn_" + jetBinStr] = new TEfficiency("AlphaTDy2_TurnOn_" + jetBinStr,"#alpha_{T} turn on (" + jetBinStr + ");PF #alpha_{T};Efficiency",  alphaTBins,alphaTMin,alphaTMax);

      hist1DEff["HTStd_TurnOn_" + jetBinStr] = new TEfficiency("HTStd_TurnOn_" + jetBinStr,"H_{T} turn on (" + jetBinStr + ");PF H_{T};Efficiency",  20,25,1025);
      hist1DEff["HTDyn_TurnOn_" + jetBinStr] = new TEfficiency("HTDyn_TurnOn_" + jetBinStr,"H_{T} turn on (" + jetBinStr + ");PF H_{T};Efficiency",  20,25,1025);
      hist1DEff["HTDy2_TurnOn_" + jetBinStr] = new TEfficiency("HTDy2_TurnOn_" + jetBinStr,"H_{T} turn on (" + jetBinStr + ");PF H_{T};Efficiency",  20,25,1025);
      hist1DEff["HTDy3_TurnOn_" + jetBinStr] = new TEfficiency("HTDy3_TurnOn_" + jetBinStr,"H_{T} turn on (" + jetBinStr + ");PF H_{T};Efficiency",  20,25,1025);


      hist2DEff["AlphaTStd_vs_HT_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDyn_vs_HT_" + jetBinStr] = new TEfficiency("AlphaTDyn_vs_HT_" + jetBinStr,"Dynamic #alpha_{T} vs static H_{T} (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDy2_vs_HT_" + jetBinStr] = new TEfficiency("AlphaTDy2_vs_HT_" + jetBinStr,"Dynamic #alpha_{T} vs dynamic H_{T} (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDy3_vs_HT_" + jetBinStr] = new TEfficiency("AlphaTDy3_vs_HT_" + jetBinStr,"Static #alpha_{T} vs dynamic H_{T} (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

//       hist2DOnEff["On_AlphaTStd_vs_HT_" + jetBinStr] = new TEfficiency("On_AlphaTStd_vs_HT_" + jetBinStr,"Online Static #alpha_{T} vs static H_{T} efficiency (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//       hist2DOnEff["On_AlphaTDyn_vs_HT_" + jetBinStr] = new TEfficiency("On_AlphaTDyn_vs_HT_" + jetBinStr,"Online Dynamic #alpha_{T} vs static H_{T} efficiency (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//       hist2DOnEff["On_AlphaTDy2_vs_HT_" + jetBinStr] = new TEfficiency("On_AlphaTDy2_vs_HT_" + jetBinStr,"Online Dynamic #alpha_{T} vs dynamic H_{T} efficiency (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//       hist2DOnEff["On_AlphaTDy3_vs_HT_" + jetBinStr] = new TEfficiency("On_AlphaTDy3_vs_HT_" + jetBinStr,"Online Static #alpha_{T} vs dynamic H_{T} efficiency (" + jetBinStr + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

    

      hist2DEff["AlphaTStd_vs_HT_Jet2gt50_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_Jet2gt50_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 50 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTStd_vs_HT_Jet2gt60_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_Jet2gt60_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 60 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTStd_vs_HT_Jet2gt70_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_Jet2gt70_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 70 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTStd_vs_HT_Jet2gt80_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_Jet2gt80_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 80 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTStd_vs_HT_Jet2gt90_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_Jet2gt90_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 90 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

      hist2DEff["AlphaTStd_vs_HT_Jet2gt100_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_Jet2gt100_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDyn_vs_HT_Jet2gt100_" + jetBinStr] = new TEfficiency("AlphaTDyn_vs_HT_Jet2gt100_" + jetBinStr,"Dynamic #alpha_{T} vs static H_{T} (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDy2_vs_HT_Jet2gt100_" + jetBinStr] = new TEfficiency("AlphaTDy2_vs_HT_Jet2gt100_" + jetBinStr,"Dynamic #alpha_{T} vs dynamic H_{T} (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDy3_vs_HT_Jet2gt100_" + jetBinStr] = new TEfficiency("AlphaTDy3_vs_HT_Jet2gt100_" + jetBinStr,"Static #alpha_{T} vs dynamic H_{T} (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

      hist2DEff["AlphaTStd_vs_HT_DynJet2gt100_" + jetBinStr] = new TEfficiency("AlphaTStd_vs_HT_DynJet2gt100_" + jetBinStr,"Static #alpha_{T} vs static H_{T} (" + jetBinStr + ", dynamic p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDyn_vs_HT_DynJet2gt100_" + jetBinStr] = new TEfficiency("AlphaTDyn_vs_HT_DynJet2gt100_" + jetBinStr,"Dynamic #alpha_{T} vs static H_{T} (" + jetBinStr + ", dynamic p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDy2_vs_HT_DynJet2gt100_" + jetBinStr] = new TEfficiency("AlphaTDy2_vs_HT_DynJet2gt100_" + jetBinStr,"Dynamic #alpha_{T} vs dynamic H_{T} (" + jetBinStr + ", dynamic p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["AlphaTDy3_vs_HT_DynJet2gt100_" + jetBinStr] = new TEfficiency("AlphaTDy3_vs_HT_DynJet2gt100_" + jetBinStr,"Static #alpha_{T} vs dynamic H_{T} (" + jetBinStr + ", dynamic p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);




//       hist2DRate["AlphaTStd_vs_HT_Jet2gt100_" + jetBinStr + "_Rate"] = new TH2D("AlphaTStd_vs_HT_Jet2gt100_" + jetBinStr + "_Rate","Static #alpha_{T} vs static H_{T} rate (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//       hist2DRate["AlphaTDyn_vs_HT_Jet2gt100_" + jetBinStr + "_Rate"] = new TH2D("AlphaTDyn_vs_HT_Jet2gt100_" + jetBinStr + "_Rate","Dynamic #alpha_{T} vs static H_{T} rate (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//       hist2DRate["AlphaTDy2_vs_HT_Jet2gt100_" + jetBinStr + "_Rate"] = new TH2D("AlphaTDy2_vs_HT_Jet2gt100_" + jetBinStr + "_Rate","Dynamic #alpha_{T} vs dynamic H_{T} rate (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
//       hist2DRate["AlphaTDy3_vs_HT_Jet2gt100_" + jetBinStr + "_Rate"] = new TH2D("AlphaTDy2_vs_HT_Jet2gt100_" + jetBinStr + "_Rate","Static #alpha_{T} vs dynamic H_{T} rate (" + jetBinStr + ", p_{T}^{j2} > 100 GeV);H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

    }


    //    hist2DEff["PF_" + trig + "_NJetvsHT_L1"]    = new TEfficiency("PF_" + trig + "_NJetvsHT_L1","Low H_{T} #alpha_{T} analysis efficiency;H_{T} (GeV);N_{Jet}",    nHTBINS, HTBINS, nJetBINS, JetBINS);




  





  // ********************************************************************************
  // *                                  Efficiency                                  *
  // ********************************************************************************


    sigChain->SetBranchAddress("genLeptonVeto", &genLeptonVeto);


  // ------------------------------------------------------------ 
  // UCT 
  // ------------------------------------------------------------ 

    // NOTE: Currently with HLT emulation the lx1extra GCT products are replaced with UCT so this is a labeling problem
    sigChain->SetBranchAddress("gctCen_Pt",       &uctJetPT);
    sigChain->SetBranchAddress("gctCen_Phi",      &uctJetPhi);
    sigChain->SetBranchAddress("gct_Ht",          &uctHT);
    sigChain->SetBranchAddress("gct_MhtPt",       &uctMHT);
    sigChain->SetBranchAddress("gct_MhtDivHt",    &uctMHToverHT);
    sigChain->SetBranchAddress("gct_MetPt",       &uctMET);


  // // ------------------------------------------------------------ 
  // // HLT 
  // // ------------------------------------------------------------ 
  //   // HLT CaloJet
  //   sigChain->SetBranchAddress("hltAk4Calo_Pt",        &hltCaloJetPT);
  //   sigChain->SetBranchAddress("hltAk4Calo_Px",        &hltCaloJetPx);
  //   sigChain->SetBranchAddress("hltAk4Calo_Py",        &hltCaloJetPy);
  //   //sigChain->SetBranchAddress("hltAk4Calo_Phi",       &hltCaloJetPhi);
  //   sigChain->SetBranchAddress("hltAk4Calo_Eta",       &hltCaloJetEta);
  
  //   // HLT PF 
  //   sigChain->SetBranchAddress("hltAk4PF_Pt",          &hltPFJetPT); 
  //   sigChain->SetBranchAddress("hltAk4PF_Px",          &hltPFJetPx);
  //   sigChain->SetBranchAddress("hltAk4PF_Py",          &hltPFJetPy);
  //   //sigChain->SetBranchAddress("hltAk4PF_Phi",         &hltPFJetPhi);
  //   sigChain->SetBranchAddress("hltAk4PF_Eta",         &hltPFJetEta);
    
  //   // // HLT PFNoPU 
  //   // sigChain->SetBranchAddress("hltAk4PFNoPU_Pt",      &hltPFNoPUPT); 
  //   // sigChain->SetBranchAddress("hltAk4PFNoPU_Px",      &hltPFNoPUPx);
  //   // sigChain->SetBranchAddress("hltAk4PFNoPU_Py",      &hltPFNoPUPy);
  //   // //sigChain->SetBranchAddress("hltAk4PFNoPU_Phi",     &hltPFNoPUPhi);
  //   // sigChain->SetBranchAddress("hltAk4PFNoPU_Eta",     &hltPFNoPUEta);
    
    
  // ------------------------------------------------------------ 
  // GEN
  // ------------------------------------------------------------ 

    // sigChain->SetBranchAddress("genAk4_Pt",         &genJetPT); 
    // sigChain->SetBranchAddress("genAk4_Px",         &genJetPx);
    // sigChain->SetBranchAddress("genAk4_Py",         &genJetPy);
    // //sigChain->SetBranchAddress("genAk4_Phi",      &genJetPhi);
    // sigChain->SetBranchAddress("genAk4_Eta",        &genJetEta);
    // sigChain->SetBranchAddress("genMetTrue_MetPt", &genMET);


    // TEMPORARY - REPLACE NAMES FOR MEETING IN 1.5 hours
    sigChain->SetBranchAddress("genAk4_Pt",         &pfJetPT);
    sigChain->SetBranchAddress("genAk4_Px",         &pfJetPx);
    sigChain->SetBranchAddress("genAk4_Py",         &pfJetPy);
    //sigChain->SetBranchAddress("genAk4_Phi",      &pfJetPhi);
    sigChain->SetBranchAddress("genAk4_Eta",        &pfJetEta);

    sigChain->SetBranchAddress("genAk4For_Pt",         &pfJetForPT);
    sigChain->SetBranchAddress("genAk4For_Px",         &pfJetForPx);
    sigChain->SetBranchAddress("genAk4For_Py",         &pfJetForPy);
    //sigChain->SetBranchAddress("genAk4For_Phi",      &pfJetForPhi);
    sigChain->SetBranchAddress("genAk4For_Eta",        &pfJetForEta);

#ifdef HLT_CALOJET
    // HLT CaloJet
    sigChain->SetBranchAddress("hltAk4Calo_Pt",        &caloJetPT);
    sigChain->SetBranchAddress("hltAk4Calo_Px",        &caloJetPx);
    sigChain->SetBranchAddress("hltAk4Calo_Py",        &caloJetPy);
    //sigChain->SetBranchAddress("hltAk4Calo_Phi",       &caloJetPhi);
    sigChain->SetBranchAddress("hltAk4Calo_Eta",       &caloJetEta);

    // sigChain->SetBranchAddress("hltAk4CaloFor_Pt",        &caloJetForPT);
    // sigChain->SetBranchAddress("hltAk4CaloFor_Px",        &caloJetForPx);
    // sigChain->SetBranchAddress("hltAk4CaloFor_Py",        &caloJetForPy);
    // //sigChain->SetBranchAddress("hltAk4CaloFor_Phi",       &caloJetForPhi);
    // sigChain->SetBranchAddress("hltAk4CaloFor_Eta",       &caloJetForEta);
#else
    // HLT PFJet
    sigChain->SetBranchAddress("hltAk4PF_Pt",        &caloJetPT);
    sigChain->SetBranchAddress("hltAk4PF_Px",        &caloJetPx);
    sigChain->SetBranchAddress("hltAk4PF_Py",        &caloJetPy);
    //sigChain->SetBranchAddress("hltAk4PF_Phi",       &caloJetPhi);
    sigChain->SetBranchAddress("hltAk4PF_Eta",       &caloJetEta);

    // sigChain->SetBranchAddress("hltAk4PFFor_Pt",        &caloJetForPT);
    // sigChain->SetBranchAddress("hltAk4PFFor_Px",        &caloJetForPx);
    // sigChain->SetBranchAddress("hltAk4PFFor_Py",        &caloJetForPy);
    // //sigChain->SetBranchAddress("hltAk4PFFor_Phi",       &caloJetForPhi);
    // sigChain->SetBranchAddress("hltAk4PFFor_Eta",       &caloJetForEta);
#endif

  

  // sigChain->SetBranchAddress("jetPtsPf",        &pfJetPT);
  // //  sigChain->SetBranchAddress("jetPhisPf",       &pfJetPhi);
  // sigChain->SetBranchAddress("jetPxsPf",        &pfJetPx);
  // sigChain->SetBranchAddress("jetPysPf",        &pfJetPy);
  // sigChain->SetBranchAddress("jetEtasPf",       &pfJetEta);
  // //  sigChain->SetBranchAddress("htPf",            &pfHT);
  // sigChain->SetBranchAddress("alphaTPf",        &pfAlphaT);

  // sigChain->SetBranchAddress("jetPtsCalo",      &caloJetPT);
  // sigChain->SetBranchAddress("jetPxsCalo",      &caloJetPx);
  // sigChain->SetBranchAddress("jetPysCalo",      &caloJetPy);
  // sigChain->SetBranchAddress("jetEtasCalo",     &caloJetEta);
  // //  sigChain->SetBranchAddress("jetPhisCalo",     &caloJetPhi);
  // //  sigChain->SetBranchAddress("htCalo",          &caloHT);
  // sigChain->SetBranchAddress("alphaTCalo",      &caloAlphaT);



  //  unsigned int nEvents = (unsigned int)ttTree->GetEntries();



  unsigned int signalEventLow  = 0;
  unsigned int signalEventHigh = sigNEvents;
  unsigned int nuEventLow      = 0;
  unsigned int nuEventHigh     = nuNEvents;

#ifdef RUN_ON_BATCH
  signalEventLow               = floor( eventLowFact*sigNEvents) + eventLowOffset;
  signalEventHigh              = floor(eventHighFact*sigNEvents) + eventHighOffset;

  nuEventLow                   = floor( eventLowFact*nuNEvents) + eventLowOffset;
  nuEventHigh                  = floor(eventHighFact*nuNEvents) + eventHighOffset;

#endif

  std::cout << "\n\nProcessing signal events:\n"
	    << "\tLow   = " << signalEventLow  
	    << "\tHigh  = " << signalEventHigh 
	    << "\tTotal = " << sigNEvents << "\n\n";

  std::cout << "Processing neutrino events:\n"
	    << "\tLow   = " << nuEventLow  
	    << "\tHigh  = " << nuEventHigh 
	    << "\tTotal = " << nuNEvents << "\n\n";

  std::cout << "Rate scalefactor = " << rateScaleFactor << "\n\n";


#ifdef SIGNAL

  // Loop over the tree
  //  for ( unsigned int iEvent = 0; iEvent < nEvents; ++iEvent ){
  for ( unsigned int iEvent = signalEventLow; iEvent < signalEventHigh; ++iEvent ){
    
    sigChain->GetEntry( iEvent );

    float pfAlphaTStandard(0);

    float caloAlphaTStandard(0);
    float caloAlphaTDynamic(0);
    float caloHTDynamic(0);
    std::pair<float,float> caloAlphaTHTDynamic;

    

    // ********************************************************************************
    // *                                 UCT triggers                                 *
    // ********************************************************************************

    bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) ); // DJ100
    bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 60) ); // DJ120
    bool l1Trig3 = ( (uctHT >= 160)  || (uctMET >= 70) ); // DJ120
    bool l1Trig4 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 150) && (uctMHToverHT >= 0.17)) );
    bool l1Trig5 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 140) && (uctMHToverHT >= 0.25)) );
    bool l1Trig6 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 130) && (uctMHToverHT >= 0.32)) );
    bool l1Trig7 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 120) && (uctMHToverHT >= 0.36)) );

    // bool l1Trig4 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHT >= 56)) );
    // bool l1Trig5 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 125) && (uctMHT >= 44)) );
    // bool l1Trig6 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHToverHT >= 0.4)) );
    // bool l1Trig7 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 125) && (uctMHToverHT >= 0.3)) );


    // bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) );
    // bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 125) && (uctMHToverHT >= 0.3)) );
    // bool l1Trig3 = ( (uctHT >= 200)  || (uctMET >= 60) );
    // bool l1Trig4 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHT >= 56)) );
    // bool l1Trig5 = ( (uctMET >= 70)  || ((uctHT >= 125) && (uctMHT >= 44)) );
    // bool l1Trig6 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHT >= 56)) );

    hist2D["UCTMET_vs_UCTMHT"]->Fill( uctMHT, uctMET );
    hist2D["UCTMET_vs_UCTHT"] ->Fill( uctHT,  uctMET );
    hist2D["UCTMHT_vs_UCTHT"] ->Fill( uctHT,  uctMHT );

    //    if (pfHT > 200){

    // calculate HT and jet multiplicity for new jet threshold
    int caloJetsAboveThresh(0);
    caloHT = 0;
    for(uint iJet = 0;iJet < caloJetPT->size(); ++iJet ){
      if( (*caloJetPT)[iJet] < caloJetThreshold ) break;
      caloHT += (*caloJetPT)[iJet];
      caloJetsAboveThresh++;
    }
    if ( caloJetsAboveThresh > MAX_JETS ){ 
      caloJetsAboveThresh = MAX_JETS; 
    }
    TString nCaloJetsStr = Form("%d", caloJetsAboveThresh );
    int pfJetsAboveThresh(0);
    pfHT = 0;
    for(uint iJet = 0;iJet < pfJetPT->size(); ++iJet ){
      if( (*pfJetPT)[iJet] < pfJetThreshold ) break;
      pfHT += (*pfJetPT)[iJet];
      pfJetsAboveThresh++;
    }
    if ( pfJetsAboveThresh > MAX_JETS ){
      pfJetsAboveThresh = MAX_JETS;
    }
    TString nPFJetsStr = Form("%d", pfJetsAboveThresh );

    // Monojets
    if ( (pfJetPT->size() == 1) || ((pfJetPT->size() > 1) && (*pfJetPT)[1] < pfJetThreshold) ){
      bool passesMonojet200(false), passesMonojet210(false), passesMonojet220(false);
      if ( (caloJetPT->size() == 1) || ((caloJetPT->size() > 1) && (*caloJetPT)[1] < caloJetThreshold) ){
	if ((*caloJetPT)[0] > 200 ){ passesMonojet200 = true; }
	if ((*caloJetPT)[0] > 210 ){ passesMonojet210 = true; }
	if ((*caloJetPT)[0] > 220 ){ passesMonojet220 = true; }
      }
      hist1DEff["MonoJetPT200_Eff"]->Fill( passesMonojet200, (*pfJetPT)[0]);
      hist1DEff["MonoJetPT210_Eff"]->Fill( passesMonojet210, (*pfJetPT)[0]);
      hist1DEff["MonoJetPT220_Eff"]->Fill( passesMonojet220, (*pfJetPT)[0]);


      bool passesUCTHT200(false);
      if (uctHT > 200){ passesUCTHT200 = true; }
      hist1DEff["HT_UCT200TurnOn"] ->Fill( passesUCTHT200, (*pfJetPT)[0]);

    }


    // Jet distributions
    if ( caloJetPT->size() > 0 ){
      hist1D["CaloJet1PT"] ->Fill( (*caloJetPT)[0] ); 
      hist1D["CaloJet1Eta"]->Fill( (*caloJetEta)[0] ); 
      if ( caloJetPT->size() > 1 ){
	hist1D["CaloJet2PT"] ->Fill( (*caloJetPT)[1] ); 
	hist1D["CaloJet2Eta"]->Fill( (*caloJetEta)[1] ); 
	if ( caloJetPT->size() > 2 ){
	  hist1D["CaloJet3PT"] ->Fill( (*caloJetPT)[2] ); 
	  hist1D["CaloJet3Eta"]->Fill( (*caloJetEta)[2] ); 
	  if ( caloJetPT->size() > 3 ){
	    hist1D["CaloJet4PT"] ->Fill( (*caloJetPT)[3] ); 
	    hist1D["CaloJet4Eta"]->Fill( (*caloJetEta)[3] ); 
	  }
	}
      }
    }
    if ( pfJetPT->size() > 0 ){
      hist1D["PFJet1PT"]  ->Fill( (*pfJetPT)[0] ); 
      hist1D["PFJet1Eta"] ->Fill( (*pfJetEta)[0] ); 
      if ( pfJetPT->size() > 1 ){
	hist1D["PFJet2PT"]  ->Fill( (*pfJetPT)[1] ); 
	hist1D["PFJet2Eta"] ->Fill( (*pfJetEta)[1] ); 
	if ( pfJetPT->size() > 2 ){
	  hist1D["PFJet3PT"]  ->Fill( (*pfJetPT)[2] ); 
	  hist1D["PFJet3Eta"] ->Fill( (*pfJetEta)[2] ); 
	  if ( pfJetPT->size() > 3 ){
	    hist1D["PFJet4PT"]  ->Fill( (*pfJetPT)[3] ); 
	    hist1D["PFJet4Eta"] ->Fill( (*pfJetEta)[3] ); 
	  }
	}
      }
    }



    // Calo HLT cuts
    bool passesATStandard(false), passesATDynamic(false), passesATDynamic2(false), passesATDynamic3(false);
    bool passesHTStandard(false), passesHTDynamic(false);
    bool passesJet2_50(false), passesJet2_60(false), passesJet2_70(false), passesJet2_80(false), passesJet2_90(false), passesJet2_100(false);
    bool passesDynJet2_100(false);    

    bool passesOffHT(false), passesOffAT(false), passesOffJet(false), passesOffAll(false);
    bool passesAnaBinOffAT(false), passesAnaBinOffJet(false), passesAnaBinOffAll(false);
    bool passesOffVetoes(false);



    pfAlphaTStandard    = calculateAlphaT( pfJetPT, pfJetPx, pfJetPy, pfJetThreshold );      

    caloAlphaTStandard  = calculateAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold );    
    // Dynamic AlphaT
    caloAlphaTDynamic   = calculateDynamicAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold, 
						  caloJetAlphaThreshold );
    // Dynamic HT and AlphaT
    caloAlphaTHTDynamic = calculateDynamicAlphaTHT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold,
						    caloJetAlphaThreshold, caloJetHTThreshold );

    hist1D["CaloJetAlphaT"]   ->Fill( caloAlphaTStandard ); 
    hist1D["CaloJetAlphaTDyn"]->Fill( caloAlphaTDynamic ); 

    hist1D["PFJetAlphaT"]   ->Fill( pfAlphaTStandard ); 

    float alphaTRes = 0;

    // HT Correlations
    hist2D["CaloHT_vs_CaloHTDyn2"]->Fill( caloAlphaTHTDynamic.second,  caloHT );	
    hist2D["PFHT_vs_CaloHTDyn2"]  ->Fill( caloAlphaTHTDynamic.second,  pfHT );	
    hist2D["PFHT_vs_CaloHT"]      ->Fill( caloHT, pfHT );

    // AlphaT Correlations
    hist2D["CaloAlphaT_vs_CaloAlphaTDyn"] ->Fill( caloAlphaTDynamic,          caloAlphaTStandard );	
    hist2D["PFAlphaT_vs_CaloAlphaT"]      ->Fill( caloAlphaTStandard,         pfAlphaTStandard );	
    hist2D["PFAlphaT_vs_CaloAlphaTDyn"]   ->Fill( caloAlphaTDynamic,          pfAlphaTStandard );	

    // ********************************************************************************
    // *                              Calojets 
    // ********************************************************************************

    if ( pfJetsAboveThresh >= 2 ){
      //      TString caloJetBinStr = TString(Form( "%d", caloJetsAboveThresh)) + TString("Jets");
      TString pfJetBinStr     = TString(Form( "%d", pfJetsAboveThresh)) + TString("Jets");

      
      // Forward jet veto
      bool forJetVeto(false);
      if (pfJetForPT->size() > 0){
	if ( (*pfJetForPT)[0] > pfJetThreshold){
	  forJetVeto = true;
	}
      }



      // Check individual offline cuts
      if ( pfHT             > pfJetHTThreshold)                                 { passesOffHT  = true; }
      if ( pfAlphaTStandard > pfJetAlphaTThreshold)                             { passesOffAT  = true; }
      if ( (pfJetsAboveThresh >= 2 ) && ((*pfJetPT)[1] > pfJet2PTThreshold) )   { passesOffJet = true; }
      if ( !(genLeptonVeto) && !(forJetVeto) )                                  { passesOffVetoes = true; }
      // Check full offline selection
      if ( passesOffHT && passesOffAT && passesOffJet && passesOffVetoes )      { passesOffAll = true; }

      // Vary the online second jet PT cut
      for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){

	float jet2PTCut = jet2PTCuts[ iJet2Cut ];
	TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
	TString stdStr  = "On_AlphaTStd_vs_HT_"  + jet2PTCutStr;
	TString dynStr  = "On_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
	TString dyn2Str = "On_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
	TString dyn3Str = "On_AlphaTDyn3_vs_HT_" + jet2PTCutStr;

	bool passHltSecondJet(false);

	// Second jet threshold
	if ((*caloJetPT)[1] > jet2PTCut){
	  passHltSecondJet = true;
	}
	// // Inclusive 
	// hist2DOnEff[stdStr  + "_Inclusive"]       ->Fill( passesOffAll, caloHT,                     caloAlphaTStandard );
	// hist2DOnEff[dynStr  + "_Inclusive"]       ->Fill( passesOffAll, caloHT,                     caloAlphaTDynamic );
	// hist2DOnEff[dyn2Str + "_Inclusive"]       ->Fill( passesOffAll, caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	// hist2DOnEff[dyn3Str + "_Inclusive"]       ->Fill( passesOffAll, caloAlphaTHTDynamic.second, caloAlphaTStandard );
	// // Jet binned
	// hist2DOnEff[stdStr  + "_" + pfJetBinStr]->Fill( passesOffAll, caloHT,                     caloAlphaTStandard );
	// hist2DOnEff[dynStr  + "_" + pfJetBinStr]->Fill( passesOffAll, caloHT,                     caloAlphaTDynamic );
	// hist2DOnEff[dyn2Str + "_" + pfJetBinStr]->Fill( passesOffAll, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	// hist2DOnEff[dyn3Str + "_" + pfJetBinStr]->Fill( passesOffAll, caloAlphaTHTDynamic.second, caloAlphaTStandard ); 
	
	

	// Efficiency offline given an HLT path (no L1)

	// ********************************************************************************
	// *                                 HLT triggers                                 *
	// ********************************************************************************
	// Rate balanced triggers to 20 Hz
	bool HLT1(false), HLT2(false), HLT3(false), HLT4(false), HLT5(false);

#ifdef HLT_CALOJET
	if (jet2PTCut == 100.){
	  HLT1 = ( (caloHT >= 240)  && (caloAlphaTStandard >= 0.92) );
	  HLT2 = ( (caloHT >= 280)  && (caloAlphaTStandard >= 0.71) );
	  HLT3 = ( (caloHT >= 320)  && (caloAlphaTStandard >= 0.64) );
	  HLT4 = ( (caloHT >= 360)  && (caloAlphaTStandard >= 0.60) );
	  HLT5 = ( (caloHT >= 400)  && (caloAlphaTStandard >= 0.58) );
	}
	else if(jet2PTCut == 50.){
	  HLT4 = ( (caloHT >= 360)  && (caloAlphaTStandard >= 0.78) );
	  HLT5 = ( (caloHT >= 400)  && (caloAlphaTStandard >= 0.69) );
	}
#else
	if (jet2PTCut == 100.){
	  HLT1 = ( (caloHT >= 240)  && (caloAlphaTStandard >= 0.90) );
	  HLT2 = ( (caloHT >= 280)  && (caloAlphaTStandard >= 0.64) );
	  HLT3 = ( (caloHT >= 320)  && (caloAlphaTStandard >= 0.61) );
	  HLT4 = ( (caloHT >= 360)  && (caloAlphaTStandard >= 0.57) );
	  HLT5 = ( (caloHT >= 400)  && (caloAlphaTStandard >= 0.55) );
	}
	else if(jet2PTCut == 50.){
	  HLT3 = ( (caloHT >= 320)  && (caloAlphaTStandard >= 0.79) );
	  HLT4 = ( (caloHT >= 360)  && (caloAlphaTStandard >= 0.64) );
	  HLT5 = ( (caloHT >= 400)  && (caloAlphaTStandard >= 0.59) );
	}
#endif

	if (passesOffVetoes && passesOffJet){
	  TString trigStr  = "";
	  bool    trigBool = false;
	  if ( (jet2PTCut == 100.) || (jet2PTCut == 50.) ){
	    trigStr  = "HLT1";
	    trigBool = HLT1;
	    hist2DHLTEff[trigStr + stdStr  + "_" + pfJetBinStr]->Fill( trigBool, pfHT,                     pfAlphaTStandard );
	    trigStr  = "HLT2";
	    trigBool = HLT2;
	    hist2DHLTEff[trigStr + stdStr  + "_" + pfJetBinStr]->Fill( trigBool, pfHT,                     pfAlphaTStandard );
	    trigStr  = "HLT3";
	    trigBool = HLT3;
	    hist2DHLTEff[trigStr + stdStr  + "_" + pfJetBinStr]->Fill( trigBool, pfHT,                     pfAlphaTStandard );
	    trigStr  = "HLT4";
	    trigBool = HLT4;
	    hist2DHLTEff[trigStr + stdStr  + "_" + pfJetBinStr]->Fill( trigBool, pfHT,                     pfAlphaTStandard );
	    trigStr  = "HLT5";
	    trigBool = HLT5;
	    hist2DHLTEff[trigStr + stdStr  + "_" + pfJetBinStr]->Fill( trigBool, pfHT,                     pfAlphaTStandard );
	  }
	}
      

     
    




















	// Get offline analysis bin
	for (int iHtBin = 0; iHtBin < anaHtBins.size(); ++iHtBin){
	  // PF - HT bin
	  float htLow  = anaHtBins[iHtBin].first;
	  float htHigh = anaHtBins[iHtBin].second;
	  if ( !((pfHT >= htLow) && (pfHT < htHigh)) ){ continue; }
	  TString anaHtStr = TString("HT") + Form("%1.0f", htLow ) + TString("to") + Form("%1.0f", htHigh );
	  
	  // PF - AlphaT
	  float alphaT         = anaAlphaTBins[iHtBin];
	  TString anaAlphaTLab = Form("%1.0f", alphaT );

	  // Check individual offline cuts
	  // --------------------------------------------------------------------------------
	  // NOTE: HT is already satisfied by binning
	  if ( pfAlphaTStandard > alphaT)                                        { passesAnaBinOffAT  = true; }
	  if ( (pfJetsAboveThresh >= 2 ) && ((*pfJetPT)[1] > pfJet2PTThreshold) ){ passesAnaBinOffJet = true; }
	  // Check full offline selection
	  if ( passesOffAT && passesOffJet )                                     { passesAnaBinOffAll = true; }


	  // Get jet bin
	  for (int iJetBin = 0; iJetBin < anaJetBins.size(); ++iJetBin){
	
	    int jetLow  = anaJetBins[ iJetBin ].first;
	    int jetHigh = anaJetBins[ iJetBin ].second;
	    if ( !((pfJetsAboveThresh >= jetLow) && (pfJetsAboveThresh < jetHigh)) ){ continue; }
	    TString anaJetStr = Form("%d", jetLow) + TString("to") + Form("%d", jetHigh);
	    
	    TString suffix = TString("_") + anaHtStr + TString("_") + anaJetStr;
	    
	    hist2DOnEff[stdStr  + suffix] ->Fill( passesAnaBinOffAll, caloHT,                     caloAlphaTStandard );
	    hist2DOnEff[dynStr  + suffix] ->Fill( passesAnaBinOffAll, caloHT,                     caloAlphaTDynamic );
	    hist2DOnEff[dyn2Str + suffix] ->Fill( passesAnaBinOffAll, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	    hist2DOnEff[dyn3Str + suffix] ->Fill( passesAnaBinOffAll, caloAlphaTHTDynamic.second, caloAlphaTStandard );
	    


	    // Make candidate L1 trigger plots
	    
	    if ( passesAnaBinOffAll){

	      TString trigStr  = "";
	      bool    trigBool = false;

	      trigStr  = "NoL1";
	      trigBool = passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig1";
	      trigBool = l1Trig1 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig2";
	      trigBool = l1Trig2 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig3";
	      trigBool = l1Trig3 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig4";
	      trigBool = l1Trig4 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig5";
	      trigBool = l1Trig5 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig6";
	      trigBool = l1Trig6 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig7";
	      trigBool = l1Trig7 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      
	    }// End signal region requirement
	    
	    
	    
	  }
	}
	
	
	//	} // Second jet requirement
      }

      hist2D["CaloAlphaT_vs_CaloAlphaTDyn_" + pfJetBinStr]->Fill( caloAlphaTDynamic,          caloAlphaTStandard );	

    }

    // ********************************************************************************
    // *                              PF Jets
    // ********************************************************************************
    if ( pfJetsAboveThresh >= 2 ){
      TString pfJetBinStr     = TString(Form( "%d", pfJetsAboveThresh)) + TString("Jets");
    
      // Check HLT triggers
      if ( caloAlphaTStandard > caloJetAlphaThreshold ){        passesATStandard = true; }
      if ( caloAlphaTDynamic  > caloJetAlphaThreshold ){        passesATDynamic  = true; }
      if ( caloAlphaTHTDynamic.first > caloJetAlphaThreshold ){ passesATDynamic2 = true; }

      if (caloHT > 200){
	passesHTStandard = true;
      }
      if (caloAlphaTHTDynamic.second > caloJetHTThreshold ){
	passesHTDynamic = true;
      }

      // Pf-binned
      hist2D["PFAlphaT_vs_CaloAlphaT_"      + pfJetBinStr]    ->Fill( caloAlphaTStandard,         pfAlphaTStandard );	
      hist2D["PFAlphaT_vs_CaloAlphaTDyn_"   + pfJetBinStr]    ->Fill( caloAlphaTDynamic,          pfAlphaTStandard );	

      hist2DEff["AlphaTStd_vs_HT_" + pfJetBinStr]->Fill( passesHTStandard&&passesATStandard,  pfHT, pfAlphaTStandard );
      hist2DEff["AlphaTDyn_vs_HT_" + pfJetBinStr]->Fill( passesHTStandard&&passesATDynamic,   pfHT, pfAlphaTStandard );
      hist2DEff["AlphaTDy2_vs_HT_" + pfJetBinStr]->Fill( passesHTDynamic &&passesATDynamic,   pfHT, pfAlphaTStandard );
      hist2DEff["AlphaTDy3_vs_HT_" + pfJetBinStr]->Fill( passesHTDynamic &&passesATStandard,  pfHT, pfAlphaTStandard );
      
      if ( caloJetsAboveThresh >= 2 ){

	if ( (*caloJetPT)[1] > 50 ){
	  passesJet2_50 = true;
	  if ( (*caloJetPT)[1] > 60 ){
	    passesJet2_60 = true;
	    if ( (*caloJetPT)[1] > 70 ){
	      passesJet2_70 = true;
	      if ( (*caloJetPT)[1] > 80 ){
		passesJet2_80 = true;
		if ( (*caloJetPT)[1] > 90 ){
		  passesJet2_90 = true;
		  if ( (*caloJetPT)[1] > 100 ){
		    passesJet2_100 = true;
		  }
		}
	      }
	    }
	  }
	}
	if ( ((*caloJetPT)[1] > 100.) || ( (caloHT <= 200) && ((*caloJetPT)[1] >= caloHT/2.)) ){
	  passesDynJet2_100 = true;
	}
      }



      if ( (*pfJetPT)[1] > 100 ){

	hist2DEff["AlphaTStd_vs_HT_Jet2gt50_" + pfJetBinStr]->Fill( passesJet2_50 && passesATStandard, pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTStd_vs_HT_Jet2gt60_" + pfJetBinStr]->Fill( passesJet2_60 && passesATStandard, pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTStd_vs_HT_Jet2gt70_" + pfJetBinStr]->Fill( passesJet2_70 && passesATStandard, pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTStd_vs_HT_Jet2gt80_" + pfJetBinStr]->Fill( passesJet2_80 && passesATStandard, pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTStd_vs_HT_Jet2gt90_" + pfJetBinStr]->Fill( passesJet2_90 && passesATStandard, pfHT, pfAlphaTStandard );

	hist2DEff["AlphaTStd_vs_HT_Jet2gt100_" + pfJetBinStr]->Fill( passesJet2_100 && passesATStandard, pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTDyn_vs_HT_Jet2gt100_" + pfJetBinStr]->Fill( passesJet2_100 && passesATDynamic,  pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTDy2_vs_HT_Jet2gt100_" + pfJetBinStr]->Fill( passesJet2_100 && passesATDynamic,  pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTDy3_vs_HT_Jet2gt100_" + pfJetBinStr]->Fill( passesJet2_100 && passesATStandard, pfHT, pfAlphaTStandard );

	hist2DEff["AlphaTStd_vs_HT_DynJet2gt100_" + pfJetBinStr]->Fill( passesDynJet2_100 && passesATStandard, pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTDyn_vs_HT_DynJet2gt100_" + pfJetBinStr]->Fill( passesDynJet2_100 && passesATDynamic,  pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTDy2_vs_HT_DynJet2gt100_" + pfJetBinStr]->Fill( passesDynJet2_100 && passesATDynamic,  pfHT, pfAlphaTStandard );
	hist2DEff["AlphaTDy3_vs_HT_DynJet2gt100_" + pfJetBinStr]->Fill( passesDynJet2_100 && passesATStandard, pfHT, pfAlphaTStandard );
      }

      // Ensure we are in the alphaT plateau
      if (pfHT > 300){
	hist1DEff["AlphaTStd_TurnOn_" + pfJetBinStr]->Fill( passesHTStandard && passesATStandard, pfAlphaTStandard );
	hist1DEff["AlphaTDyn_TurnOn_" + pfJetBinStr]->Fill( passesHTStandard && passesATDynamic,  pfAlphaTStandard );
	hist1DEff["AlphaTDy2_TurnOn_" + pfJetBinStr]->Fill( passesHTDynamic  && passesATDynamic,   pfAlphaTStandard );
      }

      // Ensure we are in the HT plateau
      if (pfAlphaTStandard > 0.7){
	hist1DEff["HTStd_TurnOn_" + pfJetBinStr]->Fill( passesHTStandard && passesATStandard, pfHT );
	hist1DEff["HTDyn_TurnOn_" + pfJetBinStr]->Fill( passesHTStandard && passesATDynamic,  pfHT );
	hist1DEff["HTDy2_TurnOn_" + pfJetBinStr]->Fill( passesHTDynamic  && passesATDynamic,  pfHT );
	hist1DEff["HTDy3_TurnOn_" + pfJetBinStr]->Fill( passesHTDynamic  && passesATStandard, pfHT );
      }


    } // End 2-PF jet requirement
    

#ifdef TEST
    if ( iEvent > maxEvents ){ break; } // Exit early
#endif


    if ( !(iEvent % 10000) ){ std::cout << "Event " << iEvent << "\n";}

  } // End event loop
#endif





  // ********************************************************************************
  // *                                      Rate                                    *
  // ********************************************************************************


  // ------------------------------------------------------------ 
  // UCT 
  // ------------------------------------------------------------ 

    // NOTE: Currently with HLT emulation the lx1extra GCT products are replaced with UCT so this is a labeling problem
    rateChain->SetBranchAddress("gctCen_Pt",       &uctJetPT);
    rateChain->SetBranchAddress("gctCen_Phi",      &uctJetPhi);
    rateChain->SetBranchAddress("gct_Ht",          &uctHT);
    rateChain->SetBranchAddress("gct_MhtPt",       &uctMHT);
    rateChain->SetBranchAddress("gct_MhtDivHt",    &uctMHToverHT);
    rateChain->SetBranchAddress("gct_MetPt",       &uctMET);


  // ------------------------------------------------------------ 
  // HLT 
  // ------------------------------------------------------------ 
    // // HLT CaloJet
    // rateChain->SetBranchAddress("hltAk4Calo_Pt",        &hltCaloJetPT);
    // rateChain->SetBranchAddress("hltAk4Calo_Px",        &hltCaloJetPx);
    // rateChain->SetBranchAddress("hltAk4Calo_Py",        &hltCaloJetPy);
    // //rateChain->SetBranchAddress("hltAk4Calo_Phi",       &hltCaloJetPhi);
    // rateChain->SetBranchAddress("hltAk4Calo_Eta",       &hltCaloJetEta);
  
    // // HLT PF 
    // rateChain->SetBranchAddress("hltAk4PF_Pt",          &hltPFJetPT); 
    // rateChain->SetBranchAddress("hltAk4PF_Px",          &hltPFJetPx);
    // rateChain->SetBranchAddress("hltAk4PF_Py",          &hltPFJetPy);
    // //rateChain->SetBranchAddress("hltAk4PF_Phi",         &hltPFJetPhi);
    // rateChain->SetBranchAddress("hltAk4PF_Eta",         &hltPFJetEta);
    
    // // // HLT PFNoPU 
    // // rateChain->SetBranchAddress("hltAk4PFNoPU_Pt",      &hltPFNoPUPT); 
    // // rateChain->SetBranchAddress("hltAk4PFNoPU_Px",      &hltPFNoPUPx);
    // // rateChain->SetBranchAddress("hltAk4PFNoPU_Py",      &hltPFNoPUPy);
    // // //rateChain->SetBranchAddress("hltAk4PFNoPU_Phi",     &hltPFNoPUPhi);
    // // rateChain->SetBranchAddress("hltAk4PFNoPU_Eta",     &hltPFNoPUEta);
    
    
  // ------------------------------------------------------------ 
  // GEN
  // ------------------------------------------------------------ 

    // rateChain->SetBranchAddress("genAk4_Pt",         &genJetPT); 
    // rateChain->SetBranchAddress("genAk4_Px",         &genJetPx);
    // rateChain->SetBranchAddress("genAk4_Py",         &genJetPy);
    // //rateChain->SetBranchAddress("genAk4_Phi",      &genJetPhi);
    // rateChain->SetBranchAddress("genAk4_Eta",        &genJetEta);

    // rateChain->SetBranchAddress("genMetTrue_MetPt", &genMET);
  

    // TEMPORARY - REPLACE NAMES FOR MEETING IN 1.5 hours
    rateChain->SetBranchAddress("genAk4_Pt",         &pfJetPT); 
    rateChain->SetBranchAddress("genAk4_Px",         &pfJetPx);
    rateChain->SetBranchAddress("genAk4_Py",         &pfJetPy);
    //rateChain->SetBranchAddress("genAk4_Phi",      &pfJetPhi);
    rateChain->SetBranchAddress("genAk4_Eta",        &pfJetEta);

#ifdef HLT_CALOJET
    // HLT CaloJet
    rateChain->SetBranchAddress("hltAk4Calo_Pt",        &caloJetPT);
    rateChain->SetBranchAddress("hltAk4Calo_Px",        &caloJetPx);
    rateChain->SetBranchAddress("hltAk4Calo_Py",        &caloJetPy);
    //rateChain->SetBranchAddress("hltAk4Calo_Phi",       &caloJetPhi);
    rateChain->SetBranchAddress("hltAk4Calo_Eta",       &caloJetEta);


    rateChain->SetBranchAddress("hltAk4PF_Pt",        &hltPFJetPT);
    rateChain->SetBranchAddress("hltAk4PF_Px",        &hltPFJetPx);
    rateChain->SetBranchAddress("hltAk4PF_Py",        &hltPFJetPy);
    //rateChain->SetBranchAddress("hltAk4PF_Phi",       &hltPFJetPhi);
    rateChain->SetBranchAddress("hltAk4PF_Eta",       &hltPFJetEta);

#else
    // HLT PFJet
    rateChain->SetBranchAddress("hltAk4PF_Pt",        &caloJetPT);
    rateChain->SetBranchAddress("hltAk4PF_Px",        &caloJetPx);
    rateChain->SetBranchAddress("hltAk4PF_Py",        &caloJetPy);
    //rateChain->SetBranchAddress("hltAk4PF_Phi",       &caloJetPhi);
    rateChain->SetBranchAddress("hltAk4PF_Eta",       &caloJetEta);
#endif




  // rateChain->SetBranchAddress("jetPtsPf",        &pfJetPT);
  // //  rateChain->SetBranchAddress("jetPhisPf",       &pfJetPhi);
  // rateChain->SetBranchAddress("jetPxsPf",        &pfJetPx);
  // rateChain->SetBranchAddress("jetPysPf",        &pfJetPy);
  // rateChain->SetBranchAddress("jetEtasPf",       &pfJetEta);
  // //  rateChain->SetBranchAddress("htPf",            &pfHT);
  // rateChain->SetBranchAddress("alphaTPf",        &pfAlphaT);

  // rateChain->SetBranchAddress("jetPtsCalo",      &caloJetPT);
  // rateChain->SetBranchAddress("jetPxsCalo",      &caloJetPx);
  // rateChain->SetBranchAddress("jetPysCalo",      &caloJetPy);
  // rateChain->SetBranchAddress("jetEtasCalo",     &caloJetEta);
  // //  rateChain->SetBranchAddress("jetPhisCalo",     &caloJetPhi);
  // //  rateChain->SetBranchAddress("htCalo",          &caloHT);
  // rateChain->SetBranchAddress("alphaTCalo",      &caloAlphaT);




  //  nEvents = (unsigned int)nuChain->GetEntries();


  dynamicRate dynamicAlphaT    = dynamicRate( HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax );
  dynamicRate dynamicAlphaTHT  = dynamicRate( HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax );
  dynamicRate standardAlphaTHT = dynamicRate( HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax );

#ifdef NEUTRINO

  // Loop over the tree
  //  for ( unsigned int iEvent = 0; iEvent < nEvents; ++iEvent ){
  for ( unsigned int iEvent = nuEventLow; iEvent < nuEventHigh; ++iEvent ){

    rateChain->GetEntry( iEvent );
    


    // ********************************************************************************
    // *                                 UCT triggers                                 *
    // ********************************************************************************

    bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) ); // DJ100
    bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 60) ); // DJ120
    bool l1Trig3 = ( (uctHT >= 160)  || (uctMET >= 70) ); // DJ120
    bool l1Trig4 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 150) && (uctMHToverHT >= 0.17)) );
    bool l1Trig5 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 140) && (uctMHToverHT >= 0.25)) );
    bool l1Trig6 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 130) && (uctMHToverHT >= 0.32)) );
    bool l1Trig7 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 120) && (uctMHToverHT >= 0.36)) );

    // bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) );
    // bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 60) );
    // bool l1Trig3 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHT >= 56)) );
    // bool l1Trig4 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 125) && (uctMHT >= 44)) );
    // bool l1Trig5 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHToverHT >= 0.4)) );
    // bool l1Trig6 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 125) && (uctMHToverHT >= 0.3)) );
    // bool l1Trig7 = ( (uctHT >= 175) );
    
    // OLD NAMING SCHEME
    // bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) );
    // bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 125) && (uctMHToverHT >= 0.3)) );
    // bool l1Trig3 = ( (uctHT >= 200)  || (uctMET >= 60) );
    // bool l1Trig4 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHT >= 56)) );
    // bool l1Trig5 = ( (uctMET >= 70)  || ((uctHT >= 125) && (uctMHT >= 44)) );
    // bool l1Trig6 = ( (uctHT >= 200)  || ((uctHT >= 112) && (uctMHT >= 56)) );


    // Validate rates
    // --------------------------------------------------------------------------------
    
    if ( (uctMET >= 70) ){
      hist1DRate["Trig1_UCTHT_Rate"]->Fill( uctHT );
    }
    if ( (uctMET >= 60) ){
      hist1DRate["Trig2_UCTHT_Rate"]->Fill( uctHT );
    }
    if ( ((uctHT >= 112) && (uctMHT >= 56)) ){
      hist1DRate["Trig3_UCTHT_Rate"]->Fill( uctHT );
    }
    if ( (uctMET >= 70) || ((uctHT >= 125) && (uctMHT >= 44)) ){
      hist1DRate["Trig4_UCTHT_Rate"]->Fill( uctHT );
    }
    if ( ((uctHT >= 112) && (uctMHToverHT >= 0.4)) ){
      hist1DRate["Trig5_UCTHT_Rate"]->Fill( uctHT );
    }
    if ( (uctMET >= 70) || ((uctHT >= 125) && (uctMHToverHT >= 0.3)) ){
      hist1DRate["Trig6_UCTHT_Rate"]->Fill( uctHT );
    }

    // Double and quad jet PT
    float sJet(0), dJet(0), qJet(0);
    if ( uctJetPT->size() >= 1){
      sJet = (*uctJetPT)[0];
      if ( uctJetPT->size() > 1){
	dJet = (*uctJetPT)[1];
	if ( uctJetPT->size() > 3){
	  qJet = (*uctJetPT)[3];
	}
      }
    }
    if ( (uctHT >= 175) || (dJet >= 100.) || (qJet >= 60) ){
	hist1DRate["TrigX_UCTHT_Rate"]->Fill( uctMET );
    }
    if ( (uctHT >= 200) || (dJet >= 120.) || (qJet >= 60) ){
	hist1DRate["TrigY_UCTHT_Rate"]->Fill( uctMET );
    }
    hist1DRate["SingleJet_Rate"]->Fill( sJet );
    hist1DRate["DoubleJet_Rate"]->Fill( dJet );
    hist1DRate["QuadJet_Rate"]  ->Fill( qJet );




    hist2DRate["UCTMHT_vs_UCTHT_Rate"] ->Fill( uctHT, uctMHT );
    hist2DRate["UCTMET_vs_UCTHT_Rate"] ->Fill( uctHT, uctMET ); 
    if (uctMHToverHT > 0.3){
      hist2DRate["UCTMET_vs_UCTHT_MHTdivHT0p3_Rate"]->Fill( uctHT, uctMET );
    }


    // Validate rates - HLT
    // --------------------------------------------------------------------------------
    
    double HLTcaloHT(0), HLTpfHT(0);
    for(uint iJet = 0;iJet < caloJetPT->size(); ++iJet ){
      if( (*caloJetPT)[iJet] < 40. ) break;
      HLTcaloHT    += (*caloJetPT)[iJet];
    }
    for(uint iJet = 0;iJet < hltPFJetPT->size(); ++iJet ){
      if( (*hltPFJetPT)[iJet] < 40. ) break;
      HLTpfHT    += (*hltPFJetPT)[iJet];
    }

    histHLTRate["HLTCaloHTRate"] ->Fill( HLTcaloHT );
    
    if (uctHT >= 200.){
      histHLTRate["L1HTT200_HLTCaloHTRate"]->Fill( HLTcaloHT );
    }
    if ( (uctHT >= 200.) && (HLTcaloHT >= 200.) ){
      histHLTRate["L1HTT200_HLTCaloHT200_HLTPFHTRate"]->Fill(HLTpfHT);
    }





    float caloAlphaTStandard(0);
    float caloAlphaTDynamic(0);
    float pfAlphaTDynamic(0);
    std::pair<float,float> caloAlphaTHTDynamic;

    
    // calculate HT for new jet threshold
    caloHT = 0;
    float caloHTJet2Cut(0);
    float caloMHTx(0), caloMHTy(0), caloMHT(0);
    
    // Calojets
    int caloJetsAboveThresh(0);
    for(uint iJet = 0;iJet < caloJetPT->size(); ++iJet ){
      if( (*caloJetPT)[iJet] < caloJetThreshold ) break;
      caloHT    += (*caloJetPT)[iJet];

      if ( (caloJetPT->size() > 1) && ((*caloJetPT)[iJet] > secondJetThreshold) ){
	caloMHTx  += (*caloJetPx)[iJet];
	caloMHTy  += (*caloJetPy)[iJet];
	caloHTJet2Cut += (*caloJetPT)[iJet];
      }

      caloJetsAboveThresh++;
    }
    if ( caloJetsAboveThresh > MAX_JETS ){
      caloJetsAboveThresh = MAX_JETS;
    }
    TString nCaloJetsStr = Form("%d", caloJetsAboveThresh );
    caloMHT = sqrt( caloMHTx*caloMHTx + caloMHTy*caloMHTy ); 

    // PF jets
    pfHT = 0;
    int pfJetsAboveThresh(0);
    for(uint iJet = 0;iJet < pfJetPT->size(); ++iJet ){
      if( (*pfJetPT)[iJet] < pfJetThreshold ) break;
      pfHT += (*pfJetPT)[iJet];
      pfJetsAboveThresh++;
    }
    if ( pfJetsAboveThresh > MAX_JETS ){
      pfJetsAboveThresh = MAX_JETS;
    }
    TString nPFJetsStr = Form("%d", pfJetsAboveThresh );


    // ********************************************************************************
    // *                                      MHT                                     *
    // ********************************************************************************


    if ( caloJetsAboveThresh >= 2 ){
      if ( (*caloJetPT)[1] > secondJetThreshold ){
	hist2DRate["MHT_vs_HT_Rate"]->Fill( caloHTJet2Cut, caloMHT );
      }
    }

    // ********************************************************************************
    // *                                     AlphaT                                   *
    // ********************************************************************************
    caloAlphaTStandard  = calculateAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold );
    caloAlphaTDynamic   = calculateDynamicAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold,
						  caloJetAlphaThreshold );
    

    pfAlphaTDynamic   = calculateDynamicAlphaT( pfJetPT, pfJetPx, pfJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold,
						  caloJetAlphaThreshold );


    // Dynamic HT and AlphaT                                                                                                              
    caloAlphaTHTDynamic = calculateDynamicAlphaTHT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold,
						    caloJetAlphaThreshold, caloJetHTThreshold );



    // ********************************************************************************
    // *                              Calojets 
    // ********************************************************************************
    if ( caloJetsAboveThresh >= 2 ){
      TString caloJetBinStr = TString(Form( "%d", caloJetsAboveThresh)) + TString("Jets");

      // Vary the online second jet PT cut
      for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){

	float jet2PTCut = jet2PTCuts[ iJet2Cut ];

	if ((*caloJetPT)[1] > jet2PTCut){
	TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
	TString stdStr  = "OnRate_AlphaTStd_vs_HT_"  + jet2PTCutStr;
	TString dynStr  = "OnRate_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
	TString dyn2Str = "OnRate_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
	TString dyn3Str = "OnRate_AlphaTDyn3_vs_HT_" + jet2PTCutStr;

	// Inclusive 
	hist2DRate[stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	hist2DRate[dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	hist2DRate[dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	hist2DRate[dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	// Jet binned
	hist2DRate[stdStr  + "_" + caloJetBinStr]->Fill( caloHT,                     caloAlphaTStandard );
	hist2DRate[dynStr  + "_" + caloJetBinStr]->Fill( caloHT,                     caloAlphaTDynamic );
	hist2DRate[dyn2Str + "_" + caloJetBinStr]->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	hist2DRate[dyn3Str + "_" + caloJetBinStr]->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard ); 
	




	// Make candidate L1 trigger plots
	TString trigStr = "";
	
	trigStr = "NoL1";
	hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );

	if ( l1Trig1){
	  trigStr = "Trig1";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}
	if ( l1Trig2){
	  trigStr = "Trig2";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}
	if ( l1Trig3){
	  trigStr = "Trig3";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}
	if ( l1Trig4){
	  trigStr = "Trig4";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}
	if ( l1Trig5){
	  trigStr = "Trig5";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}
	if ( l1Trig6){
	  trigStr = "Trig6";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}
	if ( l1Trig7){
	  trigStr = "Trig7";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	}






	}
      }

      //      hist2D["CaloAlphaT_vs_CaloAlphaTDyn_" + caloJetBinStr]->Fill( caloAlphaTDynamic,          caloAlphaTStandard );	

    }

    if ( caloJetsAboveThresh >= 2 ){  
      // Fill standard rate
      hist2DRate["AlphaTStd_vs_HT_Rate"] ->Fill( caloHT, caloAlphaTStandard );

      // Second jet threshold

//       if ( (*caloJetPT)[1] > 50. ){
// 	hist2DRate["AlphaTStd_vs_HT_Jet2PT50_Rate"]->Fill( caloHT, caloAlphaTStandard );
// 	if ( (*caloJetPT)[1] > 60. ){
// 	  hist2DRate["AlphaTStd_vs_HT_Jet2PT60_Rate"]->Fill( caloHT, caloAlphaTStandard );
// 	  if ( (*caloJetPT)[1] > 70. ){
// 	    hist2DRate["AlphaTStd_vs_HT_Jet2PT70_Rate"]->Fill( caloHT, caloAlphaTStandard );
// 	    if ( (*caloJetPT)[1] > 80. ){
// 	      hist2DRate["AlphaTStd_vs_HT_Jet2PT80_Rate"]->Fill( caloHT, caloAlphaTStandard );
// 	      if ( (*caloJetPT)[1] > 90. ){
// 		hist2DRate["AlphaTStd_vs_HT_Jet2PT90_Rate"]->Fill( caloHT, caloAlphaTStandard );
// 		if ( (*caloJetPT)[1] > 100. ){
// 		  hist2DRate["AlphaTStd_vs_HT_Jet2PT100_Rate"]->Fill( caloHT, caloAlphaTStandard );
// 		}
// 	      }
// 	    }
// 	  }
// 	}
//       }


      // Dynamic second jet threshold
      if ( ((*caloJetPT)[1] > 100.) || ( (caloHT <= 200) && ((*caloJetPT)[1] >= caloHT/2.)) ){
	hist2DRate["AlphaTStd_vs_HT_DynJet2PT_Rate"]->Fill( caloHT, caloAlphaTStandard );
      }
    
    }



    // Monojets
    if (caloJetPT->size() == 1){
      hist1DRate["MonoJetPT_Rate"]->Fill((*caloJetPT)[0]);
    }
    else if( (*caloJetPT)[1] < caloJetThreshold ){
      hist1DRate["MonoJetPT_Rate"]->Fill((*caloJetPT)[0]);
    }





//     // Calo HLT cuts
//     bool passesHLTStandard(false), passesHLTDynamic(false), passesHLTDynamic2(false), passesHLTDynamic3(false);
//     bool passesJet2_100(false);    



//     if ( caloJetsAboveThresh >= 2 ){

//       // Check HLT triggers 
//       if (caloHT > 200){
//         if ( caloAlphaTStandard > caloJetAlphaThreshold ){passesHLTStandard = true; }
//         if ( caloAlphaTDynamic  > caloJetAlphaThreshold ){passesHLTDynamic  = true; }
//       }
//       if (caloAlphaTHTDynamic.second > caloJetHTThreshold ){
//         if ( caloAlphaTHTDynamic.first > caloJetAlphaThreshold ){ passesHLTDynamic2 = true; }
//         if ( caloAlphaT                > caloJetAlphaThreshold ){ passesHLTDynamic3 = true; }
//       }
//       if ( (*caloJetPT)[1] > 100 ){
// 	passesJet2_100 = true;
//       }
//     }

//     hist2DRate["AlphaTStd_vs_HT_Jet2gt100_Rate"]->Fill( passesJet2_100 && passesHLTStandard, pfHT, pfAlphaTStandard );
//     hist2DRate["AlphaTDyn_vs_HT_Jet2gt100_Rate"]->Fill( passesJet2_100 && passesHLTDynamic,  pfHT, pfAlphaTStandard );
//     hist2DRate["AlphaTDy2_vs_HT_Jet2gt100_Rate"]->Fill( passesJet2_100 && passesHLTDynamic2, pfHT, pfAlphaTStandard );
//     hist2DRate["AlphaTDy3_vs_HT_Jet2gt100_Rate"]->Fill( passesJet2_100 && passesHLTDynamic3, pfHT, pfAlphaTStandard );








    // Rates as function of lead and second jet thresholds
    if (caloJetsAboveThresh > 1){
      float HT12 = (*caloJetPT)[0] + (*caloJetPT)[1];
      
      if ( caloHT > 200 ){

	// Standard alphaT, HT
	for (uint iCut = 0; iCut < alphaTCuts.size(); ++iCut){
	  float alphaTCut = alphaTCuts[ iCut ];
	  TString aTStr = Form("%1.2f", alphaTCut );
	  aTStr = aTStr.ReplaceAll(".","p");

	  if ( caloAlphaTStandard > alphaTCut ){
	    hist2DRate["Jet2PT_vs_Jet1PT_alphaTgt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[1]);
	    hist2DRate[  "HT12_vs_Jet1PT_alphaTgt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], HT12 );

	    // Jet binned
	    hist2DRate["Jet2PT_vs_Jet1PT_alphaTgt" + aTStr + "_" + nCaloJetsStr + "Jets_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[1]);
	    hist2DRate[  "HT12_vs_Jet1PT_alphaTgt" + aTStr + "_" + nCaloJetsStr + "Jets_Rate"]->Fill((*caloJetPT)[0], HT12 );

	    // WW scattering
	    if (caloJetsAboveThresh > 3){
	      hist2DRate["Jet4PT_vs_Jet1PT_alphaTgt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[3]);
	    }
	  }
	  else{	break; }
	}

	// Dynamic alphaT
	for (uint iCut = 0; iCut < alphaTCuts.size(); ++iCut){
	  float alphaTCut = alphaTCuts[ iCut ];
	  TString aTStr = Form("%1.2f", alphaTCut );
	  aTStr = aTStr.ReplaceAll(".","p");

	  if ( caloAlphaTDynamic > alphaTCut ){
	    hist2DRate["Jet2PT_vs_Jet1PT_alphaTDyngt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[1]);
	    hist2DRate[  "HT12_vs_Jet1PT_alphaTDyngt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], HT12 );

	    // Jet binned
	    hist2DRate["Jet2PT_vs_Jet1PT_alphaTDyngt" + aTStr + "_" + nCaloJetsStr + "Jets_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[1]);
	    hist2DRate[  "HT12_vs_Jet1PT_alphaTDyngt" + aTStr + "_" + nCaloJetsStr + "Jets_Rate"]->Fill((*caloJetPT)[0], HT12 );

	    // WW scattering
	    if (caloJetsAboveThresh > 3){
	      hist2DRate["Jet4PT_vs_Jet1PT_alphaTDyngt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[3]);
	    }
	  }
	  else{	break; }
	}

      } // HT requirement

      if (caloAlphaTHTDynamic.second > 200){

	// Dynamic alphaT, HT
	for (uint iCut = 0; iCut < alphaTCuts.size(); ++iCut){
	  float alphaTCut = alphaTCuts[ iCut ];
	  TString aTStr = Form("%1.2f", alphaTCut );
	  aTStr = aTStr.ReplaceAll(".","p");

	  if ( caloAlphaTHTDynamic.first > alphaTCut ){
	    hist2DRate["Jet2PT_vs_Jet1PT_alphaTDyn2gt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[1]);
	    hist2DRate[  "HT12_vs_Jet1PT_alphaTDyn2gt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], HT12 );

	    // Jet binned
	    hist2DRate["Jet2PT_vs_Jet1PT_alphaTDyn2gt" + aTStr + "_" + nCaloJetsStr + "Jets_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[1]);
	    hist2DRate[  "HT12_vs_Jet1PT_alphaTDyn2gt" + aTStr + "_" + nCaloJetsStr + "Jets_Rate"]->Fill((*caloJetPT)[0], HT12 );

	    // WW scattering
	    if (caloJetsAboveThresh > 3){
	      hist2DRate["Jet4PT_vs_Jet1PT_alphaTDyn2gt" + aTStr + "_Rate"]->Fill((*caloJetPT)[0], (*caloJetPT)[3]);
	    }
	  }
	  else{	break; }
	}

      }
    } // End two jet requirement
  


    // Rates as function of PF lead and second jet thresholds
    if (pfJetsAboveThresh > 1){
      float HT12 = (*pfJetPT)[0] + (*pfJetPT)[1];

      if ( pfHT > 200 ){

	float pfAlphaTStandard = calculateAlphaT( pfJetPT, pfJetPx, pfJetPy, pfJetThreshold );

	// Standard alphaT, HT
	for (uint iCut = 0; iCut < alphaTCuts.size(); ++iCut){
	  float alphaTCut = alphaTCuts[ iCut ];
	  TString aTStr = Form("%1.2f", alphaTCut );
	  aTStr = aTStr.ReplaceAll(".","p");

	  if ( pfAlphaTStandard > alphaTCut ){
	    hist2DRate["PFJet2PT_vs_Jet1PT_alphaTgt" + aTStr + "_Rate"]->Fill((*pfJetPT)[0], (*pfJetPT)[1]);
	    hist2DRate[  "PFHT12_vs_Jet1PT_alphaTgt" + aTStr + "_Rate"]->Fill((*pfJetPT)[0], HT12 );

	    // Jet binned
	    hist2DRate["PFJet2PT_vs_Jet1PT_alphaTgt" + aTStr + "_" + nPFJetsStr + "Jets_Rate"]->Fill((*pfJetPT)[0], (*pfJetPT)[1]);
	    hist2DRate[  "PFHT12_vs_Jet1PT_alphaTgt" + aTStr + "_" + nPFJetsStr + "Jets_Rate"]->Fill((*pfJetPT)[0], HT12 );
	  }
	  else{	break; }
	}

	// Dynamic alphaT
	for (uint iCut = 0; iCut < alphaTCuts.size(); ++iCut){
	  float alphaTCut = alphaTCuts[ iCut ];
	  TString aTStr = Form("%1.2f", alphaTCut );
	  aTStr = aTStr.ReplaceAll(".","p");

	  if ( pfAlphaTDynamic > alphaTCut ){
	    hist2DRate["PFJet2PT_vs_Jet1PT_alphaTDyngt" + aTStr + "_Rate"]->Fill((*pfJetPT)[0], (*pfJetPT)[1]);
	    hist2DRate[  "PFHT12_vs_Jet1PT_alphaTDyngt" + aTStr + "_Rate"]->Fill((*pfJetPT)[0], HT12 );
	  
	    // Jet binned - (Online)
	    hist2DRate["PFJet2PT_vs_Jet1PT_alphaTDyngt" + aTStr + "_" + nPFJetsStr + "Jets_Rate"]->Fill((*pfJetPT)[0], (*pfJetPT)[1]);
	    hist2DRate[  "PFHT12_vs_Jet1PT_alphaTDyngt" + aTStr + "_" + nPFJetsStr + "Jets_Rate"]->Fill((*pfJetPT)[0], HT12 );

	  }
	  else{	break; }
	}
	
	
      } // HT requirement
    }


    // calculate all possible alphaT, HT combinations
    std::vector<std::pair<float,float> > alphaTHTPairs;
    alphaTHTPairs = calculateDynamicAlphaTPairs(caloJetPT, caloJetPx, caloJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold );


    standardAlphaTHT.triggerFired( caloHT, caloAlphaTStandard );


    for (uint iPair = 0;iPair < alphaTHTPairs.size();++iPair){

      //std::cout << alphaTHTPairs[iPair].first << "\t" << alphaTHTPairs[iPair].second << "\n";

      // x = HT, y = AlphaT
      dynamicAlphaT.triggerFired(   caloHT,                      alphaTHTPairs[iPair].first );
      dynamicAlphaTHT.triggerFired( alphaTHTPairs[iPair].second, alphaTHTPairs[iPair].first );
    }
    //    std::cout << "\n";

    
    // Calculate trigger rate for event
    standardAlphaTHT.endEvent();
    dynamicAlphaT.endEvent();
    dynamicAlphaTHT.endEvent();







#ifdef TEST
    if ( iEvent > maxEvents ){ break; } // Exit early
#endif

    if ( !(iEvent % 10000) ){ std::cout << "Event " << iEvent << "\n";}

  }
#endif



  // ********************************************************************************
  // *                               Store histograms                               *
  // ********************************************************************************
  
  gStyle->SetPaintTextFormat("9.1f");

  // fOut->mkdir("Triggers");
  // fOut->mkdir("Triggers/Efficiency");
  // fOut->mkdir("Triggers/Rate");
  // fOut->mkdir("Triggers/Correlations");
  fOut->mkdir("Raw");
  fOut->mkdir("Raw/Efficiency");
  fOut->mkdir("Raw/HLTEfficiency");
  fOut->mkdir("Raw/HLTRate");
  fOut->mkdir("Raw/Rate");
  fOut->mkdir("Raw/Correlations");
  fOut->mkdir("Raw/Distributions");


  // Dynamic rates
  std::vector< std::vector< int > > stdAlphaTHTRate = standardAlphaTHT.cumulativeSum;
  std::vector< std::vector< int > > dynAlphaTRate   = dynamicAlphaT.cumulativeSum;
  std::vector< std::vector< int > > dynAlphaTHTRate = dynamicAlphaTHT.cumulativeSum;
  
  for (int xBin = 0; xBin < dynamicAlphaT.xBins; ++xBin ){
    for (int yBin = 0; yBin < dynamicAlphaT.yBins; ++yBin ){
      
      float rateStdAlphaTHT = rateScaleFactor*stdAlphaTHTRate[xBin][yBin];
      float rateDynAlphaT   = rateScaleFactor*  dynAlphaTRate[xBin][yBin];
      float rateDynAlphaTHT = rateScaleFactor*dynAlphaTHTRate[xBin][yBin];

      AlphaTStd_vs_HT_Rate->SetBinContent( xBin+1, yBin+1  , rateStdAlphaTHT );
      AlphaTDyn_vs_HT_Rate->SetBinContent( xBin+1, yBin+1  , rateDynAlphaT );
      AlphaTDy2_vs_HT_Rate->SetBinContent( xBin+1, yBin+1  , rateDynAlphaTHT );
      
    }
  }
  fOut->cd("Raw/Rate");
  AlphaTStd_vs_HT_Rate->Write();
  AlphaTDyn_vs_HT_Rate->Write();
  AlphaTDy2_vs_HT_Rate->Write();
  // fOut->cd("Triggers/Rate");
  
  // TCanvas *canv = new TCanvas();
  // canv->SetGridx();
  // canv->SetGridy();
  // canv->SetName("AlphaTStd_vs_HT_Rate_Check");
  // AlphaTStd_vs_HT_Rate->Draw("COLZTEXT");
  // canv->Write();
  // canv->SetName("AlphaTDyn_vs_HT_Rate");
  // AlphaTDyn_vs_HT_Rate->Draw("COLZTEXT");
  // canv->Write();
  // canv->SetName("AlphaTDy2_vs_HT_Rate");
  // AlphaTDy2_vs_HT_Rate->Draw("COLZTEXT");
  // canv->Write();

  // Store distributions
  for(std::map<TString, TH1*>::const_iterator itr2 = hist1D.begin(); itr2 != hist1D.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key                                                                            
    fOut->cd("Raw/Distributions");
    itr2->second->Write();
  }



  for(std::map<TString, TEfficiency*>::const_iterator itr2 = hist1DEff.begin(); itr2 != hist1DEff.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key 
    fOut->cd("Raw/Efficiency");
    itr2->second->Write();
//     fOut->cd("Triggers/Efficiency");

//     TCanvas *canv2 = new TCanvas(histoName);
//     canv2->SetGridx();
//     canv2->SetGridy();
//      TEfficiency *temp = (TEfficiency*)itr2->second->Clone();
//      temp->Draw();
// //     gPad->Update();
// //     temp->GetPaintedHistogram()->SetMaximum(1.);
// //     temp->Draw();
//     canv2->Write();
   
  }


  for(std::map<TString, TEfficiency*>::const_iterator itr2 = hist2DEff.begin(); itr2 != hist2DEff.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key 
    fOut->cd("Raw/Efficiency");
    itr2->second->Write();
//     fOut->cd("Triggers/Efficiency");

//     TCanvas *canv2 = new TCanvas(histoName);
//     canv2->SetGridx();
//     canv2->SetGridy();
//     TEfficiency *temp = (TEfficiency*)itr2->second->Clone();
//     temp->Draw("COLZTEXT");
//     gPad->Update();
//     temp->GetPaintedHistogram()->SetMaximum(1.);
//     temp->Draw("COLZTEXT");
// //     itr2->second->Draw("COLZTEXT");
// //     itr2->second->GetPaintedHistogram()->SetMaximum(1.);
// //     itr2->second->SetMinimum(0.);
// //     itr2->second->SetMaximum(1.);
// //    gPad->Update();
//     //itr2->second->GetPaintedGraph()->GetZaxis()->SetMaximum(1.);
//     //    itr2->second->GetPaintedGraph()->SetMaximum(1.);
//     //    itr2->second->Draw("COLZTEXT");
//     canv2->Write();
   
  }



  for(std::map<TString, TEfficiency*>::const_iterator itr2 = hist2DOnEff.begin(); itr2 != hist2DOnEff.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key 
    fOut->cd("Raw/Efficiency");

    itr2->second->Write();

    TEfficiency *eff      = (TEfficiency*)itr2->second     ->Clone();
    TH2* passed           = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* total            = (TH2*)eff->GetTotalHistogram() ->Clone();

    //    TEfficiency *effCumul = (TEfficiency*)eff->Clone();
    TH2* passedCumul      = (TH2*)eff->GetPassedHistogram()->Clone();
//     TH2* totalCumul       = (TH2*)effCumul->GetTotalHistogram() ->Clone();

    // Uniform
    TEfficiency *effUniCumul = (TEfficiency*)eff->Clone();
    TH2* totalUniCumul       = (TH2*)eff->GetTotalHistogram() ->Clone();

    //    reverseCumulative2D( total,  totalCumul,  1 );
    reverseCumulative2D( passed, passedCumul, 1 );
    // 
    //    fillUniform2D( totalUniCumul, passed->GetEntries() );

    // CHANGING TO NOT BE UNIFORM ANYMORE
    fillUniform2D( totalUniCumul, total->GetEntries() );


    // Extract jet2PT cut
    TString ptCutStr = histoName;
    ptCutStr = ptCutStr.Remove(0, ptCutStr.Index("Jet2gt") );
    ptCutStr.Remove( ptCutStr.Index("_"), ptCutStr.Length() );
    ptCutStr.ReplaceAll("Jet2gt","");
    float jet2PTCut = ptCutStr.Atof();
    // Cut unphysical region ht < 2*jet2PT
     clearRectangleX( passedCumul,   jet2PTCut*2 );
//     clearRectangleX( totalUniCumul, jet2PTCut*2 );



    // Have to set the TEfficiency in this order, otherwise it will reject the passed histogram and still set the new total histogram 
//     effCumul->SetTotalHistogram(  *totalCumul,  "" );
//     effCumul->SetPassedHistogram( *passedCumul, "" );
//     effCumul->SetName( histoName + "_Cumul" );
//     effCumul->Write();

    effUniCumul->SetTotalHistogram(  *totalUniCumul, "" );
    effUniCumul->SetPassedHistogram( *passedCumul,   "" );
    effUniCumul->SetName( histoName + "_UniCumul" );
    effUniCumul->Write();


//     fOut->cd("Triggers/Efficiency");
//     TCanvas *canv2 = new TCanvas(histoName);
//     canv2->SetGridx();
//     canv2->SetGridy();
//     TEfficiency *temp = (TEfficiency*)effCumul->Clone();
//     temp->Draw("COLZTEXT");
//     gPad->Update();
//     temp->GetPaintedHistogram()->SetMaximum(1.);
//     temp->Draw("COLZTEXT");
//     canv2->Write();
   
  }


  for(std::map<TString, TH1*>::const_iterator itr = hist1DRate.begin(); itr != hist1DRate.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/Rate");
    itr->second->Write();

    TH1* cumulHist = (TH1D*)itr->second->Clone();
    cumulHist->SetName( histoName + "_Cumul" );
    reverseCumulative( itr->second, cumulHist, rateScaleFactor );

    cumulHist->Write();

    // fOut->cd("Triggers/Rate");
    // TCanvas *canv2 = new TCanvas(histoName);
    // canv2->SetGridx();
    // canv2->SetGridy();
    // canv2->SetLogy();
    // cumulHist->Draw();
    // canv2->Write();
   
  }

  for(std::map<TString, TH2*>::const_iterator itr = hist2DRate.begin(); itr != hist2DRate.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/Rate");
    itr->second->Write();

    TH2* cumulHist = (TH2F*)itr->second->Clone();
    cumulHist->SetName( histoName + "_Cumul" );
    reverseCumulative2D( itr->second, cumulHist, rateScaleFactor );
    // Remove unphysical space
    if (histoName.Contains("HT12_vs_Jet1PT")){
      clearTriangle( cumulHist, -1, caloJetThreshold );
    }
    else if(histoName.Contains("Jet2PT_vs_Jet1PT")){
      clearTriangle( cumulHist, +1, 0 );
    }
    else if(histoName.Contains("MHT_vs_HT_Rate")){
      clearTriangle( cumulHist, +1, 0 );
      //      TLine *line = new TLine( secondJetThreshold*2, 0, secondJetThreshold*2, 250 );
    }

    // Cut unphysical region ht < 2*jet2PT
    if (histoName.Contains("OnRate_")){
      // Extract jet2PT cut
      TString ptCutStr = histoName;
      ptCutStr = ptCutStr.Remove(0, ptCutStr.Index("Jet2gt") );
      ptCutStr.Remove( ptCutStr.Index("_"), ptCutStr.Length() );
      ptCutStr.ReplaceAll("Jet2gt","");
      float jet2PTCut = ptCutStr.Atof();
      clearRectangleX( cumulHist, jet2PTCut*2 );
    }

    cumulHist->Write();


    // fOut->cd("Triggers/Rate");
    // TCanvas *canv2 = new TCanvas(histoName);
    // canv2->SetGridx();
    // canv2->SetGridy();
    // cumulHist->Draw("COLZ");
    // canv2->Write();
   
  }



  for(std::map<TString, TH2*>::const_iterator itr2 = hist2D.begin(); itr2 != hist2D.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key 
    fOut->cd("Raw/Correlations");
    itr2->second->Write();
    // fOut->cd("Triggers/Correlations");

    // TCanvas *canv2 = new TCanvas(histoName);
    // canv2->SetGridx();
    // canv2->SetGridy();
    // itr2->second->Draw("COLZTEXT");
    // canv2->Write();
   
  }

  // ********************************************************************************
  // *                                     HLT                                      *
  // ********************************************************************************
  for(std::map<TString, TH1*>::const_iterator itr = histHLTRate.begin(); itr != histHLTRate.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/HLTRate");
    itr->second->Write();

    TH2* cumulHist = (TH2F*)itr->second->Clone();
    cumulHist->SetName( histoName + "_Cumul" );
    reverseCumulative( itr->second, cumulHist, rateScaleFactor );
    cumulHist->Write();
   
  }


  for(std::map<TString, TEfficiency*>::const_iterator itr2 = hist2DHLTEff.begin(); itr2 != hist2DHLTEff.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key 
    fOut->cd("Raw/HLTEfficiency");
    itr2->second->Write();

    TH2* cumulHist = (TH2F*)itr2->second->Clone();
    cumulHist->SetName( histoName + "_Cumul" );


    // Uniform
    TEfficiency *effUniCumul = (TEfficiency*)itr2->second->Clone();

    TH2* passed              = (TH2*)itr2->second->GetPassedHistogram()->Clone();
    TH2* passedYCumul        = (TH2*)itr2->second->GetPassedHistogram()->Clone();
    TH2* total               = (TH2*)itr2->second->GetTotalHistogram() ->Clone();
    TH2* totalYCumul         = (TH2*)itr2->second->GetTotalHistogram() ->Clone();

    reverseCumulativeY( passed, passedYCumul, 1. );
    fillUniformY( total, totalYCumul );

    effUniCumul->SetTotalHistogram(  *totalYCumul, "" );
    effUniCumul->SetPassedHistogram( *passedYCumul,   "" );
    effUniCumul->SetName( histoName + "_UniCumul" );
    effUniCumul->Write();

  }


//   // ********************************************************************************
//   // *                               Store histograms                               *
//   // ********************************************************************************
  

//   fOut->mkdir("Efficiency");
//   fOut->mkdir("Rate");
//   fOut->cd("Rate");


//   std::map<TString, TH2*>::const_iterator itr2;
//   for(itr2 = hist2D.begin(); itr2 != hist2D.end(); ++itr2){

//     TString histoName = itr2->first; //Extract the histogram key 

//     hist2D[ histoName ]->Write();

// 	reverseCumulative2D( hist2D[ histoName ], hist2D[ histoName ], rateScaleFactor );
// 	TString seedLabel = histoName;
// 	seedLabel.ReplaceAll("NewRate_","");
// 	//	std::cout << "SeedLabel = " << seedLabel << "\n";
// 	bestCutForRate[seedLabel] = getEfficiencyBins( hist2D[ histoName ], maxRate );
// 	std::cout << "Getting efficiency from: " << histoName << "\n";


    
//   }




  exit(0);





  fOut->Close();
}
//Clears region above(+) or below(-) diagonal
void clearTriangle( TH1* histogram, int sign, float extra ){

  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();

  // Calculate the extra strip needed to be cut
  float yBinWidth = histogram->GetXaxis()->GetBinWidth(0);
  int   extraBins = extra/yBinWidth;

  if (sign > 0){
    for ( int iBinX = 1; iBinX <= nBinsX; ++iBinX ){     
      for ( int iBinY = iBinX + 1 - extraBins; iBinY < nBinsY; ++iBinY ){     
	histogram->SetBinContent( iBinX, iBinY, 0 );
	//	std::cout << "1\n";
      }
    }
  }
  else{
    for ( int iBinX = 1; iBinX <= nBinsX; ++iBinX ){     
      for ( int iBinY = iBinX - 1 + extraBins; iBinY > 0; --iBinY ){     
	if (iBinY <= 0){ continue;}
	histogram->SetBinContent( iBinX, iBinY, 0 );
	//	std::cout << iBinX << "\t" << iBinY <<  " \n";
      }
    }
  }

}


void clearRectangleX( TH1* histogram, float cut ){

  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  int xExtremeBin     = 0;

  // Find highest xBin below cut
  for ( int iBinX = 1; iBinX <= nBinsX; ++iBinX ){     

    float lowEdge = histogram->GetXaxis()->GetBinLowEdge(iBinX);
    if ( lowEdge >= cut ){
      xExtremeBin = iBinX - 1;
      break;
    }
  }

  // Remove rectangular region
  for ( int iBinX = 0; iBinX <= xExtremeBin; ++iBinX ){     
    for ( int iBinY = 0; iBinY <= nBinsY; ++iBinY ){     
      histogram->SetBinContent( iBinX, iBinY, 0 );
    }
  }

}


// Creates a normalised, reverse cumaltive of the input histogram 'histogram', returning by modifying 'rCumulHist'. Scales result by specified 'scale'.
void reverseCumulative( TH1* histogram, TH1* rCumulHist, double scale){

//   int nBins          = histogram->GetNbinsX();
//   double cumulBin    = histogram->GetBinContent( nBins + 1 ); // Initialise with overflow bin                         
//   //  double scaleFactor = scale/histogram->GetEntries();

//   for ( int iBin = nBins; iBin > 0; --iBin ){                                                                   
//     cumulBin += histogram->GetBinContent( iBin );                                                               
//     rCumulHist->SetBinContent( iBin, cumulBin );                                                                
//   }

// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin

  int nBinsX          = histogram->GetNbinsX();
  int xBinOver        = nBinsX + 1;

  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){     
      double integral = histogram->Integral( iBinX, xBinOver );
      rCumulHist->SetBinContent( iBinX, integral );
  }      
  
  rCumulHist->Scale( scale );
}


void reverseCumulativeY( TH2* histogram, TH2* rCumulHist, double scale){
// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin

  if ( histogram->GetEntries() == 0 ){    return;  }

  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  int xBinOver        = nBinsX + 1;
  int yBinOver        = nBinsY + 1;

  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){     
    double integral = 0;
    for ( int iBinY = nBinsY; iBinY > 0; --iBinY ){     
      integral = histogram->Integral( iBinX, iBinX, iBinY, yBinOver ); //      integral = histogram->Integral( iBinX, iBinX+1, iBinY, yBinOver );
      rCumulHist->SetBinContent( iBinX, iBinY, integral );
    }
  }      
  
  rCumulHist->Scale( scale );
}


void fillUniformY( TH2* histogram, TH2* rCumulHist){

  if ( histogram->GetEntries() == 0 ){    return;  }

  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  int xBinOver        = nBinsX + 1;
  int yBinOver        = nBinsY + 1;

  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){     
    double integral = histogram->Integral( iBinX, iBinX, 1, yBinOver ); //     double integral = histogram->Integral( iBinX, iBinX+1, iBinY, yBinOver );
      for ( int iBinY = nBinsY; iBinY > 0; --iBinY ){     
	rCumulHist->SetBinContent( iBinX, iBinY, integral );
      }
  }      

}


void reverseCumulative2D( TH2* histogram, TH2* rCumulHist, double scale){

  if ( histogram->GetEntries() == 0 ){
    return;
  }
  
// bin = 0;       underflow bin
// bin = 1;       first bin with low-edge xlow INCLUDED
// bin = nbins;   last bin with upper-edge xup EXCLUDED
// bin = nbins+1; overflow bin

  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  int xBinOver        = nBinsX + 1;
  int yBinOver        = nBinsY + 1;

  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){     
    for ( int iBinY = nBinsY; iBinY > 0; --iBinY ){     
      double integral = histogram->Integral( iBinX, xBinOver, iBinY, yBinOver );
      rCumulHist->SetBinContent( iBinX, iBinY, integral );
    }
  }      

  rCumulHist->Scale( scale );

}

void fillUniform2D( TH2* histogram, double value){

  if ( histogram->GetEntries() == 0 ){    return;  }

  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  int xBinOver        = nBinsX + 1;
  int yBinOver        = nBinsY + 1;

  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){     
    for ( int iBinY = nBinsY; iBinY > 0; --iBinY ){     
      double integral = histogram->Integral( iBinX, xBinOver, iBinY, yBinOver );
      histogram->SetBinContent( iBinX, iBinY, value );
    }
  }      

}


// Return the bins for seed1 vs seed2 correlations for which the rate is below a given threshold
std::vector< std::pair< int, int> > getEfficiencyBins( TH2* histogram, double maxRate ){

  std::vector< std::pair< int, int> > selectedBins;

  // Find bin with minimum cut for acceptable rate
  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  
  // Scan each y column for constant x to find minimum threshold that is below cut
  for ( int iBinX = 0; iBinX <= nBinsX; ++iBinX ){     

    for ( int iBinY = 0; iBinY <= nBinsY; ++iBinY ){
      float currRate = histogram->GetBinContent( iBinX, iBinY );
    

      std::cout << currRate << "\t" << iBinX << "\t" << iBinY << "\n";
//       if ( (iBinY == 0) && (currRate == 0) ){ break; }

//       if ( currRate < maxRate ){
// 	//	iBinY -= 1; // previous bin was last below acceptable rate

// 	std::cout << "Selected: " << currRate << "\t" << iBinX << "\t" << iBinY << "\n\n";                                        
// 	selectedBins.push_back( std::make_pair(iBinX, iBinY) );                                                                   
// 	break;                                                 
//       }

            if ( fabs( currRate - maxRate ) < maxRate/2.){
	      //	std::cout << currRate << "\t" << iBinX << "\t" << iBinY << "\n";
	      if ( currRate < maxRate ){
	  	  std::cout << "Selected: " << currRate << "\t" << iBinX << "\t" << iBinY << "\n\n";
		  selectedBins.push_back( std::make_pair(iBinX, iBinY) );
		  break;
	      }
	    }
    }
  }
  
  return selectedBins;
}

bool passesAlphaTSelection( std::vector<float> *jetPT, int &jetsAbove50 ){

    bool passesAlphaT(true);
    for (unsigned int iJet = 0; iJet < jetPT->size(); ++iJet ){
      
      float jetPt = (*jetPT)[ iJet ];

      // Leading two jet cut
      if ( (iJet+1 <= 2) && ( jetPt < 100. ) ){
	passesAlphaT = false;
	break;
      }

      if ( jetPt > 50. ){ jetsAbove50++; }
      else{ break; }

    }
    // Require at least two jets
    if ( (passesAlphaT) && ( jetsAbove50 < 2 ) ){
      passesAlphaT = false;
    }

    return passesAlphaT;
}



// Check whether event passes low-HT alpha analysis cuts
inline bool passesAlphaTFullSelection( std::vector<float> *jetPT, float HT, float AlphaT, int &jetsAbove50 ){
  // Offline HT alphaT cuts
  bool ht200to275(false), ht275to325(false), ht325to375(false);
  if ( ((HT >= 200) && (HT < 275)) && (AlphaT > 0.65) ){ ht200to275 = true; }
  if ( ((HT >= 275) && (HT < 325)) && (AlphaT > 0.60) ){ ht275to325 = true; }
  if ( ((HT >= 325) && (HT < 375)) && (AlphaT > 0.55) ){ ht325to375 = true; }



  for (unsigned int iJet = 0; iJet < jetPT->size(); ++iJet ){
      
      float jetPt = (*jetPT)[ iJet ];

      // Leading two jet cut
      if ( (iJet+1 <= 2) && ( jetPt < 100. ) ){
	return false;
      }
      if ( jetPt > 50. ){ jetsAbove50++; }
      else{ break; }

  }
  
  // Count number of analysis jets
  if ( jetsAbove50 > MAX_JETS ){
    jetsAbove50 = MAX_JETS;
  }


  if ( (*jetPT)[1] >= 100  ){
    if ( ht200to275 || ht275to325 || ht325to375 ){
      return true;
    }
  }
  return false;
}










// Returns all HT alphaT pairs
std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold, uint maxJets, float dynamicJetThreshold){

  float currentHT(0);
  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;
  //  std::vector<float> jetPhiNew;
  std::vector<std::pair<float,float> > alphaTHTPairs;

  uint jetLimit = std::min( (uint)jetPT->size(), maxJets );

  for(uint iJet = 0;iJet < jetLimit; ++iJet ){

    if( (*jetPT)[iJet] < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( (*jetPT)[iJet] );
    jetPyNew .push_back( (*jetPx)[iJet] );
    jetPxNew .push_back( (*jetPy)[iJet] );
    //    jetPhiNew.push_back( (*jetPhi)[iJet] );
    currentHT += (*jetPT)[iJet];

    //    if ( iJet > 0 ){
      float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 );
      alphaTHTPairs.push_back( std::make_pair( aT, currentHT ) );
      //    }
  }

  return alphaTHTPairs;
}





std::pair<float,float> calculateDynamicAlphaTHT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, 
						float jetThreshold, uint maxJets, 
						float dynamicJetThreshold, float dynamicAlphaTThreshold, float dynamicHTThreshold){
  
  float maxAlphaT(-1);
  float currentHT(0);
  
  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;
  //  std::vector<float> jetPhiNew;


  uint jetLimit = std::min( (uint)jetPT->size(), maxJets );

  for(uint iJet = 0;iJet < jetLimit; ++iJet ){

    //    if( std::abs(ijet->eta()) > etaJet_.at(0) ) continue;
    if( (*jetPT)[iJet] < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( (*jetPT)[iJet] );
    jetPyNew .push_back( (*jetPx)[iJet] );
    jetPxNew .push_back( (*jetPy)[iJet] );
    //    jetPhiNew.push_back( (*jetPhi)[iJet] );
    currentHT += (*jetPT)[iJet];

    //    if ( iJet > 0 ){
      float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 );
      if ( aT > maxAlphaT ){ maxAlphaT = aT; }
      if ( (maxAlphaT > dynamicAlphaTThreshold) && (currentHT > dynamicHTThreshold) ){
	return std::make_pair( maxAlphaT, currentHT );
      }
      //    }
  }

  return std::make_pair( -1, 0 );
}


// Ported from https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/HLTrigger/JetMET/src/HLTAlphaTFilter.cc
float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, 
			     float jetThreshold, uint maxJets, 
			     float dynamicJetThreshold, float dynamicAlphaTThreshold){

  float maxAlphaT(-1);

  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;
  //  std::vector<float> jetPhiNew;

  uint jetLimit = std::min( (uint)jetPT->size(), maxJets );

  for(uint iJet = 0;iJet < jetLimit; ++iJet ){

    //    if( std::abs(ijet->eta()) > etaJet_.at(0) ) continue;
    if( (*jetPT)[iJet] < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( (*jetPT)[iJet] );
    jetPxNew .push_back( (*jetPx)[iJet] );
    jetPyNew .push_back( (*jetPy)[iJet] );
    //    jetPhiNew.push_back( (*jetPhi)[iJet] );
    
    //    if ( iJet > 0 ){
      float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 );
      if ( aT > maxAlphaT ){ maxAlphaT = aT;}
      if ( maxAlphaT > dynamicAlphaTThreshold ){ break; }
    }
  //  }

  return maxAlphaT;
}



float calculateAlphaT( std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold ){

  // check the size of the input collection
  if (jetPT->size() <= 1){
    // empty jet collection, return AlphaT = 0
    return 0.;
  }
//   if (jetPT->size() == 1){ // Monojet!
//     return 10.;
//   }

//   if (jetPT->size() > (unsigned int) std::numeric_limits<unsigned int>::digits){
//     // too many jets, return AlphaT = a very large number
//     return std::numeric_limits<double>::max();
//   }
  
  // Momentum sums in transverse plane
  float sum_et(0), sum_px(0), sum_py(0);

  // Jet collection restricted to jet threshold
  std::vector<float> jetPTNew;
  
  for (unsigned int iJet = 0; iJet < jetPT->size(); ++iJet ){

    if ( (*jetPT)[iJet] < jetThreshold ){ break; }
    jetPTNew.push_back( (*jetPT)[iJet] );

    sum_et += (*jetPT)[iJet];    
    sum_px += (*jetPx)[iJet];
    sum_py += (*jetPy)[iJet];
    
  }

  // check the size of the new input collection
  if (jetPTNew.size() <= 1){
    // empty jet collection, return AlphaT = 0
    return 0.;
  }
//   if (jetPTNew.size() == 1){ // Monojet!
//     return 10.;
//   }

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
    return 10.;
  }

  // Alpha_T
  return (0.5 * (sum_et - min_delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) ));  
}




//  LocalWords:  Str
