#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include <sys/stat.h>
#include <map>
#include <iostream>
#include <iomanip> 
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
#include "TEventList.h"
#include "TROOT.h"
#include <vector>
#include <math.h>

#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/interface/dynamic.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/typeDefs.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/TriggerTuple/DynAlphaT.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/TriggerTuple/branches.cpp"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/TriggerTuple/integration.cpp"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/TriggerTuple/trigger.cpp"
//#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/TriggerTuple/samples/samples_v2.cpp"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_2/src/AlphaTHLT/Macros/TriggerTuple/samples/samples_v3.cpp"

//const int MAX_EVENTS = 10000;
const int MAX_EVENTS = -1;

// FOR CALOPREFILTER
TString PF2Jet =  "dijetAvggt90";

//#ifdef PRINT

// Remove vetoes from signal selection
//#define NOVETO
#define NOLEPTONVETO

//#define ROOTFILE // Required for ntuple production
//#define TURNON
// // #define PREFILTER
// #define HIST2DTURN
// //#define HIST1DTURN

// FOR RATE
//#define MAKE_TRIGGER_TUPLE
//#define HIST2DTURNSAVE
//#define HIST1DTURNSAVE

//#define HTTRIGGERS



#define MAKE_TRIGGER_LOG

//#define DYNALPHAT
//#define DYNALPHAT_NEW

#define TRIGGERBITS

// ****************************************
// output
// ****************************************
//TString outdir = "output/test/";
//TString outdir = "output/23_04_15/";
//TString outdir = "output/23_04_15_Ntuples_Bosh/";
//TString outdir = "output/23_04_15_Job1/";
//TString outdir = "output/23_04_15_Job2/";
//TString outdir = "output/23_04_15_50nsWPEff/";
TString outdir = "output/02_06_15_Test/";

TString fileSuffix;



// Containers 
std::map< TString, TH1*> hist;

std::map< TString, TH1*> hist1D;
std::map< TString, TH2*> hist2D;
std::map< TString, TH1*> hist1DRate;
std::map< TString, TH2*> hist2DRate;
std::map< TString, TEfficiency*> hist1DEff;
std::map< TString, TEfficiency*> hist2DEff;
std::map< TString, TEfficiency*> hist2DTurn;
std::map< TString, TEfficiency*> hist2DPrefTurn;

// Dynamic triggers 
//std::map< TString, dynamicRate > dynHistRate;
//std::map< TString, dynamicRate > dynHistEff;


TString offRecoType ;
TString offJet2PT   ;
TString offHT       ;
TString offAlphaT   ;
TString offMHT      ;
TString offMET      ;
TString offForJetPT; 
TString L1HTT     ;
TString L1ETM     ;
TString hltPFJetPT;
TString hltCaloHT;
TString hltCaloAlphaT;
TString hltCaloAlphaTPrime;

std::vector<float>   jet2PTCuts;
std::vector<TString> recoTypes;

std::map<int, std::pair< TString, TString> > jetBinStrLabels;
std::map<int, std::pair< TString, TString> > HTBinStrLabels;


// HT triggers
  std::vector<int> L1Trigger;
  std::vector<int> PFHTTrigger;
  std::vector<int> CaloHTTrigger;



bool passL1HTTorETM(false);   // Level-1 seed
bool passHLTJet2PT40(false);  // Jet2 > 40 GeV
bool passHLTJet2PT(false);    // Jet2 > threshold
bool passHLTDijetAvg(false);  // Dijet avg > threshold
bool passHLTHT1(false);       // HT WP1
bool passHLTHT2(false);       // HT WP2
bool passHLTHT3(false);       // HT WP3
bool passHLTHT4(false);       // HT WP4
bool passHLTHT5(false);       // HT WP5
bool passHLTAlphaT1(false);   // AlphaT WP1
bool passHLTAlphaT2(false);   // AlphaT WP2
bool passHLTAlphaT3(false);   // AlphaT WP3
bool passHLTAlphaT4(false);   // AlphaT WP4
bool passHLTAlphaT5(false);   // AlphaT WP5

bool passCaloDijet60(false);
bool passCaloDijet70(false);

bool passCaloPrefilter1(false);
bool passCaloPrefilter2(false);
bool passCaloPrefilter3(false);
bool passCaloPrefilter4(false);
bool passCaloPrefilter5(false);
bool passCaloBrokenPrefilter1(false);
bool passCaloBrokenPrefilter2(false);
bool passCaloBrokenPrefilter3(false);
bool passCaloBrokenPrefilter4(false);
bool passCaloBrokenPrefilter5(false);

bool passOffJet2PT40(false);  // Jet2 > 40 GeV
bool passOffJet2PT(false);    // Symmetric/asymmetric jet thresholds 
bool passOffHT200(false);     // HT   > 200 GeV
bool passOffAlphaT(false);    // Corresponding AlphaT to HT cut
bool passOffLepton(false);    // Lepton veto
bool passOffForJet(false);    // Forward jet veto
bool passOffMoM1p25(false);   // MHT/MET < 1.25
bool passOffVetoes(false);

bool passOffHT0(false);
bool passOffHT1(false);
bool passOffHT2(false);
bool passOffHT3(false);
bool passOffHT4(false);
bool passOffHT5(false);

std::map< int, bool > passOffHT;

// ****************************************************************************************************
// ****************************************************************************************************

  int minHT(175), maxHT(425), stepHT(25);
  float minAlphaT(0.50), maxAlphaT(0.60), stepAlphaT(0.005);
  int minJet2(40), maxJet2(110), stepJet2(10);


// Analysis jet thresholds 
float hltCaloJetThreshold = 40; //50; 
float hltPFJetThreshold   = 40; //50; 
float offJetThreshold     = 40; //50; 
float dynamicJetThreshold = 40;
// Second jet thresholds 
float hltPreCaloSecondJetThreshold = 70;
float hltPFSecondJetThreshold      = 90;
float hltCaloSecondJetThreshold    = 90;
const float offJet1PTThreshold     = 100;
const float offJet2PTThreshold     = 100;
float MoMVeto                      = 2.0;
float ForPTVeto                    = 70;
float MinHLTMET                    = 100.;

bool bx50(false);


inline float ATThresh(float offHT){
  if      ( offHT >= 500.){ return 0.52;}
  else if ( offHT >= 400.){ return 0.52;}
  else if ( offHT >= 350.){ return 0.53;}
  else if ( offHT >= 300.){ return 0.55;}
  else if ( offHT >= 250.){ return 0.60;}
  else if ( offHT >= 200.){ return 0.65;}
  return 99999;
}

inline float ATThresh(int offHTBin){
  if      ( offHTBin == 5){ return 0.52;}
  else if ( offHTBin == 4){ return 0.52;}
  else if ( offHTBin == 3){ return 0.53;}
  else if ( offHTBin == 2){ return 0.55;}
  else if ( offHTBin == 1){ return 0.60;}
  else if ( offHTBin == 0){ return 0.65;}
  return 99999;
}



std::map<TString, bool>  passTrigger;
std::map<TString, std::pair<float,float> > HLTHTAlphaTTrigger;

trigger::trigger     HTTrigger(TString WP, TString suffix = "" ){
    TString HTStr     = WP + "HT" + suffix;
    passTrigger[ HTStr] = false;
    return trigger::trigger( HTStr,      &passTrigger[ HTStr] ); 
}
trigger::trigger AlphaTTrigger(TString WP, TString suffix = "" ){
    TString AlphaTStr = WP + "AlphaT" + suffix;
    passTrigger[ AlphaTStr] = false;
    return trigger::trigger( AlphaTStr,  &passTrigger[ AlphaTStr] );
}
trigger::trigger DynHTAlphaTTrigger(TString WP, TString suffix = "", bool label=true ){
    TString AlphaTStr = WP;
    if (label){ AlphaTStr+= "DynHTAlphaT" + suffix; }
    passTrigger[ AlphaTStr] = false;
    return trigger::trigger( AlphaTStr,  &passTrigger[ AlphaTStr] );
}
trigger::trigger DynAlphaTTrigger(TString WP, TString suffix = "", bool label=true ){
    TString AlphaTStr = WP;
    if (label){ AlphaTStr+= "DynAlphaT" + suffix; }
    passTrigger[ AlphaTStr] = false;
    return trigger::trigger( AlphaTStr,  &passTrigger[ AlphaTStr] );
}




sampleCollection samples;
std::map<TString, std::vector<TString> > sampleStrs;
TChain *chain;
TFile* fOut;

std::map<TString, float>                 fBranch;
std::map<TString, UInt_t>                uBranch;
std::map<TString, Int_t>                 iBranch;
std::map<TString, std::vector<float>* > fvBranch;
std::map<TString, std::vector<std::pair<float,float> >* > pfvBranch;

std::vector<TString> selSampleStrs;


std::vector<int> Jet2PTVec;
std::vector<int> DijetAvgPTVec;
std::vector<int> METVec;


std::vector<TString> SelectedTriggerMenues;

std::vector<TString> turnOnTriggers;
std::vector<TString> turnOnPrefilterTriggers;
std::map< TString, trigger::triggerMenu > TriggerMenu;


std::map< TString, trigger::trigger > CaloHLTTriggers;
std::map< TString, trigger::trigger > CaloBrokenHLTTriggers;

// HLT paths                                                                                                                            
std::vector<TString>    hltPathNames;
std::map<TString, bool> passHLTPath;


TString selSampleStr;




struct HTAlphaTJet2PTTrigger{
  // Jet reconstruction/filter level 
  TString Type;
  
  float HT;
  float AlphaT;
  float Jet2PT;
  
  float Rate;
  float Eff;

  int Colour;

  HTAlphaTJet2PTTrigger():Type(""),HT(-1),AlphaT(-1),Jet2PT(-1),Rate(-1),Eff(-1),Colour(kRed){};
  HTAlphaTJet2PTTrigger( TString aType, float aHT, float aAlphaT, float aJet2PT, int aColour=kRed ):Type(aType),HT(aHT),AlphaT(aAlphaT),Jet2PT(aJet2PT),Rate(-1),Eff(-1),Colour(aColour){}

};








void initialise();
void beginJob();
void endJob();
void analyse(TString sampleStr );



// ********************************************************************************
//
// ********************************************************************************

void triggerTuple(){

  mkdir(outdir.Data(),S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


  // --------------------------------------------------------------------------------
  // selSampleStr = "PU40bx50_HCAL3_HPUV_QCD";
  //  selSampleStr = "PU40bx50_HPUV_QCD";  
  // selSampleStr  = "PU40bx25_HCAL3_HPUV_QCD";
  // selSampleStr = "PU40bx25_HPUV_QCD";

  //selSampleStr = "SM";
  //selSampleStr = "PU20bx25_HCAL3_Signal";
  //selSampleStr = "PU20bx25_Signal";
  //selSampleStr = "test";
  
  //selSampleStr = "PU40bx25_HPUV_QCD";
  selSampleStr = "PU40bx25_HCAL3_HPUV_QCD";
  //selSampleStr = "PU40bx50_HCAL3_HPUV_QCD";

  //selSampleStr = "PU20bx25_HPUV_QCD";

  // --------------------------------------------------------------------------------


  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.60 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 250., 0.57 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 300., 0.55 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 350., 0.53 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 400., 0.52 );

  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 175., 0.60 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 225., 0.57 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 275., 0.55 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 325., 0.53 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 375., 0.52 );

  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 175., 0.57 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 225., 0.55 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 275., 0.53 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 325., 0.52 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 375., 0.51 );

  // // 'Job 1'
  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.60 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 250., 0.56 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 300., 0.53 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 350., 0.52 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 400., 0.51 );

  // // 'Job 2'
  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.62 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 250., 0.57 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 300., 0.54 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 350., 0.52 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 400., 0.51 );


  HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.60 );
  HLTHTAlphaTTrigger["PF2"] = std::make_pair( 250., 0.55 );
  HLTHTAlphaTTrigger["PF3"] = std::make_pair( 300., 0.53 );
  HLTHTAlphaTTrigger["PF4"] = std::make_pair( 350., 0.52 );
  HLTHTAlphaTTrigger["PF5"] = std::make_pair( 400., 0.51 );

  // 2012 threshold
  HLTHTAlphaTTrigger["2012PF1"] = std::make_pair( 200., 0.63 ); 
  //  HLTHTAlphaTTrigger["2012PF1"] = std::make_pair( 200., 0.60 );
  HLTHTAlphaTTrigger["2012PF2"] = std::make_pair( 250., 0.55 );
  HLTHTAlphaTTrigger["2012PF3"] = std::make_pair( 300., 0.53 );
  HLTHTAlphaTTrigger["2012PF4"] = std::make_pair( 350., 0.52 );
  HLTHTAlphaTTrigger["2012PF5"] = std::make_pair( 400., 0.51 );

  // Raised thresholds
  HLTHTAlphaTTrigger["RTPF1"] = std::make_pair( 200., 0.63 );
  HLTHTAlphaTTrigger["RTPF2"] = std::make_pair( 250., 0.58 );
  HLTHTAlphaTTrigger["RTPF3"] = std::make_pair( 300., 0.54 );
  HLTHTAlphaTTrigger["RTPF4"] = std::make_pair( 350., 0.53 );
  HLTHTAlphaTTrigger["RTPF5"] = std::make_pair( 400., 0.52 );


  // // 27_04_15
  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.63 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 250., 0.58 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 300., 0.54 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 350., 0.53 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 400., 0.52 );

  // HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.590 );
  // HLTHTAlphaTTrigger["PF2"] = std::make_pair( 200., 0.580 );
  // HLTHTAlphaTTrigger["PF3"] = std::make_pair( 200., 0.570 );
  // HLTHTAlphaTTrigger["PF4"] = std::make_pair( 200., 0.560 );
  // HLTHTAlphaTTrigger["PF5"] = std::make_pair( 200., 0.550 );




  // Load all samples
  loadSamples( samples, sampleStrs );

  



#ifdef RUN_ON_BATCH
  selSampleStrs = sampleStrs["BATCH"];
  outdir = "BATCH";
#else
  selSampleStrs = sampleStrs[selSampleStr];
#endif
  initialise();

  // ****************************************************************************************************
  // *                                          Define triggers                                         *
  // ****************************************************************************************************

  trigger::triggerPath HLTNoL1J2BaseSelection;
  HLTNoL1J2BaseSelection.addTrigger( trigger::trigger("HLTJet2PT40", &passHLTJet2PT40) );

  trigger::triggerPath HLTL1Selection;
  HLTL1Selection.addTrigger( trigger::trigger("L1HTTorETM",  &passL1HTTorETM)  );
  
  trigger::triggerPath HLTNoJ2BaseSelection;
  HLTNoJ2BaseSelection.addTrigger( trigger::trigger("L1HTTorETM",  &passL1HTTorETM)  );
  HLTNoJ2BaseSelection.addTrigger( trigger::trigger("HLTJet2PT40", &passHLTJet2PT40) );


  trigger::triggerPath HLTBaseSelection;
  HLTBaseSelection.addTrigger( trigger::trigger("L1HTTorETM",  &passL1HTTorETM)  );
  HLTBaseSelection.addTrigger( trigger::trigger("HLTJet2PT40", &passHLTJet2PT40) );
  HLTBaseSelection.addTrigger( trigger::trigger("HLTJet2PT",   &passHLTJet2PT)   );

  trigger::triggerPath HLTBaseDijetAvgSelection;
  HLTBaseDijetAvgSelection.addTrigger( trigger::trigger("L1HTTorETM",    &passL1HTTorETM)  );
  HLTBaseDijetAvgSelection.addTrigger( trigger::trigger("HLTJet2PT40",   &passHLTJet2PT40) );
  HLTBaseDijetAvgSelection.addTrigger( trigger::trigger("HLTDijetAvgPT", &passHLTDijetAvg) );


  // ------------------------------
  // Prefilters
  // ------------------------------
  CaloHLTTriggers["PF1"] =        trigger::trigger("CaloPrefilter1",        &passCaloPrefilter1 );
  CaloHLTTriggers["PF2"] =        trigger::trigger("CaloPrefilter2",        &passCaloPrefilter2 );
  CaloHLTTriggers["PF3"] =        trigger::trigger("CaloPrefilter3",        &passCaloPrefilter3 );
  CaloHLTTriggers["PF4"] =        trigger::trigger("CaloPrefilter4",        &passCaloPrefilter4 );
  CaloHLTTriggers["PF5"] =        trigger::trigger("CaloPrefilter5",        &passCaloPrefilter5 );
  CaloBrokenHLTTriggers["PF1"] =  trigger::trigger("CaloBrokenPrefilter1",  &passCaloBrokenPrefilter1 );
  CaloBrokenHLTTriggers["PF2"] =  trigger::trigger("CaloBrokenPrefilter2",  &passCaloBrokenPrefilter2 );
  CaloBrokenHLTTriggers["PF3"] =  trigger::trigger("CaloBrokenPrefilter3",  &passCaloBrokenPrefilter3 );
  CaloBrokenHLTTriggers["PF4"] =  trigger::trigger("CaloBrokenPrefilter4",  &passCaloBrokenPrefilter4 );
  CaloBrokenHLTTriggers["PF5"] =  trigger::trigger("CaloBrokenPrefilter5",  &passCaloBrokenPrefilter5 );



  // ************************************************************************************************************************
  // *                                                       Turn ons                                                       
  // ************************************************************************************************************************


  // *****************************************************************************
  // *                              Trigger logging                              *
  // *****************************************************************************


 // SelectedTriggerMenues.push_back("PFNewTrigger");
 // SelectedTriggerMenues.push_back("PFNewTriggerNJet");

  // SelectedTriggerMenues.push_back( "PFJ2StdTrigger");
  // SelectedTriggerMenues.push_back( "PFJ2DynTrigger");
  // SelectedTriggerMenues.push_back( "PFAveStdTrigger");
  // SelectedTriggerMenues.push_back( "PFAveDynTrigger");
  // SelectedTriggerMenues.push_back( "PF2012J2StdTrigger");
  // SelectedTriggerMenues.push_back( "PF2012J2DynTrigger");
  // SelectedTriggerMenues.push_back( "PF2012AveStdTrigger");
  // SelectedTriggerMenues.push_back( "PF2012AveDynTrigger");

  // SelectedTriggerMenues.push_back( "PFJ2StdTriggerNJet");
  // SelectedTriggerMenues.push_back( "PFJ2DynTriggerNJet");
  // SelectedTriggerMenues.push_back( "PFAveStdTriggerNJet");
  // SelectedTriggerMenues.push_back( "PFAveDynTriggerNJet");
  // SelectedTriggerMenues.push_back( "PF2012J2StdTriggerNJet");
  // SelectedTriggerMenues.push_back( "PF2012J2DynTriggerNJet");
  // SelectedTriggerMenues.push_back( "PF2012AveStdTriggerNJet");
  // SelectedTriggerMenues.push_back( "PF2012AveDynTriggerNJet");



 // SelectedTriggerMenues.push_back("PFBrokenHLTTrigger");
 // SelectedTriggerMenues.push_back("PFHLTTriggers");
  SelectedTriggerMenues.push_back("HLTTriggerBits");

  //  SelectedTriggerMenues.push_back("CumulPFNewTriggerNJet");

  // ********************************************************************************
  // Prefilters
  // ********************************************************************************
  //  turnOnPrefilterTriggers.push_back("PFJ280AlphaTTriggers");

  turnOnPrefilterTriggers.push_back("PFJ280DynAlphaTTriggers");
  turnOnPrefilterTriggers.push_back("PFDj80DynAlphaTTriggers");

  // turnOnPrefilterTriggers.push_back("PFDj80AlphaTTriggers");
  // turnOnPrefilterTriggers.push_back("PFDj80DynAlphaTTriggers");


  // ********************************************************************************
  // Final filters
  // ********************************************************************************

  // Dijet
  turnOnTriggers.push_back("PFJ2AlphaTTriggers");
  turnOnTriggers.push_back("PFJ2DynAlphaTTriggers");
  // Dijet average
  turnOnTriggers.push_back("PFDjAlphaTTriggers");
  turnOnTriggers.push_back("PFDjDynAlphaTTriggers");
  
  // turnOnTriggers.push_back("L1Triggers");
  // turnOnTriggers.push_back("PFJet2Triggers");
  // turnOnTriggers.push_back("PFHTTriggers");


  // ************************************************************************************************************************
  // ************************************************************************************************************************

  // ------------------------------
  // PF triggers
  // ------------------------------
  TriggerMenu["L1Triggers"].setName("L1Triggers");
  TriggerMenu["L1Triggers"].addPath( "L1", HLTL1Selection);



  Jet2PTVec.push_back(40);
  Jet2PTVec.push_back(50);
  Jet2PTVec.push_back(60);
  Jet2PTVec.push_back(70);
  Jet2PTVec.push_back(80);
  Jet2PTVec.push_back(90);
  Jet2PTVec.push_back(100.);

  DijetAvgPTVec.push_back(40);
  DijetAvgPTVec.push_back(50);
  DijetAvgPTVec.push_back(60);
  DijetAvgPTVec.push_back(70);
  DijetAvgPTVec.push_back(80);
  DijetAvgPTVec.push_back(90);
  DijetAvgPTVec.push_back(100);
  //  DijetAvgPTVec.push_back(120);


  

  // *****************************************************************************



  for (int i = 0; i < 5; ++i){


	TString Jet2PTStr     =  "Jet2gt90";  // 
	TString DijetAvgPTStr =  "dijetAvggt90";


	TString WP  = TString("PF") + TString(Form("%d", (i + 1)));
	TString cWP = TString("p") + TString(Form("%d", (i + 1)));
	TString HTStr          = WP + "HT";
	TString AlphaTStr      = WP + "AlphaT";
	TString DynAlphaTHTStr = WP + "DynHTAlphaT";
	
	TString RTHTStr          = "RT" + HTStr;
        TString RTAlphaTStr      = "RT" + AlphaTStr;
	TString RTDynAlphaTHTStr = "RT" + DynAlphaTHTStr;

	TString WP2012  = "2012" + WP;
	TString cWP2012 = cWP + "2012";
	TString HT2012Str          = WP2012 + "HT";
	TString AlphaT2012Str      = WP2012 + "AlphaT";
	TString DynAlphaTHT2012Str = WP2012 + "DynHTAlphaT";
	


	//	for (int dAlphaT = 0; dAlphaT < 5; ++dAlphaT){
	for (int dAlphaT = 0; dAlphaT < 5; ++dAlphaT){
	  if ((dAlphaT != 0)){ break; }

	  float aTStep = 0.01;
	  TString aTOffset = "";
	  if (dAlphaT != 0){
	    if (i > 1){ aTStep = 0.005; }
	    aTOffset = TString(Form("%1.3f", dAlphaT*aTStep ));
	    aTOffset.ReplaceAll(".","p");
	  }
	  
	  TriggerMenu["PFNewTrigger"].setName("PFNewTrigger");
	  TriggerMenu["PFNewTrigger"].addPath(    WP+Jet2PTStr+aTOffset, HLTNoJ2BaseSelection);
	  TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(cWP + "Jet2",  &passTrigger[ cWP + "Jet2" ]) );
	  TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(cWP + "HT",    &passTrigger[ cWP + "HT" ]) );
	  TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(cWP + "AT",    &passTrigger[ cWP + "AT" ]) );
	  TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(Jet2PTStr, &passTrigger[ Jet2PTStr ]) );
	  // TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr, DynAlphaTTrigger( WP ) );
	  // TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr, DynHTAlphaTTrigger( WP ) );
	  TriggerMenu["PFNewTrigger"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(DynAlphaTHTStr + aTOffset , &passTrigger[ DynAlphaTHTStr + aTOffset ] ));
	  TriggerMenu["PFNewTrigger"].setSignal(  WP+Jet2PTStr+aTOffset, passOffHT[i] );



	  // ********************************************************************************
	  // 2012
	  // ********************************************************************************
	  TriggerMenu["PF2012J2StdTrigger"].setName("PF2012J2StdTrigger");
          TriggerMenu["PF2012J2StdTrigger"].addPath(    WP2012+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012J2StdTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "Jet2",  &passTrigger[ cWP2012 + "J2Jet2" ]) );
          TriggerMenu["PF2012J2StdTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "HT",    &passTrigger[ cWP2012 + "J2HT" ]) );
          TriggerMenu["PF2012J2StdTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "AT",    &passTrigger[ cWP2012 + "J2AT" ]) );
          TriggerMenu["PF2012J2StdTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(Jet2PTStr,         &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PF2012J2StdTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(HT2012Str,         &passTrigger[ HT2012Str ] ));
	  TriggerMenu["PF2012J2StdTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(AlphaT2012Str,     &passTrigger[ AlphaT2012Str ] ));
          TriggerMenu["PF2012J2StdTrigger"].setSignal(  WP2012+Jet2PTStr, passOffHT[i] );

	  TriggerMenu["PF2012J2DynTrigger"].setName("PF2012J2DynTrigger");
          TriggerMenu["PF2012J2DynTrigger"].addPath(    WP2012+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "Jet2",   &passTrigger[ cWP2012 + "J2Jet2" ]) );
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "HT",     &passTrigger[ cWP2012 + "J2HT" ]) );
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "AT",     &passTrigger[ cWP2012 + "J2AT" ]) );
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(Jet2PTStr,          &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(DynAlphaTHT2012Str, &passTrigger[ DynAlphaTHT2012Str ] ));
          TriggerMenu["PF2012J2DynTrigger"].setSignal(  WP2012+Jet2PTStr, passOffHT[i] );


	  TriggerMenu["PF2012J2DynTrigger"].setName("PF2012J2DynTrigger");
          TriggerMenu["PF2012J2DynTrigger"].addPath(    WP2012+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "Jet2",   &passTrigger[ cWP2012 + "J2Jet2" ]) );
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "HT",     &passTrigger[ cWP2012 + "J2HT" ]) );
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "AT",     &passTrigger[ cWP2012 + "J2AT" ]) );
          TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(Jet2PTStr,          &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PF2012J2DynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(DynAlphaTHT2012Str, &passTrigger[ DynAlphaTHT2012Str ] ));
          TriggerMenu["PF2012J2DynTrigger"].setSignal(  WP2012+Jet2PTStr, passOffHT[i] );


	  TriggerMenu["PF2012AveDynTrigger"].setName("PF2012AveDynTrigger");
          TriggerMenu["PF2012AveDynTrigger"].addPath(    WP2012+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012AveDynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "Jet2",   &passTrigger[ cWP2012 + "Jet2" ]) );
          TriggerMenu["PF2012AveDynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "HT",     &passTrigger[ cWP2012 + "HT" ]) );
          TriggerMenu["PF2012AveDynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(cWP2012 + "AT",     &passTrigger[ cWP2012 + "AT" ]) );
          TriggerMenu["PF2012AveDynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(DijetAvgPTStr,      &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PF2012AveDynTrigger"].addTrigger( WP2012+Jet2PTStr, trigger::trigger(DynAlphaTHT2012Str, &passTrigger[ DynAlphaTHT2012Str ] ));
          TriggerMenu["PF2012AveDynTrigger"].setSignal(  WP2012+Jet2PTStr, passOffHT[i] );

	  // ********************************************************************************
	  // New thresholds
	  // ********************************************************************************
	  TriggerMenu["PFJ2StdTrigger"].setName("PFJ2StdTrigger");
          TriggerMenu["PFJ2StdTrigger"].addPath(    WP+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFJ2StdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "Jet2",  &passTrigger[ cWP + "J2Jet2" ]) );
          TriggerMenu["PFJ2StdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "HT",    &passTrigger[ cWP + "J2HT" ]) );
          TriggerMenu["PFJ2StdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "AT",    &passTrigger[ cWP + "J2AT" ]) );
          TriggerMenu["PFJ2StdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(Jet2PTStr,     &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PFJ2StdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(HTStr,         &passTrigger[ RTHTStr ] ));
	  TriggerMenu["PFJ2StdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(AlphaTStr,     &passTrigger[ RTAlphaTStr ] ));
          TriggerMenu["PFJ2StdTrigger"].setSignal(  WP+Jet2PTStr, passOffHT[i] );

	  TriggerMenu["PFJ2DynTrigger"].setName("PFJ2DynTrigger");
          TriggerMenu["PFJ2DynTrigger"].addPath(    WP+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFJ2DynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "Jet2",   &passTrigger[ cWP + "J2Jet2" ]) );
          TriggerMenu["PFJ2DynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "HT",     &passTrigger[ cWP + "J2HT" ]) );
          TriggerMenu["PFJ2DynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "AT",     &passTrigger[ cWP + "J2AT" ]) );
          TriggerMenu["PFJ2DynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(Jet2PTStr,      &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PFJ2DynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(DynAlphaTHTStr, &passTrigger[ RTDynAlphaTHTStr ] ));
          TriggerMenu["PFJ2DynTrigger"].setSignal(  WP+Jet2PTStr, passOffHT[i] );

	  TriggerMenu["PFAveStdTrigger"].setName("PFAveStdTrigger");
          TriggerMenu["PFAveStdTrigger"].addPath(    WP+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFAveStdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "Jet2",  &passTrigger[ cWP + "Jet2" ]) );
          TriggerMenu["PFAveStdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "HT",    &passTrigger[ cWP + "HT" ]) );
          TriggerMenu["PFAveStdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "AT",    &passTrigger[ cWP + "AT" ]) );
          TriggerMenu["PFAveStdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(DijetAvgPTStr, &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PFAveStdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(HTStr,         &passTrigger[ RTHTStr ] ));
	  TriggerMenu["PFAveStdTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(AlphaTStr,     &passTrigger[ RTAlphaTStr ] ));
          TriggerMenu["PFAveStdTrigger"].setSignal(  WP+Jet2PTStr, passOffHT[i] );

	  TriggerMenu["PFAveDynTrigger"].setName("PFAveDynTrigger");
          TriggerMenu["PFAveDynTrigger"].addPath(    WP+Jet2PTStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFAveDynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "Jet2",   &passTrigger[ cWP + "Jet2" ]) );
          TriggerMenu["PFAveDynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "HT",     &passTrigger[ cWP + "HT" ]) );
          TriggerMenu["PFAveDynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(cWP + "AT",     &passTrigger[ cWP + "AT" ]) );
          TriggerMenu["PFAveDynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(DijetAvgPTStr,  &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PFAveDynTrigger"].addTrigger( WP+Jet2PTStr, trigger::trigger(DynAlphaTHTStr, &passTrigger[ RTDynAlphaTHTStr ] ));
          TriggerMenu["PFAveDynTrigger"].setSignal(  WP+Jet2PTStr, passOffHT[i] );

	}	

	for (int iJet = -3; iJet <=3; ++iJet){
	  if (iJet == 0){ continue;}
	  // Offline selection
	  TString jetStr = jetBinStrLabels[ iJet ].first;
	  TString HTStr  = HTBinStrLabels[  i  ].first;

	  TString HTStrPlusOne = HTStr;
	  if (i < 5){
	    HTStrPlusOne = HTBinStrLabels[  i+1  ].first;
	  }
	  TString offStr        = HTStr + jetStr;
	  TString offStrPlusOne = HTStrPlusOne + jetStr;

	  // ********************************************************************************
	  // SEED OFFLINE BIN ONE ABOVE 
	  // ********************************************************************************
	  offStr = offStrPlusOne;
	  // ********************************************************************************

	  // // Cumulative trigger
	  // TriggerMenu["CumulPFNewTriggerNJet"].setName("CumulPFNewTriggerNJet");
	  // TriggerMenu["CumulPFNewTriggerNJet"].addPath(    WP+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
	  // TriggerMenu["CumulPFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger("Cumul" + WP, &passTrigger[ "Cumul" + WP ]) );
	  // TriggerMenu["CumulPFNewTriggerNJet"].setSignal(  WP+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PFNewTriggerNJet"].setName("PFNewTriggerNJet");
	  TriggerMenu["PFNewTriggerNJet"].addPath(    WP+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
	  TriggerMenu["PFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr,    trigger::trigger(cWP + "Jet2",  &passTrigger[ cWP + "Jet2" ]) );
	  TriggerMenu["PFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr,    trigger::trigger(cWP + "HT",    &passTrigger[ cWP + "HT" ]) );
	  TriggerMenu["PFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr,    trigger::trigger(cWP + "AT",    &passTrigger[ cWP + "AT" ]) );
	  TriggerMenu["PFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(Jet2PTStr+"_"+jetStr, &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, DynAlphaTTrigger( WP ) );
	  TriggerMenu["PFNewTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, DynHTAlphaTTrigger( WP ) );
	  TriggerMenu["PFNewTriggerNJet"].setSignal(  WP+Jet2PTStr+jetStr, passTrigger[offStr] );




	  // ********************************************************************************
	  // 2012
	  // ********************************************************************************
	  TriggerMenu["PF2012J2StdTriggerNJet"].setName("PF2012J2StdTrigger");
          TriggerMenu["PF2012J2StdTriggerNJet"].addPath(    WP2012+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012J2StdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "Jet2",  &passTrigger[ cWP2012 + "J2Jet2" ]) );
          TriggerMenu["PF2012J2StdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "HT",    &passTrigger[ cWP2012 + "J2HT" ]) );
          TriggerMenu["PF2012J2StdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "AT",    &passTrigger[ cWP2012 + "J2AT" ]) );
          TriggerMenu["PF2012J2StdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(Jet2PTStr,         &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PF2012J2StdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(HT2012Str,         &passTrigger[ HT2012Str ] ));
	  TriggerMenu["PF2012J2StdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(AlphaT2012Str,     &passTrigger[ AlphaT2012Str ] ));
          TriggerMenu["PF2012J2StdTriggerNJet"].setSignal(  WP2012+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PF2012J2DynTriggerNJet"].setName("PF2012J2DynTrigger");
          TriggerMenu["PF2012J2DynTriggerNJet"].addPath(    WP2012+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012J2DynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "Jet2",   &passTrigger[ cWP2012 + "J2Jet2" ]) );
          TriggerMenu["PF2012J2DynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "HT",     &passTrigger[ cWP2012 + "J2HT" ]) );
          TriggerMenu["PF2012J2DynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "AT",     &passTrigger[ cWP2012 + "J2AT" ]) );
          TriggerMenu["PF2012J2DynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(Jet2PTStr,          &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PF2012J2DynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(DynAlphaTHT2012Str, &passTrigger[ DynAlphaTHT2012Str ] ));
          TriggerMenu["PF2012J2DynTriggerNJet"].setSignal(  WP2012+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PF2012AveStdTriggerNJet"].setName("PF2012AveStdTrigger");
          TriggerMenu["PF2012AveStdTriggerNJet"].addPath(    WP2012+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012AveStdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "Jet2",  &passTrigger[ cWP2012 + "Jet2" ]) );
          TriggerMenu["PF2012AveStdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "HT",    &passTrigger[ cWP2012 + "HT" ]) );
          TriggerMenu["PF2012AveStdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "AT",    &passTrigger[ cWP2012 + "AT" ]) );
          TriggerMenu["PF2012AveStdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(DijetAvgPTStr,     &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PF2012AveStdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(HT2012Str,         &passTrigger[ HT2012Str ] ));
	  TriggerMenu["PF2012AveStdTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(AlphaT2012Str,     &passTrigger[ AlphaT2012Str ] ));
          TriggerMenu["PF2012AveStdTriggerNJet"].setSignal(  WP2012+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PF2012AveDynTriggerNJet"].setName("PF2012AveDynTrigger");
          TriggerMenu["PF2012AveDynTriggerNJet"].addPath(    WP2012+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PF2012AveDynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "Jet2",   &passTrigger[ cWP2012 + "Jet2" ]) );
          TriggerMenu["PF2012AveDynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "HT",     &passTrigger[ cWP2012 + "HT" ]) );
          TriggerMenu["PF2012AveDynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(cWP2012 + "AT",     &passTrigger[ cWP2012 + "AT" ]) );
          TriggerMenu["PF2012AveDynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(DijetAvgPTStr,      &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PF2012AveDynTriggerNJet"].addTrigger( WP2012+Jet2PTStr+jetStr, trigger::trigger(DynAlphaTHT2012Str, &passTrigger[ DynAlphaTHT2012Str ] ));
          TriggerMenu["PF2012AveDynTriggerNJet"].setSignal(  WP2012+Jet2PTStr+jetStr, passTrigger[offStr] );

	  // ********************************************************************************
	  // New thresholds
	  // ********************************************************************************
	  TriggerMenu["PFJ2StdTriggerNJet"].setName("PFJ2StdTrigger");
          TriggerMenu["PFJ2StdTriggerNJet"].addPath(    WP+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFJ2StdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "Jet2",  &passTrigger[ cWP + "J2Jet2" ]) );
          TriggerMenu["PFJ2StdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "HT",    &passTrigger[ cWP + "J2HT" ]) );
          TriggerMenu["PFJ2StdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "AT",    &passTrigger[ cWP + "J2AT" ]) );
          TriggerMenu["PFJ2StdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(Jet2PTStr,     &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PFJ2StdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(HTStr,         &passTrigger[ RTHTStr ] ));
	  TriggerMenu["PFJ2StdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(AlphaTStr,     &passTrigger[ RTAlphaTStr ] ));
          TriggerMenu["PFJ2StdTriggerNJet"].setSignal(  WP+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PFJ2DynTriggerNJet"].setName("PFJ2DynTrigger");
          TriggerMenu["PFJ2DynTriggerNJet"].addPath(    WP+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFJ2DynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "Jet2",   &passTrigger[ cWP + "J2Jet2" ]) );
          TriggerMenu["PFJ2DynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "HT",     &passTrigger[ cWP + "J2HT" ]) );
          TriggerMenu["PFJ2DynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "AT",     &passTrigger[ cWP + "J2AT" ]) );
          TriggerMenu["PFJ2DynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(Jet2PTStr,      &passTrigger[ Jet2PTStr ]) );
	  TriggerMenu["PFJ2DynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(DynAlphaTHTStr, &passTrigger[ RTDynAlphaTHTStr ] ));
          TriggerMenu["PFJ2DynTriggerNJet"].setSignal(  WP+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PFAveStdTriggerNJet"].setName("PFAveStdTrigger");
          TriggerMenu["PFAveStdTriggerNJet"].addPath(    WP+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFAveStdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "Jet2",  &passTrigger[ cWP + "Jet2" ]) );
          TriggerMenu["PFAveStdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "HT",    &passTrigger[ cWP + "HT" ]) );
          TriggerMenu["PFAveStdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "AT",    &passTrigger[ cWP + "AT" ]) );
          TriggerMenu["PFAveStdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(DijetAvgPTStr, &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PFAveStdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(HTStr,         &passTrigger[ RTHTStr ] ));
	  TriggerMenu["PFAveStdTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(AlphaTStr,     &passTrigger[ RTAlphaTStr ] ));
          TriggerMenu["PFAveStdTriggerNJet"].setSignal(  WP+Jet2PTStr+jetStr, passTrigger[offStr] );

	  TriggerMenu["PFAveDynTriggerNJet"].setName("PFAveDynTrigger");
          TriggerMenu["PFAveDynTriggerNJet"].addPath(    WP+Jet2PTStr+jetStr, HLTNoJ2BaseSelection);
          TriggerMenu["PFAveDynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "Jet2",   &passTrigger[ cWP + "Jet2" ]) );
          TriggerMenu["PFAveDynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "HT",     &passTrigger[ cWP + "HT" ]) );
          TriggerMenu["PFAveDynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(cWP + "AT",     &passTrigger[ cWP + "AT" ]) );
          TriggerMenu["PFAveDynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(DijetAvgPTStr,  &passTrigger[ DijetAvgPTStr ]) );
	  TriggerMenu["PFAveDynTriggerNJet"].addTrigger( WP+Jet2PTStr+jetStr, trigger::trigger(DynAlphaTHTStr, &passTrigger[ RTDynAlphaTHTStr ] ));
          TriggerMenu["PFAveDynTriggerNJet"].setSignal(  WP+Jet2PTStr+jetStr, passTrigger[offStr] );



	}
  }



  for (auto itr: Jet2PTVec){
    TString Jet2PTStr = TString("Jet2gt") + TString(Form("%d", (itr)));
    
    TriggerMenu["PFJet2Triggers"].setName("PFJet2Triggers");
    TriggerMenu["PFJet2Triggers"].addPath(Jet2PTStr, HLTNoJ2BaseSelection);
    TriggerMenu["PFJet2Triggers"].addTrigger( Jet2PTStr, trigger::trigger(Jet2PTStr, &passTrigger[ Jet2PTStr ]) );


    for (int i = 0; i < 5; ++i){
      TString WP = TString("PF") + TString(Form("%d", (i + 1)));

      for (int dAlphaT = 0; dAlphaT < 5; ++dAlphaT){
	float aTStep = 0.01;
	TString aTOffset = "";
	if (dAlphaT != 0){
	  if (i > 1){ aTStep = 0.005; }
	  aTOffset = TString(Form("%1.3f", dAlphaT*aTStep ));
	  aTOffset.ReplaceAll(".","p");
	}
	TString AlphaTStr = WP + "AlphaT";
	TString DynAlphaTHTStr = WP + "DynHTAlphaT";

	// TriggerMenu["PFNoL1J2AlphaTTriggers"].setName("PFNoL1J2AlphaTTriggers");
	// TriggerMenu["PFNoL1J2AlphaTTriggers"].addPath(    WP+Jet2PTStr, HLTNoL1J2BaseSelection);
	// TriggerMenu["PFNoL1J2AlphaTTriggers"].addTrigger( WP+Jet2PTStr, trigger::trigger(Jet2PTStr, &passTrigger[ Jet2PTStr ]) );
	// TriggerMenu["PFNoL1J2AlphaTTriggers"].addTrigger( WP+Jet2PTStr, AlphaTTrigger( WP ) );
	// TriggerMenu["PFNoL1J2AlphaTTriggers"].addTrigger( WP+Jet2PTStr, HTTrigger( WP ) );
	// TriggerMenu["PFNoL1J2AlphaTTriggers"].setSignal(  WP+Jet2PTStr, passOffHT[i] );
            
	// TriggerMenu["PFNoL1J2DynAlphaTTriggers"].setName("PFNoL1J2DynAlphaTTriggers");
	// TriggerMenu["PFNoL1J2DynAlphaTTriggers"].addPath(    WP+Jet2PTStr, HLTNoL1J2BaseSelection);
	// TriggerMenu["PFNoL1J2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr, trigger::trigger(Jet2PTStr, &passTrigger[ Jet2PTStr ]) );
	// TriggerMenu["PFNoL1J2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr, DynAlphaTTrigger( WP ) );
	// TriggerMenu["PFNoL1J2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr, DynHTAlphaTTrigger( WP ) );
	// TriggerMenu["PFNoL1J2DynAlphaTTriggers"].setSignal(  WP+Jet2PTStr, passOffHT[i] );
	
	TriggerMenu["PFJ2AlphaTTriggers"].setName("PFJ2AlphaTTriggers");
	TriggerMenu["PFJ2AlphaTTriggers"].addPath(    WP+Jet2PTStr+aTOffset, HLTNoJ2BaseSelection);
	TriggerMenu["PFJ2AlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(Jet2PTStr, &passTrigger[ Jet2PTStr ]) );
	TriggerMenu["PFJ2AlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(AlphaTStr + aTOffset, &passTrigger[AlphaTStr + aTOffset] ));//AlphaTTrigger( WP, aTOffset ) );
	TriggerMenu["PFJ2AlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, HTTrigger( WP ) );
	TriggerMenu["PFJ2AlphaTTriggers"].setSignal(  WP+Jet2PTStr+aTOffset, passOffHT[i] );
	
	TriggerMenu["PFJ2DynAlphaTTriggers"].setName("PFJ2DynAlphaTTriggers");
	TriggerMenu["PFJ2DynAlphaTTriggers"].addPath(    WP+Jet2PTStr+aTOffset, HLTNoJ2BaseSelection);
	TriggerMenu["PFJ2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(Jet2PTStr, &passTrigger[ Jet2PTStr ]) );
	// TriggerMenu["PFJ2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, DynAlphaTTrigger( WP, aTOffset ) );
	// TriggerMenu["PFJ2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, DynHTAlphaTTrigger( WP, aTOffset ) );
	TriggerMenu["PFJ2DynAlphaTTriggers"].addTrigger( WP+Jet2PTStr+aTOffset, trigger::trigger(DynAlphaTHTStr + aTOffset , &passTrigger[ DynAlphaTHTStr + aTOffset ] ));
	TriggerMenu["PFJ2DynAlphaTTriggers"].setSignal(  WP+Jet2PTStr+aTOffset, passOffHT[i] );
      }
    }
  }



  // Asymmetric dijet
  for (auto itr: DijetAvgPTVec){
    TString DijetAvgPTStr = TString("dijetAvggt") + TString(Form("%d", (itr)));
    
    for (int i = 0; i < 5; ++i){
      TString WP = TString("PF") + TString(Form("%d", (i + 1)));
      for (int dAlphaT = 0; dAlphaT < 5; ++dAlphaT){
	float aTStep = 0.01;
	TString aTOffset = "";
	if (dAlphaT != 0){
	  if (i > 1){ aTStep = 0.005; }
	  aTOffset = TString(Form("%1.3f", dAlphaT*aTStep ));
	  aTOffset.ReplaceAll(".","p");
	}

	TString AlphaTStr = WP + "AlphaT";
	TString DynAlphaTHTStr = WP + "DynHTAlphaT";


      // TriggerMenu["PFNoL1DjAlphaTTriggers"].setName("PFNoL1DjAlphaTTriggers");
      // TriggerMenu["PFNoL1DjAlphaTTriggers"].addPath(    WP+DijetAvgPTStr, HLTNoL1J2BaseSelection);
      // TriggerMenu["PFNoL1DjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr, trigger::trigger(DijetAvgPTStr, &passTrigger[ DijetAvgPTStr ]) );
      // TriggerMenu["PFNoL1DjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr, AlphaTTrigger( WP ) );
      // TriggerMenu["PFNoL1DjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr, HTTrigger( WP ) );
      // TriggerMenu["PFNoL1DjAlphaTTriggers"].setSignal(  WP+DijetAvgPTStr, passOffHT[i] );
            
      // TriggerMenu["PFNoL1DjDynAlphaTTriggers"].setName("PFNoL1DjDynAlphaTTriggers");
      // TriggerMenu["PFNoL1DjDynAlphaTTriggers"].addPath(    WP+DijetAvgPTStr, HLTNoL1J2BaseSelection);
      // TriggerMenu["PFNoL1DjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr, trigger::trigger(DijetAvgPTStr, &passTrigger[ DijetAvgPTStr ]) );
      // TriggerMenu["PFNoL1DjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr, DynAlphaTTrigger( WP ) );
      // TriggerMenu["PFNoL1DjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr, DynHTAlphaTTrigger( WP ) );
      // TriggerMenu["PFNoL1DjDynAlphaTTriggers"].setSignal(  WP+DijetAvgPTStr, passOffHT[i] );

      TriggerMenu["PFDjAlphaTTriggers"].setName("PFDjAlphaTTriggers");
      TriggerMenu["PFDjAlphaTTriggers"].addPath(    WP+DijetAvgPTStr+aTOffset, HLTNoJ2BaseSelection);
      TriggerMenu["PFDjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, trigger::trigger(DijetAvgPTStr, &passTrigger[ DijetAvgPTStr ]) );
      //      TriggerMenu["PFDjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, AlphaTTrigger( WP ) );
      TriggerMenu["PFDjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, trigger::trigger(AlphaTStr + aTOffset, &passTrigger[AlphaTStr + aTOffset] ));
      TriggerMenu["PFDjAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, HTTrigger( WP ) );
      TriggerMenu["PFDjAlphaTTriggers"].setSignal(  WP+DijetAvgPTStr+aTOffset, passOffHT[i] );
            
      TriggerMenu["PFDjDynAlphaTTriggers"].setName("PFDjDynAlphaTTriggers");
      TriggerMenu["PFDjDynAlphaTTriggers"].addPath(    WP+DijetAvgPTStr+aTOffset, HLTNoJ2BaseSelection);
      TriggerMenu["PFDjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, trigger::trigger(DijetAvgPTStr, &passTrigger[ DijetAvgPTStr ]) );
      // TriggerMenu["PFDjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, DynAlphaTTrigger( WP ) );
      // TriggerMenu["PFDjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, DynHTAlphaTTrigger( WP ) );
      TriggerMenu["PFDjDynAlphaTTriggers"].addTrigger( WP+DijetAvgPTStr+aTOffset, trigger::trigger(DynAlphaTHTStr + aTOffset , &passTrigger[ DynAlphaTHTStr + aTOffset ] ));
      TriggerMenu["PFDjDynAlphaTTriggers"].setSignal(  WP+DijetAvgPTStr+aTOffset, passOffHT[i] );

      }
    }

  }



  for (int i = 0; i < 5; ++i){
    TString WP = TString("PF") + TString(Form("%d", (i + 1)));

    // NoL1J2
    TriggerMenu["PFNoL1J2HTTriggers"].setName("PFNoL1J2HTTriggers");
    TriggerMenu["PFNoL1J2HTTriggers"].addPath(    WP, HLTNoL1J2BaseSelection);
    TriggerMenu["PFNoL1J2HTTriggers"].addTrigger( WP,        HTTrigger( WP ) );
    TriggerMenu["PFNoL1J2HTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoL1J2AlphaTTriggers"].setName("PFNoL1J2AlphaTTriggers");
    TriggerMenu["PFNoL1J2AlphaTTriggers"].addPath(    WP, HLTNoL1J2BaseSelection);
    TriggerMenu["PFNoL1J2AlphaTTriggers"].addTrigger( WP,         AlphaTTrigger( WP ) );
    TriggerMenu["PFNoL1J2AlphaTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoL1J2DynAlphaTTriggers"].setName("PFNoL1J2DynAlphaTTriggers");
    TriggerMenu["PFNoL1J2DynAlphaTriggers"].addPath(    WP, HLTNoL1J2BaseSelection);
    TriggerMenu["PFNoL1J2DynAlphaTTriggers"].addTrigger( WP,       DynAlphaTTrigger( WP ) );
    TriggerMenu["PFNoL1J2DynAlphaTTriggers"].setSignal(  WP, passOffHT[i] );

    // NoJ2
    TriggerMenu["PFHTTriggers"].setName("PFHTTriggers");
    TriggerMenu["PFHTTriggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFHTTriggers"].addTrigger( WP,        HTTrigger( WP ) );
    TriggerMenu["PFHTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFAlphaTTriggers"].setName("PFAlphaTTriggers");
    TriggerMenu["PFAlphaTTriggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFAlphaTTriggers"].addTrigger( WP,         AlphaTTrigger( WP ) );
    TriggerMenu["PFAlphaTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFDynAlphaTTriggers"].setName("PFDynAlphaTTriggers");
    TriggerMenu["PFDynAlphaTriggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFDynAlphaTTriggers"].addTrigger( WP,       DynAlphaTTrigger( WP ) );
    TriggerMenu["PFDynAlphaTTriggers"].setSignal(  WP, passOffHT[i] );




    TriggerMenu["PFBrokenHLTTrigger"].setName("PFBrokenHLTTrigger");
    TriggerMenu["PFBrokenHLTTrigger"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFBrokenHLTTrigger"].addTrigger( WP, CaloBrokenHLTTriggers[ WP ] );
    TriggerMenu["PFBrokenHLTTrigger"].addTrigger( WP,             HTTrigger( WP ) );
    TriggerMenu["PFBrokenHLTTrigger"].addTrigger( WP,         AlphaTTrigger( WP ) );
    TriggerMenu["PFBrokenHLTTrigger"].setSignal(  WP, passOffHT[i] );    


    TriggerMenu["PFNoPrefJ2HLTMin25Triggers"].setName("PFNoPrefJ2HLTMin25Triggers");
    TriggerMenu["PFNoPrefJ2HLTMin25Triggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFNoPrefJ2HLTMin25Triggers"].addTrigger( WP,        HTTrigger( WP, "Min25" ) );
    TriggerMenu["PFNoPrefJ2HLTMin25Triggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFNoPrefJ2HLTMin25Triggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoPrefJ2DynHLTMin25Triggers"].setName("PFNoPrefJ2DynHLTMin25Triggers");
    TriggerMenu["PFNoPrefJ2DynHLTMin25Triggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFNoPrefJ2DynHLTMin25Triggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP, "Min25") );
    TriggerMenu["PFNoPrefJ2DynHLTMin25Triggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoPrefJ2HLTTriggers"].setName("PFNoPrefJ2HLTTriggers");
    TriggerMenu["PFNoPrefJ2HLTTriggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFNoPrefJ2HLTTriggers"].addTrigger( WP,        HTTrigger( WP ) );
    TriggerMenu["PFNoPrefJ2HLTTriggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFNoPrefJ2HLTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoPrefJ2HLTTriggers"].setName("PFNoPrefJ2HLTTriggers");
    TriggerMenu["PFNoPrefJ2HLTTriggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFNoPrefJ2HLTTriggers"].addTrigger( WP,        HTTrigger( WP ) );
    TriggerMenu["PFNoPrefJ2HLTTriggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFNoPrefJ2HLTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoPrefJ2DynHLTTriggers"].setName("PFNoPrefJ2DynHLTTriggers");
    TriggerMenu["PFNoPrefJ2DynHLTTriggers"].addPath(    WP, HLTNoJ2BaseSelection);
    TriggerMenu["PFNoPrefJ2DynHLTTriggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP ) );
    TriggerMenu["PFNoPrefJ2DynHLTTriggers"].setSignal(  WP, passOffHT[i] );





    TriggerMenu["PFHLTTriggers"].setName("PFHLTTriggers");
    TriggerMenu["PFHLTTriggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFHLTTriggers"].addTrigger( WP,  CaloHLTTriggers[ WP ] );
    TriggerMenu["PFHLTTriggers"].addTrigger( WP,        HTTrigger( WP ) );
    TriggerMenu["PFHLTTriggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFHLTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFDynHLTTriggers"].setName("PFDynHLTTriggers");
    TriggerMenu["PFDynHLTTriggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFDynHLTTriggers"].addTrigger( WP,     CaloHLTTriggers[ WP ] );
    TriggerMenu["PFDynHLTTriggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP ) );
    TriggerMenu["PFDynHLTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFHLTMin25Triggers"].setName("PFHLTMin25Triggers");
    TriggerMenu["PFHLTMin25Triggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFHLTMin25Triggers"].addTrigger( WP,  CaloHLTTriggers[ WP ] );
    TriggerMenu["PFHLTMin25Triggers"].addTrigger( WP,        HTTrigger( WP, "Min25" ) );
    TriggerMenu["PFHLTMin25Triggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFHLTMin25Triggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoPrefHLTTriggers"].setName("PFNoPrefHLTTriggers");
    TriggerMenu["PFNoPrefHLTTriggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFNoPrefHLTTriggers"].addTrigger( WP,        HTTrigger( WP ) );
    TriggerMenu["PFNoPrefHLTTriggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFNoPrefHLTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFNoPrefHLTMin25Triggers"].setName("PFNoPrefHLTMin25Triggers");
    TriggerMenu["PFNoPrefHLTMin25Triggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFNoPrefHLTMin25Triggers"].addTrigger( WP,        HTTrigger( WP, "Min25" ) );
    TriggerMenu["PFNoPrefHLTMin25Triggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    TriggerMenu["PFNoPrefHLTMin25Triggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFDynHLTMin25Trigger"].setName("PFDynHLTMin25Trigger");
    TriggerMenu["PFDynHLTMin25Triggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFDynHLTMin25Triggers"].addTrigger( WP,     CaloHLTTriggers[ WP ] );
    TriggerMenu["PFDynHLTMin25Triggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP, "Min25") );
    TriggerMenu["PFDynHLTMin25Triggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFDynNoPrefHLTTriggers"].setName("PFDynNoPrefHLTTriggers");
    TriggerMenu["PFDynNoPrefHLTTriggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFDynNoPrefHLTTriggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP ) );
    TriggerMenu["PFDynNoPrefHLTTriggers"].setSignal(  WP, passOffHT[i] );

    TriggerMenu["PFDynNoPrefHLTMin25Triggers"].setName("PFDynNoPrefHLTMin25Triggers");
    TriggerMenu["PFDynNoPrefHLTMin25Triggers"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFDynNoPrefHLTMin25Triggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP, "Min25") );
    TriggerMenu["PFDynNoPrefHLTMin25Triggers"].setSignal(  WP, passOffHT[i] );

    // TriggerMenu["PFDijetAvgHLTTriggers"].setName("PFDijetAvgHLTTriggers");
    // TriggerMenu["PFDijetAvgHLTTriggers"].addPath(    WP, HLTBaseDijetAvgSelection);
    // TriggerMenu["PFDijetAvgHLTTriggers"].addTrigger( WP,  CaloHLTTriggers[ WP ] );
    // TriggerMenu["PFDijetAvgHLTTriggers"].addTrigger( WP,        HTTrigger( WP, "Min25" ) );
    // TriggerMenu["PFDijetAvgHLTTriggers"].addTrigger( WP,    AlphaTTrigger( WP ) );
    // TriggerMenu["PFDijetAvgHLTTriggers"].setSignal(  WP, passOffHT[i] );

    // TriggerMenu["PFDijetAvgDynHLTMin25Triggers"].setName("PFDijetAvgDynHLTMin25Triggers");
    // TriggerMenu["PFDijetAvgDynHLTMin25Triggers"].addPath(    WP, HLTBaseDijetAvgSelection);
    // TriggerMenu["PFDijetAvgDynHLTMin25Triggers"].addTrigger( WP,     CaloHLTTriggers[ WP ] );
    // TriggerMenu["PFDijetAvgDynHLTMin25Triggers"].addTrigger( WP,  DynHTAlphaTTrigger( WP, "Min25") );
    // TriggerMenu["PFDijetAvgDynHLTMin25Triggers"].setSignal(  WP, passOffHT[i] );



    // PREFILTER TRIGGERS
    TString Jet2PT80Str     = "Jet2gt80";
    TString DijetAvgPT80Str =  "dijetAvggt80";

      TriggerMenu["PFJ280AlphaTTriggers"].setName("PFJ280AlphaTTriggers");
      TriggerMenu["PFJ280AlphaTTriggers"].addPath(    WP+Jet2PT80Str, HLTNoJ2BaseSelection);
      TriggerMenu["PFJ280AlphaTTriggers"].addTrigger( WP+Jet2PT80Str, trigger::trigger(Jet2PT80Str, &passTrigger[ Jet2PT80Str ]) );
      TriggerMenu["PFJ280AlphaTTriggers"].addTrigger( WP+Jet2PT80Str, AlphaTTrigger( WP ) );
      TriggerMenu["PFJ280AlphaTTriggers"].setSignal(  WP+Jet2PT80Str, passOffHT[i] );
            
      TriggerMenu["PFJ280DynAlphaTTriggers"].setName("PFJ280DynAlphaTTriggers");
      TriggerMenu["PFJ280DynAlphaTTriggers"].addPath(    WP+Jet2PT80Str, HLTNoJ2BaseSelection);
      TriggerMenu["PFJ280DynAlphaTTriggers"].addTrigger( WP+Jet2PT80Str, trigger::trigger(Jet2PT80Str, &passTrigger[ Jet2PT80Str ]) );
      TriggerMenu["PFJ280DynAlphaTTriggers"].addTrigger( WP+Jet2PT80Str, DynHTAlphaTTrigger( WP ) );
      TriggerMenu["PFJ280DynAlphaTTriggers"].setSignal(  WP+Jet2PT80Str, passOffHT[i] );
 
      TriggerMenu["PFDj80AlphaTTriggers"].setName("PFDj80AlphaTTriggers");
      TriggerMenu["PFDj80AlphaTTriggers"].addPath(    WP+DijetAvgPT80Str, HLTNoJ2BaseSelection);
      TriggerMenu["PFDj80AlphaTTriggers"].addTrigger( WP+DijetAvgPT80Str, trigger::trigger(DijetAvgPT80Str, &passTrigger[ DijetAvgPT80Str ]) );
      TriggerMenu["PFDj80AlphaTTriggers"].addTrigger( WP+DijetAvgPT80Str, AlphaTTrigger( WP ) );
      TriggerMenu["PFDj80AlphaTTriggers"].setSignal(  WP+DijetAvgPT80Str, passOffHT[i] );
            

      TriggerMenu["PFDj80DynAlphaTTriggers"].setName("PFDj80DynAlphaTTriggers");
      TriggerMenu["PFDj80DynAlphaTTriggers"].addPath(    WP+DijetAvgPT80Str, HLTNoJ2BaseSelection);
      TriggerMenu["PFDj80DynAlphaTTriggers"].addTrigger( WP+DijetAvgPT80Str, trigger::trigger(DijetAvgPT80Str, &passTrigger[ DijetAvgPT80Str ]) );
      TriggerMenu["PFDj80DynAlphaTTriggers"].addTrigger( WP+DijetAvgPT80Str, DynHTAlphaTTrigger( WP ) );
      TriggerMenu["PFDj80DynAlphaTTriggers"].setSignal(  WP+DijetAvgPT80Str, passOffHT[i] );


  }



  // ------------------------------
  // Dynamic trigger scan
  // ------------------------------
#ifdef DYNALPHAT_NEW
  for (int Jet2 = minJet2; Jet2 < maxJet2; Jet2 += stepJet2){
    for (int HT = minHT; HT < maxHT; HT += stepHT){
      for (float AlphaT = minAlphaT; AlphaT < maxAlphaT; AlphaT += stepAlphaT){

	TString HTAlphaTStr = TString(Form("%d", Jet2)) + "_" + TString(Form("%d", HT)) + "_" + TString(Form("%1.2f", AlphaT));
	TString TriggerStr = "PF_DynNoPref_" + HTAlphaTStr; 
	// + TString(Form("%d", HT)) + "_AlphaT" + TString(Form("%1.2f", AlphaT));
	TriggerMenu[TriggerStr].setName(TriggerStr); 
	TriggerMenu[TriggerStr].addPath( TriggerStr, HLTNoJ2BaseSelection); 
	TriggerMenu[TriggerStr].addTrigger( TriggerStr, DynHTAlphaTTrigger( HTAlphaTStr, "", false ) ); 

      }
    }
  }
#endif

  // HLT trigger bits
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


  TriggerMenu["HLTTriggerBits"].setName("HLTTriggerBits");
  for (auto itr : hltPathNames){ 
    TriggerMenu["HLTTriggerBits"].newPath(itr); 
    TriggerMenu["HLTTriggerBits"].addTrigger( itr, trigger::trigger( itr, &passHLTPath[itr]) );
  }


  // ****************************************************************************************************
  // ****************************************************************************************************

  for (uint i = 0; i < selSampleStrs.size(); ++i){
      TString sample   = selSampleStrs[i];
      if ( sample.Contains( "PU40bx50" ) ){ bx50 = true; std::cout << "Using low-lumi L1 thresholds\n"; }

#ifdef ROOTFILE
      // Output files 
      TString fileName = outdir + fileSuffix + sample + ".root";
      fOut = new TFile( fileName ,"recreate");
#endif
      
      beginJob();
      analyse( sample );
      endJob();
  }

  // ****************************************************************************************************
  // *                                        Print menu results                                        *
  // ****************************************************************************************************

#ifdef PRINT
  // Verbose
  for( auto itr : SelectedTriggerMenues ){
    std::cout << "\n\n" << itr.first << "\n";
    TriggerMenu[itr].printPaths(false, false);
  }

  // Summary
  for( auto itr : SelectedTriggerMenues ){
    std::cout << "\n\n" << itr.first << "\n";
    TriggerMenu[itr].printPaths(false, true); //itr.second.printPaths(false, true);
  }
#endif

#ifdef MAKE_TRIGGER_LOG
  // Save data
  std::cout << "\n\nSaving files\n";
  //  for( auto itr : TriggerMenu ){ itr.second.savePaths( outdir + selSampleStr + "_" + fileSuffix ); }
  for( auto itr : SelectedTriggerMenues ){ TriggerMenu[itr].savePaths( outdir + selSampleStr + "_" + fileSuffix ); }
#endif





  exit(0);
}




void initialise(){

  fileSuffix = "";
#ifdef NOVETO
  fileSuffix += "Noveto_";
#endif


    // Choose offline reconstruction type
    offRecoType = "genAk4";
  
    offJet2PT   = offRecoType + "_Pt";
    offHT       = offRecoType + "_HT40";
    offAlphaT   = offRecoType + "_AlphaT40";
    offMHT      = offRecoType + "_MhtPT40";
    offMET      = "genMetCalo_MetPt";
    offForJetPT = offRecoType + "For_MaxPt";


    L1HTT     = "gct_Ht";
    L1ETM     = "gct_MetPt";
    
    hltPFJetPT = "hltAk4PF_Pt";


    hltCaloHT          = "hltAk4Calo_HT40";
    hltCaloAlphaT      = "hltAk4Calo_AlphaT40";
    hltCaloAlphaTPrime = "hltAk4Calo_AlphaTPrime40";

  jet2PTCuts.push_back(40.);
  jet2PTCuts.push_back(50.); 
  jet2PTCuts.push_back(60.);
  jet2PTCuts.push_back(70.);
  jet2PTCuts.push_back(80.);
  jet2PTCuts.push_back(90.);
  jet2PTCuts.push_back(100.);

  recoTypes.push_back("PF");
  //  recoTypes.push_back("Calo");


  // Analysis binning
  jetBinStrLabels[1]  = std::make_pair( "eq2j", " N_{jet} = 2 ");
  jetBinStrLabels[2]  = std::make_pair( "eq3j", " N_{jet} = 3 ");
  jetBinStrLabels[3]  = std::make_pair( "ge4j", " N_{jet} #geq 4 ");
  jetBinStrLabels[-1] = std::make_pair( "eq2a", " N_{jet}^{asym} = 2 ");
  jetBinStrLabels[-2] = std::make_pair( "eq3a", " N_{jet}^{asym} = 3 ");
  jetBinStrLabels[-3] = std::make_pair( "ge4a", " N_{jet}^{asym} #geq 4 ");

  HTBinStrLabels[0] = std::make_pair( "HT0", " #alpha_{T} > 0.65, 200 #leq H_{T} < 250 GeV");
  HTBinStrLabels[1] = std::make_pair( "HT1", " #alpha_{T} > 0.60, 250 #leq H_{T} < 300 GeV");
  HTBinStrLabels[2] = std::make_pair( "HT2", " #alpha_{T} > 0.55, 300 #leq H_{T} < 350 GeV");
  HTBinStrLabels[3] = std::make_pair( "HT3", " #alpha_{T} > 0.53, 350 #leq H_{T} < 400 GeV");
  HTBinStrLabels[4] = std::make_pair( "HT4", " #alpha_{T} > 0.52, 400 #leq H_{T} < 500 GeV");
  HTBinStrLabels[5] = std::make_pair( "HT5", " #alpha_{T} > 0.52, H_{T} #geq 500 GeV");


}

void beginJob(){






  int   HTPrefBins(16);
  float HTPrefMin(100), HTPrefMax(500);
  int alphaTPrefBins(40);
  float alphaTPrefMin(0.5), alphaTPrefMax(0.7);


  int   HTBins(40);
  float HTMin(0), HTMax(1000);
  int   MoHBins(20);
  float MoHMin(0.), MoHMax(1.);
  int alphaTBins(100), alphaTBinsSmall(30);
  float alphaTMin(0.5), alphaTMax(1.0);
  int jetBins(120);
  float jetMin(0), jetMax(600);

  double htAnaBinArr[]  = {200.,250.,300.,350.,400.,600.,800.,1000.};
  const Int_t htAnaBins = sizeof(htAnaBinArr)/sizeof(double) - 1; 


  double alphaTAnaBinArr[] = {0.50,0.51,0.52,0.53,0.54,0.55,0.56,0.57,0.58,0.59,
			   0.60,0.61,0.62,0.63,0.64,0.65,0.66,0.67,0.68,0.69,
			      0.70};
  const Int_t alphaTAnaBins = sizeof(alphaTAnaBinArr)/sizeof(double) - 1;




  // ************************************************************
  // Turn ons
  // ************************************************************
#ifdef TURNON


#ifdef PREFILTER
  // ------------------------------------------------------------ 
  // Prefilter 
  // ------------------------------------------------------------ 
  for (auto itr:turnOnPrefilterTriggers){
    TString menuName = itr;
    std::vector<TString> paths = TriggerMenu[menuName].getPathNames();

    // Loop through paths in menues 
    for (auto itr2:paths){
      TString pathName = itr2;
      TString prefix = "Pref_" + menuName + "_" + pathName + "_";

      for ( const auto &iJetBin : jetBinStrLabels ) {
	std::pair< TString, TString> jetBinStrLabel = iJetBin.second;
	TString jetStr = jetBinStrLabel.first + "_";
	TString jetLab = jetBinStrLabel.second;

	TString binStr = prefix + jetStr;
	TString binLab = " " + jetLab; 


	for (auto itr: Jet2PTVec){
	  TString Jet2PTStr = TString("Jet2gt") + TString(Form("%d", (itr)));
  
	  hist2DPrefTurn[binStr + Jet2PTStr]     = new TEfficiency( binStr + Jet2PTStr, binLab + ";Calo H_{T} (GeV);Calo #alpha_{T}", 
								    HTPrefBins, HTPrefMin, HTPrefMax, 
								    alphaTPrefBins, alphaTPrefMin, alphaTPrefMax);
	  // hist2DPrefTurn[binStr + "Dyn_" + Jet2PTStr] = new TEfficiency( binStr + "Dyn_" + Jet2PTStr, binLab + 
	  // 								 ";Calo H_{T} (GeV);Calo #alpha_{T}^{Dyn}", 
	  // 								 HTPrefBins, HTPrefMin, HTPrefMax, 
	  // 								 alphaTPrefBins, alphaTPrefMin, alphaTPrefMax);
	  hist2DPrefTurn[binStr + "Prime_" + Jet2PTStr] = new TEfficiency( binStr + "Prime_" + Jet2PTStr, binLab + 
									   ";Calo H_{T} (GeV);Calo #alpha_{T}'", 
									   HTPrefBins, HTPrefMin, HTPrefMax, 
									   alphaTPrefBins, alphaTPrefMin, alphaTPrefMax);
  
  
	}	
	for (auto itr: DijetAvgPTVec){
	  TString DijetAvgPTStr = TString("dijetAvggt") + TString(Form("%d", (itr)));

	  hist2DPrefTurn[binStr + DijetAvgPTStr]     = new TEfficiency( binStr + DijetAvgPTStr, binLab + ";Calo H_{T} (GeV);Calo #alpha_{T}",
                                                                    HTPrefBins, HTPrefMin, HTPrefMax, alphaTPrefBins, alphaTPrefMin, alphaTPrefMax);

	  hist2DPrefTurn[binStr + "Prime_" + DijetAvgPTStr] = new TEfficiency( binStr + "Prime_" + DijetAvgPTStr, 
									       binLab + ";Calo H_{T} (GeV);Calo #alpha_{T}",
                                                                    HTPrefBins, HTPrefMin, HTPrefMax, alphaTPrefBins, alphaTPrefMin, alphaTPrefMax);

        }



      }

    }
  }
 #endif


#ifdef HTTRIGGERS

  bool passedHT125(false);
  bool passedHT175(false);
  

  // L1 HT thresholds
  L1Trigger.push_back(0);
  L1Trigger.push_back(125);
  L1Trigger.push_back(175);

  // PF HT thresholds
  PFHTTrigger.push_back(200);  
  PFHTTrigger.push_back(250);
  PFHTTrigger.push_back(300);
  PFHTTrigger.push_back(350);
  PFHTTrigger.push_back(400);  
  PFHTTrigger.push_back(475);

  // Calo HT thresholds
  CaloHTTrigger.push_back(-9999);  
  CaloHTTrigger.push_back(-100);  
  CaloHTTrigger.push_back(-90);
  CaloHTTrigger.push_back(-80);
  CaloHTTrigger.push_back(-70);
  CaloHTTrigger.push_back(-60);  
  CaloHTTrigger.push_back(-50);
  CaloHTTrigger.push_back(-40);
  CaloHTTrigger.push_back(-30);
  CaloHTTrigger.push_back(-20);
  CaloHTTrigger.push_back(-10);
  CaloHTTrigger.push_back(0);

  
  for (uint iL1 = 0; iL1 < L1Trigger.size(); ++iL1){
    int L1HT = L1Trigger[iL1];
    TString L1Str = "L1HT" + TString(Form("%d",L1HT));
    
    for (uint iPF = 0; iPF < PFHTTrigger.size(); ++iPF){
      int PFHT = PFHTTrigger[iPF];
      TString PFHTStr = "PFHT" + TString(Form("%d", PFHT ));
     
      for (uint iCalo = 0; iCalo < CaloHTTrigger.size(); ++iCalo){
	int CaloHT = PFHT + CaloHTTrigger[iCalo];
	if (CaloHT < 0){ CaloHT = 0; }
	TString CaloHTStr = "CaloHT" + TString(Form("%d", CaloHT ));


	hist1DEff[L1Str + "_" + PFHTStr + "_" + CaloHTStr + "_GenHTTurnOn"]   = new TEfficiency( L1Str + "_" + PFHTStr + "_" + 
												 CaloHTStr + "_GenHTTurnOn",
												 ";Gen H_{T} (GeV);Efficiency", 60,0,600);

	std::cout << L1Str + "_" + PFHTStr + "_" + CaloHTStr + "_GenHTTurnOn" << "\n";
      }
    }
  }



  hist2DTurn["Calo_vs_Gen_HT125"] = new TEfficiency( "Calo_vs_Gen_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",
						     40,0,400, 40,0,400);

  hist2DTurn["Calo_vs_Gen_PFHT200_HT125"] = new TEfficiency( "Calo_vs_Gen_PFHT200_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_Gen_PFHT250_HT125"] = new TEfficiency( "Calo_vs_Gen_PFHT250_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_Gen_PFHT300_HT125"] = new TEfficiency( "Calo_vs_Gen_PFHT300_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_Gen_PFHT350_HT125"] = new TEfficiency( "Calo_vs_Gen_PFHT350_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_Gen_PFHT400_HT125"] = new TEfficiency( "Calo_vs_Gen_PFHT400_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_Gen_PFHT475_HT125"] = new TEfficiency( "Calo_vs_Gen_PFHT475_HT125", ";Gen H_{T} (GeV);Calo H_{T} (GeV)",40,0,400, 40,0,400);


  hist2DTurn["Calo_vs_PF_HT125"] = new TEfficiency( "Calo_vs_PF_HT125", ";PF H_{T} (GeV);Calo H_{T} (GeV)",
                                                    40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_PF_HT175"] = new TEfficiency( "Calo_vs_PF_HT175", ";PF H_{T} (GeV);Calo H_{T} (GeV)",
                                                    40,0,400, 40,0,400);


  hist2DTurn["Calo_vs_PF_HT125"] = new TEfficiency( "Calo_vs_PF_HT125", ";PF H_{T} (GeV);Calo H_{T} (GeV)",
						    40,0,400, 40,0,400);
  hist2DTurn["Calo_vs_PF_HT175"] = new TEfficiency( "Calo_vs_PF_HT175", ";PF H_{T} (GeV);Calo H_{T} (GeV)",
						    40,0,400, 40,0,400);

#endif


  for (auto itr:turnOnTriggers){
        
    TString menuName = itr;
    std::vector<TString> paths = TriggerMenu[menuName].getPathNames();

    for (auto itr2:paths){
      TString pathName = itr2;
      TString prefix = menuName + "_" + pathName + "_";

      // ********************************************************************************
      // 2D - AlphaT vs HT
      // ********************************************************************************


#ifdef HIST2DTURN
      hist2DTurn[prefix + "NoSel"]  = new TEfficiency( prefix + "NoSel", ";Gen H_{T} (GeV);Gen #alpha_{T}",
						      htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);
      hist2DTurn[prefix + "J2"]     = new TEfficiency( prefix + "J2",    ";Gen H_{T} (GeV);Gen #alpha_{T}",
						      htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);
      // hist2DTurn[prefix + "J2L1"]   = new TEfficiency( prefix + "J2L1",  ";Gen H_{T} (GeV);Gen #alpha_{T}",
      // 					         htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);

      
      // // Lepton - binned
      // hist2DTurn[prefix + "eq0Lep_" + "J2"] = new TEfficiency( prefix + "eq0Lep_" + "J2", "0 leptons;Gen H_{T} (GeV);Gen #alpha_{T}",
      // 							       htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);
      // hist2DTurn[prefix + "eq1Mu_" + "J2"]  = new TEfficiency( prefix + "eq1Mu_" + "J2", "1 Muon;Gen H_{T} (GeV);Gen #alpha_{T}",
      // 							       htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);
      // hist2DTurn[prefix + "ge2Mu_" + "J2"]  = new TEfficiency( prefix + "ge2Mu_" + "J2", "#geq 2 Muons;Gen H_{T} (GeV);Gen #alpha_{T}",
      // 							       htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);
      // hist2DTurn[prefix + "eq1Ele_" + "J2"] = new TEfficiency( prefix + "eq1Ele_" + "J2", "1 Electron;Gen H_{T} (GeV);Gen #alpha_{T}",
      // 							       htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);
      // hist2DTurn[prefix + "ge2Ele_" + "J2"] = new TEfficiency( prefix + "ge2Ele_" + "J2", "#geq 2 Electron;Gen H_{T} (GeV);Gen #alpha_{T}",
      // 							       htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);


      for ( const auto &iJetBin : jetBinStrLabels ) {
	std::pair< TString, TString> jetBinStrLabel = iJetBin.second;
	TString jetStr = jetBinStrLabel.first + "_";
	TString jetLab = jetBinStrLabel.second;

	TString binStr = prefix + jetStr;
	TString binLab = " " + jetLab; 

	hist2DTurn[binStr + "J2"]     = new TEfficiency( binStr + "J2", binLab + ";Gen H_{T} (GeV);Gen #alpha_{T}",
							 htAnaBins, htAnaBinArr, alphaTAnaBins, alphaTAnaBinArr);

      }
#endif



      // ********************************************************************************
      // 1D
      // ********************************************************************************
#ifdef HIST1DTURN
      hist1DEff[ prefix + "HT_TurnOn" ]            = new TEfficiency( prefix + "HT_TurnOn",     ";Gen H_{T} (GeV);Efficiency", 
								      HTBins,HTMin,HTMax);
      hist1DEff[ prefix + "HT_TurnOn_OffSel" ]     = new TEfficiency( prefix + "HT_TurnOn_OffSel", ";Gen H_{T} (GeV);Efficiency", 
								      HTBins,HTMin,HTMax);

      hist1DEff[ prefix + "HT_TurnOn_Alpha0p70" ]     = new TEfficiency( prefix + "HT_TurnOn_Alpha0p70", ";Gen H_{T} (GeV);Efficiency", 
								      HTBins,HTMin,HTMax);
      hist1DEff[ prefix + "HT_TurnOn_Alpha0p65" ]     = new TEfficiency( prefix + "HT_TurnOn_Alpha0p65", ";Gen H_{T} (GeV);Efficiency", 
								      HTBins,HTMin,HTMax);
      hist1DEff[ prefix + "HT_TurnOn_Alpha0p60" ]     = new TEfficiency( prefix + "HT_TurnOn_Alpha0p60", ";Gen H_{T} (GeV);Efficiency", 
								      HTBins,HTMin,HTMax);

      hist1DEff[prefix + "HT_TurnOn_Alpha0p70_Ana" ] = new TEfficiency( prefix + "HT_TurnOn_Alpha0p70_Ana", ";Gen H_{T} (GeV);Efficiency", 
									htAnaBins, htAnaBinArr);
      hist1DEff[prefix + "HT_TurnOn_Alpha0p65_Ana" ] = new TEfficiency( prefix + "HT_TurnOn_Alpha0p65_Ana", ";Gen H_{T} (GeV);Efficiency", 
									htAnaBins, htAnaBinArr);
      hist1DEff[prefix + "HT_TurnOn_Alpha0p60_Ana" ] = new TEfficiency( prefix + "HT_TurnOn_Alpha0p60_Ana", ";Gen H_{T} (GeV);Efficiency", 
									htAnaBins, htAnaBinArr);


      hist1DEff[ prefix + "AlphaT_TurnOn" ]        = new TEfficiency( prefix + "AlphaT_TurnOn", ";Gen #alpha_{T};Efficiency", 
								      alphaTBinsSmall,alphaTMin,alphaTMax);
      hist1DEff[ prefix + "AlphaT_TurnOn_OffSel" ] = new TEfficiency( prefix + "AlphaT_TurnOn_OffSel", ";Gen #alpha_{T};Efficiency", 
								      alphaTBinsSmall,alphaTMin,alphaTMax);

      hist1DEff[ prefix + "Jet2_TurnOn" ]          = new TEfficiency( prefix + "Jet2_TurnOn",   ";Gen p_{T}^{j2} (GeV);Efficiency", 
								      jetBins,jetMin,jetMax);
      hist1DEff[ prefix + "Jet2_TurnOn_OffSel" ]   = new TEfficiency( prefix + "Jet2_TurnOn_OffSel",   ";Gen p_{T}^{j2} (GeV);Efficiency", 
								      jetBins,jetMin,jetMax);
#endif      
    }
  }
#endif

#ifdef MAKE_TRIGGER_TUPLE
  for (uint iType = 0; iType < recoTypes.size(); ++iType){
    TString recoType = recoTypes[iType];

    for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
      float jet2PTCut = jet2PTCuts[ iJet2Cut ];
      TString jet2PTCutStr     = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
      TString jet2PTCutLab     = "p_{T}^{j2} > " + TString(Form("%1.0f", jet2PTCut ));
      TString dijetAvgPTCutLab = "p_{T}^{avg} > " + TString(Form("%1.0f", jet2PTCut ));
      
      TString stdStr = recoType + "_AlphaTStd_vs_HT_"  + jet2PTCutStr;
      TString aveStr = recoType + "_AlphaTStd_vs_HT_"  + jet2PTCutStr;
      TString dynStr = recoType + "_AlphaTDyn_vs_HT_"  + jet2PTCutStr;

      TString cStdStr = "Calo_AlphaTStd_vs_HT_"  + jet2PTCutStr;
      TString cDynStr = "Calo_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
      TString cPriStr = "Calo_AlphaTPri_vs_HT_"  + jet2PTCutStr;

      TString stdLab = "#alpha_{T}^{static} vs H_{T}";
      TString dynLab = "#alpha_{T}^{dynamic} vs H_{T}^{dynamic}";
      TString priLab = "#alpha_{T}' vs H_{T}";

      //#ifdef DYNALPHAT
      // Dynamic rate triggers
      //dynHistRate[dynStr]           = dynamicRate( HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax );
      hist2D[dynStr] = new TH2D(dynStr, recoType + "jet " + dynLab + " rate (" + 
				jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
				HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2D[cDynStr] = new TH2D(cDynStr, "Calojet " + dynLab + " rate (" + 
				 jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
				 HTBins,HTMin,HTMax,  alphaTPrefBins,alphaTPrefMin,alphaTPrefMax);

      //#endif
      
      
      hist2D[stdStr]            = new TH2D(stdStr, recoType + "jet " + stdLab + " rate (" + 
					   jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
					   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2D[cStdStr]           = new TH2D(cStdStr, "Calojet " + stdLab + " rate (" + 
					   jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
					   HTBins,HTMin,HTMax,  alphaTPrefBins,alphaTPrefMin,alphaTPrefMax);
      hist2D[cPriStr]           = new TH2D(cPriStr, "Calojet " + priLab + " rate (" + 
					   jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}'",  
					   HTBins,HTMin,HTMax,  alphaTPrefBins,alphaTPrefMin,alphaTPrefMax);

      // Calo rate after final trigger decision
      for (int i = 0; i < 5; ++i){
	TString WP = TString("PF") + TString(Form("%d", (i + 1)));
	TString DynAlphaTHTStr = WP + "DynHTAlphaT";
	TString cfPriStr = cPriStr + "_" + DynAlphaTHTStr;
	hist2D[cfPriStr]          = new TH2D(cfPriStr, "Calojet " + priLab + " rate (" + 
					     jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}'",  
					     HTBins,HTMin,HTMax,  alphaTPrefBins,alphaTPrefMin,alphaTPrefMax);
      }


      hist2D["Asym_" + stdStr]  = new TH2D("Asym_" + stdStr, recoType + "jet " + stdLab + " rate (" +
					   dijetAvgPTCutLab + ");H_{T} (GeV);#alpha_{T}",
					   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2D["Asym_" + cStdStr] = new TH2D("Asym_" + cStdStr, "Calojet " + stdLab + " rate (" +
					   dijetAvgPTCutLab + ");H_{T} (GeV);#alpha_{T}",
					   HTBins,HTMin,HTMax,  alphaTPrefBins,alphaTPrefMin,alphaTPrefMax);
      hist2D["Asym_" + cPriStr] = new TH2D("Asym_" + cPriStr, "Calojet " + priLab + " rate (" +
					   dijetAvgPTCutLab + ");H_{T} (GeV);#alpha_{T}'",
					   HTBins,HTMin,HTMax,  alphaTPrefBins,alphaTPrefMin,alphaTPrefMax);
      



      hist2DEff[ stdStr ] = new TEfficiency( stdStr, recoType + "jet " + stdLab + " efficiency (" + 
					     jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",
					     HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["NoL1_" + stdStr ] = new TEfficiency( "NoL1_" + stdStr, "NoL1 " + recoType + 
						      "jet " + stdLab + " efficiency (" + 
						      jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",
						      HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["Asym_" + stdStr ] = new TEfficiency( "Asym_" + stdStr, "Asymmetric " + recoType + 
						      "jet " + stdLab + " efficiency (" + 
						      dijetAvgPTCutLab + ");H_{T} (GeV);#alpha_{T}",
						      HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);

      // ------------------------------------------------------------------------------------------------------------------------ 
      // Analysis binned 
      // ------------------------------------------------------------------------------------------------------------------------
      for ( const auto &iJetBin : jetBinStrLabels ) {
	std::pair< TString, TString> jetBinStrLabel = iJetBin.second;
	TString jetStr = jetBinStrLabel.first;
	TString jetLab = jetBinStrLabel.second;

	for ( const auto &iHTBin : HTBinStrLabels ) {
	  std::pair< TString, TString> HTBinStrLabel = iHTBin.second;
	  TString HTStr = HTBinStrLabel.first;
	  TString HTLab = HTBinStrLabel.second;
  
	  TString binStr = "_" + HTStr + "_" + jetStr;
	  TString binLab = " " + HTLab + ", " + jetLab; 

#ifdef DYNALPHAT
	  // Dynamic rate triggers
	  dynHistEff[dynStr + binStr]  = dynamicRate( HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax );
	  // hist2DEff[ dynStr + binStr ] = new TEfficiency( dynStr + binStr, recoType + "jet " + dynLab + " efficiency (" + 
	  // 						    jet2PTCutLab + ")  - " + binLab + ";H_{T} (GeV);#alpha_{T}",
	  // 						  HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
#endif

	  hist2DEff[ stdStr + binStr ] = new TEfficiency( stdStr + binStr, recoType + "jet " + stdLab + " efficiency (" + 
							    jet2PTCutLab + ")  - " + binLab + ";H_{T} (GeV);#alpha_{T}",
							    HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
	  hist2DEff[ "NoL1_" + stdStr + binStr ] = new TEfficiency( "NoL1_" + stdStr + binStr, "NoL1 " + recoType + 
								    "jet " + stdLab + " efficiency (" + 
								    jet2PTCutLab + ")  - " + binLab + ";H_{T} (GeV);#alpha_{T}",
								    HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
	  hist2DEff[ "Asym_" + stdStr + binStr ] = new TEfficiency( "Asym_" + stdStr + binStr, "Asymmetric " + recoType + 
								    "jet " + stdLab + " efficiency (" + 
								    dijetAvgPTCutLab + ")  - " + binLab + ";H_{T} (GeV);#alpha_{T}",
								    HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);


	  
	}
      }



    }
  }
#endif


  return;
}
void endJob(){

#ifdef ROOTFILE
  fOut->mkdir("Rate");
  fOut->mkdir("Rate/Differential");
  fOut->mkdir("Efficiency");
  fOut->mkdir("Efficiency/Differential");
  fOut->mkdir("Efficiency/NoL1");
  fOut->mkdir("Efficiency/NoL1/Differential");
  fOut->mkdir("TurnOn");
  fOut->mkdir("TurnOn/Cumulative");
  fOut->mkdir("TurnOn2D");
  fOut->mkdir("TurnOn2D/LeptonBinned");
  fOut->mkdir("TurnOn2D/NJetBinned");
  fOut->mkdir("TurnOn2DPref");
  fOut->mkdir("TurnOn2DPref/NJetBinned");

#ifdef DYNALPHAT_NEW
  for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
    float jet2PTCut = jet2PTCuts[ iJet2Cut ];
    TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
    
    TString dynStr = "PF_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
    TH2D* histogram = (TH2D*)hist2D[dynStr];
    if (histogram == NULL){continue;}

    for (int Jet2 = minJet2; Jet2 < maxJet2; Jet2 += stepJet2){
      if (jet2PTCut == Jet2){
	for (int HT = minHT; HT < maxHT; HT += stepHT){
	  for (float AlphaT = minAlphaT; AlphaT < maxAlphaT; AlphaT += stepAlphaT){

	    TString HTAlphaTStr = TString(Form("%d", Jet2)) + "_" + TString(Form("%d", HT)) + "_" + TString(Form("%1.2f", AlphaT));
	    TString TriggerStr = "PF_DynNoPref_" + HTAlphaTStr;

	    // Fill histogram
	    histogram->Fill( HT + 1E-3, AlphaT + 1E-3, TriggerMenu[TriggerStr].getFinalPathRate() );
	    //	    std::cout << HTAlphaTStr << "\t" << TriggerMenu[TriggerStr].getFinalPathRate() << "\n";
	    
	  }
	}
      }
    }
  
  }
#endif


#ifdef DYNALPHAT
  // Dynamic rate - Convert to TH2 and feed to hist2D (Needs to be before hist2D routine)
  // ********************************************************************************
  for(std::map<TString, dynamicRate>::iterator itr = dynHistRate.begin(); itr != dynHistRate.end(); ++itr){
    dynamicRate dynamicHist = itr->second;
    TString     histoName   = itr->first;

    TH2D* cumulHist        = (TH2D*)hist2D[histoName];
    convertDynamicToTH2(  cumulHist,   dynamicHist, true );

  }

  // Dynamic efficiency - Convert to TEfficiency and feed to hist2DEff (Needs to be before hist2DEff routine) 
  // ********************************************************************************
  for(std::map<TString, dynamicRate>::const_iterator itr = dynHistEff.begin(); itr != dynHistEff.end(); ++itr){
    //dynamicRate dynamicHist = itr->second; 
    TString     histoName   = itr->first;
   
    dynamicRate dynPassed = dynHistEff[ histoName ];

    histoName = histoName.ReplaceAll("Dyn","Std");
    std::cout << histoName << "\n";
    TEfficiency *eff      = hist2DEff[ histoName ];
    TH2D* passed          = (TH2D*)eff->GetPassedHistogram()->Clone();
    //TH2D* total           = (TH2D*)eff->GetTotalHistogram() ->Clone() ;

    convertDynamicToTH2(  passed,   dynPassed );
    // fillUniform2D(        total,    dynPassed.timesSignalFired );
    // total->SetEntries( dynPassed.timesFired );
    // eff->SetTotalHistogram(  *total,  "" );
    // eff->SetPassedHistogram( *passed, "" );



    fOut->cd("Efficiency");
    if (dynPassed.timesSignalFired > 0){ passed->Scale( 1./dynPassed.timesSignalFired); }
    //if ( dynPassed.timesSignalFired == 0 ){ continue;}
    passed->Write();

  }
#endif



  std::cout << "\thist2D\n";
  for(std::map<TString, TH2*>::const_iterator itr = hist2D.begin(); itr != hist2D.end(); ++itr){
  
    TString histoName = itr->first; //Extract the histogram key 
    if (!histoName.Contains("AlphaTDyn")){
      fOut->cd("Rate/Differential");
      itr->second->Write();
    }

    fOut->cd("Rate");
    TH2* cumulHist = (TH2F*)itr->second->Clone();
    cumulHist->SetName( histoName + "_Cumul" );
    if (!histoName.Contains("AlphaTDyn")){
      reverseCumulative2D( itr->second, cumulHist, 1. );
    }
    // Cut unphysical region ht < 2*jet2PT 
    // if (histoName.Contains("OnRate_")){
    //   // Extract jet2PT cut 
    //   TString ptCutStr = histoName;
    //   ptCutStr = ptCutStr.Remove(0, ptCutStr.Index("Jet2gt") );
    //   ptCutStr.Remove( ptCutStr.Index("_"), ptCutStr.Length() );
    //   ptCutStr.ReplaceAll("Jet2gt","");
    //   float jet2PTCut = ptCutStr.Atof();
    //   clearRectangleX( cumulHist, jet2PTCut*2 );
    // }

    cumulHist->Write();

    //    delete itr->second;
    // delete cumulHist;

    itr->second->Delete();
  }
  hist2D.clear(); // Cleanup


  std::cout << "\thist2DEff\n";
  for ( const auto &itr : hist2DEff ) {


    TString histoName = itr.first; //Extract the histogram key 
    TString relDir = "";
    if (histoName.Contains("NoL1_")){ relDir = "NoL1/"; }
    fOut->cd("Efficiency/" + relDir + "Differential");
    itr.second->Write();

    fOut->cd("Efficiency/" + relDir);
    TEfficiency *eff      = (TEfficiency*)itr.second       ->Clone();
    TH2* passed           = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* total            = (TH2*)eff->GetTotalHistogram() ->Clone();
    TH2* passedCumul      = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* totalUniCumul    = (TH2*)eff->GetTotalHistogram() ->Clone();

    // Uniform 
    TEfficiency *effUniCumul = (TEfficiency*)eff->Clone();
    reverseCumulative2D( passed, passedCumul, 1 );

    // CHANGING TO NOT BE NORMALISED ANYMORE 
    fillUniform2D( totalUniCumul, total->GetEntries() );

    effUniCumul->SetTotalHistogram(  *totalUniCumul, "" );
    effUniCumul->SetPassedHistogram( *passedCumul,   "" );
    effUniCumul->SetName( histoName + "_UniCumul" );
    effUniCumul->Write();

  }



  std::cout << "\thist1DEff\n";
  for ( const auto &itr : hist1DEff ) {

    TString histoName = itr.first; //Extract the histogram key 
    fOut->cd("TurnOn");
    itr.second->Write();


    // Make cumulative efficiency 
    fOut->cd("TurnOn/Cumulative"); 
    TEfficiency* cumulHist = (TEfficiency*)itr.second->Clone( histoName + "_Cumul");
    //cumulHist->GetYaxis()->SetTitle("Cumulative efficiency");

    TH1* passed        = (TH1*)cumulHist->GetPassedHistogram()->Clone();
    TH1* passedCumul   = (TH1*)cumulHist->GetPassedHistogram()->Clone();
    TH1* total         = (TH1*)cumulHist->GetTotalHistogram() ->Clone();
    TH1* totalCumul    = (TH1*)cumulHist->GetTotalHistogram() ->Clone();
    reverseCumulative( passed, passedCumul, 1. );
    reverseCumulative( total,  totalCumul, 1. );
 
    cumulHist->SetTotalHistogram(  *totalCumul,  "" );
    cumulHist->SetPassedHistogram( *passedCumul, "" );
    cumulHist->Write();





  }


#ifdef HIST2DTURNSAVE
  std::cout << "\thist2DTurn\n";
  for ( const auto &itr : hist2DTurn ) {

    TString histoName = itr.first; //Extract the histogram key 
    fOut->cd("TurnOn2D");

    if (histoName.Contains("Ele") || histoName.Contains("Mu") || histoName.Contains("0Lep")){ fOut->cd("TurnOn2D/LeptonBinned"); }
    else if (histoName.Contains("_ge") || histoName.Contains("_eq"))                        { fOut->cd("TurnOn2D/NJetBinned"); }
    itr.second->Write();


    TEfficiency* cumulHist = (TEfficiency*)itr.second->Clone( histoName + "_Cumul");
    TH2* passed        = (TH2*)cumulHist->GetPassedHistogram()->Clone();
    TH2* passedCumul   = (TH2*)cumulHist->GetPassedHistogram()->Clone();
    TH2* total         = (TH2*)cumulHist->GetTotalHistogram() ->Clone();
    TH2* totalCumul    = (TH2*)cumulHist->GetTotalHistogram() ->Clone();
    
    if (histoName.Contains("Calo_vs_Gen")){
      reverseCumulative( passed, passedCumul, 1. );
      reverseCumulative( total,  totalCumul,  1. );
    }
    else{

      overflowX( passed );
      overflowX( total );
      reverseCumulativeY( passed, passedCumul, 1. );
      reverseCumulativeY( total,  totalCumul,  1. );
    } 
    

    cumulHist->SetTotalHistogram(  *totalCumul,  "" );
    cumulHist->SetPassedHistogram( *passedCumul, "" );
    cumulHist->Write();
  }
#endif





#ifdef PREFILTER
  std::cout << "\thist2DPrefTurn\n";
  for ( const auto &itr : hist2DPrefTurn ) {

    TString histoName = itr.first; //Extract the histogram key 
    fOut->cd("TurnOn2DPref");

    if (histoName.Contains("_ge") || histoName.Contains("_eq"))                        { fOut->cd("TurnOn2DPref/NJetBinned"); }
    itr.second->Write();

    TEfficiency* cumulHist = (TEfficiency*)itr.second->Clone( histoName + "_Cumul");
    TH2* passed        = (TH2*)cumulHist->GetPassedHistogram()->Clone();
    TH2* passedCumul   = (TH2*)cumulHist->GetPassedHistogram()->Clone();
    TH2* total         = (TH2*)cumulHist->GetTotalHistogram() ->Clone();
    TH2* totalCumul    = (TH2*)cumulHist->GetTotalHistogram() ->Clone();

    // Uniform                                                                                                                              
    reverseCumulative2D( passed, passedCumul, 1 );
    fillUniform2D( totalCumul, total->GetEntries() );
 

    cumulHist->SetTotalHistogram(  *totalCumul,  "" );
    cumulHist->SetPassedHistogram( *passedCumul, "" );
    cumulHist->Write();
  }
#endif




  fOut->Close();
#endif

  return;
}



void makeTriggerTuple(float weight, bool rate){
 


  for (uint iType = 0; iType < recoTypes.size(); ++iType){
    TString recoType = recoTypes[iType];

#ifdef DYNALPHAT
    // calculate all possible dynamic alphaT, HT combinations 
    std::vector<std::pair<float,float> > alphaTHTPairs;
    alphaTHTPairs = calculateDynamicAlphaTPairs( fvBranch["hltAk4" + recoType + "_Pt"], 
    						 fvBranch["hltAk4" + recoType + "_Px"], 
    						 fvBranch["hltAk4" + recoType + "_Py"],
    						 dynamicJetThreshold );
#endif


    for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
      float jet2PTCut = jet2PTCuts[ iJet2Cut ];
      TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
            
      TString stdStr = recoType + "_AlphaTStd_vs_HT_"  + jet2PTCutStr;
      TString dynStr = recoType + "_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
      TString cStdStr = "Calo_AlphaTStd_vs_HT_"  + jet2PTCutStr;
      TString cDynStr = "Calo_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
      TString cPriStr = "Calo_AlphaTPri_vs_HT_"  + jet2PTCutStr;

      // ************************************************************************************************** 
      // *                                              Rate                                              * 
      // ************************************************************************************************** 
      if (rate){
	if ( passL1HTTorETM ){

	  // ********************************************************************************
	  // Calo
	  // ********************************************************************************
	  if ( fBranch["hltAk4CaloSecond_Pt"] > jet2PTCut ){ 
	  
	    hist2D[cStdStr]->Fill( fBranch["hltAk4Calo_HT40"] , fBranch["hltAk4Calo_AlphaT40"], weight);
	    hist2D[cPriStr]->Fill( fBranch["hltAk4Calo_HT40"] , fBranch["hltAk4Calo_AlphaTPrime40"], weight);
	    

	    // Calo rate after final trigger decision
	    // **************************************************
	    for (int i = 0; i < 5; ++i){
	      TString WP = TString("PF") + TString(Form("%d", (i + 1)));
	      TString DynAlphaTHTStr = WP + "DynHTAlphaT";
	      if (passTrigger[DynAlphaTHTStr] && passTrigger[PF2Jet] ){
		TString cfPriStr = cPriStr + "_" + DynAlphaTHTStr;
		hist2D[cfPriStr]->Fill( fBranch["hltAk4Calo_HT40"] , fBranch["hltAk4Calo_AlphaT40"], weight);
	      }
	    }


	    // For validation
	    //dynHistRate[dynStr].triggerFired( fBranch["hltAk4Calo_HT40"] , fBranch["hltAk4Calo_AlphaT40"] );
#ifdef DYNALPHAT
	    // Fill dynamic triggers
	    for ( const auto itr : alphaTHTPairs ){ dynHistRate[dynStr].triggerFired( itr.second, itr.first ); }
#endif
	    
	  } // End 2Jet
	  if ( fBranch["hltAk4CaloDijetAvg_Pt"] > jet2PTCut ){
	    hist2D["Asym_" + cStdStr]->Fill( fBranch["hltAk4Calo_HT40"] , fBranch["hltAk4Calo_AlphaT40"],      weight);
	    hist2D["Asym_" + cPriStr]->Fill( fBranch["hltAk4Calo_HT40"] , fBranch["hltAk4Calo_AlphaTPrime40"], weight);


	  } // End dijetAvg
	  
	  // ********************************************************************************
	  // PF
	  // ********************************************************************************
	  if ( fBranch["hltAk4" + recoType + "Second_Pt"] > jet2PTCut ){ 
	  
	    hist2D[stdStr]->Fill( fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"], weight);
	    
	    // For validation
	    //dynHistRate[dynStr].triggerFired( fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"] );
#ifdef DYNALPHAT
	    // Fill dynamic triggers
	    for ( const auto itr : alphaTHTPairs ){ dynHistRate[dynStr].triggerFired( itr.second, itr.first ); }
#endif
	    
	  } // End 2Jet
	  if ( fBranch["hltAk4" + recoType + "DijetAvg_Pt"] > jet2PTCut ){
	    hist2D["Asym_" + stdStr]  ->Fill( fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"], weight);
	  } // End dijetAvg
	}

	continue;
      }

      //      if (recoType == "Calo"){ continue; }

      // ************************************************************************************************** 
      // *                                           Efficiency                                           * 
      // ************************************************************************************************** 



      if ( (fBranch["genAk4Lead_Pt"] > 100.) && (iBranch["genAk4_HTBin40"] > -1) && (iBranch["genAk4_NJetBin40"] != 0) && passOffVetoes ){

	// Determine offline bin
	TString jetStr = jetBinStrLabels[ iBranch["genAk4_NJetBin40"] ].first;
	TString HTStr  = HTBinStrLabels[  iBranch["genAk4_HTBin40"]   ].first;
	TString binStr = "_" + HTStr + "_" + jetStr;


	bool passJet2PT     = ( fBranch["hltAk4" + recoType + "Second_Pt"]   > jet2PTCut );
	bool passDijetAvgPT = ( fBranch["hltAk4" + recoType + "DijetAvg_Pt"] > jet2PTCut );
	bool passes      = ( passL1HTTorETM && passJet2PT );
	bool passesDijet = ( passL1HTTorETM && passDijetAvgPT );



	if (fBranch["genAk4Second_Pt"] > 40.){

#ifdef DYNALPHAT
	  // Dynamic efficiency
	  for ( const auto itr : alphaTHTPairs ){
	    if (passes){ dynHistEff[dynStr + binStr].triggerFired( itr.second, itr.first ); }
	    dynHistEff[dynStr + binStr].signalFired();
	  }
#endif
  
	  hist2DEff[stdStr]          ->Fill( passes, fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"]);
	  hist2DEff[stdStr + binStr] ->Fill( passes, fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"]);
	  hist2DEff["Asym_" + stdStr]         ->Fill( passesDijet, 
						      fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"]);
	  hist2DEff["Asym_" + stdStr + binStr]->Fill( passesDijet, 
						      fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"]);

	  if ( passL1HTTorETM ){ // L1 in denominator
	    hist2DEff["NoL1_" + stdStr]         ->Fill( passes, 
							fBranch["hltAk4" + recoType + "_HT40"] , 
							fBranch["hltAk4" + recoType + "_AlphaT40"]);
	    hist2DEff["NoL1_" + stdStr + binStr]->Fill( passes, 
							fBranch["hltAk4" + recoType + "_HT40"] , 
							fBranch["hltAk4" + recoType + "_AlphaT40"]);

	  }
	} 

      }


    } // End jet2
  }

#ifdef DYNALPHAT
  // Compute dynamic triggers
  for ( auto &itr : dynHistRate ){ itr.second.endEvent(); }
  for ( auto &itr : dynHistEff ) { itr.second.endEvent(); }
#endif  


} // End MakeTriggerTuple




// ************************************************************
// Turn ons
// ************************************************************
void makeTurnOns(){


#ifdef HTTRIGGERS



  for (uint iL1 = 0; iL1 < L1Trigger.size(); ++iL1){
    int L1HT = L1Trigger[iL1];
    TString L1Str = "L1HT" + TString(Form("%d",L1HT));
    // Make 125 turnon exclusive
    if ((L1HT == 125) && (fBranch["gct_Ht"] >= 175)){ continue; }

    for (uint iPF = 0; iPF < PFHTTrigger.size(); ++iPF){
      int PFHT = PFHTTrigger[iPF];
      TString PFHTStr = "PFHT" + TString(Form("%d", PFHT ));
     
      for (uint iCalo = 0; iCalo < CaloHTTrigger.size(); ++iCalo){
	int CaloHT = PFHT + CaloHTTrigger[iCalo];
	if (CaloHT < 0){ CaloHT = 0; }
	TString CaloHTStr = "CaloHT" + TString(Form("%d", CaloHT ));
	
	bool passedL1   = (fBranch["gct_Ht"]          >= L1HT );
	bool passedCalo = (fBranch["hltAk4Calo_HT40"] >= CaloHT );
	bool passedPF   = (fBranch["hltAk4PF_HT40"]   >= PFHT );
	bool passedHT   = passedL1&&passedCalo&&passedPF;

	hist1DEff[L1Str + "_" + PFHTStr + "_" + CaloHTStr + "_GenHTTurnOn"]   ->Fill( passedHT, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
      }
    }
  }

  return;


  bool passedHT125(false);
  bool passedHT175(false);

  passedHT125 = ( fBranch["gct_Ht"] >= 125.);
  passedHT175 = ( fBranch["gct_Ht"] >= 175.);


  hist2DTurn["Calo_vs_PF_HT125"] ->Fill( passedHT125, fBranch["hltAk4PF_HT40"], fBranch["hltAk4Calo_HT40"] );  
  hist2DTurn["Calo_vs_PF_HT175"] ->Fill( passedHT175, fBranch["hltAk4PF_HT40"], fBranch["hltAk4Calo_HT40"] );

  bool passedPFHT200 = ( fBranch["hltAk4PF_HT40"] >= 200.);
  bool passedPFHT250 = ( fBranch["hltAk4PF_HT40"] >= 250.);
  bool passedPFHT300 = ( fBranch["hltAk4PF_HT40"] >= 300.);
  bool passedPFHT350 = ( fBranch["hltAk4PF_HT40"] >= 350.);
  bool passedPFHT400 = ( fBranch["hltAk4PF_HT40"] >= 400.);
  bool passedPFHT475 = ( fBranch["hltAk4PF_HT40"] >= 475.);

  bool passedCaloHT125 = ( fBranch["hltAk4Calo_HT40"] >= 125.);
  bool passedCaloHT150 = ( fBranch["hltAk4Calo_HT40"] >= 150.);
  bool passedCaloHT175 = ( fBranch["hltAk4Calo_HT40"] >= 175.);
  bool passedCaloHT200 = ( fBranch["hltAk4Calo_HT40"] >= 200.);





    if ( passedHT125 ){


    hist2DTurn["Calo_vs_Gen_PFHT200_HT125"] ->Fill( passedPFHT200, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
    hist2DTurn["Calo_vs_Gen_PFHT250_HT125"] ->Fill( passedPFHT250, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
    hist2DTurn["Calo_vs_Gen_PFHT300_HT125"] ->Fill( passedPFHT300, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
    hist2DTurn["Calo_vs_Gen_PFHT350_HT125"] ->Fill( passedPFHT350, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
    hist2DTurn["Calo_vs_Gen_PFHT400_HT125"] ->Fill( passedPFHT400, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
    hist2DTurn["Calo_vs_Gen_PFHT475_HT125"] ->Fill( passedPFHT475, fBranch["genAk4_HT40"], fBranch["hltAk4Calo_HT40"] );
    }




#endif


  // ------------------------------------------------------------
  // Prefilter
  // ------------------------------------------------------------
#ifdef PREFILTER
  for (auto itr:turnOnPrefilterTriggers){
    TString menuName = itr;
    std::vector<TString> paths = TriggerMenu[menuName].getPathNames();

    // Loop through paths in menues
    for (auto itr2:paths){
      TString pathName = itr2;
      TString prefix = "Pref_" + menuName + "_" + pathName + "_";
      // Final decision of trigger path
      bool passed = TriggerMenu[menuName].getPath(pathName).passedAllTriggers();

      if (passOffVetoes){
	if ( ( iBranch["genAk4_NJetBin40"] != 0 ) ){ //Require jet selection
	  // Determine offline bin 
	  TString jetStr = jetBinStrLabels[ iBranch["genAk4_NJetBin40"] ].first + "_";
	  TString binStr = prefix + jetStr;

	  if (passed){ // Passed final filter
	    if (( iBranch["genAk4_HTBin40"] != -1 ) && (passOffVetoes) && ( iBranch["genAk4_NJetBin40"] != 0 ) ){  // Offline selection 
		
		// Dijet
		for (auto itr: Jet2PTVec){
		  TString Jet2PTStr = TString("Jet2gt") + TString(Form("%d", (itr)));
		  bool passCaloJet2 = false;
		  if (fBranch["hltAk4CaloSecond_Pt"]   > itr){ passCaloJet2 = true; }
		  hist2DPrefTurn[binStr + Jet2PTStr]           ->Fill( passCaloJet2, fBranch["hltAk4Calo_HT40"], fBranch["hltAk4Calo_AlphaT40"] );
		  //	      hist2DPrefTurn[binStr + "Dyn_" + Jet2PTStr]  ->Fill( passCaloJet2, fBranch["hltAk4Calo_HT40"], fBranch["hltAk4Calo_AlphaT40"] );
		  hist2DPrefTurn[binStr + "Prime_" + Jet2PTStr]->Fill( passCaloJet2, fBranch["hltAk4Calo_HT40"], fBranch["hltAk4Calo_AlphaTPrime40"] );
		  
		}
		
		// Dijet average
		for (auto itr: DijetAvgPTVec){
		  TString DijetAvgPTStr = TString("dijetAvggt") + TString(Form("%d", (itr)));
		  bool passCaloJet2 = false;
		  if (fBranch["hltAk4CaloDijetAvg_Pt"] > itr){ passCaloJet2 = true; }
		  hist2DPrefTurn[binStr + DijetAvgPTStr]           ->Fill( passCaloJet2, fBranch["hltAk4Calo_HT40"], fBranch["hltAk4Calo_AlphaT40"] );
		  hist2DPrefTurn[binStr + "Prime_" + DijetAvgPTStr]->Fill( passCaloJet2, fBranch["hltAk4Calo_HT40"], fBranch["hltAk4Calo_AlphaT40"] );
		}
		
		
		
	    } // End offline selection
	  }
	}
      }
    }
  }

#endif

  // ------------------------------------------------------------
  // Final filter
  // ------------------------------------------------------------
  for (auto itr:turnOnTriggers){
    TString menuName = itr;
    std::vector<TString> paths = TriggerMenu[menuName].getPathNames();

    // Loop through paths in menues
    for (auto itr2:paths){
      TString pathName = itr2;
      TString prefix = menuName + "_" + pathName + "_";
      // Final decision of trigger path
      bool passed = TriggerMenu[menuName].getPath(pathName).passedAllTriggers();

      
      // ******************************************************************************** 
      // 2D - AlphaT vs HT 
      // ******************************************************************************** 

      // Lepton bin
      TString lepStr = "eq0Lep_";
      if (iBranch["nIsoMuons"] > 0){
	if   (iBranch["nIsoMuons"] == 1)    { lepStr = "eq1Mu_"; }  else { lepStr = "ge2Mu_"; } 
      }
      else if (iBranch["nIsoElectrons"] > 0){
	if   (iBranch["nIsoElectrons"] == 1){ lepStr = "eq1Ele_"; } else { lepStr = "ge2Ele_"; } 
      }

      
      // Final trigger
      // --------------------
#ifdef HIST2DTURN
          

      if (passOffVetoes){
	hist2DTurn[prefix + "NoSel"] ->Fill( passed, fBranch["genAk4_HT40"], fBranch["genAk4_AlphaT40"] );

	if ( ( iBranch["genAk4_NJetBin40"] != 0 ) ){ //Require jet selection
	  // Determine offline bin 
	  TString jetStr = jetBinStrLabels[ iBranch["genAk4_NJetBin40"] ].first + "_";
	  TString binStr = prefix + jetStr;

	  hist2DTurn[prefix + "J2"]         ->Fill( passed, fBranch["genAk4_HT40"], fBranch["genAk4_AlphaT40"] );
	  hist2DTurn[binStr + "J2"]         ->Fill( passed, fBranch["genAk4_HT40"], fBranch["genAk4_AlphaT40"] );
	  
	  // Lepton binned
	  //hist2DTurn[prefix + lepStr + "J2"]->Fill( passed, fBranch["genAk4_HT40"], fBranch["genAk4_AlphaT40"] );
	  
	}
      }
#endif
#ifdef HIST1DTURN
      // ******************************************************************************** 
      // 1D 
      // ******************************************************************************** 
      if (passOffVetoes){
	// No selections
	hist1DEff[ prefix + "AlphaT_TurnOn" ]->Fill( passed, fBranch["genAk4_AlphaT40"] );
	hist1DEff[ prefix + "HT_TurnOn" ]    ->Fill( passed, fBranch["genAk4_HT40"] );
	hist1DEff[ prefix + "Jet2_TurnOn" ]  ->Fill( passed, fBranch["genAk4Second_Pt"] );
	
	// Loose offline selections
	if ((( fBranch["genAk4_HT40"] > 900) || (iBranch["genAk4_HTBin40"] > -1 )) ){
	  hist1DEff[ prefix + "AlphaT_TurnOn_OffSel" ]->Fill( passed, fBranch["genAk4_AlphaT40"] );
	}
	if (( (iBranch["genAk4_HTBin40"] > -1 )|| ( fBranch["genAk4_HT40"] > 900)) ){
	  hist1DEff[ prefix + "HT_TurnOn_OffSel" ]    ->Fill( passed, fBranch["genAk4_HT40"] );
	}
        if (( (iBranch["genAk4_HTBin40"] > -1 ) && ( iBranch["genAk4_NJetBin40"] != 0 ) ) ){
	  hist1DEff[ prefix + "Jet2_TurnOn_OffSel" ]  ->Fill( passed, fBranch["genAk4Second_Pt"] );
	}
	
      	//	0.70, 0.65, and 0.60
	if (fBranch["genAk4_AlphaT40"] > 0.70){
	  hist1DEff[ prefix +"HT_TurnOn_Alpha0p70"]->Fill( passed, fBranch["genAk4_HT40"] );
	  hist1DEff[ prefix +"HT_TurnOn_Alpha0p70_Ana"]->Fill( passed, fBranch["genAk4_HT40"] );
	}
	if (fBranch["genAk4_AlphaT40"] > 0.65){
	  hist1DEff[ prefix +"HT_TurnOn_Alpha0p65"]->Fill( passed, fBranch["genAk4HT40"] );
	  hist1DEff[ prefix +"HT_TurnOn_Alpha0p65_Ana"]->Fill( passed, fBranch["genAk4HT40"] );
	}
	if (fBranch["genAk4_AlphaT40"] > 0.60){
	  hist1DEff[ prefix +"HT_TurnOn_Alpha0p60"]->Fill( passed, fBranch["genAk4_HT40"] );
	  hist1DEff[ prefix +"HT_TurnOn_Alpha0p60_Ana"]->Fill( passed, fBranch["genAk4_HT40"] );
	}

      }
#endif
    }

  }
}

void analyse( TString sampleStr ){


  float jet2Thresh     = 90;
  float dijetAvgThresh = 90;



  std::cout << "\n\n--------------------------------------------------------------------------------\n\n";
  std::cout << "Processing: " << sampleStr << "\n\n";
  sample selectedSample = samples.samples[ sampleStr ]; 
  

  bool rate = false;
  if (sampleStr.Contains("QCD")){ rate = true; } 

    
  // Setup chain
  chain = selectedSample.chain;

   
  // Relax lepton veto for DY
  TString leptonVeto = "genLeptonVeto";
  if (selectedSample.name.Contains("DY_")){ leptonVeto = "genElectronVeto"; }
  
  
  // Load branches
  loadBranches( chain, fBranch, uBranch, iBranch, fvBranch, pfvBranch, leptonVeto, hltPathNames );

   
  int iEvent(0);     
  Int_t nEvents = (Int_t)chain->GetEntries();
  float weight  = (float)chain->GetWeight();

  // Set QCD rate weights
  // ------------------------------------------------------------
  for ( auto &itr : TriggerMenu ){ itr.second.setWeight( weight ); }
  //  for ( auto &itr : dynHistRate )          { itr.second.setWeight( weight ); }




    // ------------------------------------------------------------
    // Offline selections
    // ------------------------------------------------------------
  
    float jetThresh(40.);
    
    
  
    trigger::triggerPath offSelection;
    offSelection.addTrigger( trigger::trigger("OffJet2PT40", &passOffJet2PT40) );
    offSelection.addTrigger( trigger::trigger("OffJet2PT",   &passOffJet2PT)   );
    offSelection.addTrigger( trigger::trigger("OffHT200",    &passOffHT200)    );
    offSelection.addTrigger( trigger::trigger("OffAlphaT",   &passOffAlphaT)   );
    offSelection.addTrigger( trigger::trigger("OffLepton",   &passOffLepton)   );
    offSelection.addTrigger( trigger::trigger("offForJet",   &passOffForJet)   );
    //offSelection.addTrigger( trigger::trigger("OffMoM1p25",  &passOffMoM1p25)  );


  // ------------------------------------------------------------

  // ------------------------------------------------------------
  // Online selections
  // ------------------------------------------------------------

  // 25ns
  float l1HTTThresh = 175;
  float l1ETMThresh = 70;
  // Use lower L1 thresholds for low lumi
  if (bx50){
    l1HTTThresh = 125;
    l1ETMThresh = 60;
  }








  


  // ------------------------------------------------------------
  
  // --------------------------------------------------------------------------------
  // Loop over events
  // --------------------------------------------------------------------------------
  int eventLow  = 0;
  int eventHigh = nEvents;
#ifdef RUN_ON_BATCH
  eventLow               = floor( eventLowFact*nEvents) + eventLowOffset;
  eventHigh              = floor(eventHighFact*nEvents) + eventHighOffset;
#endif



  for( iEvent = eventLow; iEvent < eventHigh; ++iEvent ){
    if ( !(iEvent % 100000) ){ std::cout << "Event " << iEvent << "\n";}
    if ( ( MAX_EVENTS != -1) && ( iEvent > MAX_EVENTS ) ){ break; }
    chain->GetEntry(iEvent);



    // Reset variables
    // ----------------------------------------
    passOffJet2PT40 = false;
    passOffJet2PT   = false;
    passOffHT200    = false;   
    passOffAlphaT   = false;      
    passOffLepton   = false;  
    passOffForJet   = false;
    passOffMoM1p25  = false;  

    // Online
    passL1HTTorETM  = false;
    passHLTJet2PT40 = false;
    passHLTJet2PT   = false;
    passHLTDijetAvg = false;

    // Offline
    for (auto &itr: passOffHT ){ itr.second = false; }

    // ================================================================================
    // Perform trigger checks here
    // ================================================================================

    // ----------------------------------------
    // Offline
    // ----------------------------------------
    if (( fvBranch[ offJet2PT ]->size() > 1) && ( (*fvBranch[ offJet2PT ])[1] > jetThresh )){ passOffJet2PT40   = true;   // Jet2
      if ( (*fvBranch[ offJet2PT ])[0] > 100. )                                             { passOffJet2PT     = true; } // Jet2
    }
    if ( fBranch[ offHT ] > 200. )                                                          { passOffHT200      = true;   // HT
      if ( fBranch[ offAlphaT ] > ATThresh( fBranch[ offHT ] ) )                            { passOffAlphaT     = true; } // AlphaT
    }
    if ( !uBranch[ leptonVeto ])                                                            { passOffLepton     = true; } // Lepton
    if ( fBranch[ offForJetPT ] > jetThresh )                                               { passOffForJet     = true; } // ForJet
    //if ( (fBranch[offMET] > 0)&&( fBranch[offMHT]/fBranch[offMET] < 1.25) )                 { passOffMoM1p25    = true; } // 

#ifdef NOLEPTONVETO
    passOffLepton     = true;
#endif
    passOffVetoes = (passOffLepton && passOffForJet ); //&& passOffMoM1p25);
#ifdef NOVETO
    passOffVetoes = true;
#endif 


    // Reset triggers
    for( std::map<TString, bool >::iterator itr = passTrigger.begin(); itr != passTrigger.end(); ++itr){ 
      itr->second = false;
    }



    // Offline selections
    if ( iBranch["genAk4_HTBin40"] != -1 ){
      if ( !((iBranch["genAk4_HTBin40"] == 4) && (fBranch["genAk4_AlphaT40"] < 0.52)) ){
	if (passOffVetoes){
	  if ( iBranch["genAk4_NJetBin40"] > 0 ){  // Dijet selection
	    if ( iBranch["genAk4_HTBin40"] > -1 ){ passOffHT[ iBranch["genAk4_HTBin40"] ] = true; }
	  }
	  
	  // Pass offline selection in NJet bin
	  TString jetStr = jetBinStrLabels[ iBranch["genAk4_NJetBin40"] ].first;
	  TString HTStr  = HTBinStrLabels[  iBranch["genAk4_HTBin40"]   ].first;
	  TString offStr = HTStr + jetStr;
	  passTrigger[offStr] = true;
	}
      }
    }
        
    // ----------------------------------------
    // Online
    // ----------------------------------------
    if ( (fBranch[L1HTT] > l1HTTThresh ) || (fBranch[L1ETM] > l1ETMThresh ) )                   { passL1HTTorETM  = true; } // L1-seed
    if ( (fvBranch[ hltPFJetPT ]->size() > 1 ) && ( (*fvBranch[ hltPFJetPT ])[1] > jetThresh ) ){ passHLTJet2PT40 = true;   // Jet2
      if ( ( (*fvBranch[ hltPFJetPT ])[1] > jet2Thresh ) )                                      { passHLTJet2PT   = true; } // Jet2
      if ( 0.5*( (*fvBranch[hltPFJetPT])[0] + (*fvBranch[hltPFJetPT])[1] ) > dijetAvgThresh )   { passHLTDijetAvg = true; } // DijetAvg
    }



    passCaloDijet70    = (fBranch["hltAk4CaloSecond_Pt"] > 70);
    passCaloDijet60    = (fBranch["hltAk4CaloSecond_Pt"] > 60);

    passCaloPrefilter1 = ( passCaloDijet70 && (fBranch[hltCaloAlphaTPrime] > 0.54)  && (fBranch[hltCaloHT] > 150) );
    passCaloPrefilter2 = ( passCaloDijet70 && (fBranch[hltCaloAlphaTPrime] > 0.535) && (fBranch[hltCaloHT] > 200) );
    passCaloPrefilter3 = ( passCaloDijet70 && (fBranch[hltCaloAlphaTPrime] > 0.525) && (fBranch[hltCaloHT] > 250) );
    passCaloPrefilter4 = ( passCaloDijet70 && (fBranch[hltCaloAlphaTPrime] > 0.52)  && (fBranch[hltCaloHT] > 350) );
    passCaloPrefilter5 = ( passCaloDijet70 && (fBranch[hltCaloAlphaTPrime] > 0.51)  && (fBranch[hltCaloHT] > 375) );
    passCaloBrokenPrefilter1 = ( passCaloDijet70 && (fBranch[hltCaloAlphaT] > 0.54)  && (fBranch[hltCaloHT] > 150) );
    passCaloBrokenPrefilter2 = ( passCaloDijet70 && (fBranch[hltCaloAlphaT] > 0.535) && (fBranch[hltCaloHT] > 200) );
    passCaloBrokenPrefilter3 = ( passCaloDijet70 && (fBranch[hltCaloAlphaT] > 0.525) && (fBranch[hltCaloHT] > 250) );
    passCaloBrokenPrefilter4 = ( passCaloDijet70 && (fBranch[hltCaloAlphaT] > 0.52)  && (fBranch[hltCaloHT] > 350) );
    passCaloBrokenPrefilter5 = ( passCaloDijet70 && (fBranch[hltCaloAlphaT] > 0.51)  && (fBranch[hltCaloHT] > 375) );

    

    TString recoType = "PF";

    // calculate all possible dynamic alphaT, HT combinations 
    std::vector<std::pair<float,float> > alphaTHTPairs;
    alphaTHTPairs = calculateDynamicAlphaTPairs( fvBranch["hltAk4" + recoType + "_Pt"], 
    						 fvBranch["hltAk4" + recoType + "_Px"], 
    						 fvBranch["hltAk4" + recoType + "_Py"],
    						 dynamicJetThreshold );






    // Check calo dijet
    // ****************************************
    for (auto itr: Jet2PTVec){
      TString cJet2PTStr = TString("cJet2gt") + TString(Form("%d", (itr)));
      if (fBranch["hltAk4CaloSecond_Pt"] > itr){ passTrigger[ cJet2PTStr ] = true; }
      else{ break; }
    }


    // ********************************************************************************
    // 2012 thresholds
    // ********************************************************************************

    // Dijet prefilter
    // passTrigger[ "p12012J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]      > 70.);
    // passTrigger[ "p12012J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 150.); 
    // passTrigger[ "p12012J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.54 ); 
    passTrigger[ "p12012J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p12012J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 175.); 
    passTrigger[ "p12012J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.59 ); 


    passTrigger[ "p22012J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p22012J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 225.); 
    passTrigger[ "p22012J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.53 ); 

    passTrigger[ "p32012J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p32012J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 275.); 
    passTrigger[ "p32012J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.525 ); 

    passTrigger[ "p42012J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p42012J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 325.); 
    passTrigger[ "p42012J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.515 ); 

    passTrigger[ "p52012J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p52012J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 375.); 
    passTrigger[ "p52012J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 

    // Asymmetric dijet prefilter
    // passTrigger[ "p12012Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    // passTrigger[ "p12012HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 150.); 
    // passTrigger[ "p12012AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.54 ); 
    passTrigger[ "p12012Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p12012HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 175.); 
    passTrigger[ "p12012AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.59 ); 

		        
    passTrigger[ "p22012Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p22012HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 225.); 
    passTrigger[ "p22012AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.53 ); 
		        
    passTrigger[ "p32012Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p32012HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 275.); 
    passTrigger[ "p32012AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.525 ); 
		        
    passTrigger[ "p42012Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p42012HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 325.); 
    passTrigger[ "p42012AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.515 ); 
		        
    passTrigger[ "p52012Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p52012HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 375.); 
    passTrigger[ "p52012AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 


    // ********************************************************************************
    // New thresholds
    // ********************************************************************************

    // Dijet prefilter
    passTrigger[ "p1J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p1J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 175.); 
    passTrigger[ "p1J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.59 ); 

    passTrigger[ "p2J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p2J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 225.); 
    passTrigger[ "p2J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.55 ); 

    passTrigger[ "p3J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p3J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 250.); 
    passTrigger[ "p3J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.53 ); 

    passTrigger[ "p4J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p4J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 300.); 
    passTrigger[ "p4J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 

    passTrigger[ "p5J2Jet2" ] = (fBranch["hltAk4CaloSecond_Pt"]    > 70.);
    passTrigger[ "p5J2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 325.); 
    passTrigger[ "p5J2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 

    // Asymmetric dijet prefilter
    passTrigger[ "p1Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p1HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 175.); 
    passTrigger[ "p1AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.59 ); 

    passTrigger[ "p2Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 225.); 
    passTrigger[ "p2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.55 ); 

    passTrigger[ "p3Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p3HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 250.); 
    passTrigger[ "p3AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.53 ); 

    passTrigger[ "p4Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p4HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 300.); 
    passTrigger[ "p4AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 

    passTrigger[ "p5Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 70.);
    passTrigger[ "p5HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 325.); 
    passTrigger[ "p5AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 



  
  // passTrigger[ "p1Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 40.);
  //   passTrigger[ "p1HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 150.); 
  //   passTrigger[ "p1AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.55 ); 

  //   passTrigger[ "p2Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 60.);
  //   passTrigger[ "p2HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 200.); 
  //   passTrigger[ "p2AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.52 ); 

  //   passTrigger[ "p3Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 40.);
  //   passTrigger[ "p3HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 225.); 
  //   passTrigger[ "p3AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.52 ); 

  //   passTrigger[ "p4Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 40.);
  //   passTrigger[ "p4HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 275.); 
  //   passTrigger[ "p4AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.51 ); 

  //   passTrigger[ "p5Jet2" ] = (fBranch["hltAk4CaloDijetAvg_Pt"]    > 40.);
  //   passTrigger[ "p5HT" ]   = (fBranch["hltAk4Calo_HT40"]          > 300.); 
  //   passTrigger[ "p5AT" ]   = (fBranch["hltAk4Calo_AlphaTPrime40"] > 0.505 ); 
  



    passTrigger[ "CaloHT175" ]   = ( fBranch["hltAk4Calo_HT40"] > 175. );
    passTrigger[ "CaloHT200" ]   = ( fBranch["hltAk4Calo_HT40"] > 200. );
    passTrigger[ "CaloATP0p52" ] = ( fBranch["hltAk4Calo_AlphaTPrime40"] > 0.52 );
    passTrigger[ "CaloATP0p53" ] = ( fBranch["hltAk4Calo_AlphaTPrime40"] > 0.53 );


    // Check dijet
    // ****************************************
    for (auto itr: Jet2PTVec){
      TString Jet2PTStr = TString("Jet2gt") + TString(Form("%d", (itr)));
      if (fBranch["hltAk4" + recoType + "Second_Pt"] > itr){ passTrigger[ Jet2PTStr ] = true; }
      else{ break; }
    }

    // Check dijet average
    // ****************************************
    for (auto itr: DijetAvgPTVec){
      TString DijetAvgPTStr = TString("dijetAvggt") + TString(Form("%d", (itr)));
      if (fBranch["hltAk4" + recoType + "DijetAvg_Pt"] > itr){ passTrigger[ DijetAvgPTStr ] = true; }
      else{ break; }
    }


    // Check HLT AlphaT triggers
    for( std::map<TString, std::pair<float,float> >::iterator itr = HLTHTAlphaTTrigger.begin(); itr != HLTHTAlphaTTrigger.end(); ++itr){ 
      TString WP                      = itr->first;
      std::pair<float,float> HTAlphaT = itr->second;

      TString HTStr     = WP + "HT"; 
      TString AlphaTStr = WP + "AlphaT";
      TString DynAlphaTStr   = WP + "DynAlphaT";
      TString DynAlphaTHTStr = WP + "DynHTAlphaT";
      TString hltHT(""), hltAlphaT(""), hltJet2PT("");
      if ( WP.Contains("PF") ) { hltHT = "hltAk4PF_HT40"; hltAlphaT = "hltAk4PF_AlphaT40"; hltJet2PT = "hltAk4PF_Second_Pt";  }
      else                     { std::cout << "Error: Unknown HLT reconstruction.\n";  continue; }//exit(0); }

      if ( fBranch[ hltHT ]        > HLTHTAlphaTTrigger[WP].first  )                          { passTrigger[ HTStr ]        = true; }
      if ( fBranch[ hltHT ]        > HLTHTAlphaTTrigger[WP].first - 25  )                     { passTrigger[ HTStr + "Min25" ] = true; }
      if ( fBranch[ hltAlphaT ]    > HLTHTAlphaTTrigger[WP].second )                          { passTrigger[ AlphaTStr ]    = true; 
     
	// Scanning higher AlphaT thresholds
	for (int dAlphaT = 1; dAlphaT < 5; ++dAlphaT){
	  float aTStep = 0.005;
	  TString aTOffset = "";
	  if (WP.Contains("1")){ aTStep = 0.01; }
	  if ( fBranch[ hltAlphaT ]    > HLTHTAlphaTTrigger[WP].second + dAlphaT*aTStep ){ // AlphaT
	    aTOffset = TString(Form("%1.3f", dAlphaT*aTStep )); aTOffset.ReplaceAll(".","p");
	    passTrigger[ AlphaTStr + aTOffset ] = true;
	  }
	  else{ break; }
	}	      
      	
      }



      
      // Dynamic triggers
      for ( const auto &itr: alphaTHTPairs ){ 
	if ( itr.first > HLTHTAlphaTTrigger[WP].second ){ // AlphaT max
	  passTrigger[ DynAlphaTStr ] = true; 
	}

      	if ( itr.second > HLTHTAlphaTTrigger[WP].first - 25 ){ // HT
      	  if ( itr.first > HLTHTAlphaTTrigger[WP].second ){ // AlphaT
      	    passTrigger[ DynAlphaTHTStr + "Min25" ] = true; 
	    
      	    if ( itr.second > HLTHTAlphaTTrigger[WP].first ){ // HT
      	      passTrigger[ DynAlphaTHTStr ] = true; 
	      
	      // Scanning higher AlphaT thresholds
	      for (int dAlphaT = 1; dAlphaT < 5; ++dAlphaT){
		float aTStep = 0.005;
		TString aTOffset = "";
		if (WP.Contains("1")){ aTStep = 0.01; }
		if ( itr.first > HLTHTAlphaTTrigger[WP].second + dAlphaT*aTStep ){ // AlphaT
		  aTOffset = TString(Form("%1.3f", dAlphaT*aTStep )); aTOffset.ReplaceAll(".","p");
		  passTrigger[ DynAlphaTHTStr + aTOffset ] = true;
		}
		else{ break; }
		
	      }	      
      	      //break;
      	    }
      	  }
      	}
      }

#ifdef ALPHATDYN_NEW
      for (int Jet2 = minJet2; Jet2 < maxJet2; Jet2 += stepJet2){

	if ( (fvBranch[ hltPFJetPT ]->size() < 2 ) || ((*fvBranch[ hltPFJetPT ])[1] < Jet2) ){ break; } // Jet2
	else{

	  for ( const auto &itr: alphaTHTPairs ){ 
	    for (int HT = minHT; HT < maxHT; HT += stepHT){
	      if ( itr.second > HT ){  // HT

		for (float AlphaT = minAlphaT; AlphaT < maxAlphaT; AlphaT += stepAlphaT){
		  if ( itr.first > AlphaT ){  // AlphaT

		    TString HTAlphaTStr = TString(Form("%d", Jet2)) + "_" + TString(Form("%d", HT)) + "_" + TString(Form("%1.2f", AlphaT));
		    passTrigger[ HTAlphaTStr] = true;
		    
		  } 
		}
	      }
	    }
	  }
	}
      }
#endif

    } // End HLTHTAlphaTTrigger

    // ----------------------------------------
    // Trigger bits
    // ----------------------------------------
    for (auto itr : hltPathNames){ passHLTPath[itr] = uBranch[itr]; }




    // ----------------------------------------
    // End check
    // ----------------------------------------

    offSelection.checkEvent();
    // Check menues
    for ( auto &itr : TriggerMenu ){ itr.second.checkEvent(); }






    // ********************************************************************************
    // Additional function calls
    // ********************************************************************************
#ifdef MAKE_TRIGGER_TUPLE
    makeTriggerTuple( chain->GetWeight(), rate );
#endif
#ifdef TURNON
    makeTurnOns();
#endif
    // ********************************************************************************

    // ================================================================================

  } // End loop


  // ----------------------------------------
  // Print logs
  // ----------------------------------------

  std::cout << "\n\n--------------------------------------------------------------------------------\n\n";
  std::cout << "Sample = " << selectedSample.name << "\n\n";

  offSelection.print();

#ifdef MAKE_TRIGGER_LOG
  if (!selectedSample.name.Contains("QCD")){
    // Save data 
    std::cout << "\n\nSaving files\n";
    //for( auto &itr : TriggerMenu ){ 
    for( auto itr : SelectedTriggerMenues ){
      TriggerMenu[itr].savePaths( outdir + selSampleStr + "_" + selectedSample.name + "_" + fileSuffix ); 
      TriggerMenu[itr].clearPaths();

      // itr.second.savePaths( outdir + selSampleStr + "_" + selectedSample.name + "_" + fileSuffix ); 
      // itr.second.clearPaths();
 
    }
  }
#endif

#ifdef ALPHATDYN_NEW
//   // Dynamic trigger scan
//   for (int Jet2 = minJet2; Jet2 < maxJet2; Jet2 += stepJet2){
//     for (int HT = minHT; HT < maxHT; HT += stepHT){
//       for (float AlphaT = minAlphaT; AlphaT < maxAlphaT; AlphaT += stepAlphaT){
// 	TString HTAlphaTStr = TString(Form("%d", Jet2)) + "_" + TString(Form("%d", HT)) + "_" + TString(Form("%1.2f", AlphaT));
// 	TString TriggerStr = "PF_DynNoPref_" + HTAlphaTStr;
// #ifdef MAKE_TRIGGER_LOG
// 	TriggerMenu[TriggerStr].savePaths( outdir + "DYN_" + selSampleStr + "_" + selectedSample.name + "_" + fileSuffix );
// #endif
// 	TriggerMenu[TriggerStr].clearPaths();
//       }
//     }
//   }
#endif

  

}







//  LocalWords:  Str
