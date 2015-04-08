#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
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

#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/interface/dynamic.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/typeDefs.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/TriggerTuple/DynAlphaT.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/TriggerTuple/branches.cpp"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/TriggerTuple/integration.cpp"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/TriggerTuple/trigger.cpp"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_4_0_pre9/src/AlphaTHLT/Macros/TriggerTuple/samples/samples_v1.cpp"

const int MAX_EVENTS = 100000;
//const int MAX_EVENTS = -1;

// #define MAKE_TRIGGER_TUPLE
// #define DYNALPHAT

//TString outdir = "output/05_04_15/";
//TString outdir = "output/06_04_15/";
TString outdir = "output/08_04_15/";
//TString outdir = "output/test/";




// Containers 
std::map< TString, TH1*> hist;

std::map< TString, TH1*> hist1D;
std::map< TString, TH2*> hist2D;
std::map< TString, TH1*> hist1DRate;
std::map< TString, TH2*> hist2DRate;
std::map< TString, TEfficiency*> hist1DEff;
std::map< TString, TEfficiency*> hist2DEff;

// Dynamic triggers 
std::map< TString, dynamicRate > dynHistRate;
std::map< TString, dynamicRate > dynHistEff;



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
// *                                           Tree skimmer                                           *
// ****************************************************************************************************

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

// int getAnalysisBin(float offHT, float offAlphaT, float offJet1PT, float offJet2PT){

//   int bin(0), asymSign(1);

//   // Second jet threshold
//   if ( (offJet1PT > offJet1PTThreshold) && (offJet2PT > offJetThreshold) ){
//     if (offJet2PT < offJet2PTThreshold){ asymSign = -1; } // Asymmetric bin
//   }
//   else{ return 0; }

//   // HT-AlphaT thresholds
//   if      ( (offHT >= 200.) && (offHT < 250.)   ){ if (offAlphaT > 0.65) bin = 1; }
//   else if ( (offHT >= 250.) && (offHT < 300.)   ){ if (offAlphaT > 0.60) bin = 2; }
//   else if ( (offHT >= 300.) && (offHT < 350.)   ){ if (offAlphaT > 0.55) bin = 3; }
//   else if ( (offHT >= 350.) && (offHT < 400.)   ){ if (offAlphaT > 0.53) bin = 4; }
//   else if ( (offHT >= 400.) && (offHT < 500.)   ){ if (offAlphaT > 0.52) bin = 5; }
//   else if ( (offHT >= 500.) && (offHT < 99999.) ){ if (offAlphaT > 0.51) bin = 6; }

//   // Encode asymmetric bin information
//   bin *= asymSign;

//   return bin;

// }

inline float ATThresh(float offHT){
  if      ( offHT >= 500.){ return 0.51;}
  else if ( offHT >= 400.){ return 0.52;}
  else if ( offHT >= 350.){ return 0.53;}
  else if ( offHT >= 300.){ return 0.55;}
  else if ( offHT >= 250.){ return 0.60;}
  else if ( offHT >= 200.){ return 0.65;}
  return 99999;
}




  std::map<TString, bool>  passTrigger;
  std::map<TString, std::pair<float,float> > HLTHTAlphaTTrigger;

trigger     HTTrigger(TString WP, TString suffix = "" ){
    TString HTStr     = WP + "HT" + suffix;
    passTrigger[ HTStr] = false;
    return trigger( HTStr,      &passTrigger[ HTStr] ); 
  }
  trigger AlphaTTrigger(TString WP ){
    TString AlphaTStr = WP + "AlphaT";
    passTrigger[ AlphaTStr] = false;
    return trigger( AlphaTStr,  &passTrigger[ AlphaTStr] );
  }
trigger DynHTAlphaTTrigger(TString WP, TString suffix = "" ){
    TString AlphaTStr = WP + "DynHTAlphaT" + suffix;
    passTrigger[ AlphaTStr] = false;
    return trigger( AlphaTStr,  &passTrigger[ AlphaTStr] );
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




std::map< TString, triggerPath > HLTTriggerBits;


std::map< TString, triggerMenu > TriggerMenu;

triggerMenu PFBrokenHLTTriggers;

std::map< TString, triggerPath > PFHLTTriggers;
std::map< TString, triggerPath > PFHLTMin25Triggers;
std::map< TString, triggerPath > PFDynHLTTriggers;
std::map< TString, triggerPath > PFDynHLTMin25Triggers;
std::map< TString, triggerPath > PFDijetAvgHLTTriggers;
std::map< TString, triggerPath > PFDijetAvgDynHLTMin25Triggers;

std::map< TString, trigger > CaloHLTTriggers;
std::map< TString, trigger > CaloBrokenHLTTriggers;

// HLT paths                                                                                                                            
std::vector<TString>    hltPathNames;
std::map<TString, bool> passHLTPath;



void initialise();
void beginJob();
void endJob();
void analyse(TString sampleStr );



// ********************************************************************************
//
// ********************************************************************************

void triggerTuple(){



  // Load all samples
  loadSamples( samples, sampleStrs );



#ifdef BATCH
  selSampleStrs = sampleStrs["BATCH"];
#else
  //selSampleStrs = sampleStrs["PU40bx50_HCAL3_HPUV_QCD"];
  //selSampleStrs = sampleStrs["PU40bx50_HPUV_QCD"];  
  //selSampleStrs = sampleStrs["PU40bx25_HCAL3_HPUV_QCD"];
  //selSampleStrs = sampleStrs["PU40bx25_HPUV_QCD"];


  //selSampleStrs = sampleStrs["SM"];
  //selSampleStrs = sampleStrs["PU20bx25_HCAL3_Signal"];
  //selSampleStrs = sampleStrs["PU20bx25_Signal"];
    selSampleStrs = sampleStrs["test"];
    // selSampleStrs = sampleStrs["PU40bx25_HCAL3_HPUV_QCD"];
#endif
  initialise();

  // ****************************************************************************************************
  // *                                          Define triggers                                         *
  // ****************************************************************************************************
  

  triggerPath HLTBaseSelection;
  HLTBaseSelection.addTrigger( trigger("L1HTTorETM",  &passL1HTTorETM)  );
  HLTBaseSelection.addTrigger( trigger("HLTJet2PT40", &passHLTJet2PT40) );
  HLTBaseSelection.addTrigger( trigger("HLTJet2PT",   &passHLTJet2PT)   );

  triggerPath HLTBaseDijetAvgSelection;
  HLTBaseDijetAvgSelection.addTrigger( trigger("L1HTTorETM",    &passL1HTTorETM)  );
  HLTBaseDijetAvgSelection.addTrigger( trigger("HLTJet2PT40",   &passHLTJet2PT40) );
  HLTBaseDijetAvgSelection.addTrigger( trigger("HLTDijetAvgPT", &passHLTDijetAvg) );


  // ------------------------------
  // Prefilters
  // ------------------------------
  CaloHLTTriggers["PF1"] =        trigger("CaloPrefilter1",        &passCaloPrefilter1 );
  CaloHLTTriggers["PF2"] =        trigger("CaloPrefilter2",        &passCaloPrefilter2 );
  CaloHLTTriggers["PF3"] =        trigger("CaloPrefilter3",        &passCaloPrefilter3 );
  CaloHLTTriggers["PF4"] =        trigger("CaloPrefilter4",        &passCaloPrefilter4 );
  CaloHLTTriggers["PF5"] =        trigger("CaloPrefilter5",        &passCaloPrefilter5 );
  CaloBrokenHLTTriggers["PF1"] =  trigger("CaloBrokenPrefilter1",  &passCaloBrokenPrefilter1 );
  CaloBrokenHLTTriggers["PF2"] =  trigger("CaloBrokenPrefilter2",  &passCaloBrokenPrefilter2 );
  CaloBrokenHLTTriggers["PF3"] =  trigger("CaloBrokenPrefilter3",  &passCaloBrokenPrefilter3 );
  CaloBrokenHLTTriggers["PF4"] =  trigger("CaloBrokenPrefilter4",  &passCaloBrokenPrefilter4 );
  CaloBrokenHLTTriggers["PF5"] =  trigger("CaloBrokenPrefilter5",  &passCaloBrokenPrefilter5 );

  // ------------------------------
  // PF triggers
  // ------------------------------
  for (int i = 0; i < 5; ++i){
    TString WP = TString("PF") + TString(Form("%d", (i + 1)));

    TriggerMenu["PFBrokenHLTTrigger"].addPath(    WP, HLTBaseSelection);
    TriggerMenu["PFBrokenHLTTrigger"].addTrigger( WP, CaloBrokenHLTTriggers[ WP ] );
    TriggerMenu["PFBrokenHLTTrigger"].addTrigger( WP,            HTTrigger( WP ) );
    TriggerMenu["PFBrokenHLTTrigger"].addTrigger( WP,        AlphaTTrigger( WP ) );
    TriggerMenu["PFBrokenHLTTrigger"].setSignal(  WP, passOffHT[i] );    


    PFBrokenHLTTriggers.addPath(    WP, HLTBaseSelection);
    PFBrokenHLTTriggers.addTrigger( WP, CaloBrokenHLTTriggers[ WP ] );
    PFBrokenHLTTriggers.addTrigger( WP,            HTTrigger( WP ) );
    PFBrokenHLTTriggers.addTrigger( WP,        AlphaTTrigger( WP ) );
    PFBrokenHLTTriggers.setSignal( WP, passOffHT[i] );    


    // PFBrokenHLTTriggers[ WP ] = HLTBaseSelection;
    // PFBrokenHLTTriggers[ WP ].addTrigger( CaloBrokenHLTTriggers[ WP ] );
    // PFBrokenHLTTriggers[ WP ].addTrigger(             HTTrigger( WP ) );
    // PFBrokenHLTTriggers[ WP ].addTrigger(         AlphaTTrigger( WP ) );
    // PFBrokenHLTTriggers[ WP ].setSignal( passOffHT[i] );    

    PFHLTTriggers[ WP ] = HLTBaseSelection;
    PFHLTTriggers[ WP ].addTrigger( CaloHLTTriggers[ WP ] );
    PFHLTTriggers[ WP ].addTrigger(       HTTrigger( WP ) );
    PFHLTTriggers[ WP ].addTrigger(   AlphaTTrigger( WP ) );
    PFHLTTriggers[ WP ].setSignal( passOffHT[i] );

    PFDynHLTTriggers[ WP ] = HLTBaseSelection;
    PFDynHLTTriggers[ WP ].addTrigger(    CaloHLTTriggers[ WP ] );
    PFDynHLTTriggers[ WP ].addTrigger( DynHTAlphaTTrigger( WP ) );
    PFDynHLTTriggers[ WP ].setSignal( passOffHT[i] );

    PFHLTMin25Triggers[ WP ] = HLTBaseSelection;
    PFHLTMin25Triggers[ WP ].addTrigger( CaloHLTTriggers[ WP ] );
    PFHLTMin25Triggers[ WP ].addTrigger(       HTTrigger( WP, "Min25" ) );
    PFHLTMin25Triggers[ WP ].addTrigger(   AlphaTTrigger( WP ) );
    PFHLTMin25Triggers[ WP ].setSignal( passOffHT[i] );

    PFDynHLTMin25Triggers[ WP ] = HLTBaseSelection;
    PFDynHLTMin25Triggers[ WP ].addTrigger(    CaloHLTTriggers[ WP ] );
    PFDynHLTMin25Triggers[ WP ].addTrigger( DynHTAlphaTTrigger( WP, "Min25") );
    PFDynHLTMin25Triggers[ WP ].setSignal( passOffHT[i] );

    PFDijetAvgHLTTriggers[ WP ] = HLTBaseDijetAvgSelection;
    PFDijetAvgHLTTriggers[ WP ].addTrigger( CaloHLTTriggers[ WP ] );
    PFDijetAvgHLTTriggers[ WP ].addTrigger(       HTTrigger( WP, "Min25" ) );
    PFDijetAvgHLTTriggers[ WP ].addTrigger(   AlphaTTrigger( WP ) );
    PFDijetAvgHLTTriggers[ WP ].setSignal( passOffHT[i] );

    PFDijetAvgDynHLTMin25Triggers[ WP ] = HLTBaseDijetAvgSelection;
    PFDijetAvgDynHLTMin25Triggers[ WP ].addTrigger(    CaloHLTTriggers[ WP ] );
    PFDijetAvgDynHLTMin25Triggers[ WP ].addTrigger( DynHTAlphaTTrigger( WP, "Min25") );
    PFDijetAvgDynHLTMin25Triggers[ WP ].setSignal( passOffHT[i] );

  }

  // //PFHLTTriggers[ "PF0" ].setSignal( passOffHT0 );
  // PFHLTTriggers[ "PF1" ].setSignal( passOffHT1 );
  // PFHLTTriggers[ "PF2" ].setSignal( passOffHT2 );
  // PFHLTTriggers[ "PF3" ].setSignal( passOffHT3 );
  // PFHLTTriggers[ "PF4" ].setSignal( passOffHT4 );
  // PFHLTTriggers[ "PF5" ].setSignal( passOffHT5 );




  // HLT trigger bits
  hltPathNames.push_back("HLT_PFHT200_DiPFJet90_PFAlphaT0p57_v1");
  hltPathNames.push_back("HLT_PFHT250_DiPFJet90_PFAlphaT0p55_v1");
  hltPathNames.push_back("HLT_PFHT300_DiPFJet90_PFAlphaT0p53_v1");
  hltPathNames.push_back("HLT_PFHT350_DiPFJet90_PFAlphaT0p52_v1");
  hltPathNames.push_back("HLT_PFHT400_DiPFJet90_PFAlphaT0p51_v1");

  hltPathNames.push_back("HLT_PFHT900_v1");
  hltPathNames.push_back("HLT_PFHT350_PFMET120_NoiseCleaned_v1");
  hltPathNames.push_back("HLT_PFMET170_NoiseCleaned_v1");
  hltPathNames.push_back("HLT_PFMET120_NoiseCleaned_BTagCSV07_v1");

  hltPathNames.push_back("HLT_HT100_v1");
  hltPathNames.push_back("HLT_HT200_v1");
  hltPathNames.push_back("HLT_HT250_v1");
  hltPathNames.push_back("HLT_HT300_v1");
  hltPathNames.push_back("HLT_HT350_v1");
  hltPathNames.push_back("HLT_HT400_v1");

  hltPathNames.push_back("HLT_PFHT100_v1");
  hltPathNames.push_back("HLT_PFHT350_v1");
  hltPathNames.push_back("HLT_PFHT600_v1");
  hltPathNames.push_back("HLT_PFHT650_v1");

  hltPathNames.push_back("HLT_Rsq0p36_v1");
  hltPathNames.push_back("HLT_RsqMR260_Rsq0p09_MR200_v1");
  hltPathNames.push_back("HLT_RsqMR260_Rsq0p09_MR200_4jet_v1");
  hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_v1");
  hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_4jet_v1");

  for (auto itr : hltPathNames){ HLTTriggerBits[itr].addTrigger( trigger( itr,    &passHLTPath[ itr ])  ); };



  // ****************************************************************************************************
  // ****************************************************************************************************

  for (uint i = 0; i < selSampleStrs.size(); ++i){
      TString sample   = selSampleStrs[i];
      if ( sample.Contains( "PU40bx50" ) ){ bx50 = true; std::cout << "Using low-lumi L1 thresholds\n"; }

      // Output files 
      TString fileName = outdir + sample + ".root";
      fOut = new TFile( fileName ,"recreate");
      
      
      beginJob();
      analyse( sample );
      endJob();
  }

  // ****************************************************************************************************
  // ****************************************************************************************************



  // CONVERT TO A TRIGGER MENU

  std::cout << "\n\nPFBrokenHLTTriggers\n";
  //  for( auto itr : PFBrokenHLTTriggers ) { itr.second.print(); }
  TriggerMenu["PFBrokenHLTTrigger"].printPaths();

  PFBrokenHLTTriggers.printPaths();
  // std::cout << "\n\nPFHLTTriggers\n";
  // for( auto itr : PFHLTTriggers)        { itr.second.print(); }
  // std::cout << "\n\nPFDynHLTTriggers\n";
  // for( auto itr : PFDynHLTTriggers)     { itr.second.print(); }
  // std::cout << "\n\nPFHLTMin25Triggers\n";
  // for( auto itr : PFHLTMin25Triggers)   { itr.second.print(); }
  // std::cout << "\n\nPFDynHLTTriggers\n";
  // for( auto itr : PFDynHLTMin25Triggers){ itr.second.print(); }
  // std::cout << "\n\nPFDijetAvgHLTTriggers\n";
  // for( auto itr : PFDijetAvgHLTTriggers){ itr.second.print(); }
  // std::cout << "\n\nPFDijetAvgDynHLTTriggers\n";
  // for( auto itr : PFDijetAvgDynHLTMin25Triggers){ itr.second.print(); }

  // std::cout << "\nTriggerBits\n";
  // for( auto itr : HLTTriggerBits){ itr.second.print(false); }






  exit(0);
}



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


void initialise(){
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
  //recoTypes.push_back("Calo");


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
  HTBinStrLabels[5] = std::make_pair( "HT5", " #alpha_{T} > 0.51, H_{T} #geq 500 GeV");


}

void beginJob(){






#ifdef MAKE_TRIGGER_TUPLE
  int   HTBins(40);
  float HTMin(0), HTMax(1000);
  int   MoHBins(20);
  float MoHMin(0.), MoHMax(1.);
  int alphaTBins(120);
  float alphaTMin(0.4), alphaTMax(1.0);









  for (uint iType = 0; iType < recoTypes.size(); ++iType){
    TString recoType = recoTypes[iType];

    for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
      float jet2PTCut = jet2PTCuts[ iJet2Cut ];
      TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
      TString jet2PTCutLab = "p_{T}^{j2} > " + TString(Form("%1.0f", jet2PTCut ));
      
      
      TString stdStr = recoType + "_AlphaTStd_vs_HT_"  + jet2PTCutStr;
      TString dynStr = recoType + "_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
      TString stdLab = "#alpha_{T}^{static} vs H_{T}^{static}";
      TString dynLab = "#alpha_{T}^{dynamic} vs H_{T}^{dynamic}";


#ifdef DYNALPHAT
      // Dynamic rate triggers
      dynHistRate[dynStr]           = dynamicRate( HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax );
      hist2D[dynStr] = new TH2D(dynStr, recoType + "jet " + dynLab + " rate (" + 
				jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
				HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
#endif
      
      
      hist2D[stdStr] = new TH2D(stdStr, recoType + "jet " + stdLab + " rate (" + 
				jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
				HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);




      hist2DEff[ stdStr ] = new TEfficiency( stdStr, recoType + "jet " + stdLab + " efficiency (" + 
					     jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",
					     HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["NoL1_" + stdStr ] = new TEfficiency( "NoL1_" + stdStr, "NoL1 " + recoType + 
						      "jet " + stdLab + " efficiency (" + 
						      jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",
						      HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
      hist2DEff["Asym_" + stdStr ] = new TEfficiency( "Asym_" + stdStr, "Asymmetric " + recoType + 
						      "jet " + stdLab + " efficiency (" + 
						      jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",
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
	  hist2DEff[ dynStr + binStr ] = new TEfficiency( dynStr + binStr, recoType + "jet " + dynLab + " efficiency (" + 
							    jet2PTCutLab + ")  - " + binLab + ";H_{T} (GeV);#alpha_{T}",
							  HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
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
								    jet2PTCutLab + ")  - " + binLab + ";H_{T} (GeV);#alpha_{T}",
								    HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);


	  
	}
      }



    }
  }
#endif


  return;
}
void endJob(){

#ifdef MAKE_TRIGGER_TUPLE
  fOut->mkdir("Rate");
  fOut->mkdir("Rate/Differential");
  fOut->mkdir("Efficiency");
  fOut->mkdir("Efficiency/Differential");
  fOut->mkdir("Efficiency/NoL1");
  fOut->mkdir("Efficiency/NoL1/Differential");


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
    if ( dynPassed.timesSignalFired == 0 ){ continue;}
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

      // ************************************************************************************************** 
      // *                                              Rate                                              * 
      // ************************************************************************************************** 
      if (rate){
	if ( passL1HTTorETM ){
	  if ( fBranch["hltAk4" + recoType + "Second_Pt"] > jet2PTCut ){ 
	    
	    hist2D[stdStr]->Fill( fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"], weight);
	    
	    // For validation
	    //dynHistRate[dynStr].triggerFired( fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"] );
#ifdef DYNALPHAT
	    // Fill dynamic triggers
	    for ( const auto itr : alphaTHTPairs ){ dynHistRate[dynStr].triggerFired( itr.second, itr.first ); }
#endif
	    
	  }
	}

	continue;
      }

      // ************************************************************************************************** 
      // *                                           Efficiency                                           * 
      // ************************************************************************************************** 



      if ( (fBranch["genAk4Lead_Pt"] > 100.) && (iBranch["genAk4_HTBin40"] > -1) && (iBranch["genAk4_NJetBin40"] != 0) && passOffVetoes ){

	// Determine offline bin
	TString jetStr = jetBinStrLabels[ iBranch["genAk4_NJetBin40"] ].first;
	TString HTStr  = HTBinStrLabels[  iBranch["genAk4_HTBin40"]   ].first;
	TString binStr = "_" + HTStr + "_" + jetStr;


	bool passJet2PT = ( fBranch["hltAk4" + recoType + "Second_Pt"] > jet2PTCut );
	bool passes = ( passL1HTTorETM && passJet2PT );

	if (fBranch["genAk4Second_Pt"] > 100.){ // Normal dijet

#ifdef DYNALPHAT
	  //if ( (iBranch["genAk4_NJetBin40"] != 1)||(iBranch["genAk4_HTBin40"] != 0) ){ continue; } 

	  // Dynamic efficiency
	  for ( const auto itr : alphaTHTPairs ){
	    if (passes){ dynHistEff[dynStr + binStr].triggerFired( itr.second, itr.first ); }
	    dynHistEff[dynStr + binStr].signalFired();
	  }
#endif
  
	  hist2DEff[stdStr]         ->Fill( passes, fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"]);
	  hist2DEff[stdStr + binStr]->Fill( passes, fBranch["hltAk4" + recoType + "_HT40"] , fBranch["hltAk4" + recoType + "_AlphaT40"]);

	  if ( passL1HTTorETM ){ // L1 in denominator
	    hist2DEff["NoL1_" + stdStr]         ->Fill( passes, 
							fBranch["hltAk4" + recoType + "_HT40"] , 
							fBranch["hltAk4" + recoType + "_AlphaT40"]);
	    hist2DEff["NoL1_" + stdStr + binStr]->Fill( passes, 
							fBranch["hltAk4" + recoType + "_HT40"] , 
							fBranch["hltAk4" + recoType + "_AlphaT40"]);





	  }
	}
	else{ //Asymmetric dijet

	  hist2DEff["Asym_" + stdStr]         ->Fill( passes, 
						      fBranch["hltAk4" + recoType + "_HT40"] , 
						      fBranch["hltAk4" + recoType + "_AlphaT40"]);
          hist2DEff["Asym_" + stdStr + binStr]->Fill( passes, 
						      fBranch["hltAk4" + recoType + "_HT40"] , 
						      fBranch["hltAk4" + recoType + "_AlphaT40"]);
	  
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


void analyse( TString sampleStr ){


  float jet2Thresh     = 90;
  float dijetAvgThresh = 90;



  std::cout << "\n\n--------------------------------------------------------------------------------\n\n";
  std::cout << "Processing: " << sampleStr << "\n\n";
  sample selectedSample = samples.samples[ sampleStr ]; //samples["QCD40bx25_PT800to1000"] ;
  

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
    PFBrokenHLTTriggers.setWeight( weight );

    //    for ( auto &itr : PFBrokenHLTTriggers ) { itr.second.setWeight( weight ); }
    for ( auto &itr : PFHLTTriggers )       { itr.second.setWeight( weight ); }
    for ( auto &itr : PFDynHLTTriggers)     { itr.second.setWeight( weight ); }
    for ( auto &itr : PFHLTMin25Triggers)   { itr.second.setWeight( weight ); }
    for ( auto &itr : PFDynHLTMin25Triggers){ itr.second.setWeight( weight ); }
    for ( auto &itr : PFDijetAvgHLTTriggers){ itr.second.setWeight( weight ); }
    for ( auto &itr : PFDijetAvgDynHLTMin25Triggers){ itr.second.setWeight( weight ); }

    for ( auto &itr : HLTTriggerBits)        { itr.second.setWeight( weight ); }
    for ( auto &itr : dynHistRate )          { itr.second.setWeight( weight ); }




    // ------------------------------------------------------------
    // Offline selections
    // ------------------------------------------------------------
  
    float jetThresh(40.);
    
    
  
    triggerPath offSelection;
    offSelection.addTrigger( trigger("OffJet2PT40", &passOffJet2PT40) );
    offSelection.addTrigger( trigger("OffJet2PT",   &passOffJet2PT)   );
    offSelection.addTrigger( trigger("OffHT200",    &passOffHT200)    );
    offSelection.addTrigger( trigger("OffAlphaT",   &passOffAlphaT)   );
    offSelection.addTrigger( trigger("OffLepton",   &passOffLepton)   );
    offSelection.addTrigger( trigger("offForJet",   &passOffForJet)   );
    offSelection.addTrigger( trigger("OffMoM1p25",  &passOffMoM1p25)  );


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




  HLTHTAlphaTTrigger["PF1"] = std::make_pair( 200., 0.57 );
  HLTHTAlphaTTrigger["PF2"] = std::make_pair( 250., 0.55 );
  HLTHTAlphaTTrigger["PF3"] = std::make_pair( 300., 0.53 );
  HLTHTAlphaTTrigger["PF4"] = std::make_pair( 350., 0.52 );
  HLTHTAlphaTTrigger["PF5"] = std::make_pair( 400., 0.51 );




  


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
    if ( (fBranch[offMET] > 0)&&( fBranch[offMHT]/fBranch[offMET] < 1.25) )                 { passOffMoM1p25    = true; } // 
    passOffVetoes = (passOffLepton && passOffForJet && passOffMoM1p25);

    // Offline selections
    if ( iBranch["genAk4_HTBin40"] != -1 ){

      if (passOffVetoes){
	if ( iBranch["genAk4_NJetBin40"] > 0 ){  // Dijet selection

	  if ( iBranch["genAk4_HTBin40"] > -1 ){
	    passOffHT[ iBranch["genAk4_HTBin40"] ] = true;
	  }

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


#ifdef MAKE_TRIGGER_TUPLE
    makeTriggerTuple( chain->GetWeight(), rate );
#endif


    // Reset triggers
    for( std::map<TString, bool >::iterator itr = passTrigger.begin(); itr != passTrigger.end(); ++itr){ 
      itr->second = false;
    }





    // Check HLT AlphaT triggers
    for( std::map<TString, std::pair<float,float> >::iterator itr = HLTHTAlphaTTrigger.begin(); itr != HLTHTAlphaTTrigger.end(); ++itr){ 
      TString WP                      = itr->first;
      std::pair<float,float> HTAlphaT = itr->second;

      TString HTStr     = WP + "HT"; 
      TString AlphaTStr = WP + "AlphaT";
      TString DynAlphaTStr = WP + "DynHTAlphaT";
      TString hltHT(""), hltAlphaT("");
      if ( WP.Contains("PF") ) { hltHT = "hltAk4PF_HT40"; hltAlphaT = "hltAk4PF_AlphaT40";  }
      else                     { std::cout << "Error: Unknown HLT reconstruction.\n"; exit(0); }

      if ( fBranch[ hltHT ]        > HLTHTAlphaTTrigger[WP].first  )                          { passTrigger[ HTStr ]        = true; }
      if ( fBranch[ hltHT ]        > HLTHTAlphaTTrigger[WP].first - 25  )                     { passTrigger[ HTStr + "Min25" ] = true; }
      if ( fBranch[ hltAlphaT ]    > HLTHTAlphaTTrigger[WP].second )                          { passTrigger[ AlphaTStr ]    = true; }

      // Dynamic
      for ( const auto &itr: alphaTHTPairs ){ 
	if ( itr.first > HLTHTAlphaTTrigger[WP].second - 25 ){ // HT
	  if ( itr.second > HLTHTAlphaTTrigger[WP].first ){ // AlphaT
	    passTrigger[ DynAlphaTStr + "Min25" ] = true; 
	    
	    if ( itr.first > HLTHTAlphaTTrigger[WP].second ){ // HT
	      passTrigger[ DynAlphaTStr ] = true; 
	      break;
	    }
	  }
	}
      }
      

    }


    // Trigger bits
    for (auto itr : hltPathNames){ passHLTPath[itr] = uBranch[itr]; }
    //HLTTriggerBits[itr].addTrigger( trigger( itr,    &passHLTPath[ itr ])  ); };



    // ----------------------------------------
    // End check
    // ----------------------------------------

    offSelection.checkEvent();
    // PF trigger menu
    //for (auto &itr : PFBrokenHLTTriggers)  { itr.second.checkEvent(); }
    for ( auto &itr : TriggerMenu ){ itr.second.checkEvent(); }

    PFBrokenHLTTriggers.checkEvent();
    for (auto &itr : PFHLTTriggers)        { itr.second.checkEvent(); }
    for (auto &itr : PFDynHLTTriggers)     { itr.second.checkEvent(); }
    for (auto &itr : PFHLTMin25Triggers)   { itr.second.checkEvent(); }
    for (auto &itr : PFDynHLTMin25Triggers){ itr.second.checkEvent(); }
    for (auto &itr : PFDijetAvgHLTTriggers){ itr.second.checkEvent(); }
    for (auto &itr : PFDijetAvgDynHLTMin25Triggers){ itr.second.checkEvent(); }
    for (auto &itr : HLTTriggerBits)       { itr.second.checkEvent(); }
    // ================================================================================

  } // End loop


  // ----------------------------------------
  // Print logs
  // ----------------------------------------

  std::cout << "\n\n--------------------------------------------------------------------------------\n\n";
  std::cout << "Sample = " << selectedSample.name << "\n\n";

  offSelection.print();

  

  

}






