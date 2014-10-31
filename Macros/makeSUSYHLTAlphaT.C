#define TEST
#define SIGNAL
//#define NEUTRINO
//#define HLT_CALOJET

//#define MATCHING


// ********************************************************************************
// *                                  Selections                                  *
// ********************************************************************************

float caloJetThreshold = 40; //50;
float genJetThreshold  = 40; //50;
float secondJetThreshold = 90;


// Dynamic alphaT variables
int maxCaloJet              = 15;
float caloJetDynThreshold   = caloJetThreshold;  
float caloJetAlphaThreshold = 0.50;  


float hltPFSecondJetThreshold = 90;


// offline selection for measuring calojet efficiency 
float genJet2PTThreshold    = 100;

const int MAX_JETS = 4;

// Jet matching
float maxDeltaR = 0.2;

// ********************************************************************************

unsigned int maxEvents(100000);


const int N_50NS_BUNCHES = 1368;
const int N_25NS_BUNCHES = 2508;
const int N_MAX_BUNCHES  = 3564;
const double CROSSING_TIME  = 25E-9;
const int LHC_FREQUENCY  = 11246;


const double INST_LUMI_25NS = 1.4E34;
const double INST_LUMI_50NS = 0.7E34;
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

#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_2_0_pre8/src/AlphaTHLT/Macros/AlphaT.h"
#include "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_2_0_pre8/src/AlphaTHLT/Macros/JetMatch.h"


// Function prototypes
void reverseCumulative( TH1* histogram, TH1* rCumulHist, double scale);
void reverseCumulative2D( TH2* histogram, TH2* rCumulHist, double scale);
void fillUniform2D( TH2* histogram, double value);
void clearTriangle( TH1* histogram, int sign, float extra );
void clearRectangleX( TH1* histogram, float cut );

void fillUniformY( TH2* histogram, TH2* rCumulHist);
void reverseCumulativeY( TH2* histogram, TH2* rCumulHist, double scale);




struct sample{
  TString Name;
  TChain* Chain;
  double  CrossSection;

  sample(TString name, TString branch, TString files, double xs):Name(name),CrossSection(xs){
    Chain = new TChain( branch );
    Chain->Add( files );
  }
};

void makePTHatBinnedHist( std::map<TString, TH1*> &hist, TString name ){
  
  TString ptStr = name + "_PTHat";

  hist[ptStr] = new TH1D( ptStr, name + " #hat{p}_{T} binned;Entries;Fired", 11, 0, 11);
  hist[ptStr]->GetXaxis()->SetBinLabel( 1,  "Not fired");
  hist[ptStr]->GetXaxis()->SetBinLabel( 2,  "Fired");
  hist[ptStr]->GetXaxis()->SetBinLabel( 3,  "QCD30to50");
  hist[ptStr]->GetXaxis()->SetBinLabel( 4,  "QCD50to80");
  hist[ptStr]->GetXaxis()->SetBinLabel( 5,  "QCD80to120");
  hist[ptStr]->GetXaxis()->SetBinLabel( 6,  "QCD120to170");
  hist[ptStr]->GetXaxis()->SetBinLabel( 7,  "QCD170to300");
  hist[ptStr]->GetXaxis()->SetBinLabel( 8,  "QCD300to470");
  hist[ptStr]->GetXaxis()->SetBinLabel( 9,  "QCD470to600");
  hist[ptStr]->GetXaxis()->SetBinLabel( 10, "QCD600to800");
  hist[ptStr]->GetXaxis()->SetBinLabel( 11, "QCD800to1000");
  
}



void makeSUSYHLTAlphaT(){

  TChain* rateChain = new TChain("MakeTrees/Ntuple");
  TChain* sigChain  = new TChain("MakeTrees/Ntuple");


  // ------------------------------------------------------------------------------------------------------------------------
  // Samples
  // ------------------------------------------------------------------------------------------------------------------------
  TString branch    = "MakeTrees/Ntuple";

  TString sampleDir = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/17Oct14/";  // pre8 with correct JEC, 25ns
  //  TString sampleDir = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/22Oct14_50ns/";  // pre8 with correct JEC (needs updating with latest GT), 50ns

  // Automatically determine bunch spacing
  double instLumi(0);
  int    nBunches(0);
  if (sampleDir.Contains("50ns")){ 
    std::cout << "\n\nSample BX = 50ns\n\n";
    instLumi = INST_LUMI_50NS;
    nBunches = N_50NS_BUNCHES;
  }
  else{
    std::cout << "\n\nSample BX = 25ns\n\n";
    instLumi = INST_LUMI_25NS; 
    nBunches = N_25NS_BUNCHES;
  }

  const int ZB_XSECTION    = LHC_FREQUENCY*nBunches;

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
  sample T2tt_300_200 = sample("T2tt_300_200", branch, sampleDir + 
  			       "T2tt_300_200_Fall13/T2tt_300_200_Fall13*.root", 1);


  sample T1bbbb_2J_mGl_1000_mLSP_900 = sample("T1bbbb_2J_mGl_1000_mLSP_900", branch, sampleDir +
					      "T1bbbb_2J_mGl-1000_mLSP-900/T1bbbb_2J_mGl-1000_mLSP-900*.root", 1);

  sample T1tttt_2J_mGl_1200_mLSP_800 = sample("T1tttt_2J_mGl_1200_mLSP_800", branch, sampleDir +
					      "T1tttt_2J_mGl-1200_mLSP-800/T1tttt_2J_mGl-1200_mLSP-800*.root", 1);

  sample T2tt_2J_mStop_425_mLSP_325  = sample("T2tt_2J_mStop_425_mLSP_325", branch, sampleDir +
					     "T2tt_2J_mStop-425_mLSP-325/T2tt_2J_mStop-425_mLSP-325*.root", 1);

  sample T2tt_2J_mStop_500_mLSP_325  = sample("T2tt_2J_mStop_500_mLSP_325", branch, sampleDir +
					      "T2tt_2J_mStop-500_mLSP-325/T2tt_2J_mStop-500_mLSP-325*.root", 1);
  
  sample T2tt_2J_mStop_850_mLSP_100  = sample("T2tt_2J_mStop_850_mLSP_100", branch, sampleDir +
					      "T2tt_2J_mStop-850_mLSP-100/T2tt_2J_mStop-850_mLSP-100*.root", 1);


  // ------------------------------------------------------------------------------------------------------------------------
  sample selectedSample =  TTBar; //QCD800to1000; //QCD30to50; //T2tt_2J_mStop_850_mLSP_100; //QCD800to1000; //T2tt_500_250; //QCD30to50;  //test; //QCD30to50; // T2tt_500_250; //T2cc_250_210; //DYJets; //NuGun; //DYJets; //TTBar; //DYJets;

  // Label QCD ptHat bins
  int samplePTHat = 0;
  if      (selectedSample.Name == "QCD30to50")   { samplePTHat = 2;}
  else if (selectedSample.Name == "QCD50to80")   { samplePTHat = 3;}
  else if (selectedSample.Name == "QCD80to120")  { samplePTHat = 4;}
  else if (selectedSample.Name == "QCD120to170") { samplePTHat = 5;}
  else if (selectedSample.Name == "QCD170to300") { samplePTHat = 6;}
  else if (selectedSample.Name == "QCD300to470") { samplePTHat = 7;}
  else if (selectedSample.Name == "QCD470to600") { samplePTHat = 8;}
  else if (selectedSample.Name == "QCD600to800") { samplePTHat = 9;}
  else if (selectedSample.Name == "QCD800to1000"){ samplePTHat = 10;}


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

  unsigned int nuNEvents  = (unsigned int)rateChain ->GetEntries(); 
  unsigned int sigNEvents = (unsigned int)sigChain  ->GetEntries();

  // New scale factor calculation for QCD samples
  double rateScaleFactor = (sampleXS * instLumi)/nuNEvents;
  if (selectedSample.Name == "NuGun"){
    rateScaleFactor = double(ZB_XSECTION)/nuNEvents;
  }

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
  std::map< TString, TEfficiency*> hist2DOverlap;
  std::map< TString, TH1*>         histHLTRate;
  std::map< TString, TH2*>         histHLTRate2D;

  std::map< TString, TH1*>         histMatch;
  std::map< TString, TH2*>         histMatch2D;


  std::vector<float> *genJetPT  = new std::vector<float>();
  std::vector<float> *genJetPx  = new std::vector<float>();
  std::vector<float> *genJetPy  = new std::vector<float>();
  std::vector<float> *genJetEta = new std::vector<float>();
  std::vector<float> *genJetPhi = new std::vector<float>();
  std::vector<float> *genJetForPT  = new std::vector<float>();
  std::vector<float> *genJetForPx  = new std::vector<float>();
  std::vector<float> *genJetForPy  = new std::vector<float>();
  std::vector<float> *genJetForEta = new std::vector<float>();
  //  std::vector<float> *genJetForPhi = new std::vector<float>();

  float pfHT(0), pfMET(0), pfMHT(0); // pfAlphaT(0),
  std::vector<float> *caloJetPT  = new std::vector<float>();
  std::vector<float> *caloJetPx  = new std::vector<float>();
  std::vector<float> *caloJetPy  = new std::vector<float>();
  std::vector<float> *caloJetEta = new std::vector<float>();
  std::vector<float> *caloJetPhi = new std::vector<float>();
  float caloHT(0), caloMET(0), caloMHT(0); //, caloAlphaT(0);

  // std::vector<float> *hltPFJetPT  = new std::vector<float>();
  // std::vector<float> *hltPFJetPx  = new std::vector<float>();
  // std::vector<float> *hltPFJetPy  = new std::vector<float>();
  // std::vector<float> *hltPFJetEta = new std::vector<float>();
  // //  std::vector<float> *hltPFJetPhi = new std::vector<float>();
  // float hltPFHT(0), hltPFAlphaT(0);

  UInt_t NVTX(0);
  UChar_t genLeptonVeto(false);

  // Trigger bits
  std::vector<TString>       hltPathNames;
  std::map<TString, UChar_t> hltPathFired;


  float genMET;
  float hltCaloMET;


  // Store HLT paths
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
    
    // hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_NoL1_v1");
    // hltPathNames.push_back("HLT_RsqMR260_Rsq0p09_MR200_4jet_NoL1_v1");
    // hltPathNames.push_back("HLT_Rsq0p36_NoL1_v1");
    // hltPathNames.push_back("HLT_RsqMR300_Rsq0p09_MR200_v1");
    // hltPathNames.push_back("HLT_RsqMR260_Rsq0p09_MR200_4jet_v1");
    // hltPathNames.push_back("HLT_Rsq0p36_v1");


  // Trigger bits 
  for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
    TString path = hltPathNames[ iPath ];
    hltPathFired[ path ] = false;
  }


  // ------------------------------------------------------------
  // UCT
  // ------------------------------------------------------------
  std::vector<float> *uctJetPT  = new std::vector<float>();
  std::vector<float> *uctJetPhi = new std::vector<float>();
  float uctHT(0), uctMHT(0), uctMHToverHT(0), uctMET(0);


  // ------------------------------------------------------------
  // HLT
  // ------------------------------------------------------------
  // HLT CaloJet
  std::vector<float> *hltCaloJetPT  = new std::vector<float>();
  std::vector<float> *hltCaloJetPx  = new std::vector<float>();
  std::vector<float> *hltCaloJetPy  = new std::vector<float>();
  std::vector<float> *hltCaloJetEta = new std::vector<float>();
  //  std::vector<float> *hltCaloJetPhi = new std::vector<float>();
  float hltCaloHT(0); //, hltCaloAlphaT(0);
  // // HLT PF
  // std::vector<float> *hltPFJetPT  = new std::vector<float>();
  // std::vector<float> *hltPFJetPx  = new std::vector<float>();
  // std::vector<float> *hltPFJetPy  = new std::vector<float>();
  // std::vector<float> *hltPFJetEta = new std::vector<float>();
  //  std::vector<float> *hltPFJetPhi = new std::vector<float>();
  //  float hltPFJetHT(0), hltPFJetAlphaT(0);
  // // HLT PFNoPU
  // // std::vector<float> *hltPFNoPUPT  = new std::vector<float>();
  // // std::vector<float> *hltPFNoPUPx  = new std::vector<float>();
  // // std::vector<float> *hltPFNoPUPy  = new std::vector<float>();
  // // std::vector<float> *hltPFNoPUEta = new std::vector<float>();
  // // //  std::vector<float> *hltPFNoPUPhi = new std::vector<float>();
  // // float hltPFNoPUHT(0), hltPFNoPUAlphaT(0);

  // ------------------------------------------------------------
  // GEN
  // ------------------------------------------------------------
  // std::vector<float> *genJetPT  = new std::vector<float>();
  // std::vector<float> *genJetPx  = new std::vector<float>();
  // std::vector<float> *genJetPy  = new std::vector<float>();
  // std::vector<float> *genJetEta = new std::vector<float>();
  //  std::vector<float> *genJetPhi = new std::vector<float>();
  //  float genHT(0), genAlphaT(0), genMET(0);




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
  jet2PTCuts.push_back(0.);
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


  
  int   HTBins(40);
  float HTMin(0);
  float HTMax(1000);

  // int   MHTBins(50);
  // float MHTMin(0);
  // float MHTMax(250);

  int alphaTBins(60);
  float alphaTMin(0.4);
  float alphaTMax(1.0);

  // int   PTBins(50);
  // float PTMin(5);
  // float PTMax(305);

  // int   HT12Bins(100);
  // float HT12Min(5);
  // float HT12Max(505);


  
  // Trigger bits 
  for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
    TString path = hltPathNames[ iPath ];
    // histHLTRate[ path ]                  = new TH1D( path, path + ";Entries;Fired", 2, 0, 2.);
    // histHLTRate[ path  + "_AND_HTT175" ] = new TH1D( path + "_AND_HTT175", path + ";Entries;Fired", 2, 0, 2.);
    // histHLTRate[ path  + "_AND_Trig1" ]  = new TH1D( path + "_AND_Trig1", path + ";Entries;Fired", 2, 0, 2.);


    // Get pt hat breakdown of rate of the triggers
    makePTHatBinnedHist( histHLTRate, path );
    makePTHatBinnedHist( histHLTRate, path + "_AND_HTT175");
    makePTHatBinnedHist( histHLTRate, path + "_AND_Trig1");
    
    
    // Add alphaT with second jet threshold

    if ( path.Contains("L1HTT175OrETM70_v1") && !path.Contains("AlphaT0p5_") ){
      for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
      	TString newPath = path;
	
	float jet2PTCut = jet2PTCuts[ iJet2Cut ];
	TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
	newPath = newPath.ReplaceAll( "L1HTT175OrETM70_v1", jet2PTCutStr + TString("_L1HTT175OrETM70_v1") );

	makePTHatBinnedHist( histHLTRate, newPath );
      }
    }

  }
  


    // histHLTRate["HLT_PFHT900_v1"]                      = new TH1D("HLT_PFHT900_v1","HLT_PFHT900_v1;Entries;Fired", 2, 0, 2.); 
    // histHLTRate["HLT_PFHT350_PFMET120_NoiseCleaned_v1"]= new TH1D("HLT_PFHT350_PFMET120_NoiseCleaned_v1",
    // 								  "HLT_PFHT350_PFMET120_NoiseCleaned_v1;Entries;Fired", 2, 0, 2.); 
    // histHLTRate["HLT_HT200_AlphaT0p4_v1"]              = new TH1D("HLT_HT200_AlphaT0p4_v1",  
    // 								  "HLT_HT200_AlphaT0p4_v1;Entries;Fired", 2, 0, 2.); 
    // histHLTRate["HLT_HT200_PFAlphaT0p4_v1"]            = new TH1D("HLT_HT200_PFAlphaT0p4_v1",
    // 								  "HLT_HT200_PFAlphaT0p4_v1;Entries;Fired", 2, 0, 2.); 
    // histHLTRate["HLT_PFHT900_v1_AND_HTT200"]           = new TH1D("HLT_PFHT900_v1_AND_HTT200","HLT_PFHT900_v1 && HTT200;Entries;Fired", 2, 0, 2.); 


  hist1DRate["L1HTT"]       = new TH1D("L1HTT",       "L1HTT;H_{T} (GeV);Fired",         2000, 0, 2000);
  hist1DRate["L1MET"]       = new TH1D("L1MET",       "L1MET;#slash{E}_{T} (GeV);Fired", 2000, 0, 2000);
  hist1DRate["L1SingleJet"] = new TH1D("L1SingleJet", "L1 Single jet;p_{T} (GeV);Fired", 500,  0, 500);
  hist1DRate["L1DoubleJet"] = new TH1D("L1DoubleJet", "L1 Double jet;p_{T} (GeV);Fired", 500,  0, 500);
  hist1DRate["L1ThirdJet"]  = new TH1D("L1ThirdJet",  "L1 Third jet;p_{T} (GeV);Fired",  500,  0, 500);
  hist1DRate["L1QuadJet"]   = new TH1D("L1QuadJet",   "L1 Quad jet;p_{T} (GeV);Fired",   500,  0, 500);


    hist1DRate["NoL1_CaloHTRate_NoL1_Inclusive"]   = new TH1D("NoL1_CaloHTRate_NoL1_Inclusive"," Calojet H_{T} rate inclusive (NoL1);H_{T} (GeV)",    80,0,2000);
    hist1DRate["Trig1_CaloHTRate_Trig1_Inclusive"] = new TH1D("Trig1_CaloHTRate_Trig1_Inclusive"," Calojet H_{T} rate inclusive (Trig1) ;H_{T} (GeV)", 80,0,2000);
    hist1DRate["CaloHTRate_HTT175_Inclusive"] = new TH1D("CaloHTRate_HT175_Inclusive"," Calojet H_{T} rate inclusive (HTT175) ;H_{T} (GeV)", 80,0,2000);



    hist1DEff["HTT200_PFHT_TurnOn"] = new TEfficiency("HTT200_PFHT_TurnOn", "PF H_{T} turn on - HTT200;H_{T} (GeV);Efficiency", 2000, 0, 2000);

  // ********************************************************************************
  // Calojet efficiency for offline selection
  // ********************************************************************************
  
  // Vary the second jet PT cut
  for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
    float jet2PTCut = jet2PTCuts[ iJet2Cut ];
    TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
    TString jet2PTCutLab = "p_{T}^{j2} > " + TString(Form("%1.0f", jet2PTCut ));

    TString stdStr      = "On_AlphaTStd_vs_HT_"      + jet2PTCutStr;
    TString stdStrRate  = "OnRate_AlphaTStd_vs_HT_"  + jet2PTCutStr;
    TString dynStr      = "On_AlphaTDyn_vs_HT_"      + jet2PTCutStr;
    TString dynStrRate  = "OnRate_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
    // TString dyn2Str = "On_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
    // TString dyn3Str = "On_AlphaTDyn3_vs_HT_" + jet2PTCutStr;


    // Rate plots vs NVTX
    histHLTRate["LegacyAlphaT_HT200_" + jet2PTCutStr + "_vs_NVTX"] = new TH1D("LegacyAlphaT_HT200" + jet2PTCutStr + "_vs_NVTX",
									      "Fires Legacy #alpha_{T} HT200 vs NVTX;N_{VTX};Entries",
									      14, 0, 70.);
    histHLTRate["LegacyAlphaT_HT250_" + jet2PTCutStr + "_vs_NVTX"] = new TH1D("LegacyAlphaT_HT250" + jet2PTCutStr + "_vs_NVTX",
									      "Fires Legacy #alpha_{T} HT250 vs NVTX;N_{VTX};Entries",
									      14, 0, 70.);
    histHLTRate["LegacyAlphaT_HT300_" + jet2PTCutStr + "_vs_NVTX"] = new TH1D("LegacyAlphaT_HT300" + jet2PTCutStr + "_vs_NVTX",
									      "Fires Legacy #alpha_{T} HT300 vs NVTX;N_{VTX};Entries",
									      14, 0, 70.);
    histHLTRate["LegacyAlphaT_HT350_" + jet2PTCutStr + "_vs_NVTX"] = new TH1D("LegacyAlphaT_HT350" + jet2PTCutStr + "_vs_NVTX",
									      "Fires Legacy #alpha_{T} HT350 vs NVTX;N_{VTX};Entries",
									      14, 0, 70.);
    histHLTRate["LegacyAlphaT_HT400_" + jet2PTCutStr + "_vs_NVTX"] = new TH1D("LegacyAlphaT_HT400" + jet2PTCutStr + "_vs_NVTX",
									      "Fires Legacy #alpha_{T} HT400 vs NVTX;N_{VTX};Entries",
									      14, 0, 70.);
    


    // Rates
    // histHLTRate["LegacyAlphaT_"      + jet2PTCutStr] = new TH1D("LegacyAlphaT_" + jet2PTCutStr,
    // 								"Fires Legacy #alpha_{T};Entries;Fired", 2, 0, 2.);
    // histHLTRate["LegacyAlphaTExcSM_" + jet2PTCutStr] = new TH1D("LegacyAlphaTExcSM_" + jet2PTCutStr,
    // 								"Fires Legacy #alpha_{T} - Excluding SUSY menu;Entries;Fired", 2, 0, 2.);


    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr] = new TH2D("SM_vs_LegacyAlphaT_" + jet2PTCutStr,
								   "SUSY hadronic menu overlap with #alpha_{T} triggers", 7, 0, 7, 2, 0, 2);
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 1, "Not fired");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 2, "Fired");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 3, "H_{T} = 200");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 4, "H_{T} = 250");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 5, "H_{T} = 300");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 6, "H_{T} = 350");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 7, "H_{T} = 400");

    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetYaxis()->SetBinLabel( 1, "Not fired");
    histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->GetYaxis()->SetBinLabel( 2, "Fired");


    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr] = new TH2D("LegacyAlphaT_vs_SM_" + jet2PTCutStr,
								   "#alpha_{T} menu overlap with SUSY hadronic triggers", 6, 0, 6, 2, 0, 2);
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 1, "Not fired");
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 2, "Fired");
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 3, "PFMET170");
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 4, "PFMET120_Btag");
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 5, "PFHT350_PFMet120");
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetXaxis()->SetBinLabel( 6, "PFHT900");

    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetYaxis()->SetBinLabel( 1, "Not fired");
    histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->GetYaxis()->SetBinLabel( 2, "Fired");



    histHLTRate2D["HLT_AlphaT_vs_HT_" + jet2PTCutStr] = new TH2D("HLT_AlphaT_vs_HT_" + jet2PTCutStr,"HLT #alpha_{T}^{static} vs H_{T}^{static} differential rate (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


    histHLTRate2D["HLT_AlphaT_vs_HT_" + jet2PTCutStr + "_" + selectedSample.Name] = new TH2D("HLT_AlphaT_vs_HT_" + jet2PTCutStr + "_" + selectedSample.Name,"HLT #alpha_{T}^{static} vs H_{T}^{static} differential rate (" + jet2PTCutLab + ") " + selectedSample.Name + ";H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    



    // ------------------------------------------------------------------------------------------------------------------------
    // Inclusive
    // ------------------------------------------------------------------------------------------------------------------------
    hist2DOnEff[stdStr  + "_Inclusive"] = new TEfficiency(stdStr  + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DOnEff[dynStr  + "_Inclusive"] = new TEfficiency(dynStr  + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DOnEff[dyn2Str + "_Inclusive"] = new TEfficiency(dyn2Str + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DOnEff[dyn3Str + "_Inclusive"] = new TEfficiency(dyn3Str + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


    // Make candidate L1 trigger plots
    for (uint trigN = 1; trigN <= 1; ++trigN ){
      TString trigStr = TString("Trig") + Form("%d", trigN);

      hist2DRate[trigStr + stdStrRate  + "_Inclusive"] = new TH2D(trigStr + stdStrRate  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DRate[trigStr + dynStrRate  + "_Inclusive"] = new TH2D(trigStr + dynStr  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DRate[trigStr + dyn2Str + "_Inclusive"] = new TH2D(trigStr + dyn2Str + "_Inclusive",trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DRate[trigStr + dyn3Str + "_Inclusive"] = new TH2D(trigStr + dyn3Str + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);



      hist2DRate[trigStr + "_CaloMHT_vs_CaloHT_" + jet2PTCutStr]  = new TH2D(trigStr + "_CaloMHT_vs_CaloHT_" + jet2PTCutStr,trigStr + " Calojet #slash{H}_{T} vs vs H_{T} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#slash{H}_{T} (GeV)",   HTBins,HTMin,HTMax,  HTBins,HTMin,HTMax);
      hist2DRate[trigStr + "_CaloMET_vs_CaloHT_" + jet2PTCutStr]  = new TH2D(trigStr + "_CaloMET_vs_CaloHT_" + jet2PTCutStr,trigStr + " Calojet #slash{E}_{T} vs vs H_{T} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#slash{E}_{T} (GeV)",   HTBins,HTMin,HTMax,  HTBins,HTMin,HTMax);



      // Excluding overlap rate with SUSY hadronic menu
      if (trigN == 1){
	hist2DRate[trigStr + "SM" + stdStrRate  + "_Inclusive"] = new TH2D(trigStr + "SM" + stdStrRate  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ") - Excluding SUSY hadronic menu overlap;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      }

    }
    hist2DRate["NoL1" + stdStrRate  + "_Inclusive"] = new TH2D("NoL1" + stdStrRate  + "_Inclusive","No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DRate["NoL1" + dynStrRate  + "_Inclusive"] = new TH2D("NoL1" + dynStrRate  + "_Inclusive","No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DRate["NoL1" + dyn2Str + "_Inclusive"] = new TH2D("NoL1" + dyn2Str + "_Inclusive","No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DRate["NoL1" + dyn3Str + "_Inclusive"] = new TH2D("NoL1" + dyn3Str + "_Inclusive","No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // ------------------------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------------

    
    // ------------------------------------------------------------------------------------------------------------------------
    // NJet binned
    // ------------------------------------------------------------------------------------------------------------------------
    for (int jetMult = 2; jetMult <= MAX_JETS; ++jetMult){
      TString jetBinStr = TString(Form( "%d", jetMult)) + TString("Jets");
      
      hist2DOnEff[stdStr  + "_" + jetBinStr] = new TEfficiency(stdStr  + "_" + jetBinStr,"Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " 
							       + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   
							       HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff[dynStr  + "_" + jetBinStr] = new TEfficiency(dynStr  + "_" + jetBinStr,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency "
      // 							       + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", 
      // 							       HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff[dyn2Str + "_" + jetBinStr] = new TEfficiency(dyn2Str + "_" + jetBinStr,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " 
      // 							       + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", 
      // 							       HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff[dyn3Str + "_" + jetBinStr] = new TEfficiency(dyn3Str + "_" + jetBinStr,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " 
      // 							       + jetBinStr + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  
      // 							       HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


      // Make CaloHLT -> PFHLT trigger plots 
      for (uint trigN = 1; trigN <= 5; ++trigN ){
	TString trigStr = TString("HLT") + Form("%d", trigN);

	hist2DHLTEff[trigStr + stdStr  + "_" + jetBinStr] = new TEfficiency(trigStr + stdStr  + "_" + jetBinStr,"HLT Calo #alpha_{T}^{static} vs H_{T}^{static} efficiency given HLT PF " + trigStr + " " + jetBinStr + " (" + jet2PTCutLab + ");HLT Calo H_{T} (GeV);HLT Calo #alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	hist2DHLTEff[trigStr + dynStr  + "_" + jetBinStr] = new TEfficiency(trigStr + dynStr  + "_" + jetBinStr,"HLT Calo #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency given HLT PF " + trigStr + " " + jetBinStr + " (" + jet2PTCutLab + ");HLT Calo H_{T} (GeV);HLT Calo #alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  
      }

    }
    // ------------------------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------------

    // ------------------------------------------------------------------------------------------------------------------------      
    // Analysis binned
    // ------------------------------------------------------------------------------------------------------------------------
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
	  // hist2DOnEff[dyn2Str + suffix] = new TEfficiency(dyn2Str + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  // hist2DOnEff[dyn3Str + suffix] = new TEfficiency(dyn3Str + suffix,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  

	  
    TString effOverlap = jet2PTCutStr + suffix;
    hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap] = new TEfficiency("SM_vs_LegacyAlphaT_" + effOverlap,
									"SUSY hadronic menu overlap with #alpha_{T} triggers " + label + 
									" (" + jet2PTCutLab + ")", 7, 0, 7, 2, 0, 2);
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 1, "Not fired");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 2, "Fired");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 3, "H_{T} = 200");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 4, "H_{T} = 250");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 5, "H_{T} = 300");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 6, "H_{T} = 350");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetXaxis()->SetBinLabel( 7, "H_{T} = 400");

    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetYaxis()->SetBinLabel( 1, "Not fired");
    // hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->GetYaxis()->SetBinLabel( 2, "Fired");


    hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap] = new TEfficiency("LegacyAlphaT_vs_SM_" + effOverlap,
									"#alpha_{T} menu overlap with SUSY hadronic triggers " + label + 
									" (" + jet2PTCutLab + ")", 6, 0, 6, 2, 0, 2);
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetXaxis()->SetBinLabel( 1, "Not fired");
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetXaxis()->SetBinLabel( 2, "Fired");
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetXaxis()->SetBinLabel( 3, "PFMET170");
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetXaxis()->SetBinLabel( 4, "PFMET120_Btag");
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetXaxis()->SetBinLabel( 5, "PFHT350_PFMet120");
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetXaxis()->SetBinLabel( 6, "PFHT900");

    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetYaxis()->SetBinLabel( 1, "Not fired");
    // hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->GetYaxis()->SetBinLabel( 2, "Fired");



    // Make candidate L1 trigger plots
    for (uint trigN = 1; trigN <= 1; ++trigN ){
      TString trigStr = TString("Trig") + Form("%d", trigN);

	hist2DOnEff[trigStr + stdStr  + suffix] = new TEfficiency(trigStr + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	if (trigN == 1){
	  hist2DOnEff[trigStr + "SM" + stdStr  + suffix] = new TEfficiency(trigStr + "SM" + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ") - Including SUSY hadronic menu overlap;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  hist2DOnEff[trigStr + "SMOnly" + stdStr  + suffix] = new TEfficiency(trigStr + "SM" + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ") - SUSY hadronic menu only;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	}

	hist2DOnEff[trigStr + dynStr  + suffix] = new TEfficiency(trigStr + dynStr  + suffix,trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	// hist2DOnEff[trigStr + dyn2Str + suffix] = new TEfficiency(trigStr + dyn2Str + suffix,trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	// hist2DOnEff[trigStr + dyn3Str + suffix] = new TEfficiency(trigStr + dyn3Str + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

    }
      
      hist2DOnEff["NoL1" + stdStr  + suffix] = new TEfficiency("NoL1" + stdStr  + suffix,"No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      hist2DOnEff["NoL1" + dynStr  + suffix] = new TEfficiency("NoL1" + dynStr  + suffix,"No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff["NoL1" + dyn2Str + suffix] = new TEfficiency("NoL1" + dyn2Str + suffix,"No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff["NoL1" + dyn3Str + suffix] = new TEfficiency("NoL1" + dyn3Str + suffix,"No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


      }
    } // End analysis bin
    
    // ------------------------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------------

  } // End second jet cut







  // ************************************************************************************************************************
  // *                                                    Matching study                                                    *
  // ************************************************************************************************************************
#ifdef MATCHING  

  std::vector<TString> jetTypes;
  jetTypes.push_back("Inclusive");
  jetTypes.push_back("Matched");
  jetTypes.push_back("Pileup");
  jetTypes.push_back("dRMatched");

  std::vector<TString> jetRanks;
  jetRanks.push_back("AllJets");
  jetRanks.push_back("Jet1");
  jetRanks.push_back("Jet2");

  // ************************************************************
  // Jet ranked distributions
  // ************************************************************
  for (uint iRank = 0; iRank < jetRanks.size(); ++iRank){
    TString rank = jetRanks[ iRank ];

    
    // ************************************************************
    // Jet type distributions
    // ************************************************************
    
      for (uint iType = 0; iType < jetTypes.size(); ++iType){
	TString jet = jetTypes[ iType ];
	

      histMatch[jet + "_" + rank + "_JetPT"]       = new TH1D(jet + "_" + rank + "_JetPT",  
							      jet + " " + rank + " - Jet P_{T};P_{T} (GeV);Entries", 
							      56, 22.5, 302.5);
      histMatch[jet + "_" + rank + "_JetEta"]      = new TH1D(jet + "_" + rank + "_JetEta",
							      jet + " " + rank + " - Jet #eta;#eta;Entries",      
							      31, -3.1, 3.1);
      
      // Gen rank matching
      histMatch[jet + "_" + rank + "_JetResponse"] = new TH1D(jet + "_" + rank + "_JetResponse", 
							      jet + " " + rank + " - rank matched Jet response;P_{T}^{HLT}/P_{T}^{GEN};Entries", 
							      200,-0.01, 3.99);
      histMatch[jet + "_" + rank + "_DeltaPTRel"]  = new TH1D(jet + "_" + rank + "_DeltaPTRel",  
							      jet + " " + rank + " - rank matched #DeltaP_{T}^{Rel};#DeltaP_{T}^{Rel};Entries", 
							      30,-1.0,3.0);
      histMatch[jet + "_" + rank + "_DeltaEta"]    = new TH1D(jet + "_" + rank + "_DeltaEta", 
							      jet + " " + rank + " - rank matched #Delta#eta;Entries", 
							      121,-0.605,0.605);
      histMatch[jet + "_" + rank + "_DeltaPhi"]    = new TH1D(jet + "_" + rank + "_DeltaPhi",  
							      jet + " " + rank + " - rank matched #Delta#phi;Entries", 
							      121,-0.605,0.605);
      histMatch[jet + "_" + rank + "_DeltaR"]      = new TH1D(jet + "_" + rank + "_DeltaR",   
							      jet + " " + rank + " - rank matched #DeltaR;Entries", 
							      100,-0.005,0.995);
      
      
      histMatch2D[jet + "_" + rank + "_HLTPT_vs_GENPT"]       = new TH2D(jet + "_" + rank + "_HLTPT_vs_GENPT", 
									 jet + " " + rank + " - rank matched HLT P_{T} vs Gen P_{T};Gen p_{T} (GeV);HLT P_{T} (GeV)", 
									 56,22.5,302.5,  56,22.5,302.5);
      histMatch2D[jet + "_" + rank + "_DeltaPTRel_vs_GENPT"]  = new TH2D(jet + "_" + rank + "_DeltaPTRel_vs_GENPT", 
									 jet + " " + rank + " - rank matched #DeltaP_{T}^{Rel} vs Gen P_{T};Gen p_{T} (GeV);#DeltaP_{T}^{Rel}", 
									 56,22.5,302.5, 30,-1.0,3.0);
      histMatch2D[jet + "_" + rank + "_DeltaPTRel_vs_NVTX"]   = new TH2D(jet + "_" + rank + "_Matched_DeltaPTRel_vs_NVTX",
									 jet + " " + rank + " - rank matched #DeltaP_{T}^{Rel} vs Gen P_{T};Gen p_{T} (GeV);#DeltaP_{T}^{Rel}", 
									 16,0,80, 30,-1.0,3.0);
   
      
      // Inter-jet HLT matching 
      //  ADD : deltaEta, Pt, phi, jet(1,2) in bins of pileup jets in lead two jets vs GEN of the same quantity


    } // End jet rank
    // ************************************************************

  } // End jet types



  // ********************************************************************************
  // Resolution of analysis quantities
  // ********************************************************************************

  histMatch2D["NPileup_vs_NVTX"]       = new TH2D("NPileup_vs_NVTX",
						  "Number of pileup analysis jets vs N_{VTX};N_{VTX};N_{Pileup}",
						  60,20.5,80.5, 7,-0.5,6.5);
  
  
  histMatch2D["HLTHT_vs_GENHT"]         = new TH2D("HLTHT_vs_GENHT", 
						   "HLT H_{T} vs Gen H_{T};Gen H_{T} (GeV);HLT H_{T} (GeV)",
						  40,0,1000,  40,0,1000);
  histMatch2D["HLTMHT_vs_GENMHT"]       = new TH2D("HLTMHT_vs_GENMHT", 
						   "HLT #slash{H}_{T} vs Gen #slash{H}_{T};Gen #slash{H}_{T} (GeV);HLT #slash{H}_{T} (GeV)", 
						   40,0,400,  40,0,400); 
  histMatch2D["HLTAlphaT_vs_GENAlphaT"] = new TH2D("HLTAlphaT_vs_GENAlphaT", 
						   "HLT #alpha_{T} vs Gen #alpha_{T};Gen #alpha_{T};HLT #alpha_{T}", 
						   230,-1.3,1.0,  230,-1.3,1.0); 
  

  histMatch2D["DeltaHT_vs_HT"]         = new TH2D("DeltaHT_vs_HT",
						  "#DeltaH_{T} vs H_{T};H_{T};(H_{T}^{HLT} - H_{T}^{Gen})/H_{T}^{Gen};",
						  40,0,1000, 40,-2.0,2.0);


  // ****************************************
  // Resolutions
  // ****************************************
  histMatch2D["DeltaHT_vs_GENHT"]     = new TH2D("DeltaHT_vs_GENHT", 
						 "#DeltaH_{T} vs Gen H_{T};Gen H_{T} (GeV);(HLT H_{T} - Gen H_{T})/Gen H_{T}",
						 40,0,1000, 33,-1.3,3.0);
  histMatch2D["DeltaMHT_vs_GENMHT"]    = new TH2D("DeltaMHT_vs_GENMHT", 
						 "#Delta#slash{H}_{T} vs Gen #slash{H}_{T};Gen #slash{H}_{T} (GeV);(HLT #slash{H}_{T} - Gen #slash{H}_{T})/Gen #slash{H}_{T}",
						 40,0,400, 33,-1.3,3.0);
  histMatch2D["DeltaAlphaT_vs_GENAlphaT"] = new TH2D("DeltaAlphaT_vs_GENAlphaT",
						     "#Delta#alpha_{T} vs Gen #alpha_{T};Gen #alpha_{T};(#alpha_{T}^{HLT} - #alpha_{T}^{Gen})/#alpha_{T}^{Gen};",
						     8,0.2,1.0, 230,-1.3,1.0);


  
  histMatch2D["DeltaHT_vs_NPileup"]   = new TH2D("DeltaHT_vs_NPileup", 
						 "#DeltaH_{T} vs Number of pileup analysis jets;N_{Pileup};(HLT H_{T} - Gen H_{T})/Gen H_{T}",
						 7,-0.5,6.5, 33,-1.3,3.0);
  histMatch2D["DeltaMHT_vs_NPileup"]  = new TH2D("DeltaMHT_vs_NPileup", 
						 "#Delta#slash{H}_{T} vs Number of pileup analysis jets;N_{Pileup};(HLT #slash{H}_{T} - Gen #slash{H}_{T})/Gen #slash{H}_{T}",
						 7,-0.5,6.5, 33,-1.3,3.0);
  histMatch2D["DeltaAlphaT_vs_NPileup"] = new TH2D("DeltaAlphaT_vs_NPileup",
						   "#Delta#alpha_{T} vs Number of pileup analysis jets;N_{Pileup};(#alpha_{T}^{HLT} - #alpha_{T}^{Gen})/#alpha_{T}^{Gen}",
						   7,-0.5,6.5, 230,-1.3,1.0);
  

  // // Performance distributions binned by NVTX
  // histMatch2D["DeltaAlphaT_vs_NVTX"] = new TH2D("DeltaAlphaT_vs_NVTX",
  // 						"#Delta#alpha_{T} vs N_{VTX};(#alpha_{T}^{HLT} - #alpha_{T}^{Gen})/#alpha_{T}^{Gen};N_{VTX};",
  // 						12,20,80, 40,-2.0,2.0);
  // histMatch2D["DeltaHT_vs_NVTX"]     = new TH2D("DeltaHT_vs_NVTX",
  // 						"#DeltaH_{T} vs N_{VTX};(H_{T}^{HLT} - H_{T}^{Gen})/H_{T}^{Gen};N_{VTX};",
  // 						12,20,80, 40,-2.0,2.0);


  






  // ************************************************************
  // Turn ons
  // ************************************************************

  // ****************************************
  // HT
  // ****************************************
  std::vector<float> HTValues;
  HTValues.push_back( 200 ); HTValues.push_back( 250 );   
  HTValues.push_back( 300 ); HTValues.push_back( 350 ); HTValues.push_back( 400 );

  for (uint iHT = 0; iHT < HTValues.size(); ++iHT){
    float HTValue = HTValues[iHT];
    TString HTStr = Form("%1.0f", HTValue );

    hist1DEff["HLTCaloHT_HLTPFHT" + HTStr + "_TurnOn"] = new TEfficiency("HLTCalo_HLTPFHT" + HTStr + "_TurnOn",
									 "HLTCalo H_{T} turn on (HLTPF H_{T} > " + HTStr + ");Efficiency", 
									 40,0,800);
    hist1DEff["HLTCaloHT_GenHT" + HTStr + "_TurnOn"]   = new TEfficiency("HLTCalo_GenHT" + HTStr + "_TurnOn",
									 "HLTCalo H_{T} turn on (Gen H_{T} > " + HTStr + ");Efficiency",  
									 40,0,800);
    hist1DEff["HLTPFHT_GenHT" + HTStr + "_TurnOn"]     = new TEfficiency("HLTPF_GenHT" + HTStr + "_TurnOn",
									 "HLTPF H_{T} turn on (Gen H_{T} > " + HTStr + ");Efficiency",  
									 40,0,800);
  }

  // ****************************************
  // Jet PT
  // ****************************************
  std::vector<float> PTValues;
  PTValues.push_back( 40 ); PTValues.push_back( 50 ); PTValues.push_back( 60 ); PTValues.push_back( 70 );
  PTValues.push_back( 80 ); PTValues.push_back( 90 ); PTValues.push_back( 100 );

  for (uint iPT = 0; iPT < PTValues.size(); ++iPT){
    float PTValue = PTValues[iPT];
    TString PTStr = Form("%1.0f", PTValue );

    // hist1DEff["HLTCaloPT_HLTPFPT" + HTStr + "_TurnOn"] = new TEfficiency("HLTCalo_HLTPFHT" + HTStr + "_TurnOn",
    // 									 "HLTCalo H_{T} turn on (HLTPF H_{T} > " + HTStr + ");Efficiency", 
    // 									 40,0,800);
    // hist1DEff["HLTCaloPT_GenPT" + HTStr + "_TurnOn"]   = new TEfficiency("HLTCalo_GenHT" + HTStr + "_TurnOn",
    // 									 "HLTCalo H_{T} turn on (Gen H_{T} > " + HTStr + ");Efficiency",  
    // 									 40,0,800);
    // hist1DEff["HLTPFPT_GenPT" + HTStr + "_TurnOn"]     = new TEfficiency("HLTPF_GenHT" + HTStr + "_TurnOn",
    // 									 "HLTPF H_{T} turn on (Gen H_{T} > " + HTStr + ");Efficiency",  
    // 									 40,0,800);
  }

  
  // ************************************************************
  // Correlations
  // ************************************************************
  histMatch2D["HLTCaloHT_vs_HLTPFHT"]  = new TH2D("HLTCaloHT_vs_HLTPFHT", 
						  "HLTCalo H_{T} vs HLTPF H_{T};H_{T}^{HLTPF} (GeV);H_{T}^{HLTCalo} (GeV)",    
						  40,0,1000, 40,0,1000 );
  histMatch2D["HLTCaloHT_vs_GenHT"]    = new TH2D("HLTCaloHT_vs_GenHT",   
						  "HLTCalo H_{T} vs Gen H_{T};H_{T}^{Gen} (GeV);H_{T}^{HLTCalo} (GeV)",   
						  40,0,1000, 40,0,1000 );
  histMatch2D["HLTPFHT_vs_GenHT"]      = new TH2D("HLTPFHT_vs_GenHT",     
						  "HLTPF H_{T} vs Gen H_{T};H_{T}^{Gen} (GeV);H_{T}^{HLTPF} (GeV)",    
						  40,0,1000, 40,0,1000 );


  // AlphaT Correlations 
  histMatch2D["HLTCaloAlphaT_vs_HLTPFAlphaT"]      = new TH2D("HLTCaloAlphaT_vs_HLTPFAlphaT",  
							 "HLTCalo #alpha_{T} vs HLTPF #alpha_{T};#alpha_{T}^{HLTPF};#alpha_{T}^{HLTCalo}",
							 alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax); 
  histMatch2D["HLTCaloAlphaT_vs_GenAlphaT"]        = new TH2D("HLTCaloAlphaT_vs_GenAlphaT",
							 "HLTPF #alpha_{T} vs Gen #alpha_{T};#alpha_{T}^{Gen};#alpha_{T}^{HLTCalo}",    
							 alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
  histMatch2D["HLTPFAlphaT_vs_GenAlphaT"]          = new TH2D("HLTPFAlphaT_vs_GenAlphaT",
							 "HLTPF #alpha_{T} vs Gen #alpha_{T};#alpha_{T}^{Gen};#alpha_{T}^{HLTPF}",    
							 alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);

  histMatch2D["HLTCaloAlphaT_vs_HLTCaloAlphaTDyn"] = new TH2D("CaloAlphaT_vs_CaloAlphaTDyn",   
							 "HLTCalo #alpha_{T} vs dynamic HLTCalo #alpha_{T};#alpha_{T}^{Dynamic HLTCalo};#alpha_{T}^{HLTCalo}",    
							 alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);
  histMatch2D["HLTCaloAlphaTDyn_vs_HLTPFAlphaT"]   = new TH2D("HLTCaloAlphaTDyn_vs_HLTPFAlphaT",  
							 "Dynamic HLTCalo #alpha_{T} vs HLTPF #alpha_{T};#alpha_{T}^{HLTPF};#alpha_{T}^{Dynamic HLTCalo}",
							 alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax); 
  histMatch2D["HLTCaloAlphaTDyn_vs_GenAlphaT"]     = new TH2D("HLTCaloAlphaTDyn_vs_GenAlphaT",
							 "Dynamic HLTPF #alpha_{T} vs Gen #alpha_{T};#alpha_{T}^{Gen};#alpha_{T}^{Dynamic HLTCalo}",    
							 alphaTBins,alphaTMin,alphaTMax, alphaTBins,alphaTMin,alphaTMax);

  histMatch2D["GenAlphaT_vs_GenHT"]                = new TH2D("GenAlphaT_vs_GenHT",  
							 "Gen #alpha_{T} vs Gen H_{T};Gen H_{T} (GeV);#alpha_{T}^{Gen}",
							 40,0,1000, alphaTBins,alphaTMin,alphaTMax); 



#endif





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
    // // HLT CaloJet
    // sigChain->SetBranchAddress("hltAk4Calo_Pt",        &hltCaloJetPT);
    // sigChain->SetBranchAddress("hltAk4Calo_Px",        &hltCaloJetPx);
    // sigChain->SetBranchAddress("hltAk4Calo_Py",        &hltCaloJetPy);
    // //sigChain->SetBranchAddress("hltAk4Calo_Phi",       &hltCaloJetPhi);
    // sigChain->SetBranchAddress("hltAk4Calo_Eta",       &hltCaloJetEta);
  
    // // HLT PF 
    // sigChain->SetBranchAddress("hltAk4PF_Pt",          &hltPFJetPT); 
    // sigChain->SetBranchAddress("hltAk4PF_Px",          &hltPFJetPx);
    // sigChain->SetBranchAddress("hltAk4PF_Py",          &hltPFJetPy);
    // //sigChain->SetBranchAddress("hltAk4PF_Phi",         &hltPFJetPhi);
    // sigChain->SetBranchAddress("hltAk4PF_Eta",         &hltPFJetEta);
   
    
    // // ------------------------------------------------------------ 
    // // GEN
    // // ------------------------------------------------------------ 
    // sigChain->SetBranchAddress("genAk4_Pt",         &genJetPT); 
    // sigChain->SetBranchAddress("genAk4_Px",         &genJetPx);
    // sigChain->SetBranchAddress("genAk4_Py",         &genJetPy);
    // //sigChain->SetBranchAddress("genAk4_Phi",      &genJetPhi);
    // sigChain->SetBranchAddress("genAk4_Eta",        &genJetEta);
    // sigChain->SetBranchAddress("genMetTrue_MetPt", &genMET);


    // TEMPORARY - REPLACE NAMES FOR MEETING IN 1.5 hours
    sigChain->SetBranchAddress("genAk4_Pt",         &genJetPT);
    sigChain->SetBranchAddress("genAk4_Px",         &genJetPx);
    sigChain->SetBranchAddress("genAk4_Py",         &genJetPy);
    //sigChain->SetBranchAddress("genAk4_Phi",      &genJetPhi);
    sigChain->SetBranchAddress("genAk4_Eta",        &genJetEta);

    sigChain->SetBranchAddress("genAk4For_Pt",         &genJetForPT);
    sigChain->SetBranchAddress("genAk4For_Px",         &genJetForPx);
    sigChain->SetBranchAddress("genAk4For_Py",         &genJetForPy);
    //sigChain->SetBranchAddress("genAk4For_Phi",      &genJetForPhi);
    sigChain->SetBranchAddress("genAk4For_Eta",        &genJetForEta);

    sigChain->SetBranchAddress("genMetCalo_MetPt",     &pfMET);


    sigChain->SetBranchAddress("hltAk4Calo_Pt",        &hltCaloJetPT);
    sigChain->SetBranchAddress("hltAk4Calo_Px",        &hltCaloJetPx);
    sigChain->SetBranchAddress("hltAk4Calo_Py",        &hltCaloJetPy);
    //sigChain->SetBranchAddress("hltAk4Calo_Phi",     &hltCaloJetPhi);
    sigChain->SetBranchAddress("hltAk4Calo_Eta",       &hltCaloJetEta);

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


    // Trigger bits
    for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
      TString path = hltPathNames[ iPath ];
      hltPathFired[ path ] = false;
      sigChain->SetBranchAddress( path, &hltPathFired[ path ]);
    }




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

  std::cout << "Processing signal events:\n"
	    << "\tLow   = " << signalEventLow  
	    << "\tHigh  = " << signalEventHigh 
	    << "\tTotal = " << sigNEvents << "\n\n";

  std::cout << "Processing neutrino events:\n"
	    << "\tLow   = " << nuEventLow  
	    << "\tHigh  = " << nuEventHigh 
	    << "\tTotal = " << nuNEvents << "\n\n";

  std::cout << "Sample xs (pb)         = " << selectedSample.CrossSection     << "\n";
  std::cout << "Sample xs (cms2)       = " << sampleXS                        << "\n";
  std::cout << "Rate scalefactor       = " << rateScaleFactor                 << "\n";
  std::cout << "Rate scalefactor error = " << rateScaleFactor*sqrt(nuNEvents) << "\n\n";


#ifdef SIGNAL

  // Loop over the tree
  for ( unsigned int iEvent = signalEventLow; iEvent < signalEventHigh; ++iEvent ){
    
    sigChain->GetEntry( iEvent );

    float pfAlphaTStandard(0);
    float caloAlphaTStandard(0);
    float caloAlphaTDynamic(0);
    float hltCaloAlphaTStandard(0);
    float hltCaloAlphaTDynamic(0);

    // float caloHTDynamic(0);
    //    std::pair<float,float> caloAlphaTHTDynamic;
    

    // ********************************************************************************
    // *                                 UCT triggers                                 *
    // ********************************************************************************

    bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) ); // DJ100
    // bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 60) ); // DJ120
    // bool l1Trig3 = ( (uctHT >= 160)  || (uctMET >= 70) ); // DJ120
    // bool l1Trig4 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 150) && (uctMHToverHT >= 0.17)) );
    // bool l1Trig5 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 140) && (uctMHToverHT >= 0.25)) );
    // bool l1Trig6 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 130) && (uctMHToverHT >= 0.32)) );
    // bool l1Trig7 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 120) && (uctMHToverHT >= 0.36)) );


    // HLTCalo
    int hltCaloJetsAboveThresh(0);
    hltCaloHT = 0;
    for(uint iJet = 0;iJet < hltCaloJetPT->size(); ++iJet ){
      if( (*hltCaloJetPT)[iJet] < caloJetThreshold ) break;
      hltCaloHT += (*hltCaloJetPT)[iJet];
      hltCaloJetsAboveThresh++;
    }
    if ( hltCaloJetsAboveThresh > MAX_JETS ){
      hltCaloJetsAboveThresh = MAX_JETS;
    }

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
    int genJetsAboveThresh(0);
    float pfMHTX(0), pfMHTY(0);
    pfHT  = 0;
    pfMHT = 0;
    for(uint iJet = 0;iJet < genJetPT->size(); ++iJet ){
      if( (*genJetPT)[iJet] < genJetThreshold ) break;
      pfHT += (*genJetPT)[iJet];
      pfMHTX += (*genJetPx)[iJet];
      pfMHTY += (*genJetPy)[iJet];
      genJetsAboveThresh++;
    }
    pfMHT = sqrt( pfMHTX*pfMHTX + pfMHTY*pfMHTY);
    if ( genJetsAboveThresh > MAX_JETS ){
      genJetsAboveThresh = MAX_JETS;
    }
    TString genJetBinStr     = TString(Form( "%d", genJetsAboveThresh)) + TString("Jets");


    // Calo HLT cuts
    // bool passesATStandard(false), passesATDynamic(false), passesATDynamic2(false), passesATDynamic3(false);
    // bool passesHTStandard(false), passesHTDynamic(false);
    //    bool passesJet2_50(false), passesJet2_60(false), passesJet2_70(false), passesJet2_80(false), passesJet2_90(false), passesJet2_100(false);
    //    bool passesDynJet2_100(false);    

    //    bool passesOffHT(false), passesOffAT(false), passesOffJet(false), passesOffAll(false);
    bool passesOffJet(false);
    bool passesAnaBinOffAT(false), passesAnaBinOffJet(false), passesAnaBinOffAll(false);
    bool passesOffVetoes(false);


    pfAlphaTStandard       = calculateAlphaT( genJetPT, genJetPx, genJetPy, genJetThreshold );      

    caloAlphaTStandard     = calculateAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold );    

    hltCaloAlphaTStandard  = calculateAlphaT( hltCaloJetPT, hltCaloJetPx, hltCaloJetPy, caloJetThreshold );    


    // Dynamic AlphaT
    hltCaloAlphaTDynamic   = calculateDynamicAlphaT( hltCaloJetPT, hltCaloJetPx, hltCaloJetPy, 
						     maxCaloJet, caloJetDynThreshold, 
						     caloJetAlphaThreshold );
  


    // ********************************************************************************
    // *                              Calojets 
    // ********************************************************************************
    if ( genJetsAboveThresh >= 2 ){
      
      // Forward jet veto
      bool forJetVeto(false);
      if (genJetForPT->size() > 0){
	if ( (*genJetForPT)[0] > genJetThreshold){
	  forJetVeto = true;
	}
      }

      // MHT over MET cleaning
      bool mhtOverMetVeto(false);
      float MHToverMET(0);
      if (pfMET > 0){
	MHToverMET = pfMHT/pfMET;
	if (MHToverMET > 1.25){ mhtOverMetVeto = true; }
      }

      // Check individual offline cuts
      if ( (genJetsAboveThresh >= 2 ) && ((*genJetPT)[1] > genJet2PTThreshold) )          { passesOffJet    = true; }
      if ( !(genLeptonVeto) && !(forJetVeto) && !(mhtOverMetVeto) )                    { passesOffVetoes = true; }
      // Check full offline selection
      //      if ( passesOffHT && passesOffAT && passesOffJet && passesOffVetoes )     { passesOffAll    = true; }
      // Vary the online second jet PT cut
      for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){

	float jet2PTCut = jet2PTCuts[ iJet2Cut ];
	TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
	TString stdStr  = "On_AlphaTStd_vs_HT_"  + jet2PTCutStr;
	TString dynStr  = "On_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
	// TString dyn2Str = "On_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
	// TString dyn3Str = "On_AlphaTDyn3_vs_HT_" + jet2PTCutStr;

	bool passHltSecondJet(false);
	bool passHltPreFilterSecondJet(false);

	// Second jet threshold
	if ((*caloJetPT)[1] > jet2PTCut){
	    passHltSecondJet = true;
	}
	if (hltCaloJetsAboveThresh >= 2){	
	  if ((*hltCaloJetPT)[1] > jet2PTCut){
	    passHltPreFilterSecondJet = true;
	  }
	}

	// Efficiency to PFHLT given CaloHLT prefilter
	// ********************************************************************************
	// *                                 HLT triggers                                 *
	// ********************************************************************************

	// Check event passes L1 and PFHLT triggers
	bool HLT1(false), HLT2(false), HLT3(false), HLT4(false), HLT5(false);
	if (l1Trig1 && ( (*caloJetPT)[1] > hltPFSecondJetThreshold) ){ // PFHLT triggers
	    HLT1 = ( (caloHT > 200)  && (caloAlphaTStandard > 0.57) );
	    HLT2 = ( (caloHT > 250)  && (caloAlphaTStandard > 0.55) );
	    HLT3 = ( (caloHT > 300)  && (caloAlphaTStandard > 0.53) );
	    HLT4 = ( (caloHT > 350)  && (caloAlphaTStandard > 0.52) );
	    HLT5 = ( (caloHT > 400)  && (caloAlphaTStandard > 0.51) );
	}

      
	// ********************************************************************************
	// *                                   Prefilter                                  *
	// ********************************************************************************
	if ( (*caloJetPT)[1] > hltPFSecondJetThreshold ){
	  TString trigStr  = "";
	  bool    trigBool = false;
	  
	  if ( passHltPreFilterSecondJet ){

	    trigStr  = "HLT1";
	    trigBool = HLT1;
	    hist2DHLTEff[trigStr + stdStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTStandard );
	    hist2DHLTEff[trigStr + dynStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTDynamic );
	    trigStr  = "HLT2";
	    trigBool = HLT2;
	    hist2DHLTEff[trigStr + stdStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTStandard );
	    hist2DHLTEff[trigStr + dynStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTDynamic );
	    trigStr  = "HLT3";
	    trigBool = HLT3;
	    hist2DHLTEff[trigStr + stdStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTStandard );
	    hist2DHLTEff[trigStr + dynStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTDynamic );
	    trigStr  = "HLT4";
	    trigBool = HLT4;
	    hist2DHLTEff[trigStr + stdStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTStandard );
	    hist2DHLTEff[trigStr + dynStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTDynamic );
	    trigStr  = "HLT5";
	    trigBool = HLT5;
	    hist2DHLTEff[trigStr + stdStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTStandard );
	    hist2DHLTEff[trigStr + dynStr  + "_" + genJetBinStr]->Fill( trigBool, hltCaloHT, hltCaloAlphaTDynamic );
  
	  }
	}
      

	// Get offline analysis bin
	for (uint iHtBin = 0; iHtBin < anaHtBins.size(); ++iHtBin){
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
	  if ( (genJetsAboveThresh >= 2 ) && ((*genJetPT)[1] > genJet2PTThreshold) ){ passesAnaBinOffJet = true; }
	  // Check full offline selection
	  if ( passesAnaBinOffAT && passesAnaBinOffJet )                         { passesAnaBinOffAll = true; }

	  // Get jet bin
	  for (uint iJetBin = 0; iJetBin < anaJetBins.size(); ++iJetBin){
	
	    int jetLow  = anaJetBins[ iJetBin ].first;
	    int jetHigh = anaJetBins[ iJetBin ].second;
	    if ( !((genJetsAboveThresh >= jetLow) && (genJetsAboveThresh < jetHigh)) ){ continue; }
	    TString anaJetStr = Form("%d", jetLow) + TString("to") + Form("%d", jetHigh);
	    
	    TString suffix = TString("_") + anaHtStr + TString("_") + anaJetStr;
	    
	    hist2DOnEff[stdStr  + suffix] ->Fill( passesAnaBinOffAll, caloHT,                     caloAlphaTStandard );
	    // hist2DOnEff[dynStr  + suffix] ->Fill( passesAnaBinOffAll, caloHT,                     caloAlphaTDynamic );
	    // hist2DOnEff[dyn2Str + suffix] ->Fill( passesAnaBinOffAll, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	    // hist2DOnEff[dyn3Str + suffix] ->Fill( passesAnaBinOffAll, caloAlphaTHTDynamic.second, caloAlphaTStandard );
	    

	    // Make candidate L1 trigger plots
    	    if ( passesAnaBinOffAll){

	      TString trigStr  = "";
	      bool    trigBool = false;

	      trigStr  = "NoL1";
	      trigBool = passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      //	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig1";
	      trigBool = l1Trig1 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      //	      hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      // ******************************************************************************** 
	      // Emulate AlphaT HLTrigger 
	      // ******************************************************************************** 
	      bool legacyAlphaTHT200(false), legacyAlphaTHT250(false), legacyAlphaTHT300(false), legacyAlphaTHT350(false), legacyAlphaTHT400(false); 
	      if (trigBool){
		if ( (caloHT > 200) && (caloAlphaTStandard > 0.57) ){ legacyAlphaTHT200 = true; }
		if ( (caloHT > 250) && (caloAlphaTStandard > 0.55) ){ legacyAlphaTHT250 = true; }
		if ( (caloHT > 300) && (caloAlphaTStandard > 0.53) ){ legacyAlphaTHT300 = true; }
		if ( (caloHT > 350) && (caloAlphaTStandard > 0.52) ){ legacyAlphaTHT350 = true; }
		if ( (caloHT > 400) && (caloAlphaTStandard > 0.51) ){ legacyAlphaTHT400 = true; }
	      }
 	      bool firesLegacyAlphaT = ( legacyAlphaTHT200 || legacyAlphaTHT250 ||
					 legacyAlphaTHT300 || legacyAlphaTHT350 ||
					 legacyAlphaTHT400 );
	      // Not including BTAG
	      bool firesSUSYMenu = ((hltPathFired["HLT_PFMET170_NoiseCleaned_v1"])           || 
				    //(hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_v1"]) ||
				    (hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_v1"] )  ||
				    (hltPathFired["HLT_PFHT900_v1"]) );
	      // Overlap efficiency with SUSY hadronic menu
	      hist2DOnEff[trigStr + "SM" + stdStr + suffix]->Fill( (firesSUSYMenu||trigBool), caloHT, caloAlphaTStandard );
	      // Efficiency of SUSY hadronic menu alone
	      hist2DOnEff[trigStr + "SMOnly" + stdStr + suffix]->Fill( firesSUSYMenu, caloHT, caloAlphaTStandard );
	      
	      
	      // ********************************************************************************
	      // *                            Overlaps between menues                           *
	      // ********************************************************************************
	      TString effOverlap = jet2PTCutStr + suffix;



      hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap] ->Fill( !firesLegacyAlphaT&&!firesSUSYMenu, 0.5,   0.5);
      hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap] ->Fill(  firesLegacyAlphaT&&!firesSUSYMenu, 1.5, 0.5);
      hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap] ->Fill( !firesLegacyAlphaT&&firesSUSYMenu,  0.5,   1.5);
      hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap] ->Fill(  firesLegacyAlphaT&&firesSUSYMenu,  1.5, 1.5);

      hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap] ->Fill( !firesLegacyAlphaT&&!firesSUSYMenu, 0.5,   0.5);
      hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap] ->Fill(  firesLegacyAlphaT&&!firesSUSYMenu, 0.5,   1.5);
      hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap] ->Fill( !firesLegacyAlphaT&& firesSUSYMenu, 1.5, 0.5);
      hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap] ->Fill(  firesLegacyAlphaT&&firesSUSYMenu,  1.5, 1.5);

	  // HT200
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT200 && !firesSUSYMenu, 2.5, 0.5);
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT200 &&  firesSUSYMenu, 2.5, 1.5);
	  // HT250
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT250 && !firesSUSYMenu, 3.5, 0.5);
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT250 &&  firesSUSYMenu, 3.5, 1.5);
	  // HT300
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT300 && !firesSUSYMenu, 4.5, 0.5);
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT300 &&  firesSUSYMenu, 4.5, 1.5);
	  // HT350
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT350 && !firesSUSYMenu, 5.5, 0.5);
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT350 &&  firesSUSYMenu, 5.5, 1.5);
	  // HT400
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT400 && !firesSUSYMenu, 6.5, 0.5);
	  hist2DOverlap["SM_vs_LegacyAlphaT_" + effOverlap]->Fill( legacyAlphaTHT400 &&  firesSUSYMenu, 6.5, 1.5);

	
	  // HLT_PFMET170_NoiseCleaned
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFMET170_NoiseCleaned_v1"] && !firesLegacyAlphaT, 2.5, 0.5);
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFMET170_NoiseCleaned_v1"] &&  firesLegacyAlphaT, 2.5, 1.5);
	  // HLT_PFMET120_NoiseCleaned_BTagCSV07
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_v1"] && !firesLegacyAlphaT, 3.5, 0.5);
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_v1"] &&  firesLegacyAlphaT, 3.5, 1.5);
	  // HLT_PFHT350_PFMET120_NoiseCleaned
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_v1"] && !firesLegacyAlphaT, 4.5, 0.5);
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_v1"] &&  firesLegacyAlphaT, 4.5, 1.5);
	  // HLT_PFHT900
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFHT900_v1"] && !firesLegacyAlphaT, 5.5, 0.5);
	  hist2DOverlap["LegacyAlphaT_vs_SM_" + effOverlap]->Fill( hltPathFired["HLT_PFHT900_v1"] &&  firesLegacyAlphaT, 5.5, 1.5);

	



	      
	    }// End signal region requirement
	    
	    
	  } // End analysis jet bin
	} // End analysis HT bin

      } // End HLT second jet cut
    } // End offline second jet cut



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


    // // ------------------------------------------------------------ 
    // // HLT 
    // // ------------------------------------------------------------ 
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
   
    
    // // ------------------------------------------------------------ 
    // // GEN
    // // ------------------------------------------------------------ 
    // rateChain->SetBranchAddress("genAk4_Pt",         &genJetPT); 
    // rateChain->SetBranchAddress("genAk4_Px",         &genJetPx);
    // rateChain->SetBranchAddress("genAk4_Py",         &genJetPy);
    // //rateChain->SetBranchAddress("genAk4_Phi",      &genJetPhi);
    // rateChain->SetBranchAddress("genAk4_Eta",        &genJetEta);


    // TEMPORARY - REPLACE NAMES FOR MEETING IN 1.5 hours
    rateChain->SetBranchAddress("genAk4_Pt",         &genJetPT); 
    rateChain->SetBranchAddress("genAk4_Px",         &genJetPx);
    rateChain->SetBranchAddress("genAk4_Py",         &genJetPy);
    rateChain->SetBranchAddress("genAk4_Phi",      &genJetPhi);
    rateChain->SetBranchAddress("genAk4_Eta",        &genJetEta);

#ifdef HLT_CALOJET
    // HLT CaloJet
    rateChain->SetBranchAddress("hltAk4Calo_Pt",        &caloJetPT);
    rateChain->SetBranchAddress("hltAk4Calo_Px",        &caloJetPx);
    rateChain->SetBranchAddress("hltAk4Calo_Py",        &caloJetPy);
    rateChain->SetBranchAddress("hltAk4Calo_Phi",       &caloJetPhi);
    rateChain->SetBranchAddress("hltAk4Calo_Eta",       &caloJetEta);


    // rateChain->SetBranchAddress("hltAk4PF_Pt",        &hltPFJetPT);
    // rateChain->SetBranchAddress("hltAk4PF_Px",        &hltPFJetPx);
    // rateChain->SetBranchAddress("hltAk4PF_Py",        &hltPFJetPy);
    // //rateChain->SetBranchAddress("hltAk4PF_Phi",       &hltPFJetPhi);
    // rateChain->SetBranchAddress("hltAk4PF_Eta",       &hltPFJetEta);

#else
    // HLT PFJet
    rateChain->SetBranchAddress("hltAk4PF_Pt",        &caloJetPT);
    rateChain->SetBranchAddress("hltAk4PF_Px",        &caloJetPx);
    rateChain->SetBranchAddress("hltAk4PF_Py",        &caloJetPy);
    rateChain->SetBranchAddress("hltAk4PF_Phi",       &caloJetPhi);
    rateChain->SetBranchAddress("hltAk4PF_Eta",       &caloJetEta);
#endif


    rateChain->SetBranchAddress("genMetCaloAndNonPrompt_MetPt", &genMET);
    rateChain->SetBranchAddress("hltMetCalo_MetPT",             &hltCaloMET);


    // Trigger bits
    for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
      TString path = hltPathNames[ iPath ];
      hltPathFired[ path ] = false;
      rateChain->SetBranchAddress( path, &hltPathFired[ path ]);
    }

    rateChain->SetBranchAddress("NVTX", &NVTX);






#ifdef NEUTRINO

    // Loop over the tree
    for ( unsigned int iEvent = nuEventLow; iEvent < nuEventHigh; ++iEvent ){



    rateChain->GetEntry( iEvent );

    // ****************************************************************************************************
    // *                                        Check L1 triggers                                         *
    // ****************************************************************************************************
    bool l1Trig1 = ( (uctHT >= 175)  || (uctMET >= 70) ); // DJ100
    // bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 60) ); // DJ120
    // bool l1Trig3 = ( (uctHT >= 160)  || (uctMET >= 70) ); // DJ120
    // bool l1Trig4 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 150) && (uctMHToverHT >= 0.17)) );
    // bool l1Trig5 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 140) && (uctMHToverHT >= 0.25)) );
    // bool l1Trig6 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 130) && (uctMHToverHT >= 0.32)) );
    // bool l1Trig7 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 120) && (uctMHToverHT >= 0.36)) );
    
    bool l1HTT175 = (uctHT  >= 175);
    bool l1ETM70  = (uctMET >= 70);
    // ****************************************************************************************************
    // *                                         Check HLTriggers                                         *
    // ****************************************************************************************************

    bool firesSUSYMenu                          = ((hltPathFired["HLT_PFMET170_NoiseCleaned_v1"])           ||
						   (hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_v1"]) ||
						   (hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_v1"] )  ||
						   (hltPathFired["HLT_PFHT900_v1"]) );


    // ************************************************************
    // Trigger bits 
    // ************************************************************
    for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
      TString path = hltPathNames[ iPath ];
      //histHLTRate[path] ->Fill( hltPathFired[ path ] );
      histHLTRate[path + "_PTHat"]->Fill( hltPathFired[ path ] );
      histHLTRate[path + "_PTHat"]->Fill( samplePTHat, hltPathFired[ path ] );


      if ( l1HTT175 ){
	//histHLTRate[ path + "_AND_HTT175" ]     ->Fill( hltPathFired[ path ] );
	histHLTRate[ path + "_AND_HTT175_PTHat"]->Fill( hltPathFired[ path ] );
	histHLTRate[ path + "_AND_HTT175_PTHat"]->Fill( samplePTHat, hltPathFired[ path ] );
      }
      if ( l1Trig1 ){
	//	histHLTRate[ path + "_AND_Trig1"]      ->Fill( hltPathFired[ path ] );
	histHLTRate[ path + "_AND_Trig1_PTHat"]->Fill( hltPathFired[ path ] );
        histHLTRate[ path + "_AND_Trig1_PTHat"]->Fill( samplePTHat, hltPathFired[ path ] );
      }
    }

    hist1DRate["L1HTT"]     ->Fill( uctHT );
    hist1DRate["L1MET"]     ->Fill( uctMET );
    if ( uctJetPT->size() >= 1){
      hist1DRate["L1SingleJet"]->Fill( (*uctJetPT)[0] );
      if ( uctJetPT->size() > 1){
        hist1DRate["L1DoubleJet"]->Fill( (*uctJetPT)[1] );
	if ( uctJetPT->size() > 2){
	  hist1DRate["L1ThirdJet"]->Fill( (*uctJetPT)[2] );
	  if ( uctJetPT->size() > 3){
	    hist1DRate["L1QuadJet"]  ->Fill( (*uctJetPT)[3] );
	  }
	}
      }
    }


    // ********************************************************************************
    // *                              Calojets 
    // ********************************************************************************
    float caloAlphaTStandard(0), caloAlphaTDynamic(0);
    caloHT        = 0;
    caloMHT       = 0;
    caloAlphaTStandard  = calculateAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold );
    // Dynamic AlphaT
    caloAlphaTDynamic   = calculateDynamicAlphaT( caloJetPT, caloJetPx, caloJetPy, //caloJetThreshold, 
						  maxCaloJet, caloJetDynThreshold, 
    						  caloJetAlphaThreshold );
    

    // Calojets
    int caloJetsAboveThresh(0);
    float caloMHTx(0), caloMHTy(0);
    for(uint iJet = 0;iJet < caloJetPT->size(); ++iJet ){
      if( (*caloJetPT)[iJet] < caloJetThreshold ) break;
      caloHT    += (*caloJetPT)[iJet];
      caloMHTx      += (*caloJetPx)[iJet];
      caloMHTy      += (*caloJetPy)[iJet];

      caloJetsAboveThresh++;
    }
    caloMHT = sqrt( caloMHTx*caloMHTx + caloMHTy*caloMHTy);
    if ( caloJetsAboveThresh > MAX_JETS ){
      caloJetsAboveThresh = MAX_JETS;
    }

    int genJetsAboveThresh(0);
    float pfMHTX(0), pfMHTY(0);
    pfHT = 0;
    for(uint iJet = 0;iJet < genJetPT->size(); ++iJet ){
      if( (*genJetPT)[iJet] < genJetThreshold ) break;
      pfHT   += (*genJetPT)[iJet];
      pfMHTX += (*genJetPx)[iJet];
      pfMHTY += (*genJetPy)[iJet];
      genJetsAboveThresh++;
    }
    pfMHT = sqrt( pfMHTX*pfMHTX + pfMHTY*pfMHTY);






    
    // Calo HT rate
    // --------------------------------------------------------------------------------
    hist1DRate["NoL1_CaloHTRate_NoL1_Inclusive"]    ->Fill( caloHT );
    if (l1Trig1){
      hist1DRate["Trig1_CaloHTRate_Trig1_Inclusive"]->Fill( caloHT ); 
    }
    if ( l1HTT175 ){
      hist1DRate["CaloHTRate_HTT175_Inclusive"]->Fill( caloHT );
    }

    hist1DEff["HTT200_PFHT_TurnOn"] ->Fill( (uctHT >= 200), caloHT );



    // ****************************************************************************************************
    // Vary the online second jet PT cut
    // ****************************************************************************************************
    for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){
      bool legacyAlphaTHT200(false), legacyAlphaTHT250(false), legacyAlphaTHT300(false), legacyAlphaTHT350(false), legacyAlphaTHT400(false);
      bool firesLegacyAlphaT(false);

      float jet2PTCut = jet2PTCuts[ iJet2Cut ];
      TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
	  
      if ( caloJetsAboveThresh >= 2 ){
	if ((*caloJetPT)[1] > jet2PTCut){


	  TString stdStr  = "OnRate_AlphaTStd_vs_HT_"  + jet2PTCutStr;
	  TString dynStr  = "OnRate_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
	  // TString dyn2Str = "OnRate_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
	  // TString dyn3Str = "OnRate_AlphaTDyn3_vs_HT_" + jet2PTCutStr;
	  

	  // ********************************************************************************
	  // Emulate AlphaT HLTrigger
	  // ********************************************************************************
	  if (l1Trig1){
	    if ( (caloHT > 200) && (caloAlphaTStandard > 0.57) ){ legacyAlphaTHT200 = true; }
	    if ( (caloHT > 250) && (caloAlphaTStandard > 0.55) ){ legacyAlphaTHT250 = true; }
	    if ( (caloHT > 300) && (caloAlphaTStandard > 0.53) ){ legacyAlphaTHT300 = true; }
	    if ( (caloHT > 350) && (caloAlphaTStandard > 0.52) ){ legacyAlphaTHT350 = true; }
	    if ( (caloHT > 400) && (caloAlphaTStandard > 0.51) ){ legacyAlphaTHT400 = true; }
	  }
	  firesLegacyAlphaT = ( legacyAlphaTHT200 || legacyAlphaTHT250 || legacyAlphaTHT300 || legacyAlphaTHT350 || legacyAlphaTHT400 );
	  

	  // Differential rate distribution 
	  histHLTRate2D["HLT_AlphaT_vs_HT_" + jet2PTCutStr]                            ->Fill(caloHT, caloAlphaTStandard );
	  histHLTRate2D["HLT_AlphaT_vs_HT_" + jet2PTCutStr + "_" + selectedSample.Name]->Fill(caloHT, caloAlphaTStandard );

	  // HLT trigger bits
	  TString hltJet2 = jet2PTCutStr + TString("_L1HTT175OrETM70_v1_PTHat");
	  histHLTRate[ "HLT_HT200_AlphaT0p57_" + hltJet2]->Fill( hltPathFired["HLT_HT200_AlphaT0p57_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT200_AlphaT0p57_" + hltJet2]->Fill( samplePTHat, hltPathFired["HLT_HT200_AlphaT0p57_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT250_AlphaT0p55_" + hltJet2]->Fill( hltPathFired["HLT_HT250_AlphaT0p55_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT250_AlphaT0p55_" + hltJet2]->Fill( samplePTHat, hltPathFired["HLT_HT250_AlphaT0p55_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT300_AlphaT0p53_" + hltJet2]->Fill( hltPathFired["HLT_HT300_AlphaT0p53_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT300_AlphaT0p53_" + hltJet2]->Fill( samplePTHat, hltPathFired["HLT_HT300_AlphaT0p53_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT350_AlphaT0p52_" + hltJet2]->Fill( hltPathFired["HLT_HT350_AlphaT0p52_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT350_AlphaT0p52_" + hltJet2]->Fill( samplePTHat, hltPathFired["HLT_HT350_AlphaT0p52_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT400_AlphaT0p51_" + hltJet2]->Fill( hltPathFired["HLT_HT400_AlphaT0p51_L1HTT175OrETM70_v1"] );
	  histHLTRate[ "HLT_HT400_AlphaT0p51_" + hltJet2]->Fill( samplePTHat, hltPathFired["HLT_HT400_AlphaT0p51_L1HTT175OrETM70_v1"] );

	  
	  // Make candidate L1 trigger plots
	  TString trigStr = "";
	  trigStr = "NoL1";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );

	  if ( l1Trig1){
	    trigStr = "Trig1";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]           ->Fill( caloHT,                     caloAlphaTStandard );
	    hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );

	    
	    hist2DRate[trigStr + "_CaloMHT_vs_CaloHT_" + jet2PTCutStr]->Fill( caloHT, caloMHT );
	    hist2DRate[trigStr + "_CaloMET_vs_CaloHT_" + jet2PTCutStr]->Fill( caloHT, hltCaloMET );

	    // Exclude overlap rate with SUSY hadronic menu - assuming SUSY menu has same HTT175||MET70 trigger
	    if ( !( firesSUSYMenu ) ){
	      hist2DRate[trigStr + "SM" + stdStr + "_Inclusive"]->Fill( caloHT,                     caloAlphaTStandard );
	    }

	  }
	  // if ( l1Trig2){
	  //   trigStr = "Trig2";
	  //   hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  //   // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  //   // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  //   // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  // }

	} // End Second jet PT requirement

      } // End two HLT jet requirement



      // ****************************************************************************************************
      // Fill trigger overlaps
      // ****************************************************************************************************
      histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr] ->Fill( firesLegacyAlphaT, firesSUSYMenu);
      histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr] ->Fill( firesSUSYMenu,     firesLegacyAlphaT);
	  
      // Breakdown of triggers fired
      if ( firesLegacyAlphaT ){
	    
	if (legacyAlphaTHT200){
	  histHLTRate["LegacyAlphaT_HT200_" + jet2PTCutStr + "_vs_NVTX"]->Fill( NVTX );
	  histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->Fill( 2.5, firesSUSYMenu);
	}
	if (legacyAlphaTHT250){
	  histHLTRate["LegacyAlphaT_HT250_" + jet2PTCutStr + "_vs_NVTX"]->Fill( NVTX );
	  histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->Fill( 3.5, firesSUSYMenu);
	}
	if (legacyAlphaTHT300){
	  histHLTRate["LegacyAlphaT_HT300_" + jet2PTCutStr + "_vs_NVTX"]->Fill( NVTX );
	  histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->Fill( 4.5, firesSUSYMenu);
	}
	if (legacyAlphaTHT350){
	  histHLTRate["LegacyAlphaT_HT350_" + jet2PTCutStr + "_vs_NVTX"]->Fill( NVTX );
	  histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->Fill( 5.5, firesSUSYMenu);
	}
	if (legacyAlphaTHT400){
	  histHLTRate["LegacyAlphaT_HT400_" + jet2PTCutStr + "_vs_NVTX"]->Fill( NVTX );
	  histHLTRate2D["SM_vs_LegacyAlphaT_" + jet2PTCutStr]->Fill( 6.5, firesSUSYMenu);
	}

      }
      if ( firesSUSYMenu ){
	
	if (hltPathFired["HLT_PFMET170_NoiseCleaned_v1"]){
	  histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->Fill( 2.5, firesLegacyAlphaT);
	}
	if (hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_v1"]){
	  histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->Fill( 3.5, firesLegacyAlphaT);
	}
	if (hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_v1"] ){
	  histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->Fill( 4.5, firesLegacyAlphaT);
	}
	if (hltPathFired["HLT_PFHT900_v1"]){
	  histHLTRate2D["LegacyAlphaT_vs_SM_" + jet2PTCutStr]->Fill( 5.5, firesLegacyAlphaT);
	}
	
      }


    } // 2nd Jet cut












    // ****************************************************************************************************
    // *                                           Jet matching                                           *
    // ****************************************************************************************************

#ifdef MATCHING  

    // std::vector<TLorentzVector> genVec, hltVec;

    // for ( uint iJet = 0; iJet < genJetPT->size(); ++iJet ){
    //   float jetPt  = (*genJetPT)[iJet];
    //   float jetEta = (*genJetEta)[iJet];
    //   float jetPhi = (*genJetPhi)[iJet];
      
    //   TLorentzVector genJet;
    //   genJet.SetPtEtaPhiM( jetPt, jetEta, jetPhi, 0 );
    //   genVec.push_back( genJet );
    // }
    // for ( uint iJet = 0; iJet < caloJetPT->size(); ++iJet ){
    //   float jetPt  = (*caloJetPT)[iJet];
    //   float jetEta = (*caloJetEta)[iJet];
    //   float jetPhi = (*caloJetPhi)[iJet];
      
    //   TLorentzVector hltJet;
    //   hltJet.SetPtEtaPhiM( jetPt, jetEta, jetPhi, 0 );
    //   hltVec.push_back( hltJet );
    // }

    // Perform matching
    std::vector<pair_info> pairs = make_pairs( genJetPT, genJetEta, genJetPhi, caloJetPT, caloJetEta, caloJetPhi, 20., 20. );
    std::sort(pairs.begin(), pairs.end(), sortDR);
    std::vector<int> HLTMatchedIndex = analyse_pairs(pairs, caloJetPT->size(), maxDeltaR);

    // Number of pileup jets that pass analysis jet threshold
    int nPileup(0);

    // Count number of analysis pileup jets
    for(unsigned int iHLT = 0; iHLT < HLTMatchedIndex.size(); ++iHLT ) {
      float hltPt  = (*caloJetPT)[iHLT];
      int iGEN = HLTMatchedIndex[iHLT];
      if ( iGEN == -1 ){ 
	if (hltPt > caloJetThreshold ){
	  nPileup++; // Count analysis jets attributed to pileup
	}
      }
    }
    

    // Iterate through HLT jets
    for(unsigned int iHLT = 0; iHLT < HLTMatchedIndex.size(); ++iHLT ) {

      float hltPt  = (*caloJetPT)[iHLT];
      float hltEta = (*caloJetEta)[iHLT];
      float hltPhi = (*caloJetPhi)[iHLT];
      // dR matched
      float dRGenPt  = 0;
      float dRGenEta = 0;
      float dRGenPhi = 0;
      // rank matched
      float rGenPt  = 0;
      float rGenEta = 0;
      float rGenPhi = 0;
      // Pileup or matched to gen
      TString jetID = "";
      // Jet leading/sub-leading rank - Jet1 or Jet2
      TString rank = "";
      
      // Find corresponding GEN jet match if it exists
      int iGEN = HLTMatchedIndex[iHLT];
      bool matched(true);
      if ( iGEN == -1 ){ matched = false; }



      histMatch["Inclusive_AllJets_JetPT"] ->Fill( hltPt );
      histMatch["Inclusive_AllJets_JetEta"]->Fill( hltEta );
      if (iHLT < 2){
	rank = TString("Jet") + Form("%d", (iHLT + 1) );

	histMatch["Inclusive_" + rank + "_JetPT"] ->Fill( hltPt );
	histMatch["Inclusive_" + rank + "_JetEta"]->Fill( hltEta );
      }
      if (matched){
	jetID = "Matched";
	dRGenPt  = (*genJetPT)[iGEN];
	dRGenEta = (*genJetEta)[iGEN];
	dRGenPhi = (*genJetPhi)[iGEN];

	histMatch["Matched_AllJets_JetPT"] ->Fill( hltPt );
	histMatch["Matched_AllJets_JetEta"]->Fill( hltEta );
	if (iHLT < 2){
	  histMatch["Matched_" + rank + "_JetPT"] ->Fill( hltPt );
	  histMatch["Matched_" + rank + "_JetEta"]->Fill( hltEta );
	}
      }
      else{
	// Make assumption non-matched jets are attributed to pileup      
	jetID = "Pileup";

	histMatch["Pileup_AllJets_JetPT"] ->Fill( hltPt );
	histMatch["Pileup_AllJets_JetEta"]->Fill( hltEta );
	if (iHLT < 2){
	  histMatch["Pileup_" + rank + "_JetPT"] ->Fill( hltPt );
	  histMatch["Pileup_" + rank + "_JetEta"]->Fill( hltEta );
	}
      }

      
      // ********************************************************************************
      // Rank matching
      // ********************************************************************************
      if (iHLT < genJetPT->size()){

	rGenPt  = (*genJetPT)[iHLT];
	rGenEta = (*genJetEta)[iHLT];
	rGenPhi = (*genJetPhi)[iHLT];

	TLorentzVector rGenJet;
	rGenJet.SetPtEtaPhiM( rGenPt, rGenEta, rGenPhi, 0 );
	TLorentzVector hltJet;
	hltJet.SetPtEtaPhiM( hltPt, hltEta, hltPhi, 0 );
	
	double rDeltaPt     = (hltPt - rGenPt);
	double rDeltaPtRel  = (hltPt - rGenPt)/rGenPt;
	double rJetResponse = hltPt/rGenPt; 
	double rDeltaEta    = (rGenEta - hltEta); 
	double rDeltaPhi    = rGenJet.DeltaPhi( hltJet );
	double rDeltaR      = sqrt( rDeltaEta*rDeltaEta + rDeltaPhi*rDeltaPhi );

	
	// Iterate through jet types
	for (uint iJet = 0; iJet < 2; ++iJet){
	  TString jet = "";
	  if (iJet == 0){ jet = "Inclusive"; }
	  else          { jet = jetID; }
	  
	  histMatch[  jet + "_AllJets_JetResponse"]        ->Fill(rJetResponse);
	  histMatch[  jet + "_AllJets_DeltaPTRel"]         ->Fill(rDeltaPtRel);
	  histMatch[  jet + "_AllJets_DeltaEta"]           ->Fill(rDeltaEta);
	  histMatch[  jet + "_AllJets_DeltaPhi"]           ->Fill(rDeltaPhi);
	  histMatch[  jet + "_AllJets_DeltaR"]             ->Fill(rDeltaR);
	  histMatch2D[jet + "_AllJets_HLTPT_vs_GENPT"]     ->Fill(rGenPt, hltPt);
	  histMatch2D[jet + "_AllJets_DeltaPTRel_vs_GENPT"]->Fill(rGenPt, rDeltaPtRel);
	  histMatch2D[jet + "_AllJets_DeltaPTRel_vs_NVTX"] ->Fill(NVTX,   rDeltaPtRel);

	  if (iHLT < 2){
	    histMatch[  jet + "_" + rank + "_JetResponse"]        ->Fill(rJetResponse);
	    histMatch[  jet + "_" + rank + "_DeltaPTRel"]         ->Fill(rDeltaPtRel);
	    histMatch[  jet + "_" + rank + "_DeltaEta"]           ->Fill(rDeltaEta);
	    histMatch[  jet + "_" + rank + "_DeltaPhi"]           ->Fill(rDeltaPhi);
	    histMatch[  jet + "_" + rank + "_DeltaR"]             ->Fill(rDeltaR);
	    histMatch2D[jet + "_" + rank + "_HLTPT_vs_GENPT"]     ->Fill(rGenPt, hltPt);
	    histMatch2D[jet + "_" + rank + "_DeltaPTRel_vs_GENPT"]->Fill(rGenPt, rDeltaPtRel);
	    histMatch2D[jet + "_" + rank + "_DeltaPTRel_vs_NVTX"] ->Fill(NVTX,   rDeltaPtRel);
	  }
	
	}
      
	// ********************************************************************************
	// dR matching
	// ********************************************************************************
	if (matched){
	  // Make dR distribution, correlations

	  TLorentzVector dRGenJet;
	  dRGenJet.SetPtEtaPhiM( dRGenPt, dRGenEta, dRGenPhi, 0 );
	
	  double dRDeltaPt     = (hltPt - dRGenPt);
	  double dRDeltaPtRel  = (hltPt - dRGenPt)/dRGenPt;
	  double dRJetResponse = hltPt/dRGenPt; 
	  double dRDeltaEta    = dRGenEta - hltEta; 
	  double dRDeltaPhi    = dRGenJet.DeltaPhi( hltJet );
	  double dRDeltaR      = sqrt( dRDeltaEta*dRDeltaEta + dRDeltaPhi*dRDeltaPhi );
	  
	
	  histMatch[  "dRMatched_AllJets_JetResponse"]        ->Fill(rJetResponse);
          histMatch[  "dRMatched_AllJets_DeltaPTRel"]         ->Fill(rDeltaPtRel);
          histMatch[  "dRMatched_AllJets_DeltaEta"]           ->Fill(rDeltaEta);
          histMatch[  "dRMatched_AllJets_DeltaPhi"]           ->Fill(rDeltaPhi);
          histMatch[  "dRMatched_AllJets_DeltaR"]             ->Fill(rDeltaR);
          histMatch2D["dRMatched_AllJets_HLTPT_vs_GENPT"]     ->Fill(rGenPt, hltPt);
          histMatch2D["dRMatched_AllJets_DeltaPTRel_vs_GENPT"]->Fill(rGenPt, rDeltaPtRel);
          histMatch2D["dRMatched_AllJets_DeltaPTRel_vs_NVTX"] ->Fill(NVTX,   rDeltaPtRel);

          if (iHLT < 2){
            histMatch[  "dRMatched_" + rank + "_JetResponse"]        ->Fill(rJetResponse);
            histMatch[  "dRMatched_" + rank + "_DeltaPTRel"]         ->Fill(rDeltaPtRel);
            histMatch[  "dRMatched_" + rank + "_DeltaEta"]           ->Fill(rDeltaEta);
            histMatch[  "dRMatched_" + rank + "_DeltaPhi"]           ->Fill(rDeltaPhi);
            histMatch[  "dRMatched_" + rank + "_DeltaR"]             ->Fill(rDeltaR);
            histMatch2D["dRMatched_" + rank + "_HLTPT_vs_GENPT"]     ->Fill(rGenPt, hltPt);
            histMatch2D["dRMatched_" + rank + "_DeltaPTRel_vs_GENPT"]->Fill(rGenPt, rDeltaPtRel);
            histMatch2D["dRMatched_" + rank + "_DeltaPTRel_vs_NVTX"] ->Fill(NVTX,   rDeltaPtRel);
          }
	} // End dRMatching



	
	// ********************************************************************************
	// Energy sums matching
	// ********************************************************************************
	
	// By default make events with no offline HT/MHT have value -1.2
	float deltaHT  = -1.2;
	float deltaMHT = -1.2;
	if (pfHT  > 0){ deltaHT  = (caloHT  - pfHT  )/pfHT; }
	if (pfMHT > 0){ deltaMHT = (caloMHT - pfMHT )/pfMHT;}
	
	float pfAlphaTStandard = calculateAlphaT( genJetPT, genJetPx, genJetPy, genJetThreshold );
	float deltaAlphaT = -1.2;
	if (pfAlphaTStandard > 0){ deltaAlphaT = ( caloAlphaTStandard - pfAlphaTStandard )/pfAlphaTStandard; }
	

	histMatch2D["NPileup_vs_NVTX"]       ->Fill( NVTX,  nPileup );
	histMatch2D["HLTHT_vs_GENHT"]        ->Fill( pfHT,  caloHT );
	histMatch2D["HLTMHT_vs_GENMHT"]      ->Fill( pfMHT, caloMHT );
	histMatch2D["HLTAlphaT_vs_GENAlphaT"]->Fill( pfAlphaTStandard, caloAlphaTStandard );

	histMatch2D["DeltaHT_vs_GENHT"]        ->Fill( pfHT,    deltaHT  );
	histMatch2D["DeltaMHT_vs_GENMHT"]      ->Fill( pfMHT,   deltaMHT );
	histMatch2D["DeltaAlphaT_vs_GENAlphaT"]->Fill( pfAlphaTStandard, deltaAlphaT );
	histMatch2D["DeltaHT_vs_NPileup"]      ->Fill( nPileup, deltaHT  );
	histMatch2D["DeltaMHT_vs_NPileup"]     ->Fill( nPileup, deltaMHT );
	histMatch2D["DeltaAlphaT_vs_NPileup"]  ->Fill( nPileup, deltaAlphaT );

	


      }
    
    }










    histMatch2D["HLTCaloHT_vs_HLTPFHT"] ->Fill( caloHT, hltCaloHT );
    histMatch2D["HLTCaloHT_vs_GenHT"]   ->Fill( pfHT, hltCaloHT );
    histMatch2D["HLTPFHT_vs_GenHT"]     ->Fill( pfHT, caloHT );

    // AlphaT Correlations 
    // histMatch2D["HLTCaloAlphaT_vs_HLTPFAlphaT"]    ->Fill( caloAlphaT, hltCaloAlphaT );
    // histMatch2D["HLTCaloAlphaT_vs_GenAlphaT"]        ->Fill( pfAlphaT,   hltCaloAlphaT );
    // histMatch2D["HLTPFAlphaT_vs_GenAlphaT"]          ->Fill( pfAlphaT,   caloAlphaT );
    // histMatch2D["HLTCaloAlphaT_vs_HLTCaloAlphaTDyn"] ->Fill( hltCaloAlphaTDynamic, hltCaloAlphaT );
    // histMatch2D["HLTCaloAlphaTDyn_vs_HLTPFAlphaT"]   ->Fill( caloAlphaT,           hltCaloAlphaTDynamic );
    // histMatch2D["HLTCaloAlphaTDyn_vs_GenAlphaT"]     ->Fill( pfAlphaT,             hltCaloAlphaTDynamic );
    
    //    histMatch2D["GenAlphaT_vs_GenHT"]

#endif




#ifdef TEST
    if ( iEvent > maxEvents ){ break; } // Exit early
#endif

    if ( !(iEvent % 10000) ){ std::cout << "Event " << iEvent << "\n";}

}
#endif



  // ********************************************************************************
  // *                               Store histograms                               *
  // ********************************************************************************


  std::cout << "Storing histograms\n";
  
  fOut->mkdir("Raw");
  fOut->mkdir("Raw/Efficiency");
  fOut->mkdir("Raw/PrefilterEfficiency");
  fOut->mkdir("Raw/PrefilterEfficiencyRaw");
  //  fOut->mkdir("Raw/HLTEfficiency");
  fOut->mkdir("Raw/HLTOverlap");
  fOut->mkdir("Raw/HLTOverlapRaw");
  fOut->mkdir("Raw/HLTRateRaw");
  fOut->mkdir("Raw/HLTRate");
  fOut->mkdir("Raw/Rate");
  fOut->mkdir("Raw/RateDifferential");
  // fOut->mkdir("Raw/Correlations");
  // fOut->mkdir("Raw/Distributions");






  std::cout << "\thist1DEff\n";
  for(std::map<TString, TEfficiency*>::const_iterator itr2 = hist1DEff.begin(); itr2 != hist1DEff.end(); ++itr2){

    TString histoName = itr2->first; //Extract the histogram key 
    fOut->cd("Raw/Efficiency");
    itr2->second->Write();
   
  }


  std::cout << "\thist2DEff\n";
  for(std::map<TString, TEfficiency*>::const_iterator itr2 = hist2DOnEff.begin(); itr2 != hist2DOnEff.end(); ++itr2){

    fOut->cd("Raw/Efficiency");
    TString histoName = itr2->first; //Extract the histogram key 
    itr2->second->Write();

    TEfficiency *eff      = (TEfficiency*)itr2->second     ->Clone();
    TH2* passed           = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* total            = (TH2*)eff->GetTotalHistogram() ->Clone();
    TH2* passedCumul      = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* totalUniCumul    = (TH2*)eff->GetTotalHistogram() ->Clone();

    // Uniform
    TEfficiency *effUniCumul = (TEfficiency*)eff->Clone();

    reverseCumulative2D( passed, passedCumul, 1 );

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


    effUniCumul->SetTotalHistogram(  *totalUniCumul, "" );
    effUniCumul->SetPassedHistogram( *passedCumul,   "" );
    effUniCumul->SetName( histoName + "_UniCumul" );
    effUniCumul->Write();

   
  }


  std::cout << "\thist1DRate\n";
  for(std::map<TString, TH1*>::const_iterator itr = hist1DRate.begin(); itr != hist1DRate.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/Rate");
    itr->second->Write();

    // Make rate distribution
    TH1* rateHist = (TH1D*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->GetYaxis()->SetTitle("Rate (Hz)");
    reverseCumulative( itr->second, rateHist, rateScaleFactor );
    rateHist->Write();
   
  }


  std::cout << "\thist2DRate\n";
  for(std::map<TString, TH2*>::const_iterator itr = hist2DRate.begin(); itr != hist2DRate.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/RateDifferential");
    itr->second->Write();

    fOut->cd("Raw/Rate");
    TH2* cumulHist = (TH2F*)itr->second->Clone();
    cumulHist->SetName( histoName + "_Cumul" );
    reverseCumulative2D( itr->second, cumulHist, rateScaleFactor );

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


   
  }



  // ********************************************************************************
  // *                                     HLT                                      *
  // ********************************************************************************
  std::cout << "\thistHLTRate\n";
  for(std::map<TString, TH1*>::const_iterator itr = histHLTRate.begin(); itr != histHLTRate.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/HLTRateRaw");
    itr->second->Write();

    fOut->cd("Raw/HLTRate");
    TH1* rateHist = (TH1F*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->GetYaxis()->SetTitle("Rate (Hz)");
    rateHist->Scale( rateScaleFactor );
    rateHist->Write();
   
  }
  std::cout << "\thistHLTRate2D\n";
  for(std::map<TString, TH2*>::const_iterator itr = histHLTRate2D.begin(); itr != histHLTRate2D.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/HLTRateRaw");
    itr->second->Write();

    fOut->cd("Raw/HLTRate");
    TH2* rateHist = (TH2F*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->Scale( rateScaleFactor );
    rateHist->Write();
   
  }

  std::cout << "\thist2DOverlap\n";
  for(std::map<TString, TEfficiency*>::const_iterator itr = hist2DOverlap.begin(); itr != hist2DOverlap.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key                                                                                       
    fOut->cd("Raw/HLTOverlapRaw");
    itr->second->Write();

    fOut->cd("Raw/HLTOverlap");
    TCanvas *c = new TCanvas(histoName);
    TEfficiency* histogramEff = (TEfficiency*)itr->second->Clone();
    histogramEff->Draw("COLZTEXTE");
    gPad->Update();
    histogramEff->GetPaintedHistogram()->SetMaximum(1.0);


    if (histoName.Contains("SM_vs_LegacyAlphaT_") ){
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 1, "Not fired");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 2, "Fired");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 3, "H_{T} = 200");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 4, "H_{T} = 250");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 5, "H_{T} = 300");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 6, "H_{T} = 350");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 7, "H_{T} = 400");
      
      histogramEff->GetPaintedHistogram()->GetYaxis()->SetBinLabel( 1, "Not fired");
      histogramEff->GetPaintedHistogram()->GetYaxis()->SetBinLabel( 2, "Fired");
    }
    else if(histoName.Contains("LegacyAlphaT_vs_SM_") ){
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 1, "Not fired");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 2, "Fired");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 3, "PFMET170");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 4, "PFMET120_Btag");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 5, "PFHT350_PFMet120");
      histogramEff->GetPaintedHistogram()->GetXaxis()->SetBinLabel( 6, "PFHT900");
      
      histogramEff->GetPaintedHistogram()->GetYaxis()->SetBinLabel( 1, "Not fired");
      histogramEff->GetPaintedHistogram()->GetYaxis()->SetBinLabel( 2, "Fired");
    }
    
    histogramEff->GetPaintedHistogram()->Draw("COLZTEXTE");
    c->Write();

  }



  std::cout << "\thist2DHLTEff\n";
  for(std::map<TString, TEfficiency*>::const_iterator itr = hist2DHLTEff.begin(); itr != hist2DHLTEff.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/PrefilterEfficiencyRaw");
    //itr->second->Write();

    TEfficiency *eff      = (TEfficiency*)itr->second      ->Clone();
    TH2* passed           = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* passedCumul      = (TH2*)eff->GetPassedHistogram()->Clone();
    TH2* passedUni        = (TH2*)eff->GetPassedHistogram()->Clone();

    // Uniform
    TEfficiency *effUniCumul = (TEfficiency*)eff->Clone();
    reverseCumulative2D( passed, passedCumul, 1 );

    // The numerator contains all the information, the denominator means little, we only consider entries in the numberator!!!
    fillUniform2D( passedUni, passed->GetEntries() );

    effUniCumul->SetTotalHistogram(  *passedUni,   "" );
    effUniCumul->SetPassedHistogram( *passedCumul, "" );
    effUniCumul->SetName( histoName + "_UniCumul" );
    effUniCumul->Write();

    fOut->cd("Raw/PrefilterEfficiency");
    TCanvas *c = new TCanvas(histoName);
    TEfficiency* histogramEff = (TEfficiency*)effUniCumul->Clone();
    histogramEff->Draw("COLZTEXTE");
    gPad->Update();
    histogramEff->GetPaintedHistogram()->SetMaximum(1.0);
    histogramEff->GetPaintedHistogram()->GetXaxis()->SetRangeUser(0,     500);
    histogramEff->GetPaintedHistogram()->GetYaxis()->SetRangeUser(0.41, 0.61);
    histogramEff->GetPaintedHistogram()->Draw("COLZTEXTE");
    c->Write();


  } 

  // ********************************************************************************
  // *                                   Matching                                   *
  // ********************************************************************************

#ifdef MATCHING
  fOut->mkdir("Raw/MatchRaw");
  fOut->mkdir("Raw/Match");

  std::cout << "\thistMatch\n";
  for(std::map<TString, TH1*>::const_iterator itr = histMatch.begin(); itr != histMatch.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/MatchRaw");
    itr->second->Write();

    fOut->cd("Raw/Match");
    TH1* rateHist = (TH1F*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->GetYaxis()->SetTitle("Rate (Hz)");
    rateHist->Scale( rateScaleFactor );
    rateHist->Write();
   
  }
  std::cout << "\thistMatch2D\n";
  for(std::map<TString, TH2*>::const_iterator itr = histMatch2D.begin(); itr != histMatch2D.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/MatchRaw");
    itr->second->Write();

    fOut->cd("Raw/Match");
    TH2* rateHist = (TH2F*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->Scale( rateScaleFactor );
    rateHist->Write();
   
  }
#endif

  exit(0);



    

  fOut->Close();
}


// ****************************************************************************************************
// ****************************************************************************************************
// ****************************************************************************************************

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
  //  int xBinOver        = nBinsX + 1;
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
  //  int xBinOver        = nBinsX + 1;
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
  // int xBinOver        = nBinsX + 1;
  // int yBinOver        = nBinsY + 1;

  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){     
    for ( int iBinY = nBinsY; iBinY > 0; --iBinY ){     
      //double integral = histogram->Integral( iBinX, xBinOver, iBinY, yBinOver );
      histogram->SetBinContent( iBinX, iBinY, value );
    }
  }      

}


