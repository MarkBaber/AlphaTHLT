//#define TEST
//#define SIGNAL
#define NEUTRINO
#define HLT_CALOJET


float caloJetThreshold = 50; //50;
float pfJetThreshold   = 50; //50;
float secondJetThreshold = 100;

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
float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx, 
			     std::vector<float> *jetPy, //float jetThreshold, 
			     uint maxJets, 
			     float dynamicJetThreshold, float dynamicAlphaTThreshold);
std::pair<float,float> calculateDynamicAlphaTHT(std::vector<float> *jetPT, std::vector<float> *jetPx, 
						std::vector<float> *jetPy, //float jetThreshold, 
						uint maxJets,
                                                float dynamicJetThreshold, float dynamicAlphaTThreshold, float dynamicHTThreshold);

std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs(std::vector<float> *jetPT, std::vector<float> *jetPx, 
								 std::vector<float> *jetPy, //float jetThreshold, 
								 uint maxJets, float dynamicJetThreshold);


float const  PI        = TMath::Pi();




// ********************************************************************************
// *                                  Selections                                  *
// ********************************************************************************


const int MAX_JETS = 4;



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
  // TString sampleDir = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/11Oct14/";  // Now includes forward jets, leptons, lepton veto, trigger bits
  TString sampleDir = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/17Oct14/";  // pre8 with correct JEC

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


  sample test = sample("test", branch, "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_2_0_pre6/src/AlphaTHLT/MakeTree/test/QCD_Pt-30to50_Tune4C_13TeV_pythia8.root", 161500000);

  // ------------------------------------------------------------------------------------------------------------------------
  sample selectedSample =  QCD30to50; //T2tt_500_250; //QCD30to50;  //test; //QCD30to50; // T2tt_500_250; //T2cc_250_210; //DYJets; //NuGun; //DYJets; //TTBar; //DYJets;

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


  unsigned int nuNEvents  = (unsigned int)rateChain ->GetEntries(); 
  unsigned int sigNEvents = (unsigned int)sigChain  ->GetEntries();
  //  double  rateScaleFactor = double(ZB_XSECTION)/nuNEvents;

  // New scale factor calculation for QCD samples
  double rateScaleFactor = (sampleXS * INST_LUMI_25NS)/nuNEvents;
  if (selectedSample.Name == "NuGun"){
    rateScaleFactor = double(ZB_XSECTION)/nuNEvents;
  }


    // rate_mc_simple = (xs * ilumi * counts)/nevt
    // The equation for the error of your rate is:
    // rateerr_mc = xs * ilumi * ((math.sqrt(counts + ((counts)**2)/nevt))/nevt)
    //   rateerr_mc_simple = ((xs * ilumi)/nevt)*math.sqrt(counts) 
    //   So if you use several samples (like QCD binned in Pt), you 




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
  std::map< TString, TH2*>         histHLTRate2D;

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

  float pfHT(0), pfMET(0), pfMHT(0); // pfAlphaT(0),
  std::vector<float> *caloJetPT  = new std::vector<float>();
  std::vector<float> *caloJetPx  = new std::vector<float>();
  std::vector<float> *caloJetPy  = new std::vector<float>();
  std::vector<float> *caloJetEta = new std::vector<float>();
  //  std::vector<float> *caloJetPhi = new std::vector<float>();
  float caloHT(0); //, caloAlphaT(0);


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
    histHLTRate[ path ]                  = new TH1D( path, path + ";Entries;Fired", 2, 0, 2.);
    histHLTRate[ path  + "_AND_HTT175" ] = new TH1D( path + "_AND_HTT175", path + ";Entries;Fired", 2, 0, 2.);
    histHLTRate[ path  + "_AND_Trig1" ]  = new TH1D( path + "_AND_Trig1", path + ";Entries;Fired", 2, 0, 2.);


    // Get pt hat breakdown of rate of the triggers
    makePTHatBinnedHist( histHLTRate, path );
    makePTHatBinnedHist( histHLTRate, path + "_AND_HTT175");
    makePTHatBinnedHist( histHLTRate, path + "_AND_Trig1");
    

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

    TString stdStr      = "On_AlphaTStd_vs_HT_"  + jet2PTCutStr;
    TString stdStrRate  = "OnRate_AlphaTStd_vs_HT_"  + jet2PTCutStr;
    // TString dynStr  = "On_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
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


    
    histHLTRate2D["HLT_AlphaT_vs_HT"] = new TH2D("HLT_AlphaT_vs_HT","HLT #alpha_{T}^{static} vs H_{T}^{static} differential rate;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


    



    // ------------------------------------------------------------------------------------------------------------------------
    // Inclusive
    // ------------------------------------------------------------------------------------------------------------------------
    hist2DOnEff[stdStr  + "_Inclusive"] = new TEfficiency(stdStr  + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DOnEff[dynStr  + "_Inclusive"] = new TEfficiency(dynStr  + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DOnEff[dyn2Str + "_Inclusive"] = new TEfficiency(dyn2Str + "_Inclusive","Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DOnEff[dyn3Str + "_Inclusive"] = new TEfficiency(dyn3Str + "_Inclusive","Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    hist2DRate[stdStrRate  + "_Inclusive"] = new TH2D(stdStr  + "_InclusiveRate","Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DRate[dynStr  + "_Inclusive"] = new TH2D(dynStr  + "_InclusiveRate","Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DRate[dyn2Str + "_Inclusive"] = new TH2D(dyn2Str + "_InclusiveRate","Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DRate[dyn3Str + "_Inclusive"] = new TH2D(dyn3Str + "_InclusiveRate","Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

    // Make candidate L1 trigger plots
    for (uint trigN = 1; trigN <= 7; ++trigN ){
      TString trigStr = TString("Trig") + Form("%d", trigN);

      hist2DRate[trigStr + stdStrRate  + "_Inclusive"] = new TH2D(trigStr + stdStrRate  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // Excluding overlap rate with SUSY hadronic menu
      if (trigN == 1){
	hist2DRate[trigStr + "SM" + stdStrRate  + "_Inclusive"] = new TH2D(trigStr + "SM" + stdStrRate  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ") - Excluding SUSY hadronic menu overlap;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      }

      // hist2DRate[trigStr + dynStr  + "_Inclusive"] = new TH2D(trigStr + dynStr  + "_Inclusive",trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DRate[trigStr + dyn2Str + "_Inclusive"] = new TH2D(trigStr + dyn2Str + "_Inclusive",trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DRate[trigStr + dyn3Str + "_Inclusive"] = new TH2D(trigStr + dyn3Str + "_Inclusive",trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    }
    hist2DRate["NoL1" + stdStrRate  + "_Inclusive"] = new TH2D("NoL1" + stdStr  + "_Inclusive","No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
    // hist2DRate["NoL1" + dynStr  + "_Inclusive"] = new TH2D("NoL1" + dynStr  + "_Inclusive","No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} rate inclusive (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
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


      // Make candidate HLTtrigger plots                                                                                                            
      for (uint trigN = 1; trigN <= 5; ++trigN ){
	TString trigStr = TString("HLT") + Form("%d", trigN);

	if ( (jet2PTCut == 100.) || (jet2PTCut == 50.) ){
	  hist2DHLTEff[trigStr + stdStr  + "_" + jetBinStr] = new TEfficiency(trigStr + stdStr  + "_" + jetBinStr,"Gen #alpha_{T}^{static} vs H_{T}^{static} efficiency "  + jetBinStr + " (" + jet2PTCutLab + ");Gen H_{T} (GeV);Gen #alpha_{T}",   20.,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	}
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
	  // hist2DOnEff[dynStr  + suffix] = new TEfficiency(dynStr  + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  // hist2DOnEff[dyn2Str + suffix] = new TEfficiency(dyn2Str + suffix,"Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  // hist2DOnEff[dyn3Str + suffix] = new TEfficiency(dyn3Str + suffix,"Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  
  
    // Make candidate L1 trigger plots
    for (uint trigN = 1; trigN <= 7; ++trigN ){
      TString trigStr = TString("Trig") + Form("%d", trigN);

	hist2DOnEff[trigStr + stdStr  + suffix] = new TEfficiency(trigStr + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	if (trigN == 1){
	  hist2DOnEff[trigStr + "SM" + stdStr  + suffix] = new TEfficiency(trigStr + "SM" + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ") - Including SUSY hadronic menu overlap;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	  hist2DOnEff[trigStr + "SMOnly" + stdStr  + suffix] = new TEfficiency(trigStr + "SM" + stdStr  + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ") - SUSY hadronic menu only;H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	}

	// hist2DOnEff[trigStr + dynStr  + suffix] = new TEfficiency(trigStr + dynStr  + suffix,trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	// hist2DOnEff[trigStr + dyn2Str + suffix] = new TEfficiency(trigStr + dyn2Str + suffix,trigStr + " Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
	// hist2DOnEff[trigStr + dyn3Str + suffix] = new TEfficiency(trigStr + dyn3Str + suffix,trigStr + " Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);

    }
      
      hist2DOnEff["NoL1" + stdStr  + suffix] = new TEfficiency("NoL1" + stdStr  + suffix,"No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",   HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff["NoL1" + dynStr  + suffix] = new TEfficiency("NoL1" + dynStr  + suffix,"No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{static} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff["NoL1" + dyn2Str + suffix] = new TEfficiency("NoL1" + dyn2Str + suffix,"No L1 - Calojet #alpha_{T}^{dynamic} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}", HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);
      // hist2DOnEff["NoL1" + dyn3Str + suffix] = new TEfficiency("NoL1" + dyn3Str + suffix,"No L1 - Calojet #alpha_{T}^{static} vs H_{T}^{dynamic} efficiency " + label + " (" + jet2PTCutLab + ");H_{T} (GeV);#alpha_{T}",  HTBins,HTMin,HTMax,  alphaTBins,alphaTMin,alphaTMax);


      }
    } // End analysis bin
    
    // ------------------------------------------------------------------------------------------------------------------------
    // ------------------------------------------------------------------------------------------------------------------------

  } // End second jet cut




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

    sigChain->SetBranchAddress("genMetCalo_MetPt",     &pfMET);

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

  std::cout << "Sample xs (pb)         = " << selectedSample.CrossSection     << "\n";
  std::cout << "Sample xs (cms2)       = " << sampleXS                        << "\n";
  std::cout << "Rate scalefactor       = " << rateScaleFactor                 << "\n";
  std::cout << "Rate scalefactor error = " << rateScaleFactor*sqrt(nuNEvents) << "\n\n";


#ifdef SIGNAL

  // Loop over the tree
  //  for ( unsigned int iEvent = 0; iEvent < nEvents; ++iEvent ){
  for ( unsigned int iEvent = signalEventLow; iEvent < signalEventHigh; ++iEvent ){
    
    sigChain->GetEntry( iEvent );



    float pfAlphaTStandard(0);

    float caloAlphaTStandard(0);
    // float caloAlphaTDynamic(0);
    // float caloHTDynamic(0);
    //    std::pair<float,float> caloAlphaTHTDynamic;
    

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
    float pfMHTX(0), pfMHTY(0);
    pfHT = 0;
    for(uint iJet = 0;iJet < pfJetPT->size(); ++iJet ){
      if( (*pfJetPT)[iJet] < pfJetThreshold ) break;
      pfHT += (*pfJetPT)[iJet];
      pfMHTX += (*pfJetPx)[iJet];
      pfMHTY += (*pfJetPy)[iJet];
      pfJetsAboveThresh++;
    }
    pfMHT = sqrt( pfMHTX*pfMHTX + pfMHTY*pfMHTY);
    if ( pfJetsAboveThresh > MAX_JETS ){
      pfJetsAboveThresh = MAX_JETS;
    }
    //    TString nPFJetsStr = Form("%d", pfJetsAboveThresh );
    TString pfJetBinStr     = TString(Form( "%d", pfJetsAboveThresh)) + TString("Jets");


    // Calo HLT cuts
    // bool passesATStandard(false), passesATDynamic(false), passesATDynamic2(false), passesATDynamic3(false);
    // bool passesHTStandard(false), passesHTDynamic(false);
    //    bool passesJet2_50(false), passesJet2_60(false), passesJet2_70(false), passesJet2_80(false), passesJet2_90(false), passesJet2_100(false);
    //    bool passesDynJet2_100(false);    

    //    bool passesOffHT(false), passesOffAT(false), passesOffJet(false), passesOffAll(false);
    bool passesOffJet(false);
    bool passesAnaBinOffAT(false), passesAnaBinOffJet(false), passesAnaBinOffAll(false);
    bool passesOffVetoes(false);


    pfAlphaTStandard    = calculateAlphaT( pfJetPT, pfJetPx, pfJetPy, pfJetThreshold );      

    caloAlphaTStandard  = calculateAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold );    
    // // Dynamic AlphaT
    // caloAlphaTDynamic   = calculateDynamicAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold, maxCaloJet, caloJetDynThreshold, 
    // 						  caloJetAlphaThreshold );




    // ********************************************************************************
    // *                              Calojets 
    // ********************************************************************************
    if ( pfJetsAboveThresh >= 2 ){

      
      // Forward jet veto
      bool forJetVeto(false);
      if (pfJetForPT->size() > 0){
	if ( (*pfJetForPT)[0] > pfJetThreshold){
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
      // if ( pfHT             > pfJetHTThreshold)                                     { passesOffHT     = true; }
      // if ( pfAlphaTStandard > pfJetAlphaTThreshold)                                 { passesOffAT     = true; }
      if ( (pfJetsAboveThresh >= 2 ) && ((*pfJetPT)[1] > pfJet2PTThreshold) )          { passesOffJet    = true; }
      if ( !(genLeptonVeto) && !(forJetVeto) && !(mhtOverMetVeto) )                    { passesOffVetoes = true; }
      // Check full offline selection
      //      if ( passesOffHT && passesOffAT && passesOffJet && passesOffVetoes )     { passesOffAll    = true; }
      // Vary the online second jet PT cut
      for (uint iJet2Cut = 0; iJet2Cut < jet2PTCuts.size(); ++iJet2Cut){

	float jet2PTCut = jet2PTCuts[ iJet2Cut ];
	TString jet2PTCutStr = "Jet2gt" + TString(Form("%1.0f", jet2PTCut ));
	TString stdStr  = "On_AlphaTStd_vs_HT_"  + jet2PTCutStr;
	// TString dynStr  = "On_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
	// TString dyn2Str = "On_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
	// TString dyn3Str = "On_AlphaTDyn3_vs_HT_" + jet2PTCutStr;

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
	  if ( (pfJetsAboveThresh >= 2 ) && ((*pfJetPT)[1] > pfJet2PTThreshold) ){ passesAnaBinOffJet = true; }
	  // Check full offline selection
	  if ( passesAnaBinOffAT && passesAnaBinOffJet )                         { passesAnaBinOffAll = true; }

	  // Get jet bin
	  for (uint iJetBin = 0; iJetBin < anaJetBins.size(); ++iJetBin){
	
	    int jetLow  = anaJetBins[ iJetBin ].first;
	    int jetHigh = anaJetBins[ iJetBin ].second;
	    if ( !((pfJetsAboveThresh >= jetLow) && (pfJetsAboveThresh < jetHigh)) ){ continue; }
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
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig1";
	      trigBool = l1Trig1 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );

	      bool susyMenu = ((hltPathFired["HLT_PFMET170_NoiseCleaned_NoL1_v1"])           || 
			       (hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_NoL1_v1"]) ||
			       (hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_NoL1_v1"] )  ||
			       (hltPathFired["HLT_PFHT900_NoL1_v1"]) );
	      // Overlap efficiency with SUSY hadronic menu
	      hist2DOnEff[trigStr + "SM" + stdStr + suffix]->Fill( (susyMenu||trigBool), caloHT, caloAlphaTStandard );
	      // Efficiency of SUSY hadronic menu alone
	      hist2DOnEff[trigStr + "SMOnly" + stdStr + suffix]->Fill( susyMenu, caloHT, caloAlphaTStandard );
	      
	      



	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig2";
	      trigBool = l1Trig2 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig3";
	      trigBool = l1Trig3 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig4";
	      trigBool = l1Trig4 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig5";
	      trigBool = l1Trig5 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig6";
	      trigBool = l1Trig6 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      trigStr  = "Trig7";
	      trigBool = l1Trig7 && passHltSecondJet;
	      hist2DOnEff[trigStr + stdStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTStandard );
	      // hist2DOnEff[trigStr + dynStr  + suffix]->Fill( trigBool, caloHT,                     caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn2Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTDynamic );
	      // hist2DOnEff[trigStr + dyn3Str + suffix]->Fill( trigBool, caloAlphaTHTDynamic.second, caloAlphaTStandard );

	      
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
    //rateChain->SetBranchAddress("hltAk4PF_Phi",       &caloJetPhi);
    rateChain->SetBranchAddress("hltAk4PF_Eta",       &caloJetEta);
#endif



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
    bool l1Trig2 = ( (uctHT >= 200)  || (uctMET >= 60) ); // DJ120
    bool l1Trig3 = ( (uctHT >= 160)  || (uctMET >= 70) ); // DJ120
    bool l1Trig4 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 150) && (uctMHToverHT >= 0.17)) );
    bool l1Trig5 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 140) && (uctMHToverHT >= 0.25)) );
    bool l1Trig6 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 130) && (uctMHToverHT >= 0.32)) );
    bool l1Trig7 = ( (uctHT >= 200)  || (uctMET >= 70) || ((uctHT >= 120) && (uctMHToverHT >= 0.36)) );
    
    bool l1HTT175 = (uctHT  >= 175);
    bool l1ETM70  = (uctMET >= 70);
    // ****************************************************************************************************
    // *                                         Check HLTriggers                                         *
    // ****************************************************************************************************
    // bool HLT_PFMET170_NoiseCleaned_v1           = (l1ETM70 && hltPathFired["HLT_PFMET170_NoiseCleaned_NoL1_v1"]);
    // bool HLT_PFMET120_NoiseCleaned_BTagCSV07_v1 = (l1ETM70 && hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_NoL1_v1"]);
    // bool HLT_PFHT350_PFMET120_NoiseCleaned_v1   = (l1HT175 && hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_NoL1_v1"]);
    // bool HLT_PFHT900_v1                         = (l1HT175 && hltPathFired["HLT_PFHT900_NoL1_v1"]);

    bool firesSUSYMenu                          = ((hltPathFired["HLT_PFMET170_NoiseCleaned_v1"])           ||
						   (hltPathFired["HLT_PFMET120_NoiseCleaned_BTagCSV07_v1"]) ||
						   (hltPathFired["HLT_PFHT350_PFMET120_NoiseCleaned_v1"] )  ||
						   (hltPathFired["HLT_PFHT900_v1"]) );


    // ************************************************************
    // Trigger bits 
    // ************************************************************
    for (uint iPath = 0; iPath < hltPathNames.size(); ++iPath){
      TString path = hltPathNames[ iPath ];
      histHLTRate[path] ->Fill( hltPathFired[ path ] );
      histHLTRate[path + "_PTHat"]->Fill( hltPathFired[ path ] );
      histHLTRate[path + "_PTHat"]->Fill( samplePTHat, hltPathFired[ path ] );


      if ( l1HTT175 ){
	histHLTRate[ path + "_AND_HTT175" ]     ->Fill( hltPathFired[ path ] );
	histHLTRate[ path + "_AND_HTT175_PTHat"]->Fill( hltPathFired[ path ] );
	histHLTRate[ path + "_AND_HTT175_PTHat"]->Fill( samplePTHat, hltPathFired[ path ] );
      }
      if ( l1Trig1 ){
	histHLTRate[ path + "_AND_Trig1"]      ->Fill( hltPathFired[ path ] );
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
    float caloAlphaTStandard(0);
    caloHT = 0;
    caloAlphaTStandard  = calculateAlphaT( caloJetPT, caloJetPx, caloJetPy, caloJetThreshold );

    // Calojets
    int caloJetsAboveThresh(0);
    for(uint iJet = 0;iJet < caloJetPT->size(); ++iJet ){
      if( (*caloJetPT)[iJet] < caloJetThreshold ) break;
      caloHT    += (*caloJetPT)[iJet];

      // if ( (caloJetPT->size() > 1) && ((*caloJetPT)[iJet] > secondJetThreshold) ){
      // 	caloMHTx  += (*caloJetPx)[iJet];
      // 	caloMHTy  += (*caloJetPy)[iJet];
      // 	caloHTJet2Cut += (*caloJetPT)[iJet];
      // }

      caloJetsAboveThresh++;
    }
    if ( caloJetsAboveThresh > MAX_JETS ){
      caloJetsAboveThresh = MAX_JETS;
    }
    TString nCaloJetsStr = Form("%d", caloJetsAboveThresh );





    
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
	  // TString dynStr  = "OnRate_AlphaTDyn_vs_HT_"  + jet2PTCutStr;
	  // TString dyn2Str = "OnRate_AlphaTDyn2_vs_HT_" + jet2PTCutStr;
	  // TString dyn3Str = "OnRate_AlphaTDyn3_vs_HT_" + jet2PTCutStr;
	  
	  // Inclusive 
	  hist2DRate[stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  // hist2DRate[dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  // hist2DRate[dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  // hist2DRate[dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );

	  // ********************************************************************************
	  // Emulate AlphaT HLTrigger
	  // ********************************************************************************
	  if (l1Trig1){
	    if ( (caloHT >= 200) && (caloAlphaTStandard >= 0.57) ){ legacyAlphaTHT200 = true; }
	    if ( (caloHT >= 250) && (caloAlphaTStandard >= 0.55) ){ legacyAlphaTHT250 = true; }
	    if ( (caloHT >= 300) && (caloAlphaTStandard >= 0.53) ){ legacyAlphaTHT300 = true; }
	    if ( (caloHT >= 350) && (caloAlphaTStandard >= 0.52) ){ legacyAlphaTHT350 = true; }
	    if ( (caloHT >= 400) && (caloAlphaTStandard >= 0.51) ){ legacyAlphaTHT400 = true; }
	  }
	  firesLegacyAlphaT = ( legacyAlphaTHT200 || legacyAlphaTHT250 || 
				legacyAlphaTHT300 || legacyAlphaTHT350 ||
				legacyAlphaTHT400 );
	  
	  
	  // Make candidate L1 trigger plots
	  TString trigStr = "";
	  trigStr = "NoL1";
	  hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	  // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	  // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	  // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );

	  if ( l1Trig1){
	    trigStr = "Trig1";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]           ->Fill( caloHT,                     caloAlphaTStandard );


	    // Exclude overlap rate with SUSY hadronic menu - assuming SUSY menu has same HTT175||MET70 trigger
	    if ( !( firesSUSYMenu ) ){
	      hist2DRate[trigStr + "SM" + stdStr + "_Inclusive"]->Fill( caloHT,                     caloAlphaTStandard );
	    }

	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }
	  if ( l1Trig2){
	    trigStr = "Trig2";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }
	  if ( l1Trig3){
	    trigStr = "Trig3";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }
	  if ( l1Trig4){
	    trigStr = "Trig4";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }
	  if ( l1Trig5){
	    trigStr = "Trig5";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }
	  if ( l1Trig6){
	    trigStr = "Trig6";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }
	  if ( l1Trig7){
	    trigStr = "Trig7";
	    hist2DRate[trigStr + stdStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTStandard );
	    // hist2DRate[trigStr + dynStr  + "_Inclusive"]       ->Fill( caloHT,                     caloAlphaTDynamic );
	    // hist2DRate[trigStr + dyn2Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTDynamic ); 
	    // hist2DRate[trigStr + dyn3Str + "_Inclusive"]       ->Fill( caloAlphaTHTDynamic.second, caloAlphaTStandard );
	  }




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

    

    histHLTRate2D["HLT_AlphaT_vs_HT"]->Fill(caloHT, caloAlphaTStandard );





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
  fOut->mkdir("Raw/HLTEfficiency");
  fOut->mkdir("Raw/HLTRate");
  fOut->mkdir("Raw/Rate");
  fOut->mkdir("Raw/Correlations");
  fOut->mkdir("Raw/Distributions");



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
    fOut->cd("Raw/Rate");
    itr->second->Write();

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
    fOut->cd("Raw/HLTRate");
    itr->second->Write();

    TH1* rateHist = (TH1F*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->GetYaxis()->SetTitle("Rate (Hz)");
    rateHist->Scale( rateScaleFactor );
    rateHist->Write();
   
  }
  std::cout << "\thistHLTRate2D\n";
  for(std::map<TString, TH2*>::const_iterator itr = histHLTRate2D.begin(); itr != histHLTRate2D.end(); ++itr){

    TString histoName = itr->first; //Extract the histogram key 
    fOut->cd("Raw/HLTRate");
    itr->second->Write();

    TH2* rateHist = (TH2F*)itr->second->Clone();
    rateHist->SetName( histoName + "_TrueRate" );
    rateHist->Scale( rateScaleFactor );
    rateHist->Write();
   
  }
 


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
std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs(std::vector<float> *jetPT, std::vector<float> *jetPx, 
								 std::vector<float> *jetPy,// float jetThreshold, 
								 uint maxJets, float dynamicJetThreshold){

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





std::pair<float,float> calculateDynamicAlphaTHT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy,//float jetThreshold,
 						uint maxJets, 
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
float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx, 
			     std::vector<float> *jetPy, //float jetThreshold, 
			     uint maxJets, 
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
