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
#include "TMath.h"


#include "typeDefs.h"
#include "ttreeDraw.C"


//#define MHT_VS_MoM
#define MHT_VS_MET
//#define MHT_VS_HT
//#define RESOLUTIONS

int readTree(){

  std::cout << "\n\nLoading chains\n";
  TString dir      = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/";  
  TString outDir   = "/home/hep/mb1512/public_html/plots/31_01_15/";
  TString filename = "MHT_vs_MET.root";
  TFile *f = new TFile(outDir + filename,"RECREATE");

  // Switches
  // -------------------------------------------------- 
  bool runRate       = true;
  bool runBackground = true;
  bool runSignal     = true;

  // Define samples
  // -------------------------------------------------- 
  sample QCD          = sample("QCD",          dir + "12Nov14QCDReweighted_New/*.root",                            runRate);
  sample ttBar        = sample("ttBar",        dir + "/12Nov14_New/TT_Tune4C_13TeV-pythia8-tauola/*.root",         runBackground );
  sample DY           = sample("DY",           dir + "/12Nov14_New/DYJetsToLL_M-50_13TeV-madgraph-pythia8/*.root", runBackground );
  sample t2tt_425_325 = sample("t2tt_425_325", dir + "/16Jan15_20PU25ns_New/T2tt_2J_mStop-425_mLSP-325/*.root",    runSignal );
  sample t2tt_500_325 = sample("t2tt_500_325", dir + "/16Jan15_20PU25ns_New/T2tt_2J_mStop-500_mLSP-325/*.root",    runSignal );
  sample t2tt_650_325 = sample("t2tt_650_325", dir + "/16Jan15_20PU25ns_New/T2tt_2J_mStop-650_mLSP-325/*.root",    runSignal );
  sample t2tt_850_100 = sample("t2tt_850_100", dir + "/16Jan15_20PU25ns_New/T2tt_2J_mStop-850_mLSP-100/*.root",    runSignal );
  sample t2qq_600_550 = sample("t2qq_600_550", dir + "/16Jan15_20PU25ns_New/T2qq_2J_mStop-600_mLSP-550/*.root",    runSignal );
  sample t2bb_600_580 = sample("t2bb_600_580", dir + "/16Jan15_20PU25ns_New/T2bb_2J_mStop-600_mLSP-580/*.root",    runSignal );

  sampleCollection rate; 
  rate.push_back( QCD );
  sampleCollection background; 
  background.push_back( ttBar );
  background.push_back( DY );
  sampleCollection signal; 
  signal.push_back( t2tt_425_325 );
  signal.push_back( t2tt_500_325 );
  signal.push_back( t2tt_650_325 );
  signal.push_back( t2tt_850_100 );
  signal.push_back( t2qq_600_550 );
  signal.push_back( t2bb_600_580 );

  // Binning
  // -------------------------------------------------- 
  int   HTBins(36);     float HTMin(100),     HTMax(1000);
  int   MHTBins(180);   float MHTMin(0),      MHTMax(1800);
  int   METBins(180);   float METMin(0),      METMax(1800);
  int   etaBins(60);    float etaMin(0.0),    etaMax(3.0);
  int   momBins(100);   float momMin(0),      momMax(10);
  int   ForBins(80);    float ForMin(0),      ForMax(400);
  int   alphaTBins(60); float alphaTMin(0.4), alphaTMax(1.0);
  int   DPhiBins(32);   float DPhiMin(0),     DPhiMax(3.2);

 // Histogram templates
 // --------------------------------------------------
 TH2D* at_vs_ht         = new TH2D("at_vs_ht",  ";HLT H_{T} (GeV);HLT #alpha_{T} (GeV)", 
				   HTBins,HTMin,HTMax, alphaTBins,alphaTMin,alphaTMax);
 TH2D* for_vs_mom       = new TH2D("for_vs_mom",";HLT #slash{H}_{T}/#slash{E}_{T};HLT leading forward jet p_{T} (GeV)", 
				   momBins,momMin,momMax, ForBins,ForMin,ForMax);
 TH2D* jet2_vs_ht       = new TH2D("jet2_vs_ht",";HLT H_{T} (GeV);HLT second leading jet p_{T} (GeV)", 
				   HTBins,HTMin,HTMax, ForBins,ForMin,ForMax);
 TH2D* for_vs_eta       = new TH2D("for_vs_eta",";HLT leading jet #eta;HLT leading forward jet p_{T} (GeV)", 
				   etaBins,etaMin,etaMax, ForBins,ForMin,ForMax);
 TH2D* mom_vs_dPhiMM    = new TH2D("mom_vs_dPhiMM",";HLT #Delta#phi(#slash{H}_{T}, #slash{E}_{T});HLT #slash{H}_{T}/#slash{E}_{T}", 
				   DPhiBins,DPhiMin,DPhiMax, momBins,momMin,momMax);
 TH2D* dPhiMM_vs_dPhiMM = new TH2D("dPhiMM_vs_dPhiMM",";HLT #Delta(#slash{H}_{T}, #slash{E}_{T});Gen #Delta(#slash{H}_{T}, #slash{E}_{T})",
				   DPhiBins,DPhiMin,DPhiMax, DPhiBins,DPhiMin,DPhiMax);
 TH2D* mht_vs_met       = new TH2D("mht_vs_met",";HLT #slash{E}_{T} (GeV);HLT #slash{H}_{T} (GeV)",
				   MHTBins,MHTMin,MHTMax, MHTBins,MHTMin,MHTMax);
 TH2D* mht_vs_mom       = new TH2D("mht_vs_mom",";HLT #slash{H}_{T}/#slash{E}_{T};HLT #slash{H}_{T} (GeV)",
				   momBins,momMin,momMax, MHTBins,MHTMin,MHTMax);
 TH2D* mht_vs_ht        = new TH2D("mht_vs_ht",";HLT #slash{H}_{T} (GeV);HLT {H}_{T} (GeV)",
				   HTBins,HTMin,HTMax, MHTBins,MHTMin,MHTMax);


 // Selections
 // --------------------------------------------------
 TString L1        = "((gct_Ht>=175)||(gct_MetPt>=70))";
 TString HTLoose40 = "((hltAk4PF_AlphaT40>0.0) &&(hltAk4PF_HT40>150)&&(hltAk4PF_Pt[1]>40))";
 TString WPLoose40 = "((hltAk4PF_AlphaT40>0.50)&&(hltAk4PF_HT40>150)&&(hltAk4PF_Pt[1]>40))";
 TString WPLoose70 = "((hltAk4PF_AlphaT40>0.50)&&(hltAk4PF_HT40>150)&&(hltAk4PF_Pt[1]>70))";
 TString WPLoose80 = "((hltAk4PF_AlphaT40>0.50)&&(hltAk4PF_HT40>150)&&(hltAk4PF_Pt[1]>80))";
 TString WPLoose90 = "((hltAk4PF_AlphaT40>0.50)&&(hltAk4PF_HT40>150)&&(hltAk4PF_Pt[1]>90))";
 TString WP1       = "((hltAk4PF_AlphaT40>0.60)&&(hltAk4PF_HT40>200)&&(hltAk4PF_Pt[1]>90))";
 TString WP5       = "((hltAk4PF_AlphaT40>0.52)&&(hltAk4PF_HT40>350)&&(hltAk4PF_Pt[1]>90))";

 TString vetoes   = "(!(genLeptonVeto && (genAk4_MhtPT40/genMetCaloAndNonPrompt > 1.25)))";
 TString offGen1  = "(   genAk4_Pt[1]>100 &&    genAk4_HT40>200 &&    genAk4_AlphaT40>0.65 &&    genAk4For_Pt[0]<40)&&" + vetoes;
 TString offGen5  = "(   genAk4_Pt[1]>100 &&    genAk4_HT40>400 &&    genAk4_AlphaT40>0.53 &&    genAk4For_Pt[0]<40)&&" + vetoes;
 TString offRECO1 = "(recoAk4PF_Pt[1]>100 && recoAk4PF_HT40>200 && recoAk4PF_AlphaT40>0.65 && recoAk4PFFor_Pt[0]<40)&&" + vetoes;
 TString offRECO5 = "(recoAk4PF_Pt[1]>100 && recoAk4PF_HT40>400 && recoAk4PF_AlphaT40>0.53 && recoAk4PFFor_Pt[0]<40)&&" + vetoes;
 // Asymmetric selection
 TString offAGen1   = offGen1;
 offAGen1.ReplaceAll( "genAk4_Pt[1]>100","(genAk4_Pt[1]>100||(genAk4_Pt[1]>40&&genAk4_Pt[1]<=100))");
 TString offAGen5   = offGen5;
 offAGen5.ReplaceAll( "genAk4_Pt[1]>100","(genAk4_Pt[1]>100||(genAk4_Pt[1]>40&&genAk4_Pt[1]<=100))");
 TString offARECO1  = offRECO1;
 offARECO1.ReplaceAll("recoAk4PF_Pt[1]>100","(recoAk4PF_Pt[1]>100||(recoAk4PF_Pt[1]>40&&recoAk4PF_Pt[1]<=100))");
 TString offARECO5  = offRECO5;
 offARECO5.ReplaceAll("recoAk4PF_Pt[1]>100","(recoAk4PF_Pt[1]>100||(recoAk4PF_Pt[1]>40&&recoAk4PF_Pt[1]<=100))");

 
 TString GEN2Jet = "(genAk4_NJet40 == 2)";
 TString GEN3Jet = "(genAk4_NJet40 == 3)";
 TString GEN4Jet = "(genAk4_NJet40 >= 4)";

 TString RECO2Jet = "(recoAk4PF_NJet40 == 2)";
 TString RECO3Jet = "(recoAk4PF_NJet40 == 3)";
 TString RECO4Jet = "(recoAk4PF_NJet40 >= 4)";

 TString HLT2Jet = "(hltAk4PF_NJet40 == 2)";
 TString HLT3Jet = "(hltAk4PF_NJet40 == 3)";
 TString HLT4Jet = "(hltAk4PF_NJet40 >= 4)";

 TString AND = "&&";

 // ***********************************************************************************************
 // *                                        Produce plots                                        *
 // ***********************************************************************************************
 std::cout << "Making plots\n";


 cutCollection rateCuts;
 // rateCuts.push_back( std::make_pair( "NoCut", TString("") ) );
 rateCuts.push_back( std::make_pair( "L1WP1", TString(L1+AND+WP1) ) );
 rateCuts.push_back( std::make_pair( "L1WP5", TString(L1+AND+WP5) ) );
 rateCuts.push_back( std::make_pair( "HTLoose40", TString(L1+AND+HTLoose40) ) );
 rateCuts.push_back( std::make_pair( "WPLoose40", TString(L1+AND+WPLoose40) ) );
 rateCuts.push_back( std::make_pair( "WPLoose70", TString(L1+AND+WPLoose70) ) );
 rateCuts.push_back( std::make_pair( "WPLoose80", TString(L1+AND+WPLoose80) ) );
 rateCuts.push_back( std::make_pair( "WPLoose90", TString(L1+AND+WPLoose90) ) );

 cutCollection effCuts;
 // effCuts.push_back( std::make_pair( "NoCut", TString("") ) );
 effCuts.push_back( std::make_pair( "L1Gen1",  TString(L1+AND+offGen1) ) );
 effCuts.push_back( std::make_pair( "L1Gen5",  TString(L1+AND+offGen5) ) );
 effCuts.push_back( std::make_pair( "L1AGen1", TString(L1+AND+offAGen1) ) );
 effCuts.push_back( std::make_pair( "L1AGen5", TString(L1+AND+offAGen5) ) );


 cutCollection effRECOCuts;
 // effCuts.push_back( std::make_pair( "NoCut", TString("") ) );
 effCuts.push_back( std::make_pair( "L1RECO1",  TString(L1+AND+offRECO1) ) );
 effCuts.push_back( std::make_pair( "L1RECO5",  TString(L1+AND+offRECO5) ) );
 effCuts.push_back( std::make_pair( "L1ARECO1", TString(L1+AND+offARECO1) ) );
 effCuts.push_back( std::make_pair( "L1ARECO5", TString(L1+AND+offARECO5) ) );



#ifdef RESOLUTIONS
 // Resolutions
 // -------------------- 
 make2DDistributions( rate,       dMht_vs_mht, "hltAk4PF_MhtPT40/genAk4_MhtPT40:hltAk4PF_MhtPT40",    rateCuts);
 make2DDistributions( background, dMht_vs_mht, "hltAk4PF_MhtPT40/genAk4_MhtPT40:hltAk4PF_MhtPT40",    effCuts);
 make2DDistributions( signal,     dMht_vs_mht, "hltAk4PF_MhtPT40/genAk4_MhtPT40:hltAk4PF_MhtPT40",    effCuts);
 make2DDistributions( signal,     dMht_vs_mht, "hltAk4PF_MhtPT40/recoAk4PF_MhtPT40:hltAk4PF_MhtPT40", effRECOCuts, false, "RECO");

 make2DDistributions( rate,       dMet_vs_met, "hltMetCalo_MetPT/genAk4_MhtPT40:hltMetCalo_MetPT",    rateCuts);
 make2DDistributions( background, dMet_vs_met, "hltMetCalo_MetPT/genAk4_MhtPT40:hltMetCalo_MetPT",    effCuts);
 make2DDistributions( signal,     dMet_vs_met, "hltMetCalo_MetPT/genAk4_MhtPT40:hltMetCalo_MetPT",    effCuts);
 make2DDistributions( signal,     dMet_vs_met, "hltMetCalo_MetPT/recoAk4PF_MhtPT40:hltMetCalo_MetPT", effRECOCuts, false, "RECO");
#endif

 // // MoM vs dPhiMM 
 // // -------------------- 
 // make2DDistributions( rate,  mom_vs_dPhiMM, "Rate_MoM_vs_DPhiMM",  "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", rateCuts, true); 
 // make2DDistributions( eff,   mom_vs_dPhiMM, "Eff_MoM_vs_DPhiMM",   "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", effCuts,  true); 
 // make2DDistributions( ttBar, mom_vs_dPhiMM, "TTBar_MoM_vs_DPhiMM", "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", effCuts,  true);
 // make2DDistributions( DY,    mom_vs_dPhiMM, "DY_MoM_vs_DPhiMM",    "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", effCuts,  true);


 // // AT vs HT
 // // --------------------
 // make2DDistributions( rate,  at_vs_ht, "Rate_HT_vs_AlphaT",  "hltAk4PF_AlphaT40:hltAk4PF_HT40", rateCuts); 
 // make2DDistributions( eff,   at_vs_ht, "Eff_HT_vs_AlphaT",   "hltAk4PF_AlphaT40:hltAk4PF_HT40", effCuts); 
 // make2DDistributions( ttBar, at_vs_ht, "TTBar_HT_vs_AlphaT", "hltAk4PF_AlphaT40:hltAk4PF_HT40", effCuts);
 // make2DDistributions( DY,    at_vs_ht, "DY_HT_vs_AlphaT",    "hltAk4PF_AlphaT40:hltAk4PF_HT40", effCuts);

#ifdef FOR_VS_MOM
 // For vs MoM
 // --------------------
 make2DDistributions( rate,       for_vs_mom, "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", rateCuts,    true);
 make2DDistributions( background, for_vs_mom, "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", effCuts,     true);
 make2DDistributions( signal,     for_vs_mom, "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", effCuts,     true);
 make2DDistributions( signal,     for_vs_mom, "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", effRECOCuts, true, "RECO");
#endif

#ifdef MHT_VS_HT
 // MHT vs HT                                                                                                                                        
 // -------------------- 
 make2DDistributions( rate,       mht_vs_ht, "hltAk4PF_MhtPT40:hltAk4PF_HT40", rateCuts);
 make2DDistributions( background, mht_vs_ht, "hltAk4PF_MhtPT40:hltAk4PF_HT40", effCuts);
 make2DDistributions( signal,     mht_vs_ht, "hltAk4PF_MhtPT40:hltAk4PF_HT40", effCuts);
 make2DDistributions( signal,     mht_vs_ht, "hltAk4PF_MhtPT40:hltAk4PF_HT40", effRECOCuts, false, "RECO");
#endif


 // make2DDistributions( signal,     mht_vs_met, "hltAk4PF_MhtPT40:1.05*hltMetCalo_MetPT", effCuts);

#ifdef MHT_VS_MET
 // MHT vs MET
 // --------------------
 make2DDistributions( rate,       mht_vs_met, "hltAk4PF_MhtPT40:hltMetCalo_MetPT", rateCuts);
 make2DDistributions( background, mht_vs_met, "hltAk4PF_MhtPT40:hltMetCalo_MetPT", effCuts);
 make2DDistributions( signal,     mht_vs_met, "hltAk4PF_MhtPT40:hltMetCalo_MetPT", effCuts);
 make2DDistributions( signal,     mht_vs_met, "hltAk4PF_MhtPT40:hltMetCalo_MetPT", effRECOCuts, false, "RECO");
#endif

#ifdef MHT_VS_MoM
 // MHT vs MoM
 // --------------------
 make2DDistributions( rate,       mht_vs_mom, "hltAk4PF_MhtPT40:hltAk4PF_MhtPT40/hltMetCalo_MetPT", rateCuts);
 make2DDistributions( background, mht_vs_mom, "hltAk4PF_MhtPT40:hltAk4PF_MhtPT40/hltMetCalo_MetPT", effCuts);
 make2DDistributions( signal,     mht_vs_mom, "hltAk4PF_MhtPT40:hltAk4PF_MhtPT40/hltMetCalo_MetPT", effCuts);
 make2DDistributions( signal,     mht_vs_mom, "hltAk4PF_MhtPT40:hltAk4PF_MhtPT40/hltMetCalo_MetPT", effRECOCuts, false, "RECO");
#endif


 // // Jet2 PT vs HT
 // // -------------------- 
 // make2DDistributions( rate,  jet2_vs_ht, "Rate_Jet2_vs_HT",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", rateCuts); 
 // make2DDistributions( eff,   jet2_vs_ht, "Eff_Jet2_vs_HT",   "hltAk4PF_Pt[1]:hltAk4PF_HT40", effCuts); 
 // make2DDistributions( ttBar, jet2_vs_ht, "TTBar_Jet2_vs_HT", "hltAk4PF_Pt[1]:hltAk4PF_HT40", effCuts);
 // make2DDistributions( DY,    jet2_vs_ht, "DY_Jet2_vs_HT",    "hltAk4PF_Pt[1]:hltAk4PF_HT40", effCuts);

 // // For vs eta
 // // -------------------- 
 // make2DDistributions( rate,  for_vs_eta, "Rate_For_vs_Eta",  "hltAk4PFFor_Pt[0]:abs(hltAk4PF_Eta[0])", rateCuts, true); 
 // make2DDistributions( eff,   for_vs_eta, "Eff_For_vs_Eta",   "hltAk4PFFor_Pt[0]:abs(hltAk4PF_Eta[0])", effCuts,  true); 
 // make2DDistributions( ttBar, for_vs_eta, "TTBar_For_vs_Eta", "hltAk4PFFor_Pt[0]:abs(hltAk4PF_Eta[0])", effCuts,  true);
 // make2DDistributions( DY,    for_vs_eta, "DY_For_vs_Eta",    "hltAk4PFFor_Pt[0]:abs(hltAk4PF_Eta[0])", effCuts,  true);

 // // Jet2 PT
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)");
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1_NJet2",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)&&"+HLT2Jet);
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1_NJet3",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)&&"+HLT3Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1_NJet4",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)&&"+HLT4Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)");
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5_NJet2",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)&&"+HLT2Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5_NJet3",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)&&"+HLT3Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5_NJet4",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)&&"+HLT4Jet));


 return 0;


 // // Correlations
 //  make2D( rate, dPhiMM_vs_dPhiMM, "Rate_DPhiMM_vs_DPhiMM_NoCut", "genMetCaloMht40_DeltaPhi:hltMetCaloPFMht40_DeltaPhi", "");
 //  make2D(  eff, dPhiMM_vs_dPhiMM, "Eff_DPhiMM_vs_DPhiMM_NoCut",  "genMetCaloMht40_DeltaPhi:hltMetCaloPFMht40_DeltaPhi", "");

 // // MoM vs dPhiMM
 // // --------------------
 // // make2D( rate, mom_vs_dPhiMM,   "Rate_MoM_vs_DPhiMM_NoCut",       "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // // 	 "",          true);
 // make2D( rate, mom_vs_dPhiMM,   "Rate_MoM_vs_DPhiMM_L1WP1",       "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	 L1+"&&"+WP1, true); 
 // make2D( rate, mom_vs_dPhiMM,   "Rate_MoM_vs_DPhiMM_L1WP5",       "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	 L1+"&&"+WP5, true);

 // // make2D(  eff, mom_vs_dPhiMM,    "Eff_MoM_vs_DPhiMM_NoCut",       "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // // 	  "",              true);
 // make2D(  eff, mom_vs_dPhiMM,    "Eff_MoM_vs_DPhiMM_L1offGen1",   "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	  L1+"&&"+offGen1, true);
 // make2D(  eff, mom_vs_dPhiMM,    "Eff_MoM_vs_DPhiMM_L1offGen5",   "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	  L1+"&&"+offGen5, true);
 
 // // make2D(  ttBar, mom_vs_dPhiMM,  "TTBar_MoM_vs_DPhiMM_NoCut",     "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // // 	  "",              true);
 // make2D(  ttBar, mom_vs_dPhiMM,  "TTBar_MoM_vs_DPhiMM_L1offGen1", "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	  L1+"&&"+offGen1, true);
 // make2D(  ttBar, mom_vs_dPhiMM,  "TTBar_MoM_vs_DPhiMM_L1offGen5", "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	  L1+"&&"+offGen5, true);

 // // make2D(  DY,   mom_vs_dPhiMM,   "DY_MoM_vs_DPhiMM_NoCut",        "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // // 	  "",              true);
 // make2D(  DY,   mom_vs_dPhiMM,   "DY_MoM_vs_DPhiMM_L1offGen1",    "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	  L1+"&&"+offGen1, true);
 // make2D(  DY,   mom_vs_dPhiMM,   "DY_MoM_vs_DPhiMM_L1offGen5",    "hltAk4PF_MhtPT40/hltMetCalo_MetPT:hltMetCaloPFMht40_DeltaPhi", 
 // 	  L1+"&&"+offGen5, true);



 // "hltAk4PFFor_Pt[0]:hltAk4PF_Eta[0]"


 // // // AT vs HT
 // // make2D( rate, at_vs_ht,   "Rate_HT_vs_AlphaT_NoCut",     "hltAk4PF_AlphaT40:hltAk4PF_HT40", "");
 // // make2D( rate, at_vs_ht,   "Rate_HT_vs_AlphaT_L1",        "hltAk4PF_AlphaT40:hltAk4PF_HT40", L1);
 // // make2D(  eff, at_vs_ht,   "Eff_HT_vs_AlphaT_NoCut",      "hltAk4PF_AlphaT40:hltAk4PF_HT40", "");
 // // make2D(  eff, at_vs_ht,   "Eff_HT_vs_AlphaT_L1offGen1",  "hltAk4PF_AlphaT40:hltAk4PF_HT40", L1+"&&"+offGen1);
 // // make2D(  eff, at_vs_ht,   "Eff_HT_vs_AlphaT_L1offGen5",  "hltAk4PF_AlphaT40:hltAk4PF_HT40", L1+"&&"+offGen5);

 // // // For vs MoM
 // // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_NoCut",       "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", "");
 // // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1",          "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1);
 // // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP1",       "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP1);
 // // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP5",       "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP5);

 // // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_NoCut",        "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", "");
 // // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1offGen1",    "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offGen1);
 // // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1offGen5",    "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offGen5);
 // // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP1",        "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP1);
 // // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP5",        "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP5);

 // // // Jet2 PT
 // // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_NoCut",       "hltAk4PF_Pt[1]:hltAk4PF_HT40", "");
 // // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1",          "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1);
 // // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1",       "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)");
 // // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5",       "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)");

 // // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_NoCut",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", "");
 // // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1offGen1",    "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offGen1);
 // // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1offGen5",    "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offGen5);
 // // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP1",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)");
 // // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP5",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)");

 // // For vs MoM
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP1",        "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP1);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP1_NJet2",  "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP1+"&&"+HLT2Jet);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP1_NJet3",  "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP1+"&&"+HLT3Jet);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP1_NJet4",  "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP1+"&&"+HLT4Jet);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP5",        "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP5);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP5_NJet2",  "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP5+"&&"+HLT2Jet);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP5_NJet3",  "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP5+"&&"+HLT3Jet);
 // make2D( rate, for_vs_mom, "Rate_For_vs_MoM_L1WP5_NJet4",  "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+WP5+"&&"+HLT4Jet);


 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP1offRECO1",       "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO1);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP1offRECO1_NJet2", "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO1+"&&"+RECO2Jet);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP1offRECO1_NJet3", "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO1+"&&"+RECO3Jet);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP1offRECO1_NJet4", "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO1+"&&"+RECO4Jet);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP5offRECO5",       "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO5);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP5offRECO5_NJet2", "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO5+"&&"+RECO2Jet);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP5offRECO5_NJet3", "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO5+"&&"+RECO3Jet);
 // make2D(  eff, for_vs_mom, "Eff_For_vs_MoM_L1WP5offRECO5_NJet4", "hltAk4PFFor_Pt[0]:hltAk4PF_MhtPT40/hltMetCalo_MetPT", L1+"&&"+offRECO5+"&&"+RECO4Jet);

 // // Jet2 PT
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)");
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1_NJet2",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)&&"+HLT2Jet);
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1_NJet3",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)&&"+HLT3Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP1_NJet4",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.60)&&"+HLT4Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5",        "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)");
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5_NJet2",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)&&"+HLT2Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5_NJet3",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)&&"+HLT3Jet));
 // make2D( rate, jet2_vs_ht, "Rate_Jet2_vs_HT_L1WP5_NJet4",  "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&(hltAk4PF_AlphaT40>0.51)&&"+HLT4Jet));

 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP1offRECO1",       "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO1);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP1offRECO1_NJet2", "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO1+"&&"+RECO2Jet);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP1offRECO1_NJet3", "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO1+"&&"+RECO3Jet);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP1offRECO1_NJet4", "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO1+"&&"+RECO4Jet);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP5offRECO5",       "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO5);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP5offRECO5_NJet2", "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO5+"&&"+RECO2Jet);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP5offRECO5_NJet3", "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO5+"&&"+RECO3Jet);
 // make2D(  eff, jet2_vs_ht, "Eff_Jet2_vs_HT_L1WP5offRECO5_NJet4", "hltAk4PF_Pt[1]:hltAk4PF_HT40", L1+"&&"+offRECO5+"&&"+RECO4Jet);
 



 delete at_vs_ht;
 delete for_vs_mom;
 delete jet2_vs_ht;
 delete mom_vs_dPhiMM;
 delete dPhiMM_vs_dPhiMM;

 exit(0);
 f->Close();

 // --------------------------------------------------------------------------------


















}
