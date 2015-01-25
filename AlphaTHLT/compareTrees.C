#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TString.h>
#include <TStyle.h>
#include <TCut.h>
#include <TMath.h>
#include <TPad.h>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>

#include <TCanvas.h>
#include <TPaveStats.h>
#include <THStack.h>
#include <TLegend.h>


#include <vector>
#include <iostream>
#include <math.h>

#include "treeFunctions.cpp"


struct corrContainer{

  TString var1;
  TString var2;
  TString cut;
  corrContainer():var1(""),var2(""),cut(""){};
  corrContainer(TString v1, TString v2, TString c):var1(v1),var2(v2),cut(c){};

};


void compareTrees(){


  // ----------------------------------------
  // Input
  // ----------------------------------------
  TString branch    = "MakeTrees/Ntuple";
  TString sampleDir = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/";

  // TString files1 = sampleDir + "01Nov14/QCD_Pt-30to50_Tune4C_13TeV_pythia8/QCD_Pt-30to50_Tune4C_13TeV_pythia8_*.root";
  // TString files2 = sampleDir + "06Nov14_50ns/QCD_Pt-30to50_Tune4C_13TeV_pythia8/QCD_Pt-30to50_Tune4C_13TeV_pythia8_*.root";


  TString QCD30to50    = "QCD_Pt-30to50_Tune4C_13TeV_pythia8/QCD_Pt-30to50_Tune4C_13TeV_pythia8_*.root";
  TString QCD50to80    = "QCD_Pt-50to80_Tune4C_13TeV_pythia8/QCD_Pt-50to80_Tune4C_13TeV_pythia8_*.root";
  TString QCD80to120   = "QCD_Pt-80to120_Tune4C_13TeV_pythia8/QCD_Pt-80to120_Tune4C_13TeV_pythia8_*.root";
  TString QCD120to170  = "QCD_Pt-120to170_Tune4C_13TeV_pythia8/QCD_Pt-120to170_Tune4C_13TeV_pythia8_*.root";
  TString QCD170to300  = "QCD_Pt-170to300_Tune4C_13TeV_pythia8/QCD_Pt-170to300_Tune4C_13TeV_pythia8_*.root";
  TString QCD300to470  = "QCD_Pt-300to470_Tune4C_13TeV_pythia8/QCD_Pt-300to470_Tune4C_13TeV_pythia8_*.root";
  TString QCD470to600  = "QCD_Pt-470to600_Tune4C_13TeV_pythia8/QCD_Pt-470to600_Tune4C_13TeV_pythia8_*.root";
  TString QCD600to800  = "QCD_Pt-600to800_Tune4C_13TeV_pythia8/QCD_Pt-600to800_Tune4C_13TeV_pythia8_*.root";
  TString QCD800to1000 = "QCD_Pt-800to1000_Tune4C_13TeV_pythia8/QCD_Pt-800to1000_Tune4C_13TeV_pythia8_*.root";

  TString selectedSample = QCD600to800;

  TString preFix25ns  = sampleDir + "01Nov14/";  
  TString preFix50ns  = sampleDir + "06Nov14_50ns/";  
  TString hcalFix25ns = sampleDir + "12Nov14/";
  TString hcalFix50ns = sampleDir + "12Nov14_50ns/";



  // TString files1 = sampleDir + "01Nov14/QCD_Pt-30to50_Tune4C_13TeV_pythia8/QCD_Pt-30to50_Tune4C_13TeV_pythia8_1.root";
  // TString files2 = sampleDir + "06Nov14_50ns/QCD_Pt-30to50_Tune4C_13TeV_pythia8/QCD_Pt-30to50_Tune4C_13TeV_pythia8_1.root";
  // TString label1 = "25ns";
  // TString label2 = "50ns";

  TString files1 = hcalFix25ns + selectedSample;
  TString files2 = hcalFix50ns + selectedSample;

  TString label1 = "25ns"; 
  TString label2 = "50ns";


  TString compBaseTitle    = label1 + TString(" vs ") + label2;
  TString compBaseFilename = label1 + TString("_vs_") + label2 + TString("_");


  std::vector<corrContainer> correlationPlots;

  correlationPlots.push_back( corrContainer("hltAk4Calo_Phi",           "hltAk4Calo_Eta",               "hltAk4Calo_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4Calo_Pt",            "hltAk4Calo_Eta",               "hltAk4Calo_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4Calo_HT40",          "hltAk4Calo_AlphaT40",          "") );
  correlationPlots.push_back( corrContainer("hltAk4CaloNoFastJet_Phi",  "hltAk4CaloNoFastJet_Eta",      "hltAk4CaloNoFastJet_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4CaloNoFastJet_Pt",  "hltAk4CaloNoFastJet_Eta",       "hltAk4CaloNoFastJet_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4PF_Phi",             "hltAk4PF_Eta",                 "hltAk4PF_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4PF_Pt",              "hltAk4PF_Eta",                 "hltAk4PF_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4PF_AlphaT40",        "hltAk4PF_HT40",                "") );
  correlationPlots.push_back( corrContainer("genAk4_Phi",               "genAk4_Eta",                   "genAk4_Pt > 40") );
  correlationPlots.push_back( corrContainer("genAk4_Pt",                "genAk4_Eta",                   "genAk4_Pt > 40") );
  correlationPlots.push_back( corrContainer("genAk4_AlphaT40",          "genAk4_HT40",                  "") );


  correlationPlots.push_back( corrContainer("hltAk4Calo_Pt",           "genAk4_Pt",  "hltAk4Calo_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4Calo_Eta",          "genAk4_Eta", "hltAk4Calo_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4Calo_Phi",          "genAk4_Eta", "hltAk4Calo_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4CaloNoFastJet_Pt",  "genAk4_Pt",  "hltAk4CaloNoFastJet_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4CaloNoFastJet_Eta", "genAk4_Eta", "hltAk4CaloNoFastJet_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4CaloNoFastJet_Phi", "genAk4_Eta", "hltAk4CaloNoFastJet_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4PF_Pt",             "genAk4_Pt",  "hltAk4PF_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4PF_Eta",            "genAk4_Eta", "hltAk4PF_Pt > 40") );
  correlationPlots.push_back( corrContainer("hltAk4PF_Phi",            "genAk4_Eta", "hltAk4PF_Pt > 40") );

  correlationPlots.push_back( corrContainer("hltAk4Calo_HT40",         "genAk4_HT40",     "") );
  correlationPlots.push_back( corrContainer("hltAk4Calo_AlphaT40",     "genAk4_AlphaT40", "") );
  correlationPlots.push_back( corrContainer("hltAk4PF_HT40",           "genAk4_HT40",     "") );
  correlationPlots.push_back( corrContainer("hltAk4PF_AlphaT40",       "genAk4_AlphaT40", "") );

  correlationPlots.push_back( corrContainer("hltMhtCalo_MetPT",        "genMetCalo_MetPt", "") );
  correlationPlots.push_back( corrContainer("hltMhtPF_MetPT",          "genMetCalo_MetPt", "") );




  // ----------------------------------------
  // Output
  // ----------------------------------------
  TString plotDirectory = "/home/hep/mb1512/SUSY/UCTHLT/CMSSW_7_2_1_patch2/src/AlphaTHLT/plotTest/";


  // Load the tchains
  TChain *chain1 = new TChain( branch, "Chain1" );
  chain1->Add( files1 );
  TChain *chain2 = new TChain( branch, "Chain2" );
  chain2->Add( files2 );
  if ( (chain1 == NULL) || (chain2 == NULL) ){
    std::cout << "Error, problem with input files.\n";
    return;
  }

  // Loop through the branches
  std::vector<TString> branchNames = getBranchNames( chain1 );

  for (uint iBranch = 0; iBranch < branchNames.size(); ++iBranch ){
    TString branchName = branchNames[ iBranch ];

    std::cout << "Branch: " << branchName << "\n";
    TString compTitle    = compBaseTitle    + TString(" - ") + branchName;
    TString compFilename = compBaseFilename + branchName;
    TString h1Name       =  branchName + "1";
    TString h2Name       =  branchName + "2";

    chain1->Draw( branchName );
    TH1* h1 = (TH1*)chain1->GetHistogram()->Clone();
    h1->SetName( h1Name );
    
    // Set up the same binning as the first histogram
    TH1F* h2 = (TH1F*)h1->Clone();
    h2->SetName( h2Name );
    chain2->Project( h2Name, branchName );

    if ( (h1 == NULL) || (h2 == NULL) ){ 
      std::cout << "Error, branch '" << branchName << "' not present in both files.\n" 
    		<< "\tFile1 = " << (h1 == NULL) << " File2 = " << (h2 == NULL) << "\n";
      continue;
    }

    // set line colours 
    h1->SetLineColor(kRed);
    h2->SetLineColor(kBlue);

    // overlay the plots 
    TCanvas* canv      = new TCanvas();
    THStack* histStack = new THStack("hs", compTitle);
    histStack->Add( h1 );
    histStack->Add( h2 );
    histStack->Draw("ehist nostack");
    canv->Update(); // Generate stats boxes

    // Statsbox wizardry 
    // TPaveStats *stat1 = (TPaveStats*)h1->GetListOfFunctions()->FindObject("stats");
    // TPaveStats *stat2 = (TPaveStats*)h2->GetListOfFunctions()->FindObject("stats");
    // // stat1->SetTextColor(kRed);
    // // stat2->SetTextColor(kBlue);
    // // stat1->SetX1NDC(.6); stat1->SetX2NDC(.8);  stat1->SetY1NDC(.8); stat1->SetY2NDC(.99);
    // // stat2->SetX1NDC(.8); stat2->SetX2NDC(.99); stat2->SetY1NDC(.8); stat2->SetY2NDC(.99);
    // canv->Modified();


    // Add a legend 
    TLegend* Legend = new TLegend(0.7, 0.7,0.84, 0.8);
    Legend->SetFillColor(kWhite);
    Legend->AddEntry( h1, label1, "l");
    Legend->AddEntry( h2, label2, "l");
    Legend->Draw();

    // Draw the statsbox on top 
    // stat1->Draw();
    // stat2->Draw();

    canv->SaveAs(plotDirectory + compFilename + ".png");
    canv->SaveAs(plotDirectory + compFilename + ".pdf");




    // ************************************************************ 
    // *                       Ratio plots                        * 
    // ************************************************************ 

    TH1* tempHist       = (TH1*)h1->Clone();
    TString titleRatio  = label1 + TString(" over ") + label2;
    tempHist->Sumw2();
    tempHist->SetTitle(titleRatio);
    tempHist->GetYaxis()->SetTitle("Ratio");
    tempHist->SetLineColor( kBlack );
    tempHist->SetMarkerStyle(20);
    tempHist->SetMarkerSize(0.5);
    tempHist->Divide( h2 );
    tempHist->SetStats(0);
    tempHist->Draw("pe");
    canv->SaveAs(plotDirectory + compFilename + "_Ratio.png");
    canv->SaveAs(plotDirectory + compFilename + "_Ratio.pdf");


    // ************************************************************************************************** 
    // *                                           Log plots                                            * 
    // ************************************************************************************************** 

    canv->SetLogy();
    histStack->Draw("nostack");

    canv->SaveAs(plotDirectory + compFilename + "_Log.png");
    canv->SaveAs(plotDirectory + compFilename + "_Log.pdf");



    // DEBUGGING
    //break;

  }


  // DEBUGGING
  exit(0);
  // DEBUGGING


    // ************************************************************ 
    // *                    Correlation plots                     * 
    // ************************************************************ 

    // Produce defined correlations
    for (uint iCorr = 0; iCorr < correlationPlots.size(); ++iCorr){
      corrContainer correlation = correlationPlots[ iCorr ];
      TString hist2D1Name     = label1 + TString("_") + correlation.var1 + TString("_") + correlation.var2;
      TString hist2D2Name     = label2 + TString("_") + correlation.var1 + TString("_") + correlation.var2;
      TString hist2DRatioName = TString("Ratio_") + correlation.var1 + TString("_") + correlation.var2;
      std::cout << correlation.var1 << "\t" << correlation.var2 << "\t" << correlation.cut << "\n";
      TCanvas* canv      = new TCanvas();
      TH2* hist2D1;
      TH2* hist2D2;
      TH2* hist2DRatio;

      // File 1
      // --------------------
      canv->SetLogz(0);

      chain1->Draw( correlation.var1 + ":" + correlation.var2, correlation.cut, "COLZ" );
      hist2D1     = (TH2*)chain1->GetHistogram()->Clone();
      hist2D2     = (TH2*)hist2D1->Clone();
      hist2DRatio = (TH2*)hist2D1->Clone();

      hist2D1->SetName( hist2D1Name );
      hist2D1->SetTitle( label1 + " " + hist2D1->GetTitle() );
      hist2D1->Draw("COLZ");

      canv->SaveAs(plotDirectory + hist2D1Name + "_Corr.png");
      canv->SaveAs(plotDirectory + hist2D1Name + "_Corr.pdf");

      canv->SetLogz(1);
      canv->SaveAs(plotDirectory + hist2D1Name + "_Corr_Log.png");
      canv->SaveAs(plotDirectory + hist2D1Name + "_Corr_Log.pdf");

      // File 2
      // --------------------
      canv->SetLogz(0);

      chain2->Project( hist2D2Name, correlation.var1 + ":" + correlation.var2, correlation.cut, "COLZ" );
      hist2D2->SetName(  hist2D2Name );
      hist2D2->SetTitle( label2 + " " + hist2D2->GetTitle() );
      hist2D2->Draw("COLZ");
      //      chain2->Draw( correlation.var1 + ":" + correlation.var2, correlation.cut, "COLZ" );

      canv->SaveAs(plotDirectory + hist2D2Name + "_Corr.png");
      canv->SaveAs(plotDirectory + hist2D2Name + "_Corr.pdf");

      canv->SetLogz(1);
      canv->SaveAs(plotDirectory + hist2D2Name + "_Corr_Log.png");
      canv->SaveAs(plotDirectory + hist2D2Name + "_Corr_Log.pdf");


      // Ratio
      // --------------------
      canv->SetLogz(0);
      hist2DRatio->SetName( hist2DRatioName );
      hist2DRatio->SetTitle( TString("Ratio ") + hist2DRatio->GetTitle() );
      hist2DRatio->Divide( hist2D2 );
      hist2DRatio->Draw("COLZ");
      canv->SaveAs(plotDirectory + hist2DRatioName + "_Corr_Ratio.png");
      canv->SaveAs(plotDirectory + hist2DRatioName + "_Corr_Ratio.pdf");


    }



}


