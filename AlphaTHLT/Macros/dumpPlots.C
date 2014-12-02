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
#include <TGraph.h>
#include <TH1.h>
#include <TH2.h>
#include <TEfficiency.h>
#include <TSpline.h>
#include <TLine.h>

#include <vector>
#include <iostream>
#include <math.h>


bool plotExists(TFile *f, TString rootDir, TString plotName);
void getPFWorkingPoint( TString effHistName, float &PFHT, float &PFAlphaT );
void reverseCumulative2D( TH2* histogram, TH2* rCumulHist );
void makePlot( TFile *f, TString outputDir, TString rootDir, TString plotName );
void makeCumulativePlot( TFile *f, TString outputDir, TString rootDir, TString plotName );
void makeJet12RatePlot( TFile *f, TString outputDir, TString rootDir, TString plotName );
std::vector< std::pair< double, double> > getEfficiencyBins( TH2* histogram, double maxRate );
void makeContourEffPlot( TFile *fRate, TFile *fEff, TString outputDir, 
                         TString rawRateDir, TString rateHistName, 
                         TString rawEffDir,  TString effHistName, 
                         std::vector< std::pair<float, float> > maxRatePair, 
                         TString extra, bool th2, 
                         double lineX, double lineY );





// Flags
#define PREFILTER
#define PFTRIGGER


void dumpPlots(){

  gStyle->SetOptStat(0);
  gStyle->SetTitleStyle(0);
  gStyle->SetPaintTextFormat("9.2f");

  // ***********************************************************************************************
  // *                                        Configuration                                        *
  // ***********************************************************************************************

  TString directory = "/home/mark/Dropbox/Ntuples/28_11_14/Samples/";

  // Samples
  // ----------
  // "T1tttt_2J_mGl_1200_mLSP_800_0";  "T2tt_2J_mStop_850_mLSP_100_0"; "T1bbbb_2J_mGl_1000_mLSP_900_0.root";
  // "T2tt_2J_mStop_425_mLSP_325_0";
  TString sample    = "DY";
  //TString sample    = "TTBar";


  TString rateFilename   = directory + "QCD.root";
  TString effFilename    = directory + sample + ".root";
  TString outputDir      = directory + sample + "/"; 
  TString outputFilename = directory + sample + "_Output.root";
  TString rateOutputDir = "Rate/"; 

  TString rawRateDir   = "Raw/Rate/";
  TString rawEffDir    = "Raw/PrefilterEfficiencyRaw/";
  TString rawPFEffDir  = "Raw/Efficiency/";
  
  std::vector< std::pair<float, float> > maxPreFilterRatePairs;
  maxPreFilterRatePairs.push_back( std::make_pair( 0.5, 250) );
  maxPreFilterRatePairs.push_back( std::make_pair( 250, 500) );
  maxPreFilterRatePairs.push_back( std::make_pair( 500, 1000) );
  
  std::vector< std::pair<float, float> > maxPFRatePairs;
  maxPFRatePairs.push_back( std::make_pair( 0.5, 10) );
  maxPFRatePairs.push_back( std::make_pair( 10, 20) );
  maxPFRatePairs.push_back( std::make_pair( 20, 30) );

  std::vector<TString> jetBins;
  jetBins.push_back("Inclusive");
  jetBins.push_back("2Jets");
  jetBins.push_back("3Jets");
  jetBins.push_back("4Jets");
  jetBins.push_back("5Jets");
  jetBins.push_back("6Jets");

  std::vector<TString> l1Triggers;
  l1Triggers.push_back("Trig1");
  l1Triggers.push_back("Trig2");

  std::vector<TString> alphaTTypes;
  alphaTTypes.push_back("Std");
  alphaTTypes.push_back("Dyn");
  alphaTTypes.push_back("MoH");
  //  alphaTTypes.push_back("NFJ");
  alphaTTypes.push_back("Pri");

  std::vector<TString> alphaTThresholds;
  alphaTThresholds.push_back("gt0p60");
  alphaTThresholds.push_back("gt0p65");
  alphaTThresholds.push_back("gt0p70");
  alphaTThresholds.push_back("gt0p75");
  alphaTThresholds.push_back("gt0p80");

  std::vector<TString> jet2PTThresholds;
  jet2PTThresholds.push_back("gt0");
  // jet2PTThresholds.push_back("gt50");
  jet2PTThresholds.push_back("gt60");
  jet2PTThresholds.push_back("gt70");
  jet2PTThresholds.push_back("gt80");
  jet2PTThresholds.push_back("gt90");
  jet2PTThresholds.push_back("gt100");

  std::vector<TString> anaJetBins;
  anaJetBins.push_back("2to3");
  anaJetBins.push_back("3to4");
  anaJetBins.push_back("4to9");
  
  std::vector<TString> anaHtBins;
  anaHtBins.push_back("HT200to300");
  anaHtBins.push_back("HT300to400");
  anaHtBins.push_back("HT400to500");
  anaHtBins.push_back("HT500to9999");
  // anaHtBins.push_back("HT200to275");
  // anaHtBins.push_back("HT275to325");
  // anaHtBins.push_back("HT325to375");
  // anaHtBins.push_back("HT375to9999");

  TFile *fRate = new TFile( rateFilename,   "OPEN");
  TFile *fEff  = new TFile( effFilename,    "OPEN");
  TFile *o     = new TFile( outputFilename, "RECREATE" );
  std::cout << "Opening file: " << effFilename << "\nWriting results to file: " << outputFilename << "\n\n";

  // --------------------------------------------------------------------------------


  // ********************************************************************************
  // *                                 Level-1 Loop                                 *
  // ********************************************************************************
  for (uint iL1 = 0; iL1 < l1Triggers.size(); ++iL1){
    TString l1Trig    = l1Triggers[ iL1 ];
    TString l1TrigEff = TString("L1") + Form("%d", iL1 + 1);


#ifdef PREFILTER
    // ****************************************************************************************************
    // *                                          Calo Prefilter                                          *
    // ****************************************************************************************************

    for (uint iType = 0;iType < alphaTTypes.size();++iType){
      TString aTType = alphaTTypes[iType];

      for (uint iJet2PT = 0; iJet2PT < jet2PTThresholds.size(); ++iJet2PT){
	TString jet2PTStr = jet2PTThresholds[ iJet2PT ];


	TString rateCaloHistName = l1Trig + "CaloOnRate_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_Inclusive_Cumul";
	TString ratePFHistName   = l1Trig + "PFOnRate_AlphaT"   + aTType + "_vs_HT_Jet2" + jet2PTStr + "_Inclusive_Cumul";

	// ------------------------------------------------------------
	// Sanity check, plot rate bands on rate plot
	// ------------------------------------------------------------
	// Calo
	makeContourEffPlot( fRate, fRate, outputDir, rawRateDir, rateCaloHistName, rawRateDir, rateCaloHistName, maxPreFilterRatePairs, jet2PTStr, true, -1, -1 );

	// PF
	if (aTType == "Std"){
	  makeContourEffPlot( fRate, fRate, outputDir, rawRateDir, ratePFHistName, rawRateDir, ratePFHistName, maxPFRatePairs,  jet2PTStr, true, -1, -1 );
	}

	//      continue;

	// ------------------------------------------------------------
	// Individual calo AlphaT HTxAlphaT triggers
	// ------------------------------------------------------------
	for (uint iJet = 0;iJet < jetBins.size(); ++iJet){
	  TString jetBin = jetBins[iJet];

      	  if ( ( jetBin == "2Jets" ) || ( jetBin == "3Jets" ) || ( jetBin == "4Jets" ) ){ 
	
      	    // Skip irrelevent bins
      	    if ( (sample.Contains("TTBar")) &&  ( jetBin == "2Jets" ) ){ continue; }
      	    if ( (sample.Contains("DY"))    && !( jetBin == "2Jets" ) ){ continue; }

	    // Get PF trigger working points
	    float PFHT(0), PFAlphaT(0);


      	    TString effHistName  = "";

	    for (uint l1On = 0; l1On <= 1; ++l1On ){
	      TString l1RequiredStr = "";
	      if (l1On == 1){
		l1RequiredStr = "_NoReqL1";
	      }
	      
	      for (uint hltN = 1; hltN <= 5; ++hltN ){
		TString hltStr = TString("HLT") + Form("%d", hltN);

		effHistName  = l1TrigEff + "_" + hltStr + "On_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_" + jetBin + 
		               l1RequiredStr + "_UniCumul";
		getPFWorkingPoint( effHistName, PFHT, PFAlphaT );
		std::cout << effHistName << "\n";
		makeContourEffPlot( fRate, fEff, outputDir, rawRateDir, rateCaloHistName, rawEffDir, effHistName, maxPreFilterRatePairs,  jet2PTStr, false, PFHT, PFAlphaT );

	      } // End HLT loop
	    } // End L1 loop
    
      	  }


	} // End njet loop

	// ------------------------------------------------------------
	// Analysis binned plots
	// ------------------------------------------------------------
	for (uint iJet = 0;iJet < anaJetBins.size(); ++iJet){
	  TString jetBin = anaJetBins[ iJet ];

	  // Restrict to Static AlphaT
	  if (aTType != "Std"){ continue; }

	  // Skip irrelevent bins
	  if ( (sample.Contains("TTBar")) &&  ( jetBin == "2to3" ) ){ continue; }
	  if ( (sample.Contains("DY"))    && !( jetBin == "2to3" ) ){ continue; }

	  for (uint iHT = 0;iHT < anaHtBins.size(); ++iHT){
	    TString htBin = anaHtBins[ iHT ];
  
	    TString effAnaPFHistName  = "Trig1PFOn_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_" + htBin + "_" + jetBin + "_UniCumul";
	    std::cout << ratePFHistName << "\n" << effAnaPFHistName << "\n";
	    makeContourEffPlot( fRate, fEff, outputDir, rawRateDir, ratePFHistName, rawPFEffDir, effAnaPFHistName, maxPFRatePairs, jet2PTStr, false, -1, -1 );

	  }
	} // End analysis binned


      }
    } // End AlphaT types


#endif


  } // End L1 triggers


  exit(0);

  // // ************************************************************
  // // Efficiency plots
  // // ************************************************************
  // for (uint iType = 0;iType < alphaTTypes.size();++iType){
    
  //   TString aTType = alphaTTypes[iType];
  //   for (uint iJet = 2;iJet < 7;++iJet){
      
  //     TString jetMult = Form("%d", iJet);
  //     TString plotName    = "AlphaT" + aTType + "_vs_HT_" + jetMult + "Jets";
  //     TString plotNameCut = "AlphaT" + aTType + "_vs_HT_Jet2gt100_" + jetMult + "Jets";


  //     makePlot(           f, outputDir, effDir, plotName );
  //     makeCumulativePlot( f, outputDir, "Raw/Efficiency/", plotName );
      
      
  //     makePlot(           f, outputDir, effDir, plotNameCut );
  //     makeCumulativePlot( f, outputDir, "Raw/Efficiency/", plotNameCut );

  //   }

  // }

  // // ************************************************************
  // // Rate plots
  // // ************************************************************
  // // "Jet2PT_vs_Jet1PT_alphaTDyngt0p65_4Jets_Rate"


  // for (uint iType = 0;iType < alphaTTypes.size();++iType){
  //   TString aTType = alphaTTypes[iType];
  //   if ( aTType == "Std" ){ aTType = ""; }
  //   else if (aTType.Contains("Dy3")){ continue; }
  //   else if ( !(aTType.Contains("Dyn")) ){ aTType.ReplaceAll("Dy","Dyn"); }

  //   for (uint iThresh = 0;iThresh < alphaTThresholds.size();++iThresh){
  //     TString aTThresh = alphaTThresholds[ iThresh ];
  //     TString plotName = "Jet2PT_vs_Jet1PT_alphaT" + aTType + aTThresh + "_Rate";
    
  //     std::cout <<  plotName << "\n";
  //     makeJet12RatePlot( f, rateOutputDir, rateDir,  plotName );
    
  //   }
  // }


}
















bool plotExists(TFile *f, TString rootDir, TString plotName){
  return f->Get( rootDir + plotName);
}

// Extract the working point of the PF HLT trigger - Needs to be updated as it changes
void getPFWorkingPoint( TString effHistName, float &PFHT, float &PFAlphaT ){
  PFHT     = 0;
  PFAlphaT = 0;

  if      ( effHistName.Contains("HLT1") ){ PFHT = 200; PFAlphaT = 0.57; }
  else if ( effHistName.Contains("HLT2") ){ PFHT = 250; PFAlphaT = 0.55; }
  else if ( effHistName.Contains("HLT3") ){ PFHT = 300; PFAlphaT = 0.53; }
  else if ( effHistName.Contains("HLT4") ){ PFHT = 350; PFAlphaT = 0.52; }
  else if ( effHistName.Contains("HLT5") ){ PFHT = 400; PFAlphaT = 0.51; }

  // Convert AlphaT to MoH, approximating AlphaTPrime = AlphaT
  if ( effHistName.Contains("MoH") ){
    float MoH = sqrt( 1 - 0.25/(PFAlphaT*PFAlphaT));
      PFAlphaT = MoH;
  }
  return;
}



void reverseCumulative2D( TH2* histogram, TH2* rCumulHist ){

  if ( histogram->GetEntries() == 0 ){
    return;
  }
 
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
}

void makePlot( TFile *f, TString outputDir, TString rootDir, TString plotName ){
  if ( !(plotExists( f, rootDir, plotName ) ) ){ return; }

  f->Get( rootDir + plotName )->SaveAs( outputDir + plotName + ".pdf");
}
void makeCumulativePlot( TFile *f, TString outputDir, TString rootDir, TString plotName ){
  if ( !(plotExists( f, rootDir, plotName ) ) ){ return; }

  TEfficiency *eff      = (TEfficiency*)f->Get( rootDir + plotName )->Clone();
  TH2* passed           = (TH2*)eff->GetPassedHistogram()->Clone();
  TH2* total            = (TH2*)eff->GetTotalHistogram() ->Clone();

  TEfficiency *effCumul = (TEfficiency*)eff->Clone();
  TH2* passedCumul      = (TH2*)effCumul->GetPassedHistogram()->Clone();
  TH2* totalCumul       = (TH2*)effCumul->GetTotalHistogram() ->Clone();

  reverseCumulative2D( total,  totalCumul );
  reverseCumulative2D( passed, passedCumul );

  // Have to set the TEfficiency in this order, otherwise it will reject the passed histogram and still set the new total histogram
  effCumul->SetTotalHistogram(  *totalCumul,  ""  );
  effCumul->SetPassedHistogram( *passedCumul, "" );

  TCanvas *canv    = new TCanvas();
  effCumul->Draw("COLZTEXT");
  canv->SaveAs(outputDir + plotName + "_Cumul.pdf");

  delete canv;
  delete eff;
  delete effCumul;
  delete passedCumul;
  delete totalCumul;
}

void makeJet12RatePlot( TFile *f, TString outputDir, TString rootDir, TString plotName ){

  TCanvas *canv = (TCanvas*)f->Get(rootDir + plotName);
  TH2D* hist    = (TH2D*)canv->GetPrimitive(plotName + "_Cumul" );
  canv->Clear();
  canv->Draw();

  hist->GetXaxis()->SetRangeUser(95, 150);
  hist->GetYaxis()->SetRangeUser(50, 150);
  hist->Draw("COLZTEXT");
  canv->SaveAs( outputDir + plotName + ".pdf");

}



// Return the bins for seed1 vs seed2 correlations for which the rate is below a given threshold 
std::vector< std::pair< double, double> > getEfficiencyBins( TH2* histogram, double maxRate ){
  
  // Points to return
  std::vector< std::pair< double, double> > selectedValues; 
  
  // Find bin with minimum cut for acceptable rate 
  int nBinsX          = histogram->GetNbinsX();
  int nBinsY          = histogram->GetNbinsY();
  
  // Scan each y column for constant x to find minimum threshold that is below cut 
  for ( int iBinX = 1; iBinX <= nBinsX; ++iBinX ){
    for ( int iBinY = 1; iBinY <= nBinsY; ++iBinY ){
      float currRate = histogram->GetBinContent( iBinX, iBinY );
      //            std::cout << currRate << "\t" << iBinX << "\t" << iBinY << "\n";

      if ( (iBinY == 1) && (currRate == 0) ){ break; } // No data

      // std::cout << "( " << iBinX << ", " << iBinY << ") = " << currRate << "\t" 
      // 		<< "( " << histogram->GetXaxis()->GetBinCenter( iBinX ) 
      // 		<< ", " << histogram->GetYaxis()->GetBinCenter( iBinY ) << ")\n";

      // Find the first bin that is below the acceptable rate
      if ( currRate < maxRate ){                             

	// Find which bin is closer to the acceptable rate
	if ( iBinY > 1 ){
	  double rateAbove = histogram->GetBinContent( iBinX, iBinY + 1);
          if ( fabs(rateAbove - maxRate) < fabs(currRate - maxRate) ){
	    iBinY = iBinY + 1; // The true rate transition point occurs in the bin above
	  }
	}
                                                   
	//	selectedBins.push_back( std::make_pair(iBinX, iBinY) );

	double xValue = histogram->GetXaxis()->GetBinLowEdge( iBinX );
	double yValue = histogram->GetYaxis()->GetBinLowEdge( iBinY );


	selectedValues.push_back( std::make_pair( xValue, yValue ) );
	break;
      }
 
    }
  }


  //  return selectedBins;
  return selectedValues;
}





  

void makeContourEffPlot( TFile *fRate, TFile *fEff, TString outputDir, 
			 TString rawRateDir, TString rateHistName,
			 TString rawEffDir,  TString effHistName, 
			 std::vector< std::pair<float, float> > maxRatePair, 
			 TString extra, bool th2,
			 double lineX, double lineY ){

  double maxYValue = 1.;


  TCanvas *c = new TCanvas(effHistName);
  c->SetGridx();
  c->SetGridy();
  TLegend *legend = new TLegend( 0.7,0.65, 0.85, 0.85);
  legend->SetFillColor( 0 );

  TH2D* histogramRate       = (TH2D*)fRate->Get( rawRateDir + rateHistName );
  if ( histogramRate == NULL ){ 
    std::cout << "\nError, could not find histogram: " << rawRateDir << rateHistName << "\n"; 
    exit(0);
  }

  TEfficiency* histogramEff;
  TH2D* histogramEff2D; 

  // Determine plot range
  float xLow(0), xHigh(0), yLow(0), yHigh(0), maximum(0);
  xLow  = 100;
  xHigh = 450;
  if ( !effHistName.Contains("MoH") ){
    yLow  = 0.49;
    yHigh = 0.6;
    if (th2){ yLow  = 0.5; yHigh = 0.75; }
  }
  else{ yLow  = 0.0; yHigh = 1.0;  }

  if ( rateHistName.Contains("CaloOnRate") ){ maximum = 1500.; }
  else                                      { maximum = 100.;  }

  // Draw histogram to be overlayed over
  if (th2){
    histogramEff2D = (TH2D*)fEff->Get( rawEffDir  + effHistName );
    if ( histogramEff2D == NULL ){ 
      std::cout << "\nError, could not find histogram: " << rawEffDir << effHistName << "\n"; 
      exit(0);
    }
    histogramEff2D->Draw("COLZTEXT"); 
    gPad->Update();

    // Set ranges
    histogramEff2D->GetXaxis()->SetRangeUser( xLow, xHigh );
    histogramEff2D->GetYaxis()->SetRangeUser( yLow, yHigh );
    histogramEff2D->SetMaximum( maximum );
    
    histogramEff2D->Draw("COLZTEXT");

  }
  else{
    histogramEff = (TEfficiency*)fEff->Get( rawEffDir  + effHistName );
    if ( histogramEff == NULL ){ 
      std::cout << "\nError, could not find histogram: " << rawEffDir << effHistName << "\n"; 
      exit(0);
    }
    histogramEff->Draw("COLZTEXT"); 
    gPad->Update();

    // Set ranges
    histogramEff->GetPaintedHistogram()->GetXaxis()->SetRangeUser( xLow, xHigh );
    if ( !effHistName.Contains("MoH") ){
      histogramEff->GetPaintedHistogram()->GetYaxis()->SetRangeUser( yLow, yHigh );
    }
    histogramEff->GetPaintedHistogram()->SetMinimum(0.);
    histogramEff->GetPaintedHistogram()->SetMaximum(1.);
    histogramEff->GetPaintedHistogram()->Draw("COLZTEXT");

  }


  // ********************************************************************************
  // Draw isocontours
  // ********************************************************************************
  for (uint iRate = 0; iRate < maxRatePair.size(); ++iRate ){
    
    // Get points bordering allowed rate
    float maxRateLow  = maxRatePair[ iRate ].first;
    float maxRateHigh = maxRatePair[ iRate ].second;
    std::vector< std::pair< double, double> > ratePairsMin = getEfficiencyBins( histogramRate, maxRateLow  );
    std::vector< std::pair< double, double> > ratePairsMax = getEfficiencyBins( histogramRate, maxRateHigh );

    int nPairsMin = ratePairsMin.size();
    int nPairsMax = ratePairsMax.size();
    int nPairsTotal = nPairsMin + nPairsMax + 2;

    //    std::cout << "MinPairs = " << nPairsMin << "\t" << nPairsMax << "\n";

    if ( nPairsMin > 0 ){
      double *xx = new double[nPairsTotal];  
      double *yy = new double[nPairsTotal];
      for (uint iPair = 0;iPair < ratePairsMin.size(); ++iPair ){

  	if (iPair == 0){
	  xx[ 0 ] = ratePairsMin[ iPair ].first;
	  yy[ 0 ] = maxYValue;
	}

	xx[ iPair+1 ] = ratePairsMin[ iPair ].first;
	yy[ iPair+1 ] = ratePairsMin[ iPair ].second;
      }

      for (int iPair = ratePairsMax.size() -1;iPair >-1; --iPair ){

	int index = nPairsTotal - 2 - iPair;

	xx[ index ] = ratePairsMax[ iPair ].first;
	yy[ index ] = ratePairsMax[ iPair ].second;


  	if (iPair == 0 ){
	  xx[ index + 1] = ratePairsMax[ iPair ].first;
	  yy[ index + 1] = maxYValue;
	}

      }


      TString maxRateStr = Form("%d", int(maxRateLow)) + TString(" < Rate < ") + Form("%d", int(maxRateHigh)) + TString(" Hz");
      
      int colour = iRate + 1;
      if (colour == 3){ colour = kYellow; }

      TGraph* graph       = new TGraph( nPairsTotal, xx, yy);
      TGraph* graphHollow = new TGraph( nPairsTotal, xx, yy);
      graph->SetFillColor( colour );
      graph->SetFillStyle( 3001 );
      graphHollow->SetFillColor( colour );
      graphHollow->SetFillStyle( 0 );

      // Draw fill
      graph      ->Draw("FL");
      // Draw hollow outline
      graphHollow->Draw("FL");

      legend ->AddEntry( graph, maxRateStr);
      delete[] xx;
      delete[] yy;

    }
  }

  // Draw lines
  if (!th2){
    if (lineX != -1){
      TLine *xLine = new TLine( lineX, yLow, lineX, lineY ); xLine->SetLineStyle(7); xLine->SetLineWidth(3); xLine->Draw();
    }
    if (lineY != -1){
      TLine *yLine = new TLine( xLow, lineY, lineX, lineY ); yLine->SetLineStyle(7); yLine->SetLineWidth(3); yLine->Draw();
    }
  }

  legend->Draw();
  c->Update();
  legend->Draw();
  c->SaveAs( outputDir + effHistName + extra + "_Contour.pdf" );
  c->SaveAs( outputDir + effHistName + extra + "_Contour.png" );
  c->Write();

}


