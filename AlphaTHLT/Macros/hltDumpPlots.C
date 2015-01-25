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

#include <vector>
#include <iostream>
#include <math.h>


bool plotExists(TFile *f, TString rootDir, TString plotName){
  return f->Get( rootDir + plotName);
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

	double xValue = histogram->GetXaxis()->GetBinCenter( iBinX );
	double yValue = histogram->GetYaxis()->GetBinCenter( iBinY );


	selectedValues.push_back( std::make_pair( xValue, yValue ) );
	break;
      }
 
    }
  }


  //  return selectedBins;
  return selectedValues;
}





  

void makeContourEffPlot( TFile *f, TFile *fRate, 
			 TString outputDir, 
			 TString rawRateDir, TString rateHistName,
			 TString rawEffDir,  TString effHistName, 
			 std::vector< std::pair<float, float> > maxRatePair ){

  double maxYValue = 1.;

  // DEBUGGING
  //  if (rateHistName != "OnRate_AlphaTStd_vs_HT_Jet2gt50_Inclusive_Cumul" ){ return;}

  TCanvas *c = new TCanvas(effHistName);
  c->SetGridx();
  c->SetGridy();
  TLegend *legend = new TLegend( 0.7,0.65, 0.85, 0.85);
  legend->SetFillColor( 0 );


  TH2D* histogramRate       =        (TH2D*)fRate->Get( rawRateDir + rateHistName );
  TEfficiency* histogramEff = (TEfficiency*)f->Get( rawEffDir  + effHistName );
  //  TH2D* histogramEff2D = (TH2D*)histogramEff->CreateHistogram();

  // Draw histogram
  histogramEff->Draw("COLZ");
    //  histogramEff2D->Draw("COLZ");
  gPad->Update();
  histogramEff->GetPaintedHistogram()->SetMaximum(1.0);
  histogramEff->GetPaintedHistogram()->Draw("COLZ");

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
	  //	  std::cout << "0" << "\t" << xx[ 0 ] << "\t" << yy[ 0 ] << "\n";
	}

	xx[ iPair+1 ] = ratePairsMin[ iPair ].first;
	yy[ iPair+1 ] = ratePairsMin[ iPair ].second;

	//	std::cout << iPair+1 << "\t" << xx[ iPair+1 ] << "\t" << yy[ iPair+1 ] << "\n";
      }

      //      std::cout << "MAXXX\n\n";
      for (int iPair = ratePairsMax.size() -1;iPair >-1; --iPair ){

	int index = nPairsTotal - 2 - iPair;

	xx[ index ] = ratePairsMax[ iPair ].first;
	yy[ index ] = ratePairsMax[ iPair ].second;

	//	std::cout << index << "\t" << xx[ index ] << "\t" << yy[ index ] << "\n";


  	if (iPair == 0 ){
	  xx[ index + 1] = ratePairsMax[ iPair ].first;
	  yy[ index + 1] = maxYValue;
	  //	  std::cout << index + 2 << "\t" << xx[ index + 2] << "\t" << yy[ index + 2] << "\n";
	}

      }

      //      std::cout  << "\n";
      

      for  (int iPair = 0;iPair <nPairsTotal; ++iPair ){

	//	std::cout << iPair << "\t" << xx[ iPair ] << "\t" << yy[ iPair ] << "\n";
      }

      TString maxRateStr = Form("%d", int(maxRateLow)) + TString(" < Rate < ") + Form("%d", int(maxRateHigh)) + TString(" Hz");
      //      TSpline3 *spline = new TSpline3( maxRateStr,xx,yy,nPairs+nPairsMax);


      
      int colour = iRate + 1;
      if (colour == 3){ colour = kYellow; }

      // spline  ->SetLineColor( colour );
      // spline  ->SetLineWidth(2);
      // spline  ->SetLineStyle(7);
      // spline  ->SetMarkerSize( 0.8 );
      // spline  ->SetMarkerStyle( 21 );
      // spline  ->SetMarkerColor( colour );
      
      // spline  ->Draw("CPSAME");


      TGraph* graph       = new TGraph( nPairsTotal, xx, yy);
      TGraph* graphHollow = new TGraph( nPairsTotal, xx, yy);
      graph->SetFillColor( colour );
      graph->SetFillStyle( 3001 );
      graphHollow->SetFillColor( colour );
      graphHollow->SetFillStyle( 0 );

      graph      ->Draw("FL");
      graphHollow->Draw("FL");

      // gr->SetFillStyle( 3001 );
      // gr->Draw("FL");
      legend ->AddEntry( graph, maxRateStr);
      delete[] xx;
      delete[] yy;

    }
  }

  //  legend->SetHeader("The Legend Title");
  legend->Draw();
  c->Update();
  legend->Draw();
  c->SaveAs( outputDir + effHistName + "_Contour.pdf" );
  c->Write();

}








void hltDumpPlots(){

  gStyle->SetOptStat(0);


  // ********************************************************************************
  // Configuration
  // ********************************************************************************

  TString sample     = "TTBar"; //"TTBar";
  TString rateSample = "QCD";
  //TString sample = "DY";

	//  TString directory = "04_08_14/";
  TString directory = "01_10_14/";

  TString filename     = directory + sample + ".root";
  TString filenameRate = directory + rateSample + ".root";
  TString outputDir = directory + sample + "/"; 
  TString rateOutputDir = "Rate/"; 
  
  TString effDir  = "Triggers/Efficiency/";
  TString rateDir = "Triggers/Rate/";


  std::vector<TString> jetBins;
  // jetBins.push_back("Inclusive");
  // jetBins.push_back("2Jets");
  // jetBins.push_back("3Jets");
  // jetBins.push_back("4Jets");
  // jetBins.push_back("5Jets");
  // jetBins.push_back("6Jets");

  std::vector<TString> alphaTTypes;
  alphaTTypes.push_back("Std");
  alphaTTypes.push_back("Dyn");
  alphaTTypes.push_back("Dy2");
  alphaTTypes.push_back("Dy3");

  std::vector<TString> alphaTThresholds;
  alphaTThresholds.push_back("gt0p60");
  alphaTThresholds.push_back("gt0p65");
  alphaTThresholds.push_back("gt0p70");
  alphaTThresholds.push_back("gt0p75");
  alphaTThresholds.push_back("gt0p80");

  std::vector<TString> jet2PTThresholds;
  jet2PTThresholds.push_back("gt50");
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
  // anaHtBins.push_back("HT200to275");
  // anaHtBins.push_back("HT275to325");
  // anaHtBins.push_back("HT325to375");
  // anaHtBins.push_back("HT375to9999");
  anaHtBins.push_back("HT200to300");
  anaHtBins.push_back("HT300to400");
  anaHtBins.push_back("HT400to500");
  anaHtBins.push_back("HT500to9999");


  // ********************************************************************************
 

  TFile *f = new TFile( filename,"OPEN");
  TFile *fRate = new TFile( filenameRate,"OPEN");
  if ( !(f->IsOpen()) ){
    std::cout << "Error: The file '" << filename << "' could not be opened.\n";
    return;
  }
  if ( !(fRate->IsOpen()) ){
    std::cout << "Error: The file '" << filenameRate << "' could not be opened.\n";
    return;
  }

  std::cout << filename << "\n";

  TString rawRateDir   = "Raw/Rate/";
  TString rawEffDir    = "Raw/Efficiency/";
  
  std::vector< std::pair<float, float> > maxRatePairs;
  // maxRatePairs.push_back( std::make_pair( 1, 10) );
  // maxRatePairs.push_back( std::make_pair( 10, 20) );
  // maxRatePairs.push_back( std::make_pair( 20, 50) );
  // maxRatePairs.push_back( std::make_pair( 50, 100) );

  maxRatePairs.push_back( std::make_pair( 0.5, 10) );
  maxRatePairs.push_back( std::make_pair( 10,15) );
  maxRatePairs.push_back( std::make_pair( 15,20) );


  std::vector<int> maxRates;
  //    maxRates.push_back( 1 );
  //  maxRates.push_back( 5 );
  maxRates.push_back( 10 );
  maxRates.push_back( 15 );
  //  maxRates.push_back( 20 );
  // maxRates.push_back( 100 );
  // maxRates.push_back( 1000 );
  
  TFile *o = new TFile( directory + sample + "_Graphs.root", "RECREATE" );
      
  gStyle->SetTitleStyle(0);

  for (uint iType = 0;iType < alphaTTypes.size();++iType){
    
    TString aTType = alphaTTypes[iType];

    if (!aTType.Contains("Dyn")){ 
      aTType = aTType.ReplaceAll("Dy","Dyn"); 
    }

    for (uint iJet2PT = 0; iJet2PT < jet2PTThresholds.size(); ++iJet2PT){

      TString jet2PTStr = jet2PTThresholds[ iJet2PT ];

      for (uint iTrig = 1;iTrig <= 6; ++iTrig){
	TString trig = TString("Trig") + Form("%d", iTrig);

	// TString rateHistName = "OnRate_AlphaTStd_vs_HT_Jet2gt100_Inclusive_Cumul";
	// TString effHistName  = "On_AlphaTStd_vs_HT_Jet2gt100_Inclusive_UniCumul";

	TString rateHistName = trig + "OnRate_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_Inclusive_Cumul";

	for (uint iJet = 0;iJet < jetBins.size(); ++iJet){
	  TString jetBin = jetBins[iJet];
	  
	  TString effHistName  = "On_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_" + jetBin + "_UniCumul";
	
	  
	  // TH2D* histogramRate       =        (TH2D*)f->Get( rawRateDir + rateHistName );
	  // TEfficiency* histogramEff = (TEfficiency*)f->Get( rawEffDir  + effHistName );
	  
	  makeContourEffPlot( f, fRate, outputDir, rawRateDir, rateHistName, rawEffDir, effHistName, maxRatePairs );

	}

	// Analysis binned plots
	for (uint iJet = 0;iJet < anaJetBins.size(); ++iJet){
	  TString jetBin = anaJetBins[ iJet ];

	  for (uint iHT = 0;iHT < anaHtBins.size(); ++iHT){
	    TString htBin = anaHtBins[ iHT ];
	    
	    TString effAnaHistName    = trig + "On_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_" + htBin + "_" + jetBin  + "_UniCumul";
	    //TString effAnaHistName         = "On_AlphaT" + aTType + "_vs_HT_Jet2" + jet2PTStr + "_" + htBin + "_" + jetBin + "_UniCumul";

	    makeContourEffPlot( f, fRate, outputDir, rawRateDir, rateHistName, rawEffDir, effAnaHistName, maxRatePairs );

	  }
	} // End analysis jet binning

      }
    }
  }


  // // Isocontours for each rate
  // std::vector< std::pair< double, double> > ratePairs10 = getEfficiencyBins( histogramRate, 10. );
  // std::vector< std::pair< double, double> > ratePairs20 = getEfficiencyBins( histogramRate, 20. );
  // std::vector< std::pair< double, double> > ratePairs30 = getEfficiencyBins( histogramRate, 30. );


    


  exit(0);

  // ************************************************************
  // Efficiency plots
  // ************************************************************
  for (uint iType = 0;iType < alphaTTypes.size();++iType){
    
    TString aTType = alphaTTypes[iType];
    for (uint iJet = 2;iJet < 7;++iJet){
      
      TString jetMult = Form("%d", iJet);
      TString plotName    = "AlphaT" + aTType + "_vs_HT_" + jetMult + "Jets";
      TString plotNameCut = "AlphaT" + aTType + "_vs_HT_Jet2gt100_" + jetMult + "Jets";


      makePlot(           f, outputDir, effDir, plotName );
      makeCumulativePlot( f, outputDir, "Raw/Efficiency/", plotName );
      
      
      makePlot(           f, outputDir, effDir, plotNameCut );
      makeCumulativePlot( f, outputDir, "Raw/Efficiency/", plotNameCut );

    }

  }

  // ************************************************************
  // Rate plots
  // ************************************************************
  // "Jet2PT_vs_Jet1PT_alphaTDyngt0p65_4Jets_Rate"


  for (uint iType = 0;iType < alphaTTypes.size();++iType){
    TString aTType = alphaTTypes[iType];
    if ( aTType == "Std" ){ aTType = ""; }
    else if (aTType.Contains("Dy3")){ continue; }
    else if ( !(aTType.Contains("Dyn")) ){ aTType.ReplaceAll("Dy","Dyn"); }

    for (uint iThresh = 0;iThresh < alphaTThresholds.size();++iThresh){
      TString aTThresh = alphaTThresholds[ iThresh ];
      TString plotName = "Jet2PT_vs_Jet1PT_alphaT" + aTType + aTThresh + "_Rate";
    
      std::cout <<  plotName << "\n";
      makeJet12RatePlot( f, rateOutputDir, rateDir,  plotName );
    
    }
  }


}
