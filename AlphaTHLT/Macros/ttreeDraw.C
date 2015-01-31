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
#include <iostream>

#include "typeDefs.h"



void cumulative2D( TH2* histogram, TH2* rCumulHist, double scale){

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
      double integral = histogram->Integral( 0, iBinX, 0, iBinY );
      rCumulHist->SetBinContent( iBinX, iBinY, integral );
    }
  }

  rCumulHist->Scale( scale );

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



void makeDiff2D(TChain *t, TH2D* h, TString title, TString plot, TString cut){
  
  h->SetName( title );
  t->Project( title, plot, cut );

  h->Write();

}

void makeCumul2D(TChain *t, TH2D* h, TString title, TString plot, TString cut){

  title += "_Cumul";

  h->SetName( title );
  t->Project( title, plot, cut );

  reverseCumulative2D( h, h, 1.0);
  h->Write();

}

void make2D(TChain *t, TH2D* hist, TString title, TString plot, TString cut, bool veto = false){

  std::cout << "\tProducing: " << title << "\n";

  TH2D* h = (TH2D*)hist->Clone();

  // Make differential distribution
  h->SetName( title );
  t->Project( title, plot, cut );
  h->Write();

  // Make cumulative distribution (or efficiency if relevant)
  TH2D* hCumul = (TH2D*)h->Clone();
  hCumul->Clear();
  hCumul->SetName( title + "_Cumul");
  if (!(title.Contains("QCD"))){
    if (veto){ cumulative2D( h, hCumul, 1/h->GetEntries());        } // Veto selected region
    else     { reverseCumulative2D( h, hCumul, 1/h->GetEntries()); }
    hCumul->SetMinimum( 0. );
    hCumul->SetMaximum( 1. );
  }
  else{
    if (veto){ cumulative2D( h, hCumul, 1.0);        }// Veto selected region
    else     { reverseCumulative2D( h, hCumul, 1.0); }
  }
  hCumul->Write();

}






// Iterate over input cutCollection and produce a plot for each
void make2DDistributions(TChain *t, TH2D* hist, TString title, TString plot, cutCollection cuts, bool veto = false){

  for (uint iCut = 0; iCut < cuts.size(); ++iCut){
    TString suffix = cuts[ iCut ].first;
    TString cut    = cuts[ iCut ].second;

    // Remove muon veto for DY to improve stats
    if (title.Contains("DY_") ){ cut.ReplaceAll("genLeptonVeto", "genElectronVeto"); }
    
    make2D( t, hist, title + "_" + suffix, plot, cut, veto);
  }
  
}


void make2DDistributions(sampleCollection samples, TH2D* hist, TString plot, cutCollection cuts, bool veto = false,  TString suffix = ""){

  for (uint iSample = 0; iSample < samples.size(); ++iSample){

    sample currentSample = samples[ iSample ];
    if (!currentSample.process){continue;}

    TString title  = currentSample.name + "_" + hist->GetName();
    if (suffix != ""){ title += "_" + suffix; }

    for (uint iCut = 0; iCut < cuts.size(); ++iCut){
      TString cutLabel = cuts[ iCut ].first;
      TString cut      = cuts[ iCut ].second;

      // Remove muon veto for DY to improve stats
      if ( currentSample.name == "DY" ){ cut.ReplaceAll("genLeptonVeto", "genElectronVeto"); }
      
      make2D( currentSample.chain, hist, title + "_" + cutLabel, plot, cut, veto);
    }
  }
  
}
