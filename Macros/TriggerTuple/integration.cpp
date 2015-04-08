#ifndef INTEGRATION
#define INTEGRATION


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
  for ( int iBinX = nBinsX; iBinX > 0; --iBinX ){
    for ( int iBinY = nBinsY; iBinY > 0; --iBinY ){
      histogram->SetBinContent( iBinX, iBinY, value );
    }
  }

}


#endif
