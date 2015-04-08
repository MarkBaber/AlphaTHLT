#ifndef ALPHAT
#define ALPHAT


#include <vector>
#include "TLorentzVector.h"




float calculateAlphaT( std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold ){

  uint maxJets     = 15;
  float maxAlphaT  = 1.0; // restrict alphaT to smaller range
  float overAlphaT = 1.1; // restrict 'overflow' alphaT to smaller range

  // check the size of the input collection
  if (jetPT->size() <= 1){
    // empty jet collection, return AlphaT = 0
    return 0.;
  }
  
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

  // Keep event if it contains many jets
  if (jetPTNew.size() >= maxJets){
    return overAlphaT;
  }

  // check the size of the new input collection
  if (jetPTNew.size() <= 1){
    // empty jet collection, return AlphaT = 0
    return 0.;
  }

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
    
    delta_sum_et = std::fabs(delta_sum_et);
    if (delta_sum_et < min_delta_sum_et) {
	min_delta_sum_et = delta_sum_et;
    }
    
  }
  
  // Return a large value of alphaT
  if ( (sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py)) <= 0 ){ 
    return overAlphaT;
  }

  // Alpha_T
  float alphaT = (0.5 * (sum_et - min_delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) )); 
  if (alphaT > maxAlphaT){alphaT = maxAlphaT;}
  return alphaT;

}







std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs( std::vector<float> *jetPT, std::vector<float> *jetPX, 
								  std::vector<float> *jetPY, float dynamicJetThreshold){

  float currentHT(0);
  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;
  std::vector<std::pair<float,float> > alphaTHTPairs;

  for (unsigned int iJet = 0; iJet < jetPT->size(); ++iJet ){

    float jetPt  = (*jetPT)[iJet];
    float jetPx  = (*jetPX)[iJet];
    float jetPy  = (*jetPY)[iJet];
    if( jetPt < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( jetPt );
    jetPyNew .push_back( jetPx );
    jetPxNew .push_back( jetPy );
    currentHT += jetPt;

    float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 );
    alphaTHTPairs.push_back( std::make_pair( aT, currentHT ) );

  }

  return alphaTHTPairs;
}

#endif
