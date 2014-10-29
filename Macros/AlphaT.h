#ifndef ALPHAT
#define ALPHAT

#include <vector>



/* float calculateAlphaT( std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold ); */
/* float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx,  */
/* 			     std::vector<float> *jetPy, //float jetThreshold,  */
/* 			     uint maxJets,  */
/* 			     float dynamicJetThreshold, float dynamicAlphaTThreshold); */
/* std::pair<float,float> calculateDynamicAlphaTHT(std::vector<float> *jetPT, std::vector<float> *jetPx,  */
/* 						std::vector<float> *jetPy, //float jetThreshold,  */
/* 						uint maxJets, */
/*                                                 float dynamicJetThreshold, float dynamicAlphaTThreshold, float dynamicHTThreshold); */

/* std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs(std::vector<float> *jetPT, std::vector<float> *jetPx,  */
/* 								 std::vector<float> *jetPy, //float jetThreshold,  */
/* 								 uint maxJets, float dynamicJetThreshold); */







float calculateAlphaT( std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy, float jetThreshold ){

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





// Returns all HT alphaT pairs
std::vector<std::pair<float,float> > calculateDynamicAlphaTPairs(std::vector<float> *jetPT, std::vector<float> *jetPx, 
								 std::vector<float> *jetPy,// float jetThreshold, 
								 uint maxJets, float dynamicJetThreshold){

  float currentHT(0);
  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;
  std::vector<std::pair<float,float> > alphaTHTPairs;

  uint jetLimit = std::min( (uint)jetPT->size(), maxJets );

  for(uint iJet = 0;iJet < jetLimit; ++iJet ){

    if( (*jetPT)[iJet] < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( (*jetPT)[iJet] );
    jetPyNew .push_back( (*jetPx)[iJet] );
    jetPxNew .push_back( (*jetPy)[iJet] );

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



float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy,
			     uint maxJets, float dynamicJetThreshold, float dynamicAlphaTThreshold){

   std::vector<std::pair<float,float> > alphaTHTPairs = calculateDynamicAlphaTPairs( jetPT, jetPx, jetPy, maxJets, dynamicJetThreshold);

   float maxAlphaT(-1); 
   for (uint iPair = 0; iPair < alphaTHTPairs.size(); ++iPair){
     float curAlphaT = alphaTHTPairs[ iPair ].first;
     if (curAlphaT > maxAlphaT){
       maxAlphaT = curAlphaT;
     }
   }
   return maxAlphaT; 
} 



// Ported from https://github.com/cms-sw/cmssw/blob/CMSSW_7_2_X/HLTrigger/JetMET/src/HLTAlphaTFilter.cc
/* float calculateDynamicAlphaT(std::vector<float> *jetPT, std::vector<float> *jetPx,  */
/* 			     std::vector<float> *jetPy, //float jetThreshold,  */
/* 			     uint maxJets,  */
/* 			     float dynamicJetThreshold, float dynamicAlphaTThreshold){ */

/*   float maxAlphaT(-1); */

/*   std::vector<float> jetPTNew; */
/*   std::vector<float> jetPxNew; */
/*   std::vector<float> jetPyNew; */
/*   //  std::vector<float> jetPhiNew; */

/*   uint jetLimit = std::min( (uint)jetPT->size(), maxJets ); */

/*   for(uint iJet = 0;iJet < jetLimit; ++iJet ){ */

/*     //    if( std::abs(ijet->eta()) > etaJet_.at(0) ) continue; */
/*     if( (*jetPT)[iJet] < dynamicJetThreshold ) break; */

/*     // Add more jets to the alphaT collection */
/*     jetPTNew .push_back( (*jetPT)[iJet] ); */
/*     jetPxNew .push_back( (*jetPx)[iJet] ); */
/*     jetPyNew .push_back( (*jetPy)[iJet] ); */
    
/*     //    if ( iJet > 0 ){ */
/*       float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 ); */
/*       if ( aT > maxAlphaT ){ maxAlphaT = aT;} */
/*       if ( maxAlphaT > dynamicAlphaTThreshold ){ break; } */
/*     } */
/*   //  } */

/*   return maxAlphaT; */
/* } */




#endif
