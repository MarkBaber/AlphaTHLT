#ifndef ALPHAT
#define ALPHAT


#include <vector>
#include "TLorentzVector.h"


// If true, testJet is in tagJet hemisphere 
inline bool inTagJetHemisphere( const TLorentzVector &testJet, const TLorentzVector &tagJet ){
  float dotProd = tagJet.Vect().Dot( testJet.Vect() );
  return ( dotProd >= 0 );
}


float calculateAlphaTPrime( float MHT, float HT ){
  float MHTOverHT(0);
  float maxAlphaT  = 1.0; // restrict alphaT to smaller range 
  if (HT > 0){
    MHTOverHT = MHT/HT;
  }
  float alphaTPrime = 0.5*( 1/sqrt( 1 - MHTOverHT*MHTOverHT )  );
  if (alphaTPrime > maxAlphaT){ alphaTPrime = maxAlphaT; }

  return alphaTPrime;
}



float calculateHTBetaT( std::vector< TLorentzVector > jetVec, int nMinTag, bool scalar ){

  //  uint maxJets = 15;

  // check the size of the input collection
  if (jetVec.size() <= 1){
    // empty jet collection, return BetaT = 0
    return 0.;
  }

  // Calculate deltaHT in hemispheres
  double delta_sum_et = 0.;
  TLorentzVector delta_sum_et_Vec;

  if ( nMinTag != -1){
    
    TLorentzVector tagJet = jetVec[ nMinTag ];
    
    // Loop through jets
    for (uint iJet = 0;iJet < jetVec.size(); ++iJet){
      // if (iJet == nMinTag){continue;}
      TLorentzVector testJet = jetVec[ iJet ];

      if (scalar){ // DeltaHT calculated from scalar sum of jet PT

	if ( inTagJetHemisphere( testJet, tagJet ) ){
	  delta_sum_et += testJet.Pt();
	} // In tag hemisphere
	else{
	  delta_sum_et -= testJet.Pt();
	} // In opposite hemisphere
      }
      else{ // DeltaHT calculated from vector sum of jet PT
	
	if ( inTagJetHemisphere( testJet, tagJet ) ){
          delta_sum_et_Vec += testJet;
        } // In tag hemisphere                                                                                                              
        else{
          delta_sum_et_Vec -= testJet;
        } // In opposite hemisphere  
      }
                 
    }
  }


  if (!scalar){
    delta_sum_et = delta_sum_et_Vec.Pt();
  }
  delta_sum_et = std::fabs(delta_sum_et);
  

  if ( delta_sum_et == 99999.){
    delta_sum_et = 0 ;
  }

  // Beta_T
  return delta_sum_et; 
}




float calculateHTAlphaT( std::vector< TLorentzVector > jetVec, bool scalar ){

  //  uint maxJets = 15;

  // check the size of the input collection
  if (jetVec.size() <= 1){
    // empty jet collection, return BetaT = 0
    return 0.;
  }

  double min_delta_sum_et = 99999.;
  if (scalar){ 
    for (unsigned int i = 0; i < (1U << (jetVec.size() - 1)); i++) { //@@ iterate through different combinations
      double delta_sum_et = 0.;
    
      for (unsigned int j = 0; j < jetVec.size(); ++j) { //@@ iterate through jets
	if (i & (1U << j))
	  delta_sum_et -= jetVec[j].Pt();
	else
	  delta_sum_et += jetVec[j].Pt();
      }
    
      delta_sum_et = std::fabs(delta_sum_et);
      if (delta_sum_et < min_delta_sum_et) {
	min_delta_sum_et = delta_sum_et;
      }
    
    }
  }
  else{

    for (unsigned int i = 0; i < (1U << (jetVec.size() - 1)); i++) { //@@ iterate through different combinations
      double delta_sum_et = 0.;
      TLorentzVector delta_sum_et_Vec;

      for (unsigned int j = 0; j < jetVec.size(); ++j) { //@@ iterate through jets
	TLorentzVector testJet = jetVec[ j ];

	if (i & (1U << j))
	  delta_sum_et_Vec -= testJet;
	else
	  delta_sum_et_Vec += testJet;
      }
    
      delta_sum_et = delta_sum_et_Vec.Pt();
      if (delta_sum_et < min_delta_sum_et) {
	min_delta_sum_et = delta_sum_et;
      }
    
    }
  }
                 
  min_delta_sum_et = std::fabs(min_delta_sum_et);
  
  if ( min_delta_sum_et == 99999.){
    min_delta_sum_et = 0 ;
  }

  // DeltaHT
  return min_delta_sum_et;
}




float calculateAlphaTVec( std::vector< TLorentzVector > jetVec ){

  uint maxJets = 15;
  float maxAlphaT  = 1.0; // restrict alphaT to smaller range 
  float overAlphaT = 1.1; // restrict 'overflow' alphaT to smaller range 

  // check the size of the input collection
  if (jetVec.size() <= 1){
    // empty jet collection, return AlphaT = 0
    return 0.;
  }
  // Keep event if it contains many jets
  if (jetVec.size() >= maxJets){
    return overAlphaT;
  }

  
  // Momentum sums in transverse plane
  float sum_et(0), sum_px(0), sum_py(0);

  for (unsigned int iJet = 0; iJet < jetVec.size(); ++iJet ){

    sum_et += jetVec[iJet].Pt(); 
    sum_px += jetVec[iJet].Px(); 
    sum_py += jetVec[iJet].Py(); 
    
  }

  double min_delta_sum_et = calculateHTAlphaT( jetVec, false );
  
  // Return a large value of alphaT
  if ( (sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py)) <= 0 ){ 
    return maxAlphaT;
  }

  // Alpha_T
  float alphaT = (0.5 * (sum_et - min_delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) ));  
  if (alphaT > maxAlphaT){alphaT = maxAlphaT;}
  return alphaT;

}





// ----------------------------------------
// CalculateBetaT
// ----------------------------------------
// Calculate BetaT for input jet collection jetVec and index of selected tagJet nMinTag. Scalar switches between 
// HT determination from scalar and vector addition of jet PT.
//
// Forms pseudojets from jets in hemispheres. The hemisphere boundary is defined as the surface normal to the tag jet vector.
float calculateBetaT( std::vector< TLorentzVector > jetVec, uint nMinTag, bool scalar ){

  float maxBetaT  = 1.0; // restrict alphaT to smaller range  
  float overBetaT = 1.1; // restrict 'overflow' alphaT to smaller range 

  // check the size of the input collection
  if (jetVec.size() <= 1){
    // empty jet collection, return BetaT = 0
    return 0.;
  }
  
  // Momentum sums in transverse plane
  float sum_et(0), sum_px(0), sum_py(0);
  for (unsigned int iJet = 0; iJet < jetVec.size(); ++iJet ){
    sum_et += jetVec[iJet].Pt();    
    sum_px += jetVec[iJet].Px();
    sum_py += jetVec[iJet].Py();
  }

  // Calculate deltaHT in hemispheres
  double delta_sum_et = 0.;
  TLorentzVector delta_sum_et_Vec;
  
  
  TLorentzVector tagJet = jetVec[ nMinTag ];
  
  // Loop through other jets
  for (uint iJet = 0;iJet < jetVec.size(); ++iJet){
    //if (iJet == nMinTag){continue;}
    TLorentzVector testJet = jetVec[ iJet ];
    
    if (scalar){ // DeltaHT calculated from scalar sum of jet PT
      
      if ( inTagJetHemisphere( testJet, tagJet ) ){
	delta_sum_et += testJet.Pt();
      } // In tag hemisphere
      else{
	delta_sum_et -= testJet.Pt();
      } // In opposite hemisphere
    }
    else{ // DeltaHT calculated from vector sum of jet PT
      
      if ( inTagJetHemisphere( testJet, tagJet ) ){
	delta_sum_et_Vec += testJet;
      } // In tag hemisphere                                                                                                              
      else{
	delta_sum_et_Vec -= testJet;
      } // In opposite hemisphere  
    }
    //      std::cout << inTagHemi << "\n";
  }


  if (!scalar){
    delta_sum_et = delta_sum_et_Vec.Pt();
  }
  delta_sum_et = std::fabs(delta_sum_et);
  
  // Return a large value of BetaT
  if ( (sum_et*sum_et - (sum_px*sum_px + sum_py*sum_py)) <= 0 ){ 
    return overBetaT;
  }

  // Beta_T
  float betaT = (0.5 * (sum_et - delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) )); 
  if (betaT > maxBetaT){betaT = maxBetaT;}
  return betaT;
}





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





// Returns all HT alphaT pairs
//std::vector<std::pair<float,float> > 
void calculateDynamicAlphaTPairs( const std::vector<const reco::Candidate*>& jetColl, float dynamicJetThreshold,
				  std::vector<float>& alphaTVec, std::vector<float>& HTVec){

  float currentHT(0);
  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;
  //std::vector<std::pair<float,float> > alphaTHTPairs;

  alphaTVec.clear();
  HTVec.clear();

  bool firstJet(true);
  for ( std::vector<const reco::Candidate*>::const_iterator itr = jetColl.begin(); itr != jetColl.end(); ++itr ){

    float jetPt = (*itr)->pt();
    float jetPx = (*itr)->px();
    float jetPy = (*itr)->py();
    if( jetPt < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( jetPt );
    jetPyNew .push_back( jetPx );
    jetPxNew .push_back( jetPy );
    currentHT += jetPt;
    
    if (firstJet){firstJet = false; continue;} // Require two jets
    float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 );
  
    alphaTVec.push_back( aT );
    HTVec.push_back( currentHT );
    //alphaTHTPairs.push_back( std::make_pair( aT, currentHT ) );

  }

  //return alphaTHTPairs;
}




std::pair<float,float> calculateDynamicAlphaTHT(std::vector<float> *jetPT, std::vector<float> *jetPx, std::vector<float> *jetPy,//float jetThreshold,
 						uint maxJets, 
						float dynamicJetThreshold, float dynamicAlphaTThreshold, float dynamicHTThreshold){
  
  float maxAlphaT(-1);
  float currentHT(0);
  

  std::vector<float> jetPTNew;
  std::vector<float> jetPxNew;
  std::vector<float> jetPyNew;

  uint jetLimit = std::min( (uint)jetPT->size(), maxJets );

  for(uint iJet = 0;iJet < jetLimit; ++iJet ){
    if( (*jetPT)[iJet] < dynamicJetThreshold ) break;

    // Add more jets to the alphaT collection
    jetPTNew .push_back( (*jetPT)[iJet] );
    jetPyNew .push_back( (*jetPx)[iJet] );
    jetPxNew .push_back( (*jetPy)[iJet] );
    currentHT += (*jetPT)[iJet];

      float aT = calculateAlphaT( &jetPTNew, &jetPxNew, &jetPyNew, 0 );
      if ( aT > maxAlphaT ){ maxAlphaT = aT; }
      if ( (maxAlphaT > dynamicAlphaTThreshold) && (currentHT > dynamicHTThreshold) ){
	return std::make_pair( maxAlphaT, currentHT );
      }
  }

  return std::make_pair( -1, 0 );
}







// ****************************************************************************************************
// *                                              AlphaT                                              *
// ****************************************************************************************************

std::pair<float, float> calculateAlphaTHT(const std::vector<const reco::Candidate*>& jets, float jetThreshold){
     
  float maxAlphaT  = 1.0; // restrict alphaT to smaller range 
  //  float overAlphaT = 1.1; // restrict 'overflow' alphaT to smaller range 

  // Momentum sums in transverse plane
  float sum_et(0), sum_px(0), sum_py(0);

  // check the size of the input collection
  if (jets.size() <= 1){
    return std::make_pair(0., sum_et);
  }

  // Jet collection restricted to jet threshold
  std::vector<float> jetPTNew;

  for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){

    if ( jets.at(iJet)->pt() < jetThreshold ){ break; }
    jetPTNew.push_back( jets.at(iJet)->pt() );

    sum_et += jets.at(iJet)->pt();
    sum_px += jets.at(iJet)->px();
    sum_py += jets.at(iJet)->py();

  }
  // check the size of the new input collection 
  if (jetPTNew.size() <= 1){
    // empty jet collection, return AlphaT = 0 
    return std::make_pair( 0., sum_et);
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
    delta_sum_et = std::abs(delta_sum_et);
    if (delta_sum_et < min_delta_sum_et) {
      min_delta_sum_et = delta_sum_et;
    }

  }

  // Return a large value of alphaT 
  if ( (sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py)) <= 0 ){
    return std::make_pair(11.,  sum_et);
  }

  // Alpha_T 
  float alphaT = 0.5 * (sum_et - min_delta_sum_et) / sqrt( sum_et*sum_et - (sum_px*sum_px+sum_py*sum_py) );
  if (alphaT > maxAlphaT){ alphaT = maxAlphaT; }
  return std::make_pair( alphaT, sum_et );

}







//float calculateBiasedDeltaPhi( std::vector<float> *JetPT, std::vector<float> *JetEta, std::vector<float> *JetPhi, double jetThreshold){
float calculateBiasedDeltaPhi( const std::vector<const reco::Candidate*>& jetColl, int& tagIndex, double jetThreshold){

    std::vector< TLorentzVector > jetVec;
    TLorentzVector jetMHTVec;

    // Require at least two jets
    if ( jetColl.size() <= 1){
      return -1;
    } 

    for ( std::vector<const reco::Candidate*>::const_iterator itr = jetColl.begin(); itr != jetColl.end(); ++itr ){

      float jetPt  = (*itr)->pt();
      float jetEta = (*itr)->eta();
      float jetPhi = (*itr)->phi();
   
      TLorentzVector jet;
      jet.SetPtEtaPhiM( jetPt, jetEta, jetPhi, 0 );
      
      // Analysis-jets                                                                                                                    
      // ----------------------------------------                                                                                         
      if (jetPt < jetThreshold){ break; }
      //HT      += jetPt;
      jetVec.push_back( jet );
      jetMHTVec += jet;
    }
    jetMHTVec = -jetMHTVec;
    //    float MHT    = jetMHTVec.Pt();

    float minDeltaPhi(9999);
    int   nMinTag(-1); // Index of tagjet that minimises deltaPhi( MHT, jet )
    for (uint iJet = 0;iJet < jetVec.size(); ++iJet){
      TLorentzVector tagJet = jetVec[ iJet ];
      
      // MHT vector excluding the tag jet: Remove tag jet from vector sum
      TLorentzVector jetMHTPrimeVec = (jetMHTVec + tagJet);
      
      // Find deltaPhi between reduced MHT and tag jet
      float deltaPhi = fabs( jetMHTPrimeVec.DeltaPhi( tagJet ) );

      if (deltaPhi < minDeltaPhi){
	minDeltaPhi = deltaPhi;
	nMinTag     = iJet;
      }
    }
    if (minDeltaPhi == 9999){minDeltaPhi = -1;}

    float biasedDPhi = minDeltaPhi; 
    tagIndex         = nMinTag;

    return biasedDPhi;
}


// ****************************************************************************************************
// *                                              AlphaT                                              *
// ****************************************************************************************************

void calculateBDPhiAlphaBetaT(const std::vector<const reco::Candidate*>& jets, float jetThreshold,
			      float& biasedDPhi, int& biasedDPhiIndex,
			      float& alphaTVector,
			      float& betaTScalar, float& betaTVector){
  
  biasedDPhiIndex = -1;
     
  // Jet collection restricted to jet threshold
  std::vector<TLorentzVector> jetVec;
  for (unsigned int iJet = 0; iJet < jets.size(); ++iJet ){
    float jetPt  = jets.at(iJet)->pt();
    float jetEta = jets.at(iJet)->eta();
    float jetPhi = jets.at(iJet)->phi();
    if ( jetPt < jetThreshold ){ break; }
    TLorentzVector jet; jet.SetPtEtaPhiM( jetPt, jetEta, jetPhi, 0 ); jetVec.push_back( jet );
  }

  biasedDPhi  = calculateBiasedDeltaPhi( jets, biasedDPhiIndex, jetThreshold);
  if (biasedDPhi != -1){
    betaTScalar  = calculateBetaT( jetVec, biasedDPhiIndex, true );
    betaTVector  = calculateBetaT( jetVec, biasedDPhiIndex, false );
    alphaTVector = calculateAlphaTVec( jetVec );
  }
  else{
    betaTScalar  = 0;
    betaTVector  = 0;
    alphaTVector = 0;
  }

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
