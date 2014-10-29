#ifndef JET_MATCH
#define JET_MATCH

#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <sstream>
#include "TLorentzVector.h"

// Takes two lists of Jets as TLorentzVectors and returns the best match, by minimising delta, for every jet (with a maximum deltaR)

class pair_info {

 public:
  pair_info(int id1, int id2, double dr);
  int ID1();
  int ID2();
  double DR();
  void Print();

 private:
  int mId1;
  int mId2;
  double mDr;
};


inline double getDeltaR(const float &jet1Eta, const float &jet1Phi, const float &jet2Eta, const float &jet2Phi){
  double deta = jet1Eta - jet2Eta;
  double dphi = TVector2::Phi_mpi_pi( jet1Phi - jet2Phi );
  return double(TMath::Sqrt( deta*deta + dphi*dphi ));
}



// Constructor
pair_info::pair_info(int id1, int id2, double dr) {
  mId1 = id1;
  mId2 = id2;
  mDr = dr;
}

// Getters
double pair_info::DR() {
  return mDr;
}

int pair_info::ID1() {
  return mId1;
}

int pair_info::ID2() {
  return mId2;
}

void pair_info::Print() {

  std::cout << "ID1 = " << this->ID1() << ", ID2 = " << this->ID2() << ", DeltaR = " << this->DR() << std::endl;
  return;

}

bool sortDR(pair_info vec1, pair_info vec2) {

  return (vec1.DR() < vec2.DR());
  
}

std::vector<pair_info> make_pairs(const std::vector<float> *jet1PT, const std::vector<float> *jet1Eta, const std::vector<float> *jet1Phi, 
				  const std::vector<float> *jet2PT, const std::vector<float> *jet2Eta, const std::vector<float> *jet2Phi,
				  const float &minJet1PT, const float &minJet2PT ){

  std::vector<pair_info> pinfo;

  // Note: Assumes jet collections are PT-ordered
  for (unsigned int i = 0; i < jet1Eta->size(); ++i){
    if ( (*jet1PT)[i] < minJet1PT){ break; }

    for (unsigned int j = 0; j < jet2Eta->size(); ++j){
      if ( (*jet2PT)[i] < minJet2PT){ break; }      

      pinfo.push_back(  pair_info( i, j, getDeltaR( (*jet1Eta)[i], (*jet1Phi)[i], (*jet2Eta)[j], (*jet2Phi)[j] ) ) );
    }
  }

  return pinfo;

}

bool in_array(std::vector<int> array, int value) {

  for(unsigned int i=0; i<array.size(); i++) {
    if(array[i] == value) {
      return true;
    }
  }

  return false;

}


std::vector<int> analyse_pairs(std::vector<pair_info> & pairs, int reco_size, double DRmax) {

  // Initialise matched array with no match, "-1"
  std::vector<int> reco_matched_index(reco_size, -1 );

  std::vector<int> gen_matched;
  std::vector<int> reco_matched;

  for(unsigned int i = 0; i < pairs.size(); i++) {

    if( !( in_array( gen_matched, pairs[i].ID1()) ) && !(in_array(reco_matched, pairs[i].ID2())) && (pairs[i].DR() < DRmax) ){
      gen_matched.push_back(pairs[i].ID1());
      reco_matched.push_back(pairs[i].ID2());
      reco_matched_index[ pairs[i].ID2() ] = pairs[i].ID1();
    }

  }

  return reco_matched_index;

}



#endif
