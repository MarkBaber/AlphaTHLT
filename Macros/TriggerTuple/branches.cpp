#ifndef BRANCHES
#define BRANCHES





template<typename T>
void getBranches(TChain *aChain, std::vector<TString>branchNames, std::map<TString, T> &branchContainer){
  for (uint iName = 0; iName < branchNames.size(); ++iName){
    TString branchName = branchNames[iName];
    aChain->SetBranchAddress( branchName.Data(), &branchContainer[branchName]); aChain->SetBranchStatus( branchName.Data(), 1);
  }
}




void loadBranches( TChain *chain, 
		   std::map<TString, float>                 &fBranch, 
		   std::map<TString, UInt_t>                &uBranch,
		   std::map<TString, Int_t>                 &iBranch,
		   std::map<TString, std::vector<float>* >  &fvBranch, 
		   std::map<TString, std::vector<std::pair<float,float> >* > &pfvBranch,
		   TString leptonVeto,
		   std::vector<TString> hltPathNames){

    // Reconstruction types
    // ----------------------------------------
    std::vector<TString> recoTypes;
    recoTypes.push_back( "hltAk4PF" );
    recoTypes.push_back( "hltAk4Calo" );
    recoTypes.push_back( "genAk4" );
    //recoTypes.push_back( "recoAk4PF" );


    // Setup branches
    // --------------------------------------------------------------------------------
    
    std::vector<TString>   fBranches;
    std::vector<TString>   uBranches;
    std::vector<TString>   iBranches;
    std::vector<TString>  fvBranches;
    std::vector<TString> pfvBranches;

    // Trigger bits
    for ( auto itr : hltPathNames ){ uBranches.push_back( itr ); }

    // L1
    fBranches.push_back("gct_Ht");
    fBranches.push_back("gct_MetPt");
    uBranches.push_back("L1HTT175OrETM70");

    // MET
    fBranches.push_back("hltMetCalo_MetPT");
    fBranches.push_back("genMetCalo_MetPt");
    
    // Vetoes
    uBranches.push_back("genLeptonVeto");  
    uBranches.push_back("genElectronVeto");
    //uBranches.push_back("genPhotonVeto");
    
    // General
    for (uint iType = 0; iType < recoTypes.size(); ++iType){
      TString recoType = recoTypes[iType];
      
      // float
      fBranches.push_back( recoType + "Lead_Pt" );
      fBranches.push_back( recoType + "Second_Pt" );
      fBranches.push_back( recoType + "DijetAvg_Pt" );
      fBranches.push_back( recoType + "_AlphaT40" );
      fBranches.push_back( recoType + "_AlphaTPrime40" );
      fBranches.push_back( recoType + "_HT40" );
      fBranches.push_back( recoType + "_MhtPT40" );
      //fBranches.push_back( recoType + "_BiasedDPhi" );
      if ( recoType != "hltAk4Calo" ){
	fBranches.push_back( recoType + "For_MaxPt" );
	fBranches.push_back( recoType + "For_MhtPT40" );
	pfvBranches.push_back( recoType + "_DynamicAlphaTHT40" );      
      }

      
      // UInt_t
      //uBranches.push_back( recoType + "_NJet40" );
      //uBranches.push_back( recoType + "_BiasedDPhiIndex" );
      
      // Int_t
      iBranches.push_back( recoType + "_NJetBin40" );
      iBranches.push_back( recoType + "_HTBin40" );

      // vfloat
      fvBranches.push_back( recoType + "_Pt" );
      fvBranches.push_back( recoType + "_Px" );
      fvBranches.push_back( recoType + "_Py" );
      // fvBranches.push_back( recoType + "_Eta" );
      // fvBranches.push_back( recoType + "_Phi" );
      // fvBranches.push_back( recoType + "For_Pt" );



    }
    
    // Load branches
    chain->SetBranchStatus("*",0);
    getBranches(chain,   fBranches,   fBranch);   
    getBranches(chain,   iBranches,   iBranch);
    getBranches(chain,   uBranches,   uBranch);
    getBranches(chain,  fvBranches,  fvBranch);
    //getBranches(chain, pfvBranches, pfvBranch);
 


}



#endif
