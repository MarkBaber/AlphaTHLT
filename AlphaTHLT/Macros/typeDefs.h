#ifndef TYPEDEFS
#define TYPEDEFS



struct sample{
  TChain *chain;
  TString name;
  sample(TString aName, TString aFiles){
    chain = new TChain("MakeTrees/Ntuple"); chain->Add( aFiles );
    name = aName;
  }
};

typedef std::vector< sample >                      sampleCollection;
typedef std::vector< std::pair<TString, TString> >    cutCollection;


#endif
