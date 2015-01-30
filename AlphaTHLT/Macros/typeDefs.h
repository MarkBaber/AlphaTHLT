#ifndef TYPEDEFS
#define TYPEDEFS



struct sample{
  TChain *chain;
  TString name;
  bool    process;
  sample(TString aName, TString aFiles, bool aProcess){
    chain = new TChain("MakeTrees/Ntuple"); chain->Add( aFiles );
    name = aName;
    process = aProcess;
  }
};

typedef std::vector< sample >                      sampleCollection;
typedef std::vector< std::pair<TString, TString> >    cutCollection;


#endif
