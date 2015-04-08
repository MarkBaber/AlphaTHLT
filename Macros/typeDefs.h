#ifndef TYPEDEFS
#define TYPEDEFS
#include <map>


struct sample{
  TChain *chain;
  TString name;
  bool    process;
  sample():name(""),process(false){};
  sample(TString aName, TString aFiles, bool aProcess){
    if ( aFiles.Contains("HPUV") ){ chain = new TChain("Ntuple"); }
    else                          { chain = new TChain("MakeTrees/Ntuple"); }
    chain->Add( aFiles );
    name = aName;
    process = aProcess;
  }
};

struct sampleCollection{


  std::map<TString, sample> samples;
  sampleCollection(){};
  void addSample(TString aName, TString aFiles, bool aProcess){
    samples[aName] = sample( aName, aFiles, aProcess);
  }

};


struct effCut{

  TString label;
  TString numCut;
  TString denCut;
  effCut(TString aLabel, TString aNumCut, TString aDenCut):label(aLabel),numCut(aNumCut),denCut(aDenCut){};
};



//typedef std::vector< sample >                      sampleCollection;
typedef std::vector< std::pair<TString, TString> >    cutCollection;
typedef std::vector< effCut >                      effCutCollection;


#endif
