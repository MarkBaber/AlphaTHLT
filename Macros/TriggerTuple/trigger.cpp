#ifndef TRIGGER
#define TRIGGER

#include <sstream>

// Trigger object - passed boolean value points to 
struct trigger{
  TString name;
  float curWeight;

  // Number passed entire sample
  bool *passed;
  int nPassed;
  int nFailed;
  float rate;
  // Number passed given signal selection
  bool *sigPassed;
  int nPassedGivenSig;
  int nFailedGivenSig;

 
  trigger():name(""),curWeight(1.),passed(NULL),nPassed(0),nFailed(0),rate(0),sigPassed(NULL),nPassedGivenSig(0),nFailedGivenSig(0){};
  trigger(TString aName, bool* aPassed, bool* aSigPassed=NULL):name(aName),curWeight(1.),passed(aPassed),nPassed(0),nFailed(0),rate(0),sigPassed(aSigPassed),nPassedGivenSig(0),nFailedGivenSig(0){} //if (aSigPassed == NULL){ *aSigPassed = 0; }};

  void setWeight(float weight)    { curWeight = weight; }
  void setSignal(bool &aSigPassed){ sigPassed = &aSigPassed; }

  float sampleEfficiency(){
    if (nPassed + nFailed == 0) return 0.;
    return float(nPassed)/(nPassed + nFailed);
  }
  float signalEfficiency(){
    if (nPassedGivenSig + nFailedGivenSig == 0) return 0.;
    return float(nPassedGivenSig)/(nPassedGivenSig + nFailedGivenSig);
  }



  float getRateError(){
    float error = 0;
    if ( nPassed > 0 ){ error = rate/sqrt( nPassed ); }
    return error;
  }

  TString setPrecisionStr(float n, int precision){
    std::ostringstream out;
    out << std::setprecision(precision) << std::fixed << n;
    return TString(out.str());
  }


  bool checkTrigger(){
    // Events in sample passing selection
    if (*passed == 1)   { nPassed++; rate += curWeight; }
    else                { nFailed++;}
    // Events in signal region passing selection
    if ((sigPassed != NULL) && (*sigPassed == 1) ){ 

      if (*passed == 1)   { nPassedGivenSig++;}
      else                { nFailedGivenSig++;}
    }

    return *passed;
  };
};


typedef std::vector< trigger > triggerCollection;
struct triggerPath{

  triggerCollection triggers;
  

  void setWeight(float weight){
    for (uint iTrig = 0; iTrig < triggers.size(); ++iTrig){ triggers[iTrig].setWeight(weight); }
  }

  void setSignal(bool &aSignal){
    for (uint iTrig = 0; iTrig < triggers.size(); ++iTrig){ triggers[iTrig].setSignal( aSignal ); };
  }


  void addTrigger( trigger aTrigger ){ triggers.push_back( aTrigger ); };

  

  void checkEvent(){ // Check triggers, halt when a trigger requirement fails
    for (uint iTrig = 0; iTrig < triggers.size(); ++iTrig){
      if ( triggers[iTrig].checkTrigger() == false ) {	break; }
    }
  };

  void print(bool header = true){

    int titleWidth  = 39;
    int columnWidth = 13;
    int rateWidth   = 26;
    float cumulSampleEff = 1.;
    float cumulSignalEff = 1.;

    ios init(NULL); init.copyfmt(std::cout);

    if (header){
    std::cout << "\n";
    std::cout << std::setw(titleWidth)  << "Trigger";
    std::cout << std::setw(columnWidth) << "Passed" 
	      << std::setw(columnWidth) << "Failed" 
              << std::setw(rateWidth)   << "Rate (Hz)"
	      << std::setw(columnWidth) << "EffSample" 
	      << std::setw(columnWidth) << "CumulEff" 
	      << std::setw(columnWidth) << "EffSignal" 
	      << std::setw(columnWidth) << "CumulEff" << "\n\n";
    }
    for (uint iTrig = 0; iTrig < triggers.size(); ++iTrig){
      trigger currTrig = triggers[iTrig];
      cumulSampleEff *= currTrig.sampleEfficiency();
      cumulSignalEff *= currTrig.signalEfficiency();

      std::cout << std::setw(titleWidth) 
		<< currTrig.name          
		<< std::setw(columnWidth) << currTrig.nPassed       
		<< std::setw(columnWidth) << currTrig.nFailed  
		<< std::setw(rateWidth)   
		<< currTrig.setPrecisionStr(currTrig.rate, 2) + " +/- " + currTrig.setPrecisionStr(currTrig.getRateError(), 2)
 		<< std::setw(columnWidth) << currTrig.setPrecisionStr(currTrig.sampleEfficiency(), 2)
		<< std::setw(columnWidth) << currTrig.setPrecisionStr(cumulSampleEff, 2)
 		<< std::setw(columnWidth) << currTrig.setPrecisionStr(currTrig.signalEfficiency(), 2)
		<< std::setw(columnWidth) << currTrig.setPrecisionStr(cumulSignalEff, 2)
		<< "\n";

    }

    std::cout.copyfmt(init);

    // for( triggerCollection::iterator itr = triggers.begin(); itr != triggers.end(); ++itr){ 
    //   std::cout << itr.name    << "\t";
    //   std::cout << itr.nPassed << "\t" << itr.nFailed "\n";
    // }
  };

};





// ****************************************************************************************************
// *                                           Trigger menu                                           *
// ****************************************************************************************************

struct triggerMenu{

  TString name;
  bool rate;
  std::map< TString, triggerPath > path;
  triggerMenu(){};
  triggerMenu(TString aName):name(aName){};
  
  void printPaths(bool title=false){
    for( auto itr : path){ 
      itr.second.print(title); 
      std::cout << "\n";
    };
  }

  void newPath(TString aPathname){
    path[aPathname] = triggerPath();
  }

  void addPath(TString aPathname, triggerPath aPath){
    path[aPathname] = aPath;
  }
  void addTrigger(TString aPathname, trigger aTrigger){
    path[aPathname].addTrigger( aTrigger );
  }


  void setWeight(float weight){ for ( auto &itr : path ) { itr.second.setWeight( weight ); } }
  void setSignal(TString aPathname, bool &aSignal){
    path[aPathname].setSignal( aSignal );
  }



  void checkEvent(){
    for (auto &itr : path)  { itr.second.checkEvent(); } 
  }


  


};
#endif
