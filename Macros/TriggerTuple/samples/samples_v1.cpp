#ifndef SAMPLES_v1
#define SAMPLES_v1

#include "../../typeDefs.h"

void loadSamples( sampleCollection &samples, std::map<TString, std::vector<TString> > &sampleStrs ){
  TString dir      = "/vols/ssd00/cms/mbaber/AlphaT/Trigger/";
  bool runRate = true;




  // PU40bx50
  // --------------------
  samples.addSample( "PU40bx50_QCD30to50",    dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD30to50.root",    runRate);
  samples.addSample( "PU40bx50_QCD50to80",    dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD50to80.root",    runRate);
  samples.addSample( "PU40bx50_QCD80to120",   dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD80to120.root",   runRate);
  samples.addSample( "PU40bx50_QCD120to170",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD120to170.root",  runRate);
  samples.addSample( "PU40bx50_QCD170to300",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD170to300.root",  runRate);
  samples.addSample( "PU40bx50_QCD300to470",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD300to470.root",  runRate);
  samples.addSample( "PU40bx50_QCD470to600",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD470to600.root",  runRate);
  samples.addSample( "PU40bx50_QCD600to800",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD600to800.root",  runRate);
  samples.addSample( "PU40bx50_QCD800to1000", dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD800to1000.root", runRate);
  samples.addSample( "PU40bx50_HPUV_QCD30to50",    dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD30to50_HPUV.root",    runRate);
  samples.addSample( "PU40bx50_HPUV_QCD50to80",    dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD50to80_HPUV.root",    runRate);
  samples.addSample( "PU40bx50_HPUV_QCD80to120",   dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD80to120_HPUV.root",   runRate);
  samples.addSample( "PU40bx50_HPUV_QCD120to170",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD120to170_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HPUV_QCD170to300",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD170to300_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HPUV_QCD300to470",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD300to470_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HPUV_QCD470to600",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD470to600_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HPUV_QCD600to800",  dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD600to800_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HPUV_QCD800to1000", dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/QCD800to1000_HPUV.root", runRate);
  samples.addSample( "PU40bx50_TTbar",        dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/TTbar.root",  runRate);
  samples.addSample( "PU40bx50_DYJets",       dir + "25Mar15_MCRUN2_72_V4A_74X_PU40bx50/DYJets.root", runRate);
  // PU40bx50 - HCAL3
  // --------------------
  samples.addSample( "PU40bx50_HCAL3_QCD30to50",    dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD30to50.root",    runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD50to80",    dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD50to80.root",    runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD80to120",   dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD80to120.root",   runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD120to170",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD120to170.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD170to300",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD170to300.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD300to470",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD300to470.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD470to600",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD470to600.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD600to800",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD600to800.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_QCD800to1000", dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD800to1000.root", runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD30to50",    dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD30to50_HPUV.root",    runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD50to80",    dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD50to80_HPUV.root",    runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD80to120",   dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD80to120_HPUV.root",   runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD120to170",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD120to170_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD170to300",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD170to300_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD300to470",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD300to470_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD470to600",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD470to600_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD600to800",  dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD600to800_HPUV.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_HPUV_QCD800to1000", dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/QCD800to1000_HPUV.root", runRate);
  samples.addSample( "PU40bx50_HCAL3_TTbar",        dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/TTbar.root",  runRate);
  samples.addSample( "PU40bx50_HCAL3_DYJets",       dir + "26Mar15_MCRUN2_72_V4A_74X_PU40bx50_HCAL3/DYJets.root", runRate);

  // PU40bx25
  // --------------------
  samples.addSample( "PU40bx25_QCD30to50",    dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD30to50.root",    runRate);
  samples.addSample( "PU40bx25_QCD50to80",    dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD50to80.root",    runRate);
  samples.addSample( "PU40bx25_QCD80to120",   dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD80to120.root",   runRate);
  samples.addSample( "PU40bx25_QCD120to170",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD120to170.root",  runRate);
  samples.addSample( "PU40bx25_QCD170to300",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD170to300.root",  runRate);
  samples.addSample( "PU40bx25_QCD300to470",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD300to470.root",  runRate);
  samples.addSample( "PU40bx25_QCD470to600",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD470to600.root",  runRate);
  samples.addSample( "PU40bx25_QCD600to800",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD600to800.root",  runRate);
  samples.addSample( "PU40bx25_QCD800to1000", dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD800to1000.root", runRate);
  samples.addSample( "PU40bx25_HPUV_QCD30to50",    dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD30to50_HPUV.root",    runRate);
  samples.addSample( "PU40bx25_HPUV_QCD50to80",    dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD50to80_HPUV.root",    runRate);
  samples.addSample( "PU40bx25_HPUV_QCD80to120",   dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD80to120_HPUV.root",   runRate);
  samples.addSample( "PU40bx25_HPUV_QCD120to170",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD120to170_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HPUV_QCD170to300",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD170to300_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HPUV_QCD300to470",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD300to470_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HPUV_QCD470to600",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD470to600_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HPUV_QCD600to800",  dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD600to800_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HPUV_QCD800to1000", dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/QCD800to1000_HPUV.root", runRate);
  samples.addSample( "PU40bx25_TTbar",        dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/TTbar.root",  runRate);
  samples.addSample( "PU40bx25_DYJets",       dir + "22Mar15_MCRUN2_72_V3A_74X_PU40bx25/DYJets.root", runRate);
  // PU40bx25 - HCAL3
  // --------------------
  samples.addSample( "PU40bx25_HCAL3_QCD30to50",    dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD30to50.root",    runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD50to80",    dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD50to80.root",    runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD80to120",   dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD80to120.root",   runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD120to170",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD120to170.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD170to300",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD170to300.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD300to470",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD300to470.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD470to600",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD470to600.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD600to800",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD600to800.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_QCD800to1000", dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD800to1000.root", runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD30to50",    dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD30to50_HPUV.root",    runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD50to80",    dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD50to80_HPUV.root",    runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD80to120",   dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD80to120_HPUV.root",   runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD120to170",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD120to170_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD170to300",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD170to300_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD300to470",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD300to470_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD470to600",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD470to600_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD600to800",  dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD600to800_HPUV.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_HPUV_QCD800to1000", dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/QCD800to1000_HPUV.root", runRate);
  samples.addSample( "PU40bx25_HCAL3_TTbar",        dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/TTbar.root",  runRate);
  samples.addSample( "PU40bx25_HCAL3_DYJets",       dir + "26Mar15_MCRUN2_72_V3A_74X_PU40bx25_HCAL3/DYJets.root", runRate);


  // SM
  sampleStrs["SM"].push_back( "PU40bx25_TTbar" );
  sampleStrs["SM"].push_back( "PU40bx25_DYJets" );
  sampleStrs["SM"].push_back( "PU40bx50_TTbar" );
  sampleStrs["SM"].push_back( "PU40bx50_DYJets" );
  sampleStrs["SM"].push_back( "PU40bx25_HCAL3_TTbar" );
  sampleStrs["SM"].push_back( "PU40bx25_HCAL3_DYJets" );
  sampleStrs["SM"].push_back( "PU40bx50_HCAL3_TTbar" );
  sampleStrs["SM"].push_back( "PU40bx50_HCAL3_DYJets" );
  sampleStrs["PU40bx25_SM"].push_back( "PU40bx25_TTbar" );
  sampleStrs["PU40bx25_SM"].push_back( "PU40bx25_DYJets" );
  sampleStrs["PU40bx50_SM"].push_back( "PU40bx50_TTbar" );
  sampleStrs["PU40bx50_SM"].push_back( "PU40bx50_DYJets" );
  sampleStrs["PU40bx25_HCAL3_SM"].push_back( "PU40bx25_HCAL3_TTbar" );
  sampleStrs["PU40bx25_HCAL3_SM"].push_back( "PU40bx25_HCAL3_DYJets" );
  sampleStrs["PU40bx50_HCAL3_SM"].push_back( "PU40bx50_HCAL3_TTbar" );
  sampleStrs["PU40bx50_HCAL3_SM"].push_back( "PU40bx50_HCAL3_DYJets" );


  // QCD
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD30to50");   
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD50to80");   
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD80to120");  
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD120to170"); 
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD170to300"); 
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD300to470"); 
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD470to600"); 
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD600to800"); 
  sampleStrs["PU40bx25_HCAL3_QCD"].push_back( "PU40bx25_HCAL3_QCD800to1000");
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD30to50");   
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD50to80");   
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD80to120");  
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD120to170"); 
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD170to300"); 
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD300to470"); 
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD470to600"); 
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD600to800"); 
  sampleStrs["PU40bx25_HCAL3_HPUV_QCD"].push_back( "PU40bx25_HCAL3_HPUV_QCD800to1000");

  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD30to50");   
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD50to80");   
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD80to120");  
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD120to170"); 
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD170to300"); 
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD300to470"); 
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD470to600"); 
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD600to800"); 
  sampleStrs["PU40bx50_HCAL3_QCD"].push_back( "PU40bx50_HCAL3_QCD800to1000");
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD30to50");   
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD50to80");   
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD80to120");  
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD120to170"); 
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD170to300"); 
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD300to470"); 
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD470to600"); 
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD600to800"); 
  sampleStrs["PU40bx50_HCAL3_HPUV_QCD"].push_back( "PU40bx50_HCAL3_HPUV_QCD800to1000");


  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD30to50");   
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD50to80");   
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD80to120");  
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD120to170"); 
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD170to300"); 
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD300to470"); 
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD470to600"); 
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD600to800"); 
  sampleStrs["PU40bx25_QCD"].push_back( "PU40bx25_QCD800to1000");
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD30to50");   
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD50to80");   
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD80to120");  
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD120to170"); 
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD170to300"); 
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD300to470"); 
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD470to600"); 
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD600to800"); 
  sampleStrs["PU40bx25_HPUV_QCD"].push_back( "PU40bx25_HPUV_QCD800to1000");
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD30to50");   
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD50to80");   
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD80to120");  
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD120to170"); 
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD170to300"); 
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD300to470"); 
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD470to600"); 
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD600to800"); 
  sampleStrs["PU40bx50_QCD"].push_back( "PU40bx50_QCD800to1000");
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD30to50");   
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD50to80");   
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD80to120");  
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD120to170"); 
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD170to300"); 
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD300to470"); 
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD470to600"); 
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD600to800"); 
  sampleStrs["PU40bx50_HPUV_QCD"].push_back( "PU40bx50_HPUV_QCD800to1000");



  // ****************************************************************************************************
  // Signal samples
  // ****************************************************************************************************



  samples.addSample( "PU20bx25_T1bbbb_2J_mGl_1000_mLSP_900", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T1bbbb_2J_mGl_1000_mLSP_900/*.root", runRate);
  samples.addSample( "PU20bx25_T1bbbb_2J_mGl_1500_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T1bbbb_2J_mGl_1500_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_T1qqqq_2J_mGl_1000_mLSP_800", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T1qqqq_2J_mGl_1000_mLSP_800/*.root", runRate);
  samples.addSample( "PU20bx25_T1qqqq_2J_mGl_1400_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T1qqqq_2J_mGl_1400_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_T1tttt_2J_mGl_1200_mLSP_800", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T1tttt_2J_mGl_1200_mLSP_800/*.root", runRate);
  samples.addSample( "PU20bx25_T1tttt_2J_mGl_1500_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T1tttt_2J_mGl_1500_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_T2bb_2J_mStop_600_mLSP_580", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2bb_2J_mStop_600_mLSP_580/*.root", runRate);
  samples.addSample( "PU20bx25_T2bb_2J_mStop_900_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2bb_2J_mStop_900_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_T2qq_2J_mStop_1200_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2qq_2J_mStop_1200_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_T2qq_2J_mStop_600_mLSP_550", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2qq_2J_mStop_600_mLSP_550/*.root", runRate);
  samples.addSample( "PU20bx25_T2tt_2J_mStop_425_mLSP_325", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2tt_2J_mStop_425_mLSP_325/*.root", runRate);
  samples.addSample( "PU20bx25_T2tt_2J_mStop_500_mLSP_325", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2tt_2J_mStop_500_mLSP_325/*.root", runRate);
  samples.addSample( "PU20bx25_T2tt_2J_mStop_650_mLSP_325", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2tt_2J_mStop_650_mLSP_325/*.root", runRate);
  samples.addSample( "PU20bx25_T2tt_2J_mStop_850_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25/T2tt_2J_mStop_850_mLSP_100/*.root", runRate);
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T1bbbb_2J_mGl_1000_mLSP_900" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T1bbbb_2J_mGl_1500_mLSP_100" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T1qqqq_2J_mGl_1000_mLSP_800" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T1qqqq_2J_mGl_1400_mLSP_100" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T1tttt_2J_mGl_1200_mLSP_800" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T1tttt_2J_mGl_1500_mLSP_100" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2bb_2J_mStop_600_mLSP_580" ); 
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2bb_2J_mStop_900_mLSP_100" ); 
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2qq_2J_mStop_1200_mLSP_100" );
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2qq_2J_mStop_600_mLSP_550" ); 
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2tt_2J_mStop_425_mLSP_325" ); 
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2tt_2J_mStop_500_mLSP_325" ); 
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2tt_2J_mStop_650_mLSP_325" ); 
  sampleStrs["PU20bx25_Signal"].push_back( "PU20bx25_T2tt_2J_mStop_850_mLSP_100" );


  samples.addSample( "PU20bx25_HCAL3_T1bbbb_2J_mGl_1000_mLSP_900", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T1bbbb_2J_mGl_1000_mLSP_900/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T1bbbb_2J_mGl_1500_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T1bbbb_2J_mGl_1500_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T1qqqq_2J_mGl_1000_mLSP_800", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T1qqqq_2J_mGl_1000_mLSP_800/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T1qqqq_2J_mGl_1400_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T1qqqq_2J_mGl_1400_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T1tttt_2J_mGl_1200_mLSP_800", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T1tttt_2J_mGl_1200_mLSP_800/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T1tttt_2J_mGl_1500_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T1tttt_2J_mGl_1500_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2bb_2J_mStop_600_mLSP_580", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2bb_2J_mStop_600_mLSP_580/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2bb_2J_mStop_900_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2bb_2J_mStop_900_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2qq_2J_mStop_1200_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2qq_2J_mStop_1200_mLSP_100/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2qq_2J_mStop_600_mLSP_550", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2qq_2J_mStop_600_mLSP_550/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2tt_2J_mStop_425_mLSP_325", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2tt_2J_mStop_425_mLSP_325/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2tt_2J_mStop_500_mLSP_325", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2tt_2J_mStop_500_mLSP_325/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2tt_2J_mStop_650_mLSP_325", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2tt_2J_mStop_650_mLSP_325/*.root", runRate);
  samples.addSample( "PU20bx25_HCAL3_T2tt_2J_mStop_850_mLSP_100", dir + "29Mar15_PHYS14_25_V1_74X_PU20BX25_HCAL3/T2tt_2J_mStop_850_mLSP_100/*.root", runRate);
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T1bbbb_2J_mGl_1000_mLSP_900" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T1bbbb_2J_mGl_1500_mLSP_100" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T1qqqq_2J_mGl_1000_mLSP_800" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T1qqqq_2J_mGl_1400_mLSP_100" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T1tttt_2J_mGl_1200_mLSP_800" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T1tttt_2J_mGl_1500_mLSP_100" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2bb_2J_mStop_600_mLSP_580" ); 
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2bb_2J_mStop_900_mLSP_100" ); 
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2qq_2J_mStop_1200_mLSP_100" );
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2qq_2J_mStop_600_mLSP_550" ); 
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2tt_2J_mStop_425_mLSP_325" ); 
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2tt_2J_mStop_500_mLSP_325" ); 
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2tt_2J_mStop_650_mLSP_325" ); 
  sampleStrs["PU20bx25_HCAL3_Signal"].push_back( "PU20bx25_HCAL3_T2tt_2J_mStop_850_mLSP_100" );


  sampleStrs["test"].push_back( "PU20bx25_HCAL3_T2tt_2J_mStop_425_mLSP_325" );
  //sampleStrs["test"].push_back( "PU40bx50_HPUV_QCD800to1000" );

}



  /*
  // Signal samples
  samples["T2bb_2J_mStop_600_mLSP_580"] = sample("T2bb_2J_mStop_600_mLSP_580", dir + "27Nov14_PU20/T2bb_2J_mStop-600_mLSP-580/T2bb_2J_mStop-600_mLSP-580*.root", runRate);
  samples["T2qq_2J_mStop_600_mLSP_550"] = sample("T2qq_2J_mStop_600_mLSP_550", dir + "27Nov14_PU20/T2qq_2J_mStop-600_mLSP-550/T2qq_2J_mStop-600_mLSP-550*.root", runRate);
  samples["T2tt_2J_mStop_425_mLSP_325"] = sample("T2tt_2J_mStop_425_mLSP_325", dir + "27Nov14_PU20/T2tt_2J_mStop-425_mLSP-325/T2tt_2J_mStop-425_mLSP-325*.root", runRate);
  samples["T2tt_2J_mStop_500_mLSP_325"] = sample("T2tt_2J_mStop_500_mLSP_325", dir + "27Nov14_PU20/T2tt_2J_mStop-500_mLSP-325/T2tt_2J_mStop-500_mLSP-325*.root", runRate);
  samples["T2tt_2J_mStop_650_mLSP_325"] = sample("T2tt_2J_mStop_650_mLSP_325", dir + "27Nov14_PU20/T2tt_2J_mStop-650_mLSP-325/T2tt_2J_mStop-650_mLSP-325*.root", runRate);
  samples["T2tt_2J_mStop_850_mLSP_100"] = sample("T2tt_2J_mStop_850_mLSP_100", dir + "27Nov14_PU20/T2tt_2J_mStop-850_mLSP-100/T2tt_2J_mStop-850_mLSP-100*.root", runRate);
  */




#endif
