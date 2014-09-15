// -*- C++ -*-
//
// Package:    HLTAnalyser
// Class:      HLTJetProducer
// 
/**\class HLTJetProducer HLTJetProducer.cc AlphaTHLT/HLTAnalyser/plugins/HLTJetProducer.cc

 Description: Extracts and stores HLT objects from specified trigger paths as four-vectors

*/
//
// Original Author:  Mark David John Baber
//         Created:  Sun, 14 Sep 2014 21:25:48 GMT
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "PhysicsTools/PatUtils/interface/TriggerHelper.h"
#include "DataFormats/PatCandidates/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObject.h"



class HLTJetProducer : public edm::EDProducer {
   public:
      explicit HLTJetProducer(const edm::ParameterSet&);
      ~HLTJetProducer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

   private:
      virtual void beginJob() override;
      virtual void produce(edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;

  // ----------member data ---------------------------
  
  edm::InputTag HLTResultsTag;
  double jetMinPT;
};


//
// constructors and destructor
//
HLTJetProducer::HLTJetProducer(const edm::ParameterSet& iConfig)
{

  HLTResultsTag = iConfig.getUntrackedParameter("HLTResults", edm::InputTag("TriggerResults","HLT"));
  jetMinPT      = iConfig.getParameter<double>("JetMinPT");

  produces <std::vector<math::PtEtaPhiELorentzVector> >("HLTCaloJet20");
  produces <std::vector<math::PtEtaPhiELorentzVector> >("HLTPFJet20");
  produces <std::vector<math::PtEtaPhiELorentzVector> >("HLTPFNoPUJet20");
  
}


HLTJetProducer::~HLTJetProducer()
{
}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
HLTJetProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  std::auto_ptr<std::vector<math::PtEtaPhiELorentzVector> > outputCaloJet20_L1Jet(   new std::vector<math::PtEtaPhiELorentzVector>() );
  std::auto_ptr<std::vector<math::PtEtaPhiELorentzVector> > outputPFJet20_L1Jet(     new std::vector<math::PtEtaPhiELorentzVector>() );
  std::auto_ptr<std::vector<math::PtEtaPhiELorentzVector> > outputPFNoPUJet20_L1Jet( new std::vector<math::PtEtaPhiELorentzVector>() );
  
  edm::Handle<edm::TriggerResults> hltresults;
  iEvent.getByLabel(HLTResultsTag, hltresults);
  
  // Get the PAT TriggerEvent
  edm::Handle< pat::TriggerEvent > triggerEvent;
  iEvent.getByLabel( "patTriggerEvent", triggerEvent );
  
  // Get a vector of all HLT paths
  std::vector<pat::TriggerPath> const* paths = triggerEvent->paths();
  
  // Find the full label of the chosen HLT path (i.e. with the version number)
  std::string full_name;
  for (unsigned i = 0; i < paths->size(); ++i) {
    std::string name = paths->at(i).name();
    //    std::cout << name << "\n";

    if ( (name.find("HLT_Calojet20_AllJets_v1")   != name.npos) || 
	 (name.find("HLT_PFjet20_AllJets_v1")     != name.npos) || 
	 (name.find("HLT_PFjetNoPU20_AllJets_v1") != name.npos)) {
      full_name = name;

      // Store in the respective jet collection
      if      (name.find("HLT_Calojet20_AllJets_v1")   != name.npos){
	// Check if trigger fired
	if (paths->at(i).wasAccept()){

	  // Get a vector of the objects used in the chosen path
	  pat::TriggerObjectRefVector pathObjects = triggerEvent->pathObjects(full_name, false);
	  for (unsigned j = 0; j < pathObjects.size(); ++j) {
	    
	    // Extract the 4-momenta of the pathObject
	    if ( pathObjects[j]->pt() < jetMinPT ){ break; }
	    math::PtEtaPhiELorentzVector tempJet;
	    tempJet.SetCoordinates( pathObjects[j]->pt(), pathObjects[j]->eta(), pathObjects[j]->phi(), pathObjects[j]->energy() );
	    
	    outputCaloJet20_L1Jet->push_back( tempJet );
	  } // End pathObject loop
	}
      }
      else if (name.find("HLT_PFjet20_AllJets_v1")     != name.npos){
	// Check if trigger fired
	if (paths->at(i).wasAccept()){

	  // Get a vector of the objects used in the chosen path
	  pat::TriggerObjectRefVector pathObjects = triggerEvent->pathObjects(full_name, false);
	  for (unsigned j = 0; j < pathObjects.size(); ++j) {
	    
	    // Extract the 4-momenta of the pathObject
	    if ( pathObjects[j]->pt() < jetMinPT ){ break; }
	    math::PtEtaPhiELorentzVector tempJet;
	    tempJet.SetCoordinates( pathObjects[j]->pt(), pathObjects[j]->eta(), pathObjects[j]->phi(), pathObjects[j]->energy() );
	    
            outputPFJet20_L1Jet->push_back( tempJet );
	  } // End pathObject loop
	}
      }
      else if (name.find("HLT_PFjetNoPU20_AllJets_v1") != name.npos){
	// Check if trigger fired
	if (paths->at(i).wasAccept()){

	  // Get a vector of the objects used in the chosen path
	  pat::TriggerObjectRefVector pathObjects = triggerEvent->pathObjects(full_name, false);
	  for (unsigned j = 0; j < pathObjects.size(); ++j) {
	    
	    // Extract the 4-momenta of the pathObject
	    if ( pathObjects[j]->pt() < jetMinPT ){ break; }
	    math::PtEtaPhiELorentzVector tempJet;
	    tempJet.SetCoordinates( pathObjects[j]->pt(), pathObjects[j]->eta(), pathObjects[j]->phi(), pathObjects[j]->energy() );
	    
	    outputPFNoPUJet20_L1Jet->push_back( tempJet );
	  } // End pathObject loop
	}
      }

      // // Check if trigger fired
      // if (paths->at(i).wasAccept()){

      // 	// Get a vector of the objects used in the chosen path
      // 	pat::TriggerObjectRefVector pathObjects = triggerEvent->pathObjects(full_name, false);
      // 	for (unsigned j = 0; j < pathObjects.size(); ++j) {

      // 	  // std::cout << "obj #" << j << ":\n"
      // 	  // 	    << "\tpt = "  << pathObjects[j]->pt()  << "\t eta = " << pathObjects[j]->eta() 
      // 	  // 	    << "\tphi = " << pathObjects[j]->phi() << "\te = "    << pathObjects[j]->energy() << "\n";
  

      // 	  // Extract the 4-momenta of the pathObject
      // 	  if ( pathObjects[j]->pt() < jetMinPT ){ break; }
      // 	  math::PtEtaPhiELorentzVector tempJet;
      // 	  tempJet.SetCoordinates( pathObjects[j]->pt(), pathObjects[j]->eta(), pathObjects[j]->phi(), pathObjects[j]->energy() );

      // 	  // Store in the respective jet collection
      // 	  if      (name.find("HLT_Calojet20_AllJets_v1")   != name.npos){
      // 	    outputCaloJet20_L1Jet->push_back( tempJet );
      // 	  }
      // 	  else if (name.find("HLT_PFjet20_AllJets_v1")     != name.npos){
      //       outputPFJet20_L1Jet->push_back( tempJet );
      //     }
      // 	  else if (name.find("HLT_PFjetNoPU20_AllJets_v1") != name.npos){
      //       outputPFNoPUJet20_L1Jet->push_back( tempJet );
      //     }
      // 	  else{
      // 	    //	    std::cout << "Error: No collection matched the HLT path: " << name << "\n";
      // 	  }
  
      // 	} // End pathObject loop

      // }

    }
  }

  iEvent.put( outputCaloJet20_L1Jet,   "HLTCaloJet20");
  iEvent.put( outputPFJet20_L1Jet,     "HLTPFJet20");
  iEvent.put( outputPFNoPUJet20_L1Jet, "HLTPFNoPUJet20");
 
}

// ------------ method called once each job just before starting event loop  ------------
void 
HLTJetProducer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
HLTJetProducer::endJob() {
}
 
// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTJetProducer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTJetProducer);
