/* \class GenParticleFilter
 *
 * \author Mark Baber, Imperial College London
 * 
 */
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Utilities/interface/InputTag.h"

class GenParticleFilter : public edm::EDFilter {
public:
  GenParticleFilter(const edm::ParameterSet &);
private:
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  edm::InputTag particleCands_;
  double ptMin_;
};

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandAssociation.h"
using namespace edm;
using namespace std;
using namespace reco;

GenParticleFilter::GenParticleFilter( const ParameterSet & cfg ) :
  particleCands_(cfg.getParameter<InputTag>("particleCands")),
  ptMin_(cfg.getParameter<double>("ptMin"))
{
}

bool GenParticleFilter::filter (Event & ev, const EventSetup &) {
  Handle<CandidateCollection> particleCands;
  ev.getByLabel(particleCands_, particleCands);
  unsigned int nPart = particleCands->size();
  if (nPart == 0) return false;
  for(unsigned int i = 0; i < nPart; ++ i) {
    const Candidate & particleCand = (*particleCands)[i];
    double pt = particleCand.pt();
    if (pt < ptMin_) return false;
  }
  return true;
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE( GenParticleFilter );
