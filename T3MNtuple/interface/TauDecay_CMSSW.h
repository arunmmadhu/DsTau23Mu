
#ifndef TauDecay_CMSSW_h
#define TauDecay_CMSSW_h

#include "DsTau23Mu/T3MNtuple/interface/TauDecay.h"
#include "DsTau23Mu/T3MNtuple/interface/T3MNtuple.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include <SimDataFormats/GeneratorProducts/interface/HepMCProduct.h>
#include <DataFormats/PatCandidates/interface/GenericParticle.h>
#include <DataFormats/PatCandidates/interface/PackedGenParticle.h>

//
// class declaration
//

using namespace reco;
using namespace edm;
using namespace std;




class TauDecay_CMSSW : public TauDecay {
 public:
  TauDecay_CMSSW();
  ~TauDecay_CMSSW();

  //Function to analyze the tau
  bool AnalyzeTau(const reco::GenParticle *Tau,unsigned int &JAK_ID,unsigned int &TauBitMask);
  // Functions to get results
  std::vector<const reco::GenParticle* > Get_TauDecayProducts(){return TauDecayProducts;}
  std::vector<unsigned int> Get_MotherIdx(){return MotherIdx;}
  //  void CheckForSignal(unsigned int &type,edm::Handle <edm::View<reco::GenParticle> > genHandle);
  void CheckForSignal(unsigned int &type, const  Handle<GenParticleCollection>&  genHandle);

 private:
  // recursive function to loop through tau decay products
  void Analyze(const reco::GenParticle *Particle,unsigned int midx);
  void AddPi0Info(const reco::GenParticle *Particle,unsigned int midx);
  //varibles
  std::vector<const reco::GenParticle*> TauDecayProducts;
  std::vector<unsigned int> MotherIdx;
  unsigned int JAK_ID, TauBitMask;

};
#endif
