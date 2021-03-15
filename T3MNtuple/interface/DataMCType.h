#ifndef DataMCType_h
#define DataMCType_h

#include "TString.h"

class DataMCType{
 public:
  enum Type {Data=1,
	     Ds_PhiPi=30,
	     Ds_Tau=40,
	     B_Tau=50,
	     B0_Tau=60,
	     Bp_Tau=90,
	     bbcc2mu=119,
	     bbcc3mu=120,
	     Tau_3Mu2Nu=121,
	     Ds_MuEtaMuMuGamma=122,
	     Ds_MuEtaMuMu=123,
	     Ds_MuEtaMuMuGammaPi0=124,
	     Ds_MuEtaPrimeMuMuGamma=128,
	     Ds_MuOmegaMuMuPi0=132,
	     Ds_MuPhiMuMuGamma=142,
	     unknown=999
	     // add here more MC types
  };

  DataMCType();
  ~DataMCType();

  unsigned int GetType(TString name);
  bool isSignalParticle(int pdg_id);
  void StoreType(TString t){type=t;}
  unsigned int GetType(){return GetType(type);}

 private:
  static TString type;

};
#endif
