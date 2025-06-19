#ifndef DataMCType_h
#define DataMCType_h

#include "TString.h"

class DataMCType{
 public:
  enum Type {Data=1,
	     Ds_PhiPi=30,
	     Ds_Tau=40,
	     Dp_Tau=45,
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
	     Ds_MuPhiMuMu=143,
	     Ds_MuPhiMuMu_MuMuGamma=145,
	     Ds_MuRhoMuMu=153,
	     Ds_MuOmegaMuMu=163,
	     D_MuOmegaMuMu_MuMuPi0=175,
	     bbcc3mu_type1=181,
	     bbcc3mu_type2=182,
	     bbcc3mu_type3=183,
	     bbcc3mu_type5=185,
	     bbcc3mu_type9=189,
	     DY_ll=200,
	     DY_tautau=205,
	     z2tautau_tau3mu=210,
             z2tautau_tau3mu165=251,
             z2tautau_tau3mu170=252,
             z2tautau_tau3mu185=253,
             z2tautau_tau3mu190=254,
             z2tautau_tau3mu195=255,
	     WW=510,
	     WGl=520,
	     ZGmumu=530,
	     ZGll=540,
	     ZZ4l=550,
	     WZ3l=560,
	     ZZ4l_1toInf=570,
             ZZ4l_powheg=571,
             ZZ4l_amcatnlo=572,
             QCD_15=601,
             QCD_20=602,
	     WZ3l_Min0p1=580,
	     ttbar=700,
	     unknown=999
	     // add here more MC types
  };

  DataMCType();
  ~DataMCType();

  unsigned int GetType(TString name);
  unsigned int SignalCode(unsigned int type,unsigned int JAK_ID1, unsigned int nprong1,unsigned int JAK_ID2, unsigned int nprong2);
  bool isSignalParticle(int pdg_id);
  bool DecodeDecayForSignalParticle(int pdg_id);
  void StoreType(TString t){type=t;}
  unsigned int GetType(){return GetType(type);}

 private:
  static TString type;

};
#endif
