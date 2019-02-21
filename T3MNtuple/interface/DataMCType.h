#ifndef DataMCType_h
#define DataMCType_h

#include "TString.h"

class DataMCType{
 public:
  enum Type {Data=1,
	     Ds_Tau=10,
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
