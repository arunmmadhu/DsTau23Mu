#include "DsTau23Mu/T3MNtuple/interface/DataMCType.h"
#include "DsTau23Mu/T3MNtuple/interface/PDGInfo.h"
#include <iostream>
#include <cstdlib>

TString DataMCType::type="unknown";


DataMCType::DataMCType(){
}

DataMCType::~DataMCType(){
}

unsigned int DataMCType::GetType(TString name){
  name.ToLower();
  StoreType(name);
  if(name=="data")                      return Data;
  if(name=="ds_tau")                    return Ds_Tau;
  if(name=="b_tau")                     return B_Tau;
  if(name=="ds_phipi")                  return Ds_PhiPi;
  if(name=="b0_tau")                    return B0_Tau;
  if(name=="bp_tau")                    return Bp_Tau;
  if(name=="tau_3mu2nu")                return Tau_3Mu2Nu;
  if(name=="ds_muetamumugamma")         return Ds_MuEtaMuMuGamma;
  if(name=="ds_muetamumu")              return Ds_MuEtaMuMu;
  if(name=="ds_muetamumugammapi0")      return Ds_MuEtaMuMuGammaPi0;
  if(name=="ds_muetaprimemumugamma")    return Ds_MuEtaPrimeMuMuGamma;
  if(name=="ds_muomegamumupi0")         return Ds_MuOmegaMuMuPi0;
  if(name=="ds_muphimumugamma")         return Ds_MuPhiMuMuGamma;
  if(name=="ds_muphimumu")              return Ds_MuPhiMuMu;
  if(name=="ds_muphimumu_mumugamma")    return Ds_MuPhiMuMu_MuMuGamma;
  if(name=="ds_murhomumu")              return Ds_MuRhoMuMu;
  if(name=="ds_muomegamumu")            return Ds_MuOmegaMuMu;
  if(name=="d_muomegamumu_mumupi0")     return D_MuOmegaMuMu_MuMuPi0;
  if(name=="bbcc2mu")                   return bbcc2mu;
  if(name=="bbcc3mu")                   return bbcc3mu;
  if(name=="dy_ll")                     return dy_ll;


  std::cout << "ERROR: Data/MC Type " << name << " UNKNOWN!!!! " << std::endl;
  return unknown;
}


bool DataMCType::isSignalParticle(int pdg_id){
  unsigned int pdgid=abs(pdg_id);
  if(pdgid==PDGInfo::Ds_plus || pdgid==PDGInfo::B_plus || pdgid==PDGInfo::B_0 ){
    return true;
  }
  return false; 
}
