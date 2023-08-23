#include "DsTau23Mu/T3MNtuple/interface/DataMCType.h"
#include "DsTau23Mu/T3MNtuple/interface/PDGInfo.h"
#include "DsTau23Mu/T3MNtuple/interface/TauDecay.h"
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
  if(name=="dp_tau")                    return Dp_Tau;
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
  if(name=="dy_ll")                     return DY_ll;
  if(name=="dy_tautau")                 return DY_tautau;
  if(name=="ww")                        return WW;
  if(name=="wgl")                       return WGl;
  if(name=="zgmumu")                    return ZGmumu;
  if(name=="zgll")                      return ZGll;
  if(name=="zz4l")                      return ZZ4l;
  if(name=="wz3l")                      return WZ3l;
  if(name=="zz4l_1toinf")               return ZZ4l_1toInf;
  if(name=="zz4l_powheg")               return ZZ4l_powheg;
  if(name=="zz4l_amcatnlo")             return ZZ4l_amcatnlo;
  if(name=="wz3l_min0p1")               return WZ3l_Min0p1;
  if(name=="qcd_15")                    return QCD_15;
  if(name=="qcd_20")                    return QCD_20;


  if(name=="z2tautau_tau3mu")           return z2tautau_tau3mu;
  if(name=="w_tau3mu")                  return w_tau3mu;

  std::cout << "ERROR: Data/MC Type " << name << " UNKNOWN!!!! " << std::endl;
  return unknown;
}


unsigned int DataMCType::SignalCode(unsigned int type,unsigned int JAK_ID1, unsigned int nprong1,unsigned int JAK_ID2, unsigned int nprong2){
  if(type==Data)return type;

  if(JAK_ID1==TauDecay::JAK_3MUON && nprong1==3 && (JAK_ID2==TauDecay::JAK_MUON && nprong2==1))   return type*1000 + JAK_ID1*10 + JAK_ID2;
  if(JAK_ID2==TauDecay::JAK_3MUON && nprong2==3 && (JAK_ID1==TauDecay::JAK_MUON && nprong1==1))   return type*1000 + JAK_ID2*10 + JAK_ID1;

  if(JAK_ID1==TauDecay::JAK_3MUON && nprong1==3 && (JAK_ID2==TauDecay::JAK_ELECTRON && nprong2==1))   return type*1000 + JAK_ID1*10 + JAK_ID2;
  if(JAK_ID2==TauDecay::JAK_3MUON && nprong2==3 && (JAK_ID1==TauDecay::JAK_ELECTRON && nprong1==1))   return type*1000 + JAK_ID2*10 + JAK_ID1;


  if(JAK_ID1==TauDecay::JAK_3MUON && nprong1==3 && (JAK_ID2!=TauDecay::JAK_ELECTRON && JAK_ID2!=TauDecay::JAK_MUON))
    {
      JAK_ID2 = TauDecay::JAK_PION;  //  use JAK_ID == 3 for any hadronic tau decay
      return type*1000 + JAK_ID1*10 + JAK_ID2;
    }
  if(JAK_ID2==TauDecay::JAK_3MUON && nprong2==3 && (JAK_ID1!=TauDecay::JAK_ELECTRON && JAK_ID1!=TauDecay::JAK_MUON))
    {
      JAK_ID1 = TauDecay::JAK_PION; //  use JAK_ID == 3 for any hadronic tau decay
      return type*1000 + JAK_ID2*10 + JAK_ID1;
    }

  ////////////////////////////////////////////////////////////////////////////////////
  //
  //   For the decay Z-> tau (3mu) tau(X) it returns: 210231 (IF X = electron);
  //                                                  210232 (IF X = muon    );
  //                                                  210233 (IF X = hadrons );
  //


  return type;
}



bool DataMCType::isSignalParticle(int pdg_id){
  unsigned int pdgid=abs(pdg_id);
  if(pdgid==PDGInfo::Ds_plus || pdgid==PDGInfo::B_plus || pdgid==PDGInfo::B_0  || pdgid ==PDGInfo::Z0  || pdgid == PDGIngo::W_plus || pdgid == PDGIngo::W_minus){
    return true;
  }
  return false; 
}

bool DataMCType::DecodeDecayForSignalParticle(int pdg_id){
  unsigned int pdgid=abs(pdg_id);
  if( pdgid ==PDGInfo::Z0 ){
    return true;
  }
  return false; 
}
