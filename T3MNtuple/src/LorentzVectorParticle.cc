#include "DsTau23Mu/T3MNtuple/interface/LorentzVectorParticle.h"


LorentzVectorParticle::LorentzVectorParticle():
  SimpleParticle(TMatrixT<double>(NLorentzandVertexPar,1),TMatrixTSym<double>(NLorentzandVertexPar),0,0,0)
{
}

LorentzVectorParticle::LorentzVectorParticle(const LorentzVectorParticle& other):
  SimpleParticle(other.getParMatrix(), other.getCovMatrix(), other.PDGID(), other.Charge(), other.BField())
{
}


LorentzVectorParticle::LorentzVectorParticle(TMatrixT<double> par_,TMatrixTSym<double> cov_,int pdgid_,double charge_,double b_):
  SimpleParticle(par_,cov_,pdgid_,charge_,b_)
{
}

/*LorentzVectorParticle& LorentzVectorParticle::operator=(const LorentzVectorParticle& other){
  return LorentzVectorParticle(other);
}
*/
TString LorentzVectorParticle::Name(int i){
  if(i==px)  return "px";
  if(i==py)  return "py";
  if(i==pz)  return "pz";
  if(i==m)   return "m";
  if(i==vx)  return "vx";
  if(i==vy)  return "vy";
  if(i==vz)  return "vz";
  return "invalid";
}

double LorentzVectorParticle::Parameter(int i){
  if(i==E)  return sqrt(pow(SimpleParticle::Parameter(m),2.0)+pow(SimpleParticle::Parameter(px),2.0)+pow(SimpleParticle::Parameter(py),2.0)+pow(SimpleParticle::Parameter(pz),2.0));
  if(i==p)  return sqrt(pow(SimpleParticle::Parameter(px),2.0)+pow(SimpleParticle::Parameter(py),2.0)+pow(SimpleParticle::Parameter(pz),2.0));
  if(i==pt) return sqrt(pow(SimpleParticle::Parameter(px),2.0)+pow(SimpleParticle::Parameter(py),2.0));
  return SimpleParticle::Parameter(i);
}

double LorentzVectorParticle::Covariance(int i,int j){
	if(i==E && j==E)
		return pow(Parameter(E), -2) * (
				pow(Parameter(px),2)*SimpleParticle::Covariance(px,px) + pow(Parameter(py),2)*SimpleParticle::Covariance(py,py) + pow(Parameter(pz),2)*SimpleParticle::Covariance(pz,pz) + Parameter(m)*Parameter(m)* SimpleParticle::Covariance(m,m) +
				2*Parameter(px)*Parameter(py)*SimpleParticle::Covariance(px,py) + 2*Parameter(px)*Parameter(pz)*SimpleParticle::Covariance(px,pz) + 2*Parameter(m)*Parameter(px)*SimpleParticle::Covariance(px,m) +
				2*Parameter(py)*Parameter(pz)*SimpleParticle::Covariance(py,pz) + 2*Parameter(m)*Parameter(py)*SimpleParticle::Covariance(py,m) +
				2*Parameter(m)*Parameter(pz)*SimpleParticle::Covariance(pz,m));
	if( i==E && (j==px || j==py || j==pz) )
		return -(SimpleParticle::Covariance(j,m)*Parameter(m) + SimpleParticle::Covariance(j,px)*Parameter(px) +
				 SimpleParticle::Covariance(j,py)*Parameter(py) + SimpleParticle::Covariance(j,pz)*Parameter(pz)) / Parameter(E);
	if (j==E && (i==px || i==py || i==pz))
		return Covariance(j,i); 
	return SimpleParticle::Covariance(i,j);
}  
 
 
