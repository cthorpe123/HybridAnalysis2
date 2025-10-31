#ifndef _EnergyEstimatorFuncs_h_
#define _EnergyEstimatorFuncs_h_

#include "Funcs.h"

namespace ee {

enum estimators { kMuonKin , kMuonKinWNP , kPeLEELike0Pi , kTotalEDep , kSFMethod , kMAX };
const std::vector<std::string> estimators_str = { "MuonKin" , "MuonKinWNP" , "PeLEELike0Pi"  , "TotalEDep" , "SFMethod" };
const std::vector<std::string> estimators_leg = { "CCQE-like" , "W^{2}" , "Proton-Based"  , "Calorimetric" , "SF" };
const std::vector<int> colors = {kCyan+2,kBlue,kRed,kMagenta,kGreen+1};

double T2KEnergy(const TLorentzVector plepton){
  return (Mp*Mp - (Mn - Eb)*(Mn - Eb) - plepton.M()*plepton.M() + 2*(Mn - Eb)*plepton.E())/(2*(Mn - Eb - plepton.E() + plepton.P()*plepton.Vect().CosTheta()));
}

double T2KEnergyW(const TLorentzVector plepton,const double& W){
  return (W*W - (Mn - Eb)*(Mn - Eb) - plepton.M()*plepton.M() + 2*(Mn - Eb)*plepton.E())/(2*(Mn - Eb - plepton.E() + plepton.P()*plepton.Vect().CosTheta()));
}

double ubooneEnergy(const TLorentzVector plepton,const double& W,const int& nprot){
  return (W*W - nprot*nprot*(Mn - Eb)*(Mn - Eb) - plepton.M()*plepton.M() + nprot*2*(Mn - Eb)*plepton.E())/(2*(nprot*Mn - nprot*Eb - plepton.E() + plepton.P()*plepton.Vect().CosTheta()));
}

double peleeEnergy(const TLorentzVector plepton,const TLorentzVector pproton,const int nprot){
  return plepton.E() + pproton.E() - Eb*nprot - Mp*nprot;
}

double totaledepEnergy(const TLorentzVector plepton,const TLorentzVector pproton,const int nprot,const TLorentzVector ppi,const TLorentzVector ppi0){
  return plepton.E() + pproton.E() - Eb*nprot - Mp*nprot + ppi.E() + ppi0.E();
}

double sfmethodEnergy(const TLorentzVector plepton,const TLorentzVector pproton){
  const TVector3& prot = pproton.Vect();
  const double Ep = pproton.E(); 
  double pt2 = (prot + plepton.Vect()).Perp2();
  double pL = (pow((MA + plepton.Pz() + prot.Z() - plepton.E() - Ep),2) - pt2 - MA1*MA1)/2/(MA + plepton.Pz() + prot.Z() - plepton.E() - Ep); 
  return plepton.Pz() + prot.Z() - pL;

}

double GetEnergy(const TLorentzVector plepton,const double& W,const TLorentzVector pprotons,const int& nprot,const TLorentzVector ppions, const int& npi,const TLorentzVector ppi0s,int npi0,int est){
  if(est < 0 || est >= kMAX) throw std::invalid_argument("GetEnergy: invalid estimator");

  switch(est){
    case kMuonKin: return T2KEnergy(plepton);
    case kMuonKinWNP: if(nprot > 0) return ubooneEnergy(plepton,W,nprot); else return -1;
    case kPeLEELike0Pi: if(npi == 0 && npi0 == 0) return peleeEnergy(plepton,pprotons,nprot); else return -1;
    case kTotalEDep: return totaledepEnergy(plepton,pprotons,nprot,ppions,ppi0s);
    case kSFMethod: if(nprot == 1 && npi == 0 && npi0 == 0) return sfmethodEnergy(plepton,pprotons); else return -1;
    default: return -1;
  }

  return -1;

}

std::vector<double> GetEnergyEst(const TLorentzVector plepton,const double& W,const TLorentzVector pprotons,const int& nprot,const TLorentzVector ppions, const int& npi,const TLorentzVector ppi0s,int npi0){

  std::vector<double> e;

  for(int i_e=0;i_e<kMAX;i_e++)
    e.push_back(GetEnergy(plepton,W,pprotons,nprot,ppions,npi,ppi0s,npi0,i_e));

  return e;

}

}

#endif
