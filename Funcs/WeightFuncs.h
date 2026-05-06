#ifndef _WeightFuncs_h_
#define _WeightFuncs_h_

#include "Funcs.h"
#include "BranchList.h"

namespace weight {

std::vector<double> alphas_LH_Bias;
std::vector<double> betas_LH_Bias;

std::vector<double> alphas_RH_Bias;
std::vector<double> betas_RH_Bias;

std::vector<double> alphas_Center_Gather;
std::vector<double> betas_Center_Gather;

std::vector<double> alphas_Center_Spread;
std::vector<double> betas_Center_Spread;

///////////////////////////////////////////////////////////////////////
// Generate alpha and beta parameter vectors with n uniformly-spaced
// points between the same limits as the hardcoded 5-point defaults above.
// Returns a map keyed by the same names as the const vectors.

std::map<std::string,std::vector<double>> MakeBetaParams(int n){

  auto linspace = [](double lo, double hi, int n){
    std::vector<double> v(n);
    for(int i=0;i<n;i++) v[i] = lo + i*(hi-lo)/(n-1);
    return v;
  };

  std::map<std::string,std::vector<double>> m;
  m["alphas_LH_Bias"]       = linspace(0.5, 1.3, n);
  m["betas_LH_Bias"]        = std::vector<double>(n, 1.0);
  m["alphas_RH_Bias"]       = std::vector<double>(n, 1.0);
  m["betas_RH_Bias"]        = linspace(0.5, 1.3, n);
  m["alphas_Center_Gather"] = linspace(1.0, 5.0, n);
  m["betas_Center_Gather"]  = linspace(1.0, 5.0, n);
  m["alphas_Center_Spread"] = linspace(1.0, 0.5, n);
  m["betas_Center_Spread"]  = linspace(1.0, 0.5, n);
  return m;

}

///////////////////////////////////////////////////////////////////////

double Asymmetry(const std::vector<TLorentzVector>* p_v)
{

  if(p_v->size() < 2) return 0;

  double leading = p_v->at(0).Vect().Mag()*(p_v->size()-1);
  double subleading = 0.0;
  for(size_t i_p=1;i_p<p_v->size();i_p++) subleading += p_v->at(i_p).Vect().Mag();

  return (leading - subleading)/(leading+subleading);

}

///////////////////////////////////////////////////////////////////////

double Asymmetry2(const std::vector<TLorentzVector> p1_v,const std::vector<TLorentzVector> p2_v)
{

  if(!p1_v.size() || !p2_v.size()) return 0;

  double leading = 0.0;
  for(size_t i_p=0;i_p<p1_v.size();i_p++) leading += p1_v.at(i_p).Vect().Mag(); 
  leading /= p1_v.size();

  double subleading = 0.0;
  for(size_t i_p=0;i_p<p2_v.size();i_p++) subleading += p2_v.at(i_p).Vect().Mag(); 
  subleading /= p2_v.size();

  return (leading - subleading)/(leading+subleading);

}

///////////////////////////////////////////////////////////////////////

double Asymmetry3(const std::vector<TLorentzVector> p1_v,const std::vector<TLorentzVector> p2_v)
{

  if(!p1_v.size() || !p2_v.size()) return 0;

  double leading = 0.0;
  for(size_t i_p=0;i_p<p1_v.size();i_p++) leading += p1_v.at(i_p).E() - p1_v.at(i_p).M(); 
  leading /= p1_v.size();

  double subleading = 0.0;
  for(size_t i_p=0;i_p<p2_v.size();i_p++) subleading += p2_v.at(i_p).E() - p2_v.at(i_p).M(); 
  subleading /= p2_v.size();

  return (leading - subleading)/(leading+subleading);

}

///////////////////////////////////////////////////////////////////////
// Calculate the variance in KE of a group of particles
// used to quantify how unevenly the momentum is shared 

double AsymmetryKE(const std::vector<double>& ke_v)
{

  if(!ke_v.size()) return -1;

  double sum = 0.0;
  for(const double k : ke_v) sum += k; 
  sum /= ke_v.size();

  double var = 0.0;
  for(const double k : ke_v) var += pow(k - sum,2); 
  var /= sum*ke_v.size();

  return var;

}

///////////////////////////////////////////////////////////////////////
// Calculate the total momentum vector of a group of angles, then a 
// weighted mean of angles individual paritcles make with that vector
// Use to quantify how coliniear paritcles are 

double Cone(const std::vector<TVector3>& p_v)
{

  if(p_v.size() < 2) return -1;

  TVector3 sum(0,0,0);
  double sum_m = 0.0;
  for(const auto p : p_v){
    sum += p; 
    sum_m += p.Mag();
  }

  double ang = 0.0;
  for(const auto p : p_v) ang += sum.Angle(p)*p.Mag(); 
  ang *= 1.0/sum_m;

  return 2*ang/3.142;

}

///////////////////////////////////////////////////////////////////////
// Set up all of the pointers to the different weight functions here

const int spline_pts = 10; // Number of points to use in the spline reweighting
const double max_weight = 5.0; // Max weight to apply in the spline reweighting

// Function pointers for the various reweighters 
std::map<std::string,std::vector<double>(*)()> r_m;
const std::vector<std::string> weight_func_labels = {"1p1piOpeningAngle_Center_Gather","1p1piOpeningAngle_Center_Spread","1p1piOpeningAngle_LH_Bias","1p1piOpeningAngle_RH_Bias","2pAsym_Center_Gather","2pAsym_Center_Spread","2pAsym_LH_Bias","2pAsym_RH_Bias","2pOpeningAngle_Center_Gather","2pOpeningAngle_Center_Spread","2pOpeningAngle_LH_Bias","2pOpeningAngle_RH_Bias","Extra1P","Extra2G","Extra2P","Extra2Pi","Extra3P","ExtraG","ExtraNP","ExtraP","ExtraPi","LeadPionE_Center_Gather","LeadPionE_Center_Spread","LeadPionE_LH_Bias","LeadPionE_RH_Bias","LeadProtonKEShape_Center_Gather","LeadProtonKEShape_Center_Spread","LeadProtonKEShape_LH_Bias","LeadProtonKEShape_RH_Bias","MuonAngleShape_Center_Gather","MuonAngleShape_Center_Spread","MuonAngleShape_LH_Bias","MuonAngleShape_RH_Bias","MuonMomShape_Center_Gather","MuonMomShape_Center_Spread","MuonMomShape_LH_Bias","MuonMomShape_RH_Bias","MuonProtonOpeningAngle_Center_Gather","MuonProtonOpeningAngle_Center_Spread","MuonProtonOpeningAngle_LH_Bias","MuonProtonOpeningAngle_RH_Bias"};

void SetWeightFuncs()
{
  auto params = MakeBetaParams(spline_pts);
  alphas_LH_Bias       = params.at("alphas_LH_Bias");
  betas_LH_Bias        = params.at("betas_LH_Bias");
  alphas_RH_Bias       = params.at("alphas_RH_Bias");
  betas_RH_Bias        = params.at("betas_RH_Bias");
  alphas_Center_Gather = params.at("alphas_Center_Gather");
  betas_Center_Gather  = params.at("betas_Center_Gather");
  alphas_Center_Spread = params.at("alphas_Center_Spread");
  betas_Center_Spread  = params.at("betas_Center_Spread");

  auto f_ExtraPi = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(npi_t > 0 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraPi",f_ExtraPi); 

  auto f_Extra2Pi = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(npi_t > 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra2Pi",f_Extra2Pi); 

  auto f_ExtraP = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t > 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraP",f_ExtraP); 

  auto f_Extra1P = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t == 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra1P",f_Extra1P); 

  auto f_Extra2P = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t == 2 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra2P",f_Extra2P); 

  auto f_Extra3P = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t == 3 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra3P",f_Extra3P); 

  auto f_ExtraNP = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t > 3 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraNP",f_ExtraNP); 

  auto f_ExtraG = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nsh_t > 0 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraG",f_ExtraG); 

  auto f_Extra2G = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nsh_t > 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra2G",f_Extra2G); 

  auto f_MuonAngleShape_LH_Bias = [](){ 
    double val = (muon_mom_t->CosTheta() + 1)/2.0; // Map from [-1,1] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("MuonAngleShape_LH_Bias",f_MuonAngleShape_LH_Bias);

  auto f_MuonAngleShape_RH_Bias = [](){ 
    double val = (muon_mom_t->CosTheta() + 1)/2.0; // Map from [-1,1] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("MuonAngleShape_RH_Bias",f_MuonAngleShape_RH_Bias);

  auto f_MuonAngleShape_Center_Gather = [](){ 
    double val = (muon_mom_t->CosTheta() + 1)/2.0; // Map from [-1,1] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("MuonAngleShape_Center_Gather",f_MuonAngleShape_Center_Gather);

  auto f_MuonAngleShape_Center_Spread = [](){ 
    double val = (muon_mom_t->CosTheta() + 1)/2.0; // Map from [-1,1] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("MuonAngleShape_Center_Spread",f_MuonAngleShape_Center_Spread);

  auto f_MuonMomShape_LH_Bias = [](){ 
    double val = muon_mom_t->Mag()/2.0; // Map from [0,2] to [0,1]
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("MuonMomShape_LH_Bias",f_MuonMomShape_LH_Bias);

  auto f_MuonMomShape_RH_Bias = [](){ 
    double val = muon_mom_t->Mag()/2.0; // Map from [0,2] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("MuonMomShape_RH_Bias",f_MuonMomShape_RH_Bias);

  auto f_MuonMomShape_Center_Gather = [](){ 
    double val = muon_mom_t->Mag()/2.0; // Map from [0,2] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("MuonMomShape_Center_Gather",f_MuonMomShape_Center_Gather);  

  auto f_MuonMomShape_Center_Spread = [](){ 
    double val = muon_mom_t->Mag()/2.0; // Map from [0,2] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || val > 1 || val < 0) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("MuonMomShape_Center_Spread",f_MuonMomShape_Center_Spread);

 auto f_LeadProtonKEShape_LH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t) return x;
    double val = (protons_t->at(0).E() - Mp)/0.6; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("LeadProtonKEShape_LH_Bias",f_LeadProtonKEShape_LH_Bias);

  auto f_LeadProtonKEShape_RH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t) return x;
    double val = (protons_t->at(0).E() - Mp)/0.6; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("LeadProtonKEShape_RH_Bias",f_LeadProtonKEShape_RH_Bias);

  auto f_LeadProtonKEShape_Center_Gather = [](){
    std::vector<double> x(spline_pts,1.0);  
    if(!is_signal_t) return x;
    double val = (protons_t->at(0).E() - Mp)/0.6; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("LeadProtonKEShape_Center_Gather",f_LeadProtonKEShape_Center_Gather);  

  auto f_LeadProtonKEShape_Center_Spread = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t) return x;
    double val = (protons_t->at(0).E() - Mp)/0.6; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("LeadProtonKEShape_Center_Spread",f_LeadProtonKEShape_Center_Spread);

  auto f_2pOpeningAngle_LH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(protons_t->at(1).Vect());
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("2pOpeningAngle_LH_Bias",f_2pOpeningAngle_LH_Bias);

  auto f_2pOpeningAngle_RH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(protons_t->at(1).Vect());
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("2pOpeningAngle_RH_Bias",f_2pOpeningAngle_RH_Bias);

  auto f_2pOpeningAngle_Center_Gather = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(protons_t->at(1).Vect());
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("2pOpeningAngle_Center_Gather",f_2pOpeningAngle_Center_Gather);

  auto f_2pOpeningAngle_Center_Spread= [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(protons_t->at(1).Vect());
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("2pOpeningAngle_Center_Spread",f_2pOpeningAngle_Center_Spread);
  
  auto f_2pAsym_LH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = Asymmetry3({protons_t->at(0)},{protons_t->at(1)})/0.8;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("2pAsym_LH_Bias",f_2pAsym_LH_Bias);

  auto f_2pAsym_RH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = Asymmetry3({protons_t->at(0)},{protons_t->at(1)})/0.8;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("2pAsym_RH_Bias",f_2pAsym_RH_Bias);

  auto f_2pAsym_Center_Gather = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = Asymmetry3({protons_t->at(0)},{protons_t->at(1)})/0.8;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("2pAsym_Center_Gather",f_2pAsym_Center_Gather);

  auto f_2pAsym_Center_Spread = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || nprot_t != 2) return x;
    double val = Asymmetry3({protons_t->at(0)},{protons_t->at(1)})/0.8;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("2pAsym_Center_Spread",f_2pAsym_Center_Spread);

  auto f_LeadPionE_LH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || npi_t == 0) return x;
    double val = pions_t->at(0).E()/0.5; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("LeadPionE_LH_Bias",f_LeadPionE_LH_Bias);

  auto f_LeadPionE_RH_Bias = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || npi_t == 0) return x;
    double val = pions_t->at(0).E()/0.5; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("LeadPionE_RH_Bias",f_LeadPionE_RH_Bias);

  auto f_LeadPionE_Center_Gather = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || npi_t == 0) return x;
    double val = pions_t->at(0).E()/0.5; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("LeadPionE_Center_Gather",f_LeadPionE_Center_Gather);

  auto f_LeadPionE_Center_Spread = [](){ 
    std::vector<double> x(spline_pts,1.0); 
    if(!is_signal_t || npi_t == 0) return x;
    double val = pions_t->at(0).E()/0.5; // Map from [0,0.5] to [0,1]
    if(val < 0 || val > 1) return x;
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("LeadPionE_Center_Spread",f_LeadPionE_Center_Spread);

  auto f_1p1piOpeningAngle_LH_Bias = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(pions_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("1p1piOpeningAngle_LH_Bias",f_1p1piOpeningAngle_LH_Bias);

  auto f_1p1piOpeningAngle_RH_Bias = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(pions_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("1p1piOpeningAngle_RH_Bias",f_1p1piOpeningAngle_RH_Bias);

  auto f_1p1piOpeningAngle_Center_Gather = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(pions_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("1p1piOpeningAngle_Center_Gather",f_1p1piOpeningAngle_Center_Gather);

  auto f_1p1piOpeningAngle_Center_Spread = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = 1.0/3.142*protons_t->at(0).Vect().Angle(pions_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("1p1piOpeningAngle_Center_Spread",f_1p1piOpeningAngle_Center_Spread);

  auto f_1p1piAsym_LH_Bias = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = (Asymmetry3({protons_t->at(0)},{pions_t->at(0)}) + 1) / 2;
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("1p1piAsym_LH_Bias",f_1p1piAsym_LH_Bias);

  auto f_1p1piAsym_RH_Bias = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = (Asymmetry3({protons_t->at(0)},{pions_t->at(0)}) + 1) / 2;
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("1p1piAsym_RH_Bias",f_1p1piAsym_RH_Bias);

  auto f_1p1piAsym_Center_Gather = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = (Asymmetry3({protons_t->at(0)},{pions_t->at(0)}) + 1) / 2;
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("1p1piAsym_Center_Gather",f_1p1piAsym_Center_Gather);

  auto f_1p1piAsym_Center_Spread = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t || nprot_t != 1 || npi_t != 1) return x;
    double val = (Asymmetry3({protons_t->at(0)},{pions_t->at(0)}) + 1) / 2;
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("1p1piAsym_Center_Spread",f_1p1piAsym_Center_Spread);

  auto f_MuonProtonOpeningAngle_LH_Bias = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t) return x;
    double val = 1.0/3.142*muon_mom_t->Angle(protons_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_LH_Bias.at(i),betas_LH_Bias.at(i));
    return x;
  };
  r_m.emplace("MuonProtonOpeningAngle_LH_Bias",f_MuonProtonOpeningAngle_LH_Bias);

  auto f_MuonProtonOpeningAngle_RH_Bias = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t) return x;
    double val = 1.0/3.142*muon_mom_t->Angle(protons_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_RH_Bias.at(i),betas_RH_Bias.at(i));
    return x;
  };
  r_m.emplace("MuonProtonOpeningAngle_RH_Bias",f_MuonProtonOpeningAngle_RH_Bias);

  auto f_MuonProtonOpeningAngle_Center_Gather = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t) return x;
    double val = 1.0/3.142*muon_mom_t->Angle(protons_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Gather.at(i),betas_Center_Gather.at(i));
    return x;
  };
  r_m.emplace("MuonProtonOpeningAngle_Center_Gather",f_MuonProtonOpeningAngle_Center_Gather);

  auto f_MuonProtonOpeningAngle_Center_Spread = [](){
    std::vector<double> x(spline_pts,1.0);
    if(!is_signal_t) return x;
    double val = 1.0/3.142*muon_mom_t->Angle(protons_t->at(0).Vect());
    for(int i=0;i<spline_pts;i++)
      x.at(i) *= ROOT::Math::beta_pdf(val,alphas_Center_Spread.at(i),betas_Center_Spread.at(i));
    return x;
  };
  r_m.emplace("MuonProtonOpeningAngle_Center_Spread",f_MuonProtonOpeningAngle_Center_Spread);

}

///////////////////////////////////////////////////////////////////////

void DrawBeta(int n = spline_pts){

  auto params = MakeBetaParams(n);
  const auto& a_LH = params.at("alphas_LH_Bias");
  const auto& b_LH = params.at("betas_LH_Bias");
  const auto& a_RH = params.at("alphas_RH_Bias");
  const auto& b_RH = params.at("betas_RH_Bias");
  const auto& a_CG = params.at("alphas_Center_Gather");
  const auto& b_CG = params.at("betas_Center_Gather");
  const auto& a_CS = params.at("alphas_Center_Spread");
  const auto& b_CS = params.at("betas_Center_Spread");

  std::vector<TH1D*> h_LH_Bias;
  std::vector<TH1D*> h_RH_Bias;
  std::vector<TH1D*> h_Center_Gather;
  std::vector<TH1D*> h_Center_Spread;

  int bins = 200;

  THStack* hs_LH_Bias = new THStack("hss_LH_Bias",";LH Bias;Events");
  THStack* hs_RH_Bias = new THStack("hss_RH_Bias",";RH Bias;Events");
  THStack* hs_Center_Gather = new THStack("hss_Center_Gather",";Center Gather;Events");
  THStack* hs_Center_Spread = new THStack("hss_Center_Spread",";Center Spread;Events");

  for(int i=0;i<n;i++){
      h_LH_Bias.push_back(new TH1D(("h_LH_Bias_"+std::to_string(i)).c_str(),";LH Bias;Events",bins,0.0,1.0));
      h_RH_Bias.push_back(new TH1D(("h_RH_Bias_"+std::to_string(i)).c_str(),";RH Bias;Events",bins,0.0,1.0));
      h_Center_Gather.push_back(new TH1D(("h_Center_Gather_"+std::to_string(i)).c_str(),";Center Gather;Events",bins,0.0,1.0));
      h_Center_Spread.push_back(new TH1D(("h_Center_Spread_"+std::to_string(i)).c_str(),";Center Spread;Events",bins,0.0,1.0));
      hs_LH_Bias->Add(h_LH_Bias.at(i));
      hs_RH_Bias->Add(h_RH_Bias.at(i));
      hs_Center_Gather->Add(h_Center_Gather.at(i));
      hs_Center_Spread->Add(h_Center_Spread.at(i));
  }

  TLegend* leg_LH_Bias = new TLegend(0.7,0.7,0.9,0.9);
  TLegend* leg_RH_Bias = new TLegend(0.7,0.7,0.9,0.9);
  TLegend* leg_Center_Gather = new TLegend(0.7,0.7,0.9,0.9);
  TLegend* leg_Center_Spread = new TLegend(0.7,0.7,0.9,0.9);

  for(int j=0;j<n;j++){

      leg_LH_Bias->AddEntry(h_LH_Bias.at(j),Form("#alpha=%.2f,#beta=%.2f",a_LH.at(j),b_LH.at(j)),"l");
      leg_RH_Bias->AddEntry(h_RH_Bias.at(j),Form("#alpha=%.2f,#beta=%.2f",a_RH.at(j),b_RH.at(j)),"l");
      leg_Center_Gather->AddEntry(h_Center_Gather.at(j),Form("#alpha=%.2f,#beta=%.2f",a_CG.at(j),b_CG.at(j)),"l");
      leg_Center_Spread->AddEntry(h_Center_Spread.at(j),Form("#alpha=%.2f,#beta=%.2f",a_CS.at(j),b_CS.at(j)),"l");

      for(int i=1;i<=bins;i++){
        double x = h_LH_Bias.at(0)->GetBinCenter(i);
        h_LH_Bias.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,a_LH.at(j),b_LH.at(j)));
        h_RH_Bias.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,a_RH.at(j),b_RH.at(j)));
        h_Center_Gather.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,a_CG.at(j),b_CG.at(j)));
        h_Center_Spread.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,a_CS.at(j),b_CS.at(j)));
      }

  }

  TCanvas* c1 = new TCanvas("c1","c1",800,600);

  for(int i=0;i<n;i++){
      h_LH_Bias.at(i)->SetLineColor(i+1);
      h_LH_Bias.at(i)->SetLineWidth(2);
      h_RH_Bias.at(i)->SetLineColor(i+1);
      h_RH_Bias.at(i)->SetLineWidth(2);
      h_Center_Gather.at(i)->SetLineColor(i+1);
      h_Center_Gather.at(i)->SetLineWidth(2);
      h_Center_Spread.at(i)->SetLineColor(i+1);
      h_Center_Spread.at(i)->SetLineWidth(2);
  }

  hs_LH_Bias->Draw("NOSTACK HIST");
  leg_LH_Bias->Draw();
  c1->Print("LH_Bias.png");
  c1->Clear();

  hs_RH_Bias->Draw("NOSTACK HIST");
  leg_RH_Bias->Draw();
  c1->Print("RH_Bias.png");
  c1->Clear();

  hs_Center_Gather->Draw("NOSTACK HIST");
  leg_Center_Gather->Draw();
  c1->Print("Center_Gather.png");
  c1->Clear();

  hs_Center_Spread->Draw("NOSTACK HIST");
  leg_Center_Spread->Draw();
  c1->Print("Center_Spread.png");
  c1->Clear();

}

///////////////////////////////////////////////////////////////////////
// Print the labels of the various groups of weight functions

void PrintWeightFuncLabels(){
  SetWeightFuncs();
  for(const auto &item : r_m)
    std::cout << "\"" << item.first << "\",";
  std::cout << std::endl;
}

///////////////////////////////////////////////////////////////////////

}

#endif
