#ifndef _WeightFuncs_h_
#define _WeightFuncs_h_

namespace weight {

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

}

#endif
