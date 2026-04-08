#ifndef _WeightFuncs_h_
#define _WeightFuncs_h_

namespace weight {

///////////////////////////////////////////////////////////////////////
// Calculate the variance in momentum of a group of particles
// used to quantify how unevenly the momentum is shared 

double Asymmetry(const std::vector<TVector3>& p_v)
{

  if(!p_v.size()) return -1;

  double sum = 0.0;
  for(const auto p : p_v) sum += p.Mag(); 
  sum /= p_v.size();

  double var = 0.0;
  for(const auto p : p_v) var += pow(p.Mag() - sum,2); 
  var /= p_v.size();

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

  return ang/3.142;

}

///////////////////////////////////////////////////////////////////////

}

#endif
