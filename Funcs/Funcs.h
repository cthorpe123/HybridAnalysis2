#ifndef _Funcs_h_
#define _Funcs_h_

#include <vector>

///////////////////////////////////////////////////////////////////////
// Fiducial volume definitions

const double FVxmin = 0.0;
const double FVxmax = 256.35;
const double FVymin = -115.53;
const double FVymax = 117.47;
const double FVzmin = 0.1;
const double FVzmax = 1036.9;
const double deadzmin = 675.1;
const double deadzmax = 775.1;

bool inActiveTPC(double x, double y, double z){
  if(x > FVxmax || x < FVxmin) return false;
  if(y > FVymax || y < FVymin) return false;
  if(z > FVzmax || z < FVzmin) return false;
  if(z > deadzmin && z < deadzmax) return false;
  return true;
}

const double padding = 10;
bool isContained(double x, double y, double z){
  if(x > FVxmax - padding || x < FVxmin + padding) return false;
  if(y > FVymax - padding || y < FVymin + padding) return false;
  if(z > FVzmax - padding || z < FVzmin + padding) return false;
  return true;
}

///////////////////////////////////////////////////////////////////////
// Particle thresholds/kinematics 

const double Mp = 0.935;
const double Mn = 0.938;
const double Eb = 0.03;
const double ml = 0.106;
const double mpi = 0.140;
const double mpi0 = 0.135;
const double MA = 22*Mn + 18*Mp - 0.34381;
const double MA1 = MA - Mn;

std::map<int,std::pair<double,double>> thresholds = {
  {13,{0.1,100.0}},
  {2212,{0.3,5.0}},
  {211,{0.1,5.0}},
  {111,{0.0,5.0}}
};

///////////////////////////////////////////////////////////////////////
// Take 2D hist and normalise each vertical strip to 1 

void Normalise(TH2D* h){
  for(int i_x=1;i_x<h->GetNbinsX()+1;i_x++){
    double content = 0.0;
    for(int i_y=1;i_y<h->GetNbinsY()+1;i_y++) content += h->GetBinContent(i_x,i_y);
    for(int i_y=1;i_y<h->GetNbinsY()+1;i_y++){
      h->SetBinContent(i_x,i_y,h->GetBinContent(i_x,i_y)/content);
      if(std::isnan(h->GetBinContent(i_x,i_y))) h->SetBinContent(i_x,i_y,0.0);
    }
  }
}

///////////////////////////////////////////////////////////////////////
// Take 1D hist divide by bin width 

void DivideByBinWidth(TH1D* h) {
  int NBins = h->GetXaxis()->GetNbins();
  for (int i=1;i<NBins+1;i++){
    double CurrentEntry = h->GetBinContent(i);
    double NewEntry = CurrentEntry / h->GetBinWidth(i);
    double CurrentError = h->GetBinError(i);
    double NewError = CurrentError / h->GetBinWidth(i);
    h->SetBinContent(i,NewEntry); 
    h->SetBinError(i,NewError); 
  }
}

///////////////////////////////////////////////////////////////////////
// Calculate bias/variance in reconstructed variable afo true variable 
// given their joint dist 

void GetBiasVariance(const TH2D* h_Data,TH1D*& h_Bias,TH1D*& h_Variance){

  int nbins_x = h_Data->GetNbinsX();
  int nbins_y = h_Data->GetNbinsY();

  for(int i=1;i<nbins_x+1;i++){
    double mean = 0.0;
    double events = 0.0;
    for(int j=1;j<nbins_y+1;j++){
      mean += h_Data->GetYaxis()->GetBinCenter(j)*h_Data->GetBinContent(i,j);
      events += h_Data->GetBinContent(i,j);
    }

    mean /= events;

    h_Bias->SetBinContent(i,(mean - h_Bias->GetBinCenter(i))/h_Bias->GetBinCenter(i));
    if(events == 0.0) h_Bias->SetBinContent(i,0);

    double var = 0.0; 
    for(int j=1;j<nbins_y+1;j++){
      var += (h_Data->GetYaxis()->GetBinCenter(j) - mean)*(h_Data->GetYaxis()->GetBinCenter(j) - mean)*h_Data->GetBinContent(i,j);
    }

    var /= events;    

    h_Variance->SetBinContent(i,var/mean/mean);
    if(events == 0.0) h_Variance->SetBinContent(i,0);

  }

}

///////////////////////////////////////////////////////////////////////
// TLorentzVector manipulation funcs

void SortTLorentzVector(std::vector<TLorentzVector>& p){
  std::sort(p.begin(),p.end(),[](const TLorentzVector &a, const TLorentzVector &b)
      { 
      return a.Vect().Mag() > b.Vect().Mag(); 
      });
}

TLorentzVector SumTLorentzVector(std::vector<TLorentzVector> p_v){
  TLorentzVector p_tot(0,0,0,0);
  for(TLorentzVector p : p_v) p_tot += p;
  return p_tot;
}

//


#endif
