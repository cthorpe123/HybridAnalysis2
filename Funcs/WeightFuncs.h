#ifndef _WeightFuncs_h_
#define _WeightFuncs_h_

#include "BranchList.h"

namespace weight {

const std::vector<double> alphas_LH_Bias = {0.7,0.8,1.0,1.2,1.5};
const std::vector<double> betas_LH_Bias = {1.0,1.0,1.0,1.0,1.0};

const std::vector<double> alphas_RH_Bias = {1.0,1.0,1.0,1.0,1.0};
const std::vector<double> betas_RH_Bias = {0.7,0.8,1.0,1.2,1.5};

const std::vector<double> alphas_Center_Gather = {1.0,1.5,2.0,2.5,3.0};
const std::vector<double> betas_Center_Gather = {1.0,1.5,2.0,2.5,3.0};

const std::vector<double> alphas_Center_Spread = {1.0,0.95,0.9,0.85,0.8};
const std::vector<double> betas_Center_Spread = {1.0,0.95,0.9,0.85,0.8};

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

const int spline_pts = 5; // Number of points to use in the spline reweighting
const double max_weight = 3.0; // Max weight to apply in the spline reweighting

// Function pointers for the various reweighters 
std::map<std::string,std::vector<double>(*)()> r_m;

void SetWeightFuncs()
{   
  // Declare the function pointers here
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

  auto f_MuonAngleShape = [](){ 
    double val = (muon_mom_t->CosTheta() + 1)/2.0; // Map from [-1,1] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= std::min(max_weight, ROOT::Math::beta_pdf(val,i*2.0/spline_pts,0.5));
    return x;
  };
  r_m.emplace("MuonAngleShape",f_MuonAngleShape);

  auto f_MuonAngleShape2 = [](){ 
    double val = (muon_mom_t->CosTheta() + 1)/2.0; // Map from [-1,1] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= std::min(max_weight, ROOT::Math::beta_pdf(val,i*1.5/spline_pts,0.1));
    return x;
  };
  r_m.emplace("MuonAngleShape2",f_MuonAngleShape2);

  auto f_MuonMomShape = [](){ 
    double val = muon_mom_t->Mag()/2.0; // Map from [0,2] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= std::min(max_weight, ROOT::Math::beta_pdf(val,0.5,i*2.0/spline_pts));
    return x;
  };
  r_m.emplace("MuonMomShape",f_MuonMomShape);

  auto f_MuonMomShape2 = [](){ 
    double val = muon_mom_t->Mag()/2.0; // Map from [0,2] to [0,1]
    std::vector<double> x(spline_pts,1.0); 
    for(int i=0;i<spline_pts;i++) 
      x.at(i) *= std::min(max_weight, ROOT::Math::beta_pdf(val,1.0,i*1.5/spline_pts));
    return x;
  };
  r_m.emplace("MuonMomShape2",f_MuonMomShape2);

}

///////////////////////////////////////////////////////////////////////

void DrawBeta(){

  std::vector<TH1D*> h_LH_Bias;
  std::vector<TH1D*> h_RH_Bias;
  std::vector<TH1D*> h_Center_Gather;
  std::vector<TH1D*> h_Center_Spread;

  int bins = 200;

  THStack* hs_LH_Bias = new THStack("hss_LH_Bias",";LH Bias;Events");
  THStack* hs_RH_Bias = new THStack("hss_RH_Bias",";RH Bias;Events");
  THStack* hs_Center_Gather = new THStack("hss_Center_Gather",";Center Gather;Events");
  THStack* hs_Center_Spread = new THStack("hss_Center_Spread",";Center Spread;Events");

  for(size_t i=0;i<alphas_LH_Bias.size();i++){
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

  for(size_t j=0;j<alphas_LH_Bias.size();j++){

      double alpha_LH_Bias = alphas_LH_Bias.at(j);
      double beta_LH_Bias = betas_LH_Bias.at(j);

      double alpha_RH_Bias = alphas_RH_Bias.at(j);
      double beta_RH_Bias = betas_RH_Bias.at(j);

      double alpha_Center_Gather = alphas_Center_Gather.at(j);
      double beta_Center_Gather = betas_Center_Gather.at(j);

      double alpha_Center_Spread = alphas_Center_Spread.at(j);
      double beta_Center_Spread = betas_Center_Spread.at(j);

      leg_LH_Bias->AddEntry(h_LH_Bias.at(j),Form("#alpha=%.2f,#beta=%.2f",alpha_LH_Bias,beta_LH_Bias),"l");
      leg_RH_Bias->AddEntry(h_RH_Bias.at(j),Form("#alpha=%.2f,#beta=%.2f",alpha_RH_Bias,beta_RH_Bias),"l");
      leg_Center_Gather->AddEntry(h_Center_Gather.at(j),Form("#alpha=%.2f,#beta=%.2f",alpha_Center_Gather,beta_Center_Gather),"l");
      leg_Center_Spread->AddEntry(h_Center_Spread.at(j),Form("#alpha=%.2f,#beta=%.2f",alpha_Center_Spread,beta_Center_Spread),"l");

        for(int i=1;i<=bins;i++){
        double x = h_LH_Bias.at(0)->GetBinCenter(i);
          h_LH_Bias.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,alpha_LH_Bias,beta_LH_Bias));
          h_RH_Bias.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,alpha_RH_Bias,beta_RH_Bias));
          h_Center_Gather.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,alpha_Center_Gather,beta_Center_Gather));
          h_Center_Spread.at(j)->SetBinContent(i,ROOT::Math::beta_pdf(x,alpha_Center_Spread,beta_Center_Spread));
        }

  }


  TCanvas* c1 = new TCanvas("c1","c1",800,600);

  for(size_t i=0;i<alphas_LH_Bias.size();i++){
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

}

#endif
