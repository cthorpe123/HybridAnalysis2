#ifndef _Systematics_h_
#define _Systematics_h_

namespace syst {

const int nuniv_Genie = 100;
const int nuniv_Flux = 100;
const int nuniv_Reint = 100;

enum e_syst {kGenie,kFlux,kReint,kSystMAX};
const std::vector<std::string> sys_str = {"Genie","Flux","Reint"};
const std::vector<int> sys_nuniv = {nuniv_Genie,nuniv_Flux,nuniv_Reint};

enum e_detvars {kLYAtt,kLYDown,kLYRayleigh,kWMX,kSCE,kRecomb2,kDetvarMAX};
const std::vector<std::string> detvar_str = {"LYAtt","LYDown","LYRayleigh","WMX","SCE","Recomb2"};

///////////////////////////////////////////////////////////////////////
// Multisim covariance calculator 

void CalcCovMultisim(std::string sys,const TH1D* h_CV,std::vector<TH1D*> h_Vars,TH2D*& h_Cov,TH2D*& h_FCov){

  std::string axis_title = h_CV->GetXaxis()->GetTitle();   

  std::vector<double> bins;
  for(int i=1;i<h_CV->GetNbinsX()+1;i++) bins.push_back(h_CV->GetBinLowEdge(i));
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  h_Cov = new TH2D(("Cov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);
  h_FCov = new TH2D(("FCov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);

  for(int i_bx=1;i_bx<h_CV->GetNbinsX()+1;i_bx++){
    for(int i_by=1;i_by<h_CV->GetNbinsX()+1;i_by++){
      double x = h_CV->GetBinContent(i_bx);
      double y = h_CV->GetBinContent(i_by);

      double cov = 0;
      for(int i_u=0;i_u<h_Vars.size();i_u++) cov += (h_Vars.at(i_u)->GetBinContent(i_bx) - x)*(h_Vars.at(i_u)->GetBinContent(i_by) - y);
      cov /= h_Vars.size(); 

      h_Cov->SetBinContent(i_bx,i_by,cov); 
      h_FCov->SetBinContent(i_bx,i_by,cov/x/y); 

    }
  }

}

///////////////////////////////////////////////////////////////////////
// Unisim covariance calculator 

void CalcCovUnisim(std::string sys,const TH1D* h_CV,TH1D* h_Var,TH2D*& h_Cov,TH2D*& h_FCov){

  std::string axis_title = h_CV->GetXaxis()->GetTitle();   

  std::vector<double> bins;
  for(int i=1;i<h_CV->GetNbinsX()+1;i++) bins.push_back(h_CV->GetBinLowEdge(i));
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  h_Cov = new TH2D(("Cov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);
  h_FCov = new TH2D(("FCov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);

  for(int i_bx=1;i_bx<h_CV->GetNbinsX()+1;i_bx++){
    for(int i_by=1;i_by<h_CV->GetNbinsX()+1;i_by++){
      double x = h_CV->GetBinContent(i_bx);
      double y = h_CV->GetBinContent(i_by);
      double cov = (h_Var->GetBinContent(i_bx) - x)*(h_Var->GetBinContent(i_by) - y);
      h_Cov->SetBinContent(i_bx,i_by,cov); 
      h_FCov->SetBinContent(i_bx,i_by,cov/x/y); 
    }
  }

}

}

#endif

