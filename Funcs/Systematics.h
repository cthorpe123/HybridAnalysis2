#ifndef _Systematics_h_
#define _Systematics_h_

namespace syst {

const int nuniv_Genie = 100;
const int nuniv_Flux = 100;
const int nuniv_Reint = 100;

enum e_syst {kGenie,kFlux,kReint,kSystMAX};
const std::vector<std::string> sys_str = {"Genie","Flux","Reint"};
const std::vector<int> sys_nuniv = {nuniv_Genie,nuniv_Flux,nuniv_Reint};

enum e_detvars {kLYAtt,kLYDown,kLYRayleigh,kSCE,kRecomb2/*,kWMX*/,kDetvarMAX};
const std::vector<std::string> detvar_str = {"LYAtt","LYDown","LYRayleigh","SCE","Recomb2"/*,"WMX"*/};

///////////////////////////////////////////////////////////////////////
// Mean value of a given bin in stack of multisim vars

double Mean(const std::vector<TH1D*>& h_Vars, int bin){
  double m=0;
  for(TH1D* h : h_Vars) m += h->GetBinContent(bin);
  return m /= h_Vars.size();
}

///////////////////////////////////////////////////////////////////////
// Multisim covariance calculator 

void CalcCovMultisim(std::string sys,const TH1D* h_CV,std::vector<TH1D*> h_Vars,TH2D*& h_Cov,TH2D*& h_FCov){

  std::cout << "Calculating covariance for " << sys << std::endl;
 
  std::string axis_title = h_CV->GetXaxis()->GetTitle();   

  std::vector<double> bins;
  for(int i=1;i<h_CV->GetNbinsX()+2;i++) bins.push_back(h_CV->GetBinLowEdge(i));
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  h_Cov = new TH2D(("Cov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);
  h_FCov = new TH2D(("FCov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);


  for(int i_bx=1;i_bx<h_CV->GetNbinsX()+1;i_bx++){
    for(int i_by=1;i_by<h_CV->GetNbinsX()+1;i_by++){

      //double x = h_CV->GetBinContent(i_bx);
      //double y = h_CV->GetBinContent(i_by);

      double x = Mean(h_Vars,i_bx);
      double y = Mean(h_Vars,i_by);

      double cov = 0;
      for(int i_u=0;i_u<h_Vars.size();i_u++) cov += (h_Vars.at(i_u)->GetBinContent(i_bx) - x)*(h_Vars.at(i_u)->GetBinContent(i_by) - y);
      cov /= h_Vars.size(); 

      //std::cout << i_bx << "  " << i_by << "  " << x << "  " <<  y << "  " << cov << std::endl;

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
  for(int i=1;i<h_CV->GetNbinsX()+2;i++) bins.push_back(h_CV->GetBinLowEdge(i));
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

///////////////////////////////////////////////////////////////////////
// Calculate correlation matrix

TH2D* CalcCorrelationMatrix(std::string sys,const TH2D* h_Cov){

  TH2D* h_Corr = (TH2D*)h_Cov->Clone(("Corr_"+sys).c_str());

  for(int i_bx=1;i_bx<h_Cov->GetNbinsX()+1;i_bx++)
    for(int i_by=1;i_by<h_Cov->GetNbinsX()+1;i_by++)
      h_Corr->SetBinContent(i_bx,i_by,h_Cov->GetBinContent(i_bx,i_by)/sqrt(h_Cov->GetBinContent(i_bx,i_bx))/sqrt(h_Cov->GetBinContent(i_by,i_by)));

  return h_Corr;
} 

}

#endif

