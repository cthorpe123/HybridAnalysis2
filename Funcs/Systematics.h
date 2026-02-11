#ifndef _Systematics_h_
#define _Systematics_h_

namespace syst {

// For speed while developing
const int nuniv_Genie = 100;
const int nuniv_Flux = 100;
const int nuniv_Reint = 100;

//const int nuniv_Genie = 500;
//const int nuniv_Flux = 1000;
//const int nuniv_Reint = 1000;

// Multisims
enum e_syst {kGenie,kFlux,kReint,kSystMAX};
const std::vector<std::string> sys_str = {"Genie","Flux","Reint"};
const std::vector<int> sys_nuniv = {nuniv_Genie,nuniv_Flux,nuniv_Reint};
const std::vector<int> sys_color = {kBlue-7,kMagenta-7,kRed-7};

// Unisims
enum e_unisims {kAxFFCCQEshape,kDecayAngMEC,kNormCCCOH,kNormNCCOH,kThetaDelta2NRad,kTheta_Delta2Npi,kVecFFCCQEshape,kXSecShape_CCMEC,kUnisimMAX};
const std::vector<std::string> unisims_str = {"AxFFCCQEshape","DecayAngMEC","NormCCCOH","NormNCCOH","ThetaDelta2NRad","Theta_Delta2Npi","VecFFCCQEshape","XSecShape_CCMEC"};
const int unisim_color = 3; 

// Detvars
enum e_detvars {kLYAtt,kLYDown,kLYRayleigh,kSCE,kRecomb2/*,kWMX*/,kDetvarMAX};
const std::vector<std::string> detvar_str = {"LYAtt","LYDown","LYRayleigh","SCE","Recomb2"/*,"WMX"*/};
const int detvar_color = kGreen+2; 

// Stats
enum e_stats {kMCStat,kDataStat,kStatMAX};
const std::vector<std::string> stat_str = {"MCStat","DataStat"};
const std::vector<int> stat_color = {kBlue+2,kRed+2};

// Special categories, only applied in certain situations
enum e_special {kEstDataStat,kSpecialMAX};
const std::vector<std::string> special_str = {"EstDataStat"}; 
const std::vector<int> special_color = {kRed+2};

///////////////////////////////////////////////////////////////////////
// Mean value of a given bin in stack of multisim vars

double Mean(const std::vector<TH1D*>& h_Vars, int bin){
  double m=0;
  for(TH1D* h : h_Vars) m += h->GetBinContent(bin);
  return m /= h_Vars.size();
}

///////////////////////////////////////////////////////////////////////
// Generate a 2D hist with binning and axis titles matching a 1D hist 
// along both axes 

TH2D* Make2DHist(std::string name,TH1D* h){
  std::string axis_title = h->GetXaxis()->GetTitle();   
  std::vector<double> bins;
  for(int i=1;i<h->GetNbinsX()+2;i++) bins.push_back(h->GetBinLowEdge(i));
  return new TH2D(name.c_str(),(";"+axis_title+";"+axis_title+";").c_str(),bins.size()-1,&bins[0],bins.size()-1,&bins[0]);
}

///////////////////////////////////////////////////////////////////////
// Multisim covariance calculator 

void CalcCovMultisim(std::string sys,std::vector<TH1D*> h_Vars,TH2D*& h_Cov,TH2D*& h_FCov){
 
  std::string axis_title = h_Vars.at(0)->GetXaxis()->GetTitle();   

  std::vector<double> bins;
  for(int i=1;i<h_Vars.at(0)->GetNbinsX()+2;i++) bins.push_back(h_Vars.at(0)->GetBinLowEdge(i));
  int n_bins = bins.size()-1;
  double* bins_a = &bins[0];

  h_Cov = new TH2D(("Cov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);
  h_FCov = new TH2D(("FCov_" + sys).c_str(),(";"+axis_title+";"+axis_title+";").c_str(),n_bins,bins_a,n_bins,bins_a);

  for(int i_bx=0;i_bx<h_Vars.at(0)->GetNbinsX()+2;i_bx++){
    for(int i_by=i_bx;i_by<h_Vars.at(0)->GetNbinsX()+2;i_by++){

      double x = Mean(h_Vars,i_bx);
      double y = Mean(h_Vars,i_by);

      double cov = 0;
      for(int i_u=0;i_u<h_Vars.size();i_u++)
        cov += (h_Vars.at(i_u)->GetBinContent(i_bx) - x)*(h_Vars.at(i_u)->GetBinContent(i_by) - y);
      cov /= h_Vars.size(); 
        
      h_Cov->SetBinContent(i_bx,i_by,cov); 
      h_FCov->SetBinContent(i_bx,i_by,cov/x/y); 
        
      h_Cov->SetBinContent(i_by,i_bx,cov); 
      h_FCov->SetBinContent(i_by,i_bx,cov/x/y); 

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

  for(int i_bx=0;i_bx<h_CV->GetNbinsX()+2;i_bx++){
    for(int i_by=0;i_by<h_CV->GetNbinsX()+2;i_by++){
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

///////////////////////////////////////////////////////////////////////
// Convert a 2D histogram into a tmatrix 

TMatrixDSym MakeCovMat(TH2D* h){

  TMatrixDSym m(h->GetNbinsX());
  for(int i=1;i<h->GetNbinsX()+1;i++)
    for(int j=1;j<h->GetNbinsX()+1;j++)
      m[i-1][j-1] = h->GetBinContent(i,j);

  return m;
}

///////////////////////////////////////////////////////////////////////
// Convert a 2D histogram into a tmatrix 

std::vector<TH1D*> PadUniverses(TH1D* h,TH2D* h_cov,int nuniv){

  //std::cout << "Doing PadUniverses" << std::endl;

  TMatrixDSym m_Cov = MakeCovMat(h_cov); 
  int dim = m_Cov.GetNrows();

  TDecompChol* decomp = new TDecompChol(m_Cov);
  bool valid = decomp->Decompose();
  //std::cout << "valid = " << valid << std::endl;  
 
  if(!valid) return std::vector<TH1D*>();

  TMatrixD m_decomp = decomp->GetU(); 
  TMatrixD m_decomp_t = m_decomp;

  TRandom2* r = new TRandom2();
  std::vector<TH1D*> univ;

  for(int i_u=0;i_u<nuniv;i_u++){
    TMatrixD v(dim,1);
    for(int i=0;i<dim;i++) v[i][0] = r->Gaus(0.0,1.0); 
    TMatrixD weights = m_decomp*v;
    univ.push_back((TH1D*)h->Clone((string(h->GetName())+"_"+std::to_string(i_u)).c_str()));
    for(int i=1;i<h->GetNbinsX()+1;i++) univ.back()->SetBinContent(i,(1+weights[i-1][0])*h->GetBinContent(i));
  }

  delete r;
  delete decomp;

  return univ;

}

///////////////////////////////////////////////////////////////////////
// Select the correct unisim weight 

double ChooseUnisimWeight(int u,std::map<std::string,std::vector<double>>* w){

  if(u == kXSecShape_CCMEC) return w->at("XSecShape_CCMEC_UBGenie").at(1); 
  else if(u < kUnisimMAX) return w->at(unisims_str.at(u)+"_UBGenie").at(0);
 
  return 0;
}

}

#endif

