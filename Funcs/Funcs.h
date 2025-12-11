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
  {111,{0.0,5.0}},
  {22,{0.0,5.0}}
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

///////////////////////////////////////////////////////////////////////
// Generate variably binned histogram  

void MakeBinningTemplate(std::string label,TH1D* h_data,double target_fe=0.05){

  std::cout << "Generating binning template for " << label << std::endl;

  gSystem->Exec(("mkdir -p Analysis/"+label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+label+"/rootfiles/BinningTemplate.root").c_str(),"RECREATE");

  std::vector<double> bin_edges;
  double events = 0;
  const double _EPSILON_ = 1e-10;
  for(int i=1;i<h_data->GetNbinsX()+1;i++){

    if(h_data->GetBinContent(i) < _EPSILON_) continue;

    if(!bin_edges.size()) bin_edges.push_back(h_data->GetBinLowEdge(i));

    events += h_data->GetBinContent(i);
    if(1/sqrt(events) < target_fe){
      bin_edges.push_back(h_data->GetBinLowEdge(i+1));
      events = 0;
    } 

  }

  std::cout << "Number of bins generated " << bin_edges.size()-1 << std::endl;
  std::cout << "Low edge = " << bin_edges.front() << " Hgh edge = " << bin_edges.back() << std::endl;

  TH1D* h_template = new TH1D(("h_template_"+label).c_str(),"",bin_edges.size()-1,&bin_edges[0]); 
  h_template->GetXaxis()->SetTitle(h_data->GetXaxis()->GetTitle());
  h_template->GetYaxis()->SetTitle(h_data->GetYaxis()->GetTitle());
  h_template->Write();
  f_out->Close();

}

///////////////////////////////////////////////////////////////////////
// Make histogram with underflow/overflow 1 bin histograms that can be 
// plotted alongside 

void MakeOU(std::string label,TH1D* h,TH1D*& h_underflow,TH1D*& h_overflow,std::string u_label="Underflow",std::string o_label="Overflow",TH2D* h_cov=nullptr){

  h_underflow = new TH1D(("h_underflow_"+label).c_str(),"",1,0.0,1.0);
  h_overflow = new TH1D(("h_overflow_"+label).c_str(),"",1,0.0,1.0);

  h_underflow->SetBinContent(1,h->GetBinContent(0));
  h_overflow->SetBinContent(1,h->GetBinContent(h->GetNbinsX()+1));

  if(h_cov != nullptr) h_underflow->SetBinError(1,sqrt(h_cov->GetBinContent(0,0)));
  else h_underflow->SetBinError(1,h->GetBinError(0));
  if(h_cov != nullptr) h_overflow->SetBinError(1,sqrt(h_cov->GetBinContent(h->GetNbinsX()+1,h->GetNbinsX()+1)));
  else h_overflow->SetBinError(1,h->GetBinContent(h->GetNbinsX()+1));
  
  h_underflow->GetXaxis()->SetBinLabel(1,u_label.c_str());
  h_overflow->GetXaxis()->SetBinLabel(1,o_label.c_str());

}

///////////////////////////////////////////////////////////////////////
// Get maximum of hist with error bars 

double GetMax(const TH1D* h){
  double max = 0;
  for(int i=1;i<h->GetNbinsX()+1;i++) max = std::max(max,h->GetBinContent(i)+h->GetBinError(i));
  return max;
}

///////////////////////////////////////////////////////////////////////
// Unroll a 2D dist into one large 1D dist

TH1D* Unroll2DDist(TH2D* h,std::vector<std::pair<int,int>> skip={}){

  int nbins_x = h->GetNbinsX();
  int nbins_y = h->GetNbinsY();
  std::string name = string(h->GetName())+"_unroll";
  TH1D* h_unroll = new TH1D(name.c_str(),"",nbins_x*nbins_y-skip.size(),0.5,nbins_x*nbins_y-skip.size()+0.5);

  int ctr = 1;
  for(int i=1;i<nbins_x+1;i++){
    for(int j=1;j<nbins_y+1;j++){
      if(std::find(skip.begin(),skip.end(),std::make_pair(i,j)) != skip.end()) continue;
      h_unroll->SetBinContent(ctr,h->GetBinContent(i,j));
      h_unroll->SetBinError(ctr,h->GetBinError(i,j));
      ctr++;
    }
  }

  return h_unroll;
}

///////////////////////////////////////////////////////////////////////
// Given the true dist and joint truth/reco dist for selected events
// renormalise the 2D hist to give the response 

void NormaliseResponse(TH1D* h_true,TH2D* h_true_reco){

  for(int i=0;i<h_true->GetNbinsX()+2;i++){
    double total_true = h_true->GetBinContent(i);
    if(total_true > 0){
      for(int j=0;j<h_true_reco->GetNbinsY()+2;j++){
        h_true_reco->SetBinContent(i,j,h_true_reco->GetBinContent(i,j)/total_true);
        h_true_reco->SetBinError(i,j,h_true_reco->GetBinError(i,j)/total_true);
      }
    }
  }

}

#endif
