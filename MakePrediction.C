#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"
#include "Systematics.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

void MakePrediction(){

  bool add_detvars = false;
  bool blinded = true;
  bool draw_underflow = false;
  bool draw_overflow = false;
  bool divide_by_bin_width = false;

  std::vector<std::string> vars = {"MuonMom"};

  for(size_t i_f=0;i_f<vars.size();i_f++){
    std::string label = vars.at(i_f);

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;
 
    hist::MultiChannelHistogramManager mchm(label);
    mchm.LoadTemplates();
    
    // Build the CV prediction
    TH1D* h_CV_Tot = (TH1D*)(f_in_hist->Get("Reco/CV/h_Tot"));;
    mchm.Restore(h_CV_Tot);

    std::vector<TH1D*> h_CV;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;
    for(size_t i_c=0;i_c<categories.size();i_c++){
      if(i_c == kData) continue;
      h_CV.push_back((TH1D*)f_in_hist->Get(("Reco/CV/h_"+categories.at(i_c)).c_str()));
      mchm.Restore(h_CV.back());
      if(divide_by_bin_width) DivideByBinWidth(h_CV.back());
      fill_colors.push_back(cat_colors[i_c]); 
      legs.push_back(categories.at(i_c));
    }

    TH2D* h_Cov = (TH2D*)(f_in_hist->Get("Reco/Cov/Total/Cov_Tot"));

    // Detvar prediction isn't guaranteed to match CV's exposure - use FCov and scale
    if(add_detvars){
      TH2D* h_FCov_Detvar = (TH2D*)f_in_detvar->Get("Reco/Cov/Total/FCov_Tot");
      for(int i=0;i<h_CV_Tot->GetNbinsX()+2;i++)
        for(int j=0;j<h_CV_Tot->GetNbinsX()+2;j++)
          h_FCov_Detvar->SetBinContent(i,j,h_FCov_Detvar->GetBinContent(i,j)*h_CV_Tot->GetBinContent(i)*h_CV_Tot->GetBinContent(j));
      h_Cov->Add(h_FCov_Detvar); 
    }

    mchm.Restore(h_Cov);
    for(int i=0;i<h_CV_Tot->GetNbinsX()+2;i++) h_CV_Tot->SetBinError(i,sqrt(h_Cov->GetBinContent(i,i)));
    if(divide_by_bin_width) DivideByBinWidth(h_CV_Tot);
        
    // If not blindded, load the data
    TH1D* h_Data = nullptr;
    if(!blinded){
      h_Data = (TH1D*)f_in_hist->Get("Reco/CV/h_Data");
      mchm.Restore(h_Data);
      if(divide_by_bin_width) DivideByBinWidth(h_Data);
    }

    std::string plot_dir = "Analysis/"+label+"/Plots/MakePrediction/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());
    pfs::DrawStacked(h_CV,fill_colors,legs,h_CV_Tot,h_Data,draw_overflow,draw_underflow,plot_dir+"Pred.png");

  }

}
