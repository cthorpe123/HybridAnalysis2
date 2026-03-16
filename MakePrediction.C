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

void MakePredictionMB(){

  bool add_detvars = false;
  bool blinded = true;
  bool draw_underflow = false;
  bool draw_overflow = true;
  bool divide_by_bin_width = true;

  std::vector<std::string> label_v = {"MuonMom","MuonMom2"};

  for(size_t i_f=0;i_f<label_v.size();i_f++){
    std::string label = label_v.at(i_f);

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
