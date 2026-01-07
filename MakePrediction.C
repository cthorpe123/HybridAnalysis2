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

using namespace syst;

void MakePrediction(){

  bool add_detvars = false;
  bool blinded = true;
  bool draw_underflow = false;
  bool draw_overflow = false;

  std::vector<std::string> label_v = {"MuonMom"};

  for(size_t i_f=0;i_f<label_v.size();i_f++){
    std::string label = label_v.at(i_f);

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;

    // Build the CV prediction
    TH1D* h_CV_Tot = static_cast<TH1D*>(f_in_hist->Get("Reco/CV/h_Tot"));

    std::vector<TH1D*> h_CV;
    std::vector<int> fill_colors;
    std::vector<std::string> legs;
    for(size_t i_c=0;i_c<categories.size();i_c++){
      if(i_c == kData) continue;
      h_CV.push_back(static_cast<TH1D*>(f_in_hist->Get(("Reco/CV/h_"+categories.at(i_c)).c_str())));
      fill_colors.push_back(cat_colors[i_c]); 
      legs.push_back(categories.at(i_c));
    }

    TH2D* h_Cov = static_cast<TH2D*>(f_in_hist->Get("Reco/Cov/Total/Cov_Tot"));
    for(int i=0;i<h_CV_Tot->GetNbinsX()+2;i++) h_CV_Tot->SetBinError(i,sqrt(h_Cov->GetBinContent(i,i)));

    // If not blindded, load the data
    TH1D* h_Data = !blinded ? (TH1D*)f_in_hist->Get("Reco/CV/h_Data") : nullptr; 

    std::string plot_dir = "Analysis/"+label+"/Plots/MakePrediction/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());
    pfs::DrawStacked(h_CV,fill_colors,legs,h_CV_Tot,h_Data,draw_overflow,draw_underflow,plot_dir+"Pred.png");

  }

}
