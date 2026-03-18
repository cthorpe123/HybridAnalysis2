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
#include "DetvarHistograms.h"
#include "PlotFuncs.h"

using namespace syst;

void CheckDetvars(){

  // Label and set the branches defining the selection and systematics
  std::string label = "MuonMom";
  bool draw_truth = true;
  bool draw_hist = true; // Grab the CV from the non-detvar file
  bool draw_o=false,draw_u=false;

  std::string plot_dir = "Analysis/"+label+"/Plots/CheckDetvars/";
  gSystem->Exec(("mkdir -p "+plot_dir).c_str());

  TFile* f_in = TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str());
  TFile* f_in_hist = draw_hist ? TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str()) : nullptr;

  std::vector<TH1D*> h_v;
  std::vector<int> fill_colors;
  std::vector<std::string> legs;

  h_v.push_back((TH1D*)f_in->Get("Reco/CV/h_Tot"));
  fill_colors.push_back(1);
  legs.push_back("CV");
   
  std::cout << std::endl;
  for(int i=1;i<h_v.back()->GetNbinsX()+1;i++) std::cout << h_v.back()->GetBinLowEdge(i) << " " << h_v.back()->GetBinContent(i) << std::endl;

  //for(std::string var : detvar_str){
  for(size_t i_s=0;i_s<detvar_str.size();i_s++){
    std::string var = detvar_str.at(i_s);
    h_v.push_back((TH1D*)f_in->Get(("Reco/Vars/"+var+"/h_Tot").c_str()));
    fill_colors.push_back(i_s+2);
    legs.push_back(var);
  }

  if(draw_hist){
    h_v.push_back((TH1D*)f_in_hist->Get("Reco/CV/h_Tot")->Clone("h_CV_Reco_Tot_Hist"));
    h_v.back()->SetLineStyle(2);
    std::cout << std::endl;
    for(int i=1;i<h_v.back()->GetNbinsX()+1;i++) std::cout << h_v.back()->GetBinLowEdge(i) << " " << h_v.back()->GetBinContent(i) << std::endl;
    fill_colors.push_back(1);
    legs.push_back("CV HIST");
  }

  pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,plot_dir+"Reco.png"); 

  if(draw_truth){
    h_v.clear();
    h_v.push_back((TH1D*)f_in->Get("Truth/CV/h_Signal"));
    for(size_t i_s=0;i_s<detvar_str.size();i_s++){
      std::string var = detvar_str.at(i_s);
      h_v.push_back((TH1D*)f_in->Get(("Truth/Vars/"+var+"/h_Signal").c_str()));
    } 

  if(draw_hist){
    h_v.push_back((TH1D*)f_in_hist->Get("Truth/CV/h_Signal")->Clone("h_CV_Truth_Signal_Hist"));
    h_v.back()->SetLineStyle(2);
  }

    pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,plot_dir+"Truth.png"); 

    // Draw the efficiency afo truth variable
    h_v.clear();

    TH2D* h_CV_Response_Signal = (TH2D*)f_in->Get("Response/CV/h_Signal");
    int bins = h_CV_Response_Signal->GetNbinsY()+1;
    h_v.push_back(h_CV_Response_Signal->ProjectionX("h_CV_Eff",0,bins)); 

    for(size_t i_s=0;i_s<detvar_str.size();i_s++){
      std::string var = detvar_str.at(i_s);
      TH2D* h_Var_Response_Signal = (TH2D*)f_in->Get(("Response/Vars/"+var+"/h_Signal").c_str());
      h_v.push_back(h_Var_Response_Signal->ProjectionX(("h_"+var+"_Eff").c_str(),0,bins));
    } 

    if(draw_hist){
      TH2D* h_CV_Hist_Response_Signal = (TH2D*)f_in_hist->Get("Response/CV/h_Signal");
      h_v.push_back(h_CV_Hist_Response_Signal->ProjectionX("h_CV_Hist_Eff",0,bins));
      h_v.back()->SetLineStyle(2);
    }

    pfs::DrawUnstacked(h_v,fill_colors,legs,draw_o,draw_u,plot_dir+"Efficiency.png"); 

  }

}
