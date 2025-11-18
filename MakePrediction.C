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

using namespace syst;

void MakePrediction(){

  bool blinded = true;

  std::string label = "RecoE_Test";

  TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
  TFile* f_in_detvar = TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str());

  // Build the CV prediction
  TH1D* h_CV_Tot = static_cast<TH1D*>(f_in_hist->Get("h_CV_Tot"));

  std::vector<TH1D*> h_CV;
  for(size_t i_c=0;i_c<categories.size();i_c++){
    if(i_c == kData) continue;
    h_CV.push_back(static_cast<TH1D*>(f_in_hist->Get(("h_CV_"+categories.at(i_c)).c_str())));
  }

  // Build total covariance 
  TH2D* h_Cov_Sum = static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(0)).c_str())->Clone("h_Cov_Sum"));
  h_Cov_Sum->Reset();
  for(int i_s=0;i_s<kSystMAX;i_s++)
      h_Cov_Sum->Add(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(i_s)).c_str())));

  TH2D* h_Cov_Sum_Detvar = static_cast<TH2D*>(f_in_detvar->Get(("Cov_"+detvar_str.at(0)).c_str())->Clone("h_Cov_Sum_Detvar"));
  h_Cov_Sum_Detvar->Reset();
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    h_Cov_Sum_Detvar->Add(static_cast<TH2D*>(f_in_detvar->Get(("Cov_"+detvar_str.at(i_s)).c_str())));

  h_Cov_Sum->Add(h_Cov_Sum_Detvar);

  for(int i=1;i<h_CV_Tot->GetNbinsX()+1;i++){
    h_CV_Tot->SetBinError(i,sqrt(h_Cov_Sum->GetBinContent(i,i)));
  }

  std::string plot_dir = "Analysis/"+label+"/Plots/MakePrediction/";
  gSystem->Exec(("mkdir -p "+plot_dir).c_str());

  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");
  l->SetNColumns(3);

  THStack* hs = new THStack("hs",h_CV_Tot->GetTitle());

  for(size_t i_c=0;i_c<categories.size();i_c++){
    if(i_c == kData) continue;
    h_CV.at(i_c)->SetFillColor(i_c+2);
    l->AddEntry(h_CV.at(i_c),categories.at(i_c).c_str(),"F");
    hs->Add(h_CV.at(i_c));
  }

  h_CV_Tot->SetFillStyle(3253);
  h_CV_Tot->SetFillColor(1);

  hs->Draw("HIST");
  h_CV_Tot->Draw("same e2");
  l->Draw();
  c->Print((plot_dir+"Prediction.png").c_str());
  c->Clear();

/*
  std::vector<std::vector<TH2D*>> h_Cov;
  for(size_t i_c=0;i_c<categories.size();i_c++){
    if(i_c == kData) continue;
    h_Cov.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Cov.back().push_back(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+categories.at(i_c)+"_"+sys_str.at(i_s)).c_str())));
      h_Cov_Sum->Add(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+categories.at(i_c)+"_"+sys_str.at(i_s)).c_str())));
    }
  } 
*/

/*
  std::vector<std::vector<TH2D*>> h_Cov_Detvars;
  for(size_t i_c=0;i_c<categories.size();i_c++){
    if(i_c == kData) continue;
    h_Cov.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kSystMAX;i_s++){
      h_Cov.back().push_back(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+categories.at(i_c)+"_"+sys_str.at(i_s)).c_str())));
      h_Cov_Sum->Add(static_cast<TH2D*>(f_in_hist->Get(("Cov_"+categories.at(i_c)+"_"+sys_str.at(i_s)).c_str())));
    }
  } 
  


*/



}
