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

void SysTest(){

  bool blinded = true;

  std::string label = "RecoE_Test";

  TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
  TFile* f_in_detvar = TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str());

  std::string plot_dir = "Analysis/"+label+"/Plots/SysTest/";
  gSystem->Exec(("mkdir -p "+plot_dir).c_str());
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");
  l->SetNColumns(3);

  // Print total covariance and total fractional covariance
  TH2D* h_Cov = static_cast<TH2D*>(f_in_hist->Get("Cov"));
  h_Cov->Draw("colz");
  h_Cov->SetStats(0);
  c->Print((plot_dir+"Cov.png").c_str());
  c->Clear(); 

  TH2D* h_FCov = static_cast<TH2D*>(f_in_hist->Get("FCov"));
  h_FCov->Draw("colz");
  h_FCov->SetStats(0);
  c->Print((plot_dir+"FCov.png").c_str());
  c->Clear(); 

  TH2D* h_Corr = CalcCorrelationMatrix("",h_Cov);
  h_Corr->Draw("colz");
  h_Corr->SetStats(0);
  c->Print((plot_dir+"Corr.png").c_str());
  c->Clear(); 

  // Print Genie, Flux and Reint covariance and fractional covariance
  for(int i_s=0;i_s<kSystMAX;i_s++){
    TH2D* h = static_cast<TH2D*>(f_in_hist->Get(("Cov_"+sys_str.at(i_s)).c_str()));
    h->Draw("colz");
    h->SetStats(0);
    c->Print((plot_dir+"Cov_"+sys_str.at(i_s)+".png").c_str());
    c->Clear();
    TH2D* h2 = CalcCorrelationMatrix(sys_str.at(i_s),h);
    h2->Draw("colz");
    c->Print((plot_dir+"Corr_"+sys_str.at(i_s)+".png").c_str());
    c->Clear(); 
    delete h;
    delete h2;
  }

  for(int i_s=0;i_s<kSystMAX;i_s++){
    TH2D* h = static_cast<TH2D*>(f_in_hist->Get(("FCov_"+sys_str.at(i_s)).c_str()));
    h->Draw("colz");
    h->SetStats(0);
    c->Print((plot_dir+"FCov_"+sys_str.at(i_s)+".png").c_str());
    delete h;
  }

  // Print the detector covarianc and fractional covariance  
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    TH2D* h = static_cast<TH2D*>(f_in_detvar->Get(("Cov_"+detvar_str.at(i_s)).c_str()));
    h->Draw("colz");
    h->SetStats(0);
    c->Print((plot_dir+"Cov_"+detvar_str.at(i_s)+".png").c_str());
    TH2D* h2 = CalcCorrelationMatrix(detvar_str.at(i_s),h);
    h2->Draw("colz");
    h2->SetStats(0);
    c->Print((plot_dir+"Corr_"+detvar_str.at(i_s)+".png").c_str());
    delete h;
    delete h2;
  }

  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    TH2D* h = static_cast<TH2D*>(f_in_detvar->Get(("FCov_"+detvar_str.at(i_s)).c_str()));
    h->Draw("colz");
    h->SetStats(0);
    c->Print((plot_dir+"FCov_"+detvar_str.at(i_s)+".png").c_str());
    delete h;
  }

  std::string axis_title = h_Cov->GetXaxis()->GetTitle();

  TH1D* h_FE_Tot = (TH1D*)f_in_hist->Get("h_CV_Tot")->Clone("h_FE_Tot");
  TH2D* h_Tot_Detvar = (TH2D*)f_in_detvar->Get("FCov")->Clone("FCov_Detvar");
  TH2D* h_Tot = (TH2D*)f_in_hist->Get("FCov")->Clone("FCov");

  THStack* hs_FE = new THStack("hs_FE",h_FE_Tot->GetTitle());

  // Calculate total frac error from all uncertainties
  for(int i=1;i<h_FE_Tot->GetNbinsX()+1;i++) h_FE_Tot->SetBinContent(i,sqrt(h_Tot->GetBinContent(i,i) + h_Tot_Detvar->GetBinContent(i,i)));
  h_FE_Tot->SetLineColor(1);
  h_FE_Tot->SetLineWidth(2);
  l->AddEntry(h_FE_Tot,"Total","L");
  hs_FE->Add(h_FE_Tot);

  // Calculate total frac error from all detvars combined
  TH1D* h_FE_Tot_Detvar = (TH1D*)f_in_detvar->Get("h_Detvar_CV_Tot")->Clone("h_Detvar_FE_Tot");
  for(int i=1;i<h_FE_Tot_Detvar->GetNbinsX()+1;i++) h_FE_Tot_Detvar->SetBinContent(i,sqrt(h_Tot_Detvar->GetBinContent(i,i)));
  h_FE_Tot_Detvar->SetLineColor(2);
  h_FE_Tot_Detvar->SetLineWidth(2);
  l->AddEntry(h_FE_Tot_Detvar,"Detector","L");
  hs_FE->Add(h_FE_Tot_Detvar);

  std::vector<TH1D*> h_FE;
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_FE.push_back((TH1D*)f_in_hist->Get("h_CV_Tot")->Clone(("h_FE_"+sys_str.at(i_s)).c_str()));
    TH2D* h = static_cast<TH2D*>(f_in_hist->Get(("FCov_"+sys_str.at(i_s)).c_str()));
    for(int i=1;i<h_FE.back()->GetNbinsX()+1;i++) h_FE.back()->SetBinContent(i,sqrt(h->GetBinContent(i,i)));
    h_FE.back()->SetLineColor(i_s+3);
    h_FE.back()->SetLineWidth(2);
    delete h;
    hs_FE->Add(h_FE.back());
    l->AddEntry(h_FE.back(),sys_str.at(i_s).c_str(),"L");
  }

  hs_FE->Draw("nostack HIST");
  hs_FE->GetXaxis()->SetTitle(axis_title.c_str());
  l->Draw();
  c->Print((plot_dir+"FE.png").c_str());
  c->Clear();
  l->Clear(); 

  THStack* hs_FE_Detvar = new THStack("hs_FE_Detvar",h_FE_Tot_Detvar->GetTitle());

  h_FE_Tot_Detvar->SetLineColor(1);
  hs_FE_Detvar->Add(h_FE_Tot_Detvar);
  l->AddEntry(h_FE_Tot_Detvar,"Total","L"); 

  std::vector<TH1D*> h_FE_Detvar;
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    h_FE_Detvar.push_back((TH1D*)f_in_detvar->Get("h_Detvar_CV_Tot")->Clone(("h_FE_"+detvar_str.at(i_s)).c_str()));
    TH2D* h = static_cast<TH2D*>(f_in_detvar->Get(("FCov_"+detvar_str.at(i_s)).c_str()));
    for(int i=1;i<h_FE_Detvar.back()->GetNbinsX()+1;i++) h_FE_Detvar.back()->SetBinContent(i,sqrt(h->GetBinContent(i,i)));
    delete h; 
    h_FE_Detvar.back()->SetLineColor(i_s+2);
    h_FE_Detvar.back()->SetLineWidth(2);
    l->AddEntry(h_FE_Detvar.back(),detvar_str.at(i_s).c_str(),"L");
    hs_FE_Detvar->Add(h_FE_Detvar.back());    
  }

  hs_FE_Detvar->Draw("nostack HIST");
  hs_FE_Detvar->GetXaxis()->SetTitle(axis_title.c_str());
  hs_FE_Detvar->GetYaxis()->SetTitle("Frac. Uncertainty");
  l->Draw();
  c->Print((plot_dir+"FE_Detvar.png").c_str());
  c->Clear();
  l->Clear(); 

  // Validation - check the mean of each multiverse systematic matches the CV prediction
  for(int i_s=0;i_s<kSystMAX;i_s++){

    THStack* hs = new THStack("hs","");

    TH1D* h_CV = (TH1D*)f_in_hist->Get("h_CV_Tot")->Clone(("h_CV_"+sys_str.at(i_s)).c_str());
    TH1D* h_Var_Mean = (TH1D*)f_in_hist->Get(("h_Vars_Tot_"+sys_str.at(i_s)+"_0").c_str())->Clone("h_Var_Mean"); 

    // Calculate mean of alt universes 
    for(int i_u=1;i_u<sys_nuniv.at(i_s);i_u++)
      h_Var_Mean->Add((TH1D*)f_in_hist->Get(("h_Vars_Tot_"+sys_str.at(i_s)+"_"+std::to_string(i_u)).c_str())); 
    h_Var_Mean->Scale(1.0/sys_nuniv.at(i_s));

    h_CV->SetLineColor(1);
    h_CV->SetLineWidth(2);
    hs->Add(h_CV);
    l->AddEntry(h_CV,"CV","L");

    h_Var_Mean->SetLineColor(2);
    h_Var_Mean->SetLineWidth(2);
    hs->Add(h_Var_Mean);
    l->AddEntry(h_Var_Mean,"Alt Univ. Mean","L");

    hs->Draw("HIST nostack"); 
    hs->GetXaxis()->SetTitle(h_CV->GetXaxis()->GetTitle());
    hs->GetYaxis()->SetTitle(h_CV->GetYaxis()->GetTitle());
    l->Draw();
    c->Print((plot_dir+"SysValidation_"+sys_str.at(i_s)+".png").c_str());

    c->Clear();
    l->Clear();

    delete h_CV;
    delete h_Var_Mean;

  }

  // Detvar validation - draw each variation alongside the CV
  TH1D* h_Detvar_CV = (TH1D*)f_in_detvar->Get("h_Detvar_CV_Tot");
  h_Detvar_CV->SetLineColor(1);
  h_Detvar_CV->SetLineWidth(2);

  for(int i_s=0;i_s<kDetvarMAX;i_s++){

    THStack* hs = new THStack("hs","");

    TH1D* h = static_cast<TH1D*>(f_in_detvar->Get(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)).c_str()));
    h->SetLineColor(2);  
    h->SetLineWidth(2);

    hs->Add(h_Detvar_CV); 
    hs->Add(h);

    l->AddEntry(h,detvar_str.at(i_s).c_str(),"L");
    l->AddEntry(h_Detvar_CV,"CV","L");

    hs->Draw("HIST nostack");
    hs->GetXaxis()->SetTitle(axis_title.c_str());
    l->Draw();
    c->Print((plot_dir+"SysValidation_"+detvar_str.at(i_s)+".png").c_str());
    c->Clear();
    l->Clear();

    delete hs;  

  }

}
