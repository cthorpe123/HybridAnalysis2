#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "MultiChannelHistograms.h"
#include "WeightFuncs.h"
#include "PlotFuncs.h"

void CheckCorr(){

  std::string var = "NProt";

  std::string plot_dir = "Analysis/"+var+"/Plots/CheckCorr/";
  gSystem->Exec(("mkdir -p " + plot_dir).c_str());

  TFile* f_hist = TFile::Open(("Analysis/"+var+"/rootfiles/Histograms.root").c_str());
  const double POT = ((TH1D*)f_hist->Get("Meta/POT"))->GetBinContent(1);

  TH2D* h_eff_cv = (TH2D*)f_hist->Get("Response/CV/h_Signal");
  TH1D* h_bg_cv = (TH1D*)f_hist->Get("Reco/CV/h_AllBG");
  TH1D* h_eff_proj_cv = h_eff_cv->ProjectionX((std::string(h_eff_cv->GetName())+"_px_cv").c_str());

  // Block 1 = efficiency (ProjectionX of the truth-reco signal matrix), block 2 = background
  const int n_eff = h_eff_cv->GetNbinsX();
  const int n_bg = h_bg_cv->GetNbinsX();
  const int n_tot = n_eff+n_bg;

  TH1D* h_eff_bg_cv = new TH1D("h_eff_bg_cv","",n_tot,-0.5,n_tot-0.5);
  for(int i_b=0;i_b<n_eff;i_b++) h_eff_bg_cv->SetBinContent(i_b+1,h_eff_proj_cv->GetBinContent(i_b+1));
  for(int i_b=0;i_b<n_bg;i_b++) h_eff_bg_cv->SetBinContent(n_eff+i_b+1,h_bg_cv->GetBinContent(i_b+1));

  TH2D* h_fcov = new TH2D("h_eff_bg","",n_tot,-0.5,n_tot-0.5,n_tot,-0.5,n_tot-0.5);

  for(int i_s=0;i_s<kSystMAX;i_s++){
    gSystem->Exec(("mkdir -p " + plot_dir + "/" + sys_str.at(i_s)).c_str());
    std::vector<TH1D*> h_eff_bg_v;
    std::vector<TH1D*> h_eff_bg_ratio_v;
    std::vector<int> cols;
    for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
      TH2D* h_eff = (TH2D*)f_hist->Get(("Response/Vars/"+sys_str.at(i_s)+"/h_Signal_"+std::to_string(i_u)).c_str());
      TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+sys_str.at(i_s)+"/h_AllBG_"+std::to_string(i_u)).c_str());
      TH1D* h_eff_proj = h_eff->ProjectionX((std::string(h_eff->GetName())+"_px").c_str());
      h_eff_bg_v.push_back(new TH1D(("h_eff_bg_"+std::to_string(i_u)).c_str(),"",n_tot,-0.5,n_tot-0.5));
      for(int i_b=0;i_b<n_eff;i_b++) h_eff_bg_v.back()->SetBinContent(i_b+1,h_eff_proj->GetBinContent(i_b+1));
      for(int i_b=0;i_b<n_bg;i_b++) h_eff_bg_v.back()->SetBinContent(n_eff+i_b+1,h_bg->GetBinContent(i_b+1));
      if(i_u < 5){
        TH1D* h = (TH1D*)h_eff_bg_v.back()->Clone(Form("h%i",i_u));
        h->Add(h_eff_bg_cv,-1);
        h_eff_bg_ratio_v.push_back(h);
        cols.push_back(i_u+1);
      }
    }
  
    TH2D *c,*fc;
    CalcCovMultisim(sys_str.at(i_s),h_eff_bg_v,c,fc);
    h_fcov->Add(fc);

    for(TH1D* hh : h_eff_bg_ratio_v){
      for(int i_b=1;i_b<hh->GetNbinsX()+1;i_b++){
        double err = hh->GetBinContent(i_b)/sqrt(c->GetBinContent(i_b,i_b));
        hh->SetBinContent(i_b,err);
      }
    }

    pfs::DrawUnstacked2(h_eff_bg_ratio_v,cols,std::vector<std::string>(h_eff_bg_ratio_v.size(),""),plot_dir+"/"+sys_str.at(i_s)+"/EffBG.png",false);

    TH2D* corr = CalcCorrelationMatrix("EffBG",c);
    pfs::Draw2DHist(corr,(plot_dir+"Corr_EffBG_"+sys_str.at(i_s)+".png").c_str());
    pfs::Draw2DHist(fc,(plot_dir+"Cov_EffBG_"+sys_str.at(i_s)+".png").c_str());

    for(TH1D* hh : h_eff_bg_v) delete hh;
    for(TH1D* hh : h_eff_bg_ratio_v) delete hh;
    h_eff_bg_v.clear();
    h_eff_bg_ratio_v.clear();
    delete corr;

  }

  // Unisims - group together into one systematic


  for(int i_s=0;i_s<kUnisimMAX;i_s++){

    TH2D* h_eff = (TH2D*)f_hist->Get(("Response/Vars/"+unisims_str.at(i_s)+"/h_Signal").c_str());
    TH1D* h_bg = (TH1D*)f_hist->Get(("Reco/Vars/"+unisims_str.at(i_s)+"/h_AllBG").c_str());
    TH1D* h_eff_proj = h_eff->ProjectionX((std::string(h_eff->GetName())+"_px").c_str());
    TH1D* h_eff_bg = new TH1D(("h_eff_bg_"+unisims_str.at(i_s)).c_str(),"",n_tot,-0.5,n_tot-0.5);
    for(int i_b=0;i_b<n_eff;i_b++) h_eff_bg->SetBinContent(i_b+1,h_eff_proj->GetBinContent(i_b+1));
    for(int i_b=0;i_b<n_bg;i_b++) h_eff_bg->SetBinContent(n_eff+i_b+1,h_bg->GetBinContent(i_b+1));

    TH2D *c,*fc;
    CalcCovUnisim(unisims_str.at(i_s),h_eff_bg_cv,h_eff_bg,c,fc);
    h_fcov->Add(fc);

    TH2D* corr = CalcCorrelationMatrix("EffBG",c);
    pfs::Draw2DHist(corr,(plot_dir+"Corr_EffBG_"+unisims_str.at(i_s)+".png").c_str());
    pfs::Draw2DHist(fc,(plot_dir+ "Cov_EffBG_"+unisims_str.at(i_s)+".png").c_str());

    delete c;
    delete fc;
    
  }

  pfs::Draw2DHist(h_fcov,(plot_dir+"Cov_EffBG_All.png").c_str());
  TH2D* corr = CalcCorrelationMatrix("EffBG",h_fcov);
  pfs::Draw2DHist(corr,(plot_dir+"Corr_EffBG_All.png").c_str());
  delete corr;

}