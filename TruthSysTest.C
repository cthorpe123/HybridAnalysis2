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

void TruthSysTest(){

  bool draw_underflow = false;
  bool draw_overflow = true;

  bool add_detvars = false;
  bool blinded = true;

  std::vector<std::string> label_v = {"MuonMom"};

  for(size_t i_f=0;i_f<label_v.size();i_f++){

    std::string label = label_v.at(i_f);

    TFile* f_in_hist = TFile::Open(("Analysis/"+label+"/rootfiles/Histograms.root").c_str());
    TFile* f_in_detvar = add_detvars ? TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str()) : nullptr;

    std::string plot_dir = "Analysis/"+label+"/Plots/SysTest/Truth/";
    gSystem->Exec(("mkdir -p "+plot_dir).c_str());

    std::string sub = "Cov/Truth/Total/"; 

    // Print total covariance and total fractional covariance
    TH2D* h_Cov = static_cast<TH2D*>(f_in_hist->Get((sub+"Cov_Signal").c_str()));
    pfs::Draw2DHist(h_Cov,plot_dir+"Cov.png");

    TH2D* h_FCov = static_cast<TH2D*>(f_in_hist->Get((sub+"FCov_Signal").c_str()));
    pfs::Draw2DHist(h_FCov,plot_dir+"FCov.png");

    TH2D* h_Corr = CalcCorrelationMatrix("",h_Cov);
    pfs::Draw2DHist(h_Corr,plot_dir+"Corr.png");

    // Print Genie, Flux and Reint covariance and fractional covariance
    for(int i_s=0;i_s<kSystMAX;i_s++){
      plot_dir = "Analysis/"+label+"/Plots/SysTest/Truth/";
      sub = "Cov/Truth/"+sys_str.at(i_s)+"/"; 
      TH2D* h = static_cast<TH2D*>(f_in_hist->Get((sub+"Cov_Signal").c_str()));
      pfs::Draw2DHist(h,plot_dir+"Cov_"+sys_str.at(i_s)+".png");
      TH2D* hf = static_cast<TH2D*>(f_in_hist->Get((sub+"FCov_Signal").c_str()));
      pfs::Draw2DHist(hf,plot_dir+"FCov_"+sys_str.at(i_s)+".png");
      TH2D* h2 = CalcCorrelationMatrix(sys_str.at(i_s),h);
      pfs::Draw2DHist(h2,plot_dir+"Corr_"+sys_str.at(i_s)+".png");
      delete h;
      delete hf;
      delete h2;
    }

    plot_dir = "Analysis/"+label+"/Plots/SysTest/Truth/";
    TH2D* h_Cov_MCStat = static_cast<TH2D*>(f_in_hist->Get("Cov/Truth/MCStat/Cov_Signal"));
    pfs::Draw2DHist(h_Cov_MCStat,plot_dir+"Cov_MCStat.png");

    TH2D* h_FCov_MCStat = static_cast<TH2D*>(f_in_hist->Get("Cov/Truth/MCStat/FCov_Signal"));
    pfs::Draw2DHist(h_FCov_MCStat,plot_dir+"FCov_MCStat.png");

    TH2D* h_Cov_EstDataStat = static_cast<TH2D*>(f_in_hist->Get("Cov/Truth/EstDataStat/Cov_Signal"));
    pfs::Draw2DHist(h_Cov_EstDataStat,plot_dir+"Cov_EstDataStat.png");

    TH2D* h_FCov_EstDataStat = static_cast<TH2D*>(f_in_hist->Get("Cov/Truth/EstDataStat/FCov_Signal"));
    pfs::Draw2DHist(h_FCov_EstDataStat,plot_dir+"FCov_EstDataStat.png");

    plot_dir = "Analysis/"+label+"/Plots/SysTest/Truth/";
    std::vector<TH1D*> h_FE;
    std::vector<int> colors;
    std::vector<std::string> legs;
    for(int i_s=0;i_s<kSystMAX;i_s++){
      sub = "Cov/Truth/"+sys_str.at(i_s)+"/"; 
      h_FE.push_back((TH1D*)f_in_hist->Get("CV/Truth/h_Signal")->Clone(("h_FE_"+sys_str.at(i_s)).c_str()));
      colors.push_back(sys_color.at(i_s));
      legs.push_back(sys_str.at(i_s));
      TH2D* h = static_cast<TH2D*>(f_in_hist->Get((sub+"FCov_Signal").c_str()));
      for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h->GetBinContent(i,i)));
      delete h;
    }

    h_FE.push_back((TH1D*)f_in_hist->Get("CV/Truth/h_Signal")->Clone("h_FE_MCStat"));
    for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h_FCov_MCStat->GetBinContent(i,i)));
    colors.push_back(stat_color[kMCStat]);
    legs.push_back("MCStat");

    h_FE.push_back((TH1D*)f_in_hist->Get("CV/Truth/h_Signal")->Clone("h_FE_EstDataStat"));
    for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h_FCov_EstDataStat->GetBinContent(i,i)));
    colors.push_back(special_color[kEstDataStat]);
    legs.push_back("EstDataStat");

    h_FE.push_back((TH1D*)f_in_hist->Get("CV/Truth/h_Signal")->Clone("h_FE_Total"));
    for(int i=0;i<h_FE.back()->GetNbinsX()+2;i++) h_FE.back()->SetBinContent(i,sqrt(h_FCov->GetBinContent(i,i)));
    colors.push_back(1);
    legs.push_back("Total");

    pfs::DrawUnstacked(h_FE,colors,legs,draw_overflow,draw_underflow,plot_dir+"FE.png");

  }

}
