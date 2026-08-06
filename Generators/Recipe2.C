#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "BinningFuncs.h"
#include "PlotFuncs.h"
#include "MultiChannelHistograms.h"

using namespace syst;

// Calculate the forward folded cross sections and covariances, this time
// assume the background prediction can be factorised out

void Recipe2(){

  std::vector<std::string> vars = {"Norm","Enu","MuonMom","MuonCosTheta"};
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/Recipe2/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+var+"/rootfiles/FFGenerators.root").c_str());
    TFile* f_out = new TFile(("Analysis/"+var+"/rootfiles/Recipe2.root").c_str(),"RECREATE");

    TH2D* h_cov_data_stat = (TH2D*)f_in->Get("Cov/DataStat/h_Cov"); // Cov for errors on data 
    TH2D* h_cov_bg_mc_stat = (TH2D*)f_in->Get("Cov/BGMCStat/h_Cov");

    std::vector<TH2D*> h_cov_tot; // Total cov for each generator 
    std::map<std::string,std::vector<TH2D*>> h_cov_m; // map with cov by category for each generator

    // Calculate the stat errors
    for(std::string gen : generators){

      f_out->cd();
      f_out->mkdir(gen.c_str());
      f_out->cd(gen.c_str());

      const TH1D* h_bgs_data = (TH1D*)f_in->Get("CV/BGSData");
      const TH1D* h_pred = (TH1D*)f_in->Get(("CV/"+gen).c_str());

      h_pred->Write("Pred");
      h_bgs_data->Write("BGSData");

      h_cov_tot.push_back((TH2D*)f_in->Get("Cov/DataStat/h_Cov")->Clone(("h_Cov_Tot_"+gen).c_str()));
      h_cov_tot.back()->Reset();

      h_cov_m["DataStat"].push_back(h_cov_data_stat);
      h_cov_m["DataStat"].back()->Write("Cov_DataStat");
      h_cov_tot.back()->Add(h_cov_data_stat);

      h_cov_m["BGMCStat"].push_back(h_cov_bg_mc_stat);
      h_cov_m["BGMCStat"].back()->Write("Cov_BGMCStat");
      h_cov_tot.back()->Add(h_cov_bg_mc_stat);

      // Multisims
      for(int i_s=0;i_s<kSystMAX;i_s++){
        std::string sys = sys_str.at(i_s);
        std::vector<TH1D*> h;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          h.push_back((TH1D*)f_in->Get(("Vars/"+sys+"/"+gen+"/Pred_"+std::to_string(i_u)).c_str()));
        }
        TH2D *c,*fc;
        CalcCovMultisim(gen+"_"+sys,h,c,fc);
        c->Add((TH2D*)f_in->Get(("Cov/"+sys+"/BG/Cov_BG").c_str()));
        h_cov_tot.back()->Add(c);
        h_cov_m[sys].push_back(c);
        h_cov_m[sys].back()->Write(("Cov_"+sys).c_str());
      }
      
      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string sys = unisims_str.at(i_s);
        TH1D* h = (TH1D*)f_in->Get(("Vars/"+sys+"/"+gen+"/Pred").c_str());
        TH2D *c,*fc;
        CalcCovUnisim(gen+"_"+sys,h_pred,h,c,fc);
        c->Add((TH2D*)f_in->Get(("Cov/"+sys+"/BG/Cov_BG").c_str()));
        h_cov_tot.back()->Add(c);
        h_cov_m[sys].push_back(c);
        h_cov_m[sys].back()->Write(("Cov_"+sys).c_str());
      }
      
      h_cov_tot.back()->Write("Cov_Tot");

      std::vector<TH1D*> h_fe_v;
      std::vector<std::string> legs;
      std::vector<int> cols;
      int col = 2;
      for(auto item : h_cov_m){
        pfs::Draw2DHist(item.second.back(),plot_dir+"Cov_"+item.first+"_"+gen+".png");
        h_fe_v.push_back((TH1D*)h_pred->Clone(("h_fe_"+item.first+"_"+gen).c_str()));
        legs.push_back(item.first);
        cols.push_back(col);
        MakeFEHist(h_fe_v.back(),h_pred,item.second.back());
        h_fe_v.back()->Write(("FE_"+item.first).c_str());
        col++;
      }
      h_fe_v.push_back((TH1D*)h_pred->Clone(("h_fe_tot_"+gen).c_str()));
      MakeFEHist(h_fe_v.back(),h_pred,h_cov_tot.back());
      h_fe_v.back()->Write("FE_Total");
      legs.push_back("Total");
      cols.push_back(1);
      pfs::DrawUnstacked(h_fe_v,cols,legs,draw_o,draw_u,false,false,plot_dir+"FE_"+gen+".png");

    } 

    f_in->Close();
    f_out->Close();

  }

}
