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

// Calculate the forward folded cross sections and chi2s using the "right"
// method, use the full set of responses in calculating every systematic
// Combine the background predictions with the generators
//
// Same as Recipe5, but instead of stitching the BGSData and Pred histograms
// together to get the covariance of their difference from a single call to
// CalcCovMultisim, build it directly out of the auto/cross covariance
// blocks: Cov(x-y) = Cov(x,x) - Cov(x,y) - Cov(y,x) + Cov(y,y), using
// CalcCovMultisim for the auto blocks and the new CalcCovMultisimBlock for
// the cross blocks.

void Recipe6(){

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","LeadProtonKE","ProtonKE"};
  //std::vector<std::string> vars = var_names;
  //vars.push_back("Enu");
  //vars.push_back("Norm");
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/Recipe6/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+var+"/rootfiles/FFGenerators.root").c_str());
    TFile* f_out = new TFile(("Analysis/"+var+"/rootfiles/Recipe6.root").c_str(),"RECREATE");

    TH2D* h_cov_data_stat = (TH2D*)f_in->Get("Cov/DataStat/h_Cov"); // Cov for errors on data
    TH2D* h_cov_bg_mc_stat = (TH2D*)f_in->Get("Cov/BGMCStat/h_Cov");

    std::vector<TH2D*> h_cov_diff_tot; // Total cov on the BGSData-Pred residual for each generator
    std::map<std::string,std::vector<TH2D*>> h_cov_diff_m; // map with residual cov by category for each generator

    for(std::string gen : generators){

      f_out->cd();
      f_out->mkdir(gen.c_str());
      f_out->cd(gen.c_str());

      TH1D* h_bgs_data = (TH1D*)f_in->Get("CV/BGSData");
      TH1D* h_pred = (TH1D*)f_in->Get(("CV/"+gen).c_str());

      h_cov_diff_tot.push_back(Make2DHist("h_Cov_Diff_Tot_"+gen,h_bgs_data));

      // DataStat and BGMCStat are stat-only uncertainties on the BGSData
      // side of the residual, Pred has no corresponding stat uncertainty
      // here, so these blocks carry over unchanged
      h_cov_diff_m["DataStat"].push_back((TH2D*)h_cov_data_stat->Clone(("h_Cov_DataStat_"+gen).c_str()));
      h_cov_diff_tot.back()->Add(h_cov_data_stat);

      h_cov_diff_m["BGMCStat"].push_back((TH2D*)h_cov_bg_mc_stat->Clone(("h_Cov_BGMCStat_"+gen).c_str()));
      h_cov_diff_tot.back()->Add(h_cov_bg_mc_stat);

      // Multisims
      for(int i_s=0;i_s<kSystMAX;i_s++){
        std::string sys = sys_str.at(i_s);

        std::string plot_dir_sys = plot_dir+sys+"/";
        gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());

        std::vector<TH1D*> h_bgsdata_v;
        std::vector<TH1D*> h_pred_v;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          h_bgsdata_v.push_back((TH1D*)f_in->Get(("Vars/"+sys+"/BGSData/BGSData_"+std::to_string(i_u)).c_str()));
          h_pred_v.push_back((TH1D*)f_in->Get(("Vars/"+sys+"/"+gen+"/Pred_"+std::to_string(i_u)).c_str()));
        }

        TH2D *c_xx,*fc_xx,*c_yy,*fc_yy,*c_xy,*fc_xy,*c_yx,*fc_yx;
        CalcCovMultisim(gen+"_"+sys+"_BGSData",h_bgsdata_v,c_xx,fc_xx);
        CalcCovMultisim(gen+"_"+sys+"_Pred",h_pred_v,c_yy,fc_yy);
        CalcCovMultisimBlock(gen+"_"+sys+"_XY",h_bgsdata_v,h_pred_v,c_xy,fc_xy);
        CalcCovMultisimBlock(gen+"_"+sys+"_YX",h_pred_v,h_bgsdata_v,c_yx,fc_yx);

        pfs::Draw2DHist(c_xx,plot_dir_sys+"Cov_"+sys+"_"+gen+"_xx.png");
        pfs::Draw2DHist(c_yy,plot_dir_sys+"Cov_"+sys+"_"+gen+"_yy.png");
        pfs::Draw2DHist(c_xy,plot_dir_sys+"Cov_"+sys+"_"+gen+"_xy.png");
        pfs::Draw2DHist(c_yx,plot_dir_sys+"Cov_"+sys+"_"+gen+"_yx.png");

        TH2D* c_diff = Make2DHist("h_Cov_"+sys+"_"+gen,h_bgs_data);
        for(int i=0;i<h_bgs_data->GetNbinsX()+2;i++)
          for(int j=0;j<h_bgs_data->GetNbinsX()+2;j++)
            c_diff->SetBinContent(i,j,c_xx->GetBinContent(i,j) - c_xy->GetBinContent(i,j) - c_yx->GetBinContent(i,j) + c_yy->GetBinContent(i,j));

        pfs::Draw2DHist(c_diff,plot_dir_sys+"Cov_"+sys+"_"+gen+".png");
        h_cov_diff_tot.back()->Add(c_diff);
        h_cov_diff_m[sys].push_back(c_diff);
      }

      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string sys = unisims_str.at(i_s);

        std::string plot_dir_sys = plot_dir+sys+"/";
        gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());

        TH1D* h_bgsdata_var = (TH1D*)f_in->Get(("Vars/"+sys+"/BGSData/BGSData").c_str());
        TH1D* h_pred_var = (TH1D*)f_in->Get(("Vars/"+sys+"/"+gen+"/Pred").c_str());

        TH2D *c_xx,*fc_xx,*c_yy,*fc_yy,*c_xy,*fc_xy,*c_yx,*fc_yx;
        CalcCovUnisim(gen+"_"+sys+"_BGSData",h_bgs_data,h_bgsdata_var,c_xx,fc_xx);
        CalcCovUnisim(gen+"_"+sys+"_Pred",h_pred,h_pred_var,c_yy,fc_yy);
        CalcCovUnisimBlock(gen+"_"+sys+"_XY",h_bgs_data,h_bgsdata_var,h_pred,h_pred_var,c_xy,fc_xy);
        CalcCovUnisimBlock(gen+"_"+sys+"_YX",h_pred,h_pred_var,h_bgs_data,h_bgsdata_var,c_yx,fc_yx);

        pfs::Draw2DHist(c_xx,plot_dir_sys+"Cov_"+sys+"_"+gen+"_xx.png");
        pfs::Draw2DHist(c_yy,plot_dir_sys+"Cov_"+sys+"_"+gen+"_yy.png");
        pfs::Draw2DHist(c_xy,plot_dir_sys+"Cov_"+sys+"_"+gen+"_xy.png");
        pfs::Draw2DHist(c_yx,plot_dir_sys+"Cov_"+sys+"_"+gen+"_yx.png");

        TH2D* c = Make2DHist("h_Cov_"+sys+"_"+gen,h_bgs_data);
        for(int i=0;i<h_bgs_data->GetNbinsX()+2;i++)
          for(int j=0;j<h_bgs_data->GetNbinsX()+2;j++)
            c->SetBinContent(i,j,c_xx->GetBinContent(i,j) - c_xy->GetBinContent(i,j) - c_yx->GetBinContent(i,j) + c_yy->GetBinContent(i,j));

        pfs::Draw2DHist(c,plot_dir_sys+"Cov_"+sys+"_"+gen+".png");
        h_cov_diff_tot.back()->Add(c);
        h_cov_diff_m[sys].push_back(c);
      }

      // Set the error on the Pred histogram from the diagonal of the total
      // residual covariance, then write it out alongside BGSData
      for(int i=0;i<=h_pred->GetNbinsX()+1;i++)
        h_pred->SetBinError(i,std::sqrt(h_cov_diff_tot.back()->GetBinContent(i,i)));
      h_pred->Write("Pred");
      h_bgs_data->Write("BGSData");

      h_cov_diff_tot.back()->Write("Cov_Total");
      for(auto item : h_cov_diff_m)
        item.second.back()->Write(("Cov_"+item.first).c_str());

      std::vector<TH1D*> h_fe_v;
      std::vector<std::string> legs;
      std::vector<int> cols;
      int col = 2;
      for(auto item : h_cov_diff_m){
        std::string plot_dir_sys = plot_dir+item.first+"/";
        gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());
        pfs::Draw2DHist(item.second.back(),plot_dir_sys+"Cov_"+item.first+"_"+gen+".png");
        h_fe_v.push_back((TH1D*)h_pred->Clone(("h_fe_"+item.first+"_"+gen).c_str()));
        legs.push_back(item.first);
        cols.push_back(col);
        MakeFEHist(h_fe_v.back(),h_pred,item.second.back());
        h_fe_v.back()->Write(("FE_"+item.first).c_str());
        col++;
      }
      h_fe_v.push_back((TH1D*)h_pred->Clone(("h_fe_tot_"+gen).c_str()));
      MakeFEHist(h_fe_v.back(),h_pred,h_cov_diff_tot.back());
      h_fe_v.back()->Write("FE_Total");
      legs.push_back("Total");
      cols.push_back(1);
      pfs::DrawUnstacked(h_fe_v,cols,legs,draw_o,draw_u,false,false,plot_dir+"FE_"+gen+".png");

    }

    f_in->Close();
    f_out->Close();

  }

}
