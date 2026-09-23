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
// Same as Recipe6, which builds Cov(x-y) = Cov(x,x) - Cov(x,y) - Cov(y,x) +
// Cov(y,y) directly out of the auto/cross covariance blocks for every
// generator. Here, following the approach used in Recipe4, the four blocks
// are only computed as fractional covariances using the first generator in
// the list, then rescaled using the CV predictions of the other generators
// to get their residual covariance matrices. Assume the FCov blocks from
// the first generator are valid for the others.

void Recipe7(){

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","LeadProtonKE","ProtonKE"};
  //std::vector<std::string> vars = var_names;
  //vars.push_back("Enu");
  //vars.push_back("Norm");
  std::vector<std::string> generators = {"v3.0.6","Untunedv3.0.6","NuWro","GiBUU"};
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = AnalysisDir()+"/"+var+"/Plots/Recipe7/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile* f_in = TFile::Open((AnalysisDir()+"/"+var+"/rootfiles/FFGenerators.root").c_str());
    TFile* f_out = new TFile((AnalysisDir()+"/"+var+"/rootfiles/Recipe7.root").c_str(),"RECREATE");

    TH2D* h_cov_data_stat = (TH2D*)f_in->Get("Cov/DataStat/h_Cov"); // Cov for errors on data
    TH2D* h_cov_bg_mc_stat = (TH2D*)f_in->Get("Cov/BGMCStat/h_Cov");

    std::vector<TH2D*> h_cov_diff_tot; // Total cov on the BGSData-Pred residual for each generator
    std::map<std::string,std::vector<TH2D*>> h_cov_diff_m; // map with residual cov by category for each generator

    // Calculate the fractional residual covariance blocks using the first generator
    std::string gen_ref = generators.at(0);
    TH1D* h_bgs_data_ref = (TH1D*)f_in->Get("CV/BGSData");
    TH1D* h_pred_ref = (TH1D*)f_in->Get(("CV/"+gen_ref).c_str());

    std::map<std::string,TH2D*> h_fcov_xx,h_fcov_yy,h_fcov_xy,h_fcov_yx;

    // Multisims
    for(int i_s=0;i_s<kSystMAX;i_s++){
      std::string sys = sys_str.at(i_s);

      std::vector<TH1D*> h_bgsdata_v;
      std::vector<TH1D*> h_pred_v;
      for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
        h_bgsdata_v.push_back((TH1D*)f_in->Get(("Vars/"+sys+"/BGSData/BGSData_"+std::to_string(i_u)).c_str()));
        h_pred_v.push_back((TH1D*)f_in->Get(("Vars/"+sys+"/"+gen_ref+"/Pred_"+std::to_string(i_u)).c_str()));
      }

      TH2D *c_xx,*fc_xx,*c_yy,*fc_yy,*c_xy,*fc_xy,*c_yx,*fc_yx;
      CalcCovMultisim(gen_ref+"_"+sys+"_BGSData",h_bgsdata_v,c_xx,fc_xx);
      CalcCovMultisim(gen_ref+"_"+sys+"_Pred",h_pred_v,c_yy,fc_yy);
      CalcCovMultisimBlock(gen_ref+"_"+sys+"_XY",h_bgsdata_v,h_pred_v,c_xy,fc_xy);
      CalcCovMultisimBlock(gen_ref+"_"+sys+"_YX",h_pred_v,h_bgsdata_v,c_yx,fc_yx);

      std::string plot_dir_sys = plot_dir+sys+"/";
      gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());
      pfs::Draw2DHist(fc_xx,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_xx.png");
      pfs::Draw2DHist(fc_yy,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_yy.png");
      pfs::Draw2DHist(fc_xy,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_xy.png");
      pfs::Draw2DHist(fc_yx,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_yx.png");

      h_fcov_xx[sys] = fc_xx;
      h_fcov_yy[sys] = fc_yy;
      h_fcov_xy[sys] = fc_xy;
      h_fcov_yx[sys] = fc_yx;
    }

    // Unisims
    for(int i_s=0;i_s<kUnisimMAX;i_s++){
      std::string sys = unisims_str.at(i_s);

      TH1D* h_bgsdata_var = (TH1D*)f_in->Get(("Vars/"+sys+"/BGSData/BGSData").c_str());
      TH1D* h_pred_var = (TH1D*)f_in->Get(("Vars/"+sys+"/"+gen_ref+"/Pred").c_str());

      TH2D *c_xx,*fc_xx,*c_yy,*fc_yy,*c_xy,*fc_xy,*c_yx,*fc_yx;
      CalcCovUnisim(gen_ref+"_"+sys+"_BGSData",h_bgs_data_ref,h_bgsdata_var,c_xx,fc_xx);
      CalcCovUnisim(gen_ref+"_"+sys+"_Pred",h_pred_ref,h_pred_var,c_yy,fc_yy);
      CalcCovUnisimBlock(gen_ref+"_"+sys+"_XY",h_bgs_data_ref,h_bgsdata_var,h_pred_ref,h_pred_var,c_xy,fc_xy);
      CalcCovUnisimBlock(gen_ref+"_"+sys+"_YX",h_pred_ref,h_pred_var,h_bgs_data_ref,h_bgsdata_var,c_yx,fc_yx);

      std::string plot_dir_sys = plot_dir+sys+"/";
      gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());
      pfs::Draw2DHist(fc_xx,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_xx.png");
      pfs::Draw2DHist(fc_yy,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_yy.png");
      pfs::Draw2DHist(fc_xy,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_xy.png");
      pfs::Draw2DHist(fc_yx,plot_dir_sys+"FCov_Ref_"+sys+"_"+gen_ref+"_yx.png");

      h_fcov_xx[sys] = fc_xx;
      h_fcov_yy[sys] = fc_yy;
      h_fcov_xy[sys] = fc_xy;
      h_fcov_yx[sys] = fc_yx;
    }

    std::vector<std::string> all_sys;
    for(int i_s=0;i_s<kSystMAX;i_s++) all_sys.push_back(sys_str.at(i_s));
    for(int i_s=0;i_s<kUnisimMAX;i_s++) all_sys.push_back(unisims_str.at(i_s));

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

      // Rescale the reference fractional covariance blocks using the CV
      // BGSData/Pred for this generator to get its residual covariance
      for(const std::string& sys : all_sys){
        TH2D* fc_xx = h_fcov_xx.at(sys);
        TH2D* fc_yy = h_fcov_yy.at(sys);
        TH2D* fc_xy = h_fcov_xy.at(sys);
        TH2D* fc_yx = h_fcov_yx.at(sys);

        std::string plot_dir_sys = plot_dir+sys+"/";
        gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());

        TH2D* c_xx = Make2DHist("h_Cov_"+sys+"_"+gen+"_xx",h_bgs_data);
        TH2D* c_yy = Make2DHist("h_Cov_"+sys+"_"+gen+"_yy",h_bgs_data);
        TH2D* c_xy = Make2DHist("h_Cov_"+sys+"_"+gen+"_xy",h_bgs_data);
        TH2D* c_yx = Make2DHist("h_Cov_"+sys+"_"+gen+"_yx",h_bgs_data);
        TH2D* c_diff = Make2DHist("h_Cov_"+sys+"_"+gen,h_bgs_data);
        for(int i=0;i<h_bgs_data->GetNbinsX()+2;i++){
          for(int j=0;j<h_bgs_data->GetNbinsX()+2;j++){
            double v_xx = fc_xx->GetBinContent(i,j)*h_bgs_data->GetBinContent(i)*h_bgs_data->GetBinContent(j);
            double v_yy = fc_yy->GetBinContent(i,j)*h_pred->GetBinContent(i)*h_pred->GetBinContent(j);
            double v_xy = fc_xy->GetBinContent(i,j)*h_bgs_data->GetBinContent(i)*h_pred->GetBinContent(j);
            double v_yx = fc_yx->GetBinContent(i,j)*h_pred->GetBinContent(i)*h_bgs_data->GetBinContent(j);
            c_xx->SetBinContent(i,j,v_xx);
            c_yy->SetBinContent(i,j,v_yy);
            c_xy->SetBinContent(i,j,v_xy);
            c_yx->SetBinContent(i,j,v_yx);
            c_diff->SetBinContent(i,j,v_xx - v_xy - v_yx + v_yy);
          }
        }

        pfs::Draw2DHist(c_xx,plot_dir_sys+"Cov_"+sys+"_"+gen+"_xx.png");
        pfs::Draw2DHist(c_yy,plot_dir_sys+"Cov_"+sys+"_"+gen+"_yy.png");
        pfs::Draw2DHist(c_xy,plot_dir_sys+"Cov_"+sys+"_"+gen+"_xy.png");
        pfs::Draw2DHist(c_yx,plot_dir_sys+"Cov_"+sys+"_"+gen+"_yx.png");
        pfs::Draw2DHist(c_diff,plot_dir_sys+"Cov_"+sys+"_"+gen+".png");
        h_cov_diff_tot.back()->Add(c_diff);
        h_cov_diff_m[sys].push_back(c_diff);
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
