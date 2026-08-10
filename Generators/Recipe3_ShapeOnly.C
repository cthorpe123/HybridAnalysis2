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
// assume the background prediction can be factorised out, and that the
// response can be handled using only the first generator

void Recipe3_ShapeOnly(){

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","LeadProtonKE","ProtonKE"};
  //std::vector<std::string> vars = var_names;
  vars.push_back("Enu");
  vars.push_back("Norm");
  std::vector<std::string> generators = {"Untunedv3.0.6","v3.0.6","NuWro","GiBUU"};
  bool add_detvars = false;
  bool draw_o = false;
  bool draw_u = false;

  for(const std::string& var : vars){
    std::cout << var << std::endl;

    std::string plot_dir = "Analysis/"+var+"/Plots/Recipe3_ShapeOnly/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+var+"/rootfiles/FFGenerators.root").c_str());
    TFile* f_out = new TFile(("Analysis/"+var+"/rootfiles/Recipe3_ShapeOnly.root").c_str(),"RECREATE");

    TH1D* h_bgs_data = (TH1D*)f_in->Get("CV/BGSData");
    double integral = IntegralWithOU(h_bgs_data);
    h_bgs_data->Scale(1.0/integral);

    TH2D* h_cov_data_stat = (TH2D*)f_in->Get("Cov/DataStat/h_Cov"); // Cov for errors on data 
    TH2D* h_cov_bg_mc_stat = (TH2D*)f_in->Get("Cov/BGMCStat/h_Cov");
    h_cov_data_stat->Scale(1.0/integral/integral);
    h_cov_bg_mc_stat->Scale(1.0/integral/integral);

    std::vector<TH2D*> h_cov_tot; // Total cov for each generator 
    std::map<std::string,std::vector<TH2D*>> h_cov_m; // map with cov by category for each generator

    // Calculate the fractional covariance from smearing using the first generator
    std::string gen_ref = generators.at(0);
    std::map<std::string,TH2D*> h_fcov_ref;
    TH1D* h_pred_ref = (TH1D*)f_in->Get(("CV/"+gen_ref).c_str())->Clone("h_pred_ref");
    h_pred_ref->Scale(1.0/IntegralWithOU(h_pred_ref));

    for(int i_s=0;i_s<kSystMAX;i_s++){
        std::string sys = sys_str.at(i_s);
        std::vector<TH1D*> h;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          h.push_back((TH1D*)f_in->Get(("Vars/"+sys+"/"+gen_ref+"/Pred_"+std::to_string(i_u)).c_str()));
          h.back()->Scale(1.0/IntegralWithOU(h.back()));
        }
        TH2D *c,*fc;
        CalcCovMultisim(gen_ref+"_"+sys,h,c,fc);
        h_fcov_ref[sys] = fc;
      }

      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string sys = unisims_str.at(i_s);
        TH1D* h = (TH1D*)f_in->Get(("Vars/"+sys+"/"+gen_ref+"/Pred").c_str());
        h->Scale(1.0/IntegralWithOU(h));
        TH2D *c,*fc;
        CalcCovUnisim(gen_ref+"_"+sys,h_pred_ref,h,c,fc);
        h_fcov_ref[sys] = fc;
      }

    // Calculate the stat errors
    for(std::string gen : generators){

      f_out->cd();
      f_out->mkdir(gen.c_str());
      f_out->cd(gen.c_str());

      TH1D* h_pred = (TH1D*)f_in->Get(("CV/"+gen).c_str());
      double integral_pred = IntegralWithOU(h_pred);
      h_pred->Scale(1.0/integral_pred);

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
        TH2D* c = (TH2D*)h_fcov_ref.at(sys)->Clone((sys+"_"+gen).c_str());
        for(int i_b=0;i_b<c->GetNbinsX()+2;i_b++)
          for(int j_b=0;j_b<c->GetNbinsX()+2;j_b++)
            c->SetBinContent(i_b,j_b,c->GetBinContent(i_b,j_b)*h_pred->GetBinContent(i_b)*h_pred->GetBinContent(j_b));
        
        TH2D* c_bg = (TH2D*)f_in->Get(("Cov/"+sys+"/BGSData/Cov_BGSData").c_str());
        c_bg->Scale(1.0/integral/integral);
        c->Add(c_bg);
        h_cov_tot.back()->Add(c);
        h_cov_m[sys].push_back(c);
        h_cov_m[sys].back()->Write(("Cov_"+sys).c_str());
      }

      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string sys = unisims_str.at(i_s);
        TH2D* c = (TH2D*)h_fcov_ref.at(sys)->Clone((sys+"_"+gen).c_str());
        for(int i_b=0;i_b<c->GetNbinsX()+2;i_b++)
          for(int j_b=0;j_b<c->GetNbinsX()+2;j_b++)
            c->SetBinContent(i_b,j_b,c->GetBinContent(i_b,j_b)*h_pred->GetBinContent(i_b)*h_pred->GetBinContent(j_b));
        TH2D* c_bg = (TH2D*)f_in->Get(("Cov/"+sys+"/BGSData/Cov_BGSData").c_str());
        c_bg->Scale(1.0/integral/integral);
        c->Add(c_bg);
        h_cov_tot.back()->Add(c);
        h_cov_m[sys].push_back(c);
        h_cov_m[sys].back()->Write(("Cov_"+sys).c_str());
      }      

      h_cov_tot.back()->Write("Cov_Total");

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
