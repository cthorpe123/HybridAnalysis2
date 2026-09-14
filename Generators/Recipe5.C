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

// Stitch two histograms together into a single histogram: bin 1 of the
// output holds h1's bin 1, bin 2 holds h1's bin 2, and so on, followed by
// h2's bins in the same way. This is a purely numerical concatenation,
// not a physical x axis, so the output just uses unit-width bins (0,1,2,...)
// numbered in that order. Each input's underflow/overflow bins are carried
// over too, as extra bins flanking its block of real bins.
TH1D* StitchHistograms(TH1D* h1,TH1D* h2,std::string name){

  int n1 = h1->GetNbinsX();
  int n2 = h2->GetNbinsX();

  TH1D* h_out = new TH1D(name.c_str(),"",n1+n2+2,0.5,n1+n2+2+0.5);

  //std::cout << "n1=" << n1 << " n2=" << n2 << std::endl;

  for(int i=0;i<n1+2;i++){
    //std::cout << i << " " << i << std::endl;
    h_out->SetBinContent(i,h1->GetBinContent(i));
    h_out->SetBinError(i,h1->GetBinError(i));
  }
  for(int i=0;i<n2+2;i++){
    //std::cout << i << " " << n1+i+2 << std::endl;
    h_out->SetBinContent(n1+i+2,h2->GetBinContent(i));
    h_out->SetBinError(n1+i+2,h2->GetBinError(i));
  }

  return h_out;
}

// Calculate the forward folded cross sections and chi2s using the "right"
// method, use the full set of responses in calculating every systematic
// Combine the background predictions with the generators

void Recipe5(){

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

    std::string plot_dir = "Analysis/"+var+"/Plots/Recipe5/";
    gSystem->Exec(("mkdir -p " + plot_dir).c_str());

    TFile* f_in = TFile::Open(("Analysis/"+var+"/rootfiles/FFGenerators.root").c_str());
    TFile* f_out = new TFile(("Analysis/"+var+"/rootfiles/Recipe5.root").c_str(),"RECREATE");

    TH2D* h_cov_data_stat = (TH2D*)f_in->Get("Cov/DataStat/h_Cov"); // Cov for errors on data 
    TH2D* h_cov_bg_mc_stat = (TH2D*)f_in->Get("Cov/BGMCStat/h_Cov");

    std::vector<TH2D*> h_fcov_tot_stitch; // Total fractional cov on the stitched (BGSData+Pred) histogram for each generator
    std::map<std::string,std::vector<TH2D*>> h_fcov_m_stitch; // map with stitched fractional cov by category for each generator

    for(std::string gen : generators){

      f_out->cd();
      f_out->mkdir(gen.c_str());
      f_out->cd(gen.c_str());

      const TH1D* h_bgs_data = (TH1D*)f_in->Get("CV/BGSData");
      const TH1D* h_pred = (TH1D*)f_in->Get(("CV/"+gen).c_str());

      h_pred->Write("Pred");
      h_bgs_data->Write("BGSData");

      // Size this to match the stitched (BGSData+Pred) histograms used
      // throughout this loop, rather than the plain BGSData binning
      int n1 = h_bgs_data->GetNbinsX();
      int n_stitch = n1+h_pred->GetNbinsX()+2;
      h_fcov_tot_stitch.push_back(new TH2D(("h_FCov_Tot_Stitch_"+gen).c_str(),"",n_stitch,0.5,n_stitch+0.5,n_stitch,0.5,n_stitch+0.5));

      // DataStat and BGMCStat are stat-only uncertainties on the BGSData
      // part of the stitched vector, so their binning matches h_bgs_data
      // rather than the full stitched histogram. Build a stitched-size
      // clone of h_fcov_tot_stitch and fill in only the diagonal of the
      // first (BGSData) block, leaving the rest (including the Pred block)
      // zero. The diagonal is normalised by h_bgs_data to turn the absolute
      // DataStat/BGMCStat covariance into a fractional one.
      TH2D* h_fcov_data_stat_stitch = (TH2D*)h_fcov_tot_stitch.back()->Clone(("h_FCov_DataStat_Stitch_"+gen).c_str());
      for(int i=0;i<=n1+1;i++){
        double x = h_bgs_data->GetBinContent(i);
        if(std::abs(x) > 0)
          h_fcov_data_stat_stitch->SetBinContent(i,i,h_cov_data_stat->GetBinContent(i,i)/x/x);
      }
      h_fcov_m_stitch["DataStat"].push_back(h_fcov_data_stat_stitch);
      h_fcov_tot_stitch.back()->Add(h_fcov_data_stat_stitch);

      TH2D* h_fcov_bg_mc_stat_stitch = (TH2D*)h_fcov_tot_stitch.back()->Clone(("h_FCov_BGMCStat_Stitch_"+gen).c_str());
      for(int i=0;i<=n1+1;i++){
        double x = h_bgs_data->GetBinContent(i);
        if(std::abs(x) > 0)
          h_fcov_bg_mc_stat_stitch->SetBinContent(i,i,h_cov_bg_mc_stat->GetBinContent(i,i)/x/x);
      }
      h_fcov_m_stitch["BGMCStat"].push_back(h_fcov_bg_mc_stat_stitch);
      h_fcov_tot_stitch.back()->Add(h_fcov_bg_mc_stat_stitch);

      // Multisims
      for(int i_s=0;i_s<kSystMAX;i_s++){
        std::string sys = sys_str.at(i_s);

        std::string plot_dir_sys = plot_dir+sys+"/";
        gSystem->Exec(("mkdir -p " + plot_dir_sys).c_str());

        std::vector<TH1D*> h;
        for(int i_u=0;i_u<sys_nuniv.at(i_s);i_u++){
          TH1D* h_bgsdata_i = (TH1D*)f_in->Get(("Vars/"+sys+"/BGSData/BGSData_"+std::to_string(i_u)).c_str());
          TH1D* h_pred_i = (TH1D*)f_in->Get(("Vars/"+sys+"/"+gen+"/Pred_"+std::to_string(i_u)).c_str());
          TH1D* h_stitch = StitchHistograms(h_bgsdata_i,h_pred_i,"h_stitch_"+gen+"_"+sys+"_"+std::to_string(i_u));
          h.push_back(h_stitch);
        }
        TH2D *c,*fc;
        CalcCovMultisim(gen+"_"+sys,h,c,fc);
        pfs::Draw2DHist(c,plot_dir+"Cov_"+sys+"_"+gen+".png");
        h_fcov_tot_stitch.back()->Add(fc);
        h_fcov_m_stitch[sys].push_back(fc);
      }

      // Unisims
      for(int i_s=0;i_s<kUnisimMAX;i_s++){
        std::string sys = unisims_str.at(i_s);

        TH1D* h_bgsdata_var = (TH1D*)f_in->Get(("Vars/"+sys+"/BGSData/BGSData").c_str());
        TH1D* h_pred_var = (TH1D*)f_in->Get(("Vars/"+sys+"/"+gen+"/Pred").c_str());

        TH1D* h_stitch_cv = StitchHistograms((TH1D*)h_bgs_data,(TH1D*)h_pred,"h_stitch_cv_"+gen+"_"+sys);
        TH1D* h_stitch_var = StitchHistograms(h_bgsdata_var,h_pred_var,"h_stitch_"+gen+"_"+sys);

        TH2D *c,*fc;
        CalcCovUnisim(gen+"_"+sys,h_stitch_cv,h_stitch_var,c,fc);
        pfs::Draw2DHist(c,plot_dir+"Cov_"+sys+"_"+gen+".png");
        h_fcov_tot_stitch.back()->Add(fc);
        h_fcov_m_stitch[sys].push_back(fc);
      }

      // Convert the fractional covariances back into absolute covariances.
      // Rather than re-stitching h_bgs_data and h_pred into a single CV
      // histogram, look up each stitched bin's CV directly from whichever
      // of the two it came from: bins 0..n1+1 are h_bgs_data's (including
      // its under/overflow), the rest are h_pred's (see StitchHistograms).
      auto StitchCV = [&](int i){ return i<=n1+1 ? h_bgs_data->GetBinContent(i) : h_pred->GetBinContent(i-(n1+2)); };

      std::vector<TH2D*> h_cov_tot_stitch;
      h_cov_tot_stitch.push_back((TH2D*)h_fcov_tot_stitch.back()->Clone(("h_Cov_Tot_Stitch_"+gen).c_str()));
      for(int i=0;i<n_stitch+2;i++)
        for(int j=0;j<n_stitch+2;j++)
          h_cov_tot_stitch.back()->SetBinContent(i,j,h_fcov_tot_stitch.back()->GetBinContent(i,j)*StitchCV(i)*StitchCV(j));

      std::map<std::string,std::vector<TH2D*>> h_cov_m_stitch;
      for(auto item : h_fcov_m_stitch){
        h_cov_m_stitch[item.first].push_back((TH2D*)item.second.back()->Clone(("h_Cov_"+item.first+"_Stitch_"+gen).c_str()));
        for(int i=0;i<n_stitch+2;i++)
          for(int j=0;j<n_stitch+2;j++)
            h_cov_m_stitch[item.first].back()->SetBinContent(i,j,item.second.back()->GetBinContent(i,j)*StitchCV(i)*StitchCV(j));
      }

      // Use the (now absolute) stitched covariance matrices to build the
      // covariance matrix of the residual z = h_bgs_data-h_pred. With x,y
      // the BGSData/Pred blocks of the stitched vector,
      // Cov(z_i,z_j) = Cov(x_i,x_j) - Cov(x_i,y_j) - Cov(y_i,x_j) + Cov(y_i,y_j)
      auto DiffCov = [&](const TH2D* h_cov_stitch,int i,int j){
        return h_cov_stitch->GetBinContent(i,j)
             - h_cov_stitch->GetBinContent(i,n1+2+j)
             - h_cov_stitch->GetBinContent(n1+2+i,j)
             + h_cov_stitch->GetBinContent(n1+2+i,n1+2+j);
      };

      // Use h_bgs_data's physical binning (rather than the unit-width
      // binning used for the stitched matrices above), since this matrix
      // is indexed the same way as h_bgs_data/h_pred themselves
      TH2D* h_cov_diff_tot = Make2DHist("h_Cov_Diff_Tot_"+gen,(TH1D*)h_bgs_data);
      for(int i=0;i<=n1+1;i++)
        for(int j=0;j<=n1+1;j++)
          h_cov_diff_tot->SetBinContent(i,j,DiffCov(h_cov_tot_stitch.back(),i,j));
      h_cov_diff_tot->Write("Cov_Total");

      std::map<std::string,TH2D*> h_cov_diff_m;
      for(auto item : h_cov_m_stitch){
        TH2D* h_cov_diff = Make2DHist("h_Cov_Diff_"+item.first+"_"+gen,(TH1D*)h_bgs_data);
        for(int i=0;i<=n1+1;i++)
          for(int j=0;j<=n1+1;j++)
            h_cov_diff->SetBinContent(i,j,DiffCov(item.second.back(),i,j));
        h_cov_diff_m[item.first] = h_cov_diff;
        h_cov_diff->Write(("Cov_"+item.first).c_str());
      }



      // h_cov_diff_tot/h_cov_diff_m are absolute covariances on the same
      // binning as h_pred, so MakeFEHist (which divides sqrt(diag) by
      // h_val) can be used directly, same as Recipe1-4
      std::vector<TH1D*> h_fe_v;
      std::vector<std::string> legs;
      std::vector<int> cols;
      int col = 2;
      for(auto item : h_cov_diff_m){
        pfs::Draw2DHist(item.second,plot_dir+"Cov_"+item.first+"_"+gen+".png");
        h_fe_v.push_back((TH1D*)h_pred->Clone(("h_fe_"+item.first+"_"+gen).c_str()));
        legs.push_back(item.first);
        cols.push_back(col);
        MakeFEHist(h_fe_v.back(),h_pred,item.second);
        h_fe_v.back()->Write(("FE_"+item.first).c_str());
        col++;
      }
      h_fe_v.push_back((TH1D*)h_pred->Clone(("h_fe_tot_"+gen).c_str()));
      MakeFEHist(h_fe_v.back(),h_pred,h_cov_diff_tot);
      h_fe_v.back()->Write("FE_Total");
      legs.push_back("Total");
      cols.push_back(1);
      pfs::DrawUnstacked(h_fe_v,cols,legs,draw_o,draw_u,false,false,plot_dir+"FE_"+gen+".png");
    }
      
    f_in->Close();
    f_out->Close();

  }

}
