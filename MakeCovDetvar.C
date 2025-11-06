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

void MakeCovDetvar(){

  const double scale = 1.0;

  const double data_POT = 1.332E+20;
  const double data_Trig = 31582916.0;

  // Label and set the branches defining the selection and systematics
  std::string label = "RecoW";
  double* var = &W_lt; 
  bool* sel = &sel_h8;

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";

  std::string file_CV = "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_cv_surprise_reco2_hist.root";
  double pot_CV = 5.42073e+20;
  double weight_CV = data_POT/pot_CV;

  std::vector<std::string> files_Vars = {
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lya_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyd_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyr_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_WMX_surprise_TEST_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_SCE_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_recomb2_surprise_reco2_hist.root" 
  };

  std::vector<double> pot_Vars = {
    1.11917e+20,
    1.19862e+20,
    1.17663e+20,
    1.19177e+20,
    1.19752e+20,
    1.14542e+20
  };

  std::vector<double> weight_Vars;
  for(size_t i_f=0;i_f<pot_Vars.size();i_f++) weight_Vars.push_back(data_POT/pot_Vars.at(i_f));

  // Giant contained for all of the systematics
  std::vector<std::vector<TH1D*>> h_Vars;
  std::vector<TH1D*> h_CV;

  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_CV.push_back(new TH1D(("h_Detvar_CV_"+categories.at(i_c)).c_str(),";Reco W (GeV);Events",50,0.95,5.0));
    h_Vars.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      h_Vars.back().push_back(new TH1D(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)).c_str(),";Reco W (GeV);Events",50,0.95,5.0));
    }
  }


  // Fill the CV histogram
  TFile* f_in_CV = nullptr;
  TTree* t_in_CV = nullptr;
  LoadTreeFiltered(in_dir+file_CV,f_in_CV,t_in_CV,false,false);
  for(int ievent=0;ievent<t_in_CV->GetEntries();ievent++){
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_CV->GetEntries() << std::endl;
    t_in_CV->GetEntry(ievent);
    if(!*sel) continue;
    h_CV.at(category)->Fill(*var,weight_CV);
  }
  f_in_CV->Close();

  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    TFile* f_in_Vars = nullptr;
    TTree* t_in_Vars = nullptr;
    LoadTreeFiltered(in_dir+files_Vars.at(i_s),f_in_Vars,t_in_Vars,false,false);
    for(int ievent=0;ievent<t_in_Vars->GetEntries();ievent++){
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_Vars->GetEntries() << std::endl;
      t_in_Vars->GetEntry(ievent);
      if(!*sel) continue;
      h_Vars.at(category).at(i_s)->Fill(*var,weight_Vars.at(i_s));
    }
   f_in_Vars->Close();
  }

  
  gSystem->Exec(("mkdir -p rootfiles/"+label).c_str());
  TFile* f_out = TFile::Open(("rootfiles/"+label+"/Detvars.root").c_str(),"RECREATE");

  for(size_t i_c=0;i_c<categories.size();i_c++)
    h_CV.at(i_c)->Write();

  for(size_t i_c=0;i_c<categories.size();i_c++)
    for(int i_s=0;i_s<kDetvarMAX;i_s++)
        h_Vars.at(i_c).at(i_s)->Write();        

  std::vector<std::vector<TH2D*>> h_Cov; 
  std::vector<std::vector<TH2D*>> h_FCov; 

  for(size_t i_c=0;i_c<categories.size();i_c++){
     h_Cov.push_back(std::vector<TH2D*>());
     h_FCov.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      h_Cov.back().push_back(nullptr);
      h_FCov.back().push_back(nullptr);
      CalcCovUnisim(categories.at(i_c)+"_"+detvar_str.at(i_s),h_CV.at(i_c),h_Vars.at(i_c).at(i_s),h_Cov.back().back(),h_FCov.back().back());
      h_Cov.back().back()->Write();
      h_FCov.back().back()->Write();
    }
  }       

  f_out->Close();


}
