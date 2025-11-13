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
  std::string label = "TrueE";
  std::string axis_title = ";True E (GeV);Events";
  int nbins = 50;
  double low = 0.1;
  double high = 2.5;

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";

  std::string file_CV = "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_cv_surprise_reco2_hist.root";
  double pot_CV = 5.42073e+20;
  double weight_CV = data_POT/pot_CV;

  std::vector<std::string> files_Vars = {
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lya_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyd_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyr_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_SCE_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_recomb2_surprise_reco2_hist.root" 
    //"Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_WMX_surprise_TEST_reco2_hist.root",
  };

  std::vector<double> pot_Vars = {
    1.11917e+20,
    1.19862e+20,
    1.17663e+20,
    1.19752e+20,
    1.14542e+20
    //1.19177e+20,
  };

  std::vector<double> weight_Vars;
  for(size_t i_f=0;i_f<pot_Vars.size();i_f++) weight_Vars.push_back(data_POT/pot_Vars.at(i_f));

 
  TH1D* h_CV_Tot =  new TH1D("h_Detvar_CV_Tot",axis_title.c_str(),nbins,low,high);
  std::vector<TH1D*> h_CV;

  std::vector<TH1D*> h_Vars_Tot;
  std::vector<std::vector<TH1D*>> h_Vars;
  
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    h_Vars_Tot.push_back(new TH1D(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)).c_str(),axis_title.c_str(),nbins,low,high));
 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_CV.push_back(new TH1D(("h_Detvar_CV_"+categories.at(i_c)).c_str(),axis_title.c_str(),nbins,low,high));
    h_Vars.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      h_Vars.back().push_back(new TH1D(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)).c_str(),axis_title.c_str(),nbins,low,high));
    }
  }


  // Fill the CV histogram
  TFile* f_in_CV = nullptr;
  TTree* t_in_CV = nullptr;
  LoadTreeFiltered(in_dir+file_CV,f_in_CV,t_in_CV,false,false);
  for(int ievent=0;ievent<t_in_CV->GetEntries();ievent++){
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_CV->GetEntries() << std::endl;
    t_in_CV->GetEntry(ievent);
    if(is_data || is_ext || is_dirt) continue;
    double var = nu_e;
    if(!is_data) h_CV_Tot->Fill(var,weight_CV);
    h_CV.at(category)->Fill(var,weight_CV);
  }
  f_in_CV->Close();

  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    TFile* f_in_Vars = nullptr;
    TTree* t_in_Vars = nullptr;
    LoadTreeFiltered(in_dir+files_Vars.at(i_s),f_in_Vars,t_in_Vars,false,false);
    for(int ievent=0;ievent<t_in_Vars->GetEntries();ievent++){
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_Vars->GetEntries() << std::endl;
      t_in_Vars->GetEntry(ievent);
      if(is_data || is_ext || is_dirt) continue;
      double var = nu_e;
      h_Vars_Tot.at(i_s)->Fill(var,weight_Vars.at(i_s));
      h_Vars.at(category).at(i_s)->Fill(var,weight_Vars.at(i_s));
    }
   f_in_Vars->Close();
  }

  
  gSystem->Exec(("mkdir -p Analysis/"+label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str(),"RECREATE");

  h_CV_Tot->Write();
 
  for(size_t i_c=0;i_c<categories.size();i_c++)
    h_CV.at(i_c)->Write();

  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    h_Vars_Tot.at(i_s)->Write();        

  for(size_t i_c=0;i_c<categories.size();i_c++)
    for(int i_s=0;i_s<kDetvarMAX;i_s++)
      h_Vars.at(i_c).at(i_s)->Write();        

  std::vector<TH2D*> h_Cov_Tot; 
  std::vector<TH2D*> h_FCov_Tot; 
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    h_Cov_Tot.push_back(nullptr);
    h_FCov_Tot.push_back(nullptr);
    CalcCovUnisim(detvar_str.at(i_s),h_CV_Tot,h_Vars_Tot.at(i_s),h_Cov_Tot.back(),h_FCov_Tot.back());
    h_Cov_Tot.back()->Write();
    h_FCov_Tot.back()->Write();
  }
  
  TH2D* h_Cov = static_cast<TH2D*>(h_Cov_Tot.at(0)->Clone("Cov"));
  TH2D* h_FCov = static_cast<TH2D*>(h_FCov_Tot.at(0)->Clone("FCov"));
  h_Cov->Reset();
  h_FCov->Reset();
  for(int i_s=0;i_s<kSystMAX;i_s++){
    h_Cov->Add(h_Cov_Tot.at(i_s));    
    h_FCov->Add(h_FCov_Tot.at(i_s));    
  }
  h_Cov->Write();
  h_FCov->Write();

  std::vector<std::vector<TH2D*>> h_Cov_Cat; 
  std::vector<std::vector<TH2D*>> h_FCov_Cat; 

  for(size_t i_c=0;i_c<categories.size();i_c++){
     h_Cov_Cat.push_back(std::vector<TH2D*>());
     h_FCov_Cat.push_back(std::vector<TH2D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      h_Cov_Cat.back().push_back(nullptr);
      h_FCov_Cat.back().push_back(nullptr);
      CalcCovUnisim(categories.at(i_c)+"_"+detvar_str.at(i_s),h_CV.at(i_c),h_Vars.at(i_c).at(i_s),h_Cov_Cat.back().back(),h_FCov_Cat.back().back());
      h_Cov_Cat.back().back()->Write();
      h_FCov_Cat.back().back()->Write();
    }
  }       

  f_out->Close();


}
