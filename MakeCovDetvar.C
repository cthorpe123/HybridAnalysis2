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
  std::string label = "RecoE_Test";
  std::string axis_title = ";Reco E (GeV);Events";
  std::vector<double> bins = {0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.90,0.95,1.00,1.1,1.2,1.3,1.4,1.5,1.75,2.0,2.25,2.5};

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";

  std::string file_CV = "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_cv_surprise_reco2_hist.root";
  double pot_CV = 5.42073e+20;
  double weight_CV = data_POT/pot_CV;

  std::string file_Dirt = "../Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root";
  double pot_Dirt = 3.06E+20;

  std::string file_EXT = "../Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root";
  double trig_EXT = 88445969.0;

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

  double weight_Dirt = data_POT/pot_Dirt;
  double weight_EXT = data_Trig/trig_EXT;

  std::vector<double> weight_Vars;
  for(size_t i_f=0;i_f<pot_Vars.size();i_f++) weight_Vars.push_back(data_POT/pot_Vars.at(i_f));

  TH1D* h_CV_Tot =  new TH1D("h_Detvar_CV_Tot",axis_title.c_str(),bins.size()-1,&bins[0]);
  std::vector<TH1D*> h_CV;

  std::vector<TH1D*> h_Vars_Tot;
  std::vector<std::vector<TH1D*>> h_Vars;
  
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    h_Vars_Tot.push_back(new TH1D(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    h_CV.push_back(new TH1D(("h_Detvar_CV_"+categories.at(i_c)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
    h_Vars.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      h_Vars.back().push_back(new TH1D(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
    }
  }


  // Fill the CV histogram
  std::cout << "Analysing CV" << std::endl;
  TFile* f_in_CV = nullptr;
  TTree* t_in_CV = nullptr;
  LoadTreeFiltered(in_dir+file_CV,f_in_CV,t_in_CV,false,false);
  for(int ievent=0;ievent<t_in_CV->GetEntries();ievent++){
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_CV->GetEntries() << std::endl;
    t_in_CV->GetEntry(ievent);
    if(!sel_h8) continue;
    double var = est_nu_e_h8->at(ee::kMuonKin);
    if(!is_data) h_CV_Tot->Fill(var,weight_CV);
    h_CV.at(category)->Fill(var,weight_CV);
  }
  f_in_CV->Close();

  std::cout << "Analysing Dirt" << std::endl;
  TFile* f_in_Dirt = nullptr;
  TTree* t_in_Dirt = nullptr;
  LoadTreeFiltered(in_dir+file_Dirt,f_in_Dirt,t_in_Dirt,false,false);
  for(int ievent=0;ievent<t_in_Dirt->GetEntries();ievent++){
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_Dirt->GetEntries() << std::endl;
    t_in_Dirt->GetEntry(ievent);
    if(!sel_h8) continue;
    double var = est_nu_e_h8->at(ee::kMuonKin);
    if(!is_data) h_CV_Tot->Fill(var,weight_Dirt);
    h_CV.at(category)->Fill(var,weight_Dirt);
  }
  f_in_Dirt->Close();

  std::cout << "Analysing EXT" << std::endl;
  TFile* f_in_EXT = nullptr;
  TTree* t_in_EXT = nullptr;
  LoadTreeFiltered(in_dir+file_EXT,f_in_EXT,t_in_EXT,true,false);
  for(int ievent=0;ievent<t_in_EXT->GetEntries();ievent++){
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_EXT->GetEntries() << std::endl;
    t_in_EXT->GetEntry(ievent);
    if(!sel_h8) continue;
    double var = est_nu_e_h8->at(ee::kMuonKin);
    if(!is_data) h_CV_Tot->Fill(var,weight_EXT);
    h_CV.at(category)->Fill(var,weight_EXT);
  }
  f_in_EXT->Close();


  std::cout << "Analysing Variations" << std::endl;
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    TFile* f_in_Vars = nullptr;
    TTree* t_in_Vars = nullptr;
    LoadTreeFiltered(in_dir+files_Vars.at(i_s),f_in_Vars,t_in_Vars,false,false);
    for(int ievent=0;ievent<t_in_Vars->GetEntries();ievent++){
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_Vars->GetEntries() << std::endl;
      t_in_Vars->GetEntry(ievent);
      if(!sel_h8) continue;
      double var = est_nu_e_h8->at(ee::kMuonKin);
      h_Vars_Tot.at(i_s)->Fill(var,weight_Vars.at(i_s));
      h_Vars.at(category).at(i_s)->Fill(var,weight_Vars.at(i_s));
    }
   f_in_Vars->Close();
  }

  
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    h_Vars_Tot.at(i_s)->Add(h_CV.at(kEXT));
    h_Vars_Tot.at(i_s)->Add(h_CV.at(kDirt));
    h_Vars.at(kEXT).at(i_s)->Add(h_CV.at(kEXT));
    h_Vars.at(kDirt).at(i_s)->Add(h_CV.at(kDirt));
  }
  
  gSystem->Exec(("mkdir -p Analysis/"+label+"/rootfiles/").c_str());
  TFile* f_out = TFile::Open(("Analysis/"+label+"/rootfiles/Detvars.root").c_str(),"RECREATE");

  DivideByBinWidth(h_CV_Tot);
  h_CV_Tot->Write();
 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    DivideByBinWidth(h_CV.at(i_c));
    h_CV.at(i_c)->Write();
  }

  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    DivideByBinWidth(h_Vars_Tot.at(i_s));
    h_Vars_Tot.at(i_s)->Write();      
  }  

  for(size_t i_c=0;i_c<categories.size();i_c++){
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      DivideByBinWidth(h_Vars.at(i_c).at(i_s));
      h_Vars.at(i_c).at(i_s)->Write();  
    }
  }      

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
