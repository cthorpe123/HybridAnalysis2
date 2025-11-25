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

  // Label and set the branches defining the selection and systematics
  std::string label = "RecoE_Test";
  std::string axis_title = ";Reco E (GeV);Events";

  TFile* f_tp = TFile::Open(("Analysis/"+label+"/rootfiles/BinningTemplate.root").c_str());
  TH1D* h_tp = (TH1D*)f_tp->Get("h_template");
  h_tp->SetDirectory(0);
  f_tp->Close();

  h_tp->GetXaxis()->SetTitle("Reco E (GeV)");
  h_tp->GetYaxis()->SetTitle("Events/GeV");

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";

  std::vector<std::string> files_CV = {
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_cv_surprise_reco2_hist.root",
    "../Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "../Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  std::vector<std::string> files_Vars = {
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lya_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyd_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyr_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_SCE_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_recomb2_surprise_reco2_hist.root" 
      //"Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_WMX_surprise_TEST_reco2_hist.root",
  };

  //TH1D* h_CV_Tot =  new TH1D("h_Detvar_CV_Tot",axis_title.c_str(),bins.size()-1,&bins[0]);
  TH1D* h_CV_Tot = (TH1D*)h_tp->Clone("h_Detvar_CV_Tot");
  std::vector<TH1D*> h_CV;

  std::vector<TH1D*> h_Vars_Tot;
  std::vector<std::vector<TH1D*>> h_Vars;
  
  for(int i_s=0;i_s<kDetvarMAX;i_s++)
    //h_Vars_Tot.push_back(new TH1D(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
    h_Vars_Tot.push_back((TH1D*)h_tp->Clone(("h_Detvar_Vars_Tot_"+detvar_str.at(i_s)).c_str()));
 
  for(size_t i_c=0;i_c<categories.size();i_c++){
    //h_CV.push_back(new TH1D(("h_Detvar_CV_"+categories.at(i_c)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
    h_CV.push_back((TH1D*)h_tp->Clone(("h_Detvar_CV_"+categories.at(i_c)).c_str()));
    h_Vars.push_back(std::vector<TH1D*>());
    for(int i_s=0;i_s<kDetvarMAX;i_s++){
      //h_Vars.back().push_back(new TH1D(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)).c_str(),axis_title.c_str(),bins.size()-1,&bins[0]));
      h_Vars.back().push_back((TH1D*)h_tp->Clone(("h_Detvar_Vars_"+categories.at(i_c)+"_"+detvar_str.at(i_s)).c_str()));
    }
  }

  
  for(int i_f=0;i_f<files_CV.size();i_f++){
    std::string file = files_CV.at(i_f);
    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(in_dir+file,f_in,t_in,is_overlay,load_syst);
    for(int ievent=0;ievent<t_in->GetEntries();ievent++){
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);
      if(!sel_h8) continue;
      double var = est_nu_e_h8->at(ee::kMuonKin);
      if(!is_data) h_CV_Tot->Fill(var,POT_weight);
      h_CV.at(category)->Fill(var,POT_weight);
    }
    f_in->Close();
  }


  std::cout << "Analysing Variations" << std::endl;
  for(int i_s=0;i_s<kDetvarMAX;i_s++){
    TFile* f_in_Vars = nullptr;
    TTree* t_in_Vars = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(in_dir+files_Vars.at(i_s),f_in_Vars,t_in_Vars,is_overlay,load_syst);
    for(int ievent=0;ievent<t_in_Vars->GetEntries();ievent++){
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in_Vars->GetEntries() << std::endl;
      t_in_Vars->GetEntry(ievent);
      if(!sel_h8) continue;
      double var = est_nu_e_h8->at(ee::kMuonKin);
      h_Vars_Tot.at(detvar_univ)->Fill(var,POT_weight);
      h_Vars.at(category).at(detvar_univ)->Fill(var,POT_weight);
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
