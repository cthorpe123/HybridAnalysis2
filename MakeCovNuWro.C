#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "Systematics.h"
#include "EnergyEstimatorFuncs.h"
#include "MultiChannelHistograms.h"
#include "WeightFuncs.h"

void MakeCovNuWro(){

  bool add_nuwro_fd = true;

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run45_v10_04_07_23_BNB_nuwro_overlay_surprise_reco2_hist_4c.root"
  };

  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  std::vector<std::string> vars = {"MuonMom"};
  //std::vector<std::string> vars = var_names;
  std::vector<std::string> int_vars = {"NProt","NPi","NSh","NPi0"};

  std::map<std::string,hist::MultiChannelHistogramManager> h_m;
  for(std::string var : vars){
    h_m.emplace(var,hist::MultiChannelHistogramManager(var,true));
    h_m.at(var).SetTrueChannelList(channels_t);
    h_m.at(var).SetRecoChannelList(channels_r);
    h_m.at(var).KeepAll();
    if(!in_vec(int_vars,var)) h_m.at(var).LoadTemplates();
    else h_m.at(var).SetTemplates("",4,-0.5,3.5,4,-0.5,3.5); 
    h_m.at(var).MakeHM();
  } 

  weight::SetWeightFuncs();
  for(std::string var : vars){
    for(const auto& wf_label : weight::r_m){ 
      for(int i=0;i<weight::spline_pts;i++){
        h_m.at(var).AddSpecialUniv(wf_label.first+"_"+std::to_string(i));
      }
    }
  }

  if(add_nuwro_fd){
    for(std::string var : vars)
      h_m.at(var).AddSpecialUniv("NuWro_0");
  }

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);
      
      std::string channel_t = "All";
      std::string channel_h8 = "All";

      // Set the spec universe weights if being used
      std::map<std::string,std::vector<double>> weight_m;
      for(const auto &w : weight::r_m){ 
        weight_m[w.first] = w.second();
      } 

      if(add_nuwro_fd) weight_m["NuWro"] = {1.0};
      
      // Fill the histograms
      for(const auto &item : h_m){
        std::string var = item.first;
        if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
        if(vars_h8->find(var) == vars_h8->end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
        const double& t = vars_t->at(var);
        const double& r = vars_h8->at(var);
        h_m.at(var).FillHistograms2D(is_signal_t,sel_h8,t,r,load_syst,channel_t,channel_h8);
        for(const auto &w : weight_m){ 
          for(int i=0;i<w.second.size();i++){
            h_m.at(var).FillSpecialHistograms2D(w.first+"_"+std::to_string(i),is_signal_t,sel_h8,t,r,w.second.at(i),channel_t,channel_h8);
          }
        }
          
      }

    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var).Write();
  }

}

