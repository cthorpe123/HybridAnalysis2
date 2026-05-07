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

void MakeCov(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root"
  };

  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  std::vector<std::string> vars = {"2pAsym","2shwOpenAngle","2shwAsym"};

  std::map<std::string,hist::MultiChannelHistogramManager> h_m;
  for(std::string var : vars){
    h_m.emplace(var,hist::MultiChannelHistogramManager(var,true));
    h_m.at(var).SetTrueChannelList(channels_t);
    h_m.at(var).SetRecoChannelList(channels_r);
    h_m.at(var).KeepAll();
    h_m.at(var).LoadTemplates();
    h_m.at(var).MakeHM();
  } 

  for(std::string var : vars){
    for(const std::string wf_label : weight::weight_func_labels){ 
      for(int i=0;i<weight::spline_pts;i++){
        h_m.at(var).AddSpecialUniv(wf_label+"_"+std::to_string(i));
      }
    }
  }

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);
      
      std::string channel_t = "All";
      std::string channel_h8 = "All";

      for(const auto &item : h_m){
        std::string var = item.first;
        if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
        if(vars_h8->find(var) == vars_h8->end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
        const double& t = vars_t->at(var);
        const double& r = vars_h8->at(var);
        h_m.at(var).FillHistograms2D(is_signal_t,sel_h8,t,r,load_syst,channel_t,channel_h8);
        for(const auto &w : weight::r_m){ 
          for(int i=0;i<weight::spline_pts;i++)
            h_m.at(var).FillSpecialHistograms2D(w.first+"_"+std::to_string(i),is_signal_t,sel_h8,t,r,weight_funcs_m->at(w.first).at(i),channel_t,channel_h8);
        }
      }

    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var).Write();
  }

}

