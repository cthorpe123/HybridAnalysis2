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

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/retupled/";
  std::vector<std::string> files_v = {
    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root",
    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_on_metapatch_retuple_retuple_hist.root",
    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_off_metapatch_retuple_retuple_hist.root",
    
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4c.root",
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_on_surprise_reco2_hist_4c.root",
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root",

    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_on_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root",
    
    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_on_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root"
  };

  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  //std::vector<std::string> vars = {"MuonMom","MuonCosTheta","LeadProtonKE","ProtonKE"};
  std::vector<std::string> vars = var_names;
  std::vector<std::string> int_vars = {"NProt","NPi","NSh","NPi0"};

  std::map<std::string,hist::MultiChannelHistogramManager> h_m;
  for(std::string var : vars){
    h_m.emplace(var,hist::MultiChannelHistogramManager(var,true));
    h_m.at(var).SetTrueChannelList(channels_t);
    h_m.at(var).SetRecoChannelList(channels_r);
    h_m.at(var).KeepAll();
    if(!in_vec(int_vars,var)) h_m.at(var).LoadTemplates();
    else if(var != "NProt") h_m.at(var).SetTemplates("",3,-0.5,2.5,3,-0.5,2.5); 
    else h_m.at(var).SetTemplates("",3,0.5,3.5,3,0.5,3.5); 
    h_m.at(var).MakeHM();
  } 
  
  h_m.emplace("Norm",hist::MultiChannelHistogramManager("Norm",true));
  h_m.at("Norm").SetTrueChannelList(channels_t);
  h_m.at("Norm").SetRecoChannelList(channels_r);
  h_m.at("Norm").KeepAll();
  h_m.at("Norm").SetTemplates("",1,0,1,1,0,1); 
  h_m.at("Norm").MakeHM();

  h_m.emplace("Enu",hist::MultiChannelHistogramManager("Enu",true));
  h_m.at("Enu").SetTrueChannelList(channels_t);
  h_m.at("Enu").SetRecoChannelList(channels_r);
  h_m.at("Enu").KeepAll();
  h_m.at("Enu").LoadTemplates();
  h_m.at("Enu").MakeHM();
  

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 20000) break;
      if (ievent % 1000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);
      
      std::string channel_t = "All";
      std::string channel_h8 = "All";

      vars_t->emplace("Enu",nu_e);
      vars_h8->emplace("Enu",nu_e);

      vars_t->emplace("Norm",0.5);
      vars_h8->emplace("Norm",0.5);

      // Fill the histograms
      for(const auto &item : h_m){
        std::string var = item.first;
        if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
        if(vars_h8->find(var) == vars_h8->end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
        const double& t = vars_t->at(var);
        const double& r = vars_h8->at(var);
        h_m.at(var).FillHistograms2D(is_signal_t,sel_h8,t,r,load_syst,channel_t,channel_h8);
      }

    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var).Write();
  }

}

