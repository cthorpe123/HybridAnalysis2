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
#include "DetvarHistograms.h"
#include "MultiChannelHistograms.h"

using namespace syst;

void MakeCovDetvar(){

  //std::vector<std::string> channels_t = {"1p","2p","3p"}; 
  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta"};
  //std::vector<std::string> vars = var_names;
  std::vector<std::string> int_vars = {"NProt","NPi","NSh","NPi0"};

  std::map<std::string,hist::MultiChannelHistogramManager> h_m;
  for(std::string var : vars){
    h_m.emplace(var,hist::MultiChannelHistogramManager(var,true));
    h_m.at(var).DetvarMode();
    h_m.at(var).SetTrueChannelList(channels_t);
    h_m.at(var).SetRecoChannelList(channels_r);
    h_m.at(var).KeepAll();
    if(!in_vec(int_vars,var)) h_m.at(var).LoadTemplates();
    else h_m.at(var).SetTemplates("",4,-0.5,3.5,4,-0.5,3.5); 
    h_m.at(var).MakeHM();
  } 

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/retupled/";
  std::vector<std::string> files = {
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_cv_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lya_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lyd_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lyr_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_recomb2_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_sce_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_WMX_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_WMYZ_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root",

    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_cv_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lya_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lyd_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lyr_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_recomb2_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_sce_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_WMX_surprise_reco2_hist_5.root",
    "run5_detvar/Filtered_Merged_checkout_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_WMYZ_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root" 
  };

  for(int i_f=0;i_f<files.size();i_f++){
    std::string file = files.at(i_f);
    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(in_dir+file,f_in,t_in,is_overlay,load_syst);
    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if (ievent % 1000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);

      //if(run > 25400 && run < 25600) continue;

      std::string channel_t = "All";
      std::string channel_h8 = "All";

      for(const auto &item : h_m){
        std::string var = item.first;
        if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
        if(vars_h8->find(var) == vars_h8->end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
        const double& t = vars_t->at(var);
        const double& r = vars_h8->at(var);
        h_m.at(var).FillHistograms2D(is_signal_t,sel_h8,t,r,load_syst,channel_t,channel_h8);
      }

    }

    f_in->Close();

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var).Write();
  }

}
