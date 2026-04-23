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

void MakeCov(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root",
    "run4c/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root",
    "run4c/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
    "run4c/Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root",
    "run4d/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root",
    "run5/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root"
  };

  std::vector<std::string> channels_t = {"1p0pi","2p0pi","1p1pi","2p1pi"};
  std::vector<std::string> channels_r = {"All"};

  std::vector<std::string> vars = {"MuonMom","MuonCosTheta","NProt","NPi","NSh","ProtonKE","PionE","PiZeroE","W"};
  for(int i_e=0;i_e<ee::kMAX;i_e++)
    vars.push_back(ee::estimators_str.at(i_e));

  std::map<std::string,hist::MultiChannelHistogramManager> h_m;
  for(std::string var : vars){
    h_m.emplace(var,hist::MultiChannelHistogramManager(var,true));
    h_m.at(var).SetTrueChannelList(channels_t);
    h_m.at(var).SetRecoChannelList(channels_r);
    h_m.at(var).KeepAll();
    h_m.at(var).LoadTemplates();
    h_m.at(var).MakeHM();
  } 

  const int spline_pts = 10;

  // Function pointers for the various reweighters 
  std::map<std::string,std::vector<double>(*)()> r_m;

  // Declare the function pointers here
  auto f_ExtraPi = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(npi_t > 0 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraPi",f_ExtraPi); 

  auto f_Extra2Pi = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(npi_t > 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra2Pi",f_Extra2Pi); 
  
  auto f_ExtraP = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t > 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraP",f_ExtraP); 

  auto f_Extra1P = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t == 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra1P",f_Extra1P); 

  auto f_Extra2P = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t == 2 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra2P",f_Extra2P); 

  auto f_Extra3P = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t == 3 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra3P",f_Extra3P); 

  auto f_ExtraNP = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nprot_t > 3 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraNP",f_ExtraNP); 

  auto f_ExtraG = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nsh_t > 0 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("ExtraG",f_ExtraG); 

  auto f_Extra2G = [](){ std::vector<double> x; for(int i=0;i<spline_pts;i++) x.push_back(nsh_t > 1 ? i*5.0/spline_pts : 1.0); return x; };
  r_m.emplace("Extra2G",f_Extra2G); 

  for(std::string var : vars){
    std::cout << var << std::endl;
    for(const auto &item : r_m){ 
      std::cout << item.first << std::endl;
      for(int i=0;i<spline_pts;i++){
        h_m.at(var).AddSpecialUniv(item.first+"_"+std::to_string(i));
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
      
      std::string ch_t = "All";        
      if(nprot_t == 1 && npi_t == 0) ch_t = "1p0pi";
      if(nprot_t == 2 && npi_t == 0) ch_t = "2p0pi";
      if(nprot_t == 1 && npi_t == 1) ch_t = "1p1pi";
      if(nprot_t == 2 && npi_t == 1) ch_t = "2p1pi";

      std::string ch_h8 = "All";

      std::map<std::string,double> vars_t = {
        {"MuonMom",muon_mom_t->Mag()},
        {"MuonCosTheta",muon_mom_t->CosTheta()},
        {"NProt",nprot_t},
        {"NPi",npi_t},
        {"NSh",nsh_t},
        {"ProtonKE",proton_p4_t->E()-nprot_t*Mp},
        {"PionE",pion_p4_t->E()},
        {"PiZeroE",gamma_p4_t->E()}, 
        {"W",W_t}
      };

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        vars_t[ee::estimators_str.at(i_e)] = est_nu_e_t->at(i_e);

      std::map<std::string,double> vars_h8 = {
        {"MuonMom",muon_mom_h8->Mag()},
        {"MuonCosTheta",muon_mom_h8->CosTheta()},
        {"NProt",nprot_h8},
        {"NPi",npi_h8},
        {"NSh",nsh_h8},
        {"ProtonKE",proton_p4_h8->E()-nprot_h8*Mp},
        {"PionE",pion_p4_h8->E()},
        {"PiZeroE",gamma_p4_h8->E()},
        {"W",W_h8}
      };

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        vars_h8[ee::estimators_str.at(i_e)] = est_nu_e_h8->at(i_e);

      std::map<std::string,std::vector<double>> w_m; 
      for(const auto &item : r_m) w_m[item.first] = item.second(); 

      for(const auto &item : h_m){
        std::string var = item.first;
        if(vars_t.find(var) == vars_t.end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
        if(vars_h8.find(var) == vars_h8.end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
        const double& t = vars_t.at(var);
        const double& r = vars_h8.at(var);
        h_m.at(var).FillHistograms2D(is_signal_t,sel_h8,t,r,load_syst,ch_t,ch_h8);
        for(const auto &w : w_m){ 
          for(int i=0;i<spline_pts;i++)
            h_m.at(var).FillSpecialHistograms2D(w.first+"_"+std::to_string(i),is_signal_t,sel_h8,t,r,w.second.at(i),ch_t,ch_h8);
        }
      }

    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    h_m.at(var).Write();
  }

}

