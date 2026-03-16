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
#include "BinningFuncs.h"

using namespace syst;
using namespace binning;

// Tune binning so all bins have data FE of this value

void MakeBinning(){

  std::vector<std::string> channels = {"1p","2p","3p","Other"};

  std::string vars;
  std::map<std::string,TH1D*> h_m;
  
  h_m["MuonMom"] = new TH1D("h_MuonMom",";Reco Muon Momentum (GeV);Events/GeV",10000,0.0,2.0);
  h_m["MuonCosTheta"] = new TH1D("h_MuonCosTheta",";Reco Muon Cos(#theta);Events/Unit",10000,-1.0,1.0);
  h_m["NProt"] = new TH1D("h_NProt",";Reco N Protons;Events",10000,0.5,4.5);
  h_m["ProtonKE"] = new TH1D("h_ProtonKE",";Reco Proton KE (GeV);Events/GeV",10000,0.01,1.0);
  h_m["PionE"] = new TH1D("h_PionE",";Reco Pion E (GeV);Events/GeV",10000,0.01,1.0);
  h_m["PiZeroE"] = new TH1D("h_PiZeroE",";Reco #pi^{0} E (GeV);Events/GeV",10000,0.01,1.0);

  std::map<std::string,std::map<std::string,TH1D*>> h_true_m;
  for(const auto &item : h_m){
    h_true_m[item.first] = std::map<std::string,TH1D*>(); 
    for(std::string ch :channels) h_true_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_"+ch).c_str());
  }

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4d.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
    "Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_5.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
    "Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root"
  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 100000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);
        
      std::string ch_t = "Other";
      if(is_signal_t && nprot_t == 1) ch_t = "1p";
      else if (is_signal_t && nprot_t == 2) ch_t = "2p";

      std::map<std::string,double> vars_t = {
        {"MuonMom",muon_mom_t->Mag()},
        {"MuonCosTheta",muon_mom_t->CosTheta()},
        {"NProt",nprot_t-0.5},
        {"ProtonKE",proton_p4_t->E()-nprot_t*Mp},
        {"PionE",pion_p4_t->E()},
        {"PiZeroE",gamma_p4_t->E()}
      };

      std::map<std::string,double> vars_h8 = {
        {"MuonMom",muon_mom_h8->Mag()},
        {"MuonCosTheta",muon_mom_h8->CosTheta()},
        {"NProt",nprot_h8-0.5},
        {"ProtonKE",proton_p4_h8->E()-nprot_h8*Mp},
        {"PionE",pion_p4_h8->E()},
        {"PiZeroE",gamma_p4_h8->E()}
      };


      if(is_signal_t){
        for(const auto &item : h_true_m){
          std::string var = item.first;
          if(vars_t.find(var) == vars_t.end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
          h_true_m.at(var).at(ch_t)->Fill(vars_t.at(var),POT_weight);
        }
      }

      if(is_signal_t && sel_h8){
        for(const auto &item : h_m){
          std::string var = item.first;
          if(vars_h8.find(var) == vars_h8.end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
          h_m.at(var)->Fill(vars_h8.at(var),POT_weight);
        }
      }

    }

  }


  for(const auto &item : h_m){
    std::string var = item.first;
    MakeMultiChannelTemplate(var,h_true_m.at(var),true);
    MakeBinningTemplate(var,h_m.at(var),false);
  }


}
