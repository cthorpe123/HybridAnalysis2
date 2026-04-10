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

  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  std::string vars;
  std::map<std::string,TH1D*> h_m;
  std::map<std::string,std::map<std::string,TH1D*>> h_reco_m,h_true_m;
  
  h_m["MuonMom"] = new TH1D("h_MuonMom",";Reco Muon Momentum (GeV);Events/GeV",10000,0.0,2.0);
   
  for(const auto &item : h_m){
    h_reco_m[item.first] = std::map<std::string,TH1D*>(); 
    h_true_m[item.first] = std::map<std::string,TH1D*>(); 
    for(std::string ch :channels_r) h_reco_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_Reco_"+ch).c_str());
    for(std::string ch :channels_t) h_true_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_True_"+ch).c_str());
  }

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "run4b/Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root"
  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 200000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);
      
      std::string ch_t = "All";        
      std::string ch_h8 = "All";

      std::map<std::string,double> vars_t = {
        {"MuonMom",muon_mom_t->Mag()}
      };

      std::map<std::string,double> vars_h8 = {
        {"MuonMom",muon_mom_h8->Mag()}
      };

        if(is_signal_t && in_vec(channels_t,ch_t)){
        for(const auto &item : h_true_m){
          std::string var = item.first;
          if(vars_t.find(var) == vars_t.end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
          h_true_m.at(var).at(ch_t)->Fill(vars_t.at(var),POT_weight);
        }
      }

      if(is_signal_t && sel_h8 && in_vec(channels_r,ch_h8)){
        for(const auto &item : h_reco_m){
          std::string var = item.first;
          if(vars_h8.find(var) == vars_h8.end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
          h_reco_m.at(var).at(ch_h8)->Fill(vars_h8.at(var),POT_weight);
        }
      }

    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    MakeMultiChannelTemplate(var,h_reco_m.at(var),false);
    MakeMultiChannelTemplate(var,h_true_m.at(var),true);
  }


}
