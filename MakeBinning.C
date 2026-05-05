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
#include "WeightFuncs.h"

using namespace syst;
using namespace binning;

// Tune binning so all bins have data FE of this value

void MakeBinning(){

  std::vector<std::string> channels_t = {"All"};
  std::vector<std::string> channels_r = {"All"};

  std::string vars;
  std::map<std::string,TH1D*> h_m;
  std::map<std::string,std::map<std::string,TH1D*>> h_reco_m,h_true_m;

  h_m["2pOpeningAngle"] = new TH1D("h_2pOpeningAngle",";2p Opening Angle;Events",10000,0.0,180);
  h_m["2pAsym"] = new TH1D("h_2pAsym",";2p Asym;Events",10000,0.0,1.0);

  for(const auto &item : h_m){
    h_reco_m[item.first] = std::map<std::string,TH1D*>(); 
    h_true_m[item.first] = std::map<std::string,TH1D*>(); 
    for(std::string ch :channels_r) h_reco_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_Reco_"+ch).c_str());
    for(std::string ch :channels_t) h_true_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_True_"+ch).c_str());
  }

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root"
  };

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

      //std::cout << "Event " << ievent << std::endl;

      std::string channel_t = "All";
      std::string channel_h8 = "All";

      std::map<std::string,double> vars_t = {
        {"MuonMom",muon_mom_t->Mag()},
        {"MuonCosTheta",muon_mom_t->CosTheta()},
        {"LeadProtonKE",-1},
        {"2pOpeningAngle",-1},
        {"2pAsym",-1},
        {"NProt",nprot_t},
        {"NPi",npi_t},
        {"NSh",nsh_t},
        {"ProtonKE",proton_p4_t->E()-nprot_t*Mp},
        {"PionE",pion_p4_t->E()},
        {"PiZeroE",gamma_p4_t->E()}, 
        {"W",W_t},
        {"Channel",ch_t}
      };

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        vars_t[ee::estimators_str.at(i_e)] = est_nu_e_t->at(i_e);

      std::map<std::string,double> vars_h8 = {
        {"MuonMom",muon_mom_h8->Mag()},
        {"MuonCosTheta",muon_mom_h8->CosTheta()},
        {"LeadProtonKE",-1},
        {"2pOpeningAngle",-1},
        {"2pAsym",-1},
        {"NProt",nprot_h8},
        {"NPi",npi_h8},
        {"NSh",nsh_h8},
        {"ProtonKE",proton_p4_h8->E()-nprot_h8*Mp},
        {"PionE",pion_p4_h8->E()},
        {"PiZeroE",gamma_p4_h8->E()},
        {"W",W_h8},
        {"Channel",ch_h8}
      };

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        vars_h8[ee::estimators_str.at(i_e)] = est_nu_e_h8->at(i_e);

      if(is_signal_t){
        vars_t.at("LeadProtonKE") = protons_t->at(0).E() - Mp;
        if(nprot_t == 2){
          vars_t.at("2pOpeningAngle") = 180/3.142*protons_t->at(0).Vect().Angle(protons_t->at(1).Vect());
          vars_t.at("2pAsym") = weight::Asymmetry3({protons_t->at(0)},{protons_t->at(1)});
          }
        }
      

      if(sel_h8){
        vars_h8.at("LeadProtonKE") = protons_h8->at(0).E() - Mp;
        if(nprot_h8 == 2){
          vars_h8.at("2pOpeningAngle") = 180/3.142*protons_h8->at(0).Vect().Angle(protons_h8->at(1).Vect());
          vars_h8.at("2pAsym") = weight::Asymmetry3({protons_h8->at(0)},{protons_h8->at(1)});
        }
      }

      if(is_signal_t){
        for(const auto &item : h_true_m){
          std::string var = item.first;
          if(vars_t.find(var) == vars_t.end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
          if(in_vec(channels_t,channel_t))
           h_true_m.at(var).at(channel_t)->Fill(vars_t.at(var),POT_weight);
        }
      }

      if(is_signal_t && sel_h8){
        for(const auto &item : h_reco_m){
          std::string var = item.first;
          if(vars_h8.find(var) == vars_h8.end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
          if(in_vec(channels_r,channel_h8))
            h_reco_m.at(var).at(channel_h8)->Fill(vars_h8.at(var),POT_weight);
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