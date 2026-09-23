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

  std::map<std::string,TH1D*> h_m;
  std::map<std::string,std::map<std::string,TH1D*>> h_reco_m,h_true_m;

  h_m["Enu"] = new TH1D("h_Enu",";Neutrino Energy [GeV];Events/GeV",10000, 0,3);
  h_m["MuonMom"]                = new TH1D("h_MuonMom",                ";Muon Momentum [GeV];Events/GeV",              10000, 0,   3);
  h_m["MuonCosTheta"]           = new TH1D("h_MuonCosTheta",           ";Muon cos#theta;Events/Unit",                   10000,-1,   1);
  h_m["LeadProtonKE"]           = new TH1D("h_LeadProtonKE",           ";Leading Proton KE [GeV];Events/GeV",          10000, 0,   1);
  h_m["ProtonKE"]               = new TH1D("h_ProtonKE",               ";Total Proton KE [GeV];Events/GeV",            10000, 0,   2);
  h_m["LeadPionE"]              = new TH1D("h_LeadPionE",              ";Leading Pion Energy [GeV];Events/GeV",        10000, 0,   1.5);
  h_m["1p1piOpeningAngle"]      = new TH1D("h_1p1piOpeningAngle",      ";1p1#pi Opening Angle [deg];Events/GeV",       10000, 0, 180);
  h_m["1p1piAsym"]              = new TH1D("h_1p1piAsym",              ";1p1#pi Asymmetry;Events/Unit",                 10000,-1,   1);
  h_m["MuonProtonOpeningAngle"] = new TH1D("h_MuonProtonOpeningAngle", ";Muon-Proton Opening Angle [deg];Events",  10000, 0, 180);
  h_m["2pOpeningAngle"]         = new TH1D("h_2pOpeningAngle",         ";2p Opening Angle [deg];Events",           10000, 0, 180);
  h_m["2pAsym"]                 = new TH1D("h_2pAsym",                 ";2p Asymmetry;Events/Unit",                     10000,-1,   1);
  h_m["2shwOpenAngle"]          = new TH1D("h_2shwOpenAngle",          ";2-Shower Opening Angle [deg];Events/Unit",     10000, 0, 180);
  h_m["2shwAsym"]               = new TH1D("h_2shwAsym",               ";2-Shower Asymmetry;Events/Unit",               10000,-1,   1);
  h_m["PionE"]                  = new TH1D("h_PionE",                  ";Total Pion Energy [GeV];Events/GeV",          10000, 0,   1.5);
  h_m["PiZeroE"]                = new TH1D("h_PiZeroE",                ";Total #pi^{0} Energy [GeV];Events/GeV",       10000, 0,   1.5);
  h_m["W"]                      = new TH1D("h_W",                      ";W [GeV];Events/GeV",                          10000, 1,   4);
  h_m["Channel"]                = new TH1D("h_Channel",                ";Channel;Events",                             27, 0,  27);
  h_m["MuonKin"]                = new TH1D("h_MuonKin",                ";E_{#nu} Muon Kin. [GeV];Events/GeV",          10000, 0,   4);
  h_m["MuonKinWNP"]             = new TH1D("h_MuonKinWNP",             ";E_{#nu} Muon Kin. + NP [GeV];Events/GeV",     10000, 0,   4);
  h_m["PeLEELike0Pi"]           = new TH1D("h_PeLEELike0Pi",           ";E_{#nu} PeLEE-Like 0#pi [GeV];Events/GeV",    10000, 0,   4);
  h_m["TotalEDep"]              = new TH1D("h_TotalEDep",              ";E_{#nu} Total E Dep. [GeV];Events/GeV",       10000, 0,   4);
  h_m["SFMethod"]               = new TH1D("h_SFMethod",               ";E_{#nu} SF Method [GeV];Events/GeV",          10000, 0,   4);
 
  for(const auto &item : h_m){
    h_reco_m[item.first] = std::map<std::string,TH1D*>(); 
    h_true_m[item.first] = std::map<std::string,TH1D*>(); 
    for(std::string ch :channels_r) h_reco_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_Reco_"+ch).c_str());
    for(std::string ch :channels_t) h_true_m.at(item.first)[ch] = (TH1D*)h_m.at(item.first)->Clone((item.first+"_True_"+ch).c_str());
  }

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/retupled/";
  std::vector<std::string> files_v = {

    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root",
    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_beam_off_metapatch_retuple_retuple_hist.root",

    "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4c.root",
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4c.root",
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4c.root",

    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root",

    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_5.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_5.root"

  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if (ievent % 5000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);

      std::string channel_t = "All";
      std::string channel_h8 = "All";

      vars_t->emplace("Enu",nu_e);
      vars_h8->emplace("Enu",nu_e);

      if(std::isnan(weightSplineTimesTune) || std::isinf(weightSplineTimesTune)) continue;

      if(is_signal_t){
        for(const auto &item : h_true_m){
          std::string var = item.first;
          if(vars_t->find(var) == vars_t->end()) throw std::invalid_argument("Variable " + var + " missing from true var map");
          if(in_vec(channels_t,channel_t))
            h_true_m.at(var).at(channel_t)->Fill(vars_t->at(var),POT_weight*weightSplineTimesTune);
        }
      }

      if(is_signal_t && sel_h8){
        for(const auto &item : h_reco_m){
          std::string var = item.first;
          if(vars_h8->find(var) == vars_h8->end()) throw std::invalid_argument("Variable " + var + " missing from reco var map");
          if(in_vec(channels_r,channel_h8))
            h_reco_m.at(var).at(channel_h8)->Fill(vars_h8->at(var),POT_weight*weightSplineTimesTune);
        }
      }

    }

  }

  for(const auto &item : h_m){
    std::string var = item.first;
    MakeMultiChannelTemplate(var,h_reco_m.at(var),false,0.03);
    MakeMultiChannelTemplate(var,h_true_m.at(var),true,0.03);
  }

}
