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

using namespace syst;

void MakeCovDetvar(){

  const double scale = 1.0;

  // Label and set the branches defining the selection and systematics
  std::string label = "MuonMom";
  hist::DetvarHistogramManager h(label,true);
  h.LoadTemplate();

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files = {
    "run4_detvar/Filtered_Merged_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_cv_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_lya_surprise_reco2_hist_4d.root",
    "run4_detvar/Filtered_Merged_DetVar_Run45_v10_04_07_19_BNB_nu_overlay_sce_surprise_reco2_hist_4d.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_dirt_overlay_surprise_reco2_hist_4d.root",
    "Filtered_Merged_MCC9.10_Run4acd5_v10_04_07_14_BNB_beam_off_surprise_reco2_hist_4d.root"
  };

 
  for(int i_f=0;i_f<files.size();i_f++){
    std::string file = files.at(i_f);
    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(in_dir+file,f_in,t_in,is_overlay,load_syst);
    for(int ievent=0;ievent<t_in->GetEntries();ievent++){
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

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

      const double& t = vars_t.at("MuonMom");
      const double& r = vars_h8.at("MuonMom");
      h.FillHistograms2D(is_signal_t,sel_h8,t,r);

    }
    f_in->Close();
  }

  //h.Write();

}
