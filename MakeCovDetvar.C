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
#include "Histograms.h"

using namespace syst;

void MakeCovDetvar(){

  const double scale = 1.0;

  // Label and set the branches defining the selection and systematics
  std::string label = "RecoE_0p_Test";
  hist::DetvarHistogramManager h(label);
  h.DBBW(); 

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/detvar/";
  std::vector<std::string> files = {
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_cv_surprise_reco2_hist.root",
    "../Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "../Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lya_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyd_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_lyr_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_SCE_surprise_reco2_hist.root",
    "Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_recomb2_surprise_reco2_hist.root" 
      //"Filtered_Merged_DetVar_Run45_v10_04_07_15_BNB_nu_overlay_WMX_surprise_TEST_reco2_hist.root",
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
      bool sel = in_tpc_pd && has_muon_pd && in_tpc_lt && !nprot_lt && !npi_lt && !nsh_lt;
      if(!sel) continue;
      double var = muon_mom_h8->Mag();
      double var2 = 0.5;
      h.FillHistograms(var);

    }
    f_in->Close();
  }

  h.Write();

}
