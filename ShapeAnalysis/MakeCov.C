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
#include "Histograms.h"

void MakeCov(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  for(int i_e=0;i_e<ee::kMAX;i_e++){

    std::string label = "RecoE_"+ee::estimators_str.at(i_e);
    hist::HistogramManager h(label);
    h.LoadTemplate();
    h.DBBW();

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

        if(!sel_h8) continue;
        double var = est_nu_e_h8->at(i_e);

        h.FillHistograms(var,load_syst);

      }

    }

    h.Write();

  }

}
