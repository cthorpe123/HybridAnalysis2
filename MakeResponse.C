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
#include "Response.h"

void MakeResponse(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4c.root",
    "Filtered_Merged_MCC9.10_Run4a4c4d5_v10_04_07_13_BNB_nu_overlay_surprise_reco2_hist_4d.root"
  };

  hist::Response r_W("W");
  r_W.SetTemplate(";W (GeV)",50,0.9,5.0);
  r_W.AddSpecialUniv("ExtraGamma");
  r_W.AddSpecialUniv("ExtraProt");
  r_W.AddSpecialUniv("ExtraPi");
  r_W.DiscardEmpty();

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

      double w_ExtraGamma = npi0_t > 0 ? 5.0 : 1.0;
      double w_ExtraProt = nprot_t > 1 ? 5.0 : 1.0;
      double w_ExtraPi = npi_t > 0 ? 5.0 : 1.0;

      r_W.Fill(W_t,W_h8,is_signal_t,sel_h8);
      r_W.FillSpecialUniv("ExtraGamma",W_t,W_h8,is_signal_t,sel_h8,w_ExtraGamma);
      r_W.FillSpecialUniv("ExtraProt",W_t,W_h8,is_signal_t,sel_h8,w_ExtraProt);
      r_W.FillSpecialUniv("ExtraPi",W_t,W_h8,is_signal_t,sel_h8,w_ExtraPi);

    }

  }

  r_W.Write();

}

