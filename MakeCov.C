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
#include "Histograms2.h"
//#include "Response.h"

void MakeCov(){

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  hist::HistogramManager h("MuonMom",true);
  h.LoadTemplate();
  h.KeepAll();
  h.DBBW();
  //h.ShapeOnly();
  h.AddSpecialUniv("ExtraPi");

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      if(ievent > 100000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

      double alt_weight = npi_t > 0 && is_signal_t ? 10.0 : 1.0;

      //h.FillTruthHistograms(is_signal_t,muon_mom_t->Mag(),load_syst);
      //h.FillRecoHistograms(sel_h8,muon_mom_h8->Mag(),load_syst);
      h.FillHistograms2D(is_signal_t,sel_h8,muon_mom_t->Mag(),muon_mom_h8->Mag(),load_syst);

      //h.FillSpecialTruthHistograms("ExtraPi",is_signal_t,muon_mom_t->Mag(),alt_weight);
      //h.FillSpecialRecoHistograms("ExtraPi",sel_h8,muon_mom_h8->Mag(),alt_weight);
      h.FillSpecialHistograms2D("ExtraPi",is_signal_t,sel_h8,muon_mom_t->Mag(),muon_mom_h8->Mag(),alt_weight);
 
    }

  }

  h.Write();

}

