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

  hist::HistogramManager h_MuonMom("MuonMom",true);
  h_MuonMom.LoadTemplate();
  //h_MuonMom.KeepAll();
  h_MuonMom.DBBW();

  hist::HistogramManager h_MuonCosTheta("MuonCosTheta",true);
  h_MuonCosTheta.LoadTemplate();
  //h_MuonCosTheta.KeepAll();
  h_MuonCosTheta.DBBW();

  hist::HistogramManager h_ProtonE("ProtonE",true);
  h_ProtonE.LoadTemplate();
  //h_ProtonE.KeepAll();
  h_ProtonE.DBBW();

  hist::HistogramManager h_PionE("PionE",true);
  h_PionE.LoadTemplate();
  //h_PionE.KeepAll();
  h_PionE.DBBW();

  hist::HistogramManager h_ShowerE("ShowerE",true);
  h_ShowerE.LoadTemplate();
  //h_ShowerE.KeepAll();
  h_ShowerE.DBBW();

  hist::HistogramManager h_W("W",true);
  h_W.LoadTemplate();
  //h_W.KeepAll();
  h_W.DBBW();

  std::vector<hist::HistogramManager> h_NuE;
  for(int i_e=0;i_e<ee::kMAX;i_e++){
    h_NuE.push_back(hist::HistogramManager(ee::estimators_str.at(i_e),true));
    h_NuE.back().LoadTemplate();
    h_NuE.back().DBBW();
    //h_NuE.back().KeepAll();
  }

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

      h_MuonMom.FillHistograms2D(is_signal_t,sel_h8,muon_mom_t->Mag(),muon_mom_h8->Mag(),load_syst);
      h_MuonCosTheta.FillHistograms2D(is_signal_t,sel_h8,muon_mom_t->CosTheta(),muon_mom_h8->CosTheta(),load_syst);
      h_ProtonE.FillHistograms2D(is_signal_t,sel_h8,proton_p4_t->E()-nprot_t*Mp,proton_p4_h8->E()-nprot_h8*Mp,load_syst);
      h_PionE.FillHistograms2D(is_signal_t,sel_h8,pion_p4_t->E(),pion_p4_h8->E(),load_syst);
      h_ShowerE.FillHistograms2D(is_signal_t,sel_h8,gamma_p4_t->E(),gamma_p4_h8->E(),load_syst);
      h_W.FillHistograms2D(is_signal_t,sel_h8,W_t,W_h8,load_syst);

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        h_NuE.at(i_e).FillHistograms2D(is_signal_t,sel_h8,est_nu_e_t->at(i_e),est_nu_e_h8->at(i_e),load_syst);
 
    }

  }

  h_MuonMom.Write();
  h_MuonCosTheta.Write();
  h_ProtonE.Write();
  h_PionE.Write();
  h_ShowerE.Write();
  h_W.Write();
  for(int i_e=0;i_e<ee::kMAX;i_e++) h_NuE.at(i_e).Write();

}

