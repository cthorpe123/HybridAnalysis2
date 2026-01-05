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

using namespace syst;

// Tune binning so all bins have data FE of this value

void MakeBinning(){

  TH1D* h_TrueMuonMom = new TH1D("h_MuonMom",";True Muon Momentum (GeV);Events/GeV",1000,0.0,2.0);
  TH1D* h_MuonMom = new TH1D("h_MuonMom",";Reco Muon Momentum (GeV);Events/GeV",1000,0.0,2.0);

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";

  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

      if(is_signal_t)
        h_TrueMuonMom->Fill(muon_mom_t->Mag(),POT_weight);

      if(is_signal_t && sel_h8)
        h_MuonMom->Fill(muon_mom_h8->Mag(),POT_weight);

    }

  }

  MakeBinningTemplate("MuonMom",h_TrueMuonMom,true); 
  MakeBinningTemplate("MuonMom",h_MuonMom); 

}
