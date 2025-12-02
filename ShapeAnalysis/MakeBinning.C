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

  TH1D* h_MuonMom = new TH1D("h_MuonMom",";Muon Momentum (GeV);Events",1000,0.0,2.5);
  TH1D* h_MuonCosTheta = new TH1D("h_MuonCosTheta",";Muon Cos(#theta);Events",1000,-1,1);
  TH1D* h_ProtonE = new TH1D("h_ProtonE",";Proton Energy (GeV);Events",1000,0.0,2.5);
  TH1D* h_PionE = new TH1D("h_PionE",";Charged Pion Energy (GeV);Events",1000,0.0,2.5);
  TH1D* h_PiZeroE = new TH1D("h_PiZeroE",";#pi^{0} Energy (GeV);Events",1000,0.0,2.5);
  TH1D* h_W = new TH1D("h_W",";W (GeV);Events",1000,0.0,5.0);
  
  std::vector<TH1D*> h_NuE_v; 
  for(int i_e=0;i_e<ee::kMAX;i_e++)
    h_NuE_v.push_back(new TH1D(("h_"+ee::estimators_str.at(i_e)).c_str(),";#pi^{0} Energy (GeV);Events",1000,0.0,2.5));

  const double target_fe = 0.1;
  const double scale = 1.0;

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

      if(!sel_h8) continue;

      h_MuonMom->Fill(muon_mom_h8->Mag(),POT_weight);
      h_MuonCosTheta->Fill(muon_mom_h8->CosTheta(),POT_weight);
      h_ProtonE->Fill(proton_p4_h8->E() - nprot_h8*Mp,POT_weight);
      if(npi_h8 > 0) h_PionE->Fill(pion_p4_h8->E(),POT_weight);
      if(nsh_h8 > 0) h_PiZeroE->Fill(shower_p4_h8->E(),POT_weight);
      h_W->Fill(W_h8,POT_weight);

      for(int i_e=0;i_e<ee::kMAX;i_e++)
        h_NuE_v.at(i_e)->Fill(est_nu_e_h8->at(i_e),POT_weight);

    }

  }

  MakeBinningTemplate("MuonMom",h_MuonMom); 
  MakeBinningTemplate("MuonCosTheta",h_MuonCosTheta); 
  MakeBinningTemplate("ProtonE",h_ProtonE); 
  MakeBinningTemplate("PionE",h_PionE); 
  MakeBinningTemplate("PiZeroE",h_PiZeroE); 
  MakeBinningTemplate("W",h_W);
  for(int i_e=0;i_e<ee::kMAX;i_e++) MakeBinningTemplate(ee::estimators_str.at(i_e),h_NuE_v.at(i_e));


}
