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

  TH1D* h_TrueMuonMom = new TH1D("h_TrueMuonMom",";True Muon Momentum (GeV);Events/GeV",1000,0.0,2.0);
  TH1D* h_MuonMom = new TH1D("h_MuonMom",";Reco Muon Momentum (GeV);Events/GeV",1000,0.0,2.0);

  TH1D* h_TrueMuonCosTheta = new TH1D("h_TrueMuonCosTheta",";True Muon Cos(#theta);Events/Unit",1000,-1.0,1.0);
  TH1D* h_MuonCosTheta = new TH1D("h_MuonCosTheta",";Reco Muon Cos(#theta);Events/Unit",1000,-1.0,1.0);

  TH1D* h_TrueProtonE = new TH1D("h_TrueProtonE",";True Proton Energy (GeV);Events/GeV",1000,0.0,1.0);
  TH1D* h_ProtonE = new TH1D("h_ProtonE",";Reco Proton Energy (GeV);Events/GeV",1000,0.0,1.0);

  TH1D* h_TruePionE = new TH1D("h_TruePionE",";True Pion Energy (GeV);Events/GeV",1000,0.05,1.0);
  TH1D* h_PionE = new TH1D("h_PionE",";Reco Pion Energy (GeV);Events/GeV",1000,0.05,1.0);

  TH1D* h_TrueShowerE = new TH1D("h_TrueShowerE",";True Shower Energy (GeV);Events/GeV",1000,0.05,1.0);
  TH1D* h_ShowerE = new TH1D("h_ShowerE",";Reco Shower Energy (GeV);Events/GeV",1000,0.05,1.0);

  TH1D* h_TrueW = new TH1D("h_TrueW",";True W_{vis} (GeV);Events/GeV",1000,Mp+0.05,5.0);
  TH1D* h_W = new TH1D("h_W",";Reco W_{vis} (GeV);Events/GeV",1000,Mp+0.05,5.0);

  std::vector<TH1D*> h_TrueNuE;
  std::vector<TH1D*> h_NuE;
  for(int i_e=0;i_e<ee::kMAX;i_e++){
    h_TrueNuE.push_back(new TH1D(("h_True"+ee::estimators_str.at(i_e)).c_str(),";Est Neutrino Energy (GeV);Events",1000,0.0,2.5));
    h_NuE.push_back(new TH1D(("h_"+ee::estimators_str.at(i_e)).c_str(),";Reco Neutrino Energy (GeV);Events",1000,0.0,2.5));
  }


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

      if(is_signal_t){
        h_TrueMuonMom->Fill(muon_mom_t->Mag(),POT_weight);
        h_TrueMuonCosTheta->Fill(muon_mom_t->CosTheta(),POT_weight);
        h_TrueProtonE->Fill(proton_p4_t->E()-nprot_t*Mp,POT_weight);
        h_TruePionE->Fill(pion_p4_t->E(),POT_weight);
        h_TrueShowerE->Fill(gamma_p4_t->E(),POT_weight);
        h_TrueW->Fill(W_t,POT_weight);
        for(int i_e=0;i_e<ee::kMAX;i_e++) h_TrueNuE.at(i_e)->Fill(est_nu_e_t->at(i_e),POT_weight);
      }

      if(is_signal_t && sel_h8){
        h_MuonMom->Fill(muon_mom_h8->Mag(),POT_weight);
        h_MuonCosTheta->Fill(muon_mom_h8->CosTheta(),POT_weight);
        h_ProtonE->Fill(proton_p4_h8->E()-nprot_h8*Mp,POT_weight);
        h_PionE->Fill(pion_p4_h8->E(),POT_weight);
        h_ShowerE->Fill(gamma_p4_h8->E(),POT_weight);
        h_W->Fill(W_h8,POT_weight);
        for(int i_e=0;i_e<ee::kMAX;i_e++) h_NuE.at(i_e)->Fill(est_nu_e_h8->at(i_e),POT_weight);
      }

    }

  }

  MakeBinningTemplate("MuonMom",h_TrueMuonMom,true); 
  MakeBinningTemplate("MuonMom",h_MuonMom); 
  MakeBinningTemplate("MuonCosTheta",h_TrueMuonCosTheta,true); 
  MakeBinningTemplate("MuonCosTheta",h_MuonCosTheta); 
  MakeBinningTemplate("ProtonE",h_TrueProtonE,true); 
  MakeBinningTemplate("ProtonE",h_ProtonE); 
  MakeBinningTemplate("PionE",h_TruePionE,true); 
  MakeBinningTemplate("PionE",h_PionE); 
  MakeBinningTemplate("ShowerE",h_TrueShowerE,true); 
  MakeBinningTemplate("ShowerE",h_ShowerE); 
  MakeBinningTemplate("W",h_TrueW,true); 
  MakeBinningTemplate("W",h_W); 
  for(int i_e=0;i_e<ee::kMAX;i_e++){
    MakeBinningTemplate(ee::estimators_str.at(i_e),h_TrueNuE.at(i_e),true); 
    MakeBinningTemplate(ee::estimators_str.at(i_e),h_NuE.at(i_e)); 
  }

}
