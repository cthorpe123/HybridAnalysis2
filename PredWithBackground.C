#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"
#include "EnergyEstimatorFuncs.h"

// Hybrid 7 - muon kinematics from pandora and hadrons from LT
// Hybrid 8 - muon kinematics from pandora with MCS for uncontained muons, and hadrons from LT

const std::vector<std::string> methods_str = {"PD","WC","LT","H7","H8"};
enum methds_e {kpd,kwc,klt,h7,h8,kMethMAX};

// Try using different combinations of cuts and frameworks to calculate W

void PredWithBackground(){

  const double scale = 1.0;

  const double data_POT = 1.332E+20;
  const double data_Trig = 31582916.0;

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";

  std::vector<std::string> files_v = {
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_BNB_dirt_surpise_reco2_hist.root",
    "Filtered_Merged_MCC9.10_Run4b_v10_04_07_09_Run4b_BNB_beam_off_surprise_reco2_hist.root"
  };

  std::vector<double> weights = {
    scale*data_POT/7.88166e+20,
    scale*data_POT/3.06E+20,
    scale*data_Trig/88445969.0  
  };

  std::vector<bool> is_data_v = {
    false,
    false,
    true
  };

  std::vector<std::string> categories = {"Signal","BG","Nue","Dirt","EXT"};

  //std::vector<std::vector<std::vector<TH1D*>>> h_Selected_RecoE_v_v_v;
  TH1D* h_Selected_RecoE_v_v_v[ee::kMAX][kMethMAX][categories.size()];


  for(int i_e=0;i_e<ee::kMAX;i_e++){
    for(int i_m=0;i_m<kMethMAX;i_m++){
      for(size_t i_c=0;i_c<categories.size();i_c++){
        std::string label = ee::estimators_str.at(i_e) + "_" + methods_str.at(i_m) + "_" + categories.at(i_c);
        h_Selected_RecoE_v_v_v[i_e][i_m][i_c] = new TH1D(("h_Selected_RecoE_v_v_v"+label).c_str(),";Reconstructed Energy (GeV);Events",50,0.0,3.0);
      }
    }
  }


  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    LoadTreeFiltered(file,f_in,t_in,is_data_v.at(i_f));

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
      t_in->GetEntry(ievent);

      //if(!is_signal_t) continue;
      //if(!muon_contained_t) continue;

      bool sel_h7 = has_muon_pd && sel_lt; 
      TLorentzVector plepton_pd(muon_mom_pd->X(),muon_mom_pd->Y(),muon_mom_pd->Z(),sqrt(muon_mom_pd->Mag()*muon_mom_pd->Mag() + ml*ml));
      std::vector<double> est_nu_e_h7 = ee::GetEnergyEst(plepton_pd,W_lt,*proton_p4_lt,nprot_lt,*pion_p4_lt,npi_lt,*shower_p4_lt,nsh_lt);

      bool sel_h8 = has_muon_pd && sel_lt; 
      TLorentzVector plepton_h8 = muon_contained_pd ? TLorentzVector(muon_mom_pd->X(),muon_mom_pd->Y(),muon_mom_pd->Z(),sqrt(muon_mom_pd->Mag()*muon_mom_pd->Mag() + ml*ml))
        : TLorentzVector(muon_mom_mcs_pd->X(),muon_mom_mcs_pd->Y(),muon_mom_mcs_pd->Z(),sqrt(muon_mom_mcs_pd->Mag()*muon_mom_mcs_pd->Mag() + ml*ml));
      std::vector<double> est_nu_e_h8 = ee::GetEnergyEst(plepton_h8,W_lt,*proton_p4_lt,nprot_lt,*pion_p4_lt,npi_lt,*shower_p4_lt,nsh_lt);

      std::vector<bool> sel_v = {sel_pd,sel_wc,sel_lt,sel_h7,sel_h8};    
      std::vector<bool> cont_muon_v = {muon_contained_pd,muon_contained_wc,muon_contained_lt,muon_contained_pd,muon_contained_pd};
      std::vector<double> W_v = {W_pd,W_wc,W_lt,W_lt,W_lt};    
      std::vector<std::vector<double>*> e_v = {est_nu_e_pd,est_nu_e_wc,est_nu_e_lt,&est_nu_e_h7,&est_nu_e_h8}; 

      int cat = -1;
      if(i_f == 0){
        if(is_signal_t) cat = 0;
        else if(abs(nu_pdg) == 12) cat = 2; 
        else cat = 1;
      }
      if(i_f == 1) cat = 3;  
      if(i_f == 2) cat = 4;  
      if(cat == -1){
        std::cout << "Uncategorised event!" << std::endl;
        continue;
      }

      for(int i_m=0;i_m<kMethMAX;i_m++){
        if(!sel_v.at(i_m)) continue;
        //if(!cont_muon_v.at(i)) continue;
        for(int i_e=0;i_e<ee::kMAX;i_e++){
          h_Selected_RecoE_v_v_v[i_e][i_m][cat]->Fill(e_v.at(i_m)->at(i_e),weights.at(i_f)); 
        }
      }

    }

    f_in->Close();

  }


  gSystem->Exec("mkdir -p Plots/PredWithBackground/"); 
  TLegend* l = new TLegend(0.75,0.75,0.98,0.98);
  TCanvas* c = new TCanvas("c","c");
  l->SetNColumns(3);

  for(int i_m=0;i_m<kMethMAX;i_m++){
    for(int i_e=0;i_e<ee::kMAX;i_e++){
      THStack* hs = new THStack("hs",";Reconstructed Energy;Events");
      for(size_t i_c=0;i_c<categories.size();i_c++){
        h_Selected_RecoE_v_v_v[i_e][i_m][i_c]->SetFillColor(i_c+2);
        l->AddEntry(h_Selected_RecoE_v_v_v[i_e][i_m][i_c],categories.at(i_c).c_str(),"F");
        hs->Add(h_Selected_RecoE_v_v_v[i_e][i_m][i_c]);  
      }

      hs->Draw("HIST");
      l->Draw();
      c->Print(("Plots/PredWithBackground/RecoE_"+ee::estimators_str.at(i_e)+"_"+methods_str.at(i_m)+".png").c_str());
      c->Clear();
      l->Clear();
      delete hs;
    }
  }


}
