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
#include "BinningFuncs.h"
#include "WeightFuncs.h"

using namespace syst;
using namespace binning;

// Tune binning so all bins have data FE of this value

void MakeFluxRatios(){

  TH1D* h_NuMu_CV = new TH1D("h_NuMu_CV","#nu_{#mu} True Neutrino Energy [GeV];Events/GeV",100,0,3);
  TH1D* h_NuMuBar_CV = new TH1D("h_NuMuBar_CV","#nu_{#mu} True Neutrino Energy [GeV];Events/GeV",100,0,3);

  std::vector<TH1D*> h_NuMu_FluxUniv(nuniv_Flux,nullptr);
  std::vector<TH1D*> h_NuMuBar_FluxUniv(nuniv_Flux,nullptr);
  for(int i_u=0;i_u<nuniv_Flux;i_u++){
    h_NuMu_FluxUniv.at(i_u) = new TH1D(("h_NuMu_FluxUniv_"+std::to_string(i_u)).c_str(),"#nu_{#mu} True Neutrino Energy [GeV];Events/GeV",100,0,3);
    h_NuMuBar_FluxUniv.at(i_u) = new TH1D(("h_NuMuBar_FluxUniv_"+std::to_string(i_u)).c_str(),"#nu_{#mu} True Neutrino Energy [GeV];Events/GeV",100,0,3);
  }

  std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/retupled/";
  std::vector<std::string> files_v = {

    "run4b/Filtered_Merged_checkout_MCC9.10_Run4b_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist.root"/*,
    "run4c/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4c.root",
    "run4d/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_4d.root",
    "run5/Filtered_Merged_checkout_MCC9.10_Run4acd5_v10_04_07_20_BNB_nu_overlay_retuple_retuple_hist_5.root",*/

  };

  for(int i_f=0;i_f<files_v.size();i_f++){

    std::string file = in_dir + files_v.at(i_f);

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    bool is_overlay,load_syst;
    LoadTreeFiltered(file,f_in,t_in,is_overlay,load_syst);

    for(int ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 10000) break;
      if (ievent % 5000 == 0) std::cout << "  " << ievent << " / " << t_in->GetEntries() << "\r" << std::flush;
      t_in->GetEntry(ievent);

      if(nu_pdg == 14) h_NuMu_CV->Fill(nu_e);
      if(nu_pdg == -14) h_NuMuBar_CV->Fill(nu_e);

      if(load_syst){
        for(int i_u=0;i_u<nuniv_Flux;i_u++){
          double w = weightsFlux->at(i_u)/1000.0;
          if(nu_pdg == 14) h_NuMu_FluxUniv.at(i_u)->Fill(nu_e,w);
          if(nu_pdg == -14) h_NuMuBar_FluxUniv.at(i_u)->Fill(nu_e,w);
        }
      }

    }

   f_in->Close();

  }

  // Ratio histograms (universe shape / CV shape)
  std::vector<TH1D*> h_NuMu_FluxRatio(nuniv_Flux,nullptr);
  std::vector<TH1D*> h_NuMuBar_FluxRatio(nuniv_Flux,nullptr);
  for(int i_u=0;i_u<nuniv_Flux;i_u++){
    h_NuMu_FluxRatio.at(i_u) = static_cast<TH1D*>(h_NuMu_FluxUniv.at(i_u)->Clone(("h_NuMu_FluxRatio_"+std::to_string(i_u)).c_str()));
    h_NuMu_FluxRatio.at(i_u)->Divide(h_NuMu_CV);
    h_NuMuBar_FluxRatio.at(i_u) = static_cast<TH1D*>(h_NuMuBar_FluxUniv.at(i_u)->Clone(("h_NuMuBar_FluxRatio_"+std::to_string(i_u)).c_str()));
    h_NuMuBar_FluxRatio.at(i_u)->Divide(h_NuMuBar_CV);
  }

  // Integral ratios (one entry per universe)
  TH1D* h_NuMu_IntegralRatio = new TH1D("h_NuMu_IntegralRatio","#nu_{#mu} Flux Universe Integral Ratio;Universe;Integral / CV Integral",nuniv_Flux,0,nuniv_Flux);
  TH1D* h_NuMuBar_IntegralRatio = new TH1D("h_NuMuBar_IntegralRatio","#bar{#nu}_{#mu} Flux Universe Integral Ratio;Universe;Integral / CV Integral",nuniv_Flux,0,nuniv_Flux);
  double cv_numu_integral = h_NuMu_CV->Integral();
  double cv_numubar_integral = h_NuMuBar_CV->Integral();
  for(int i_u=0;i_u<nuniv_Flux;i_u++){
    h_NuMu_IntegralRatio->SetBinContent(i_u+1, h_NuMu_FluxUniv.at(i_u)->Integral() / cv_numu_integral);
    h_NuMuBar_IntegralRatio->SetBinContent(i_u+1, h_NuMuBar_FluxUniv.at(i_u)->Integral() / cv_numubar_integral);
  }

  TFile* f_out = new TFile("FluxRatios.root","RECREATE");

  TDirectory* d_cv      = f_out->mkdir("CV");
  TDirectory* d_univ    = f_out->mkdir("Universes");
  TDirectory* d_univ_numu    = d_univ->mkdir("NuMu");
  TDirectory* d_univ_numubar = d_univ->mkdir("NuMuBar");
  TDirectory* d_ratio   = f_out->mkdir("ShapeRatios");
  TDirectory* d_ratio_numu    = d_ratio->mkdir("NuMu");
  TDirectory* d_ratio_numubar = d_ratio->mkdir("NuMuBar");
  TDirectory* d_int     = f_out->mkdir("IntegralRatios");

  d_cv->cd();
  h_NuMu_CV->Write();
  h_NuMuBar_CV->Write();

  d_univ_numu->cd();
  for(int i_u=0;i_u<nuniv_Flux;i_u++) h_NuMu_FluxUniv.at(i_u)->Write();

  d_univ_numubar->cd();
  for(int i_u=0;i_u<nuniv_Flux;i_u++) h_NuMuBar_FluxUniv.at(i_u)->Write();

  d_ratio_numu->cd();
  for(int i_u=0;i_u<nuniv_Flux;i_u++) h_NuMu_FluxRatio.at(i_u)->Write();

  d_ratio_numubar->cd();
  for(int i_u=0;i_u<nuniv_Flux;i_u++) h_NuMuBar_FluxRatio.at(i_u)->Write();

  d_int->cd();
  h_NuMu_IntegralRatio->Write();
  h_NuMuBar_IntegralRatio->Write();

  f_out->Close();

}
