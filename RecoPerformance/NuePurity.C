#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"


void NuePurity(){

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";

  const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  double pot = 7.88E+20;

  const std::string nue_file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nue_overlay_surprise_reco2_hist.root";
  double nue_pot = 1.17842e+23;

  TH1D* h_reco_electron_mom_all_pd = new TH1D("h_true_electron_mom_all_pd",";Reco Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_reco_electron_costheta_all_pd = new TH1D("h_true_electron_costheta_all_pd",";Reco Nue Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_reco_electron_mom_signal_pd = new TH1D("h_true_electron_mom_signal_pd",";Reco Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_reco_electron_costheta_signal_pd = new TH1D("h_true_electron_costheta_signal_pd",";Reco Nue Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_reco_electron_mom_all_wc = new TH1D("h_true_electron_mom_all_wc",";Reco Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_reco_electron_costheta_all_wc = new TH1D("h_true_electron_costheta_all_wc",";Reco Nue Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_reco_electron_mom_signal_wc = new TH1D("h_true_electron_mom_signal_wc",";Reco Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_reco_electron_costheta_signal_wc = new TH1D("h_true_electron_costheta_signal_wc",";Reco Nue Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_reco_electron_mom_all_lt = new TH1D("h_true_electron_mom_all_lt",";Reco Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_reco_electron_costheta_all_lt = new TH1D("h_true_electron_costheta_all_lt",";Reco Nue Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_reco_electron_mom_signal_lt = new TH1D("h_true_electron_mom_signal_lt",";Reco Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_reco_electron_costheta_signal_lt = new TH1D("h_true_electron_costheta_signal_lt",";Reco Nue Cos(#theta);Events",40,-1.0,1.0);


  std::vector<double> weights = {1.0,pot/nue_pot};
  std::vector<std::string> files = {file,nue_file};

  for(size_t i_f=0;i_f<files.size();i_f++){

    TFile* f_in = nullptr;
    TTree* t_in = nullptr;
    LoadTree(in_dir+files.at(i_f),f_in,t_in,false,false);

    for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

      //if(ievent > 50000) break;
      if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

      t_in->GetEntry(ievent);
      if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;

      if(i_f == 0 && nu_pdg == 12 && ccnc == 0) continue;
      if(i_f != 0 && !(nu_pdg == 12 && ccnc == 0)) continue;
      bool signal = nu_pdg == 12 && ccnc == 0;

      TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
      double p = plepton.Mag();
      double theta = plepton.Theta();
      if(p < 0.05) continue;

      // Wirecell
      int wc_electron = wc::SimpleNueSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT);
      if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ) && wc_electron != -1){

        TVector3 reco_plepton(reco_startMomentum[wc_electron][0],reco_startMomentum[wc_electron][1],reco_startMomentum[wc_electron][2]);

        h_reco_electron_mom_all_wc->Fill(plepton.Mag(),weights.at(i_f));
        h_reco_electron_costheta_all_wc->Fill(plepton.CosTheta(),weights.at(i_f));    

        if(signal && abs(reco_truthMatch_pdg[wc_electron]) == 11){
          h_reco_electron_mom_signal_wc->Fill(plepton.Mag(),weights.at(i_f));
          h_reco_electron_costheta_signal_wc->Fill(plepton.CosTheta(),weights.at(i_f));    
        }

      }

      // Lantern
      int lt_electron = lt::SimpleNueSelection(nShowers,showerIsSecondary,showerPID);
      if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ) && lt_electron != -1){

        double mom = sqrt(showerRecoE[lt_electron]*showerRecoE[lt_electron]/1e6 + 2*0.106*showerRecoE[lt_electron]/1e3);
        TVector3 reco_plepton(mom*showerStartDirX[lt_electron],mom*showerStartDirY[lt_electron],mom*showerStartDirZ[lt_electron]);

        h_reco_electron_mom_all_lt->Fill(plepton.Mag(),weights.at(i_f));
        h_reco_electron_costheta_all_lt->Fill(plepton.CosTheta(),weights.at(i_f));    

        if(signal && abs(showerTruePID[lt_electron]) == 11){
          h_reco_electron_mom_signal_lt->Fill(plepton.Mag(),weights.at(i_f));
          h_reco_electron_costheta_signal_lt->Fill(plepton.CosTheta(),weights.at(i_f));    
        }

      }

    }  

    f_in->Close();

  }

  gSystem->Exec("mkdir -p Plots/NueEfficiency/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);

  h_reco_electron_mom_signal_pd->Divide(h_reco_electron_mom_all_pd);
  h_reco_electron_costheta_signal_pd->Divide(h_reco_electron_costheta_all_pd);
  h_reco_electron_mom_signal_wc->Divide(h_reco_electron_mom_all_wc);
  h_reco_electron_costheta_signal_wc->Divide(h_reco_electron_costheta_all_wc);
  h_reco_electron_mom_signal_lt->Divide(h_reco_electron_mom_all_lt);
  h_reco_electron_costheta_signal_lt->Divide(h_reco_electron_costheta_all_lt);

  THStack* hs_mom = new THStack("hs_mom",";Reco Momentum (GeV);Purity");

  h_reco_electron_mom_signal_pd->SetLineColor(1);
  h_reco_electron_mom_signal_pd->SetLineWidth(2);
  hs_mom->Add(h_reco_electron_mom_signal_pd);
  l->AddEntry(h_reco_electron_mom_signal_pd,"PD","L");
  
  h_reco_electron_mom_signal_wc->SetLineColor(2);
  h_reco_electron_mom_signal_wc->SetLineWidth(2);
  hs_mom->Add(h_reco_electron_mom_signal_wc);
  l->AddEntry(h_reco_electron_mom_signal_wc,"WC","L");

  h_reco_electron_mom_signal_lt->SetLineColor(3);
  h_reco_electron_mom_signal_lt->SetLineWidth(2);
  hs_mom->Add(h_reco_electron_mom_signal_lt);
  l->AddEntry(h_reco_electron_mom_signal_lt,"LT","L");

  hs_mom->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/NueEfficiency/MomentumPurity.png");
  c->Clear();
  l->Clear();

  THStack* hs_costheta = new THStack("hs_costheta",";Reco Momentum (GeV);Purity");

  h_reco_electron_costheta_signal_pd->SetLineColor(1);
  h_reco_electron_costheta_signal_pd->SetLineWidth(2);
  hs_costheta->Add(h_reco_electron_costheta_signal_pd);
  l->AddEntry(h_reco_electron_costheta_signal_pd,"PD","L");
  
  h_reco_electron_costheta_signal_wc->SetLineColor(2);
  h_reco_electron_costheta_signal_wc->SetLineWidth(2);
  hs_costheta->Add(h_reco_electron_costheta_signal_wc);
  l->AddEntry(h_reco_electron_costheta_signal_wc,"WC","L");

  h_reco_electron_costheta_signal_lt->SetLineColor(3);
  h_reco_electron_costheta_signal_lt->SetLineWidth(2);
  hs_costheta->Add(h_reco_electron_costheta_signal_lt);
  l->AddEntry(h_reco_electron_costheta_signal_lt,"LT","L");

  hs_costheta->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/NueEfficiency/CosThetaPurity.png");
  c->Clear();
  l->Clear();
  
}
