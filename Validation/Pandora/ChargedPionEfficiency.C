#include "Funcs.h"
#include "PD_Funcs.h"

using namespace pd;

void ChargedPionEfficiency(){

  const double data_POT = 1.5e21;
  TFile* f_in = new TFile("/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("MergedNtuple"));
  const double POT = 7.88166e+20; 

  int run,sub,evt;
  std::vector<int>* mc_pdg=0;
  std::vector<double>* mc_E=0;
  std::vector<double>* mc_px=0;
  std::vector<double>* mc_py=0;
  std::vector<double>* mc_pz=0;

  float true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z;
  int nu_pdg,ccnc,interaction;
  float nu_e; 

  Float_t         reco_nu_vtx_x;
  Float_t         reco_nu_vtx_y;
  Float_t         reco_nu_vtx_z;
  std::vector<float>* trk_len_v=0;
  std::vector<float>* trk_dir_x_v=0;
  std::vector<float>* trk_dir_y_v=0;
  std::vector<float>* trk_dir_z_v=0;
  std::vector<float>* trk_range_muon_mom_v=0;
  std::vector<float>* trk_llr_pid_score_v=0;
  std::vector<float>* trk_start_x_v=0;
  std::vector<float>* trk_start_y_v=0;
  std::vector<float>* trk_start_z_v=0;
  std::vector<int>* backtracked_pdg=0;

  t_in->SetBranchAddress("run",&run);
  t_in->SetBranchAddress("sub",&sub);
  t_in->SetBranchAddress("evt",&evt);

  t_in->SetBranchAddress("mc_pdg",&mc_pdg);
  t_in->SetBranchAddress("mc_E",&mc_E);
  t_in->SetBranchAddress("mc_px",&mc_px);
  t_in->SetBranchAddress("mc_py",&mc_py);
  t_in->SetBranchAddress("mc_pz",&mc_pz);

  t_in->SetBranchAddress("true_nu_vtx_x",&true_nu_vtx_x);
  t_in->SetBranchAddress("true_nu_vtx_y",&true_nu_vtx_y);
  t_in->SetBranchAddress("true_nu_vtx_z",&true_nu_vtx_z);
  t_in->SetBranchAddress("nu_pdg",&nu_pdg);
  t_in->SetBranchAddress("interaction",&interaction);
  t_in->SetBranchAddress("ccnc",&ccnc);
  t_in->SetBranchAddress("nu_e",&nu_e); 

  t_in->SetBranchAddress("reco_nu_vtx_x", &reco_nu_vtx_x);
  t_in->SetBranchAddress("reco_nu_vtx_y", &reco_nu_vtx_y);
  t_in->SetBranchAddress("reco_nu_vtx_z", &reco_nu_vtx_z);
  t_in->SetBranchAddress("backtracked_pdg",&backtracked_pdg);
  t_in->SetBranchAddress("trk_len_v",&trk_len_v);
  t_in->SetBranchAddress("trk_range_muon_mom_v",&trk_range_muon_mom_v);
  t_in->SetBranchAddress("trk_dir_x_v",&trk_dir_x_v);
  t_in->SetBranchAddress("trk_dir_y_v",&trk_dir_y_v);
  t_in->SetBranchAddress("trk_dir_z_v",&trk_dir_z_v);
  t_in->SetBranchAddress("trk_start_x_v",&trk_start_x_v);
  t_in->SetBranchAddress("trk_start_y_v",&trk_start_y_v);
  t_in->SetBranchAddress("trk_start_z_v",&trk_start_z_v);
  t_in->SetBranchAddress("trk_llr_pid_score_v",&trk_llr_pid_score_v);

  TH1D* h_true_pion_mom = new TH1D("h_true_pion_mom",";True Pion Momentum (GeV);Events",40,0.0,2.5);
  TH1D* h_selected_true_pion_mom = new TH1D("h_selected_true_pion_mom",";True Pion Momentum (GeV);Events",40,0.0,2.5);

  TH1D* h_true_pion_costheta = new TH1D("h_true_pion_costheta",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_true_pion_costheta = new TH1D("h_selected_true_pion_costheta",";True Pion Cos(#theta);Events",40,-1.0,1.0);

  TH2D* h_selected_true_pion_mom_reco_pion_mom = new TH2D("h_selected_true_pion_mom_reco_pion_mom",";True Pion Momentum (Gev);Reco Pion Momentum (Gev);",40,0.0,2.5,40,0.0,2.5);
  TH2D* h_selected_true_pion_costheta_reco_pion_costheta = new TH2D("h_selected_true_pion_costheta_reco_pion_costheta",";True Pion Cos(#theta);Reco Pion Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);

  TH1D* h_pion_mom_error = new TH1D("h_pion_mom_error",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 1000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);

    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 
 
    double leading_pion_mom = 0.0;
    double leading_pion_costheta = 0.0;
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      TVector3 mom(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
      if(abs(mc_pdg->at(i_p)) == 211 && mom.Mag() > leading_pion_mom){
        leading_pion_mom = mom.Mag();
        leading_pion_costheta = mom.CosTheta();
      }
    } 

    // If leading proton has less than 0.3 GeV, event is not signal
    if(leading_pion_mom < 0.1) continue;

    if(!inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z)) continue;
    int muon = SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    if(muon == -1 || abs(backtracked_pdg->at(muon)) != 13) continue;
    
    h_true_pion_mom->Fill(leading_pion_mom); 
    h_true_pion_costheta->Fill(leading_pion_costheta); 

    // Single proton selection
    int pion = SimpleChargedPionSelection(trk_llr_pid_score_v,trk_len_v,muon);
    if(pion == -1 || abs(backtracked_pdg->at(pion)) != 211) continue;

    double reco_p = PionMom(trk_len_v->at(pion)); 
    double reco_costheta = TVector3(trk_dir_x_v->at(pion),trk_dir_y_v->at(pion),trk_dir_z_v->at(pion)).CosTheta();

    h_selected_true_pion_mom->Fill(leading_pion_mom);
    h_selected_true_pion_costheta->Fill(leading_pion_costheta);

    h_selected_true_pion_mom_reco_pion_mom->Fill(leading_pion_mom,reco_p); 
    h_selected_true_pion_costheta_reco_pion_costheta->Fill(leading_pion_costheta,reco_costheta);

    h_pion_mom_error->Fill((reco_p-leading_pion_mom)/leading_pion_mom);
   
  }  

  gSystem->Exec("mkdir -p Plots/PionEfficiency/");
  TCanvas* c = new TCanvas("c","c");

  h_selected_true_pion_mom->Divide(h_true_pion_mom);
  h_selected_true_pion_costheta->Divide(h_true_pion_costheta);

  h_selected_true_pion_mom->Draw("HIST");
  h_selected_true_pion_mom->SetStats(0);
  c->Print("Plots/PionEfficiency/MomentumEfficiency.png");
  c->Clear();

  h_selected_true_pion_costheta->Draw("HIST");
  h_selected_true_pion_costheta->SetStats(0);
  c->Print("Plots/PionEfficiency/CosThetaEfficiency.png");
  c->Clear();

  h_selected_true_pion_mom_reco_pion_mom->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom->SetStats(0);
  c->Print("Plots/PionEfficiency/MomentumReconstruction.png");
  c->Clear();

  Normalise(h_selected_true_pion_mom_reco_pion_mom); 
  h_selected_true_pion_mom_reco_pion_mom->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_MomentumReconstruction.png");
  c->Clear();

  h_selected_true_pion_costheta_reco_pion_costheta->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta->SetStats(0);
  c->Print("Plots/PionEfficiency/CosThetaReconstruction.png");
  c->Clear();

  Normalise(h_selected_true_pion_costheta_reco_pion_costheta);
  h_selected_true_pion_costheta_reco_pion_costheta->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_CosThetaReconstruction.png");
  c->Clear();
 
  h_pion_mom_error->Draw("HIST");
  h_pion_mom_error->SetStats(0);
  c->Print("Plots/PionEfficiency/MomentumError.png");
  c->Clear();

}
