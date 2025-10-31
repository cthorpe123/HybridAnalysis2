#include "Funcs.h"
#include "PD_Funcs.h"

using namespace pd;

void MuonEfficiency(){

  bool check_containment = true;

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
  std::vector<float>* trk_mcs_muon_mom_v=0;
  std::vector<float>* trk_llr_pid_score_v=0;
  std::vector<float>* trk_start_x_v=0;
  std::vector<float>* trk_start_y_v=0;
  std::vector<float>* trk_start_z_v=0;
  std::vector<float>* trk_end_x_v=0;
  std::vector<float>* trk_end_y_v=0;
  std::vector<float>* trk_end_z_v=0;
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
  t_in->SetBranchAddress("trk_mcs_muon_mom_v",&trk_mcs_muon_mom_v);
  t_in->SetBranchAddress("trk_dir_x_v",&trk_dir_x_v);
  t_in->SetBranchAddress("trk_dir_y_v",&trk_dir_y_v);
  t_in->SetBranchAddress("trk_dir_z_v",&trk_dir_z_v);
  t_in->SetBranchAddress("trk_start_x_v",&trk_start_x_v);
  t_in->SetBranchAddress("trk_start_y_v",&trk_start_y_v);
  t_in->SetBranchAddress("trk_start_z_v",&trk_start_z_v);
  t_in->SetBranchAddress("trk_end_x_v",&trk_end_x_v);
  t_in->SetBranchAddress("trk_end_y_v",&trk_end_y_v);
  t_in->SetBranchAddress("trk_end_z_v",&trk_end_z_v);
  t_in->SetBranchAddress("trk_llr_pid_score_v",&trk_llr_pid_score_v);

  TH1D* h_true_muon_mom = new TH1D("h_true_muon_mom",";True Muon Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_muon_mom = new TH1D("h_selected_true_muon_mom",";True Muon Momentum (GeV);Events",40,0.0,2.0);

  TH1D* h_true_muon_costheta = new TH1D("h_true_muon_costheta",";True Muon Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_true_muon_costheta = new TH1D("h_selected_true_muon_costheta",";True Muon Cos(#theta);Events",40,-1.0,1.0);

  TH2D* h_selected_true_muon_mom_reco_muon_mom = new TH2D("h_selected_true_muon_mom_reco_muon_mom",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_muon_costheta_reco_muon_costheta = new TH2D("h_selected_true_muon_costheta_reco_muon_costheta",";True Muon Cos(#theta);Reco Muon Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  
  TH1D* h_muon_mom_error = new TH1D("h_muon_mom_error",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 10000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 

    TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
    double p = plepton.Mag();
    double theta = plepton.Theta();
    if(p < 0.1) continue;

    h_true_muon_mom->Fill(p);    
    h_true_muon_costheta->Fill(cos(theta));    

    // Simple preselection and muon id
    if(!inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z)) continue;
    int muon = SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    if(muon == -1 || abs(backtracked_pdg->at(muon)) != 13) continue;
    if(check_containment && !isContained(trk_end_x_v->at(muon),trk_end_y_v->at(muon),trk_end_z_v->at(muon))) continue;

    double reco_p_range = trk_range_muon_mom_v->at(muon); 
    double reco_p_mcs = trk_mcs_muon_mom_v->at(muon); 
    double reco_costheta = TVector3(trk_dir_x_v->at(muon),trk_dir_y_v->at(muon),trk_dir_z_v->at(muon)).CosTheta();
 
    h_selected_true_muon_mom->Fill(p);
    h_selected_true_muon_costheta->Fill(cos(theta));    

    h_selected_true_muon_mom_reco_muon_mom->Fill(p,reco_p_range); 
    h_selected_true_muon_costheta_reco_muon_costheta->Fill(cos(theta),reco_costheta);
    h_muon_mom_error->Fill((reco_p_range-p)/p); 

  }  

  gSystem->Exec("mkdir -p Plots/MuonEfficiency/");
  TCanvas* c = new TCanvas("c","c");

  h_selected_true_muon_mom->Divide(h_true_muon_mom);
  h_selected_true_muon_costheta->Divide(h_true_muon_costheta);

  h_selected_true_muon_mom->Draw("HIST");
  h_selected_true_muon_mom->SetStats(0);
  c->Print("Plots/MuonEfficiency/MomentumEfficiency.png");
  c->Clear();

  h_selected_true_muon_costheta->Draw("HIST");
  h_selected_true_muon_costheta->SetStats(0);
  c->Print("Plots/MuonEfficiency/CosThetaEfficiency.png");
  c->Clear();

  h_selected_true_muon_mom_reco_muon_mom->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom->SetStats(0);
  c->Print("Plots/MuonEfficiency/MomentumReconstruction.png");
  c->Clear();
 
  Normalise(h_selected_true_muon_mom_reco_muon_mom);
  h_selected_true_muon_mom_reco_muon_mom->Draw("colz");
  h_selected_true_muon_mom_reco_muon_mom->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_MomentumReconstruction.png");
  c->Clear();

  h_selected_true_muon_costheta_reco_muon_costheta->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta->SetStats(0);
  c->Print("Plots/MuonEfficiency/CosThetaReconstruction.png");
  c->Clear();

   Normalise(h_selected_true_muon_costheta_reco_muon_costheta);
  h_selected_true_muon_costheta_reco_muon_costheta->Draw("colz");
  h_selected_true_muon_costheta_reco_muon_costheta->SetStats(0);
  c->Print("Plots/MuonEfficiency/Normalised_CosThetaReconstruction.png");
  c->Clear();
 
  h_muon_mom_error->Draw("HIST");
  h_muon_mom_error->SetStats(0);
  c->Print("Plots/MuonEfficiency/MomentumError.png");
  c->Clear();

}
