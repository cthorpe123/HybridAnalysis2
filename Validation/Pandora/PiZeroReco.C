#include "Funcs.h"
#include "PD_Funcs.h"

using namespace pd;

void PiZeroReco(){

  const double data_POT = 1.5e21;
  TFile* f_in = new TFile("/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Run4b_Overlay.root");
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
  std::vector<int>* backtracked_pdg=0;

  std::vector<float>* shr_px_v=0;
  std::vector<float>* shr_py_v=0;
  std::vector<float>* shr_pz_v=0;
  std::vector<float>* shr_energy_u_v=0;
  std::vector<float>* shr_energy_v_v=0;
  std::vector<float>* shr_energy_y_v=0;

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
  t_in->SetBranchAddress("trk_dir_x_v",&trk_dir_x_v);
  t_in->SetBranchAddress("trk_dir_y_v",&trk_dir_y_v);
  t_in->SetBranchAddress("trk_dir_z_v",&trk_dir_z_v);

  t_in->SetBranchAddress("shr_px_v",&shr_px_v);
  t_in->SetBranchAddress("shr_py_v",&shr_py_v);
  t_in->SetBranchAddress("shr_pz_v",&shr_pz_v);
  //t_in->SetBranchAddress("shr_energy_u_v",&shr_energy_u_v);
  //t_in->SetBranchAddress("shr_energy_v_v",&shr_energy_v_v);
  t_in->SetBranchAddress("shr_energy_y_v",&shr_energy_y_v);

  TH1D* h_pi0_true_mom = new TH1D("h_pi0_true_mom",";Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_selected_pi0_true_mom = new TH1D("h_selected_pi0_true_mom",";Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_pi0_true_costheta = new TH1D("h_pi0_true_costheta",";#pi^{0} Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_pi0_true_costheta = new TH1D("h_selected_pi0_true_costheta",";#pi^{0} Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_diphoton_mass = new TH1D("h_diphoton_mass",";W (GeV);Events",40,0.0,0.6);
  TH2D* h_pi0_mom = new TH2D("h_pi0_mom",";#pi^{0} Momentum (GeV);Reco #pi^{0} Momentum (GeV);Events",40,0.0,1.0,40,0.0,1.0);
  TH2D* h_pi0_costheta = new TH2D("h_pi0_costheta",";#pi^{0} Cos(#theta);Reco #pi^{0} Cos(#theta);Events",40,-1.0,1.0,40,-1.0,1.0);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 10000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);

    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 0) continue; 

    // Calculating the true invariant mass
    TVector3 pi0_mom(0,0,0);
    double leading_pi0_p = 0.0;
    int n_pi0 = 0;
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      if(mc_pdg->at(i_p) == 111){
        double  p = TVector3(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p)).Mag();
        if(p > leading_pi0_p){
          pi0_mom = TVector3(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
          leading_pi0_p = p;
          n_pi0++;
        }
      }
    }

    if(n_pi0 == 0) continue;

    if(!inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z)) continue;

    h_pi0_true_mom->Fill(pi0_mom.Mag());
    h_pi0_true_costheta->Fill(pi0_mom.CosTheta());

    bool has_good_reco_pi0 = true;
    TLorentzVector p4_tot_reco(0,0,0,0);
    int showers = 0;
    for(size_t i_shr=0;i_shr<shr_px_v->size();i_shr++){
      double E = (shr_energy_y_v->at(i_shr))/1e3/0.83; 
      if(E < 0) continue;
      showers++;       
      p4_tot_reco += TLorentzVector(shr_px_v->at(i_shr)*E,shr_py_v->at(i_shr)*E,shr_pz_v->at(i_shr)*E,E);
    } 

   if(showers < 1) continue;

   TVector3 p_reco = p4_tot_reco.Vect();  

   h_pi0_mom->Fill(pi0_mom.Mag(),TVector3(p4_tot_reco.X(),p4_tot_reco.Y(),p4_tot_reco.Z()).Mag());
   h_pi0_costheta->Fill(pi0_mom.CosTheta(),TVector3(p4_tot_reco.X(),p4_tot_reco.Y(),p4_tot_reco.Z()).CosTheta());

   h_selected_pi0_true_mom->Fill(pi0_mom.Mag());
   h_selected_pi0_true_costheta->Fill(pi0_mom.CosTheta());

   if(shr_px_v->size() != 2) continue;
  
   double W_reco = p4_tot_reco.M();
   h_diphoton_mass->Fill(W_reco);
 
  }  

  gSystem->Exec("mkdir -p Plots/PiZero/");
  TCanvas* c = new TCanvas("c","c");

  h_selected_pi0_true_mom->Divide(h_pi0_true_mom);
  h_selected_pi0_true_costheta->Divide(h_pi0_true_costheta);

  h_diphoton_mass->Draw("HIST");
  h_diphoton_mass->SetStats(0);
  c->Print("Plots/PiZero/DiphotonMass.png");
  c->Clear();
 
  h_selected_pi0_true_mom->Draw("HIST");
  h_selected_pi0_true_mom->SetStats(0);
  c->Print("Plots/PiZero/MomentumEfficiency.png");
  c->Clear();
 
  h_selected_pi0_true_costheta->Draw("HIST");
  h_selected_pi0_true_costheta->SetStats(0);
  c->Print("Plots/PiZero/CosThetaEfficiency.png");
  c->Clear();

  h_pi0_mom->Draw("colz");
  h_pi0_mom->SetStats(0);
  c->Print("Plots/PiZero/Pi0Momentum.png");
  c->Clear();

  Normalise(h_pi0_mom);
  h_pi0_mom->Draw("colz");
  h_pi0_mom->SetStats(0);
  c->Print("Plots/PiZero/Normalised_Pi0Momentum.png");
  c->Clear();

  h_pi0_costheta->Draw("colz");
  h_pi0_costheta->SetStats(0);
  c->Print("Plots/PiZero/Pi0Costheta.png");
  c->Clear();

  Normalise(h_pi0_costheta);
  h_pi0_costheta->Draw("colz");
  h_pi0_costheta->SetStats(0);
  c->Print("Plots/PiZero/Normalised_Pi0Costheta.png");
  c->Clear();

  c->Close();

}
