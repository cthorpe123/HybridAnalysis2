#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
//#include "Funcs/Funcs.h"
#include "BranchList.h"

void Validate(){

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(file,f_in,t_in,false);

  TH1D* h_vtx_x_pd_minus_wc = new TH1D("h_vtx_x_pd_minus_wc",";Vertex X PD - WC (cm);Events",100,-10,10);
  TH1D* h_vtx_y_pd_minus_wc = new TH1D("h_vtx_y_pd_minus_wc",";Vertex Y PD - WC (cm);Events",100,-10,10);
  TH1D* h_vtx_z_pd_minus_wc = new TH1D("h_vtx_z_pd_minus_wc",";Vertex Z PD - WC (cm);Events",100,-10,10);

  TH1D* h_vtx_x_pd_minus_lt = new TH1D("h_vtx_x_pd_minus_lt",";Vertex X PD - LT (cm);Events",100,-10,10);
  TH1D* h_vtx_y_pd_minus_lt = new TH1D("h_vtx_y_pd_minus_lt",";Vertex Y PD - LT (cm);Events",100,-10,10);
  TH1D* h_vtx_z_pd_minus_lt = new TH1D("h_vtx_z_pd_minus_lt",";Vertex Z PD - LT (cm);Events",100,-10,10);

  TH1D* h_nu_e_pd_minus_wc = new TH1D("h_nu_e_pd_minus_wc",";True Neutrino Energy PD - WC (GeV);Events",100,-0.5,0.5);
  TH1D* h_nu_e_pd_minus_lt = new TH1D("h_nu_e_pd_minus_lt",";True Neutrino Energy PD - LT (GeV);Events",100,-0.5,0.5);

  TH1D* h_n_prot_pd_minus_wc = new TH1D("h_n_prot_pd_minus_wc",";N Protons PD - WC;Events",11,-5.5,5.5); 
  TH1D* h_n_prot_pd_minus_lt = new TH1D("h_n_prot_pd_minus_lt",";N Protons PD - LT;Events",11,-5.5,5.5); 

  TH1D* h_muon_mom_x_pd_minus_wc = new TH1D("h_muon_mom_x_pd_minus_wc",";Muon X Momentum PD - WC (GeV);Events",100,-0.5,0.5);
  TH1D* h_muon_mom_x_pd_minus_lt = new TH1D("h_muon_mom_x_pd_minus_lt",";Muon X Momentum PD - LT (GeV);Events",100,-0.5,0.5);
  TH1D* h_muon_mom_y_pd_minus_wc = new TH1D("h_muon_mom_y_pd_minus_wc",";Muon Y Momentum PD - WC (GeV);Events",100,-0.5,0.5);
  TH1D* h_muon_mom_y_pd_minus_lt = new TH1D("h_muon_mom_y_pd_minus_lt",";Muon Y Momentum PD - LT (GeV);Events",100,-0.5,0.5);
  TH1D* h_muon_mom_z_pd_minus_wc = new TH1D("h_muon_mom_z_pd_minus_wc",";Muon Z Momentum PD - WC (GeV);Events",100,-0.5,0.5);
  TH1D* h_muon_mom_z_pd_minus_lt = new TH1D("h_muon_mom_z_pd_minus_lt",";Muon Z Momentum PD - LT (GeV);Events",100,-0.5,0.5);

  for(int ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 50000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl;
    t_in->GetEntry(ievent);

    // Check neutrino vertex is in same place event by event 
    h_vtx_x_pd_minus_wc->Fill(true_nu_vtx_x - mc_nu_pos[0]);
    h_vtx_y_pd_minus_wc->Fill(true_nu_vtx_y - mc_nu_pos[1]);
    h_vtx_z_pd_minus_wc->Fill(true_nu_vtx_z - mc_nu_pos[2]);

    h_vtx_x_pd_minus_lt->Fill(true_nu_vtx_x - trueVtxX);
    h_vtx_y_pd_minus_lt->Fill(true_nu_vtx_y - trueVtxY);
    h_vtx_z_pd_minus_lt->Fill(true_nu_vtx_z - trueVtxZ);

    // Check neutrino energy event by event
    h_nu_e_pd_minus_wc->Fill(nu_e - mc_nu_mom[3]); 
    h_nu_e_pd_minus_lt->Fill(nu_e - trueNuE); 

    // Count number of FS protons from each framework
    int n_prot_pd = 0;
    for(int i=0;i<mc_pdg->size();i++) if(mc_pdg->at(i) == 2212) n_prot_pd++;

    int n_prot_wc = 0;
    for(int i=0;i<truth_Ntrack;i++) if(truth_mother[i] == 0 && truth_pdg[i] == 2212) n_prot_wc++;

    int n_prot_lt = 0;
    for(int i=0;i<nTruePrimParts;i++) if(truePrimPartPDG[i] == 2212) n_prot_lt++;

    h_n_prot_pd_minus_wc->Fill(n_prot_pd - n_prot_wc);
    h_n_prot_pd_minus_lt->Fill(n_prot_pd - n_prot_lt);

    // Check lepton momentum
    TVector3 plepton_pd(0,0,0);
    for(int i=0;i<mc_pdg->size();i++) if(abs(mc_pdg->at(i)) == 13) plepton_pd = TVector3(mc_px->at(i),mc_py->at(i),mc_pz->at(i));
     
    TVector3 plepton_wc(0,0,0);
    for(int i=0;i<truth_Ntrack;i++) if(truth_mother[i] == 0 && abs(truth_pdg[i]) == 13) plepton_wc = TVector3(truth_startMomentum[i][0],truth_startMomentum[i][1],truth_startMomentum[i][2]);

    TVector3 plepton_lt(0,0,0);
    for(int i=0;i<nTruePrimParts;i++) if(abs(truePrimPartPDG[i]) == 13) plepton_lt = TVector3(truePrimPartPx[i],truePrimPartPy[i],truePrimPartPz[i]);

    h_muon_mom_x_pd_minus_wc->Fill(plepton_pd.X() - plepton_wc.X());
    h_muon_mom_y_pd_minus_wc->Fill(plepton_pd.Y() - plepton_wc.Y());
    h_muon_mom_z_pd_minus_wc->Fill(plepton_pd.Z() - plepton_wc.Z());
    h_muon_mom_x_pd_minus_lt->Fill(plepton_pd.X() - plepton_lt.X());
    h_muon_mom_y_pd_minus_lt->Fill(plepton_pd.Y() - plepton_lt.Y());
    h_muon_mom_z_pd_minus_lt->Fill(plepton_pd.Z() - plepton_lt.Z());

  }

  gSystem->Exec("mkdir -p Plots/Validation/"); 
  TCanvas* c = new TCanvas("c","c");

  h_vtx_x_pd_minus_wc->Draw("HIST"); 
  h_vtx_x_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/vtx_x_pd_minus_wc.png");
  c->Clear();

  h_vtx_y_pd_minus_wc->Draw("HIST"); 
  h_vtx_y_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/vtx_y_pd_minus_wc.png");
  c->Clear();

  h_vtx_z_pd_minus_wc->Draw("HIST"); 
  h_vtx_z_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/vtx_z_pd_minus_wc.png");
  c->Clear();

  h_vtx_x_pd_minus_lt->Draw("HIST"); 
  h_vtx_x_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/vtx_x_pd_minus_lt.png");
  c->Clear();

  h_vtx_y_pd_minus_lt->Draw("HIST"); 
  h_vtx_y_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/vtx_y_pd_minus_lt.png");
  c->Clear();

  h_vtx_z_pd_minus_lt->Draw("HIST"); 
  h_vtx_z_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/vtx_z_pd_minus_lt.png");
  c->Clear();

  h_nu_e_pd_minus_wc->Draw("HIST");
  h_nu_e_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/nu_e_pd_minus_wc.png");
  c->Clear();

  h_nu_e_pd_minus_lt->Draw("HIST");
  h_nu_e_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/nu_e_pd_minus_lt.png");
  c->Clear();

  h_n_prot_pd_minus_wc->Draw("HIST");
  h_n_prot_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/n_prot_pd_minus_wc.png");
  c->Clear();

  h_n_prot_pd_minus_lt->Draw("HIST");
  h_n_prot_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/n_prot_pd_minus_lt.png");
  c->Clear();
  
  h_muon_mom_x_pd_minus_wc->Draw("HIST");
  h_muon_mom_x_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/muon_mom_x_pd_minus_wc.png");
  c->Clear();
  
  h_muon_mom_y_pd_minus_wc->Draw("HIST");
  h_muon_mom_y_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/muon_mom_y_pd_minus_wc.png");
  c->Clear();

  h_muon_mom_z_pd_minus_wc->Draw("HIST");
  h_muon_mom_z_pd_minus_wc->SetStats(0);
  c->Print("Plots/Validation/muon_mom_z_pd_minus_wc.png");
  c->Clear();

  h_muon_mom_x_pd_minus_lt->Draw("HIST");
  h_muon_mom_x_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/muon_mom_x_pd_minus_lt.png");
  c->Clear();
  
  h_muon_mom_y_pd_minus_lt->Draw("HIST");
  h_muon_mom_y_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/muon_mom_y_pd_minus_lt.png");
  c->Clear();

  h_muon_mom_z_pd_minus_lt->Draw("HIST");
  h_muon_mom_z_pd_minus_lt->SetStats(0);
  c->Print("Plots/Validation/muon_mom_z_pd_minus_lt.png");
  c->Clear();

}
