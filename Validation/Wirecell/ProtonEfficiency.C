#include "Funcs.h"
#include "WC_Funcs.h"

using namespace wc;

void ProtonEfficiency(){

  const double data_POT = 1.5e21;
  TFile* f_in = new TFile("/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root");
  TTree* t_in = static_cast<TTree*>(f_in->Get("MergedNtuple"));
  const double POT = 7.88166e+20; 

  Int_t run,subrun,event;
  Int_t mc_nu_pdg;
  Int_t mc_nu_ccnc;
  Int_t mc_nu_mode;
  Float_t mc_nu_pos[4];
  Float_t mc_nu_mom[4];

  Int_t truth_Ntrack;
  Int_t truth_pdg[2000];   //[truth_Ntrack]
  Int_t truth_mother[2000];   //[truth_Ntrack]
  Float_t truth_startMomentum[2000][4];   //[truth_Ntrack]
  Float_t truth_endMomentum[2000][4];   //[truth_Ntrack]
  Int_t truth_id[2000];   //[truth_Ntrack]

  Float_t reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ;
  Int_t reco_Ntrack;
  Int_t reco_pdg[1000];   //[reco_Ntrack]
  Int_t reco_mother[1000];   //[reco_Ntrack]
  Float_t reco_startMomentum[1000][4];   //[reco_Ntrack]
  Float_t reco_endMomentum[1000][4];   //[reco_Ntrack]
  Float_t reco_endXYZT[1000][4];   //[reco_Ntrack]
  Int_t reco_id[1000];   //[reco_Ntrack]
  Int_t reco_truthMatch_pdg[1000];   //[reco_Ntrack]
  Int_t reco_truthMatch_id[1000];   //[reco_Ntrack]

  t_in->SetBranchAddress("run",&run);
  t_in->SetBranchAddress("sub",&subrun);
  t_in->SetBranchAddress("evt",&event);
 
  t_in->SetBranchAddress("mc_nu_pdg",&mc_nu_pdg); 
  t_in->SetBranchAddress("mc_nu_ccnc",&mc_nu_ccnc); 
  t_in->SetBranchAddress("mc_nu_mode",&mc_nu_mode); 
  t_in->SetBranchAddress("mc_nu_pos",mc_nu_pos); 
  t_in->SetBranchAddress("mc_nu_mom",mc_nu_mom); 
 
  t_in->SetBranchAddress("truth_Ntrack",&truth_Ntrack); 
  t_in->SetBranchAddress("truth_pdg",truth_pdg); 
  t_in->SetBranchAddress("truth_mother",truth_mother); 
  t_in->SetBranchAddress("truth_startMomentum",truth_startMomentum); 
  t_in->SetBranchAddress("truth_endMomentum",truth_endMomentum); 
  t_in->SetBranchAddress("truth_id",&truth_id);

  t_in->SetBranchAddress("reco_nuvtxX",&reco_nuvtxX); 
  t_in->SetBranchAddress("reco_nuvtxY",&reco_nuvtxY); 
  t_in->SetBranchAddress("reco_nuvtxZ",&reco_nuvtxZ); 
  t_in->SetBranchAddress("reco_Ntrack",&reco_Ntrack); 
  t_in->SetBranchAddress("reco_pdg",reco_pdg); 
  t_in->SetBranchAddress("reco_mother",reco_mother); 
  t_in->SetBranchAddress("reco_startMomentum",reco_startMomentum); 
  t_in->SetBranchAddress("reco_endMomentum",reco_endMomentum); 
  t_in->SetBranchAddress("reco_endXYZT",reco_endXYZT); 
  t_in->SetBranchAddress("reco_id",&reco_id);
  t_in->SetBranchAddress("reco_truthMatch_pdg",reco_truthMatch_pdg);
  t_in->SetBranchAddress("reco_truthMatch_id",reco_truthMatch_id);

  TH1D* h_true_proton_mom = new TH1D("h_true_proton_mom",";True Proton Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_proton_mom = new TH1D("h_selected_true_proton_mom",";True Proton Momentum (GeV);Events",40,0.0,2.0);

  TH1D* h_true_proton_costheta = new TH1D("h_true_proton_costheta",";True Proton Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_true_proton_costheta = new TH1D("h_selected_true_proton_costheta",";True Proton Cos(#theta);Events",40,-1.0,1.0);

  TH2D* h_selected_true_proton_mom_reco_proton_mom = new TH2D("h_selected_true_proton_mom_reco_proton_mom",";True Proton Momentum (Gev);Reco Proton Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_proton_costheta_reco_proton_costheta = new TH2D("h_selected_true_proton_costheta_reco_proton_costheta",";True Proton Cos(#theta);Reco Proton Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);

  TH1D* h_proton_mom_error = new TH1D("h_proton_mom_error",";(Reco - True)/True;Events",100,-1,1);
 
  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 10000) break;
    if(ievent % 50000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 
    t_in->GetEntry(ievent);   
    
    if(!inActiveTPC(mc_nu_pos[0],mc_nu_pos[1],mc_nu_pos[2])) continue;
    if(mc_nu_ccnc == 1 || mc_nu_pdg != 14) continue; 

    int true_proton = SimpleProtonSelection(truth_Ntrack,truth_pdg,truth_mother,truth_startMomentum);
    if(true_proton == -1) continue;

    TVector3 pproton(truth_startMomentum[true_proton][0],truth_startMomentum[true_proton][1],truth_startMomentum[true_proton][2]);
    double p = pproton.Mag();
    double costheta = pproton.CosTheta();

    h_true_proton_mom->Fill(p); 
    h_true_proton_costheta->Fill(costheta); 

    if(!inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ)) continue;

    int reco_proton = SimpleProtonSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum);
    if(reco_proton == -1 || abs(reco_truthMatch_pdg[reco_proton]) != 2212) continue;

    TVector3 reco_pproton(reco_startMomentum[reco_proton][0],reco_startMomentum[reco_proton][1],reco_startMomentum[reco_proton][2]);
    double reco_p = reco_pproton.Mag();
    double reco_costheta = reco_pproton.CosTheta();

    h_selected_true_proton_mom->Fill(p);
    h_selected_true_proton_costheta->Fill(costheta);

    h_selected_true_proton_mom_reco_proton_mom->Fill(p,reco_p); 
    h_selected_true_proton_costheta_reco_proton_costheta->Fill(costheta,reco_costheta);

    h_proton_mom_error->Fill((reco_p-p)/p);  

  }

  std::cout << "Before Reco" << std::endl;
  for(int i=1;i<h_true_proton_mom->GetNbinsX()+1;i++) std::cout << h_true_proton_mom->GetBinContent(i) << std::endl; 
  std::cout << "After Reco" << std::endl;
  for(int i=1;i<h_selected_true_proton_mom->GetNbinsX()+1;i++) std::cout << h_selected_true_proton_mom->GetBinContent(i) << std::endl; 

 
  gSystem->Exec("mkdir -p Plots/ProtonEfficiency/");
  TCanvas* c = new TCanvas("c","c");

  h_selected_true_proton_mom->Divide(h_true_proton_mom);
  h_selected_true_proton_costheta->Divide(h_true_proton_costheta);

  h_selected_true_proton_mom->Draw("HIST");
  h_selected_true_proton_mom->SetStats(0);
  c->Print("Plots/ProtonEfficiency/MomentumEfficiency.png");
  c->Clear();

  h_selected_true_proton_costheta->Draw("HIST");
  h_selected_true_proton_costheta->SetStats(0);
  c->Print("Plots/ProtonEfficiency/CosThetaEfficiency.png");
  c->Clear();

  h_selected_true_proton_mom_reco_proton_mom->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom->SetStats(0);
  c->Print("Plots/ProtonEfficiency/MomentumReconstruction.png");
  c->Clear();

  Normalise(h_selected_true_proton_mom_reco_proton_mom); 
  h_selected_true_proton_mom_reco_proton_mom->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_MomentumReconstruction.png");
  c->Clear();

  h_selected_true_proton_costheta_reco_proton_costheta->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta->SetStats(0);
  c->Print("Plots/ProtonEfficiency/CosThetaReconstruction.png");
  c->Clear();

  Normalise(h_selected_true_proton_costheta_reco_proton_costheta);
  h_selected_true_proton_costheta_reco_proton_costheta->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_CosThetaReconstruction.png");
  c->Clear();
 
  h_proton_mom_error->Draw("HIST");
  h_proton_mom_error->SetStats(0);
  c->Print("Plots/ProtonEfficiency/MomentumError.png");
  c->Clear();

}
