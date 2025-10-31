#include "Funcs.h"
#include "WC_Funcs.h"

using namespace wc;

void PiZeroReco(){

  const double data_POT = 1.5e21;
  TFile* f_in = new TFile("/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Run4b_Overlay.root");
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

  Float_t reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ;
  Int_t reco_Ntrack;
  Int_t reco_pdg[1000];   //[reco_Ntrack]
  Int_t reco_mother[1000];   //[reco_Ntrack]
  Float_t reco_startMomentum[1000][4];   //[reco_Ntrack]
  Float_t reco_endMomentum[1000][4];   //[reco_Ntrack]
  Int_t reco_id[1000];   //[reco_Ntrack]
  Float_t reco_startXYZT[1000][4];   //[reco_Ntrack]
  Float_t reco_endXYZT[1000][4];   //[reco_Ntrack]
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

  t_in->SetBranchAddress("reco_nuvtxX",&reco_nuvtxX); 
  t_in->SetBranchAddress("reco_nuvtxY",&reco_nuvtxY); 
  t_in->SetBranchAddress("reco_nuvtxZ",&reco_nuvtxZ); 
  t_in->SetBranchAddress("reco_Ntrack",&reco_Ntrack); 
  t_in->SetBranchAddress("reco_pdg",reco_pdg); 
  t_in->SetBranchAddress("reco_mother",reco_mother); 
  t_in->SetBranchAddress("reco_startMomentum",reco_startMomentum); 
  t_in->SetBranchAddress("reco_endMomentum",reco_endMomentum); 
  t_in->SetBranchAddress("reco_id",&reco_id);
  //t_in->SetBranchAddress("reco_startXYZT",reco_startXYZT); 
  t_in->SetBranchAddress("reco_endXYZT",reco_endXYZT); 
  t_in->SetBranchAddress("reco_truthMatch_pdg",reco_truthMatch_pdg);
  t_in->SetBranchAddress("reco_truthMatch_id",reco_truthMatch_id);

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

    if(!inActiveTPC(mc_nu_pos[0],mc_nu_pos[1],mc_nu_pos[2])) continue;
    if(mc_nu_ccnc == 0) continue; 

    int n_pi0 = 0;
    TVector3 pi0_mom(0,0,0);
    double leading_pi0_E = 0.0;
    for(size_t i_tr=0;i_tr<truth_Ntrack;i_tr++){
      if(truth_mother[i_tr] == 0 && truth_pdg[i_tr] == 111 && truth_startMomentum[i_tr][3] > leading_pi0_E){
        pi0_mom = TVector3(truth_startMomentum[i_tr][0],truth_startMomentum[i_tr][1],truth_startMomentum[i_tr][2]);
        leading_pi0_E = truth_startMomentum[i_tr][3];
        n_pi0++;
      }
    }

    if(n_pi0 == 0) continue;

    if(!inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ)) continue;

    h_pi0_true_mom->Fill(pi0_mom.Mag());
    h_pi0_true_costheta->Fill(pi0_mom.CosTheta());

    int reco_pi0 = 0;
    bool has_reco_pi0 = false;
    double leading_reco_pi0_E = 0.0;
    int leading_reco_pi0_id = -1;
    std::map<int,int> part_by_id;
    for(size_t i_tr=0;i_tr<reco_Ntrack;i_tr++){
      part_by_id[reco_id[i_tr]] = reco_pdg[i_tr];
    }

    int showers = 0;
    TLorentzVector p4_tot_reco(0,0,0,0);
    for(size_t i_tr=0;i_tr<reco_Ntrack;i_tr++){
      if(reco_pdg[i_tr] == 22 && reco_mother[i_tr] != 0 && (part_by_id.at(reco_mother[i_tr]) != 2212 && part_by_id.at(reco_mother[i_tr]) != 211 && part_by_id.at(reco_mother[i_tr]) != 11 && part_by_id.at(reco_mother[i_tr]) != 13)){
        showers++;
        p4_tot_reco += TLorentzVector(reco_startMomentum[i_tr][0],reco_startMomentum[i_tr][1],reco_startMomentum[i_tr][2],reco_startMomentum[i_tr][3]);
      }
    }

   if(showers < 1) continue;
   TVector3 p_reco = p4_tot_reco.Vect();  
   h_pi0_mom->Fill(pi0_mom.Mag(),p_reco.Mag());
   h_pi0_costheta->Fill(pi0_mom.CosTheta(),p_reco.CosTheta());

   h_selected_pi0_true_mom->Fill(pi0_mom.Mag());
   h_selected_pi0_true_costheta->Fill(pi0_mom.CosTheta());

   if(showers != 2) continue;
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
