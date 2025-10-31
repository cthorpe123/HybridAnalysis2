#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

using namespace lt;

void PionPerformance(){

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(file,f_in,t_in,false);

  TH1D* h_true_pion_mom = new TH1D("h_true_pion_mom",";True Pion Momentum (GeV);Events",40,0.0,2.5);
  TH1D* h_selected_true_pion_mom = new TH1D("h_selected_true_pion_mom",";True Pion Momentum (GeV);Events",40,0.0,2.5);

  TH1D* h_true_pion_costheta = new TH1D("h_true_pion_costheta",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_true_pion_costheta = new TH1D("h_selected_true_pion_costheta",";True Pion Cos(#theta);Events",40,-1.0,1.0);

  TH2D* h_selected_true_pion_mom_reco_pion_mom = new TH2D("h_selected_true_pion_mom_reco_pion_mom",";True Pion Momentum (Gev);Reco Pion Momentum (Gev);",40,0.0,2.5,40,0.0,2.5);
  TH2D* h_selected_true_pion_costheta_reco_pion_costheta = new TH2D("h_selected_true_pion_costheta_reco_pion_costheta",";True Pion Cos(#theta);Reco Pion Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);

  TH1D* h_pion_mom_error = new TH1D("h_pion_mom_error",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ientry=0;ientry<t_in->GetEntries();ientry++){

    t_in->GetEntry(ientry);
    //if(ientry > 20000) break;
    if(ientry % 50000 == 0) std::cout << ientry << "/" << t_in->GetEntries() << std::endl;

    if(abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    if(!inActiveTPC(trueVtxX,trueVtxY,trueVtxZ)) continue; 

    int true_pion = lt::FindLeadingTruePart(nTruePrimParts,truePrimPartPDG,truePrimPartPx,truePrimPartPy,truePrimPartPz,211);
    if(true_pion == -1) continue;
    TVector3 true_pion_mom(truePrimPartPx[true_pion],truePrimPartPy[true_pion],truePrimPartPz[true_pion]);

    // Demand event has reco'd vertex
    if(!foundVertex || !inActiveTPC(vtxX,vtxY,vtxZ)) continue;

    // PD and WC plots run muon ID first
    int reco_muon = -1; 
    double leading_reco_muon_e = -1; 
    for(int i=0;i<nTracks;i++){
      if(trackIsSecondary[i]) continue;
      if(abs(trackPID[i]) == 13 && trackRecoE[i] > leading_reco_muon_e){
        reco_muon = i;
        leading_reco_muon_e = trackRecoE[i];
      }
    } 
    if(reco_muon == -1 || abs(trackTruePID[reco_muon]) != 13) continue;

    h_true_pion_mom->Fill(true_pion_mom.Mag());
    h_true_pion_costheta->Fill(true_pion_mom.CosTheta());

    int reco_pion = -1; 
    double leading_reco_pion_e = -1; 
    for(int i=0;i<nTracks;i++){
      if(trackIsSecondary[i]) continue;
      if(abs(trackPID[i]) == 211 && trackRecoE[i] > leading_reco_pion_e){
        reco_pion = i;
        leading_reco_pion_e = trackRecoE[i];
      }
    } 
    if(reco_pion == -1 || abs(trackTruePID[reco_pion]) != 211) continue;

    double mom = sqrt(trackRecoE[reco_pion]*trackRecoE[reco_pion]/1e6 + 2*0.140*trackRecoE[reco_pion]/1e3);
    TVector3 reco_pion_mom(mom*trackStartDirX[reco_pion],mom*trackStartDirY[reco_pion],mom*trackStartDirZ[reco_pion]);

    h_selected_true_pion_mom->Fill(true_pion_mom.Mag());
    h_selected_true_pion_costheta->Fill(true_pion_mom.CosTheta());

    h_selected_true_pion_mom_reco_pion_mom->Fill(true_pion_mom.Mag(),reco_pion_mom.Mag());
    h_selected_true_pion_costheta_reco_pion_costheta->Fill(true_pion_mom.CosTheta(),reco_pion_mom.CosTheta());

    h_pion_mom_error->Fill((reco_pion_mom.Mag()-true_pion_mom.Mag())/true_pion_mom.Mag()); 

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

