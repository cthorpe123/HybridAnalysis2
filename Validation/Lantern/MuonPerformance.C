#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

using namespace lt;

void MuonPerformance(){

  bool check_containment = true;

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(file,f_in,t_in,false);

  TH1D* h_true_muon_mom = new TH1D("h_true_muon_mom",";True Muon Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_muon_mom = new TH1D("h_selected_true_muon_mom",";True Muon Momentum (GeV);Events",40,0.0,2.0);

  TH1D* h_true_muon_costheta = new TH1D("h_true_muon_costheta",";True Muon Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_true_muon_costheta = new TH1D("h_selected_true_muon_costheta",";True Muon Cos(#theta);Events",40,-1.0,1.0);

  TH2D* h_selected_true_muon_mom_reco_muon_mom = new TH2D("h_selected_true_muon_mom_reco_muon_mom",";True Muon Momentum (Gev);Reco Muon Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_muon_costheta_reco_muon_costheta = new TH2D("h_selected_true_muon_costheta_reco_muon_costheta",";True Muon Cos(#theta);Reco Muon Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  
  TH1D* h_muon_mom_error = new TH1D("h_muon_mom_error",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ientry=0;ientry<t_in->GetEntries();ientry++){

    t_in->GetEntry(ientry);
    //if(ientry > 20000) break;
    if(ientry % 50000 == 0) std::cout << ientry << "/" << t_in->GetEntries() << std::endl;

    if(abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    if(!inActiveTPC(trueVtxX,trueVtxY,trueVtxZ)) continue; 

    int true_muon = lt::FindLeadingTruePart(nTruePrimParts,truePrimPartPDG,truePrimPartPx,truePrimPartPy,truePrimPartPz,13);
    TVector3 true_muon_mom(truePrimPartPx[true_muon],truePrimPartPy[true_muon],truePrimPartPz[true_muon]);
    if(true_muon == -1) continue;

    h_true_muon_mom->Fill(true_muon_mom.Mag());
    h_true_muon_costheta->Fill(true_muon_mom.CosTheta());

    // Demand event has reco'd vertex
    if(!foundVertex || !inActiveTPC(vtxX,vtxY,vtxZ)) continue;

    int reco_muon = -1; 
    for(int i=0;i<nTracks;i++){
      if(trackIsSecondary[i]) continue;
      if(abs(trackPID[i]) == 13){
        reco_muon = i;
        break;
      }
    } 

    if(reco_muon == -1 || abs(trackTruePID[reco_muon]) != 13) continue;
    if(check_containment && !isContained(trackEndPosX[reco_muon],trackEndPosY[reco_muon],trackEndPosZ[reco_muon])) continue;

    double mom = sqrt(trackRecoE[reco_muon]*trackRecoE[reco_muon]/1e6 + 2*0.106*trackRecoE[reco_muon]/1e3);
    TVector3 reco_muon_mom(mom*trackStartDirX[reco_muon],mom*trackStartDirY[reco_muon],mom*trackStartDirZ[reco_muon]);

    h_selected_true_muon_mom->Fill(true_muon_mom.Mag());
    h_selected_true_muon_costheta->Fill(true_muon_mom.CosTheta());

    h_selected_true_muon_mom_reco_muon_mom->Fill(true_muon_mom.Mag(),reco_muon_mom.Mag());
    h_selected_true_muon_costheta_reco_muon_costheta->Fill(true_muon_mom.CosTheta(),reco_muon_mom.CosTheta());
    h_muon_mom_error->Fill((reco_muon_mom.Mag()-true_muon_mom.Mag())/true_muon_mom.Mag()); 

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

