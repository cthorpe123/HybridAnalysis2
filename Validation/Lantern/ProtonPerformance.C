#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

using namespace lt;

void ProtonPerformance(){

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(file,f_in,t_in,false,false);

  TH1D* h_true_proton_mom = new TH1D("h_true_proton_mom",";True Proton Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_proton_mom = new TH1D("h_selected_true_proton_mom",";True Proton Momentum (GeV);Events",40,0.0,2.0);

  TH1D* h_true_proton_costheta = new TH1D("h_true_proton_costheta",";True Proton Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_true_proton_costheta = new TH1D("h_selected_true_proton_costheta",";True Proton Cos(#theta);Events",40,-1.0,1.0);

  TH2D* h_selected_true_proton_mom_reco_proton_mom = new TH2D("h_selected_true_proton_mom_reco_proton_mom",";True Proton Momentum (Gev);Reco Proton Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_proton_costheta_reco_proton_costheta = new TH2D("h_selected_true_proton_costheta_reco_proton_costheta",";True Proton Cos(#theta);Reco Proton Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);

  TH1D* h_proton_mom_error = new TH1D("h_proton_mom_error",";(Reco - True)/True;Events",100,-1,1);

  int ctr = 0;

  for(Long64_t ientry=0;ientry<t_in->GetEntries();ientry++){

    t_in->GetEntry(ientry);
    if(ientry > 100) break;
    if(ientry % 50000 == 0) std::cout << ientry << "/" << t_in->GetEntries() << std::endl;

    if(abs(trueNuPDG) != 14 || trueNuCCNC != 0) continue;
    if(!inActiveTPC(trueVtxX,trueVtxY,trueVtxZ)) continue; 

    int true_proton = lt::FindLeadingTruePart(nTruePrimParts,truePrimPartPDG,truePrimPartPx,truePrimPartPy,truePrimPartPz,2212);

    if(true_proton == -1) continue;
    TVector3 true_proton_mom(truePrimPartPx[true_proton],truePrimPartPy[true_proton],truePrimPartPz[true_proton]);

    h_true_proton_mom->Fill(true_proton_mom.Mag());
    h_true_proton_costheta->Fill(true_proton_mom.CosTheta());

    // Demand event has reco'd vertex
    if(!foundVertex || !inActiveTPC(vtxX,vtxY,vtxZ)) continue;

     ctr++;

    int reco_proton = -1; 
    double leading_reco_proton_e = -1; 
    for(int i=0;i<nTracks;i++){
      if(trackIsSecondary[i]) continue;
      if(abs(trackPID[i]) == 2212 && trackRecoE[i] > leading_reco_proton_e){
        reco_proton = i;
        leading_reco_proton_e = trackRecoE[i];
      }
    } 

    if(reco_proton == -1) continue;
    if(abs(trackTruePID[reco_proton]) != 2212) continue;
    std::cout << run << " " << subrun <<  " "<< event << " " << reco_proton << std::endl;

    double mom = sqrt(trackRecoE[reco_proton]*trackRecoE[reco_proton]/1e6 + 2*0.938*trackRecoE[reco_proton]/1e3);
    TVector3 reco_proton_mom(mom*trackStartDirX[reco_proton],mom*trackStartDirY[reco_proton],mom*trackStartDirZ[reco_proton]);

    h_selected_true_proton_mom->Fill(true_proton_mom.Mag());
    h_selected_true_proton_costheta->Fill(true_proton_mom.CosTheta());

    h_selected_true_proton_mom_reco_proton_mom->Fill(true_proton_mom.Mag(),reco_proton_mom.Mag());
    h_selected_true_proton_costheta_reco_proton_costheta->Fill(true_proton_mom.CosTheta(),reco_proton_mom.CosTheta());

    h_proton_mom_error->Fill((reco_proton_mom.Mag()-true_proton_mom.Mag())/true_proton_mom.Mag());

  }

  std::cout << ctr << std::endl;

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

