#include "TTree.h"
#include "TLeafF16.h"
#pragma link C++ class TLeafF16+;
#include "Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

using namespace lt;

void PiZeroPerformance(){

  const std::string file = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/Run4b_Overlay.root";
  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(file,f_in,t_in,false);

  TH1D* h_pi0_true_mom = new TH1D("h_pi0_true_mom",";Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_selected_pi0_true_mom = new TH1D("h_selected_pi0_true_mom",";Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_pi0_true_costheta = new TH1D("h_pi0_true_costheta",";#pi^{0} Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_selected_pi0_true_costheta = new TH1D("h_selected_pi0_true_costheta",";#pi^{0} Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_diphoton_mass = new TH1D("h_diphoton_mass",";W (GeV);Events",40,0.0,0.6);
  TH2D* h_pi0_mom = new TH2D("h_pi0_mom",";#pi^{0} Momentum (GeV);Reco #pi^{0} Momentum (GeV);Events",40,0.0,1.0,40,0.0,1.0);
  TH2D* h_pi0_costheta = new TH2D("h_pi0_costheta",";#pi^{0} Cos(#theta);Reco #pi^{0} Cos(#theta);Events",40,-1.0,1.0,40,-1.0,1.0);

  for(Long64_t ientry=0;ientry<t_in->GetEntries();ientry++){

    t_in->GetEntry(ientry);
    //if(ientry > 20000) break;
    if(ientry % 50000 == 0) std::cout << ientry << "/" << t_in->GetEntries() << std::endl;

    if(trueNuCCNC != 1) continue;
    if(!inActiveTPC(trueVtxX,trueVtxY,trueVtxZ)) continue; 

    int true_pizero = lt::FindLeadingTruePart(nTruePrimParts,truePrimPartPDG,truePrimPartPx,truePrimPartPy,truePrimPartPz,111);
    if(true_pizero == -1) continue;
    TVector3 true_pizero_mom(truePrimPartPx[true_pizero],truePrimPartPy[true_pizero],truePrimPartPz[true_pizero]);

    // Demand event has reco'd vertex
    if(!foundVertex || !inActiveTPC(vtxX,vtxY,vtxZ)) continue;
 
    h_pi0_true_mom->Fill(true_pizero_mom.Mag());
    h_pi0_true_costheta->Fill(true_pizero_mom.CosTheta());  

    int nphotons = 0;
    TLorentzVector p4(0,0,0,0); 
    for(int i=0;i<nShowers;i++){
      //std::cout << "PID = " << showerPID[i] << "  True PID = " << showerTruePID[i] << std::endl;
      if(showerPID[i] == 22){
        p4 += TLorentzVector(showerRecoE[i]*showerStartDirX[i]/1e3,showerRecoE[i]*showerStartDirY[i]/1e3,showerRecoE[i]*showerStartDirZ[i]/1e3,showerRecoE[i]/1e3);
        nphotons++;
      } 
    }    
    if(nphotons < 1) continue;
    h_pi0_mom->Fill(true_pizero_mom.Mag(),p4.Vect().Mag());   
    h_pi0_costheta->Fill(true_pizero_mom.CosTheta(),p4.Vect().CosTheta());   
    h_selected_pi0_true_mom->Fill(true_pizero_mom.Mag());
    h_selected_pi0_true_costheta->Fill(true_pizero_mom.CosTheta());  

    if(nphotons != 2) continue;
    h_diphoton_mass->Fill(p4.M());
 
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

