#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"


void NueEfficiency(){

  //const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  //const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nue_overlay_surprise_reco2_hist.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  bool check_containment = false;

  TH1D* h_true_electron_mom = new TH1D("h_true_electron_mom",";True Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_true_electron_costheta = new TH1D("h_true_electron_costheta",";True Nue Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_selected_true_electron_mom_pd = new TH1D("h_selected_true_electron_mom_pd",";True Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_electron_costheta_pd = new TH1D("h_selected_true_electron_costheta_pd",";True Nue Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_electron_mom_reco_electron_mom_pd = new TH2D("h_selected_true_electron_mom_reco_electron_mom_pd",";True Nue Momentum (Gev);Reco Nue Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_electron_costheta_reco_electron_costheta_pd = new TH2D("h_selected_true_electron_costheta_reco_electron_costheta_pd",";True Nue Cos(#theta);Reco Nue Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_electron_mom_error_pd = new TH1D("h_electron_mom_error_pd",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_electron_mom_wc = new TH1D("h_selected_true_electron_mom_wc",";True Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_electron_costheta_wc = new TH1D("h_selected_true_electron_costheta_wc",";True Nue Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_electron_mom_reco_electron_mom_wc = new TH2D("h_selected_true_electron_mom_reco_electron_mom_wc",";True Nue Momentum (Gev);Reco Nue Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_electron_costheta_reco_electron_costheta_wc = new TH2D("h_selected_true_electron_costheta_reco_electron_costheta_wc",";True Nue Cos(#theta);Reco Nue Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_electron_mom_error_wc = new TH1D("h_electron_mom_error_wc",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_electron_mom_lt = new TH1D("h_selected_true_electron_mom_lt",";True Nue Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_electron_costheta_lt = new TH1D("h_selected_true_electron_costheta_lt",";True Nue Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_electron_mom_reco_electron_mom_lt = new TH2D("h_selected_true_electron_mom_reco_electron_mom_lt",";True Nue Momentum (Gev);Reco Nue Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_electron_costheta_reco_electron_costheta_lt = new TH2D("h_selected_true_electron_costheta_reco_electron_costheta_lt",";True Nue Cos(#theta);Reco Nue Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_electron_mom_error_lt = new TH1D("h_electron_mom_error_lt",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 50000) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 12) continue; 

    TVector3 plepton(mc_px->at(0),mc_py->at(0),mc_pz->at(0));
    double p = plepton.Mag();
    double theta = plepton.Theta();
    if(p < 0.05) continue;

    h_true_electron_mom->Fill(p);    
    h_true_electron_costheta->Fill(cos(theta));    


    // Wirecell
    int wc_electron = wc::SimpleNueSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT);
    if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ) && wc_electron != -1 && abs(reco_truthMatch_pdg[wc_electron]) == 11){

      TVector3 reco_plepton(reco_startMomentum[wc_electron][0],reco_startMomentum[wc_electron][1],reco_startMomentum[wc_electron][2]);

      h_selected_true_electron_mom_wc->Fill(plepton.Mag());
      h_selected_true_electron_costheta_wc->Fill(plepton.CosTheta());    

      h_selected_true_electron_mom_reco_electron_mom_wc->Fill(plepton.Mag(),reco_plepton.Mag()); 
      h_selected_true_electron_costheta_reco_electron_costheta_wc->Fill(plepton.CosTheta(),reco_plepton.CosTheta());
      h_electron_mom_error_wc->Fill((reco_plepton.Mag()-plepton.Mag())/plepton.Mag()); 

    }

    // Lantern
    int lt_electron = lt::SimpleNueSelection(nShowers,showerIsSecondary,showerPID);
    if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ) && lt_electron != -1 && abs(showerTruePID[lt_electron]) == 11){

      double mom = sqrt(showerRecoE[lt_electron]*showerRecoE[lt_electron]/1e6 + 2*0.106*showerRecoE[lt_electron]/1e3);
      TVector3 reco_electron_mom(mom*showerStartDirX[lt_electron],mom*showerStartDirY[lt_electron],mom*showerStartDirZ[lt_electron]);

      h_selected_true_electron_mom_lt->Fill(plepton.Mag());
      h_selected_true_electron_costheta_lt->Fill(plepton.CosTheta());

      h_selected_true_electron_mom_reco_electron_mom_lt->Fill(plepton.Mag(),reco_electron_mom.Mag());
      h_selected_true_electron_costheta_reco_electron_costheta_lt->Fill(plepton.CosTheta(),reco_electron_mom.CosTheta());
      h_electron_mom_error_lt->Fill((reco_electron_mom.Mag()-plepton.Mag())/plepton.Mag()); 

    }

  }  


  gSystem->Exec("mkdir -p Plots/NueEfficiency/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);

  h_selected_true_electron_mom_pd->Divide(h_true_electron_mom);
  h_selected_true_electron_costheta_pd->Divide(h_true_electron_costheta);
  h_selected_true_electron_mom_wc->Divide(h_true_electron_mom);
  h_selected_true_electron_costheta_wc->Divide(h_true_electron_costheta);
  h_selected_true_electron_mom_lt->Divide(h_true_electron_mom);
  h_selected_true_electron_costheta_lt->Divide(h_true_electron_costheta);

  THStack* hs_electron_mom_eff = new THStack("hs_electron_mom_eff",";True Nue Momentum (GeV);Efficiency");

  h_selected_true_electron_mom_pd->SetLineColor(1);
  h_selected_true_electron_mom_pd->SetLineWidth(2);
  hs_electron_mom_eff->Add(h_selected_true_electron_mom_pd);
  l->AddEntry(h_selected_true_electron_mom_pd,"PD","L");
    
  h_selected_true_electron_mom_wc->SetLineColor(2);
  h_selected_true_electron_mom_wc->SetLineWidth(2);
  hs_electron_mom_eff->Add(h_selected_true_electron_mom_wc);
  l->AddEntry(h_selected_true_electron_mom_wc,"WC","L");

  h_selected_true_electron_mom_lt->SetLineColor(3);
  h_selected_true_electron_mom_lt->SetLineWidth(2);
  hs_electron_mom_eff->Add(h_selected_true_electron_mom_lt);
  l->AddEntry(h_selected_true_electron_mom_lt,"LT","L");

  hs_electron_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/NueEfficiency/MomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_electron_costheta_eff = new THStack("hs_electron_costheta_eff",";True Nue Cos(#theta);Efficiency");

  h_selected_true_electron_costheta_pd->SetLineColor(1);
  h_selected_true_electron_costheta_pd->SetLineWidth(2);
  hs_electron_costheta_eff->Add(h_selected_true_electron_costheta_pd);
  l->AddEntry(h_selected_true_electron_costheta_pd,"PD","L");
    
  h_selected_true_electron_costheta_wc->SetLineColor(2);
  h_selected_true_electron_costheta_wc->SetLineWidth(2);
  hs_electron_costheta_eff->Add(h_selected_true_electron_costheta_wc);
  l->AddEntry(h_selected_true_electron_costheta_wc,"WC","L");

  h_selected_true_electron_costheta_lt->SetLineColor(3);
  h_selected_true_electron_costheta_lt->SetLineWidth(2);
  hs_electron_costheta_eff->Add(h_selected_true_electron_costheta_lt);
  l->AddEntry(h_selected_true_electron_costheta_lt,"LT","L");

  hs_electron_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/NueEfficiency/CosThetaEfficiency.png");
  c->Clear();
  l->Clear();

  h_selected_true_electron_mom_reco_electron_mom_pd->Draw("colz");
  h_selected_true_electron_mom_reco_electron_mom_pd->SetStats(0);
  c->Print("Plots/NueEfficiency/MomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_electron_mom_reco_electron_mom_pd);
  h_selected_true_electron_mom_reco_electron_mom_pd->Draw("colz");
  h_selected_true_electron_mom_reco_electron_mom_pd->SetStats(0);
  c->Print("Plots/NueEfficiency/Normalised_MomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_electron_costheta_reco_electron_costheta_pd->Draw("colz");
  h_selected_true_electron_costheta_reco_electron_costheta_pd->SetStats(0);
  c->Print("Plots/NueEfficiency/CosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_electron_costheta_reco_electron_costheta_pd);
  h_selected_true_electron_costheta_reco_electron_costheta_pd->Draw("colz");
  h_selected_true_electron_costheta_reco_electron_costheta_pd->SetStats(0);
  c->Print("Plots/NueEfficiency/Normalised_CosThetaReconstruction_PD.png");
  c->Clear();

  h_selected_true_electron_mom_reco_electron_mom_wc->Draw("colz");
  h_selected_true_electron_mom_reco_electron_mom_wc->SetStats(0);
  c->Print("Plots/NueEfficiency/MomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_electron_mom_reco_electron_mom_wc);
  h_selected_true_electron_mom_reco_electron_mom_wc->Draw("colz");
  h_selected_true_electron_mom_reco_electron_mom_wc->SetStats(0);
  c->Print("Plots/NueEfficiency/Normalised_MomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_electron_costheta_reco_electron_costheta_wc->Draw("colz");
  h_selected_true_electron_costheta_reco_electron_costheta_wc->SetStats(0);
  c->Print("Plots/NueEfficiency/CosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_electron_costheta_reco_electron_costheta_wc);
  h_selected_true_electron_costheta_reco_electron_costheta_wc->Draw("colz");
  h_selected_true_electron_costheta_reco_electron_costheta_wc->SetStats(0);
  c->Print("Plots/NueEfficiency/Normalised_CosThetaReconstruction_WC.png");
  c->Clear();


  h_selected_true_electron_mom_reco_electron_mom_lt->Draw("colz");
  h_selected_true_electron_mom_reco_electron_mom_lt->SetStats(0);
  c->Print("Plots/NueEfficiency/MomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_electron_mom_reco_electron_mom_lt);
  h_selected_true_electron_mom_reco_electron_mom_lt->Draw("colz");
  h_selected_true_electron_mom_reco_electron_mom_lt->SetStats(0);
  c->Print("Plots/NueEfficiency/Normalised_MomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_electron_costheta_reco_electron_costheta_lt->Draw("colz");
  h_selected_true_electron_costheta_reco_electron_costheta_lt->SetStats(0);
  c->Print("Plots/NueEfficiency/CosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_electron_costheta_reco_electron_costheta_lt);
  h_selected_true_electron_costheta_reco_electron_costheta_lt->Draw("colz");
  h_selected_true_electron_costheta_reco_electron_costheta_lt->SetStats(0);
  c->Print("Plots/NueEfficiency/Normalised_CosThetaReconstruction_LT.png");
  c->Clear();

  THStack* hs_electron_mom_err = new THStack("hs_electron_mom_err",";Nue Momentum (Reco - True)/True;Events");

  h_electron_mom_error_pd->SetLineColor(1);
  h_electron_mom_error_pd->SetLineWidth(2);
  hs_electron_mom_err->Add(h_electron_mom_error_pd);
  l->AddEntry(h_electron_mom_error_pd,"PD","L");

  h_electron_mom_error_wc->SetLineColor(2);
  h_electron_mom_error_wc->SetLineWidth(2);
  hs_electron_mom_err->Add(h_electron_mom_error_wc);
  l->AddEntry(h_electron_mom_error_wc,"WC","L");

  h_electron_mom_error_lt->SetLineColor(3);
  h_electron_mom_error_lt->SetLineWidth(2);
  hs_electron_mom_err->Add(h_electron_mom_error_lt);
  l->AddEntry(h_electron_mom_error_lt,"LT","L");
    
  hs_electron_mom_err->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/NueEfficiency/MomentumError.png");
  c->Clear();
  l->Clear();


}
