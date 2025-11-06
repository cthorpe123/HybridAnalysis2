#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

void PionEfficiency(){

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
  const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  TH1D* h_true_pion_mom_pd = new TH1D("h_true_pion_mom_pd",";True Pion Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_true_pion_costheta_pd = new TH1D("h_true_pion_costheta_pd",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_true_pion_mom_wc = new TH1D("h_true_pion_mom_wc",";True Pion Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_true_pion_costheta_wc = new TH1D("h_true_pion_costheta_wc",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH1D* h_true_pion_mom_lt = new TH1D("h_true_pion_mom_lt",";True Pion Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_true_pion_costheta_lt = new TH1D("h_true_pion_costheta_lt",";True Pion Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_selected_true_pion_mom_pd = new TH1D("h_selected_true_pion_mom_pd",";True Pion Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_selected_true_pion_costheta_pd = new TH1D("h_selected_true_pion_costheta_pd",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_pion_mom_reco_pion_mom_pd = new TH2D("h_selected_true_pion_mom_reco_pion_mom_pd",";True Pion Momentum (Gev);Reco Pion Momentum (Gev);",40,0.0,1.0,40,0.0,1.0);
  TH2D* h_selected_true_pion_costheta_reco_pion_costheta_pd = new TH2D("h_selected_true_pion_costheta_reco_pion_costheta_pd",";True Pion Cos(#theta);Reco Pion Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_pion_mom_error_pd = new TH1D("h_pion_mom_error_pd",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_pion_mom_wc = new TH1D("h_selected_true_pion_mom_wc",";True Pion Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_selected_true_pion_costheta_wc = new TH1D("h_selected_true_pion_costheta_wc",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_pion_mom_reco_pion_mom_wc = new TH2D("h_selected_true_pion_mom_reco_pion_mom_wc",";True Pion Momentum (Gev);Reco Pion Momentum (Gev);",40,0.0,1.0,40,0.0,1.0);
  TH2D* h_selected_true_pion_costheta_reco_pion_costheta_wc = new TH2D("h_selected_true_pion_costheta_reco_pion_costheta_wc",";True Pion Cos(#theta);Reco Pion Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_pion_mom_error_wc = new TH1D("h_pion_mom_error_wc",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_pion_mom_lt = new TH1D("h_selected_true_pion_mom_lt",";True Pion Momentum (GeV);Events",40,0.0,1.0);
  TH1D* h_selected_true_pion_costheta_lt = new TH1D("h_selected_true_pion_costheta_lt",";True Pion Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_pion_mom_reco_pion_mom_lt = new TH2D("h_selected_true_pion_mom_reco_pion_mom_lt",";True Pion Momentum (Gev);Reco Pion Momentum (Gev);",40,0.0,1.0,40,0.0,1.0);
  TH2D* h_selected_true_pion_costheta_reco_pion_costheta_lt = new TH2D("h_selected_true_pion_costheta_reco_pion_costheta_lt",";True Pion Cos(#theta);Reco Pion Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_pion_mom_error_lt = new TH1D("h_pion_mom_error_lt",";(Reco - True)/True;Events",100,-1,1);

  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 100) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 
    double leading_pion_mom = 0.0;
    double leading_pion_costheta = 0.0;
    TVector3 ppion(0,0,0);
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      TVector3 mom(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
      if(abs(mc_pdg->at(i_p)) == 211 && mom.Mag() > leading_pion_mom){
        ppion = mom;
        leading_pion_mom = mom.Mag();
        leading_pion_costheta = mom.CosTheta();
      }
    } 

    // If leading pion has less than 0.3 GeV, event is not signal
    if(leading_pion_mom < 0.1) continue;

    // Pandora
    int pd_muon = pd::SimpleMuonSelection(trk_llr_pid_score_v,trk_len_v);
    if(inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z) && pd_muon != -1 && abs(backtracked_pdg->at(pd_muon)) == 13){

      h_true_pion_mom_pd->Fill(leading_pion_mom); 
      h_true_pion_costheta_pd->Fill(leading_pion_costheta); 

      int pd_pion = pd::SimpleChargedPionSelection(trk_llr_pid_score_v,trk_len_v,pd_muon);
      if(pd_pion != -1 && abs(backtracked_pdg->at(pd_pion)) == 211){

        double reco_p = pd::PionMom(trk_len_v->at(pd_pion)); 
        double reco_costheta = TVector3(trk_dir_x_v->at(pd_pion),trk_dir_y_v->at(pd_pion),trk_dir_z_v->at(pd_pion)).CosTheta();

        h_selected_true_pion_mom_pd->Fill(ppion.Mag());
        h_selected_true_pion_costheta_pd->Fill(ppion.CosTheta());    

        h_selected_true_pion_mom_reco_pion_mom_pd->Fill(ppion.Mag(),reco_p); 
        h_selected_true_pion_costheta_reco_pion_costheta_pd->Fill(ppion.CosTheta(),reco_costheta);
        h_pion_mom_error_pd->Fill((reco_p-ppion.Mag())/ppion.Mag()); 

      }
    }

    // Wirecell
    int wc_muon = wc::SimpleMuonSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum,reco_endXYZT);
    if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ) && wc_muon != -1 && abs(reco_truthMatch_pdg[wc_muon]) == 13){

      h_true_pion_mom_wc->Fill(leading_pion_mom); 
      h_true_pion_costheta_wc->Fill(leading_pion_costheta); 

      int wc_pion = wc::SimplePionSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum);
      if(wc_pion != -1 && abs(reco_truthMatch_pdg[wc_pion]) == 211){

        TVector3 reco_ppion(reco_startMomentum[wc_pion][0],reco_startMomentum[wc_pion][1],reco_startMomentum[wc_pion][2]);

        h_selected_true_pion_mom_wc->Fill(ppion.Mag());
        h_selected_true_pion_costheta_wc->Fill(ppion.CosTheta());    

        h_selected_true_pion_mom_reco_pion_mom_wc->Fill(ppion.Mag(),reco_ppion.Mag()); 
        h_selected_true_pion_costheta_reco_pion_costheta_wc->Fill(ppion.CosTheta(),reco_ppion.CosTheta());
        h_pion_mom_error_wc->Fill((reco_ppion.Mag()-ppion.Mag())/ppion.Mag()); 

      }
    }

    // Lantern
    int lt_muon = lt::SimpleMuonSelection(nTracks,trackIsSecondary,trackPID);
    if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ) && lt_muon != -1 && abs(trackTruePID[lt_muon]) == 13){

      h_true_pion_mom_lt->Fill(leading_pion_mom); 
      h_true_pion_costheta_lt->Fill(leading_pion_costheta); 

      int lt_pion = -1; 
      double leading_reco_pion_e = -1; 
      for(int i=0;i<nTracks;i++){
        if(trackIsSecondary[i]) continue;
        if(abs(trackPID[i]) == 211 && trackRecoE[i] > leading_reco_pion_e){
          lt_pion = i;
          leading_reco_pion_e = trackRecoE[i];
        }
      } 

      if(lt_pion != -1 && abs(trackTruePID[lt_pion]) == 211){

        double mom = sqrt(trackRecoE[lt_pion]*trackRecoE[lt_pion]/1e6 + 2*mpi*trackRecoE[lt_pion]/1e3);
        TVector3 reco_pion_mom(mom*trackStartDirX[lt_pion],mom*trackStartDirY[lt_pion],mom*trackStartDirZ[lt_pion]);

        h_selected_true_pion_mom_lt->Fill(ppion.Mag());
        h_selected_true_pion_costheta_lt->Fill(ppion.CosTheta());

        h_selected_true_pion_mom_reco_pion_mom_lt->Fill(ppion.Mag(),reco_pion_mom.Mag());
        h_selected_true_pion_costheta_reco_pion_costheta_lt->Fill(ppion.CosTheta(),reco_pion_mom.CosTheta());
        h_pion_mom_error_lt->Fill((reco_pion_mom.Mag()-ppion.Mag())/ppion.Mag()); 

      }

    }
  }  

  gSystem->Exec("mkdir -p Plots/PionEfficiency/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);

  h_selected_true_pion_mom_pd->Divide(h_true_pion_mom_pd);
  h_selected_true_pion_costheta_pd->Divide(h_true_pion_costheta_pd);
  h_selected_true_pion_mom_wc->Divide(h_true_pion_mom_wc);
  h_selected_true_pion_costheta_wc->Divide(h_true_pion_costheta_wc);
  h_selected_true_pion_mom_lt->Divide(h_true_pion_mom_lt);
  h_selected_true_pion_costheta_lt->Divide(h_true_pion_costheta_lt);

  THStack* hs_pion_mom_eff = new THStack("hs_pion_mom_eff",";True Pion Momentum (GeV);Efficiency");

  h_selected_true_pion_mom_pd->SetLineColor(1);
  h_selected_true_pion_mom_pd->SetLineWidth(2);
  hs_pion_mom_eff->Add(h_selected_true_pion_mom_pd);
  l->AddEntry(h_selected_true_pion_mom_pd,"PD","L");

  h_selected_true_pion_mom_wc->SetLineColor(2);
  h_selected_true_pion_mom_wc->SetLineWidth(2);
  hs_pion_mom_eff->Add(h_selected_true_pion_mom_wc);
  l->AddEntry(h_selected_true_pion_mom_wc,"WC","L");

  h_selected_true_pion_mom_lt->SetLineColor(3);
  h_selected_true_pion_mom_lt->SetLineWidth(2);
  hs_pion_mom_eff->Add(h_selected_true_pion_mom_lt);
  l->AddEntry(h_selected_true_pion_mom_lt,"LT","L");

  hs_pion_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PionEfficiency/MomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_pion_costheta_eff = new THStack("hs_pion_costheta_eff",";True Pion Cos(#theta);Efficiency");

  h_selected_true_pion_costheta_pd->SetLineColor(1);
  h_selected_true_pion_costheta_pd->SetLineWidth(2);
  hs_pion_costheta_eff->Add(h_selected_true_pion_costheta_pd);
  l->AddEntry(h_selected_true_pion_costheta_pd,"PD","L");

  h_selected_true_pion_costheta_wc->SetLineColor(2);
  h_selected_true_pion_costheta_wc->SetLineWidth(2);
  hs_pion_costheta_eff->Add(h_selected_true_pion_costheta_wc);
  l->AddEntry(h_selected_true_pion_costheta_wc,"WC","L");

  h_selected_true_pion_costheta_lt->SetLineColor(3);
  h_selected_true_pion_costheta_lt->SetLineWidth(2);
  hs_pion_costheta_eff->Add(h_selected_true_pion_costheta_lt);
  l->AddEntry(h_selected_true_pion_costheta_lt,"LT","L");

  hs_pion_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PionEfficiency/CosThetaEfficiency.png");
  c->Clear();
  l->Clear();

  h_selected_true_pion_mom_reco_pion_mom_pd->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom_pd->SetStats(0);
  c->Print("Plots/PionEfficiency/MomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_pion_mom_reco_pion_mom_pd);
  h_selected_true_pion_mom_reco_pion_mom_pd->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom_pd->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_MomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_pion_costheta_reco_pion_costheta_pd->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta_pd->SetStats(0);
  c->Print("Plots/PionEfficiency/CosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_pion_costheta_reco_pion_costheta_pd);
  h_selected_true_pion_costheta_reco_pion_costheta_pd->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta_pd->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_CosThetaReconstruction_PD.png");
  c->Clear();

  h_selected_true_pion_mom_reco_pion_mom_wc->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom_wc->SetStats(0);
  c->Print("Plots/PionEfficiency/MomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_pion_mom_reco_pion_mom_wc);
  h_selected_true_pion_mom_reco_pion_mom_wc->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom_wc->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_MomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_pion_costheta_reco_pion_costheta_wc->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta_wc->SetStats(0);
  c->Print("Plots/PionEfficiency/CosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_pion_costheta_reco_pion_costheta_wc);
  h_selected_true_pion_costheta_reco_pion_costheta_wc->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta_wc->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_CosThetaReconstruction_WC.png");
  c->Clear();


  h_selected_true_pion_mom_reco_pion_mom_lt->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom_lt->SetStats(0);
  c->Print("Plots/PionEfficiency/MomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_pion_mom_reco_pion_mom_lt);
  h_selected_true_pion_mom_reco_pion_mom_lt->Draw("colz");
  h_selected_true_pion_mom_reco_pion_mom_lt->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_MomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_pion_costheta_reco_pion_costheta_lt->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta_lt->SetStats(0);
  c->Print("Plots/PionEfficiency/CosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_pion_costheta_reco_pion_costheta_lt);
  h_selected_true_pion_costheta_reco_pion_costheta_lt->Draw("colz");
  h_selected_true_pion_costheta_reco_pion_costheta_lt->SetStats(0);
  c->Print("Plots/PionEfficiency/Normalised_CosThetaReconstruction_LT.png");
  c->Clear();

  THStack* hs_pion_mom_err = new THStack("hs_pion_mom_err",";Pion Momentum (Reco - True)/True;Events");

  h_pion_mom_error_pd->SetLineColor(1);
  h_pion_mom_error_pd->SetLineWidth(2);
  hs_pion_mom_err->Add(h_pion_mom_error_pd);
  l->AddEntry(h_pion_mom_error_pd,"PD","L");

  h_pion_mom_error_wc->SetLineColor(2);
  h_pion_mom_error_wc->SetLineWidth(2);
  hs_pion_mom_err->Add(h_pion_mom_error_wc);
  l->AddEntry(h_pion_mom_error_wc,"WC","L");

  h_pion_mom_error_lt->SetLineColor(3);
  h_pion_mom_error_lt->SetLineWidth(2);
  hs_pion_mom_err->Add(h_pion_mom_error_lt);
  l->AddEntry(h_pion_mom_error_lt,"LT","L");

  hs_pion_mom_err->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/PionEfficiency/MomentumError.png");
  c->Clear();
  l->Clear();


}
