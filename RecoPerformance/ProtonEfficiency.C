#include "Funcs.h"
#include "PD_Funcs.h"
#include "WC_Funcs.h"
#include "LT_Funcs.h"
#include "BranchList.h"

void ProtonEfficiency(){

 // const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/";
 // const std::string file = "Merged_MCC9.10_Run4b_v10_04_07_09_BNB_nu_overlay_surprise_reco2_hist.root";

  const std::string in_dir = "/exp/uboone/data/users/cthorpe/DIS/Lanpandircell/test/";
  const std::string file = "Merged_larpid_patch_smart_patch_test10_full_more.root";


  TFile* f_in = nullptr;
  TTree* t_in = nullptr;
  LoadTree(in_dir+file,f_in,t_in,false,false);

  TH1D* h_true_proton_mom = new TH1D("h_true_proton_mom",";True Proton Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_true_proton_costheta = new TH1D("h_true_proton_costheta",";True Proton Cos(#theta);Events",40,-1.0,1.0);

  TH1D* h_selected_true_proton_mom_pd = new TH1D("h_selected_true_proton_mom_pd",";True Proton Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_proton_costheta_pd = new TH1D("h_selected_true_proton_costheta_pd",";True Proton Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_proton_mom_reco_proton_mom_pd = new TH2D("h_selected_true_proton_mom_reco_proton_mom_pd",";True Proton Momentum (Gev);Reco Proton Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_proton_costheta_reco_proton_costheta_pd = new TH2D("h_selected_true_proton_costheta_reco_proton_costheta_pd",";True Proton Cos(#theta);Reco Proton Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_proton_mom_error_pd = new TH1D("h_proton_mom_error_pd",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_proton_mom_wc = new TH1D("h_selected_true_proton_mom_wc",";True Proton Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_proton_costheta_wc = new TH1D("h_selected_true_proton_costheta_wc",";True Proton Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_proton_mom_reco_proton_mom_wc = new TH2D("h_selected_true_proton_mom_reco_proton_mom_wc",";True Proton Momentum (Gev);Reco Proton Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_proton_costheta_reco_proton_costheta_wc = new TH2D("h_selected_true_proton_costheta_reco_proton_costheta_wc",";True Proton Cos(#theta);Reco Proton Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_proton_mom_error_wc = new TH1D("h_proton_mom_error_wc",";(Reco - True)/True;Events",100,-1,1);


  TH1D* h_selected_true_proton_mom_lt = new TH1D("h_selected_true_proton_mom_lt",";True Proton Momentum (GeV);Events",40,0.0,2.0);
  TH1D* h_selected_true_proton_costheta_lt = new TH1D("h_selected_true_proton_costheta_lt",";True Proton Cos(#theta);Events",40,-1.0,1.0);
  TH2D* h_selected_true_proton_mom_reco_proton_mom_lt = new TH2D("h_selected_true_proton_mom_reco_proton_mom_lt",";True Proton Momentum (Gev);Reco Proton Momentum (Gev);",40,0.0,2.0,40,0.0,2.0);
  TH2D* h_selected_true_proton_costheta_reco_proton_costheta_lt = new TH2D("h_selected_true_proton_costheta_reco_proton_costheta_lt",";True Proton Cos(#theta);Reco Proton Cos(#theta);",40,-1.0,1.0,40,-1.0,1.0);
  TH1D* h_proton_mom_error_lt = new TH1D("h_proton_mom_error_lt",";(Reco - True)/True;Events",100,-1,1);
  
  for(Long64_t ievent=0;ievent<t_in->GetEntries();ievent++){

    //if(ievent > 100) break;
    if(ievent % 10000 == 0) std::cout << ievent << "/" << t_in->GetEntries() << std::endl; 

    t_in->GetEntry(ievent);
    if(!inActiveTPC(true_nu_vtx_x,true_nu_vtx_y,true_nu_vtx_z)) continue;
    if(ccnc == 1 || nu_pdg != 14) continue; 
    double leading_proton_mom = 0.0;
    double leading_proton_costheta = 0.0;
    TVector3 pproton(0,0,0);
    for(size_t i_p=0;i_p<mc_pdg->size();i_p++){
      TVector3 mom(mc_px->at(i_p),mc_py->at(i_p),mc_pz->at(i_p));
      if(abs(mc_pdg->at(i_p)) == 2212 && mom.Mag() > leading_proton_mom){
        pproton = mom;
        leading_proton_mom = mom.Mag();
        leading_proton_costheta = mom.CosTheta();
      }
    } 

    // If leading proton has less than 0.3 GeV, event is not signal
    if(leading_proton_mom < 0.3) continue;

    h_true_proton_mom->Fill(leading_proton_mom); 
    h_true_proton_costheta->Fill(leading_proton_costheta); 

    // Pandora
    int pd_proton = pd::SimpleProtonSelection(trk_llr_pid_score_v,trk_len_v);
    if(inActiveTPC(reco_nu_vtx_x,reco_nu_vtx_y,reco_nu_vtx_z) && pd_proton != -1 && abs(backtracked_pdg->at(pd_proton)) == 2212){

      double reco_p = pd::ProtonMom(trk_len_v->at(pd_proton)); 
      double reco_costheta = TVector3(trk_dir_x_v->at(pd_proton),trk_dir_y_v->at(pd_proton),trk_dir_z_v->at(pd_proton)).CosTheta();

      h_selected_true_proton_mom_pd->Fill(pproton.Mag());
      h_selected_true_proton_costheta_pd->Fill(pproton.CosTheta());    

      h_selected_true_proton_mom_reco_proton_mom_pd->Fill(pproton.Mag(),reco_p); 
      h_selected_true_proton_costheta_reco_proton_costheta_pd->Fill(pproton.CosTheta(),reco_costheta);
      h_proton_mom_error_pd->Fill((reco_p-pproton.Mag())/pproton.Mag()); 

    }

    // Wirecell
    int wc_proton = wc::SimpleProtonSelection(reco_Ntrack,reco_pdg,reco_mother,reco_startMomentum);
    if(inActiveTPC(reco_nuvtxX,reco_nuvtxY,reco_nuvtxZ) && wc_proton != -1 && abs(reco_truthMatch_pdg[wc_proton]) == 2212){

      TVector3 reco_pproton(reco_startMomentum[wc_proton][0],reco_startMomentum[wc_proton][1],reco_startMomentum[wc_proton][2]);

      h_selected_true_proton_mom_wc->Fill(pproton.Mag());
      h_selected_true_proton_costheta_wc->Fill(pproton.CosTheta());    

      h_selected_true_proton_mom_reco_proton_mom_wc->Fill(pproton.Mag(),reco_pproton.Mag()); 
      h_selected_true_proton_costheta_reco_proton_costheta_wc->Fill(pproton.CosTheta(),reco_pproton.CosTheta());
      h_proton_mom_error_wc->Fill((reco_pproton.Mag()-pproton.Mag())/pproton.Mag()); 

    }

    // Lantern
    int lt_proton = -1; 
    double leading_reco_proton_e = -1; 
    for(int i=0;i<nTracks;i++){
      if(trackIsSecondary[i]) continue;
      if(abs(trackPID[i]) == 2212 && trackRecoE[i] > leading_reco_proton_e){
        lt_proton = i;
        leading_reco_proton_e = trackRecoE[i];
      }
    } 

    if(foundVertex && inActiveTPC(vtxX,vtxY,vtxZ) && lt_proton != -1 && abs(trackTruePID[lt_proton]) == 2212){

      double mom = sqrt(trackRecoE[lt_proton]*trackRecoE[lt_proton]/1e6 + 2*Mp*trackRecoE[lt_proton]/1e3);
      TVector3 reco_proton_mom(mom*trackStartDirX[lt_proton],mom*trackStartDirY[lt_proton],mom*trackStartDirZ[lt_proton]);

      h_selected_true_proton_mom_lt->Fill(pproton.Mag());
      h_selected_true_proton_costheta_lt->Fill(pproton.CosTheta());

      h_selected_true_proton_mom_reco_proton_mom_lt->Fill(pproton.Mag(),reco_proton_mom.Mag());
      h_selected_true_proton_costheta_reco_proton_costheta_lt->Fill(pproton.CosTheta(),reco_proton_mom.CosTheta());
      h_proton_mom_error_lt->Fill((reco_proton_mom.Mag()-pproton.Mag())/pproton.Mag()); 

    }

  }  

  gSystem->Exec("mkdir -p Plots/ProtonEfficiency/");
  TCanvas* c = new TCanvas("c","c");
  TLegend* l = new TLegend(0.8,0.8,0.95,0.95);

  h_selected_true_proton_mom_pd->Divide(h_true_proton_mom);
  h_selected_true_proton_costheta_pd->Divide(h_true_proton_costheta);
  h_selected_true_proton_mom_wc->Divide(h_true_proton_mom);
  h_selected_true_proton_costheta_wc->Divide(h_true_proton_costheta);
  h_selected_true_proton_mom_lt->Divide(h_true_proton_mom);
  h_selected_true_proton_costheta_lt->Divide(h_true_proton_costheta);

  THStack* hs_proton_mom_eff = new THStack("hs_proton_mom_eff",";True Proton Momentum (GeV);Efficiency");

  h_selected_true_proton_mom_pd->SetLineColor(1);
  h_selected_true_proton_mom_pd->SetLineWidth(2);
  hs_proton_mom_eff->Add(h_selected_true_proton_mom_pd);
  l->AddEntry(h_selected_true_proton_mom_pd,"PD","L");
    
  h_selected_true_proton_mom_wc->SetLineColor(2);
  h_selected_true_proton_mom_wc->SetLineWidth(2);
  hs_proton_mom_eff->Add(h_selected_true_proton_mom_wc);
  l->AddEntry(h_selected_true_proton_mom_wc,"WC","L");

  h_selected_true_proton_mom_lt->SetLineColor(3);
  h_selected_true_proton_mom_lt->SetLineWidth(2);
  hs_proton_mom_eff->Add(h_selected_true_proton_mom_lt);
  l->AddEntry(h_selected_true_proton_mom_lt,"LT","L");

  hs_proton_mom_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/ProtonEfficiency/MomentumEfficiency.png");
  c->Clear();
  l->Clear();

  THStack* hs_proton_costheta_eff = new THStack("hs_proton_costheta_eff",";True Proton Cos(#theta);Efficiency");

  h_selected_true_proton_costheta_pd->SetLineColor(1);
  h_selected_true_proton_costheta_pd->SetLineWidth(2);
  hs_proton_costheta_eff->Add(h_selected_true_proton_costheta_pd);
  l->AddEntry(h_selected_true_proton_costheta_pd,"PD","L");
    
  h_selected_true_proton_costheta_wc->SetLineColor(2);
  h_selected_true_proton_costheta_wc->SetLineWidth(2);
  hs_proton_costheta_eff->Add(h_selected_true_proton_costheta_wc);
  l->AddEntry(h_selected_true_proton_costheta_wc,"WC","L");

  h_selected_true_proton_costheta_lt->SetLineColor(3);
  h_selected_true_proton_costheta_lt->SetLineWidth(2);
  hs_proton_costheta_eff->Add(h_selected_true_proton_costheta_lt);
  l->AddEntry(h_selected_true_proton_costheta_lt,"LT","L");

  hs_proton_costheta_eff->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/ProtonEfficiency/CosThetaEfficiency.png");
  c->Clear();
  l->Clear();

  h_selected_true_proton_mom_reco_proton_mom_pd->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom_pd->SetStats(0);
  c->Print("Plots/ProtonEfficiency/MomentumReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_proton_mom_reco_proton_mom_pd);
  h_selected_true_proton_mom_reco_proton_mom_pd->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom_pd->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_MomentumReconstruction_PD.png");
  c->Clear();

  h_selected_true_proton_costheta_reco_proton_costheta_pd->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta_pd->SetStats(0);
  c->Print("Plots/ProtonEfficiency/CosThetaReconstruction_PD.png");
  c->Clear();

  Normalise(h_selected_true_proton_costheta_reco_proton_costheta_pd);
  h_selected_true_proton_costheta_reco_proton_costheta_pd->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta_pd->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_CosThetaReconstruction_PD.png");
  c->Clear();

  h_selected_true_proton_mom_reco_proton_mom_wc->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom_wc->SetStats(0);
  c->Print("Plots/ProtonEfficiency/MomentumReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_proton_mom_reco_proton_mom_wc);
  h_selected_true_proton_mom_reco_proton_mom_wc->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom_wc->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_MomentumReconstruction_WC.png");
  c->Clear();

  h_selected_true_proton_costheta_reco_proton_costheta_wc->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta_wc->SetStats(0);
  c->Print("Plots/ProtonEfficiency/CosThetaReconstruction_WC.png");
  c->Clear();

  Normalise(h_selected_true_proton_costheta_reco_proton_costheta_wc);
  h_selected_true_proton_costheta_reco_proton_costheta_wc->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta_wc->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_CosThetaReconstruction_WC.png");
  c->Clear();


  h_selected_true_proton_mom_reco_proton_mom_lt->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom_lt->SetStats(0);
  c->Print("Plots/ProtonEfficiency/MomentumReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_proton_mom_reco_proton_mom_lt);
  h_selected_true_proton_mom_reco_proton_mom_lt->Draw("colz");
  h_selected_true_proton_mom_reco_proton_mom_lt->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_MomentumReconstruction_LT.png");
  c->Clear();

  h_selected_true_proton_costheta_reco_proton_costheta_lt->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta_lt->SetStats(0);
  c->Print("Plots/ProtonEfficiency/CosThetaReconstruction_LT.png");
  c->Clear();

  Normalise(h_selected_true_proton_costheta_reco_proton_costheta_lt);
  h_selected_true_proton_costheta_reco_proton_costheta_lt->Draw("colz");
  h_selected_true_proton_costheta_reco_proton_costheta_lt->SetStats(0);
  c->Print("Plots/ProtonEfficiency/Normalised_CosThetaReconstruction_LT.png");
  c->Clear();

  THStack* hs_proton_mom_err = new THStack("hs_proton_mom_err",";Proton Momentum (Reco - True)/True;Events");

  h_proton_mom_error_pd->SetLineColor(1);
  h_proton_mom_error_pd->SetLineWidth(2);
  hs_proton_mom_err->Add(h_proton_mom_error_pd);
  l->AddEntry(h_proton_mom_error_pd,"PD","L");

  h_proton_mom_error_wc->SetLineColor(2);
  h_proton_mom_error_wc->SetLineWidth(2);
  hs_proton_mom_err->Add(h_proton_mom_error_wc);
  l->AddEntry(h_proton_mom_error_wc,"WC","L");

  h_proton_mom_error_lt->SetLineColor(3);
  h_proton_mom_error_lt->SetLineWidth(2);
  hs_proton_mom_err->Add(h_proton_mom_error_lt);
  l->AddEntry(h_proton_mom_error_lt,"LT","L");
    
  hs_proton_mom_err->Draw("nostack HIST");
  l->Draw();
  c->Print("Plots/ProtonEfficiency/MomentumError.png");
  c->Clear();
  l->Clear();


}
